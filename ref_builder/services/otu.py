"""Manage OTU data."""

from uuid import UUID

import structlog

from ref_builder.models.isolate import IsolateName
from ref_builder.models.molecule import Molecule
from ref_builder.models.plan import (
    Plan,
    Segment,
    SegmentRule,
    extract_segment_name_from_record,
)
from ref_builder.ncbi.models import NCBIGenbank, NCBIRank, NCBITaxonomy
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.isolate import create_sequence_from_record
from ref_builder.otu.promote import (
    promote_otu_accessions,
    promote_otu_accessions_from_records,
    upgrade_outdated_sequences_in_otu,
)
from ref_builder.otu.utils import (
    assign_records_to_segments,
    create_segments_from_records,
    get_segments_max_length,
    get_segments_min_length,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.services import Service

logger = structlog.get_logger("services.otu")

UUID_STRING_LENGTH = 36


def create_plan_from_records(
    records: list[NCBIGenbank],
    length_tolerance: float,
    segments: list[Segment] | None = None,
) -> Plan | None:
    """Return a plan from a list of records representing an isolate."""
    if len(records) == 1:
        record = records[0]

        return Plan.new(
            segments=[
                Segment.new(
                    length=len(record.sequence),
                    length_tolerance=length_tolerance,
                    name=extract_segment_name_from_record(record),
                    rule=SegmentRule.REQUIRED,
                )
            ]
        )

    if len(group_genbank_records_by_isolate(records)) > 1:
        logger.warning("More than one isolate found. Cannot create plan.")
        return None

    if segments is None:
        segments = create_segments_from_records(
            records,
            rule=SegmentRule.REQUIRED,
            length_tolerance=length_tolerance,
        )

    if segments is not None:
        return Plan.new(segments=segments)

    return None


def get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
    """Return relevant molecule metadata from one or more records.

    Molecule metadata is retrieved from the first RefSeq record in the list.
    If no RefSeq record is found in the list, molecule metadata is retrieved
    from record[0].
    """
    if not records:
        raise ValueError("No records given")

    # Assign first record as benchmark to start
    representative_record = records[0]

    if not representative_record.refseq:
        for record in records:
            if record.refseq:
                # Replace representative record with first RefSeq record found
                representative_record = record
                break

    return Molecule.model_validate(
        {
            "strandedness": representative_record.strandedness.value,
            "type": representative_record.moltype.value,
            "topology": representative_record.topology.value,
        },
    )


class OTUService(Service):
    """Service for managing OTU operations."""

    def get_otu(self, identifier: str) -> OTUBuilder | None:
        """Return an OTU from the repo if identifier matches a single OTU.

        The identifier can either be a stringified UUIDv4, a NCBI Taxonomy ID,
        or an acronym associated with the OTU.

        :param identifier: a non-UUID identifier.
            Can be an integer Taxonomy ID or acronym.
        :return: the OTU or None if not found
        """
        otu_id = None

        if len(identifier) == UUID_STRING_LENGTH:
            try:
                otu_id = UUID(identifier)
            except ValueError:
                otu_id = None

        elif identifier.isnumeric():
            try:
                taxid = int(identifier)
            except ValueError:
                pass
            else:
                otu_id = self._repo.get_otu_id_by_taxid(taxid)

        else:
            otu_id = self._repo.get_otu_id_by_acronym(identifier)

        if otu_id is None:
            return None

        return self._repo.get_otu(otu_id)

    def create(
        self,
        accessions: list[str],
    ) -> OTUBuilder | None:
        """Create a new OTU from a list of accessions.

        Uses the provided accessions to generate a plan and add a first isolate.
        Derives the taxonomy ID from the accessions.

        :param accessions: accessions to build the new OTU from
        :return: the created OTU or None if creation failed
        """
        otu_logger = logger.bind(accessions=accessions)

        if not accessions:
            otu_logger.error("OTU could not be created to spec based on given data.")
            return None

        records = self.ncbi.fetch_genbank_records(accessions)

        if len(records) != len(accessions):
            otu_logger.fatal("Could not retrieve all requested accessions.")
            return None

        if len({record.source.taxid for record in records}) > 1:
            otu_logger.fatal("Not all records are from the same organism.")
            return None

        taxid = records[0].source.taxid

        binned_records = group_genbank_records_by_isolate(records)

        if len(binned_records) > 1:
            otu_logger.fatal("More than one isolate found. Cannot create plan.")
            return None

        taxonomy = self.ncbi.fetch_taxonomy_record(taxid)

        if taxonomy is None:
            otu_logger.fatal("Could not retrieve data from NCBI Taxonomy", taxid=taxid)
            return None

        if taxonomy.rank != NCBIRank.SPECIES:
            taxonomy = self.ncbi.fetch_taxonomy_record(taxonomy.species.id)

        if taxonomy is None:
            otu_logger.fatal("Could not retrieve data from NCBI Taxonomy", taxid=taxid)
            return None

        if self._repo.get_otu_id_by_taxid(taxonomy.id):
            otu_logger.error(
                f"Taxonomy ID {taxonomy.id} has already been added to this reference."
            )
            return None

        with self._repo.use_transaction():
            return self._write_otu(
                taxonomy,
                records,
                isolate_name=next(iter(binned_records.keys()))
                if binned_records
                else None,
            )

    def _write_otu(
        self,
        taxonomy: NCBITaxonomy,
        records: list[NCBIGenbank],
        isolate_name: IsolateName | None,
    ) -> OTUBuilder | None:
        """Create a new OTU from an NCBI Taxonomy record and a list of Nucleotide records.

        :param taxonomy: the NCBI taxonomy record
        :param records: the GenBank records
        :param isolate_name: the isolate name
        :return: the created OTU or None if creation failed
        """
        otu_logger = logger.bind(taxid=taxonomy.id)

        plan = create_plan_from_records(
            records,
            length_tolerance=self._repo.settings.default_segment_length_tolerance,
        )

        if plan is None:
            otu_logger.fatal("Could not create plan from records.")
            return None

        molecule = get_molecule_from_records(records)

        lineage = self.ncbi.fetch_lineage(records[0].source.taxid)

        otu = self._repo.create_otu(
            lineage=lineage,
            molecule=molecule,
            plan=plan,
        )

        isolate = self._repo.create_isolate(
            otu_id=otu.id,
            name=isolate_name,
            taxid=records[0].source.taxid,
        )

        otu.add_isolate(isolate)

        if otu.plan.monopartite:
            record = records[0]

            sequence = create_sequence_from_record(
                self._repo, otu, record, plan.segments[0].id
            )

            self._repo.link_sequence(otu.id, isolate.id, sequence.id)

            if record.refseq:
                _, old_accession = parse_refseq_comment(record.comment)

                self._repo.exclude_accessions(
                    otu.id,
                    [old_accession],
                )

        else:
            for segment_id, record in assign_records_to_segments(records, plan).items():
                sequence = create_sequence_from_record(
                    self._repo, otu, record, segment_id
                )

                self._repo.link_sequence(otu.id, isolate.id, sequence.id)

                if record.refseq:
                    _, old_accession = parse_refseq_comment(record.comment)
                    self._repo.exclude_accessions(
                        otu.id,
                        [old_accession],
                    )

        return self._repo.get_otu(otu.id)

    def exclude_accessions(
        self,
        otu_id: UUID,
        accessions: set[str],
    ) -> None:
        """Exclude accessions from future addition to an OTU.

        :param otu_id: the OTU ID
        :param accessions: accessions to exclude
        """
        otu = self._repo.get_otu(otu_id)

        if otu is None:
            logger.error("OTU not found", otu_id=str(otu_id))
            return

        original_excluded_accessions = otu.excluded_accessions.copy()

        with self._repo.use_transaction():
            excluded_accessions = self._repo.exclude_accessions(
                otu_id=otu.id, accessions=accessions
            )

        if excluded_accessions == original_excluded_accessions:
            logger.info(
                "Excluded accession list already up to date.",
                excluded_accessions=excluded_accessions,
            )
        else:
            logger.info(
                "Updated excluded accession list.",
                otu_id=str(otu.id),
                excluded_accessions=sorted(excluded_accessions),
            )

    def allow_accessions(
        self,
        otu_id: UUID,
        accessions: set[str],
    ) -> None:
        """Allow accessions for future addition to an OTU.

        This reverses the effect of exclude_accessions.

        :param otu_id: the OTU ID
        :param accessions: accessions to allow
        """
        otu = self._repo.get_otu(otu_id)

        if otu is None:
            logger.error("OTU not found", otu_id=str(otu_id))
            return

        original_excluded_accessions = otu.excluded_accessions.copy()

        with self._repo.use_transaction():
            excluded_accessions = self._repo.allow_accessions(
                otu_id=otu.id, accessions=accessions
            )

        if excluded_accessions == original_excluded_accessions:
            logger.info(
                "Excluded accession list already up to date.",
                excluded_accessions=excluded_accessions,
            )
        else:
            logger.info(
                "Updated excluded accession list.",
                otu_id=str(otu.id),
                excluded_accessions=sorted(excluded_accessions),
            )

    def update(self, otu_id: UUID, ignore_cache: bool = False) -> OTUBuilder | None:
        """Update an OTU by promoting, upgrading, and adding new isolates.

        This method performs three operations in sequence:
        1. Promote GenBank accessions to RefSeq equivalents where available
        2. Upgrade outdated sequence versions (e.g., v1 → v2)
        3. Add new isolates from newly available accessions

        :param otu_id: the OTU ID
        :param ignore_cache: whether to ignore the NCBI cache
        :return: the updated OTU or None if the OTU was not found
        """
        otu = self._repo.get_otu(otu_id)

        if otu is None:
            logger.error("OTU not found", otu_id=str(otu_id))
            return None

        log = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

        log.info("Starting comprehensive OTU update.")

        # Step 1: Promote GenBank accessions to RefSeq
        promoted_accessions = promote_otu_accessions(self._repo, otu, ignore_cache)
        if promoted_accessions:
            log.info("Promoted sequences", count=len(promoted_accessions))
            otu = self._repo.get_otu(otu.id)

        # Step 2: Upgrade outdated sequence versions
        upgraded_sequence_ids = upgrade_outdated_sequences_in_otu(
            self._repo,
            otu,
            modification_date_start=None,
            ignore_cache=ignore_cache,
        )
        if upgraded_sequence_ids:
            log.info("Upgraded sequences", count=len(upgraded_sequence_ids))
            otu = self._repo.get_otu(otu.id)

        # Step 3: Add new isolates
        accessions = self.ncbi.fetch_accessions_by_taxid(
            otu.taxid,
            sequence_min_length=get_segments_min_length(otu.plan.segments),
            sequence_max_length=get_segments_max_length(otu.plan.segments),
        )

        fetch_set = {accession.key for accession in accessions} - otu.blocked_accessions

        if fetch_set:
            log.info("Adding new isolates from NCBI.")
            log.info(
                "Fetching records.",
                count=len(fetch_set),
                fetch_list=sorted(fetch_set),
            )

            records = self.ncbi.fetch_genbank_records(fetch_set)

            if records:
                # Promote RefSeq records first
                refseq_records = [r for r in records if r.refseq]
                if refseq_records and promote_otu_accessions_from_records(
                    self._repo, otu, refseq_records
                ):
                    otu = self._repo.get_otu(otu.id)

                # Create isolates
                new_isolate_ids = []

                for isolate_name, isolate_records in group_genbank_records_by_isolate(
                    records
                ).items():
                    try:
                        isolate = self._services.isolate.create_from_records(
                            otu.id, isolate_name, list(isolate_records.values())
                        )
                        if isolate:
                            new_isolate_ids.append(isolate.id)
                    except ValueError as e:
                        log.error(
                            "Error creating isolate",
                            error=str(e),
                            isolate_name=isolate_name,
                        )

                if new_isolate_ids:
                    log.info("Added new isolates", count=len(new_isolate_ids))

        self._repo.write_otu_update_history_entry(otu.id)

        log.info("Comprehensive OTU update complete.")

        return self._repo.get_otu(otu.id)
