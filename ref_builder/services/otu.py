"""Manage OTU data."""

from uuid import UUID, uuid4

import structlog

from ref_builder.errors import OTUExistsError
from ref_builder.events.isolate import CreateIsolateData
from ref_builder.models.accession import Accession
from ref_builder.models.isolate import IsolateName
from ref_builder.models.molecule import Molecule
from ref_builder.models.otu import OTU
from ref_builder.models.sequence import Sequence
from ref_builder.ncbi.models import NCBIGenbank, NCBIRank, NCBITaxonomy
from ref_builder.ncbi.utils import (
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.otu import assign_records_to_segments
from ref_builder.plan import (
    create_plan_from_records,
    get_segments_max_length,
    get_segments_min_length,
)
from ref_builder.promote import (
    assign_segment_id_to_record,
    promote_otu_from_records,
)
from ref_builder.services import Service

logger = structlog.get_logger("services.otu")

UUID_STRING_LENGTH = 36


class OTUService(Service):
    """Service for managing OTU operations."""

    def get_otu(self, identifier: str) -> OTU | None:
        """Return an OTU from the repo if identifier matches a single OTU.

        The identifier can either be a stringified UUIDv4 or a NCBI Taxonomy ID.

        :param identifier: an OTU identifier (UUID or taxid)
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

        if otu_id is None:
            return None

        return self._repo.get_otu(otu_id)

    def create(
        self,
        accessions: list[str],
    ) -> OTU | None:
        """Create a new OTU from a list of accessions.

        Uses the provided accessions to generate a plan and add a first isolate.
        Derives the taxonomy ID from the accessions.

        :param accessions: accessions to build the new OTU from
        :return: the created OTU or None if creation failed
        """
        log = logger.bind(accessions=accessions)

        if not accessions:
            log.error("OTU could not be created to spec based on given data.")
            return None

        records = self.ncbi.fetch_genbank_records(accessions)

        if len(records) != len(accessions):
            log.fatal("Could not retrieve all requested accessions.")
            return None

        if len({record.source.taxid for record in records}) > 1:
            log.fatal("Not all records are from the same organism.")
            return None

        taxid = records[0].source.taxid

        binned_records = group_genbank_records_by_isolate(records)

        if len(binned_records) > 1:
            log.fatal("More than one isolate found. Cannot create plan.")
            return None

        taxonomy = self.ncbi.fetch_taxonomy_record(taxid)

        if taxonomy is None:
            log.fatal("Could not retrieve data from NCBI Taxonomy", taxid=taxid)
            return None

        if taxonomy.rank != NCBIRank.SPECIES:
            taxonomy = self.ncbi.fetch_taxonomy_record(taxonomy.species.id)

        if taxonomy is None:
            log.fatal("Could not retrieve data from NCBI Taxonomy", taxid=taxid)
            return None

        if otu_id := self._repo.get_otu_id_by_taxid(taxonomy.id):
            raise OTUExistsError(taxonomy.id, otu_id)

        return self._write_otu(
            taxonomy,
            records,
            isolate_name=next(iter(binned_records.keys())) if binned_records else None,
        )

    def _write_otu(
        self,
        taxonomy: NCBITaxonomy,
        records: list[NCBIGenbank],
        isolate_name: IsolateName | None,
    ) -> OTU | None:
        """Create a new OTU from an NCBI Taxonomy record and a list of Nucleotide records.

        :param taxonomy: the NCBI taxonomy record
        :param records: the GenBank records
        :param isolate_name: the isolate name
        :return: the created OTU or None if creation failed
        """
        log = logger.bind(taxid=taxonomy.id)

        plan = create_plan_from_records(
            records,
            length_tolerance=self._repo.settings.default_segment_length_tolerance,
        )

        if plan is None:
            log.fatal("Could not create plan from records.")
            return None

        molecule = _get_molecule_from_records(records)

        lineage = self.ncbi.fetch_lineage(records[0].source.taxid)

        isolate_id = uuid4()

        if plan.monopartite:
            record = records[0]
            assigned = {plan.segments[0].id: record}
        else:
            assigned = assign_records_to_segments(records, plan)

        sequences = [
            Sequence(
                accession=Accession.from_string(record.accession_version),
                definition=record.definition,
                segment=segment_id,
                sequence=record.sequence,
            )
            for segment_id, record in assigned.items()
        ]

        promoted_accessions = set()
        for record in assigned.values():
            if record.refseq:
                _, old_accession = parse_refseq_comment(record.comment)
                promoted_accessions.add(old_accession)

        isolate_data = CreateIsolateData(
            id=isolate_id,
            name=isolate_name,
            sequences=sequences,
            taxid=records[0].source.taxid,
        )

        otu = self._repo.create_otu(
            isolate=isolate_data,
            lineage=lineage,
            molecule=molecule,
            plan=plan,
            promoted_accessions=promoted_accessions,
        )

        if otu is None:
            log.fatal("Could not create OTU.")
            return None

        return otu

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

    def _promote_accessions(self, otu: OTU) -> set[str]:
        """Fetch new accessions from NCBI Nucleotide and promote accessions
        with newly added RefSeq equivalents.

        :param otu: the OTU builder instance
        :return: set of promoted accession keys
        """
        log = logger.bind(otu_id=otu.id, taxid=otu.taxid)

        log.info("Checking for promotable sequences.")

        accessions = self.ncbi.fetch_accessions_by_taxid(
            otu.taxid,
            sequence_min_length=get_segments_min_length(otu.plan.segments),
            sequence_max_length=get_segments_max_length(otu.plan.segments),
            refseq_only=True,
        )

        fetch_set = {accession.key for accession in accessions} - otu.blocked_accessions

        if fetch_set:
            records = self.ncbi.fetch_genbank_records(fetch_set)

            log.debug(
                "New accessions found. Checking for promotable records.",
                fetch_list=sorted(fetch_set),
            )

            if promoted_accessions := promote_otu_from_records(
                self._repo, otu, records
            ):
                log.info(
                    "Sequences promoted.", new_accessions=sorted(promoted_accessions)
                )

                return promoted_accessions

        log.info("Records are already up to date.")

        return set()

    def _upgrade_outdated_sequences(
        self,
        otu: OTU,
    ) -> set[str]:
        """Fetch all extant accessions in the OTU and check if the record has been
        modified since last addition. Replace the sequence if an upgrade is found.
        """
        all_server_accessions = self.ncbi.fetch_accessions_by_taxid(
            otu.taxid,
            sequence_min_length=get_segments_min_length(otu.plan.segments),
            sequence_max_length=get_segments_max_length(otu.plan.segments),
        )

        server_upgraded_accessions = {
            accession for accession in all_server_accessions if accession.version > 2
        }

        replacement_index = {}

        for accession in server_upgraded_accessions:
            if (
                accession.key in otu.accessions
                and accession not in otu.versioned_accessions
            ):
                replacement_index[accession] = otu.get_sequence(accession.key)

        if not replacement_index:
            logger.info("All sequences are up to date.")
            return set()

        logger.info(
            "Upgradable sequences found. Fetching records...",
            upgradable_accessions=[str(accession) for accession in replacement_index],
        )

        records = self.ncbi.fetch_genbank_records(
            [str(accession) for accession in replacement_index]
        )

        upgraded_accessions = set()

        for record in records:
            outdated_sequence = otu.get_sequence(record.accession)
            versioned_accession = Accession.from_string(record.accession_version)

            segment_id = assign_segment_id_to_record(record, otu.plan)
            if segment_id is None:
                logger.error(
                    "Segment does not match plan",
                    accession=record.accession_version,
                )
                continue

            if versioned_accession not in otu.versioned_accessions:
                new_sequence = Sequence(
                    accession=Accession.from_string(record.accession_version),
                    definition=record.definition,
                    segment=segment_id,
                    sequence=record.sequence,
                )
            else:
                logger.info(
                    "Retrieving existing sequence",
                    accession=record.accession,
                )
                new_sequence = otu.get_sequence(record.accession)

            try:
                self._repo.update_sequence(
                    otu.id,
                    old_accession=outdated_sequence.accession,
                    new_sequence=new_sequence,
                )
                upgraded_accessions.add(record.accession)
            except ValueError as e:
                logger.error(
                    "Update failed",
                    error=str(e),
                    old_accession=outdated_sequence.accession,
                    new_accession=record.accession_version,
                )
                continue

            otu = self._repo.get_otu(otu.id)

        if upgraded_accessions:
            logger.info(
                "Replaced sequences",
                upgraded_accessions=sorted(upgraded_accessions),
            )

        return upgraded_accessions

    def update(self, otu_id: UUID) -> OTU | None:
        """Update an OTU by promoting, upgrading, and adding new isolates.

        This method performs three operations in sequence:
        1. Promote GenBank accessions to RefSeq equivalents where available
        2. Upgrade outdated sequence versions (e.g., v1 → v2)
        3. Add new isolates from newly available accessions

        :param otu_id: the OTU ID
        :return: the updated OTU or None if the OTU was not found
        """
        otu = self._repo.get_otu(otu_id)

        if otu is None:
            logger.error("OTU not found", otu_id=str(otu_id))
            return None

        log = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

        log.info("Starting comprehensive OTU update.")

        # Step 1: Promote GenBank accessions to RefSeq
        if self._promote_accessions(otu):
            otu = self._repo.get_otu(otu.id)

        # Step 2: Upgrade outdated sequence versions
        upgraded_accessions = self._upgrade_outdated_sequences(otu)

        if upgraded_accessions:
            log.info("Upgraded sequences", count=len(upgraded_accessions))
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
                refseq_records = [r for r in records if r.refseq]
                promoted_accessions = set()
                if refseq_records:
                    promoted_accessions = promote_otu_from_records(
                        self._repo, otu, refseq_records
                    )
                    if promoted_accessions:
                        otu = self._repo.get_otu(otu.id)

                records = [
                    r
                    for r in records
                    if r.accession.split(".")[0] not in promoted_accessions
                ]

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


def _get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
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
