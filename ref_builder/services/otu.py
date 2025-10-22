"""Manage OTU data."""

from uuid import UUID

import structlog

from ref_builder.errors import InvalidInputError, PartialIDConflictError
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
from ref_builder.otu.utils import (
    assign_records_to_segments,
    create_segments_from_records,
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

        The identifier can either be a stringified UUIDv4, a truncated portion
        of a UUID, a NCBI Taxonomy ID or an acronym associated with the OTU.

        :param identifier: a non-UUID identifier.
            Can be an integer Taxonomy ID, acronym or truncated partial UUID.
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

        elif (otu_id := self._repo.get_otu_id_by_acronym(identifier)) is None:
            try:
                otu_id = self._repo.get_otu_id_by_partial(identifier)

            except (PartialIDConflictError, InvalidInputError):
                return None

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

        acronym = ""
        if taxonomy.other_names.acronym:
            acronym = taxonomy.other_names.acronym[0]

        with self._repo.use_transaction():
            return self._write_otu(
                taxonomy,
                records,
                acronym,
                isolate_name=next(iter(binned_records.keys()))
                if binned_records
                else None,
            )

    def _write_otu(
        self,
        taxonomy: NCBITaxonomy,
        records: list[NCBIGenbank],
        acronym: str,
        isolate_name: IsolateName | None,
    ) -> OTUBuilder | None:
        """Create a new OTU from an NCBI Taxonomy record and a list of Nucleotide records.

        :param taxonomy: the NCBI taxonomy record
        :param records: the GenBank records
        :param acronym: the OTU acronym
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

        otu = self._repo.create_otu(
            acronym=acronym,
            molecule=molecule,
            name=taxonomy.name,
            plan=plan,
            taxid=taxonomy.id,
        )

        isolate = self._repo.create_isolate(
            otu_id=otu.id,
            name=isolate_name,
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

                self._repo.exclude_accession(
                    otu.id,
                    old_accession,
                )

        else:
            for segment_id, record in assign_records_to_segments(records, plan).items():
                sequence = create_sequence_from_record(
                    self._repo, otu, record, segment_id
                )

                self._repo.link_sequence(otu.id, isolate.id, sequence.id)

                if record.refseq:
                    _, old_accession = parse_refseq_comment(record.comment)
                    self._repo.exclude_accession(
                        otu.id,
                        old_accession,
                    )

        return self._repo.get_otu(otu.id)
