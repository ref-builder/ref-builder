"""Manage OTU data."""

import structlog

from ref_builder.ncbi.models import NCBIGenbank, NCBIRank, NCBITaxonomy
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.isolate import create_sequence_from_record
from ref_builder.otu.utils import (
    assign_records_to_segments,
    create_plan_from_records,
    get_molecule_from_records,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.services import Service
from ref_builder.utils import IsolateName

logger = structlog.get_logger("services.otu")


class OTUService(Service):
    """Service for managing OTU operations."""

    def create(
        self,
        accessions: list[str],
        acronym: str = "",
    ) -> OTUBuilder | None:
        """Create a new OTU from a list of accessions.

        Uses the provided accessions to generate a plan and add a first isolate.
        Derives the taxonomy ID from the accessions.

        :param accessions: accessions to build the new OTU from
        :param acronym: an alternative name to use during searches
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
            otu_logger.fatal(f"Could not retrieve {taxid} from NCBI Taxonomy")
            return None

        if taxonomy.rank != NCBIRank.SPECIES:
            taxonomy = self.ncbi.fetch_taxonomy_record(taxonomy.species.id)

        if self.repo.get_otu_id_by_taxid(taxonomy.id):
            otu_logger.error(
                f"Taxonomy ID {taxonomy.id} has already been added to this reference."
            )
            return None

        if not acronym and taxonomy.other_names.acronym:
            acronym = taxonomy.other_names.acronym[0]

        with self.repo.use_transaction():
            return self._write_otu(
                taxonomy,
                records,
                acronym,
                isolate_name=next(iter(binned_records.keys())) if binned_records else None,
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
            length_tolerance=self.repo.settings.default_segment_length_tolerance,
        )

        if plan is None:
            otu_logger.fatal("Could not create plan from records.")
            return None

        molecule = get_molecule_from_records(records)

        otu = self.repo.create_otu(
            acronym=acronym,
            legacy_id=None,
            molecule=molecule,
            name=taxonomy.name,
            plan=plan,
            taxid=taxonomy.id,
        )

        isolate = self.repo.create_isolate(
            otu_id=otu.id,
            legacy_id=None,
            name=isolate_name,
        )

        otu.add_isolate(isolate)

        if otu.plan.monopartite:
            record = records[0]

            sequence = create_sequence_from_record(
                self.repo, otu, record, plan.segments[0].id
            )

            self.repo.link_sequence(otu.id, isolate.id, sequence.id)

            if record.refseq:
                _, old_accession = parse_refseq_comment(record.comment)

                self.repo.exclude_accession(
                    otu.id,
                    old_accession,
                )

        else:
            for segment_id, record in assign_records_to_segments(records, plan).items():
                sequence = create_sequence_from_record(
                    self.repo, otu, record, segment_id
                )

                self.repo.link_sequence(otu.id, isolate.id, sequence.id)

                if record.refseq:
                    _, old_accession = parse_refseq_comment(record.comment)
                    self.repo.exclude_accession(
                        otu.id,
                        old_accession,
                    )

        return self.repo.get_otu(otu.id)
