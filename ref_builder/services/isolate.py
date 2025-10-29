"""Manage isolate data."""

from uuid import UUID, uuid4

import structlog

from ref_builder.errors import PlanValidationError
from ref_builder.models.accession import Accession
from ref_builder.models.isolate import Isolate, IsolateName
from ref_builder.models.otu import OTU
from ref_builder.models.sequence import Sequence
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.ncbi.utils import (
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.otu import assign_records_to_segments
from ref_builder.promote import promote_otu_from_records
from ref_builder.services import Service
from ref_builder.utils import filter_accessions

logger = structlog.get_logger("services.isolate")


class IsolateService(Service):
    """Service for managing isolate operations."""

    def create(
        self,
        accessions: list[str],
    ) -> Isolate | None:
        """Create a new isolate from a list of accessions.

        The isolate name is extracted from GenBank record metadata. If the records
        contain multiple isolate names or no isolate name, creation will fail.

        The OTU is automatically determined by fetching the GenBank records and
        extracting the taxid, then finding an OTU that contains that taxid in its
        lineage.

        :param accessions: accessions to build the new isolate from
        :return: the created isolate or None if creation failed
        """
        log = logger.bind(accessions=accessions)

        # Fetch records to determine taxid
        fetched_records = self.ncbi.fetch_genbank_records(accessions)

        if not fetched_records:
            log.error("Failed to fetch records from NCBI")
            return None

        # Extract taxid from records
        taxids = {record.source.taxid for record in fetched_records}
        if len(taxids) > 1:
            log.error(
                "Not all records have the same taxid.",
                taxids=sorted(taxids),
            )
            return None

        taxid = fetched_records[0].source.taxid
        log = log.bind(taxid=taxid)

        # Find OTU by taxid
        otu = self._repo.get_otu_by_taxid(taxid)

        if otu is None:
            log.error("No OTU found for taxid")
            return None

        log = log.bind(otu_id=str(otu.id), otu_name=otu.name)

        # Filter out blocked accessions
        eligible_accessions = filter_accessions(
            [r.accession for r in fetched_records],
            otu.blocked_accessions,
        )

        if not eligible_accessions:
            log.error("All fetched accessions are blocked for this OTU")
            return None

        records = [r for r in fetched_records if r.accession in eligible_accessions]

        # Validate that we don't mix RefSeq and non-RefSeq sequences
        if any(record.refseq for record in records) and not all(
            record.refseq for record in records
        ):
            log.error(
                "Cannot mix RefSeq and non-RefSeq sequences in multipartite isolate."
            )
            return None

        binned_records = group_genbank_records_by_isolate(records)

        if len(binned_records) != 1:
            log.error("More than one isolate name found in requested accessions.")
            return None

        isolate_name, _ = next(iter(binned_records.items()))

        # Try to promote sequences if all records are RefSeq
        if all(record.refseq for record in records):
            log.info("Checking for promotable sequences")

            promoted_accessions = promote_otu_from_records(self._repo, otu, records)

            if promoted_accessions:
                otu = self._repo.get_otu(otu.id)
                # Find the isolate by one of the promoted accessions
                promoted_accession = next(iter(promoted_accessions))
                isolate = otu.get_isolate_by_accession(promoted_accession)

                if isolate:
                    log.info(
                        "Sequences promoted",
                        promoted_accessions=sorted(promoted_accessions),
                    )
                    return isolate

            log.warning("No promotable sequences found")

        # Create the isolate
        isolate = self._write_isolate(otu, isolate_name, records)

        if isolate:
            return isolate

        return None

    def create_from_records(
        self,
        otu_id: UUID,
        isolate_name: IsolateName | None,
        records: list[NCBIGenbank],
    ) -> Isolate | None:
        """Create a new isolate from pre-fetched GenBank records.

        Use this method when records are already fetched (e.g., in batch operations).
        For creating isolates from accessions, use create() instead.

        :param otu_id: the OTU ID to add the isolate to
        :param isolate_name: the isolate name (or None for unnamed)
        :param records: the GenBank records
        :return: the created isolate or None if creation failed
        """
        otu = self._repo.get_otu(otu_id)

        if otu is None:
            logger.error("OTU not found", otu_id=str(otu_id))
            return None

        isolate = self._write_isolate(otu, isolate_name, records)

        if isolate:
            return isolate

        return None

    def delete(self, isolate_id: UUID, message: str) -> bool:
        """Delete an isolate.

        :param isolate_id: the isolate ID to delete
        :return: True if deletion succeeded, False otherwise
        """
        isolate = self._repo.get_isolate(isolate_id)

        if isolate is None:
            logger.error("Isolate not found", isolate_id=str(isolate_id))
            return False

        otu_id = self._repo.get_otu_id_by_isolate_id(isolate_id)

        if otu_id is None:
            logger.error("OTU not found for isolate", isolate_id=str(isolate_id))
            return False

        otu = self._repo.get_otu(otu_id)

        log = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

        self._repo.delete_isolate(otu.id, isolate.id, message)

        log.info(
            "Isolate removed.",
            name=isolate.name,
            removed_isolate_id=str(isolate_id),
            current_isolate_ids=[str(isolate_id_) for isolate_id_ in otu.isolate_ids],
        )

        self._repo.get_otu(otu.id)

        return True

    def _write_isolate(
        self,
        otu: OTU,
        isolate_name: IsolateName | None,
        records: list[NCBIGenbank],
    ) -> Isolate | None:
        """Create a new isolate and its sequences in the repository.

        :param otu: the OTU to add the isolate to
        :param isolate_name: the isolate name (or None for unnamed)
        :param records: the GenBank records
        :return: the created isolate or None if creation failed
        """
        log = logger.bind(
            name=str("Unnamed" if isolate_name is None else str(isolate_name)),
            otu_name=otu.name,
            otu_id=str(otu.id),
            taxid=otu.taxid,
            accessions=[record.accession for record in records],
        )

        # Validate all records have the same taxid
        taxids = {record.source.taxid for record in records}
        if len(taxids) > 1:
            log.error(
                "Not all records have the same taxid.",
                taxids=sorted(taxids),
            )
            return None

        taxid = records[0].source.taxid

        # Validate taxid exists in OTU lineage
        lineage_taxids = {taxon.id for taxon in otu.lineage.taxa}

        if taxid not in lineage_taxids:
            log.error(
                "Taxid not found in OTU lineage.",
                taxid=taxid,
                lineage_taxids=sorted(lineage_taxids),
            )
            return None

        try:
            assigned = assign_records_to_segments(records, otu.plan)
        except PlanValidationError as e:
            log.warning(
                str(e),
                segment_names=[str(segment.name) for segment in otu.plan.segments],
            )

            return None

        isolate_id = uuid4()

        sequences = [
            Sequence(
                accession=Accession.from_string(record.accession_version),
                definition=record.definition,
                segment=segment_id,
                sequence=record.sequence,
            )
            for segment_id, record in assigned.items()
        ]

        isolate = self._repo.create_isolate(
            otu_id=otu.id,
            isolate_id=isolate_id,
            name=isolate_name,
            taxid=taxid,
            sequences=sequences,
        )

        for record in assigned.values():
            if record.refseq:
                _, old_accession = parse_refseq_comment(record.comment)

                self._repo.exclude_accessions(
                    otu.id,
                    [old_accession],
                )

        log.info("Isolate created", id=str(isolate_id))

        return isolate
