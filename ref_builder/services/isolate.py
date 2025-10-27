"""Manage isolate data."""

from uuid import UUID

import structlog

from ref_builder.errors import PlanConformationError
from ref_builder.models.isolate import IsolateName
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.isolate import create_sequence_from_record
from ref_builder.otu.promote import promote_otu_accessions_from_records
from ref_builder.otu.utils import (
    DeleteRationale,
    assign_records_to_segments,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.services import Service
from ref_builder.utils import filter_accessions

logger = structlog.get_logger("services.isolate")


class IsolateService(Service):
    """Service for managing isolate operations."""

    def create(
        self,
        accessions: list[str],
    ) -> IsolateBuilder | None:
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

        otu_logger = log.bind(otu_id=str(otu.id), otu_name=otu.name)

        # Filter out blocked accessions
        eligible_accessions = filter_accessions(
            [r.accession for r in fetched_records],
            otu.blocked_accessions,
        )

        if not eligible_accessions:
            otu_logger.error("All fetched accessions are blocked for this OTU")
            return None

        records = [r for r in fetched_records if r.accession in eligible_accessions]

        # Validate that we don't mix RefSeq and non-RefSeq sequences
        if any(record.refseq for record in records) and not all(
            record.refseq for record in records
        ):
            otu_logger.error(
                "Cannot mix RefSeq and non-RefSeq sequences in multipartite isolate."
            )
            return None

        binned_records = group_genbank_records_by_isolate(records)

        if len(binned_records) != 1:
            otu_logger.error(
                "More than one isolate name found in requested accessions."
            )
            return None

        isolate_name, _ = next(iter(binned_records.items()))

        # Check for existing isolate with same name
        if (isolate_id := otu.get_isolate_id_by_name(isolate_name)) is not None:
            otu_logger.warning(
                "Isolate name already exists in this OTU.", name=isolate_name
            )

            # Try to promote sequences if all records are RefSeq
            if all(record.refseq for record in records):
                otu_logger.info("Attempting to promote sequences to RefSeq")

                promoted_accessions = promote_otu_accessions_from_records(
                    self._repo, otu, records
                )

                if promoted_accessions:
                    otu = self._repo.get_otu(otu.id)
                    isolate = otu.get_isolate(isolate_id)
                    otu_logger.info(
                        "Sequences promoted",
                        promoted_accessions=sorted(promoted_accessions),
                    )
                    return isolate

                otu_logger.warning("No promotable sequences found")

            return None

        # Create the isolate
        with self._repo.use_transaction() as transaction:
            isolate = self._write_isolate(otu, isolate_name, records)

            if isolate:
                return isolate

            transaction.abort()

        return None

    def create_from_records(
        self,
        otu_id: UUID,
        isolate_name: IsolateName | None,
        records: list[NCBIGenbank],
    ) -> IsolateBuilder | None:
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

        with self._repo.use_transaction() as transaction:
            isolate = self._write_isolate(otu, isolate_name, records)

            if isolate:
                return isolate

            transaction.abort()

        return None

    def delete(self, otu_id: UUID, isolate_id: UUID) -> bool:
        """Delete an isolate from an OTU.

        :param otu_id: the OTU ID
        :param isolate_id: the isolate ID to delete
        :return: True if deletion succeeded, False otherwise
        """
        otu = self._repo.get_otu(otu_id)

        if otu is None:
            logger.error("OTU not found", otu_id=str(otu_id))
            return False

        otu_logger = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

        isolate = otu.get_isolate(isolate_id)

        if not isolate:
            otu_logger.error("Isolate not found.", isolate_id=str(isolate_id))
            return False

        with self._repo.use_transaction():
            self._repo.delete_isolate(
                otu.id, isolate.id, rationale=DeleteRationale.USER
            )

        otu_logger.info(
            "Isolate removed.",
            name=isolate.name,
            removed_isolate_id=str(isolate_id),
            current_isolate_ids=[str(isolate_id_) for isolate_id_ in otu.isolate_ids],
        )

        self._repo.get_otu(otu.id)

        return True

    def _write_isolate(
        self,
        otu: OTUBuilder,
        isolate_name: IsolateName | None,
        records: list[NCBIGenbank],
    ) -> IsolateBuilder | None:
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
        except PlanConformationError as e:
            log.warning(
                str(e),
                segment_names=[str(segment.name) for segment in otu.plan.segments],
            )

            return None

        isolate = self._repo.create_isolate(
            otu.id,
            name=isolate_name,
            taxid=taxid,
        )

        for segment_id, record in assigned.items():
            if (sequence := otu.get_sequence_by_accession(record.accession)) is None:
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

        log.info("Isolate created", id=str(isolate.id))

        return self._repo.get_isolate(isolate.id)


def _fetch_records(
    accessions: list | set,
    blocked_accessions: set,
    ncbi_client,
) -> list[NCBIGenbank]:
    """Fetch GenBank records from a list of accessions.

    Don't fetch accessions in ``blocked_accessions``.

    :param accessions: A list of accessions to fetch.
    :param blocked_accessions: A set of accessions to ignore.
    :param ncbi_client: NCBI client to use for fetching records.
    """
    log = logger.bind(
        requested=sorted(accessions),
        blocked=sorted(blocked_accessions),
    )

    eligible = set(accessions) - blocked_accessions

    if not eligible:
        log.error("None of the requested accessions were eligible for inclusion.")
        return []

    log.debug(
        "Fetching accessions.",
        accessions=sorted(eligible),
        count=len(eligible),
    )

    return ncbi_client.fetch_genbank_records(eligible)
