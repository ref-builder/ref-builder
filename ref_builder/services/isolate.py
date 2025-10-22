"""Manage isolate data."""

from uuid import UUID

import structlog

from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.isolate import create_sequence_from_record
from ref_builder.otu.promote import promote_otu_accessions_from_records
from ref_builder.otu.utils import (
    DeleteRationale,
    assign_records_to_segments,
    fetch_records_from_accessions,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.plan import PlanConformationError
from ref_builder.services import Service
from ref_builder.utils import IsolateName

logger = structlog.get_logger("services.isolate")


class IsolateService(Service):
    """Service for managing isolate operations."""

    def create(
        self,
        otu_id: UUID,
        accessions: list[str],
        ignore_cache: bool = False,
    ) -> IsolateBuilder | None:
        """Create a new isolate from a list of accessions.

        The isolate name is extracted from GenBank record metadata. If the records
        contain multiple isolate names or no isolate name, creation will fail.

        :param otu_id: the OTU ID to add the isolate to
        :param accessions: accessions to build the new isolate from
        :param ignore_cache: whether to ignore the NCBI cache
        :return: the created isolate or None if creation failed
        """
        otu = self._repo.get_otu(otu_id)

        if otu is None:
            logger.error("OTU not found", otu_id=str(otu_id))
            return None

        otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

        records = fetch_records_from_accessions(
            accessions, otu.blocked_accessions, self.ncbi
        )

        if not records:
            return None

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
        )

        for segment_id, record in assigned.items():
            if (sequence := otu.get_sequence_by_accession(record.accession)) is None:
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

        log.info("Isolate created", id=str(isolate.id))

        return self._repo.get_isolate(isolate.id)
