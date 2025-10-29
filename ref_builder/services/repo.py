import datetime
from collections.abc import Collection, Iterable, Iterator
from uuid import UUID

import arrow
import structlog

from ref_builder.models.otu import OTU
from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.ncbi.utils import group_genbank_records_by_isolate
from ref_builder.plan import get_segments_max_length, get_segments_min_length
from ref_builder.promote import promote_otu_from_records
from ref_builder.services import Service

logger = structlog.get_logger("services.repo")

OTU_FEEDBACK_INTERVAL = 100
"""A default interval for batch OTU feedback."""

RECORD_FETCH_CHUNK_SIZE = 500
"""A default chunk size for NCBI EFetch calls."""

UPDATE_COOLDOWN_INTERVAL_IN_DAYS = 14
"""A default chunk size for NCBI EFetch calls."""


class RepoService(Service):
    """A service for repo-level operations.

    TODO: Avoid refetching recently updated accessions.
    """

    def update(self) -> set[UUID]:
        """Update all OTUs in the repository.

        Fetch new accessions for all OTUs in the repo and create isolates as possible.

        :param start_date: exclude records edited before this date
        :param skip_recently_updated: skip OTUs updated within cooldown period
        :param ignore_cache: whether to ignore the NCBI cache
        :return: set of updated OTU IDs
        """
        log = logger.bind(
            path=str(self._repo.path),
        )

        log.info("Starting batch update")

        operation_run_timestamp = arrow.utcnow().naive

        otu_iterator = (
            otu
            for otu in self._repo.iter_otus()
            if _otu_is_cooled(
                self._repo,
                otu.id,
                timestamp_current=operation_run_timestamp,
            )
        )

        batch_fetch_index = _fetch_new_accessions(otu_iterator)

        if not batch_fetch_index:
            logger.info("OTUs are up to date.")
            return set()

        logger.info(
            "Batch fetch index contains potential new accessions.",
            otu_count=len(batch_fetch_index),
        )

        fetch_set = {
            accession
            for otu_accessions in batch_fetch_index.values()
            for accession in otu_accessions
        }

        record_index_by_accession = _fetch_new_records(fetch_set)

        if not record_index_by_accession:
            logger.info("No valid accessions found.")
            return set()

        updated_otu_ids = set()

        for taxid, accessions in batch_fetch_index.items():
            if (otu_id := self._repo.get_otu_id_by_taxid(taxid)) is None:
                logger.debug("No corresponding OTU found in this repo", taxid=taxid)
                continue

            otu_records = [
                record
                for accession in accessions
                if (record := record_index_by_accession.get(accession)) is not None
            ]

            if otu_records:
                otu = self._repo.get_otu(otu_id)

                # Promote RefSeq accessions first
                refseq_records = [r for r in otu_records if r.refseq]
                if refseq_records and promote_otu_from_records(
                    self._repo, otu, refseq_records
                ):
                    otu = self._repo.get_otu(otu_id)

                # Create isolates from records
                for isolate_name, isolate_records in group_genbank_records_by_isolate(
                    otu_records
                ).items():
                    try:
                        isolate = self._services.isolate.create_from_records(
                            otu.id, isolate_name, list(isolate_records.values())
                        )
                        if isolate:
                            updated_otu_ids.add(otu_id)
                    except ValueError as e:
                        logger.error(
                            "Error creating isolate",
                            error=str(e),
                            isolate_name=isolate_name,
                        )

            self._repo.write_otu_update_history_entry(otu_id)

        log.info("Batch update complete.", new_isolate_count=len(updated_otu_ids))

        return updated_otu_ids


def _otu_is_cooled(
    repo,
    otu_id: UUID,
    timestamp_current: datetime.datetime | None,
    cooldown: int = UPDATE_COOLDOWN_INTERVAL_IN_DAYS,
) -> bool:
    """Return the update cooldown status of an OTU.

    Return True if:
        A) The OTU has been recorded as updated within ``cooldown`` days,
        B) The OTU was created within ``cooldown`` days.
        C) The OTU was created outside and last updated outside ``cooldown`` days.
    """
    cooldown_delta = datetime.timedelta(days=cooldown)

    if timestamp_current is None:
        timestamp_current = arrow.utcnow().naive

    if (timestamp_last_updated := repo.get_otu_last_updated(otu_id)) is not None:
        return timestamp_current - timestamp_last_updated > cooldown_delta

    timestamp_created = repo.get_otu_first_created(otu_id)

    if (timestamp_current - timestamp_created) <= cooldown_delta:
        return True

    timestamp_latest = repo.get_otu_last_modified(otu_id)

    if (timestamp_current - timestamp_latest) > cooldown_delta:
        return True

    logger.debug(
        "Last modified timestamp occured within cooldown",
        timestamp_current=timestamp_current.isoformat(),
        timestamp_otu_created=timestamp_created.isoformat(),
        timestamp_otu_latest=timestamp_latest.isoformat(),
        otu_last_modified_delta=timestamp_current - timestamp_latest,
    )

    return False


def _fetch_new_accessions(
    otus: Iterable[OTU],
    ignore_cache: bool = False,
) -> dict[int, set[str]]:
    """Check OTU iterator for new accessions and return results indexed by taxid."""
    ncbi = NCBIClient(ignore_cache)

    otu_counter = 0

    taxid_accession_index = {}

    for otu in otus:
        log = logger.bind(taxid=otu.taxid, name=otu.name)

        otu_counter += 1

        if otu_counter % OTU_FEEDBACK_INTERVAL == 0:
            log.info(
                "Fetching accession updates...",
                otu_counter=otu_counter,
            )

        accessions = ncbi.fetch_accessions_by_taxid(
            otu.taxid,
            sequence_min_length=get_segments_min_length(otu.plan.segments),
            sequence_max_length=get_segments_max_length(otu.plan.segments),
        )

        accessions_to_fetch = {
            accession.key for accession in accessions
        } - otu.blocked_accessions

        if accessions_to_fetch:
            log.debug(
                "Potential accessions found.",
                accession_count=len(accessions_to_fetch),
                otu_counter=otu_counter,
            )

            taxid_accession_index[otu.taxid] = accessions_to_fetch

    return taxid_accession_index


def _iter_fetch_list(
    fetch_list: list[str], page_size: int = RECORD_FETCH_CHUNK_SIZE
) -> Iterator[list[str]]:
    """Divide a list of accessions and yield in pages."""
    page_size = max(1, page_size)
    page_count = len(fetch_list) // page_size + int(len(fetch_list) % page_size != 0)

    for iterator in range(page_count):
        yield fetch_list[iterator * page_size : (iterator + 1) * page_size]


def _fetch_new_records(
    accessions: Collection[str],
    chunk_size: int = RECORD_FETCH_CHUNK_SIZE,
    ignore_cache: bool = False,
) -> dict[str, NCBIGenbank]:
    """Download a batch of records and return in a dictionary indexed by accession."""
    log = logger.bind(
        accession_count=len(accessions),
        chunk_size=chunk_size,
        ignore_cache=ignore_cache,
    )

    if not accessions:
        return {}

    ncbi = NCBIClient(ignore_cache)

    fetch_list = list(accessions)

    page_counter = 0

    indexed_records = {}
    for fetch_list_chunk in _iter_fetch_list(fetch_list, chunk_size):
        log.info("Fetching records...", page_counter=page_counter)

        chunked_records = ncbi.fetch_genbank_records(fetch_list_chunk)

        indexed_records.update({record.accession: record for record in chunked_records})

        page_counter += 1

    if indexed_records:
        return indexed_records

    return {}
