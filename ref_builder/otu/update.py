import datetime
from abc import ABC, abstractmethod
from collections.abc import Collection, Iterable, Iterator
from pathlib import Path
from uuid import UUID

import arrow
from pydantic import RootModel
from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.promote import (
    promote_otu_accessions,
    promote_otu_accessions_from_records,
    upgrade_outdated_sequences_in_otu,
)
from ref_builder.otu.utils import (
    DeleteRationale,
    get_segments_max_length,
    get_segments_min_length,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.repo import Repo
from ref_builder.services.isolate import IsolateService

logger = get_logger("otu.update")

OTU_FEEDBACK_INTERVAL = 100
"""A default interval for batch OTU feedback."""

RECORD_FETCH_CHUNK_SIZE = 500
"""A default chunk size for NCBI EFetch calls."""

UPDATE_COOLDOWN_INTERVAL_IN_DAYS = 14


BatchFetchIndex = RootModel[dict[int, set[str]]]
"""Assists in reading and writing the fetched accessions by taxid index from file."""


class BaseBatchRecordGetter(ABC):
    """An abstract class with a .get_records() method."""

    @abstractmethod
    def get_records(self, taxid: int) -> list[NCBIGenbank]:
        """Return Genbank records corresponding to the given Taxonomy ID."""
        return NotImplemented


def comprehensive_update_otu(
    repo: Repo,
    otu: OTUBuilder,
    ignore_cache: bool = False,
) -> OTUBuilder:
    """Comprehensively update an OTU by promoting, upgrading, and adding new isolates.

    This function performs three operations in sequence:
    1. Promote GenBank accessions to RefSeq equivalents where available
    2. Upgrade outdated sequence versions (e.g., v1 → v2)
    3. Add new isolates from newly available accessions
    """
    log = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    log.info("Starting comprehensive OTU update.")

    # Step 1: Promote GenBank accessions to RefSeq
    promoted_accessions = promote_otu_accessions(repo, otu, ignore_cache)
    if promoted_accessions:
        log.info("Promoted sequences", count=len(promoted_accessions))
        otu = repo.get_otu(otu.id)

    # Step 2: Upgrade outdated sequence versions
    upgraded_sequence_ids = upgrade_outdated_sequences_in_otu(
        repo,
        otu,
        modification_date_start=None,
        ignore_cache=ignore_cache,
    )
    if upgraded_sequence_ids:
        log.info("Upgraded sequences", count=len(upgraded_sequence_ids))
        otu = repo.get_otu(otu.id)

    # Step 3: Add new isolates
    ncbi = NCBIClient(False)
    accessions = ncbi.filter_accessions(
        ncbi.fetch_accessions_by_taxid(
            otu.taxid,
            sequence_min_length=get_segments_min_length(otu.plan.segments),
            sequence_max_length=get_segments_max_length(otu.plan.segments),
        ),
    )

    fetch_set = {accession.key for accession in accessions} - otu.blocked_accessions

    if fetch_set:
        log.info("Adding new isolates from NCBI.")
        new_isolate_ids = update_otu_with_accessions(repo, otu, fetch_set, ignore_cache)

        if new_isolate_ids:
            log.info("Added new isolates", count=len(new_isolate_ids))

    repo.write_otu_update_history_entry(otu.id)

    log.info("Comprehensive OTU update complete.")

    return repo.get_otu(otu.id)


class PrecachedRecordStore(BaseBatchRecordGetter):
    """Retrieves records from an indexed dictionary of records and
    a batch fetch index set at initialization.
    """

    def __init__(
        self,
        batch_fetch_index: dict[int, set[str]],
        record_index: dict[str, NCBIGenbank],
    ) -> None:
        self.batch_fetch_index = batch_fetch_index
        self.record_index = record_index

    def get_records(self, taxid: int) -> list[NCBIGenbank]:
        accessions = self.batch_fetch_index.get(taxid, [])

        otu_records = [
            record
            for accession in accessions
            if (record := self.record_index.get(accession)) is not None
        ]

        return otu_records


class RecordFetcher(BaseBatchRecordGetter):
    """Retrieves records from NCBI Nucleotide based on a batch fetch index
    set at initialization.
    """

    def __init__(
        self, batch_fetch_index: dict[int, set[str]], ignore_cache: bool = False
    ) -> None:
        self.batch_fetch_index = batch_fetch_index
        self.ncbi = NCBIClient(ignore_cache)

    def get_records(self, taxid: int) -> list[NCBIGenbank]:
        accessions = self.batch_fetch_index.get(taxid, [])

        return self.ncbi.fetch_genbank_records(accessions)


def batch_update_repo(
    repo: Repo,
    start_date: datetime.date | None = None,
    skip_recently_updated: bool = True,
    ignore_cache: bool = False,
) -> set[UUID]:
    """Fetch new accessions for all OTUs in the repo and create isolates as possible."""
    operation_run_timestamp = arrow.utcnow().naive

    updated_otu_ids = set()

    repo_logger = logger.bind(
        path=str(repo.path),
    )
    if start_date is not None:
        repo_logger = repo_logger.bind(start_date.isoformat())

    repo_logger.info("Starting batch update...")

    if skip_recently_updated:
        otu_iterator = (
            otu
            for otu in repo.iter_otus()
            if _otu_is_cooled(
                repo,
                otu.id,
                timestamp_current=operation_run_timestamp,
            )
        )
    else:
        otu_iterator = repo.iter_otus()

    batch_fetch_index = batch_fetch_new_accessions(
        otu_iterator,
        modification_date_start=start_date,
        ignore_cache=ignore_cache,
    )

    fetch_index_cache_path = _cache_fetch_index(batch_fetch_index, repo.path / ".cache")

    repo_logger.info("Fetch index cached", fetch_index_path=fetch_index_cache_path)

    if not batch_fetch_index:
        logger.info("OTUs are up to date.")

        return updated_otu_ids

    logger.info(
        "Batch fetch index contains potential new accessions.",
        otu_count=len(batch_fetch_index),
    )

    fetch_set = {
        accession
        for otu_accessions in batch_fetch_index.values()
        for accession in otu_accessions
    }

    logger.info("Precaching records...", accession_count=len(fetch_set))

    record_index_by_accession = batch_fetch_new_records(
        fetch_set,
        chunk_size=RECORD_FETCH_CHUNK_SIZE,
        ignore_cache=ignore_cache,
    )

    if not record_index_by_accession:
        logger.info("No valid accessions found.")
        return updated_otu_ids

    record_getter = PrecachedRecordStore(batch_fetch_index, record_index_by_accession)

    for taxid, accessions in batch_fetch_index.items():
        if (otu_id := repo.get_otu_id_by_taxid(taxid)) is None:
            logger.debug("No corresponding OTU found in this repo", taxid=taxid)
            continue

        if skip_recently_updated and not _otu_is_cooled(
            repo,
            otu_id,
            timestamp_current=operation_run_timestamp,
        ):
            logger.info(
                "This OTU was updated recently. Skipping...",
                cooldown=UPDATE_COOLDOWN_INTERVAL_IN_DAYS,
                otu_id=str(otu_id),
                taxid=str(taxid),
            )
            continue

        otu_records = record_getter.get_records(taxid)

        if otu_records:
            isolate_ids = promote_and_update_otu_from_records(
                repo, repo.get_otu(otu_id), otu_records
            )

            if isolate_ids:
                updated_otu_ids.add(otu_id)

        repo.write_otu_update_history_entry(otu_id)

    repo_logger.info("Batch update complete.", new_isolate_count=len(updated_otu_ids))

    return updated_otu_ids


def batch_fetch_new_accessions(
    otus: Iterable[OTUBuilder],
    modification_date_start: datetime.date | None = None,
    modification_date_end: datetime.date | None = None,
    ignore_cache: bool = False,
) -> dict[int, set[str]]:
    """Check OTU iterator for new accessions and return results indexed by taxid."""
    ncbi = NCBIClient(ignore_cache)

    otu_counter = 0

    taxid_accession_index = {}

    for otu in otus:
        otu_logger = logger.bind(taxid=otu.taxid, name=otu.name)

        otu_counter += 1

        if otu_counter % OTU_FEEDBACK_INTERVAL == 0:
            otu_logger.info(
                "Fetching accession updates...",
                otu_counter=otu_counter,
            )

        raw_accessions = ncbi.fetch_accessions_by_taxid(
            otu.taxid,
            sequence_min_length=get_segments_min_length(otu.plan.segments),
            sequence_max_length=get_segments_max_length(otu.plan.segments),
            modification_date_start=modification_date_start,
            modification_date_end=modification_date_end,
        )

        accessions_to_fetch = {
            accession.key for accession in ncbi.filter_accessions(raw_accessions)
        } - otu.blocked_accessions

        if accessions_to_fetch:
            otu_logger.debug(
                "Potential accessions found.",
                accession_count=len(accessions_to_fetch),
                otu_counter=otu_counter,
            )

            taxid_accession_index[otu.taxid] = accessions_to_fetch

    return taxid_accession_index


def batch_fetch_new_records(
    accessions: Collection[str],
    chunk_size: int = RECORD_FETCH_CHUNK_SIZE,
    ignore_cache: bool = False,
) -> dict[str, NCBIGenbank]:
    """Download a batch of records and return in a dictionary indexed by accession."""
    fetch_logger = logger.bind(
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
    for fetch_list_chunk in iter_fetch_list(fetch_list, chunk_size):
        fetch_logger.info("Fetching records...", page_counter=page_counter)

        chunked_records = ncbi.fetch_genbank_records(fetch_list_chunk)

        indexed_records.update({record.accession: record for record in chunked_records})

        page_counter += 1

    if indexed_records:
        return indexed_records

    return {}


def update_isolate_from_records(
    repo: Repo,
    otu: OTUBuilder,
    isolate_id: UUID,
    records: list[NCBIGenbank],
) -> IsolateBuilder | None:
    """Take a list of GenBank records and replace the existing sequences,
    adding the previous sequence accessions to the excluded accessions list.
    """
    isolate = otu.get_isolate(isolate_id)

    if len(records) == 1:
        assigned = {otu.plan.segments[0].id: records[0]}
    else:
        try:
            assigned = assign_records_to_segments(records, otu.plan)
        except ValueError as e:
            logger.debug(e)
            return None

    for segment_id, record in assigned.items():
        _, accession = parse_refseq_comment(record.comment)

        if accession in isolate.accessions:
            old_sequence = isolate.get_sequence_by_accession(accession)

            # Use one transaction per sequence
            with repo.use_transaction():
                new_sequence = repo.create_sequence(
                    otu.id,
                    accession=record.accession_version,
                    definition=record.definition,
                    segment=segment_id,
                    sequence=record.sequence,
                )

                if new_sequence is None:
                    logger.error("Isolate update failed when creating new sequence.")
                    return None

                repo.replace_sequence(
                    otu.id,
                    isolate.id,
                    new_sequence.id,
                    replaced_sequence_id=old_sequence.id,
                    rationale=DeleteRationale.REFSEQ,
                )

                repo.exclude_accession(otu.id, old_sequence.accession.key)

    otu = repo.get_otu(otu.id)
    isolate = otu.get_isolate(isolate_id)

    logger.info(
        "Isolate updated",
        name=str(isolate.name),
        id=str(isolate.id),
        accessions=sorted([str(accession) for accession in isolate.accessions]),
    )

    return isolate


def promote_and_update_otu_from_records(
    repo: Repo,
    otu: OTUBuilder,
    records: list[NCBIGenbank],
) -> list[UUID]:
    """Promote new RefSeq accessions and add new isolates."""
    genbank_records, refseq_records = [], []

    for record in records:
        if record.refseq:
            refseq_records.append(record)

    if promote_otu_accessions_from_records(
        repo,
        otu=repo.get_otu(otu.id),
        records=refseq_records,
    ):
        otu = repo.get_otu(otu.id)

    new_isolate_ids = update_otu_with_records(
        repo,
        otu=otu,
        records=records,
    )

    repo.get_otu(otu.id)

    return new_isolate_ids


def update_otu_with_accessions(
    repo: Repo,
    otu: OTUBuilder,
    accessions: Collection[str],
    ignore_cache: bool = False,
) -> list[UUID]:
    """Take a list of accessions, filter for eligible accessions and
    add new sequences to the OTU.
    """
    log = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    ncbi = NCBIClient(ignore_cache)

    log.info(
        "Fetching records.",
        count=len(accessions),
        fetch_list=sorted(accessions),
    )

    records = ncbi.fetch_genbank_records(accessions)

    if records:
        return promote_and_update_otu_from_records(repo, otu, records)


def update_otu_with_records(
    repo: Repo,
    otu: OTUBuilder,
    records: list[NCBIGenbank],
) -> list[UUID]:
    """Take a list of downloaded NCBI Genbank records, filter for eligible records
    and add new sequences to the OTU.
    """
    new_isolate_ids = []
    isolate_service = IsolateService(repo, NCBIClient(False))

    for divided_records in (
        [record for record in records if record.refseq],
        [record for record in records if not record.refseq],
    ):
        otu = repo.get_otu(otu.id)

        for isolate_name, isolate_records in group_genbank_records_by_isolate(
            divided_records
        ).items():
            try:
                isolate = isolate_service.create_from_records(
                    otu.id, isolate_name, list(isolate_records.values())
                )
            except ValueError as e:
                logger.error(
                    "Error encountered while creating sequence.",
                    error_message=str(e),
                )

                isolate = None

            if isolate:
                new_isolate_ids.append(isolate.id)

    return new_isolate_ids


def iter_fetch_list(
    fetch_list: list[str], page_size=RECORD_FETCH_CHUNK_SIZE
) -> Iterator[list[str]]:
    """Divide a list of accessions and yield in pages."""
    page_size = max(1, page_size)
    page_count = len(fetch_list) // page_size + int(len(fetch_list) % page_size != 0)

    for iterator in range(page_count):
        yield fetch_list[iterator * page_size : (iterator + 1) * page_size]


def _generate_datestamp_filename() -> str:
    """Get the current UTC date and return as a a filename_safe string."""
    timestamp = arrow.utcnow().naive

    return f"{timestamp:%Y}_{timestamp:%m}_{timestamp:%d}"


def _cache_fetch_index(
    fetch_index: dict[int, set[str]],
    cache_path: Path,
) -> Path | None:
    """Write a batch fetch index to file."""
    validated_fetch_index = BatchFetchIndex.model_validate(fetch_index)

    fetch_index_path = (
        cache_path / f"fetch_index__{_generate_datestamp_filename()}.json"
    )

    with open(fetch_index_path, "w") as f:
        f.write(validated_fetch_index.model_dump_json())

    if fetch_index_path.exists():
        return fetch_index_path

    return None


def _otu_is_cooled(
    repo: Repo,
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
