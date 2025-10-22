from collections.abc import Collection
from uuid import UUID

from pydantic import ValidationError
from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.builders.sequence import SequenceBuilder
from ref_builder.otu.utils import (
    DeleteRationale,
    assign_segment_id_to_record,
    create_segments_from_records,
)
from ref_builder.plan import (
    Plan,
    SegmentRule,
)
from ref_builder.repo import Repo

logger = get_logger("otu.modify")


def exclude_accessions_from_otu(
    repo: Repo,
    otu: OTUBuilder,
    accessions: Collection[str],
) -> None:
    """Exclude accessions from future addition to an OTU."""
    original_excluded_accessions = otu.excluded_accessions.copy()

    with repo.use_transaction():
        excluded_accessions = repo.exclude_accessions(
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


def allow_accessions_into_otu(
    repo: Repo,
    otu: OTUBuilder,
    accessions: Collection[str],
) -> None:
    """Allow accessions for future addition to an OTU.

    This reverses the effect of exclude_accessions_from_otu.
    """
    original_excluded_accessions = otu.excluded_accessions.copy()

    with repo.use_transaction():
        excluded_accessions = repo.allow_accessions(
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


def delete_isolate_from_otu(repo: Repo, otu: OTUBuilder, isolate_id: UUID) -> bool:
    """Remove an isolate from a specified OTU."""
    otu_logger = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

    isolate = otu.get_isolate(isolate_id)

    if not isolate:
        otu_logger.error("Isolate not found.", isolate_id=str(isolate_id))
        return False

    with repo.use_transaction():
        repo.delete_isolate(otu.id, isolate.id, rationale=DeleteRationale.USER)

    otu_logger.info(
        "Isolate removed.",
        name=isolate.name,
        removed_isolate_id=str(isolate_id),
        current_isolate_ids=[str(isolate_id_) for isolate_id_ in otu.isolate_ids],
    )

    repo.get_otu(otu.id)

    return True


def set_plan(
    repo: Repo,
    otu: OTUBuilder,
    plan: Plan,
) -> Plan | None:
    """Set an OTU's plan."""
    log = logger.bind(name=otu.name, taxid=otu.taxid, plan=plan.model_dump())

    try:
        with repo.use_transaction():
            repo.set_plan(otu.id, plan)
    except ValueError:
        log.exception()
        return None

    return repo.get_otu(otu.id).plan


def add_segments_to_plan(
    repo: Repo,
    otu: OTUBuilder,
    rule: SegmentRule,
    accessions: Collection[str],
    ignore_cache: bool = False,
) -> set[UUID]:
    """Add new segments to a multipartite plan."""
    log = logger.bind(name=otu.name, taxid=otu.taxid, rule=str(rule))

    sequence_length_tolerance = repo.settings.default_segment_length_tolerance

    if otu.plan.monopartite and otu.plan.segments[0].name is None:
        log.error(
            "Cannot add new segments unless all existing segments are named.",
            segment_id=str(otu.plan.segments[0].id),
            segment_name=otu.plan.segments[0].name,
        )

        return set()

    client = NCBIClient(ignore_cache)

    if not (records := client.fetch_genbank_records(accessions)):
        log.error(
            "Could not fetch records associated with requested accessions.",
            accessions=sorted(accessions),
        )
        return set()

    if len(records) < len(accessions):
        log.error(
            "Not all requested accessions could be fetched.",
            requested_accessions=sorted(accessions),
            accessible_record_accessions=sorted(
                record.accession_version for record in records
            ),
        )
        return set()

    new_segments = create_segments_from_records(
        records,
        rule,
        length_tolerance=sequence_length_tolerance,
    )
    if not new_segments:
        log.warning("No segments can be added.")
        return set()

    if len(new_segments) < len(accessions):
        log.error("Could not create all new segments.")
        return set()

    new_plan = otu.plan.model_copy()
    new_plan.segments.extend(new_segments)

    set_plan(repo, otu, new_plan)

    log.info("Added new segments", ids=[str(segment.id) for segment in new_segments])

    return {segment.id for segment in new_segments}
