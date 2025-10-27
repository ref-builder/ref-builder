import datetime
from uuid import UUID

from structlog import get_logger

from ref_builder.models.accession import Accession
from ref_builder.models.plan import Plan
from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.utils import (
    get_segments_max_length,
    get_segments_min_length,
    parse_refseq_comment,
)
from ref_builder.plan import extract_segment_name_from_record_with_plan
from ref_builder.repo import Repo

logger = get_logger("otu.promote")


def assign_segment_id_to_record(
    record: NCBIGenbank,
    plan: Plan,
) -> UUID | None:
    """Assign a segment ID to a record based on the plan."""
    segment_name = extract_segment_name_from_record_with_plan(record, plan)

    if segment_name is None and plan.monopartite:
        return plan.segments[0].id

    for segment in plan.segments:
        if segment_name == segment.name:
            return segment.id

    return None


def promote_otu_accessions(
    repo: Repo, otu: OTUBuilder, ignore_cache: bool = False
) -> set[str]:
    """Fetch new accessions from NCBI Nucleotide and promote accessions
    with newly added RefSeq equivalents.
    """
    ncbi = NCBIClient(ignore_cache)

    log = logger.bind(otu_id=otu.id, taxid=otu.taxid)

    log.info("Checking for promotable sequences.")

    accessions = ncbi.fetch_accessions_by_taxid(
        otu.taxid,
        sequence_min_length=get_segments_min_length(otu.plan.segments),
        sequence_max_length=get_segments_max_length(otu.plan.segments),
        refseq_only=True,
    )
    fetch_set = {accession.key for accession in accessions} - otu.blocked_accessions

    if fetch_set:
        records = ncbi.fetch_genbank_records(fetch_set)

        log.debug(
            "New accessions found. Checking for promotable records.",
            fetch_list=sorted(fetch_set),
        )

        if promoted_accessions := promote_otu_accessions_from_records(
            repo, otu, records
        ):
            log.info("Sequences promoted.", new_accessions=sorted(promoted_accessions))

            return promoted_accessions

    log.info("Records are already up to date.")

    return set()


def promote_otu_accessions_from_records(
    repo: Repo, otu: OTUBuilder, records: list[NCBIGenbank]
) -> set[str]:
    """Promote GenBank sequences to their RefSeq equivalents.

    Takes a list of records, identifies RefSeq records that replace existing
    GenBank sequences in the OTU, and promotes them using the PromoteSequence event.
    Returns the set of promoted accessions.
    """
    otu_logger = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

    refseq_records = [record for record in records if record.refseq]

    records_by_promotable_sequence_id = {}

    for record in refseq_records:
        try:
            _, predecessor_accession = parse_refseq_comment(record.comment)

        except ValueError as e:
            logger.debug(e, accession=record.accession_version, comment=record.comment)
            continue

        if predecessor_accession in otu.accessions:
            otu_logger.debug(
                "Replaceable accession found",
                predecessor_accession=predecessor_accession,
                promoted_accession=record.accession,
            )

            predecessor_sequence = otu.get_sequence_by_accession(predecessor_accession)

            records_by_promotable_sequence_id[predecessor_sequence.id] = record

    if not records_by_promotable_sequence_id:
        otu_logger.info("No promotable sequences found.")
        return set()

    promoted_accessions = set()

    for old_sequence_id, refseq_record in records_by_promotable_sequence_id.items():
        predecessor_sequence = otu.get_sequence_by_id(old_sequence_id)

        if predecessor_sequence is None:
            logger.warning(
                "Predecessor sequence not found", sequence_id=str(old_sequence_id)
            )
            continue

        segment_id = assign_segment_id_to_record(refseq_record, otu.plan)
        if segment_id is None:
            logger.error(
                "Segment does not match plan",
                accession=refseq_record.accession_version,
            )
            continue

        versioned_accession = Accession.from_string(refseq_record.accession_version)

        with repo.use_transaction() as active_transaction:
            if versioned_accession not in otu.versioned_accessions:
                new_sequence = repo.create_sequence(
                    otu.id,
                    accession=refseq_record.accession_version,
                    definition=refseq_record.definition,
                    segment=segment_id,
                    sequence=refseq_record.sequence,
                )

                if new_sequence is None:
                    logger.error(
                        "Failed to create new sequence",
                        accession=refseq_record.accession_version,
                    )
                    active_transaction.abort()
                    continue
            else:
                logger.info(
                    "Retrieving existing RefSeq sequence",
                    accession=refseq_record.accession,
                )
                new_sequence = otu.get_sequence_by_accession(refseq_record.accession)

            try:
                repo.promote_sequence(
                    otu.id,
                    old_sequence_id=predecessor_sequence.id,
                    new_sequence_id=new_sequence.id,
                )
                promoted_accessions.add(new_sequence.accession.key)
            except ValueError as e:
                logger.error(
                    "Promotion failed",
                    error=str(e),
                    old_accession=predecessor_sequence.accession.key,
                    new_accession=refseq_record.accession_version,
                )
                active_transaction.abort()
                continue

        otu = repo.get_otu(otu.id)

    if promoted_accessions:
        otu_logger.info(
            "Sequences promoted.",
            promoted_accessions=sorted(promoted_accessions),
        )

    return promoted_accessions


def upgrade_outdated_sequences_in_otu(
    repo: Repo,
    otu: OTUBuilder,
    modification_date_start: datetime.datetime | None = None,
    ignore_cache: bool = False,
) -> set[UUID]:
    """Fetch all extant accessions in the OTU and check if the record has been
    modified since last addition. Replace the sequence if an upgrade is found.
    """
    ncbi = NCBIClient(ignore_cache)

    all_server_accessions = ncbi.fetch_accessions_by_taxid(
        otu.taxid,
        modification_date_start=modification_date_start,
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
            replacement_index[accession] = otu.get_sequence_by_accession(
                accession.key
            ).id

    if not replacement_index:
        logger.info("All sequences are up to date.")
        return set()

    logger.info(
        "Upgradable sequences found. Fetching records...",
        upgradable_accessions=[str(accession) for accession in replacement_index],
    )

    records = ncbi.fetch_genbank_records(
        [str(accession) for accession in replacement_index]
    )

    replacement_sequence_ids = set()

    for record in records:
        outmoded_sequence = otu.get_sequence_by_accession(record.accession)
        versioned_accession = Accession.from_string(record.accession_version)

        logger.info(
            "Replacing sequence...",
            sequence_id=str(outmoded_sequence.id),
            outdated_accession=str(outmoded_sequence.accession),
            new_accession=str(versioned_accession),
        )

        segment_id = assign_segment_id_to_record(record, otu.plan)
        if segment_id is None:
            logger.error(
                "Segment does not match plan",
                accession=record.accession_version,
            )
            continue

        with repo.use_transaction() as active_transaction:
            if versioned_accession not in otu.versioned_accessions:
                new_sequence = repo.create_sequence(
                    otu.id,
                    accession=record.accession_version,
                    definition=record.definition,
                    segment=segment_id,
                    sequence=record.sequence,
                )

                if new_sequence is None:
                    logger.error(
                        "Failed to create new sequence",
                        accession=record.accession_version,
                    )
                    active_transaction.abort()
                    continue
            else:
                logger.info(
                    "Retrieving existing sequence",
                    accession=record.accession,
                )
                new_sequence = otu.get_sequence_by_accession(record.accession)

            try:
                repo.update_sequence(
                    otu.id,
                    old_sequence_id=outmoded_sequence.id,
                    new_sequence_id=new_sequence.id,
                )
                replacement_sequence_ids.add(new_sequence.id)
            except ValueError as e:
                logger.error(
                    "Update failed",
                    error=str(e),
                    old_accession=outmoded_sequence.accession,
                    new_accession=record.accession_version,
                )
                active_transaction.abort()
                continue

        otu = repo.get_otu(otu.id)

    if replacement_sequence_ids:
        logger.info(
            "Replaced sequences",
            new_sequence_ids=[
                str(sequence_id) for sequence_id in replacement_sequence_ids
            ],
        )

    return replacement_sequence_ids
