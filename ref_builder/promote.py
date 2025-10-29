from uuid import UUID

from structlog import get_logger

from ref_builder.models.otu import OTU
from ref_builder.models.plan import Plan
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.ncbi.utils import parse_refseq_comment
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


def promote_otu_from_records(
    repo: Repo, otu: OTU, records: list[NCBIGenbank]
) -> set[str]:
    """Promote GenBank sequences to their RefSeq equivalents.

    Takes a list of records, identifies RefSeq records that replace existing
    GenBank sequences in the OTU, and promotes them using the PromoteSequence event.
    Returns the set of promoted accessions.
    """
    log = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

    refseq_records = [record for record in records if record.refseq]

    # Group RefSeq records by the isolate they replace
    isolate_promotion_map = {}

    for record in refseq_records:
        try:
            _, predecessor_accession = parse_refseq_comment(record.comment)

        except ValueError as e:
            logger.debug(e, accession=record.accession_version, comment=record.comment)
            continue

        if predecessor_accession in otu.accessions:
            log.debug(
                "Replaceable accession found",
                predecessor_accession=predecessor_accession,
                promoted_accession=record.accession,
            )

            isolate = otu.get_isolate_by_accession(predecessor_accession)

            if isolate is None:
                logger.warning(
                    "Isolate not found for accession", accession=predecessor_accession
                )
                continue

            if isolate.id not in isolate_promotion_map:
                isolate_promotion_map[isolate.id] = {}

            isolate_promotion_map[isolate.id][predecessor_accession] = record

    if not isolate_promotion_map:
        log.info("No promotable sequences found.")
        return set()

    promoted_accessions = set()

    # Promote each isolate
    for isolate_id, accession_records in isolate_promotion_map.items():
        isolate = otu.get_isolate(isolate_id)

        if isolate is None:
            logger.warning("Isolate not found", isolate_id=str(isolate_id))
            continue

        accession_map = {}

        for predecessor_accession, refseq_record in accession_records.items():
            segment_id = assign_segment_id_to_record(refseq_record, otu.plan)
            if segment_id is None:
                logger.error(
                    "Segment does not match plan",
                    accession=refseq_record.accession_version,
                )
                continue

            from ref_builder.models.sequence import Sequence

            sequence = None
            for seq in isolate.sequences:
                if seq.accession.key == predecessor_accession:
                    sequence = seq
                    break

            if sequence is None:
                logger.warning(
                    "Predecessor sequence not found in isolate",
                    predecessor_accession=predecessor_accession,
                    isolate_id=str(isolate_id),
                )
                continue

            accession_map[sequence.accession] = Sequence(
                accession=refseq_record.accession_version,
                definition=refseq_record.definition,
                segment=segment_id,
                sequence=refseq_record.sequence,
            )

        if not accession_map:
            logger.warning(
                "No valid promotions for isolate", isolate_id=str(isolate_id)
            )
            continue

        try:
            repo.promote_isolate(otu.id, isolate_id, accession_map)
            promoted_accessions.update(
                seq.accession.key for seq in accession_map.values()
            )
        except ValueError as e:
            logger.error(
                "Promotion failed for isolate",
                error=str(e),
                isolate_id=str(isolate_id),
            )
            continue

        otu = repo.get_otu(otu.id)

    if promoted_accessions:
        log.info(
            "Sequences promoted.",
            promoted_accessions=sorted(promoted_accessions),
        )

    return promoted_accessions
