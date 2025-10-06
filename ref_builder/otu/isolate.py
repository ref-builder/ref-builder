from uuid import UUID

from structlog import get_logger

from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.builders.sequence import SequenceBuilder
from ref_builder.otu.promote import promote_otu_accessions_from_records
from ref_builder.otu.utils import (
    assign_records_to_segments,
    fetch_records_from_accessions,
    group_genbank_records_by_isolate,
    parse_refseq_comment,
)
from ref_builder.plan import PlanConformationError
from ref_builder.repo import Repo
from ref_builder.utils import IsolateName

logger = get_logger("otu.isolate")


def add_genbank_isolate(
    repo: Repo,
    otu: OTUBuilder,
    accessions: list[str],
    ignore_cache: bool = False,
) -> IsolateBuilder | None:
    """Take a list of accessions that make up a new isolate and check that they make up
     a valid isolate before adding new isolate to the OTU.

    Download the GenBank records, categorize into an isolate bin and pass
    the isolate name and records to the add method.
    """
    otu_logger = logger.bind(taxid=otu.taxid, otu_id=str(otu.id), name=otu.name)

    records = fetch_records_from_accessions(
        accessions, otu.blocked_accessions, ignore_cache
    )

    if not records:
        return None

    if any(record.refseq for record in records) and not all(
        record.refseq for record in records
    ):
        raise ValueError(
            "Cannot mix RefSeq and non-RefSeq sequences in multipartite isolate."
        )

    binned_records = group_genbank_records_by_isolate(records)

    if len(binned_records) != 1:
        otu_logger.error("More than one isolate name found in requested accession.")
        return None

    isolate_name, isolate_records = next(iter(binned_records.items()))

    if (isolate_id := otu.get_isolate_id_by_name(isolate_name)) is not None:
        otu_logger.warning(
            "Isolate name already exists in this OTU.", name=isolate_name
        )

        if all(record.refseq for record in isolate_records.values()):
            otu_logger.info("Attempting to promote sequences to RefSeq")

            promoted_accessions = promote_otu_accessions_from_records(
                repo, otu, list(isolate_records.values())
            )

            if promoted_accessions:
                otu = repo.get_otu(otu.id)
                isolate = otu.get_isolate(isolate_id)
                otu_logger.info(
                    "Sequences promoted",
                    promoted_accessions=sorted(promoted_accessions),
                )
                return isolate

            otu_logger.warning("No promotable sequences found")

        return None

    with repo.use_transaction() as transaction:
        try:
            isolate = create_isolate(
                repo,
                otu,
                isolate_name,
                records,
            )
        except ValueError:
            otu_logger.exception()

        if isolate:
            return isolate

        transaction.abort()

    return None


def add_unnamed_isolate(
    repo: Repo,
    otu: OTUBuilder,
    accessions: list[str],
    ignore_cache: bool = False,
) -> IsolateBuilder | None:
    """Create an unnamed isolate from a list of accessions.

    Download the GenBank records and pass the isolate name and records to the add
    method.
    """
    records = fetch_records_from_accessions(
        accessions, otu.blocked_accessions, ignore_cache
    )

    if not records:
        return None

    with repo.use_transaction() as transaction:
        isolate = create_isolate(repo, otu, None, records)

        if isolate:
            return isolate

        transaction.abort()

    return None


def add_and_name_isolate(
    repo: Repo,
    otu: OTUBuilder,
    accessions: list[str],
    isolate_name: IsolateName,
    ignore_cache: bool = False,
) -> IsolateBuilder | None:
    """Take a list of accessions that make up a new isolate and a preferred
    isolate name, then add a new isolate to the OTU.

    Download the GenBank records and pass the isolate name and records to
    the add method.
    """
    records = fetch_records_from_accessions(
        accessions, otu.blocked_accessions, ignore_cache
    )

    if not records:
        return None

    with repo.use_transaction() as transaction:
        isolate = create_isolate(repo, otu, isolate_name, records)

        if isolate:
            return isolate

        transaction.abort()

    return None


def create_isolate(
    repo: Repo,
    otu: OTUBuilder,
    isolate_name: IsolateName | None,
    records: list[NCBIGenbank],
) -> IsolateBuilder | None:
    """Take a dictionary of records keyed by segment name and add them to the OTU."""
    log = get_logger(
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

    isolate = repo.create_isolate(
        otu.id,
        name=isolate_name,
    )

    for segment_id, record in assigned.items():
        if (sequence := otu.get_sequence_by_accession(record.accession)) is None:
            sequence = repo.create_sequence(
                otu.id,
                accession=record.accession_version,
                definition=record.definition,
                segment=segment_id,
                sequence=record.sequence,
            )

        repo.link_sequence(otu.id, isolate.id, sequence.id)

        if record.refseq:
            _, old_accession = parse_refseq_comment(record.comment)

            repo.exclude_accession(
                otu.id,
                old_accession,
            )

    log.info("Isolate created", id=str(isolate.id))

    return repo.get_isolate(isolate.id)


def create_sequence_from_record(
    repo: Repo,
    otu: OTUBuilder,
    record: NCBIGenbank,
    segment_id: UUID,
) -> SequenceBuilder | None:
    """Take a NCBI Nucleotide record and create a new sequence."""
    return repo.create_sequence(
        otu.id,
        accession=record.accession_version,
        definition=record.definition,
        segment=segment_id,
        sequence=record.sequence,
    )
