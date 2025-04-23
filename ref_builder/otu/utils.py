import math
import re
from collections import defaultdict
from enum import StrEnum
from uuid import UUID

import structlog

from ref_builder.models import Molecule
from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.plan import (
    Plan,
    PlanConformationError,
    Segment,
    SegmentRule,
    extract_segment_name_from_record,
)
from ref_builder.utils import (
    Accession,
    IsolateName,
    IsolateNameType,
    generate_natural_sort_key,
)

logger = structlog.get_logger()


class DeleteRationale(StrEnum):
    """Default strings delineating reasons for resource deletion."""

    USER = "Requested by user"
    REFSEQ = "Superceded by RefSeq"


class RefSeqConflictError(ValueError):
    """Raised when a potential RefSeq replacement is found."""

    def __init__(
        self,
        message: str,
        isolate_id: UUID,
        isolate_name: IsolateName,
        accessions: list[str],
    ) -> None:
        super().__init__(message)

        self.isolate_id = isolate_id

        self.isolate_name = isolate_name

        self.accessions = accessions


def get_segments_min_length(segments: list[Segment]) -> int:
    """Return the shortest minimum length from a list of segments."""
    shortest_segment = min(segments, key=lambda s: s.length)

    return math.floor(
        shortest_segment.length * (1.0 - shortest_segment.length_tolerance)
    )


def get_segments_max_length(segments: list[Segment]) -> int:
    """Return the longest maximum length from a list of segments."""
    longest_segment = max(segments, key=lambda s: s.length)

    return math.ceil(longest_segment.length * (1.0 + longest_segment.length_tolerance))


def check_sequence_length(sequence: str, segment_length: int, tolerance: float) -> bool:
    """Check if the sequence length is within acceptable segment length tolerance."""
    min_length = segment_length * (1.0 - tolerance)
    max_length = segment_length * (1.0 + tolerance)

    return min_length <= len(sequence) <= max_length


def create_plan_from_records(
    records: list[NCBIGenbank],
    length_tolerance: float,
    segments: list[Segment] | None = None,
) -> Plan | None:
    """Return a plan from a list of records representing an isolate."""
    if len(records) == 1:
        record = records[0]

        return Plan.new(
            segments=[
                Segment.new(
                    length=len(record.sequence),
                    length_tolerance=length_tolerance,
                    name=extract_segment_name_from_record(record),
                    rule=SegmentRule.REQUIRED,
                )
            ]
        )

    if len(group_genbank_records_by_isolate(records)) > 1:
        logger.warning("More than one isolate found. Cannot create plan.")
        return None

    if segments is None:
        segments = create_segments_from_records(
            records,
            rule=SegmentRule.REQUIRED,
            length_tolerance=length_tolerance,
        )

    if segments is not None:
        return Plan.new(segments=segments)

    return None


def fetch_records_from_accessions(
    accessions: list | set,
    blocked_accessions: set,
    ignore_cache: bool = False,
) -> list[NCBIGenbank]:
    """Fetch Genbank records from a list of accessions.

    Don't fetch accessions in ``blocked_accessions``.

    :param accessions: A list of accessions to fetch.
    :param blocked_accessions: A set of accessions to ignore.
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

    return NCBIClient(ignore_cache).fetch_genbank_records(eligible)


def create_segments_from_records(
    records: list[NCBIGenbank], rule: SegmentRule, length_tolerance: float
) -> list[Segment]:
    """Return a list of SegmentPlans."""
    if len(records) > 1 and not all(r.source.segment for r in records):
        raise ValueError("Segment name not found for multipartite OTU segment.")

    segments = [
        Segment.from_record(record, length_tolerance, rule) for record in records
    ]

    return sorted(segments, key=lambda s: generate_natural_sort_key(str(s.name)))


def get_molecule_from_records(records: list[NCBIGenbank]) -> Molecule:
    """Return relevant molecule metadata from one or more records.
    Molecule metadata is retrieved from the first RefSeq record in the list.
    If no RefSeq record is found in the list, molecule metadata is retrieved
    from record[0].
    """
    if not records:
        raise ValueError("No records given")

    # Assign first record as benchmark to start
    representative_record = records[0]

    if not representative_record.refseq:
        for record in records:
            if record.refseq:
                # Replace representative record with first RefSeq record found
                representative_record = record
                break

    return Molecule.model_validate(
        {
            "strandedness": representative_record.strandedness.value,
            "type": representative_record.moltype.value,
            "topology": representative_record.topology.value,
        },
    )


def group_genbank_records_by_isolate(
    records: list[NCBIGenbank],
) -> dict[IsolateName, dict[Accession, NCBIGenbank]]:
    """Indexes Genbank records by isolate name."""
    isolates = defaultdict(dict)

    for record in records:
        isolate_name = _extract_isolate_name_from_record(record)
        if isolate_name is None:
            # Assume this is a monopartite OTU and do not group.
            continue

        versioned_accession = Accession.from_string(record.accession_version)
        isolates[isolate_name][versioned_accession] = record

    return isolates


def parse_refseq_comment(comment: str) -> tuple[str, str]:
    """Parse a standard RefSeq comment."""
    if not comment:
        raise ValueError("Empty comment")

    pattern = re.compile(r"^(\w+ REFSEQ): [\w ]+. [\w ]+ (\w+).")

    if match := pattern.search(comment):
        return match.group(1), match.group(2)

    raise ValueError("Invalid RefSeq comment")


def _extract_isolate_name_from_record(record: NCBIGenbank) -> IsolateName | None:
    """Get the isolate name from a Genbank record."""
    if record.source.model_fields_set.intersection(
        {IsolateNameType.ISOLATE, IsolateNameType.STRAIN, IsolateNameType.CLONE},
    ):
        for source_type in IsolateNameType:
            if source_type in record.source.model_fields_set:
                return IsolateName(
                    type=IsolateNameType(source_type),
                    value=record.source.model_dump()[source_type],
                )

    return None
