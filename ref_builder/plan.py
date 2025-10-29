import math
from contextlib import suppress

from ref_builder.errors import PlanCreationError
from ref_builder.models.plan import (
    Plan,
    Segment,
    SegmentName,
    SegmentRule,
    extract_segment_name_from_record,
)
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.ncbi.utils import group_genbank_records_by_isolate
from ref_builder.utils import generate_natural_sort_key


def extract_segment_name_from_record_with_plan(
    record: NCBIGenbank, plan: Plan
) -> SegmentName | None:
    """Extract a segment name from a Genbank record using an OTU plan.

    If the segment name can be parsed normally (eg. DNA A), it is returned. If there is
    no delimiter or the prefix is missing, the function tries to match the segment name
    with a segment in the plan.

    If no segment name can be extracted, `None` is returned.

    :param record: A Genbank record.
    :param plan: A plan.
    :return: A segment name or `None`.
    """
    if not record.source.segment:
        return None

    if (segment_name := SegmentName.from_string(record.source.segment)) is not None:
        return segment_name

    if not plan.monopartite:
        try:
            plan_keys_and_prefixes = {
                segment.name.key: segment.name.prefix for segment in plan.segments
            }
        except AttributeError:
            raise ValueError("Multipartite plan contains unnamed segments")

        # Handle no prefix.
        with suppress(KeyError):
            return SegmentName(
                prefix=plan_keys_and_prefixes[record.source.segment],
                key=record.source.segment,
            )

    return None


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
        raise PlanCreationError("More than one isolate found. Cannot create plan.")

    if segments is None:
        segments = create_segments_from_records(
            records,
            rule=SegmentRule.REQUIRED,
            length_tolerance=length_tolerance,
        )

    return Plan.new(segments=segments)


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
