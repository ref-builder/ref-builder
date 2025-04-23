from collections import Counter
from uuid import UUID

from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.plan import (
    Plan,
    PlanConformationError,
    extract_segment_name_from_record_with_plan,
)


def assign_records_to_segments(
    records: list[NCBIGenbank], plan: Plan
) -> dict[UUID, NCBIGenbank]:
    """Assign genbank records segment IDs based on the passed plan.

    Assignment is based on segment naming. The function tries to normalize the segment
    name from the Genbank source table and match it with a segment in the plan. If a
    match to the normalized segment name is found in the plan, the record is assigned to
    the segment.

    Exceptions are raised when:

    * More than one segment is unnamed. Monopartite plans are allowed one unnamed
      segment.
    * A segment name is found in the records that doesn't exist in the plan.
    * A segment name is required by the plan but not found in the records.
    * There are duplicate segment names in the records.1

    The assigned segments are returned as a dictionary with segment IDs as records keyed
    by segment IDs.

    :param records: A list of Genbank records.
    :param plan: A plan.
    :return: A dictionary of segment IDs as keys and records as values.
    """
    seen_segment_names = Counter(
        extract_segment_name_from_record_with_plan(record, plan) for record in records
    )

    if seen_segment_names.total() > 1 and seen_segment_names[None]:
        raise PlanConformationError(
            "If a segment has no name, it must be the only segment. Only monopartite "
            "plans may have unnamed segments."
        )

    unassigned_segments = {segment.name: segment for segment in plan.segments}

    duplicate_segment_names = [
        segment_name for segment_name, count in seen_segment_names.items() if count > 1
    ]

    if duplicate_segment_names:
        # This also covers monopartite plans. Duplicate `None` segment names will be
        # caught here.
        joined = ", ".join(sorted(str(n) for n in duplicate_segment_names))
        raise PlanConformationError(
            f"Duplicate segment names found in records: {joined}."
        )

    segment_names_not_in_plan = [
        segment_name
        for segment_name in seen_segment_names
        if segment_name not in unassigned_segments
    ]

    if segment_names_not_in_plan:
        # This also covers monopartite plans. Plans are guaranteed to only one `None`
        # segment name and only when they are monopartite.
        joined = ", ".join(sorted(str(n) for n in segment_names_not_in_plan))
        raise PlanConformationError(f"Segment names not found in plan: {joined}.")

    segment_names_not_in_records = [
        segment.name
        for segment in plan.segments
        if segment.name not in seen_segment_names
    ]

    if segment_names_not_in_records:
        raise PlanConformationError(
            "Required segment names not found in records: ",
            f"{segment_names_not_in_records}",
        )

    return {
        unassigned_segments[
            extract_segment_name_from_record_with_plan(record, plan)
        ].id: record
        for record in records
    }


def assign_segment_id_to_record(
    record: NCBIGenbank,
    plan: Plan,
) -> UUID | None:
    """Return the ID of the matching segment, else None."""
    segment_name = extract_segment_name_from_record_with_plan(record, plan)

    if segment_name is None and plan.monopartite:
        return plan.segments[0].id

    for segment in plan.segments:
        if segment_name == segment.name:
            return segment.id

    return None
