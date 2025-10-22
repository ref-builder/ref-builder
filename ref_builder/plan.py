from contextlib import suppress

from ref_builder.models.plan import Plan, SegmentName
from ref_builder.ncbi.models import NCBIGenbank


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
