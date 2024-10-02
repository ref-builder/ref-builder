import re

from pydantic import UUID4, BaseModel, computed_field
from pydantic.dataclasses import dataclass

from ref_builder.models import Molecule
from ref_builder.ncbi.models import NCBIGenbank, NCBISourceMolType

SIMPLE_NAME_PATTERN = re.compile(r"([A-Za-z0-9])+")
"""Regex pattern for parsing segment name strings with no prefix."""

COMPLEX_NAME_PATTERN = re.compile(r"([A-Za-z]+)[-_ ]+([A-Za-z0-9]+)")
"""Regex pattern for parsing segment name strings consisting of a prefix and a key."""


@dataclass(frozen=True)
class SegmentName:
    """A normalized segment name. Can be used as a key."""

    prefix: str
    """The prefix of the segment name."""

    key: str
    """The identifying key portion of the segment name."""

    def __str__(self) -> str:
        """Return the segment name as a formatted string."""
        return f"{self.prefix} {self.key}"


class Segment(BaseModel):
    """The metadata of the segment."""

    id: UUID4
    """The unique id number of this segment"""

    name: str
    """The name of the segment"""

    required: bool
    """Whether this segment must be present in all additions."""

    length: int | None = None


class OTUSchema(BaseModel):
    """A schema for the intended data."""

    molecule: Molecule
    """The molecular metadata for this OTU."""

    segments: list[Segment]
    """The segments contained in this OTU."""

    @computed_field
    def multipartite(self) -> bool:
        """Is true if multiple sequences are required to form a complete isolate."""
        return len(self.segments) > 1


def determine_segment_prefix(moltype: NCBISourceMolType) -> str:
    """Return an acceptable SegmentName prefix corresponding to
    the given NCBISourceMolType.
    """
    if moltype in (
        NCBISourceMolType.GENOMIC_DNA,
        NCBISourceMolType.OTHER_DNA,
        NCBISourceMolType.UNASSIGNED_DNA,
    ):
        return "DNA"

    prefix = None

    match moltype:
        case NCBISourceMolType.GENOMIC_RNA:
            prefix = "RNA"
        case NCBISourceMolType.MRNA:
            prefix = "mRNA"
        case NCBISourceMolType.TRNA:
            prefix = "tRNA"
        case NCBISourceMolType.TRANSCRIBED_RNA:
            prefix = "RNA"
        case NCBISourceMolType.VIRAL_CRNA:
            prefix = "cRNA"
        case NCBISourceMolType.OTHER_RNA:
            prefix = "RNA"

    if prefix:
        return prefix

    raise ValueError(f"{moltype} may not be a valid NCBISourceMolType.")


def parse_segment_name(raw: str) -> SegmentName:
    """Parse a SegmentName from a raw string."""
    segment_name_parse = COMPLEX_NAME_PATTERN.fullmatch(raw)
    if segment_name_parse:
        return SegmentName(
            prefix=segment_name_parse.group(1),
            key=segment_name_parse.group(2),
        )

    raise ValueError(f"{raw} is not a valid segment name")


def get_multipartite_segment_name(record: NCBIGenbank) -> SegmentName:
    """Get a multipartite segment name from the record."""
    if SIMPLE_NAME_PATTERN.fullmatch(record.source.segment):
        return SegmentName(
            prefix=record.moltype, key=record.source.segment,
        )

    return parse_segment_name(record.source.segment)
