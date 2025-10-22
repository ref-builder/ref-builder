from pydantic import UUID4

from ref_builder.models.otu import SequenceModel


class SequenceBuilder(SequenceModel):
    """Represents a mutable sequence in a Virtool reference repository."""

    segment: UUID4 | None
    """The sequence segment."""
