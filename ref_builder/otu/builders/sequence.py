from pydantic import UUID4, BaseModel, field_serializer, field_validator

from ref_builder.utils import Accession


class SequenceBuilder(BaseModel):
    """Represents a mutable sequence in a Virtool reference repository."""

    id: UUID4
    """The sequence id."""

    accession: Accession
    """The sequence accession."""

    definition: str
    """The sequence definition."""

    legacy_id: str | None
    """A string based ID carried over from a legacy Virtool reference repository.

    It the sequence was not migrated from a legacy repository, this will be `None`.
    """

    sequence: str
    """The sequence."""

    segment: UUID4
    """The sequence segment."""

    @field_validator("accession", mode="before")
    @classmethod
    def convert_accession(cls: "SequenceBuilder", value: Accession | str) -> Accession:
        """Convert the accession to an Accession object."""
        if isinstance(value, Accession):
            return value

        if isinstance(value, str):
            return Accession.from_string(value)

        raise ValueError(f"Invalid type for accession: {type(value)}")

    @field_serializer("accession")
    @classmethod
    def serialize_accession(cls: "SequenceBuilder", accession: Accession) -> str:
        """Serialize the accession to a string."""
        return str(accession)