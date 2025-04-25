from pydantic import (
    UUID4,
    BaseModel,
    ConfigDict,
    field_serializer,
    field_validator,
)

from ref_builder.otu.models import SequenceModel
from ref_builder.utils import Accession, is_accession_key_valid, is_refseq


class SequenceBase(SequenceModel):
    """A class representing a sequence with basic validation."""

    model_config = ConfigDict(validate_assignment=True)

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

    @property
    def refseq(self) -> bool:
        """Return True if this sequence was sourced from NCBI's RefSeq database."""
        return is_refseq(self.accession.key)

    @field_validator("accession", mode="before")
    @classmethod
    def convert_accession(cls, value: Accession | str) -> Accession:
        """Convert the accession to an Accession object."""
        if isinstance(value, Accession):
            return value

        if isinstance(value, str):
            return Accession.from_string(value)

        raise ValueError(f"Invalid type for accession: {type(value)}")

    @field_serializer("accession")
    @classmethod
    def serialize_accession(cls, accession: Accession) -> str:
        """Serialize the accession to a string."""
        return str(accession)


class Sequence(SequenceBase):
    """A class representing a sequence with full validation."""

    @field_validator("accession", mode="after")
    @classmethod
    def check_accession_key_is_valid(cls, v: Accession) -> Accession:
        """Check the accession key against INSDC and NCBI RefSeq patterns."""
        if is_accession_key_valid(v.key):
            return v

        raise ValueError(f"Accession {v} does not match a valid accession pattern.")
