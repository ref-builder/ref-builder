from pydantic import (
    UUID4,
    ConfigDict,
    Field,
    field_serializer,
    field_validator,
)

from ref_builder.models.accession import Accession
from ref_builder.models.otu import SequenceModel
from ref_builder.utils import is_accession_key_valid, is_refseq


class SequenceBase(SequenceModel):
    """A class representing a sequence with basic validation."""

    model_config = ConfigDict(validate_assignment=True)

    id: UUID4
    """The sequence id."""

    accession: Accession
    """The sequence accession."""

    definition: str = Field(min_length=1)
    """The sequence definition."""

    sequence: str = Field(pattern=r"[ATGCURYKMSWBDHVN]+")
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
