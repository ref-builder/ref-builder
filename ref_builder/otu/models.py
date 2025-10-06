
from pydantic import UUID4, BaseModel, field_serializer, field_validator

from ref_builder.models import Molecule
from ref_builder.plan import Plan
from ref_builder.utils import (
    Accession,
    IsolateName,
)


class SequenceModel(BaseModel):
    """A class representing the fields of a sequence."""

    id: UUID4
    """The sequence id."""

    accession: Accession
    """The sequence accession."""

    definition: str
    """The sequence definition."""

    sequence: str
    """The sequence."""

    segment: UUID4
    """The sequence segment."""

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


class IsolateModel(BaseModel):
    """A class representing the fields of an isolate."""

    id: UUID4
    """The isolate id."""

    name: IsolateName | None
    """The isolate's name."""


class OTUModel(BaseModel):
    """A class representing the fields of an OTU."""

    id: UUID4
    """The OTU id."""

    acronym: str
    """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

    excluded_accessions: set[str]
    """A set of accessions that should not be retrieved in future fetch operations."""

    molecule: Molecule
    """The type of molecular information contained in this OTU."""

    name: str
    """The name of the OTU (eg. TMV for Tobacco mosaic virus)"""

    plan: Plan
    """The plan for the OTU."""

    taxid: int
    """The NCBI Taxonomy id for this OTU."""
