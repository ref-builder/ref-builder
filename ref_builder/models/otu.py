from uuid import UUID

from pydantic import UUID4, BaseModel, field_serializer, field_validator

from ref_builder.models.accession import Accession
from ref_builder.models.molecule import Molecule
from ref_builder.models.plan import Plan


class OTUMinimal(BaseModel):
    """A minimal representation of an OTU."""

    acronym: str
    id: UUID
    name: str
    taxid: int


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
