from uuid import UUID

from pydantic import UUID4, BaseModel, field_serializer, field_validator

from ref_builder.models.accession import Accession
from ref_builder.models.lineage import Lineage
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

    excluded_accessions: set[str]
    """A set of accessions that should not be retrieved in future fetch operations."""

    lineage: Lineage
    """The taxonomic lineage from species down to the target taxon."""

    molecule: Molecule
    """The type of molecular information contained in this OTU."""

    plan: Plan
    """The plan for the OTU."""

    @property
    def acronym(self) -> str:
        """The OTU acronym computed from the lineage."""
        return self.lineage.acronym

    @property
    def name(self) -> str:
        """The OTU name computed from the species-level taxon in the lineage."""
        return self.lineage.name

    @property
    def synonyms(self) -> set[str]:
        """All possible names derived from the OTU's lineage.

        Returns all taxon names, acronyms, and synonyms from the lineage.
        """
        return self.lineage.synonyms

    @property
    def taxid(self) -> int:
        """The NCBI Taxonomy ID for this OTU (species-level taxon from lineage)."""
        return self.lineage.taxa[-1].id


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
