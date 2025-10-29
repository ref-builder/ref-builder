from typing import Annotated

from pydantic import BaseModel, ConfigDict, Field, field_validator

from ref_builder.ncbi.models import NCBIRank


class TaxonOtherNames(BaseModel):
    """Alternate names and identifiers for a taxon."""

    model_config = ConfigDict(frozen=True)

    acronym: list[str]
    """Official acronyms from NCBI taxonomy (e.g., TMV, BBSV)."""

    synonyms: list[str]
    """Alternative valid names for the taxon (parsed from NCBI's EquivalentName field)."""


class Taxon(BaseModel):
    """A single taxon in a lineage."""

    model_config = ConfigDict(frozen=True)

    id: int
    name: str
    parent: int | None
    rank: NCBIRank
    other_names: TaxonOtherNames


class Lineage(BaseModel):
    """An OTU subspecific lineage.

    Ordered from species (root) to isolates/strains (leaves).
    """

    model_config = ConfigDict(frozen=True)

    taxa: Annotated[list[Taxon], Field(min_length=1)]

    @field_validator("taxa", mode="after")
    @classmethod
    def validate_taxa(cls, taxa: list[Taxon]) -> list[Taxon]:
        """Validate the taxa in this lineage.

        - Ensure species is the root (first position with parent=None).
        - Ensure all taxa reference valid parents within the lineage.
        """
        if taxa[0].rank != NCBIRank.SPECIES:
            raise ValueError(f"Root must be species rank, got '{taxa[0].rank}'")

        if taxa[0].parent is not None:
            raise ValueError(
                f"Species root must have parent=None, got parent={taxa[0].parent}"
            )

        taxon_ids = {taxon.id for taxon in taxa}

        for i in range(1, len(taxa)):
            if taxa[i].parent not in taxon_ids:
                raise ValueError(
                    f"Taxon {taxa[i].name} (id={taxa[i].id}) has parent={taxa[i].parent} "
                    f"which is not in the lineage"
                )

        return taxa

    @property
    def acronym(self) -> str:
        """Return first acronym found in lineage, searching from root to leaf.

        Searches taxa from species (root at taxa[0]) down to leaf (taxa[-1]).
        Returns empty string if no acronyms found.
        """
        for taxon in self.taxa:
            if taxon.other_names.acronym:
                return taxon.other_names.acronym[0]
        return ""

    @property
    def name(self) -> str:
        """Return the species name from the lineage.

        The root taxon (first position) is always species rank (enforced by validator).
        """
        return self.taxa[0].name

    @property
    def synonyms(self) -> set[str]:
        """Return all possible names derived from the lineage.

        Collects all taxon names, acronyms, and synonyms from all taxa in the lineage.
        """
        synonyms = set()

        for taxon in self.taxa:
            synonyms.add(taxon.name)
            synonyms.update(taxon.other_names.acronym)
            synonyms.update(taxon.other_names.synonyms)

        return synonyms
