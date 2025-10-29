"""Tests for Lineage models."""

import pytest
from pydantic import ValidationError

from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.ncbi.models import NCBIRank


class TestLineage:
    """Test the Lineage model."""

    def test_synonyms_single_taxon(self):
        """Test synonyms with a single species taxon."""
        lineage = Lineage(
            taxa=[
                Taxon(
                    id=1,
                    name="Tobacco mosaic virus",
                    parent=None,
                    rank=NCBIRank.SPECIES,
                    other_names=TaxonOtherNames(acronym=["TMV"], synonyms=["TMV-1"]),
                ),
            ]
        )

        assert lineage.synonyms == {"Tobacco mosaic virus", "TMV", "TMV-1"}

    def test_synonyms_multiple_taxa(self):
        """Test synonyms with multiple taxa in lineage."""
        lineage = Lineage(
            taxa=[
                Taxon(
                    id=3432891,
                    name="Tobamovirus tabaci",
                    parent=None,
                    rank=NCBIRank.SPECIES,
                    other_names=TaxonOtherNames(
                        acronym=[], synonyms=["Tobacco mosaic tobamovirus"]
                    ),
                ),
                Taxon(
                    id=12242,
                    name="Tobacco mosaic virus",
                    parent=3432891,
                    rank=NCBIRank.NO_RANK,
                    other_names=TaxonOtherNames(acronym=["TMV"], synonyms=["TMV-U1"]),
                ),
            ]
        )

        expected = {
            "Tobacco mosaic virus",
            "TMV",
            "TMV-U1",
            "Tobamovirus tabaci",
            "Tobacco mosaic tobamovirus",
        }
        assert lineage.synonyms == expected

    def test_synonyms_empty_other_names(self):
        """Test synonyms when other_names are empty."""
        lineage = Lineage(
            taxa=[
                Taxon(
                    id=1,
                    name="Test virus",
                    parent=None,
                    rank=NCBIRank.SPECIES,
                    other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                ),
            ]
        )

        assert lineage.synonyms == {"Test virus"}

    def test_synonyms_deduplication(self):
        """Test that synonyms are deduplicated across taxa."""
        lineage = Lineage(
            taxa=[
                Taxon(
                    id=2,
                    name="Virus A",
                    parent=None,
                    rank=NCBIRank.SPECIES,
                    other_names=TaxonOtherNames(acronym=["VA"], synonyms=[]),
                ),
                Taxon(
                    id=1,
                    name="Virus A",
                    parent=2,
                    rank=NCBIRank.NO_RANK,
                    other_names=TaxonOtherNames(acronym=["VA"], synonyms=["Virus A"]),
                ),
            ]
        )

        assert lineage.synonyms == {"Virus A", "VA"}

    def test_name_property(self):
        """Test that name returns the species taxon name."""
        lineage = Lineage(
            taxa=[
                Taxon(
                    id=2,
                    name="Species Y",
                    parent=None,
                    rank=NCBIRank.SPECIES,
                    other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                ),
                Taxon(
                    id=1,
                    name="Strain X",
                    parent=2,
                    rank=NCBIRank.NO_RANK,
                    other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                ),
            ]
        )

        assert lineage.name == "Species Y"

    def test_acronym_property(self):
        """Test that acronym returns first acronym from root to leaf."""
        lineage = Lineage(
            taxa=[
                Taxon(
                    id=2,
                    name="Species Y",
                    parent=None,
                    rank=NCBIRank.SPECIES,
                    other_names=TaxonOtherNames(acronym=["SY"], synonyms=[]),
                ),
                Taxon(
                    id=1,
                    name="Strain X",
                    parent=2,
                    rank=NCBIRank.NO_RANK,
                    other_names=TaxonOtherNames(acronym=["SX"], synonyms=[]),
                ),
            ]
        )

        assert lineage.acronym == "SY"

    def test_acronym_property_empty(self):
        """Test that acronym returns empty string when no acronyms exist."""
        lineage = Lineage(
            taxa=[
                Taxon(
                    id=1,
                    name="Test virus",
                    parent=None,
                    rank=NCBIRank.SPECIES,
                    other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                ),
            ]
        )

        assert lineage.acronym == ""

    def test_validation_non_species_root(self):
        """Test that validation fails if root is not species rank."""
        with pytest.raises(ValidationError, match="Root must be species rank"):
            Lineage(
                taxa=[
                    Taxon(
                        id=1,
                        name="Test",
                        parent=None,
                        rank=NCBIRank.NO_RANK,
                        other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                    ),
                ]
            )

    def test_validation_root_has_parent(self):
        """Test that validation fails if root has a parent."""
        with pytest.raises(ValidationError, match="Species root must have parent=None"):
            Lineage(
                taxa=[
                    Taxon(
                        id=1,
                        name="Test",
                        parent=999,
                        rank=NCBIRank.SPECIES,
                        other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                    ),
                ]
            )

    def test_validation_non_contiguous_taxa(self):
        """Test that validation fails if taxa reference invalid parents."""
        with pytest.raises(ValidationError, match="which is not in the lineage"):
            Lineage(
                taxa=[
                    Taxon(
                        id=2,
                        name="Parent",
                        parent=None,
                        rank=NCBIRank.SPECIES,
                        other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                    ),
                    Taxon(
                        id=1,
                        name="Child",
                        parent=999,
                        rank=NCBIRank.NO_RANK,
                        other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                    ),
                ]
            )
