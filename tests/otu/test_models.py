"""Tests for OTU models."""

import warnings
from uuid import uuid4

import pytest
from pydantic import ValidationError

from ref_builder.isolate import IsolateNameType
from ref_builder.models.accession import Accession
from ref_builder.models.isolate import IsolateName
from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.models.molecule import Molecule, MoleculeType, Strandedness, Topology
from ref_builder.models.plan import Plan, Segment, SegmentRule
from ref_builder.ncbi.models import NCBIRank
from ref_builder.otu.validators.isolate import Isolate
from ref_builder.otu.validators.otu import OTU, OTUBase
from ref_builder.otu.validators.sequence import Sequence
from ref_builder.warnings import PlanWarning


class TestSequence:
    """Test the ``Sequence`` model which is used for complete validation of sequences."""

    def test_ok(self):
        """Test that a valid sequence passes validation."""
        assert Sequence.model_validate(
            {
                "id": uuid4(),
                "accession": Accession("NC_001234", 1),
                "definition": "Test virus, complete genome",
                "segment": uuid4(),
                "sequence": "ATCGATCGATCG",
            }
        )

    def test_invalid_accession(self):
        """Test that an invalid accession fails validation."""
        bad_sequence_data = {
            "id": uuid4(),
            "accession": Accession("BADSEQUENCE", 1),
            "definition": "Test virus, complete genome",
            "segment": uuid4(),
            "sequence": "ATCGATCGATCG",
        }

        try:
            Sequence.model_validate(bad_sequence_data)
        except ValidationError as e:
            for error in e.errors():
                assert "accession" in error["loc"]
                assert (
                    "Accession BADSEQUENCE.1 does not match a valid accession pattern"
                    in error["msg"]
                )


class TestIsolate:
    """Test the ``Isolate`` model which is used for complete validation of isolates."""

    def test_accession_consistency_warning(self):
        """Test that a warning is raised if accession provenances are mixed."""
        segment_id = uuid4()

        # Create isolate with mixed RefSeq and GenBank accessions
        mixed_isolate_data = {
            "id": uuid4(),
            "name": None,
            "sequences": [
                {
                    "id": uuid4(),
                    "accession": Accession("NC_000001", 1),  # RefSeq
                    "definition": "Test virus segment A",
                    "segment": segment_id,
                    "sequence": "ATCGATCGATCG",
                },
                {
                    "id": uuid4(),
                    "accession": Accession("BD000001", 1),  # GenBank
                    "definition": "Test virus segment B",
                    "segment": segment_id,
                    "sequence": "GCTAGCTAGCTA",
                },
            ],
        }

        with warnings.catch_warnings(record=True) as warning_list:
            Isolate.model_validate(mixed_isolate_data)

        assert len(warning_list) > 0
        warning_msg = warning_list[0]

        assert warning_msg.category.__name__ == "IsolateInconsistencyWarning"
        assert (
            "Combination of RefSeq and non-RefSeq sequences found in multipartite isolate"
            in str(warning_msg.message)
        )
        assert "NC_000001" in str(warning_msg)
        assert "BD000001" in str(warning_msg)


class TestOTU:
    """Test the ``OTU`` model which is used for complete validation of OTUs."""

    otu: OTUBase

    @pytest.fixture(autouse=True)
    def _build_otu(self):
        """Build static OTU data for testing."""
        segment_id = uuid4()
        plan_id = uuid4()

        self.otu = OTUBase(
            id=uuid4(),
            acronym="TMV",
            excluded_accessions=set(),
            lineage=Lineage(
                taxa=[
                    Taxon(
                        id=12242,
                        name="Tobacco mosaic virus",
                        parent=3432891,
                        rank=NCBIRank.NO_RANK,
                        other_names=TaxonOtherNames(acronym=["TMV"], synonyms=[]),
                    ),
                    Taxon(
                        id=3432891,
                        name="Tobamovirus tabaci",
                        parent=None,
                        rank=NCBIRank.SPECIES,
                        other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                    ),
                ]
            ),
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MoleculeType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            plan=Plan(
                id=plan_id,
                segments=[
                    Segment(
                        id=segment_id,
                        length=20,
                        length_tolerance=0.1,
                        name=None,
                        rule=SegmentRule.REQUIRED,
                    )
                ],
            ),
            taxid=12242,
            isolates=[
                {
                    "id": uuid4(),
                    "name": IsolateName(type=IsolateNameType.ISOLATE, value="TMV-001"),
                    "sequences": [
                        {
                            "id": uuid4(),
                            "accession": Accession("NC_001367", 1),
                            "definition": "Tobacco mosaic virus, complete genome",
                            "segment": segment_id,
                            "sequence": "ATCGATCGATCGATCGATCG",
                        }
                    ],
                },
                {
                    "id": uuid4(),
                    "name": IsolateName(type=IsolateNameType.ISOLATE, value="TMV-002"),
                    "sequences": [
                        {
                            "id": uuid4(),
                            "accession": Accession("AF395128", 1),
                            "definition": "Tobacco mosaic virus isolate TMV-017",
                            "segment": segment_id,
                            "sequence": "GCTAGCTAGCTAGCTAGCTA",
                        }
                    ],
                },
            ],
        )

    def test_ok(self):
        """Test that a valid OTU passes validation."""
        assert OTU.model_validate(self.otu.model_dump())

    def test_synonyms(self):
        """Test that synonyms returns all names from lineage."""
        assert self.otu.synonyms == {"Tobacco mosaic virus", "TMV", "Tobamovirus tabaci"}

    def test_no_required_segments(self):
        """Test that OTU raises a warning if initialized without required segments."""
        segment_id = uuid4()
        plan_id = uuid4()

        otu_data = OTUBase(
            id=uuid4(),
            acronym="TMV",
            excluded_accessions=set(),
            lineage=Lineage(
                taxa=[
                    Taxon(
                        id=3432891,
                        name="Tobamovirus tabaci",
                        parent=None,
                        rank=NCBIRank.SPECIES,
                        other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                    ),
                ]
            ),
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MoleculeType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            plan=Plan(
                id=plan_id,
                segments=[
                    Segment(
                        id=segment_id,
                        length=20,
                        length_tolerance=0.1,
                        name=None,
                        rule=SegmentRule.RECOMMENDED,  # Not required
                    )
                ],
            ),
            taxid=3432891,
            isolates=[
                {
                    "id": uuid4(),
                    "name": IsolateName(type=IsolateNameType.ISOLATE, value="TMV-001"),
                    "sequences": [
                        {
                            "id": uuid4(),
                            "accession": Accession("NC_001367", 1),
                            "definition": "Tobacco mosaic virus, complete genome",
                            "segment": segment_id,
                            "sequence": "ATCGATCGATCGATCGATCG",
                        }
                    ],
                }
            ],
        )

        assert not otu_data.plan.required_segments

        with pytest.warns(PlanWarning):
            OTU.model_validate(otu_data.model_dump())

    def test_excluded_accessions(self):
        """Test that validation fails if the OTU includes accessions that are included
        in ``excluded_accessions``.
        """
        accession = self.otu.isolates[0].sequences[0].accession.key

        self.otu.excluded_accessions.add(accession)

        with pytest.raises(
            ValueError, match=f"Excluded accessions found in the OTU: {accession}"
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_no_isolates(self):
        """Test that validation fails if the OTU has no isolates."""
        self.otu.isolates = []

        with pytest.raises(
            ValueError, match="List should have at least 1 item after validation, not 0"
        ):
            assert OTU.model_validate(self.otu.model_dump())

    def test_unique_isolate_names(self):
        """Test that validation fails if the OTU has duplicate isolate names."""
        duplicate_name = self.otu.isolates[1].name = self.otu.isolates[0].name

        with pytest.raises(
            ValueError,
            match=f"Isolate names must be unique. Non-unique names: {duplicate_name}",
        ):
            assert OTU.model_validate(self.otu.model_dump())
