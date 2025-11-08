"""Tests for OTU model."""

from uuid import uuid4

import pytest
from pydantic import ValidationError

from ref_builder.models.accession import Accession
from ref_builder.models.isolate import IsolateName, IsolateNameType
from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.models.molecule import Molecule, MoleculeType, Strandedness, Topology
from ref_builder.models.otu import OTU
from ref_builder.models.plan import Plan, Segment, SegmentRule
from ref_builder.ncbi.models import NCBIRank
from ref_builder.warnings import PlanWarning


class TestOTU:
    """Test the ``OTU`` model which is used for complete validation of OTUs."""

    otu: OTU

    @pytest.fixture(autouse=True)
    def setup(self):
        """Build static OTU data for testing."""
        segment_id = uuid4()
        plan_id = uuid4()

        self.otu = OTU.model_validate(
            {
                "id": uuid4(),
                "acronym": "TMV",
                "excluded_accessions": set(),
                "promoted_accessions": set(),
                "lineage": Lineage(
                    taxa=[
                        Taxon(
                            id=3432891,
                            name="Tobamovirus tabaci",
                            parent=None,
                            rank=NCBIRank.SPECIES,
                            other_names=TaxonOtherNames(acronym=[], synonyms=[]),
                        ),
                        Taxon(
                            id=12242,
                            name="Tobacco mosaic virus",
                            parent=3432891,
                            rank=NCBIRank.NO_RANK,
                            other_names=TaxonOtherNames(acronym=["TMV"], synonyms=[]),
                        ),
                    ]
                ),
                "molecule": Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MoleculeType.RNA,
                    topology=Topology.LINEAR,
                ),
                "name": "Tobacco mosaic virus",
                "plan": Plan(
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
                "taxid": 12242,
                "isolates": [
                    {
                        "id": uuid4(),
                        "name": IsolateName(
                            type=IsolateNameType.ISOLATE, value="TMV-001"
                        ),
                        "taxid": 12242,
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
                        "name": IsolateName(
                            type=IsolateNameType.ISOLATE, value="TMV-002"
                        ),
                        "taxid": 12242,
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
            }
        )

    def test_ok(self):
        """Test that a valid OTU passes validation."""
        assert OTU.model_validate(self.otu.model_dump())

    def test_synonyms(self):
        """Test that synonyms returns all names from lineage."""
        assert self.otu.synonyms == {
            "Tobacco mosaic virus",
            "TMV",
            "Tobamovirus tabaci",
        }

    def test_no_required_segments(self):
        """Test that OTU raises a warning if initialized without required segments."""
        segment_id = uuid4()
        plan_id = uuid4()

        otu = OTU.model_validate(
            {
                "id": uuid4(),
                "acronym": "TMV",
                "excluded_accessions": set(),
                "promoted_accessions": set(),
                "lineage": Lineage(
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
                "molecule": Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MoleculeType.RNA,
                    topology=Topology.LINEAR,
                ),
                "name": "Tobacco mosaic virus",
                "plan": Plan(
                    id=plan_id,
                    segments=[
                        Segment(
                            id=segment_id,
                            length=20,
                            length_tolerance=0.1,
                            name=None,
                            rule=SegmentRule.RECOMMENDED,
                        )
                    ],
                ),
                "taxid": 3432891,
                "isolates": [
                    {
                        "id": uuid4(),
                        "name": IsolateName(
                            type=IsolateNameType.ISOLATE, value="TMV-001"
                        ),
                        "taxid": 3432891,
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
            }
        )

        assert not otu.plan.required_segments

        with pytest.warns(PlanWarning):
            OTU.model_validate(otu.model_dump())

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
        otu_data = self.otu.model_dump()
        otu_data["isolates"] = []

        with pytest.raises(
            ValidationError,
            match="List should have at least 1 item after validation, not 0",
        ):
            OTU.model_validate(otu_data)

    def test_isolate_taxid_not_in_lineage(self):
        """Test that validation fails if an isolate taxid is not in the OTU lineage."""
        otu_data = self.otu.model_dump()
        otu_data["isolates"][0]["taxid"] = 99999

        with pytest.raises(
            ValidationError,
            match="Isolate taxid 99999 is not found in OTU lineage",
        ):
            OTU.model_validate(otu_data)

    def test_isolate_taxid_valid_species_level(self):
        """Test that validation passes when isolate taxid is the species-level taxid."""
        otu_data = self.otu.model_dump()
        otu_data["isolates"][0]["taxid"] = 3432891

        assert OTU.model_validate(otu_data)
