from uuid import uuid4

from ref_builder.isolate import IsolateNameType
from ref_builder.models.isolate import IsolateName
from ref_builder.models.molecule import Molecule, MoleculeType, Strandedness, Topology
from ref_builder.models.plan import Plan, Segment
from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.repo import Repo
from tests.fixtures.factories import IsolateFactory, OTUFactory


class TestSequence:
    """Test properties of SequenceBuilder."""

    def test_equivalence(self, scratch_repo: Repo, mock_ncbi_client):
        """Test that the == operator works correctly."""
        taxid = mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        accessions = ["DQ178614", "DQ178613", "DQ178610", "DQ178611"]

        otu = scratch_repo.get_otu_by_taxid(taxid)

        assert otu

        for accession in accessions:
            assert otu.get_sequence_by_accession(
                accession,
            ) == otu.get_sequence_by_accession(accession)


class TestIsolate:
    """Test properties of IsolateBuilder."""

    def test_no_sequences(self):
        """Test that an isolate initializes correctly with no sequences."""
        isolate = IsolateBuilder(
            id=uuid4(),
            name=IsolateName(type=IsolateNameType.ISOLATE, value="A"),
            sequences=[],
        )

        assert isolate.sequences == []

    def test_equivalence(self, scratch_repo: Repo, mock_ncbi_client):
        """Test that the == operator works correctly."""
        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        )

        assert otu

        for isolate in otu.isolates:
            assert isolate == otu.get_isolate(isolate.id)


class TestOTU:
    """Test properties of OTUBuilder."""

    def test_no_isolates(self):
        """Test that an isolate initializes correctly with no isolates."""
        otu = OTUBuilder(
            id=uuid4(),
            acronym="TMV",
            excluded_accessions=set(),
            isolates=[],
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MoleculeType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            plan=Plan.new(
                segments=[
                    Segment.new(
                        length=6395,
                        length_tolerance=0.03,
                        name=None,
                    )
                ]
            ),
            taxid=12242,
        )

        assert otu.isolates == []

    def test_equivalence(self, scratch_repo: Repo):
        """Test that the == operator works correctly."""
        taxid = 345184

        assert scratch_repo.get_otu_by_taxid(taxid) == scratch_repo.get_otu_by_taxid(
            taxid,
        )

    def test_get_sequence_id_hierarchy(self, scratch_repo: Repo, mock_ncbi_client):
        """Test that the isolate ID can be found from a sequence ID."""
        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        )

        assert otu

        sequence = otu.get_sequence_by_accession("DQ178610")

        assert sequence

        containing_isolate_ids = otu.get_isolate_ids_containing_sequence_id(sequence.id)

        assert len(containing_isolate_ids) == 1

        isolate_id = next(iter(containing_isolate_ids))

        assert otu.get_isolate(isolate_id) is not None

    def test_check_get_sequence_by_id_integrity(self):
        """Test that OTUBuilder.get_sequence() can retrieve every sequence ID
        within each constituent isolate from its own private index.
        """
        otu_ = OTUBuilder(**OTUFactory.build().model_dump())

        for _ in range(10):
            isolate = IsolateBuilder.model_validate(
                IsolateFactory.build_on_plan(otu_.plan).model_dump()
            )
            otu_.add_isolate(isolate)

        sequence_ids_in_isolates = {
            sequence.id for isolate in otu_.isolates for sequence in isolate.sequences
        }

        for sequence_id in sequence_ids_in_isolates:
            assert otu_.get_sequence_by_id(sequence_id) is not None

    def test_check_get_sequence_by_accession_integrity(self):
        """Test that OTUBuilder.get_sequence() can retrieve every accession
        within each constituent isolate.
        """
        otu = OTUBuilder.model_validate(OTUFactory.build().model_dump())

        for _ in range(10):
            isolate = IsolateBuilder.model_validate(
                IsolateFactory.build_on_plan(otu.plan).model_dump()
            )
            otu.add_isolate(isolate)

        for isolate in otu.isolates:
            for sequence in isolate.sequences:
                assert otu.get_sequence_by_accession(sequence.accession.key) == sequence

    def test_delete_isolate(self):
        """Test that OTUBuilder.delete_isolate() does not affect other isolates."""
        otu = OTUBuilder.model_validate(OTUFactory.build().model_dump())

        for _ in range(10):
            isolate = IsolateBuilder.model_validate(
                IsolateFactory.build_on_plan(otu.plan).model_dump()
            )
            otu.add_isolate(isolate)

        initial_isolates = {isolate.id: isolate for isolate in otu.isolates}
        deleted_id = list(initial_isolates.keys())[3]

        otu.delete_isolate(deleted_id)

        assert otu.get_isolate(deleted_id) is None

        for isolate_id, isolate in initial_isolates.items():
            if isolate_id != deleted_id:
                assert otu.get_isolate(isolate_id) == isolate
