import pytest
from syrupy.assertion import SnapshotAssertion
from syrupy.filters import props

from ref_builder.models.plan import Plan
from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBISourceMolType
from ref_builder.otu import assign_records_to_segments
from ref_builder.plan import create_plan_from_records
from ref_builder.repo import Repo
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory
from tests.fixtures.mock_ncbi_client import MockNCBIClient
from tests.fixtures.utils import uuid_matcher


class TestCreatePlanFromRecords:
    """Test `create_plan_from_records` function."""

    def test_monopartite(
        self,
        ncbi_genbank_factory: NCBIGenbankFactory,
        snapshot: SnapshotAssertion,
    ):
        """Test that a monopartite plan is created from a single Genbank record."""
        record = ncbi_genbank_factory.build()

        plan = create_plan_from_records([record], 0.03)

        assert isinstance(plan, Plan)
        assert len(plan.segments) == 1
        assert plan.model_dump() == snapshot(matcher=uuid_matcher)

    def test_multipartite(
        self,
        scratch_ncbi_client: NCBIClient,
        snapshot: SnapshotAssertion,
    ):
        """Test that a multipartite plan is created from multiple Genbank records."""
        records = scratch_ncbi_client.fetch_genbank_records(
            [
                "NC_010314",
                "NC_010315",
                "NC_010316",
                "NC_010317",
                "NC_010318",
                "NC_010319",
            ]
        )

        plan = create_plan_from_records(records, length_tolerance=0.03)

        assert isinstance(plan, Plan)
        assert plan.model_dump() == snapshot(matcher=uuid_matcher)

        segment_names = [str(segment.name) for segment in plan.segments]

        # Make sure the segment names are ordered alphabetically.
        assert segment_names == sorted(segment_names)

    def test_numeric_sorting(self, scratch_ncbi_client: NCBIClient):
        """Test that segment names are sorted numerically."""
        records = scratch_ncbi_client.fetch_genbank_records(
            ["NC_010314", "NC_010315", "NC_010316"]
        )

        for record, segment_name in zip(
            records, ["RNA 3", "RNA 10", "RNA 1"], strict=True
        ):
            record.source.segment = segment_name

        plan = create_plan_from_records(records, length_tolerance=0.03)

        assert plan
        assert [str(segment.name) for segment in plan.segments] == [
            "RNA 1",
            "RNA 3",
            "RNA 10",
        ]


class TestAssignRecordsToSegments:
    """Test whether a list of records can be assigned to plan segments."""

    def test_ok(
        self,
        mock_ncbi_client: MockNCBIClient,
        ncbi_genbank_factory: NCBIGenbankFactory,
        ncbi_source_factory: NCBISourceFactory,
        scratch_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test that a list of records generated to match the OTU plan."""
        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        )

        assert otu

        records = [
            ncbi_genbank_factory.build(
                source=ncbi_source_factory.build(
                    mol_type=NCBISourceMolType.from_molecule(otu.molecule),
                    segment=str(segment.name),
                    organism=otu.name,
                    taxid=otu.taxid,
                )
            )
            for segment in otu.plan.required_segments
        ]

        assigned_records = assign_records_to_segments(records, otu.plan)

        segment_names_by_id = {
            segment.id: segment.name for segment in otu.plan.segments
        }

        assert {
            (segment_names_by_id[segment_id], record.accession, record.source.segment)
            for segment_id, record in assigned_records.items()
        } == snapshot()

    def test_names_not_in_plan(
        self,
        mock_ncbi_client: MockNCBIClient,
        ncbi_genbank_factory: NCBIGenbankFactory,
        scratch_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test that a randomly generated list of records raises a ValueError."""
        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        )

        assert otu

        records = ncbi_genbank_factory.build_from_plan(otu.plan)

        # The example OTU uses a naming pattern like "RNA B".
        records[1].source.segment = "DNA 1"

        with pytest.raises(ValueError, match="Segment names not found in plan:") as e:
            assign_records_to_segments(records, otu.plan)

        assert str(e.value) == snapshot(exclude=props("id"))
