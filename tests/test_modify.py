from uuid import UUID, uuid4

import pytest
from pydantic import ValidationError
from syrupy.assertion import SnapshotAssertion
from syrupy.filters import props

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.modify import (
    add_segments_to_plan,
    allow_accessions_into_otu,
    delete_isolate_from_otu,
    exclude_accessions_from_otu,
    rename_plan_segment,
    replace_sequence_in_otu,
    set_plan,
    set_plan_length_tolerances,
)
from ref_builder.plan import (
    Plan,
    Segment,
    SegmentName,
    SegmentRule,
)
from ref_builder.repo import Repo
from ref_builder.services.otu import OTUService
from ref_builder.utils import IsolateName, IsolateNameType
from tests.fixtures.factories import IsolateFactory


def test_exclude_accessions(scratch_repo: Repo):
    """Test accession exclusion."""
    taxid = 345184

    otu_before = scratch_repo.get_otu_by_taxid(taxid)
    assert otu_before is not None

    assert not otu_before.excluded_accessions

    with scratch_repo.lock():
        exclude_accessions_from_otu(
            scratch_repo, otu_before, accessions=["DQ178608", "DQ178609"]
        )

    otu_after = scratch_repo.get_otu_by_taxid(taxid)
    assert otu_after is not None

    assert otu_after.excluded_accessions == {"DQ178608", "DQ178609"}


def test_allow_accessions(scratch_repo: Repo):
    taxid = 345184

    otu_initial = scratch_repo.get_otu_by_taxid(taxid)
    assert otu_initial is not None

    with scratch_repo.lock():
        exclude_accessions_from_otu(
            scratch_repo,
            otu=otu_initial,
            accessions={"DQ178608", "DQ178609"},
        )

    otu_before = scratch_repo.get_otu_by_taxid(taxid)
    assert otu_before is not None

    assert otu_before.excluded_accessions == {"DQ178608", "DQ178609"}

    with scratch_repo.lock():
        allow_accessions_into_otu(
            scratch_repo,
            otu=otu_before,
            accessions={"DQ178608"},
        )

    otu_after = scratch_repo.get_otu_by_taxid(taxid)
    assert otu_after is not None

    assert otu_after.excluded_accessions == {"DQ178609"}


class TestSetPlan:
    """Test functions that make changes to an OTU plan."""

    def test_ok(self, scratch_repo: Repo):
        """Test that an OTU's plan can be replaced."""
        otu_before = scratch_repo.get_otu_by_taxid(223262)
        assert otu_before is not None

        original_plan = otu_before.plan

        assert type(original_plan) is Plan

        new_plan = Plan.new(segments=original_plan.segments)

        new_plan.segments.append(
            Segment.new(
                length=2000,
                length_tolerance=scratch_repo.settings.default_segment_length_tolerance,
                name=SegmentName(prefix="DNA", key="C"),
                rule=SegmentRule.RECOMMENDED,
            ),
        )

        new_plan.segments.append(
            Segment.new(
                length=1000,
                length_tolerance=scratch_repo.settings.default_segment_length_tolerance,
                name=SegmentName(prefix="DNA", key="Z"),
                rule=SegmentRule.OPTIONAL,
            ),
        )
        with scratch_repo.lock():
            set_plan(scratch_repo, otu_before, new_plan)

        assert type(new_plan) is Plan

        otu_after = scratch_repo.get_otu(otu_before.id)
        assert otu_after is not None

        assert len(otu_after.plan.segments) == len(otu_before.plan.segments) + 2

        assert otu_after.plan == new_plan

    def test_rename_segment_ok(self, scratch_repo: Repo):
        """Test that a given plan segment can be renamed."""
        otu_before = scratch_repo.get_otu_by_taxid(223262)
        assert otu_before is not None

        first_segment_id = otu_before.plan.segments[0].id

        new_name = SegmentName(prefix="RNA", key="TestName")

        segment_before = otu_before.plan.get_segment_by_id(first_segment_id)
        assert segment_before is not None
        assert segment_before.name != new_name

        with scratch_repo.lock():
            rename_plan_segment(
                scratch_repo,
                otu_before,
                segment_id=first_segment_id,
                segment_name=SegmentName(prefix="RNA", key="TestName"),
            )

        otu_after = scratch_repo.get_otu_by_taxid(223262)
        assert otu_after is not None

        segment_after = otu_after.plan.get_segment_by_id(first_segment_id)
        assert segment_after is not None
        assert segment_after.name == new_name

    def test_rename_segment_fail(self, scratch_repo: Repo):
        """Test that an attempt to rename a nonexistent segment does not change the OTU."""
        otu_before = scratch_repo.get_otu_by_taxid(223262)
        assert otu_before is not None

        first_segment = otu_before.plan.segments[0]

        new_name = SegmentName(prefix="RNA", key="TestName")

        segment_before = otu_before.plan.get_segment_by_id(first_segment.id)
        assert segment_before is not None
        assert segment_before.name != new_name

        assert (
            rename_plan_segment(
                scratch_repo,
                otu_before,
                segment_id=uuid4(),
                segment_name=SegmentName(prefix="RNA", key="TestName"),
            )
            is None
        )

        otu_after = scratch_repo.get_otu_by_taxid(223262)
        assert otu_after is not None

        assert otu_after.plan.model_dump_json() == otu_before.plan.model_dump_json()

        segment_after = otu_after.plan.get_segment_by_id(first_segment.id)
        assert segment_after is not None
        assert segment_after.name == first_segment.name

    @pytest.mark.parametrize(
        ("initial_accessions", "new_accessions"),
        [
            (["MF062136", "MF062137"], ["MF062138"]),
            (["MF062136"], ["MF062137", "MF062138"]),
        ],
    )
    def test_add_segments_to_plan_ok(
        self,
        precached_repo: Repo,
        initial_accessions: list[str],
        new_accessions: list[str],
        snapshot: SnapshotAssertion,
    ):
        """Test the addition of segments to an OTU plan."""
        otu_service = OTUService(precached_repo, NCBIClient(False))

        with precached_repo.lock():
            otu_before = otu_service.create(initial_accessions)

        assert otu_before is not None

        plan_before = otu_before.plan

        assert isinstance(plan_before, Plan)

        with precached_repo.lock():
            new_segment_ids = add_segments_to_plan(
                precached_repo,
                otu_before,
                rule=SegmentRule.OPTIONAL,
                accessions=new_accessions,
            )

        assert len(new_segment_ids) == len(new_accessions)

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after
        assert new_segment_ids.issubset(otu_before.plan.segment_ids)
        assert otu_after.plan.segment_ids.issubset(otu_before.plan.segment_ids)
        assert otu_after.plan.model_dump() == snapshot(exclude=props("id"))

    def test_add_segments_to_plan_fail(
        self,
        scratch_repo: Repo,
    ):
        """Test that segments cannot be added to a monopartite plan with
        a preexisting unnamed segment.
        """
        otu_before = scratch_repo.get_otu_by_taxid(96892)
        assert otu_before is not None

        assert otu_before.plan.monopartite

        with scratch_repo.lock():
            assert not add_segments_to_plan(
                scratch_repo,
                otu_before,
                rule=SegmentRule.OPTIONAL,
                accessions=["NC_010620"],
            )

    @pytest.mark.parametrize("tolerance", [0.05, 0.5, 1.0])
    def test_set_length_tolerances_ok(self, scratch_repo: Repo, tolerance: float):
        """Check that plan length tolerances can be modified by function."""
        otu_before = scratch_repo.get_otu_by_taxid(96892)
        assert otu_before is not None

        assert (
            otu_before.plan.segments[0].length_tolerance
            == scratch_repo.settings.default_segment_length_tolerance
        )

        with scratch_repo.lock():
            set_plan_length_tolerances(scratch_repo, otu_before, tolerance)

        otu_after = scratch_repo.get_otu(otu_before.id)
        assert otu_after is not None

        assert otu_after.plan.segments[0].length_tolerance == tolerance

    @pytest.mark.parametrize("bad_tolerance", [-1.0, 1.1, 100.0])
    def test_set_length_tolerances_fail(self, scratch_repo: Repo, bad_tolerance: float):
        """Check that plan length tolerances cannot be set to an invalid float value."""
        otu_before = scratch_repo.get_otu_by_taxid(96892)
        assert otu_before is not None

        assert (
            otu_before.plan.segments[0].length_tolerance
            == scratch_repo.settings.default_segment_length_tolerance
        )

        with scratch_repo.lock():
            try:
                set_plan_length_tolerances(scratch_repo, otu_before, bad_tolerance)
            except ValidationError as exc:
                for error in exc.errors():
                    assert error["type"] in ("less_than_equal", "greater_than_equal")

        otu_after = scratch_repo.get_otu(otu_before.id)
        assert otu_after is not None

        assert (
            otu_after.plan.segments[0].length_tolerance
            == otu_before.plan.segments[0].length_tolerance
        )


class TestDeleteIsolate:
    """Test isolate deletion behaviour."""

    def test_ok(self, scratch_repo):
        """Test that a given isolate can be deleted from the OTU."""
        taxid = 1169032

        otu_before = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_before is not None

        isolate_id = otu_before.get_isolate_id_by_name(
            IsolateName(type=IsolateNameType.ISOLATE, value="WMoV-6.3"),
        )

        assert type(isolate_id) is UUID

        with scratch_repo.lock():
            assert delete_isolate_from_otu(scratch_repo, otu_before, isolate_id)

        otu_after = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_after is not None

        assert isolate_id not in otu_after.isolate_ids

        assert otu_after.get_isolate(isolate_id) is None

        isolate_before = otu_before.get_isolate(isolate_id)
        assert isolate_before is not None
        assert isolate_before.accessions not in otu_after.accessions

        assert len(otu_after.isolate_ids) == len(otu_before.isolate_ids) - 1


class TestReplaceSequence:
    def test_ok(self, precached_repo: Repo):
        """Test sequence replacement and deletion."""
        otu_service = OTUService(precached_repo, NCBIClient(False))

        with precached_repo.lock():
            otu_before = otu_service.create(["MK431779"])

        assert otu_before

        sequence_before = otu_before.get_sequence_by_accession("MK431779")

        assert sequence_before

        isolate_id = next(
            iter(otu_before.get_isolate_ids_containing_sequence_id(sequence_before.id))
        )

        with precached_repo.lock():
            new_sequence = replace_sequence_in_otu(
                repo=precached_repo,
                otu=otu_before,
                new_accession="NC_003355",
                replaced_accession="MK431779",
            )

        assert new_sequence

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after

        isolate_after = otu_after.get_isolate(isolate_id)

        assert isolate_after is not None
        assert otu_after.accessions == isolate_after.accessions == {"NC_003355"}

    def test_multiple_links(self, precached_repo: Repo):
        otu_service = OTUService(precached_repo, NCBIClient(False))

        with precached_repo.lock():
            otu_before = otu_service.create(["DQ178608", "DQ178609"])

        assert otu_before
        assert otu_before == precached_repo.get_otu(otu_before.id)

        mock_isolate = IsolateFactory.build_on_plan(otu_before.plan)
        mock_sequence = mock_isolate.sequences[1]

        with precached_repo.lock(), precached_repo.use_transaction():
            sequence = precached_repo.create_sequence(
                otu_id=otu_before.id,
                accession=str(mock_sequence.accession),
                definition=mock_sequence.definition,
                segment=mock_sequence.segment,
                sequence=mock_sequence.sequence,
            )

            assert sequence

            isolate_before = precached_repo.create_isolate(
                otu_id=otu_before.id,
                name=mock_isolate.name,
            )

            assert isolate_before

            first_isolate = otu_before.isolates[0]

            precached_repo.link_sequence(
                otu_id=otu_before.id,
                isolate_id=isolate_before.id,
                sequence_id=first_isolate.sequences[0].id,
            )

            precached_repo.link_sequence(
                otu_id=otu_before.id,
                isolate_id=isolate_before.id,
                sequence_id=sequence.id,
            )

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after
        assert otu_after.accessions == {
            "DQ178608",
            "DQ178609",
            sequence.accession.key,
        }

        with precached_repo.lock():
            replace_sequence_in_otu(
                precached_repo,
                otu_after,
                new_accession="NC_038792",
                replaced_accession="DQ178608",
            )
