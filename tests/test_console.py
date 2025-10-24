import pytest
from faker import Faker
from syrupy import SnapshotAssertion

from ref_builder.console import (
    console,
    print_isolate,
    print_isolate_as_json,
    print_otu,
    print_otu_list,
)
from ref_builder.models.otu import OTUMinimal
from ref_builder.models.plan import Plan, Segment, SegmentName, SegmentRule
from ref_builder.otu.builders.isolate import IsolateBuilder
from tests.fixtures.factories import IsolateFactory, OTUMinimalFactory
from tests.fixtures.providers import AccessionProvider, SequenceProvider


class TestPrintOTUList:
    """Tests for the ``print_otu_list`` function."""

    def test_ok(
        self,
        capsys: pytest.CaptureFixture,
        otu_minimal_factory: OTUMinimalFactory,
        snapshot: SnapshotAssertion,
    ) -> None:
        """Test that listed OTUs are printed."""
        print_otu_list(otu_minimal_factory.build() for _ in range(5))

        assert capsys.readouterr().out == snapshot

    def test_empty(self) -> None:
        """Test that an empty list of OTUs is printed."""
        with console.capture() as capture:
            print_otu_list(OTUMinimal(**otu) for otu in [])

        assert capture.get() == "No OTUs found\n"


class TestPrintIsolate:
    """Test isolate console output."""

    def test_ok(self, snapshot: SnapshotAssertion):
        """Test that an isolate is printed as expected by ``print_isolate``."""
        fake = Faker(["en_US"])
        fake.add_provider(AccessionProvider)
        fake.add_provider(SequenceProvider)
        fake.seed_instance(8801)

        plan = Plan.new(
            [
                Segment(
                    id=fake.uuid4(),
                    length=1099,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "R"),
                    rule=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1074,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "M"),
                    rule=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1087,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "S"),
                    rule=SegmentRule.REQUIRED,
                ),
            ],
        )

        isolate = IsolateBuilder(**IsolateFactory.build_on_plan(plan).model_dump())

        with console.capture() as capture:
            print_isolate(isolate, plan)

        assert capture.get() == snapshot

    def test_json_ok(self, snapshot: SnapshotAssertion):
        """Test that an isolate is printed as expected by ``print_isolate_as_json``."""
        fake = Faker(["en_US"])
        fake.add_provider(AccessionProvider)
        fake.add_provider(SequenceProvider)
        fake.seed_instance(8801)

        plan = Plan.new(
            [
                Segment(
                    id=fake.uuid4(),
                    length=1099,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "R"),
                    rule=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1074,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "M"),
                    rule=SegmentRule.REQUIRED,
                ),
                Segment(
                    id=fake.uuid4(),
                    length=1087,
                    length_tolerance=0.03,
                    name=SegmentName("DNA", "S"),
                    rule=SegmentRule.REQUIRED,
                ),
            ],
        )

        isolate = IsolateBuilder(**IsolateFactory.build_on_plan(plan).model_dump())

        with console.capture() as capture:
            print_isolate_as_json(isolate)

        assert capture.get() == snapshot


class TestPrintOTU:
    """Test OTU console output."""

    def test_print_otu(self, scratch_repo):
        """Test that an OTU is printed as expected by ``print_otu``."""
        otu = scratch_repo.get_otu_by_taxid(3429802)

        with console.capture() as capture:
            print_otu(otu)

        output = capture.get()
        assert "Hostuviroid latensdahliae" in output
        assert "TAXID" in output
        assert "3429802" in output
        assert "NC_020160.1" in output
