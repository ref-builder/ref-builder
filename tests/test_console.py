import pytest

from ref_builder.console import (
    console,
    print_isolate,
    print_otu,
    print_otu_list,
)
from ref_builder.models.otu import OTUMinimal
from ref_builder.ncbi.client import NCBIClientProtocol
from ref_builder.repo import Repo
from tests.fixtures.mock_ncbi_client import MockNCBIClient


class TestPrintOTUList:
    """Tests for the ``print_otu_list`` function."""

    def test_ok(
        self,
        capsys: pytest.CaptureFixture,
        scratch_repo: Repo,
    ) -> None:
        """Test that listed OTUs are printed."""
        otus = scratch_repo.iter_minimal_otus()
        print_otu_list(otus)

        output = capsys.readouterr().out

        assert "NAME" in output
        assert "ACRONYM" in output
        assert "TAXID" in output
        assert "ID" in output

        for otu in otus:
            assert otu.name in output
            assert str(otu.taxid) in output

    def test_empty(self) -> None:
        """Test that an empty list of OTUs is printed."""
        with console.capture() as capture:
            print_otu_list(OTUMinimal(**otu) for otu in [])

        assert capture.get() == "No OTUs found\n"


class TestPrintIsolate:
    """Test isolate console output."""

    def test_ok(
        self,
        mock_ncbi_client: MockNCBIClient,
        scratch_repo: Repo,
    ):
        """Test that an isolate is printed as expected by ``print_isolate``."""
        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.abaca_bunchy_top_virus.taxid
        )

        assert otu

        isolate = otu.isolates[0]

        with console.capture() as capture:
            print_isolate(isolate, otu.plan)

        output = capture.get()

        assert str(isolate.name) in output

        assert "ACCESSION" in output
        assert "LENGTH" in output
        assert "SEGMENT" in output
        assert "DEFINITION" in output

        for sequence in isolate.sequences:
            assert str(sequence.accession) in output
            assert str(len(sequence.sequence)) in output


class TestPrintOTU:
    """Test OTU console output."""

    def test_ok(
        self,
        mock_ncbi_client: NCBIClientProtocol,
        scratch_repo: Repo,
    ):
        """Test that an OTU is printed as expected by ``print_otu``."""
        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.abaca_bunchy_top_virus.taxid
        )

        assert otu

        with console.capture() as capture:
            print_otu(otu)

        output = capture.get()

        assert otu.name in output

        # Check metadata section
        assert "ACRONYM" in output
        assert "ID" in output
        assert "TAXID" in output
        assert str(otu.taxid) in output

        # Check plan section
        assert "PLAN" in output
        assert "NAME" in output
        assert "REQUIRED" in output
        assert "LENGTH" in output
        assert "TOLERANCE" in output

        # Check all segments in plan
        for segment in otu.plan.segments:
            assert str(segment.name) in output
            assert str(segment.length) in output

        # Check isolates section
        assert "ISOLATES" in output

        # Check all isolates are present
        for isolate in otu.isolates:
            assert str(isolate.name) in output

            # Check sequences from each isolate
            for sequence in isolate.sequences:
                assert str(sequence.accession) in output
