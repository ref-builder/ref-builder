import datetime

import pytest
from syrupy.assertion import SnapshotAssertion

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.isolate import add_genbank_isolate
from ref_builder.otu.promote import promote_otu_accessions
from ref_builder.otu.update import (
    auto_update_otu,
    batch_update_repo,
)
from ref_builder.repo import Repo
from ref_builder.services.otu import OTUService


class TestPromoteOTU:
    """Test OTU accession promotion from Genbank to RefSeq."""

    def test_ok(self, empty_repo: Repo):
        """Test that RefSeq accessions can be promoted automatically."""
        ncbi_client = NCBIClient(ignore_cache=False)
        otu_service = OTUService(empty_repo, ncbi_client)

        with empty_repo.lock():
            otu = otu_service.create(["MF062136", "MF062137", "MF062138"])

            assert otu

            isolate = add_genbank_isolate(
                empty_repo, otu, ["MF062125", "MF062126", "MF062127"], ncbi_client
            )

            assert isolate

        otu_before = empty_repo.get_otu(otu.id)

        assert otu_before
        assert otu_before.accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
            "MF062136",
            "MF062137",
            "MF062138",
        }

        isolate_before = otu_before.get_isolate(isolate.id)

        assert isolate_before
        assert isolate_before.accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
        }

        with empty_repo.lock():
            promoted_accessions = promote_otu_accessions(empty_repo, otu_before)

        assert promoted_accessions == {"NC_055390", "NC_055391", "NC_055392"}

        otu_after = empty_repo.get_otu(otu.id)

        assert otu_after
        assert otu_after.isolate_ids == otu_before.isolate_ids
        assert otu_after.accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "MF062136",
            "MF062137",
            "MF062138",
        }
        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}

        isolate_after = otu_after.get_isolate(isolate_before.id)

        assert isolate_after
        assert isolate_after.accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
        }


class TestUpdateOTU:
    """Test automatic OTU update functionality."""

    def test_ok(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test automatic update behaviour."""
        otu_service = OTUService(precached_repo, NCBIClient(False))

        with precached_repo.lock():
            otu_before = otu_service.create(["NC_055390", "NC_055391", "NC_055392"])

        assert otu_before
        assert otu_before.accessions == {"NC_055390", "NC_055391", "NC_055392"}
        assert otu_before.blocked_accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "MF062125",
            "MF062126",
            "MF062127",
        }

        with precached_repo.lock():
            auto_update_otu(precached_repo, otu_before)

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after
        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}
        assert otu_after.id == otu_before.id
        assert otu_after.isolate_ids.issuperset(otu_before.isolate_ids)
        assert otu_after.accessions == {
            "MF062130",
            "MF062131",
            "MF062132",
            "MF062136",
            "MF062137",
            "MF062138",
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "OR889795",
            "OR889796",
            "OR889797",
        }
        assert {
            str(isolate.name): isolate.accessions for isolate in otu_after.isolates
        } == snapshot()

    def test_start_date_limit(self, precached_repo: Repo):
        """Test automatic update with the start date set to ``today``."""
        otu_service = OTUService(precached_repo, NCBIClient(False))

        with precached_repo.lock():
            otu_before = otu_service.create(["NC_055390", "NC_055391", "NC_055392"])

        assert otu_before
        assert otu_before.accessions == {"NC_055390", "NC_055391", "NC_055392"}

        with precached_repo.lock():
            otu_after = auto_update_otu(
                precached_repo,
                otu_before,
                start_date=datetime.date.today(),
            )

        assert otu_after
        assert {
            "MF062130",
            "MF062131",
            "MF062132",
            "MF062136",
            "MF062137",
            "MF062138",
            "OR889795",
            "OR889796",
            "OR889797",
        }.isdisjoint(otu_after.accessions)

    def test_with_refseq_replacement_ok(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test that automatic update replaces accessions superceded by RefSeq."""
        otu_service = OTUService(precached_repo, NCBIClient(False))

        with precached_repo.lock():
            otu_before = otu_service.create(["MF062125", "MF062126", "MF062127"])

        assert otu_before

        isolate_before = otu_before.isolates[0]

        assert isolate_before
        assert (
            otu_before.accessions
            == isolate_before.accessions
            == otu_before.get_isolate(isolate_before.id).accessions
            == {"MF062125", "MF062126", "MF062127"}
        )

        with precached_repo.lock():
            otu_after = auto_update_otu(precached_repo, otu_before)

        isolate_after = otu_after.get_isolate(isolate_before.id)

        assert isolate_after
        assert isolate_after.accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
        }
        assert {"MF062125", "MF062126", "MF062127"}.isdisjoint(otu_after.accessions)
        assert isolate_before.accessions != isolate_after.accessions
        assert otu_after.id == otu_before.id
        assert otu_after.isolate_ids.issuperset(otu_before.isolate_ids)
        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}
        assert {
            str(isolate.name): isolate.accessions for isolate in otu_after.isolates
        } == snapshot()


class TestBatchUpdate:
    """Test rudimentary batch update operation with a single OTU."""

    @pytest.fixture(autouse=True)
    def setup(self, precached_repo: Repo):
        """Set up mock repo for batch update tests."""
        otu_service = OTUService(precached_repo, NCBIClient(False))

        with precached_repo.lock():
            otu = otu_service.create(["MF062125", "MF062126", "MF062127"])

        self.repo = precached_repo
        self.expected_accessions = {
            "MF062130",
            "MF062131",
            "MF062132",
            "MF062136",
            "MF062137",
            "MF062138",
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "OR889795",
            "OR889796",
            "OR889797",
        }


    def test_ok(self):
        """Test that batch update works as expected."""
        with self.repo.lock():
            assert len(batch_update_repo(self.repo)) == 1

        otu_after = next(self.repo.iter_otus())

        assert otu_after.accessions == self.expected_accessions

        with self.repo.lock():
            assert len(batch_update_repo(self.repo)) == 0
