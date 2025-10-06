import pytest
from structlog.testing import capture_logs
from syrupy.assertion import SnapshotAssertion
from syrupy.filters import props

from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.create import (
    create_otu_with_taxid,
    create_otu_without_taxid,
)
from ref_builder.otu.isolate import (
    add_and_name_isolate,
    add_genbank_isolate,
    add_unnamed_isolate,
)
from ref_builder.repo import Repo
from ref_builder.services.otu import OTUService
from ref_builder.utils import IsolateName, IsolateNameType


class TestCreateOTU:
    @pytest.mark.parametrize(
        ("accessions", "expected_taxid"),
        [
            (["DQ178610", "DQ178611"], 3426695),
            (["NC_043170"], 3240630),
        ],
    )
    def test_ok(self, accessions: list[str], expected_taxid: int, precached_repo: Repo):
        with precached_repo.lock():
            otu = create_otu_without_taxid(
                precached_repo, accessions=accessions, acronym=""
            )

        assert otu
        assert otu.taxid == expected_taxid
        assert precached_repo.get_otu_by_taxid(expected_taxid) == otu

    def test_empty_repo(
        self,
        precached_repo: Repo,
        snapshot: SnapshotAssertion,
    ):
        """Test that an OTU can be created in an empty repository."""
        with precached_repo.lock():
            otu = create_otu_with_taxid(
                precached_repo,
                345184,
                ["DQ178610", "DQ178611"],
                "",
            )

        assert otu
        assert otu.model_dump() == snapshot(
            exclude=props("id", "isolates"),
        )
        assert list(precached_repo.iter_otus()) == [otu]

    def test_no_accessions(self, empty_repo: Repo):
        """Test that creating an OTU with not accession fails with the expected log."""
        with empty_repo.lock(), capture_logs() as logs:
            otu = create_otu_with_taxid(
                empty_repo,
                345184,
                [],
                "",
            )

            assert otu is None
            assert any(
                "OTU could not be created to spec" in (log.get("event") or "")
                for log in logs
            )

    def test_duplicate_taxid(self, precached_repo: Repo):
        """Test that an OTU with the same taxid cannot be created."""
        accessions = ["DQ178610", "DQ178611"]
        taxid = 345184

        with precached_repo.lock():
            otu = create_otu_with_taxid(
                precached_repo,
                taxid,
                accessions,
                "",
            )

        assert otu

        with (
            pytest.raises(
                ValueError,
                match="Taxonomy ID 345184 has already been added to this reference.",
            ),
            precached_repo.lock(),
        ):
            create_otu_with_taxid(
                precached_repo,
                taxid,
                accessions,
                "",
            )

    def test_refseq_autoexclude(self, precached_repo: Repo):
        """Test that the superceded accessions included in RefSeq metadata are
        automatically added to the OTU's excluded accessions list.
        """
        with precached_repo.lock():
            otu = create_otu_with_taxid(
                precached_repo,
                3158377,
                [
                    "NC_010314",
                    "NC_010316",
                    "NC_010315",
                    "NC_010317",
                    "NC_010318",
                    "NC_010319",
                ],
                "",
            )

        assert otu
        assert otu.excluded_accessions == {
            "EF546808",
            "EF546809",
            "EF546810",
            "EF546811",
            "EF546812",
            "EF546813",
        }

    def test_acronym(self, precached_repo: Repo) -> None:
        """Test that the acronym pulled from NCBI is correct."""
        with precached_repo.lock():
            otu = create_otu_with_taxid(
                precached_repo,
                132477,
                ["NC_013006"],
                "",
            )

        assert otu is not None
        assert otu.acronym == "KLV"

    def test_acronym_manual(self, precached_repo: Repo) -> None:
        """Test the acronym can be set manually."""
        with precached_repo.lock():
            otu = create_otu_with_taxid(
                precached_repo,
                1441799,
                ["NC_023881"],
                "FBNSV",
            )

        assert otu
        assert otu.acronym == "FBNSV"


class TestAddIsolate:
    def test_multipartite(self, precached_repo: Repo):
        otu_service = OTUService(precached_repo, NCBIClient(True))

        with precached_repo.lock():
            otu_before = otu_service.create(
                ["MF062136", "MF062137", "MF062138"],
                acronym="",
            )

        assert otu_before
        assert otu_before.accessions == {"MF062136", "MF062137", "MF062138"}
        assert len(otu_before.isolate_ids) == 1

        with precached_repo.lock():
            isolate = add_genbank_isolate(
                precached_repo, otu_before, ["MF062125", "MF062126", "MF062127"]
            )

        assert isolate

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after
        assert otu_after.accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
            "MF062136",
            "MF062137",
            "MF062138",
        }

        isolate_after = otu_after.get_isolate(isolate.id)

        assert isolate_after
        assert isolate_after.accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
        }

    def test_genbank(self, precached_repo: Repo):
        """Test that add_genbank_isolate() adds an isolate with a correctly parsed
        name.
        """
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        otu_service = OTUService(precached_repo, NCBIClient(True))

        with precached_repo.lock():
            otu_before = otu_service.create(isolate_1_accessions, acronym="")

            assert otu_before
            assert otu_before.accessions == set(isolate_1_accessions)

            isolate = add_genbank_isolate(
                precached_repo, otu_before, isolate_2_accessions
            )

            assert isolate

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after
        assert otu_after.accessions == set(isolate_1_accessions).union(
            set(isolate_2_accessions),
        )

        assert isolate == otu_after.get_isolate(isolate.id)
        assert isolate.accessions == set(isolate_2_accessions)
        assert isolate.name == IsolateName(
            IsolateNameType.ISOLATE,
            "Douglas Castle",
        )

    def test_ignore_name(self, precached_repo: Repo):
        """Test that add_unnamed_isolate() adds the isolate with a ``None`` name."""
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        otu_service = OTUService(precached_repo, NCBIClient(True))

        with precached_repo.lock():
            otu_before = otu_service.create(isolate_1_accessions, acronym="")

        assert otu_before
        assert otu_before.accessions == set(isolate_1_accessions)

        isolate_1 = otu_before.isolates[0]

        assert isolate_1
        assert isolate_1.accessions == set(isolate_1_accessions)

        with precached_repo.lock():
            isolate_2 = add_unnamed_isolate(
                precached_repo, otu_before, isolate_2_accessions
            )

        assert isolate_2

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after
        assert otu_after.isolate_ids == {isolate_1.id, isolate_2.id}

        assert isolate_2 == otu_after.get_isolate(isolate_2.id)
        assert isolate_2.name is None
        assert isolate_2.accessions == set(isolate_2_accessions)

    def test_add_and_name_isolate(self, precached_repo: Repo):
        """Test that add_and_name_isolate() creates an isolate with the correct name."""
        isolate_1_accessions = ["DQ178610", "DQ178611"]
        isolate_2_accessions = ["DQ178613", "DQ178614"]

        otu_service = OTUService(precached_repo, NCBIClient(True))

        with precached_repo.lock():
            otu = otu_service.create(isolate_1_accessions, acronym="")

        assert otu
        assert otu.accessions == set(isolate_1_accessions)

        isolate_1 = otu.isolates[0]

        assert isolate_1
        assert isolate_1 == otu.get_isolate(isolate_1.id)

        with precached_repo.lock():
            isolate_2 = add_and_name_isolate(
                precached_repo,
                otu,
                isolate_2_accessions,
                isolate_name=IsolateName(type=IsolateNameType.ISOLATE, value="dummy"),
            )

        assert isolate_2

        otu_after = precached_repo.get_otu(otu.id)

        assert otu_after
        assert otu_after.isolate_ids == {isolate_1.id, isolate_2.id}

        assert isolate_2 == otu_after.get_isolate(isolate_2.id)
        assert isolate_2.accessions == set(isolate_2_accessions)
        assert isolate_2.name == IsolateName(
            type=IsolateNameType.ISOLATE,
            value="dummy",
        )

    def test_blocked(self, precached_repo: Repo):
        """Test that an isolate cannot be added to an OTU if both its name and its
        accessions are already contained.
        """
        accessions = ["MF062136", "MF062137", "MF062138"]

        otu_service = OTUService(precached_repo, NCBIClient(True))

        with precached_repo.lock():
            otu = otu_service.create(accessions, "")

            assert otu
            assert add_genbank_isolate(precached_repo, otu, accessions) is None

    def test_name_conflict(self, precached_repo: Repo):
        """Test that an isolate cannot be added to an OTU if its name is already contained."""
        otu_service = OTUService(precached_repo, NCBIClient(True))

        with precached_repo.lock():
            otu_before = otu_service.create(["MF062136", "MF062137", "MF062138"])

        assert otu_before

        last_event_id_before_test = precached_repo.last_id

        with precached_repo.lock():
            isolate = add_and_name_isolate(
                precached_repo,
                otu_before,
                ["NC_055390", "NC_055391", "NC_055392"],
                IsolateName(type=IsolateNameType.ISOLATE, value="4342-5"),
            )

        assert isolate

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after
        assert isolate.id not in otu_after.isolate_ids
        assert precached_repo.last_id == last_event_id_before_test

    def test_plan_mismatch(self, scratch_repo: Repo):
        """Test that an isolate cannot be added to if it does not match the plan."""
        otu_before = scratch_repo.get_otu_by_taxid(345184)

        assert otu_before

        with scratch_repo.lock():
            assert (
                add_genbank_isolate(
                    scratch_repo,
                    otu_before,
                    accessions=["AB017503"],
                )
                is None
            )

        otu_after = scratch_repo.get_otu_by_taxid(345184)

        assert otu_after
        assert otu_after == otu_before
        assert "AB017503" not in otu_after.accessions

    @pytest.mark.parametrize(
        "taxid, original_accessions, refseq_accessions",
        [
            (1169032, ["AB017503"], ["NC_003355"]),
            (345184, ["DQ178608", "DQ178609"], ["NC_038792", "NC_038793"]),
        ],
    )
    def test_auto_promote_refseq(
        self,
        empty_repo: Repo,
        taxid: int,
        original_accessions: list[str],
        refseq_accessions: list[str],
    ):
        """Test that RefSeq accessions automatically promote GenBank sequences."""
        otu_service = OTUService(empty_repo, NCBIClient(True))

        with empty_repo.lock():
            otu = otu_service.create(accessions=original_accessions, acronym="")

        assert otu
        assert otu.accessions == set(original_accessions)

        original_isolate = otu.isolates[0]
        original_isolate_name = original_isolate.name

        with empty_repo.lock():
            isolate = add_genbank_isolate(empty_repo, otu, refseq_accessions)

        assert isolate is not None
        assert isolate.name == original_isolate_name
        assert isolate.accessions == set(refseq_accessions)

        otu_after = empty_repo.get_otu(otu.id)

        assert otu_after
        assert otu_after.accessions == set(refseq_accessions)
        assert set(original_accessions).issubset(otu_after.excluded_accessions)
