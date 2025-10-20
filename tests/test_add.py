import pytest
from pytest_mock import MockerFixture
from structlog.testing import capture_logs

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import (
    NCBIGenbank,
    NCBIRank,
    NCBITaxonomy,
    NCBITaxonomyOtherNames,
)
from ref_builder.otu.isolate import (
    add_and_name_isolate,
    add_genbank_isolate,
    add_unnamed_isolate,
)
from ref_builder.repo import Repo
from ref_builder.services.otu import OTUService
from ref_builder.utils import IsolateName, IsolateNameType
from tests.fixtures.factories import (
    NCBIGenbankFactory,
    NCBISourceFactory,
    NCBITaxonomyFactory,
)


def create_mocked_otu_service(
    repo: Repo,
    mocker: MockerFixture,
    taxonomy: NCBITaxonomy,
    records: list[NCBIGenbank],
) -> OTUService:
    """Create an OTUService with a mocked NCBIClient.

    :param repo: repository to use
    :param mocker: pytest-mock fixture
    :param taxonomy: taxonomy record to return from fetch_taxonomy_record
    :param records: genbank records to return from fetch_genbank_records
    :return: OTUService with mocked NCBIClient
    """
    ncbi_client = mocker.create_autospec(NCBIClient, instance=True)
    ncbi_client.fetch_genbank_records = mocker.Mock(return_value=records)
    ncbi_client.fetch_taxonomy_record = mocker.Mock(return_value=taxonomy)
    return OTUService(repo, ncbi_client)


class TestCreateOTU:
    @pytest.mark.parametrize(
        ("taxid", "segment_count"),
        [
            (3426695, 2),
            (3240630, 1),
        ],
    )
    def test_ok(
        self,
        empty_repo: Repo,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
        taxid: int,
        segment_count: int,
    ):
        """Test the the default, happy case works."""
        taxonomy = ncbi_taxonomy_factory.build(id=taxid, rank=NCBIRank.SPECIES)

        if segment_count == 1:
            records = [
                ncbi_genbank_factory.build(
                    source__taxid=taxid,
                    source__isolate="Isolate A",
                    accession="AB123456",
                )
            ]
        else:
            records = ncbi_genbank_factory.build_isolate(
                segment_count=segment_count,
                refseq=False,
                base_source=NCBISourceFactory.build(taxid=taxid, isolate="Isolate A"),
            )

        otu_service = create_mocked_otu_service(empty_repo, mocker, taxonomy, records)
        accessions = [r.accession for r in records]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu
        assert otu.taxid == taxid
        assert empty_repo.get_otu_by_taxid(taxid) == otu

    def test_no_accessions(self, empty_repo: Repo):
        """Test that creating an OTU with not accession fails with the expected log."""
        otu_service = OTUService(empty_repo, NCBIClient(False))

        with empty_repo.lock(), capture_logs() as logs:
            otu = otu_service.create([])

            assert otu is None
            assert any(
                "OTU could not be created to spec" in (log.get("event") or "")
                for log in logs
            )

    def test_duplicate_taxid(
        self,
        empty_repo: Repo,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that an OTU with the same taxid cannot be created."""
        taxid = 345184
        taxonomy = ncbi_taxonomy_factory.build(id=taxid, rank=NCBIRank.SPECIES)
        records = ncbi_genbank_factory.build_isolate(
            segment_count=2,
            refseq=False,
            base_source=NCBISourceFactory.build(taxid=taxid, isolate="Isolate A"),
        )

        otu_service = create_mocked_otu_service(empty_repo, mocker, taxonomy, records)
        accessions = [r.accession for r in records]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu
        assert otu.taxid == taxid

        with empty_repo.lock():
            duplicate_otu = otu_service.create(accessions)
            assert duplicate_otu is None

    def test_refseq_autoexclude(
        self,
        empty_repo: Repo,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that the superceded accessions included in RefSeq metadata are
        automatically added to the OTU's excluded accessions list.
        """
        taxid = 3158377
        taxonomy = ncbi_taxonomy_factory.build(id=taxid, rank=NCBIRank.SPECIES)

        # Create 6 RefSeq records with comments pointing to old accessions
        old_accessions = [f"EF546{800 + i}" for i in range(8, 14)]

        records = ncbi_genbank_factory.build_isolate(
            segment_count=6,
            refseq=True,
            base_source=NCBISourceFactory.build(taxid=taxid, isolate="Isolate A"),
        )

        # Update the records with proper RefSeq accessions and comments
        for i, (record, old_acc) in enumerate(
            zip(records, old_accessions, strict=False), start=4
        ):
            record.accession = f"NC_01031{i}"
            record.comment = (
                f"PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review. "
                f"The reference sequence is identical to {old_acc}."
            )

        otu_service = create_mocked_otu_service(empty_repo, mocker, taxonomy, records)
        accessions = [r.accession for r in records]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu
        assert otu.excluded_accessions == set(old_accessions)

    def test_acronym(
        self,
        empty_repo: Repo,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ) -> None:
        """Test that the acronym pulled from NCBI is correct."""
        other_names = NCBITaxonomyOtherNames()
        other_names.acronym = ["KLV"]

        taxonomy = NCBITaxonomy(
            id=132477,
            name="Kasba virus",
            rank=NCBIRank.SPECIES,
            other_names=other_names,
            lineage=[],
        )

        records = [
            ncbi_genbank_factory.build(source__taxid=taxonomy.id, accession="AB123456")
        ]

        otu_service = create_mocked_otu_service(empty_repo, mocker, taxonomy, records)

        with empty_repo.lock():
            otu = otu_service.create(["AB123456"])

        assert otu is not None
        assert otu.acronym == "KLV"


class TestAddIsolate:
    def test_multipartite(self, precached_repo: Repo):
        ncbi_client = NCBIClient(ignore_cache=True)
        otu_service = OTUService(precached_repo, ncbi_client)

        with precached_repo.lock():
            otu_before = otu_service.create(["MF062136", "MF062137", "MF062138"])

        assert otu_before
        assert otu_before.accessions == {"MF062136", "MF062137", "MF062138"}
        assert len(otu_before.isolate_ids) == 1

        with precached_repo.lock():
            isolate = add_genbank_isolate(
                precached_repo,
                otu_before,
                ["MF062125", "MF062126", "MF062127"],
                ncbi_client,
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

        ncbi_client = NCBIClient(ignore_cache=True)
        otu_service = OTUService(precached_repo, ncbi_client)

        with precached_repo.lock():
            otu_before = otu_service.create(isolate_1_accessions)

            assert otu_before
            assert otu_before.accessions == set(isolate_1_accessions)

            isolate = add_genbank_isolate(
                precached_repo, otu_before, isolate_2_accessions, ncbi_client
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

        ncbi_client = NCBIClient(ignore_cache=True)
        otu_service = OTUService(precached_repo, ncbi_client)

        with precached_repo.lock():
            otu_before = otu_service.create(isolate_1_accessions)

        assert otu_before
        assert otu_before.accessions == set(isolate_1_accessions)

        isolate_1 = otu_before.isolates[0]

        assert isolate_1
        assert isolate_1.accessions == set(isolate_1_accessions)

        with precached_repo.lock():
            isolate_2 = add_unnamed_isolate(
                precached_repo, otu_before, isolate_2_accessions, ncbi_client
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

        ncbi_client = NCBIClient(ignore_cache=True)
        otu_service = OTUService(precached_repo, ncbi_client)

        with precached_repo.lock():
            otu = otu_service.create(isolate_1_accessions)

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
                IsolateName(type=IsolateNameType.ISOLATE, value="dummy"),
                ncbi_client,
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

        ncbi_client = NCBIClient(ignore_cache=True)
        otu_service = OTUService(precached_repo, ncbi_client)

        with precached_repo.lock():
            otu = otu_service.create(accessions)

            assert otu
            assert (
                add_genbank_isolate(precached_repo, otu, accessions, ncbi_client)
                is None
            )

    def test_name_conflict(self, precached_repo: Repo):
        """Test that an isolate cannot be added to an OTU if its name is already contained."""
        ncbi_client = NCBIClient(ignore_cache=True)
        otu_service = OTUService(precached_repo, ncbi_client)

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
                ncbi_client,
            )

        assert isolate

        otu_after = precached_repo.get_otu(otu_before.id)

        assert otu_after
        assert isolate.id not in otu_after.isolate_ids
        assert precached_repo.last_id == last_event_id_before_test

    def test_plan_mismatch(self, scratch_repo: Repo, mock_ncbi_client):
        """Test that an isolate cannot be added to if it does not match the plan."""
        taxid = mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        otu_before = scratch_repo.get_otu_by_taxid(taxid)

        assert otu_before

        ncbi_client = NCBIClient(ignore_cache=True)

        with scratch_repo.lock():
            assert (
                add_genbank_isolate(
                    scratch_repo,
                    otu_before,
                    ["AB017503"],
                    ncbi_client,
                )
                is None
            )

        otu_after = scratch_repo.get_otu_by_taxid(taxid)

        assert otu_after
        assert otu_after == otu_before
        assert "AB017503" not in otu_after.accessions

    def test_auto_promote_refseq_monopartite(self, empty_repo: Repo):
        """Test that RefSeq accessions automatically promote GenBank sequences (monopartite)."""
        taxid = 1169032
        original_accessions = ["AB017503"]
        refseq_accessions = ["NC_003355"]

        ncbi_client = NCBIClient(ignore_cache=True)
        otu_service = OTUService(empty_repo, ncbi_client)

        with empty_repo.lock():
            otu = otu_service.create(accessions=original_accessions)

        assert otu
        assert otu.accessions == set(original_accessions)

        original_isolate = otu.isolates[0]
        original_isolate_name = original_isolate.name

        with empty_repo.lock():
            isolate = add_genbank_isolate(
                empty_repo, otu, refseq_accessions, ncbi_client
            )

        assert isolate is not None
        assert isolate.name == original_isolate_name
        assert isolate.accessions == set(refseq_accessions)

        otu_after = empty_repo.get_otu(otu.id)

        assert otu_after
        assert otu_after.accessions == set(refseq_accessions)
        assert set(original_accessions).issubset(otu_after.excluded_accessions)

    def test_auto_promote_refseq_multipartite(self, empty_repo: Repo, mock_ncbi_client):
        """Test that RefSeq accessions automatically promote GenBank sequences (multipartite)."""
        taxid = mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        original_accessions = ["DQ178608", "DQ178609"]
        refseq_accessions = ["NC_038792", "NC_038793"]

        ncbi_client = NCBIClient(ignore_cache=True)
        otu_service = OTUService(empty_repo, ncbi_client)

        with empty_repo.lock():
            otu = otu_service.create(accessions=original_accessions)

        assert otu
        assert otu.accessions == set(original_accessions)

        original_isolate = otu.isolates[0]
        original_isolate_name = original_isolate.name

        with empty_repo.lock():
            isolate = add_genbank_isolate(
                empty_repo, otu, refseq_accessions, ncbi_client
            )

        assert isolate is not None
        assert isolate.name == original_isolate_name
        assert isolate.accessions == set(refseq_accessions)

        otu_after = empty_repo.get_otu(otu.id)

        assert otu_after
        assert otu_after.accessions == set(refseq_accessions)
        assert set(original_accessions).issubset(otu_after.excluded_accessions)
