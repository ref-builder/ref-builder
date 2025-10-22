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
from ref_builder.repo import Repo
from ref_builder.services.otu import OTUService
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
