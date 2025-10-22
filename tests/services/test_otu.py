import pytest
from pytest_mock import MockerFixture

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIRank, NCBITaxonomy, NCBITaxonomyOtherNames
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.repo import Repo
from ref_builder.services.otu import OTUService
from tests.fixtures.factories import (
    NCBIGenbankFactory,
    NCBISourceFactory,
    NCBITaxonomyFactory,
)


@pytest.fixture
def otu_service(empty_repo: Repo, mocker: MockerFixture) -> OTUService:
    """Create an OTUService with a mocked NCBI client."""
    ncbi_client = mocker.create_autospec(NCBIClient, instance=True)
    return OTUService(empty_repo, ncbi_client)


class TestOTUServiceCreate:
    """Test successful OTU creation scenarios."""

    def test_ok(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test basic OTU creation with monopartite genome."""
        taxonomy = ncbi_taxonomy_factory.build(
            id=345184,
            rank=NCBIRank.SPECIES,
        )
        records = [
            ncbi_genbank_factory.build(source__taxid=taxonomy.id, accession="AB123456")
        ]

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=records)
        otu_service.ncbi.fetch_taxonomy_record = mocker.Mock(return_value=taxonomy)

        accessions = [records[0].accession]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu is not None
        assert isinstance(otu, OTUBuilder)
        assert otu.taxid == taxonomy.id
        assert otu.name == taxonomy.name

        otus = list(empty_repo.iter_otus())
        assert len(otus) == 1
        assert otus[0].id == otu.id

    def test_ok_multipartite(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test OTU creation with multipartite genome."""
        taxonomy = ncbi_taxonomy_factory.build(
            id=223262,
            rank=NCBIRank.SPECIES,
        )
        records = ncbi_genbank_factory.build_isolate(
            segment_count=3,
            refseq=False,
            base_source__taxid=taxonomy.id,
        )

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=records)
        otu_service.ncbi.fetch_taxonomy_record = mocker.Mock(return_value=taxonomy)

        accessions = [r.accession for r in records]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu is not None
        assert len(otu.plan.segments) == 3
        assert len(otu.isolates) == 1
        assert len(otu.isolates[0].sequences) == 3

    def test_acronym_from_ncbi(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test auto-fill acronym from NCBI taxonomy."""
        other_names = NCBITaxonomyOtherNames()
        other_names.acronym = ["TestAcronym"]

        taxonomy = NCBITaxonomy(
            id=12345,
            name="Test Virus",
            rank=NCBIRank.SPECIES,
            other_names=other_names,
            lineage=[],
        )

        records = [
            ncbi_genbank_factory.build(source__taxid=taxonomy.id, accession="AB123456")
        ]

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=records)
        otu_service.ncbi.fetch_taxonomy_record = mocker.Mock(return_value=taxonomy)

        accessions = [records[0].accession]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu is not None
        assert otu.acronym == "TestAcronym"

    def test_refseq_exclusion(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test RefSeq comment parsing and old accession exclusion."""
        taxonomy = ncbi_taxonomy_factory.build(rank=NCBIRank.SPECIES)
        record = ncbi_genbank_factory.build(
            source__taxid=taxonomy.id,
            accession="NC_123456",
            comment="PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review. "
            "The reference sequence is identical to DQ178610,COMPLETENESS: full length.",
        )

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=[record])
        otu_service.ncbi.fetch_taxonomy_record = mocker.Mock(return_value=taxonomy)

        accessions = [record.accession]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu is not None
        assert "DQ178610" in otu.excluded_accessions


class TestOTUServiceCreateValidation:
    """Test input validation failures."""

    def test_empty_accessions(self, otu_service: OTUService):
        """Test that empty accessions list returns None."""
        otu = otu_service.create([])

        assert otu is None

    def test_missing_accessions(
        self,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that missing accessions returns None."""
        records = [ncbi_genbank_factory.build()]

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=records)

        accessions = [records[0].accession, "MISSING123"]
        otu = otu_service.create(accessions)

        assert otu is None

    def test_multiple_organisms(
        self,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that records from different organisms returns None."""
        record_1 = ncbi_genbank_factory.build(source__taxid=12345)
        record_2 = ncbi_genbank_factory.build(source__taxid=67890)

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(
            return_value=[record_1, record_2]
        )

        accessions = [record_1.accession, record_2.accession]
        otu = otu_service.create(accessions)

        assert otu is None

    def test_multiple_isolates(
        self,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that records from multiple isolates returns None."""
        taxid = 12345
        record_1 = ncbi_genbank_factory.build(
            source__taxid=taxid,
            source__isolate="isolate1",
        )
        record_2 = ncbi_genbank_factory.build(
            source__taxid=taxid,
            source__isolate="isolate2",
        )

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(
            return_value=[record_1, record_2]
        )

        accessions = [record_1.accession, record_2.accession]
        otu = otu_service.create(accessions)

        assert otu is None


class TestOTUServiceCreateTaxonomy:
    """Test taxonomy-related edge cases."""

    def test_strain_level_enforcement(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test auto-promotion of strain-level taxid to species level."""
        strain_taxonomy = ncbi_taxonomy_factory.build(
            id=999999,
            rank=NCBIRank.ISOLATE,
        )

        species_id = strain_taxonomy.lineage[-1].id
        species_taxonomy = ncbi_taxonomy_factory.build(
            id=species_id,
            rank=NCBIRank.SPECIES,
        )

        source = NCBISourceFactory.build(taxid=strain_taxonomy.id)
        records = [ncbi_genbank_factory.build(source=source, accession="AB123456")]

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=records)
        otu_service.ncbi.fetch_taxonomy_record = mocker.Mock(
            side_effect=[strain_taxonomy, species_taxonomy]
        )

        accessions = [records[0].accession]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu is not None
        assert otu.taxid == species_id

        assert otu_service.ncbi.fetch_taxonomy_record.call_count == 2
        otu_service.ncbi.fetch_taxonomy_record.assert_any_call(strain_taxonomy.id)
        otu_service.ncbi.fetch_taxonomy_record.assert_any_call(species_id)

    def test_duplicate_taxid(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that duplicate taxid returns None."""
        taxonomy = ncbi_taxonomy_factory.build(rank=NCBIRank.SPECIES)
        records = [
            ncbi_genbank_factory.build(source__taxid=taxonomy.id, accession="AB123456")
        ]

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=records)
        otu_service.ncbi.fetch_taxonomy_record = mocker.Mock(return_value=taxonomy)

        with empty_repo.lock():
            first_otu = otu_service.create([records[0].accession])

        assert first_otu is not None

        record_2 = ncbi_genbank_factory.build(
            source__taxid=taxonomy.id, accession="AB789012"
        )
        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=[record_2])

        with empty_repo.lock():
            second_otu = otu_service.create([record_2.accession])

        assert second_otu is None

    def test_taxonomy_not_found(
        self,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that taxonomy lookup failure returns None."""
        records = [ncbi_genbank_factory.build()]

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=records)
        otu_service.ncbi.fetch_taxonomy_record = mocker.Mock(return_value=None)

        accessions = [records[0].accession]
        otu = otu_service.create(accessions)

        assert otu is None


class TestOTUServiceCreatePlan:
    """Test plan creation failures."""

    def test_plan_creation_failure(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
        mocker: MockerFixture,
        ncbi_taxonomy_factory: type[NCBITaxonomyFactory],
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that plan creation failure returns None."""
        taxonomy = ncbi_taxonomy_factory.build(rank=NCBIRank.SPECIES)
        records = [ncbi_genbank_factory.build(source__taxid=taxonomy.id)]

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=records)
        otu_service.ncbi.fetch_taxonomy_record = mocker.Mock(return_value=taxonomy)

        mocker.patch(
            "ref_builder.services.otu.create_plan_from_records",
            return_value=None,
        )

        accessions = [records[0].accession]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu is None


class TestOTUServiceExcludeAccessions:
    """Test accession exclusion."""

    def test_ok(self, scratch_repo: Repo, mock_ncbi_client):
        """Test accession exclusion."""
        taxid = mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid

        otu_before = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_before is not None

        assert not otu_before.excluded_accessions

        otu_service = OTUService(scratch_repo, mock_ncbi_client)

        with scratch_repo.lock():
            otu_service.exclude_accessions(otu_before.id, {"DQ178608", "DQ178609"})

        otu_after = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_after is not None

        assert otu_after.excluded_accessions == {"DQ178608", "DQ178609"}


class TestOTUServiceAllowAccessions:
    """Test allowing previously excluded accessions."""

    def test_ok(self, scratch_repo: Repo, mock_ncbi_client):
        taxid = mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid

        otu_initial = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_initial is not None

        otu_service = OTUService(scratch_repo, mock_ncbi_client)

        with scratch_repo.lock():
            otu_service.exclude_accessions(otu_initial.id, {"DQ178608", "DQ178609"})

        otu_before = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_before is not None

        assert otu_before.excluded_accessions == {"DQ178608", "DQ178609"}

        with scratch_repo.lock():
            otu_service.allow_accessions(otu_before.id, {"DQ178608"})

        otu_after = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_after is not None

        assert otu_after.excluded_accessions == {"DQ178609"}
