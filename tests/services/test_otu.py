import pytest
from pytest_mock import MockerFixture

from ref_builder.ncbi.client import NCBIClientProtocol
from ref_builder.ncbi.models import NCBIRank
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.repo import Repo
from ref_builder.services.cls import Services
from ref_builder.services.otu import OTUService
from tests.fixtures.factories import (
    NCBIGenbankFactory,
    NCBITaxonomyFactory,
)


@pytest.fixture
def otu_service(empty_repo: Repo, mock_ncbi_client: NCBIClientProtocol) -> OTUService:
    """Create an OTUService with a mock NCBI client."""
    services = Services(empty_repo, mock_ncbi_client)
    return services.otu


class TestOTUServiceCreate:
    """Test successful OTU creation scenarios."""

    def test_ok(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
    ):
        """Test basic OTU creation with monopartite genome."""
        # Use tobacco_mosaic_virus - monopartite with refseq
        # Note: TMV taxid 12242 is promoted to species level 3432891
        accessions = ["NC_001367"]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu is not None
        assert isinstance(otu, OTUBuilder)
        assert otu.taxid == 3432891  # Species level (auto-promoted from 12242)
        assert otu.name == "Tobamovirus tabaci"

        otus = list(empty_repo.iter_otus())
        assert len(otus) == 1
        assert otus[0].id == otu.id

    def test_ok_multipartite(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
    ):
        """Test OTU creation with multipartite genome."""
        # Use abaca_bunchy_top_virus - multipartite with 6 segments
        accessions = [
            "NC_010314",
            "NC_010315",
            "NC_010316",
            "NC_010317",
            "NC_010318",
            "NC_010319",
        ]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu is not None
        assert len(otu.plan.segments) == 6
        assert len(otu.isolates) == 1
        assert len(otu.isolates[0].sequences) == 6

    def test_refseq_exclusion(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
    ):
        """Test RefSeq comment parsing and old accession exclusion."""
        # Use TMV RefSeq which references V01408 in its comment
        accessions = ["NC_001367"]

        with empty_repo.lock():
            otu = otu_service.create(accessions)

        assert otu is not None
        assert "V01408" in otu.excluded_accessions


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

        services = Services(scratch_repo, mock_ncbi_client)

        with scratch_repo.lock():
            services.otu.exclude_accessions(otu_before.id, {"DQ178608", "DQ178609"})

        otu_after = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_after is not None

        assert otu_after.excluded_accessions == {"DQ178608", "DQ178609"}


class TestOTUServiceAllowAccessions:
    """Test allowing previously excluded accessions."""

    def test_ok(self, scratch_repo: Repo, mock_ncbi_client):
        taxid = mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid

        otu_initial = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_initial is not None

        services = Services(scratch_repo, mock_ncbi_client)

        with scratch_repo.lock():
            services.otu.exclude_accessions(otu_initial.id, {"DQ178608", "DQ178609"})

        otu_before = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_before is not None

        assert otu_before.excluded_accessions == {"DQ178608", "DQ178609"}

        with scratch_repo.lock():
            services.otu.allow_accessions(otu_before.id, {"DQ178608"})

        otu_after = scratch_repo.get_otu_by_taxid(taxid)
        assert otu_after is not None

        assert otu_after.excluded_accessions == {"DQ178609"}


class TestOTUServiceUpdate:
    """Test OTU update functionality."""

    def test_ok(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
        mocker: MockerFixture,
    ):
        """Test basic OTU update with mock data."""
        # Create an OTU first
        with empty_repo.lock():
            otu = otu_service.create(["NC_001367"])

        assert otu is not None

        # Mock the update functions to avoid real NCBI calls
        mocker.patch(
            "ref_builder.services.otu.promote_otu_accessions", return_value=set()
        )
        mocker.patch(
            "ref_builder.services.otu.upgrade_outdated_sequences_in_otu",
            return_value=set(),
        )
        mocker.patch.object(
            otu_service._ncbi, "fetch_accessions_by_taxid", return_value=[]
        )

        # Update the OTU
        with empty_repo.lock():
            updated_otu = otu_service.update(otu.id, ignore_cache=True)

        assert updated_otu is not None
        assert updated_otu.id == otu.id

    def test_otu_not_found(
        self,
        otu_service: OTUService,
        mocker: MockerFixture,
    ):
        """Test that update returns None when OTU is not found."""
        from uuid import uuid4

        result = otu_service.update(uuid4(), ignore_cache=True)

        assert result is None

    def test_promotion(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
    ):
        """Test that GenBank accessions are promoted to RefSeq during update."""
        # Create OTU with GenBank isolate V01408 (which has RefSeq equivalent NC_001367)
        with empty_repo.lock():
            otu = otu_service.create(["V01408"])

        assert otu is not None
        assert otu.accessions == {"V01408"}

        isolate_before = otu.isolates[0]
        assert isolate_before.accessions == {"V01408"}

        # Update should promote V01408 to NC_001367 and also discover other isolates
        with empty_repo.lock():
            updated_otu = otu_service.update(otu.id, ignore_cache=True)

        assert updated_otu is not None
        assert updated_otu.id == otu.id
        assert "NC_001367" in updated_otu.accessions
        assert updated_otu.excluded_accessions == {"V01408"}

        # The original isolate should be promoted
        isolate_after = updated_otu.get_isolate(isolate_before.id)
        assert isolate_after
        assert isolate_after.accessions == {"NC_001367"}

        # New isolates should have been added
        assert len(updated_otu.isolate_ids) > len(otu.isolate_ids)

    def test_add_new_isolates(
        self,
        empty_repo: Repo,
        otu_service: OTUService,
    ):
        """Test that update discovers and adds new isolates from NCBI."""
        # Create OTU with only RefSeq (TMV has 5 isolates in mock data)
        with empty_repo.lock():
            otu_before = otu_service.create(["NC_001367"])

        assert otu_before is not None
        assert otu_before.accessions == {"NC_001367"}
        assert len(otu_before.isolates) == 1
        initial_isolate_count = len(otu_before.isolate_ids)

        # Update should discover and add the GenBank isolates
        with empty_repo.lock():
            otu_after = otu_service.update(otu_before.id, ignore_cache=True)

        assert otu_after is not None
        assert otu_after.id == otu_before.id

        # Should have added new isolates (TMV has 4 GenBank isolates in mock data)
        assert len(otu_after.isolate_ids) > initial_isolate_count
        # Should now include both RefSeq and GenBank accessions
        assert "NC_001367" in otu_after.accessions
        # Check some of the GenBank isolates were added
        genbank_isolates = {"OQ953825", "HE818414", "AJ011933", "AF395128"}
        assert genbank_isolates.intersection(otu_after.accessions), (
            "Expected some GenBank isolates to be added"
        )
