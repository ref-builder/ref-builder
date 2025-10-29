import pytest
from pytest_mock import MockerFixture

from ref_builder.models.otu import OTU
from ref_builder.ncbi.client import NCBIClientProtocol
from ref_builder.ncbi.models import NCBIRank
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


class TestCreate:
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
        assert isinstance(otu, OTU)
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

    def test_empty_accessions(self, otu_service: OTUService):
        """Test that empty accessions list returns None."""
        assert otu_service.create([]) is None

    def test_missing_accessions(
        self,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
        otu_service: OTUService,
    ):
        """Test that missing accessions returns None."""
        records = [ncbi_genbank_factory.build()]

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(return_value=records)

        accessions = [records[0].accession, "MISSING123"]
        otu = otu_service.create(accessions)

        assert otu is None

    def test_multiple_organisms(
        self,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
        otu_service: OTUService,
    ):
        """Test that records from different organisms returns None."""
        record_1 = ncbi_genbank_factory.build(source__taxid=12345)
        record_2 = ncbi_genbank_factory.build(source__taxid=67890)

        otu_service.ncbi.fetch_genbank_records = mocker.Mock(
            return_value=[record_1, record_2]
        )

        assert otu_service.create([record_1.accession, record_2.accession]) is None

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


class TestSetPlan:
    """Test plan creation failures."""

    def test_failure(
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


class TestOTUServiceManageAccessions:
    """Test accession exclusion and allowing."""

    def test_ok(self, empty_repo: Repo, mock_ncbi_client):
        """Test excluding and allowing accessions."""
        services = Services(empty_repo, mock_ncbi_client)

        # Create OTU with only RefSeq accessions
        with empty_repo.lock():
            otu = services.otu.create(["DQ178610", "DQ178611"])

        assert otu is not None
        assert otu.accessions == {"DQ178610", "DQ178611"}
        assert not otu.excluded_accessions

        # Exclude GenBank isolate accessions (not yet in OTU)
        with empty_repo.lock():
            services.otu.exclude_accessions(otu.id, {"DQ178608", "DQ178609"})

        otu_excluded = empty_repo.get_otu(otu.id)
        assert otu_excluded is not None
        assert otu_excluded.excluded_accessions == {"DQ178608", "DQ178609"}

        # Allow one accession back
        with empty_repo.lock():
            services.otu.allow_accessions(otu.id, {"DQ178608"})

        otu_final = empty_repo.get_otu(otu.id)
        assert otu_final is not None
        assert otu_final.excluded_accessions == {"DQ178609"}


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
        mocker.patch.object(otu_service, "_promote_accessions", return_value=set())
        mocker.patch.object(
            otu_service, "_upgrade_outdated_sequences", return_value=set()
        )
        mocker.patch.object(
            otu_service._ncbi, "fetch_accessions_by_taxid", return_value=[]
        )

        # Update the OTU
        with empty_repo.lock():
            updated_otu = otu_service.update(otu.id)

        assert updated_otu is not None
        assert updated_otu.id == otu.id

    def test_otu_not_found(
        self,
        otu_service: OTUService,
        mocker: MockerFixture,
    ):
        """Test that update returns None when OTU is not found."""
        from uuid import uuid4

        result = otu_service.update(uuid4())

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
            updated_otu = otu_service.update(otu.id)

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
            otu_after = otu_service.update(otu_before.id)

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

    def test_upgrade_sequences(
        self, empty_repo: Repo, otu_service: OTUService, mock_ncbi_client
    ):
        """Test sequence version upgrade from .1 to .3."""
        # Block newer versions so only .1 is discoverable during creation
        with mock_ncbi_client.blocking(["NC_004452.2", "NC_004452.3"]):
            with empty_repo.lock():
                otu_before = otu_service.create(["NC_004452.1"])

        assert otu_before
        assert "NC_004452" in otu_before.accessions
        assert len(otu_before.isolates) == 1

        sequence_before = otu_before.get_sequence("NC_004452")
        assert sequence_before
        assert sequence_before.accession.version == 1

        isolate_before = otu_before.isolates[0]
        assert isolate_before.accessions == {"NC_004452"}

        # Now .2 and .3 are discoverable - update should upgrade to .3
        with empty_repo.lock():
            otu_after = otu_service.update(otu_before.id)

        assert otu_after
        assert "NC_004452" in otu_after.accessions

        # Find the sequence with the highest version
        nc_sequences = [
            seq for seq in otu_after.sequences if seq.accession.key == "NC_004452"
        ]
        assert len(nc_sequences) == 1

        sequence_after = nc_sequences[0]
        assert sequence_after.accession.version == 3

        # Verify the isolate was updated
        isolate_after = otu_after.get_isolate(isolate_before.id)
        assert isolate_after
        assert isolate_after.accessions == {"NC_004452"}
