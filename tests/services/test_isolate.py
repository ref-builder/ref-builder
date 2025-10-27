from uuid import uuid4

from pytest_mock import MockerFixture

from ref_builder.models.isolate import IsolateName, IsolateNameType
from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.repo import Repo
from ref_builder.services.cls import Services
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory
from tests.fixtures.mock_ncbi_client import MockNCBIClient


class TestIsolateServiceCreate:
    """Test successful isolate creation scenarios."""

    def test_ok(
        self,
        empty_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
    ):
        """Test basic isolate creation with monopartite genome."""
        services = Services(empty_repo, mock_ncbi_client)

        with empty_repo.lock():
            otu = services.otu.create(["NC_003355"])

        assert otu

        with empty_repo.lock():
            isolate = services.isolate.create(["MH200607"])

        assert isolate is not None
        assert isinstance(isolate, IsolateBuilder)
        assert isolate.name == IsolateName(IsolateNameType.ISOLATE, "WMoV-6.3")
        assert len(isolate.sequences) == 1

    def test_multipartite(
        self,
        empty_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
    ):
        """Test isolate creation with multipartite genome."""
        services = Services(empty_repo, mock_ncbi_client)

        with empty_repo.lock():
            otu = services.otu.create(
                [
                    "EF546808",
                    "EF546809",
                    "EF546810",
                    "EF546811",
                    "EF546812",
                    "EF546813",
                ]
            )

            assert otu

            isolate = services.isolate.create(
                [
                    "EF546802",
                    "EF546803",
                    "EF546804",
                    "EF546805",
                    "EF546806",
                    "EF546807",
                ]
            )

        assert isolate
        assert len(isolate.sequences) == 6

        name = isolate.name

        assert name
        assert name.value == "Q1108"

    def test_refseq_exclusion(
        self,
        empty_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
    ):
        """Test RefSeq promotion and exclusion of original GenBank accessions for multipartite isolate."""
        services = Services(empty_repo, mock_ncbi_client)

        with empty_repo.lock():
            otu = services.otu.create(
                [
                    "EF546808",
                    "EF546809",
                    "EF546810",
                    "EF546811",
                    "EF546812",
                    "EF546813",
                ]
            )

            assert otu

            isolate = services.isolate.create(
                [
                    "NC_010314",
                    "NC_010315",
                    "NC_010316",
                    "NC_010317",
                    "NC_010318",
                    "NC_010319",
                ]
            )

        assert isolate

        otu = empty_repo.get_otu(otu.id)

        assert otu
        assert "EF546808" in otu.excluded_accessions
        assert "EF546809" in otu.excluded_accessions
        assert "EF546810" in otu.excluded_accessions
        assert "EF546811" in otu.excluded_accessions
        assert "EF546812" in otu.excluded_accessions
        assert "EF546813" in otu.excluded_accessions


class TestIsolateServiceCreateValidation:
    """Test input validation failures."""

    def test_otu_not_found(
        self,
        empty_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that non-existent OTU returns None when taxid doesn't match any OTU."""
        services = Services(empty_repo, mock_ncbi_client)

        # Mock a record with a taxid that doesn't match any OTU
        source = NCBISourceFactory.build(taxid=999999)
        record = ncbi_genbank_factory.build(source=source, accession="AB123456")

        mocker.patch.object(
            mock_ncbi_client,
            "fetch_genbank_records",
            return_value=[record],
        )

        with empty_repo.lock():
            isolate = services.isolate.create(["AB123456"])

        assert isolate is None

    def test_no_records_fetched(
        self,
        empty_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
    ):
        """Test that failure to fetch records returns None."""
        services = Services(empty_repo, mock_ncbi_client)

        with empty_repo.lock():
            otu = services.otu.create(["NC_003355"])

            assert otu

            isolate = services.isolate.create(["MISS123"])

        assert isolate is None

    def test_mixed_refseq_nonrefseq(
        self,
        empty_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
    ):
        """Test that mixing RefSeq and non-RefSeq sequences returns None."""
        services = Services(empty_repo, mock_ncbi_client)

        with empty_repo.lock():
            otu = services.otu.create(
                [
                    "EF546808",
                    "EF546809",
                    "EF546810",
                    "EF546811",
                    "EF546812",
                    "EF546813",
                ]
            )

            assert otu

            isolate = services.isolate.create(
                [
                    "EF546803",
                    "EF546804",
                    "EF546805",
                    "NC_010314",
                    "NC_010318",
                    "NC_010319",
                ]
            )

        with empty_repo.lock():
            isolate = services.isolate.create([])

        assert isolate is None

    def test_multiple_isolate_names(
        self,
        scratch_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that multiple isolate names in records returns None."""
        services = Services(scratch_repo, mock_ncbi_client)

        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.wasabi_mottle_virus.taxid
        )

        assert otu

        source_1 = NCBISourceFactory.build(
            taxid=otu.taxid,
            isolate="Isolate-1",
        )
        source_2 = NCBISourceFactory.build(
            taxid=otu.taxid,
            isolate="Isolate-2",
        )

        record_1 = ncbi_genbank_factory.build(source=source_1, accession="AB123456")
        record_2 = ncbi_genbank_factory.build(source=source_2, accession="AB123457")

        mocker.patch.object(
            mock_ncbi_client,
            "fetch_genbank_records",
            return_value=[record_1, record_2],
        )

        with scratch_repo.lock():
            isolate = services.isolate.create([record_1.accession, record_2.accession])

        assert isolate is None

    def test_duplicate_isolate_name(
        self,
        scratch_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test that duplicate isolate name returns None."""
        services = Services(scratch_repo, mock_ncbi_client)

        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.wasabi_mottle_virus.taxid
        )

        assert otu

        # Use an existing isolate name from scratch_repo
        existing_isolate = otu.isolates[0]
        existing_name = existing_isolate.name

        source = NCBISourceFactory.build(
            taxid=otu.taxid,
            isolate=existing_name.value,
        )
        record = ncbi_genbank_factory.build(
            source=source,
            accession="AB999999",
            refseq=False,
        )

        mocker.patch.object(
            mock_ncbi_client,
            "fetch_genbank_records",
            return_value=[record],
        )

        with scratch_repo.lock():
            isolate = services.isolate.create([record.accession])

        assert isolate is None

    def test_plan_mismatch(
        self,
        scratch_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
    ):
        """Test that isolate cannot be added if it doesn't match OTU plan."""
        services = Services(scratch_repo, mock_ncbi_client)

        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        )

        assert otu

        with scratch_repo.lock():
            isolate = services.isolate.create(["NC_003355"])

        assert isolate is None

        otu_after = scratch_repo.get_otu(otu.id)
        assert otu_after
        assert "NC_003355" not in otu_after.accessions

    def test_sequence_too_short(
        self,
        scratch_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
        ncbi_source_factory: type[NCBISourceFactory],
    ):
        """Test that sequences below length tolerance are rejected."""
        services = Services(scratch_repo, mock_ncbi_client)

        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        )

        assert otu
        assert len(otu.plan.required_segments) == 2

        segment = otu.plan.required_segments[0]
        too_short_length = int(segment.length * 0.5)

        records = []
        for i, seg in enumerate(otu.plan.required_segments):
            source = ncbi_source_factory.build(
                taxid=otu.taxid,
                organism=otu.name,
                segment=str(seg.name),
            )
            sequence_length = too_short_length if i == 0 else seg.length
            records.append(
                ncbi_genbank_factory.build(
                    accession=f"AB{100000 + i}",
                    sequence="A" * sequence_length,
                    source=source,
                )
            )

        mocker.patch.object(
            mock_ncbi_client,
            "fetch_genbank_records",
            return_value=records,
        )

        with scratch_repo.lock():
            isolate = services.isolate.create([r.accession for r in records])

        assert isolate is None

        otu_after = scratch_repo.get_otu(otu.id)
        assert otu_after
        assert otu_after.isolate_ids == otu.isolate_ids

    def test_sequence_too_long(
        self,
        scratch_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
        ncbi_source_factory: type[NCBISourceFactory],
    ):
        """Test that sequences above length tolerance are rejected."""
        services = Services(scratch_repo, mock_ncbi_client)

        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid
        )

        assert otu
        assert len(otu.plan.required_segments) == 2

        segment = otu.plan.required_segments[0]
        too_long_length = int(segment.length * 2.0)

        records = []
        for i, seg in enumerate(otu.plan.required_segments):
            source = ncbi_source_factory.build(
                taxid=otu.taxid,
                organism=otu.name,
                segment=str(seg.name),
            )
            sequence_length = too_long_length if i == 0 else seg.length
            records.append(
                ncbi_genbank_factory.build(
                    accession=f"AB{200000 + i}",
                    sequence="A" * sequence_length,
                    source=source,
                )
            )

        mocker.patch.object(
            mock_ncbi_client,
            "fetch_genbank_records",
            return_value=records,
        )

        with scratch_repo.lock():
            isolate = services.isolate.create([r.accession for r in records])

        assert isolate is None

        otu_after = scratch_repo.get_otu(otu.id)
        assert otu_after
        assert otu_after.isolate_ids == otu.isolate_ids


class TestIsolateServicePromotion:
    """Test RefSeq promotion logic."""

    def test_promote_existing_isolate(
        self,
        scratch_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
        mocker: MockerFixture,
        ncbi_genbank_factory: type[NCBIGenbankFactory],
    ):
        """Test promoting sequences when adding RefSeq version of existing isolate."""
        services = Services(scratch_repo, mock_ncbi_client)

        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.wasabi_mottle_virus.taxid
        )

        # Get an existing isolate
        existing_isolate = otu.isolates[0]
        existing_name = existing_isolate.name

        # Create RefSeq record with same isolate name
        source = NCBISourceFactory.build(
            taxid=otu.taxid,
            isolate=existing_name.value,
        )
        refseq_record = ncbi_genbank_factory.build(
            source=source,
            accession="NC_999999",
            refseq=True,
            comment="PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review. "
            f"The reference sequence is identical to {existing_isolate.sequences[0].accession.key},COMPLETENESS: full length.",
        )

        mocker.patch.object(
            mock_ncbi_client,
            "fetch_genbank_records",
            return_value=[refseq_record],
        )

        # Mock promote_otu_accessions_from_records to return promoted accessions
        mocker.patch(
            "ref_builder.services.isolate.promote_otu_accessions_from_records",
            return_value={refseq_record.accession},
        )

        with scratch_repo.lock():
            isolate = services.isolate.create([refseq_record.accession])

        assert isolate is not None
        assert isolate.id == existing_isolate.id


class TestIsolateServiceDelete:
    """Test isolate deletion."""

    def test_ok(self, scratch_repo: Repo, mock_ncbi_client: MockNCBIClient):
        """Test successful isolate deletion."""
        services = Services(scratch_repo, mock_ncbi_client)

        otu_before = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.wasabi_mottle_virus.taxid
        )
        isolate_id = next(iter(otu_before.isolate_ids))
        isolate_before = otu_before.get_isolate(isolate_id)

        with scratch_repo.lock():
            result = services.isolate.delete(isolate_id)

        assert result is True

        otu_after = scratch_repo.get_otu(otu_before.id)
        assert isolate_id not in otu_after.isolate_ids
        assert otu_after.get_isolate(isolate_id) is None
        assert isolate_before.accessions not in otu_after.accessions
        assert len(otu_after.isolate_ids) == len(otu_before.isolate_ids) - 1

    def test_otu_not_found(
        self,
        empty_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
    ):
        """Test deletion with non-existent OTU returns False."""
        services = Services(empty_repo, mock_ncbi_client)

        result = services.isolate.delete(uuid4())

        assert result is False

    def test_isolate_not_found(
        self,
        scratch_repo: Repo,
        mock_ncbi_client: MockNCBIClient,
    ):
        """Test deletion with non-existent isolate returns False."""
        services = Services(scratch_repo, mock_ncbi_client)

        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.wasabi_mottle_virus.taxid
        )

        with scratch_repo.lock():
            result = services.isolate.delete(uuid4())

        assert result is False
