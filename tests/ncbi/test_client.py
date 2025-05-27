import datetime

import pytest
from syrupy import SnapshotAssertion

from ref_builder.ncbi.client import NCBIClient, TaxonLevelError
from ref_builder.utils import Accession


class TestFetchGenbank:
    @pytest.mark.ncbi()
    def test_fetch_genbank_records_from_ncbi(
        self,
        snapshot: SnapshotAssertion,
        uncached_ncbi_client: NCBIClient,
    ):
        accessions = ["NC_036587", "MT240513", "AB017504"]

        assert all(
            uncached_ncbi_client.cache.load_genbank_record(accession) is None
            for accession in accessions
        )

        for record in uncached_ncbi_client.fetch_genbank_records(accessions):
            assert record.model_dump() == snapshot

        assert all(
            uncached_ncbi_client.cache.load_genbank_record(accession)
            for accession in accessions
        )

    def test_fetch_genbank_records_from_cache(
        self,
        snapshot: SnapshotAssertion,
        scratch_ncbi_client: NCBIClient,
    ):
        for record in scratch_ncbi_client.fetch_genbank_records(
            ["NC_036587", "MT240513", "AB017504"]
        ):
            assert record.model_dump() == snapshot

    @pytest.mark.ncbi()
    def test_fetch_partially_cached_genbank_records(
        self,
        snapshot: SnapshotAssertion,
        uncached_ncbi_client: NCBIClient,
    ):
        """Test that the client can fetch records when some are cached and some aren't.

        Start by fetching one accession, then ensure it was cached. Then, fetch is again
        along with two other accessions.
        """
        assert uncached_ncbi_client.cache.load_genbank_record("NC_036587") is None
        assert uncached_ncbi_client.fetch_genbank_records(["NC_036587"])
        assert uncached_ncbi_client.cache.load_genbank_record("NC_036587")

        for record in uncached_ncbi_client.fetch_genbank_records(
            ["NC_036587", "MT240513", "AB017504"]
        ):
            assert record.model_dump() == snapshot

    def test_fetch_multi_source_record(
        self, scratch_ncbi_client: NCBIClient, snapshot: SnapshotAssertion
    ):
        """Test the parsing of a Genbank record containing multiple source features."""
        record = scratch_ncbi_client.fetch_genbank_records(["Y18556"])[0]

        assert record.source.organism == record.organism

        assert record.model_dump() == snapshot

    @pytest.mark.ncbi()
    def test_fetch_non_existent_accession(self, scratch_ncbi_client: NCBIClient):
        """Test that the client returns an empty list when the fetched accession does
        not exist.
        """
        assert scratch_ncbi_client.fetch_genbank_records(["paella"]) == []


class TestClientFetchRawGenbank:
    @pytest.mark.ncbi()
    @pytest.mark.parametrize(
        "accessions",
        [
            ["AB017504", "MH200607", "MK431779", "NC_003355"],
            ["NC_036587", "MT240513", "MT240490"],
        ],
    )
    def test_fetch_raw_via_accessions(self, accessions: list[str]):
        records = NCBIClient.fetch_unvalidated_genbank_records(accessions)

        for record in records:
            assert record.get("GBSeq_locus")
            assert record.get("GBSeq_sequence")

    @pytest.mark.ncbi()
    @pytest.mark.parametrize(
        "accessions",
        [
            ["paella", "MH200607", "MK431779", "NC_003355"],
            ["friday", "MT240513", "MT240490"],
        ],
    )
    def test_fetch_raw_via_accessions_partial(self, accessions: list[str]):
        records = NCBIClient.fetch_unvalidated_genbank_records(accessions)

        assert len(records) == len(accessions) - 1

        for record in records:
            assert record.get("GBSeq_locus", None)
            assert record.get("GBSeq_sequence", None)

    @pytest.mark.ncbi()
    def test_fetch_raw_via_accessions_fail(self):
        records = NCBIClient.fetch_unvalidated_genbank_records(
            ["friday", "paella", "111"],
        )

        assert not records


@pytest.mark.ncbi()
class TestFetchAccessionsByTaxid:
    """Test NCBI ESearch interface functionality."""

    def test_ok(self):
        """Test that the client can fetch accessions by taxid."""
        assert NCBIClient.fetch_accessions_by_taxid(1198450) == [
            "NC_038797.1",
            "NC_038796.1",
            "JQ821387.1",
            "JQ821386.1",
        ]

    def test_nonexistent(self):
        """Test that the client returns an empty list when the taxid does not exist."""
        assert NCBIClient.fetch_accessions_by_taxid(99999999) == []

    def test_giant_ok(self):
        """Test that the client returns a list of results when the accession list
        is long enough to force pagination.
        """
        assert len(NCBIClient.fetch_accessions_by_taxid(12585)) > 1000

    def test_refseq_limit(self, uncached_ncbi_client: NCBIClient):
        """Test that RefSeq-only searches return only RefSeq accessions."""

        for raw_accession in NCBIClient.fetch_accessions_by_taxid(
            438782, refseq_only=True
        ):
            assert raw_accession[:3] == "NC_"

    def test_esearch_limit_by_length(self, uncached_ncbi_client: NCBIClient):
        """Test that narrow search terms fetch a smaller subset of accessions
        than wide search terms.
        """
        taxid = 12232

        segment_length = 9591

        wide_filtered_accessions = set(
            uncached_ncbi_client.fetch_accessions_by_taxid(
                taxid,
                sequence_min_length=segment_length - 6,
                sequence_max_length=segment_length + 6,
            )
        )

        assert all(
            [
                segment_length - 6 <= len(record.sequence) <= segment_length + 6
                for record in uncached_ncbi_client.fetch_genbank_records(
                    wide_filtered_accessions
                )
            ]
        )

        narrow_filtered_accessions = set(
            NCBIClient.fetch_accessions_by_taxid(
                taxid,
                sequence_min_length=segment_length - 1,
                sequence_max_length=segment_length + 1,
            )
        )

        assert all(
            [
                segment_length - 1 <= len(record.sequence) <= segment_length + 1
                for record in uncached_ncbi_client.fetch_genbank_records(
                    narrow_filtered_accessions
                )
            ]
        )

        assert narrow_filtered_accessions.issubset(wide_filtered_accessions)

        assert wide_filtered_accessions - narrow_filtered_accessions

    def test_esearch_limit_by_date(self, uncached_ncbi_client: NCBIClient):
        """Test that a more recent modification date search term fetches a smaller subset of accessions"""
        taxid = 12232

        unfiltered_accessions = set(
            uncached_ncbi_client.fetch_accessions_by_taxid(taxid)
        )

        after_2024_accessions = set(
            NCBIClient.fetch_accessions_by_taxid(
                taxid,
                modification_date_start=datetime.date(year=2024, month=1, day=1),
            )
        )

        assert len(after_2024_accessions) < len(unfiltered_accessions)

        assert after_2024_accessions.issubset(unfiltered_accessions)

        after_2025_accessions = set(
            NCBIClient.fetch_accessions_by_taxid(
                taxid,
                modification_date_start=datetime.date(year=2025, month=1, day=1),
            )
        )

        assert len(after_2025_accessions) < len(after_2024_accessions)

        assert after_2025_accessions.issubset(after_2024_accessions)


@pytest.mark.ncbi()
class TestFetchTaxonomy:
    def test_ok(
        self,
        snapshot: SnapshotAssertion,
        uncached_ncbi_client: NCBIClient,
    ):
        """Test that the client can fetch a taxonomy record and that it is cached."""
        # Make sure the taxid is not cached.
        assert uncached_ncbi_client.cache.load_taxonomy(1198450) is None

        # Fetch the taxonomy record and check its contents.
        assert (
            uncached_ncbi_client.fetch_taxonomy_record(1198450).model_dump() == snapshot
        )

        # Make sure the taxid is now cached.
        assert uncached_ncbi_client.cache.load_taxonomy(1198450)

    def test_fetch_taxonomy_by_name(self):
        """Test that the client can fetch a taxid by name."""
        assert (
            NCBIClient.fetch_taxonomy_id_by_name("Rhynchosia golden mosaic virus")
            == 117198
        )

    def test_not_found(self, uncached_ncbi_client: NCBIClient):
        """Test that the client returns None when the taxid does not exist."""
        assert uncached_ncbi_client.fetch_taxonomy_record(99999999) is None

    def test_rank_too_high_fail(
        self,
        uncached_ncbi_client: NCBIClient,
    ):
        """Test that if the requested taxonomy rank is above species level,
        the correct error is raised.
        """
        with pytest.raises(TaxonLevelError):
            uncached_ncbi_client.fetch_taxonomy_record(190729)


@pytest.mark.ncbi()
def test_fetch_spelling():
    """Test that the client can fetch the correct spelling of a virus name."""
    assert (
        NCBIClient.fetch_spelling(name="Angelica bush stunt virus")
        == "angelica bushy stunt virus"
    )


def test_filter_accessions():
    """Test that accession filter only allows valid versioned accessions to pass."""
    assert NCBIClient.filter_accessions(
        ["friday", "paella", "111", "NC_038797.1", "NO_000001"]
    ) == {Accession("NC_038797", 1)}
