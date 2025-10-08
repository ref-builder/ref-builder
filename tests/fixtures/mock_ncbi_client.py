"""Mock NCBI client for testing without real API calls or file cache."""

from collections.abc import Collection

from structlog import get_logger

from ref_builder.ncbi.client import NCBIClientProtocol
from ref_builder.ncbi.models import NCBIGenbank, NCBITaxonomy
from ref_builder.utils import Accession
from tests.fixtures.ncbi_mock_data import MOCK_GENBANK_RECORDS, MOCK_TAXONOMY_RECORDS

logger = get_logger("tests.mock_ncbi_client")


class MockDataNotFoundError(LookupError):
    """Raised when requested data is not available in mock.

    This indicates missing test data, not a test implementation error.
    Add the missing accession/taxid to tests/fixtures/ncbi_mock_data.py
    """
    pass


class MockNCBIClient:
    """Mock NCBIClient that returns hardcoded test data.

    Implements NCBIClientProtocol for type compatibility with real NCBIClient.
    Raises clear errors when unmocked data is requested to guide test development.
    """

    def __init__(self, ignore_cache: bool = False) -> None:
        """Initialize mock client.

        Args:
            ignore_cache: Ignored in mock (for interface compatibility)
        """
        self.ignore_cache = ignore_cache
        self.cache = None  # Mock has no cache

    def fetch_genbank_records(
        self,
        accessions: Collection[str | Accession],
    ) -> list[NCBIGenbank]:
        """Fetch mocked GenBank records.

        Args:
            accessions: List of accessions to fetch

        Returns:
            List of NCBIGenbank records

        Raises:
            KeyError: If any accession is not mocked
        """
        if not accessions:
            return []

        records = []

        for accession in accessions:
            # Normalize accession to string without version
            if isinstance(accession, Accession):
                acc_key = accession.key
            else:
                try:
                    acc_key = Accession.from_string(accession).key
                except ValueError:
                    # Not a versioned accession, use as-is
                    acc_key = accession

            try:
                records.append(MOCK_GENBANK_RECORDS[acc_key])
            except KeyError:
                available = ", ".join(sorted(MOCK_GENBANK_RECORDS.keys())[:10])
                msg = (
                    f"Accession '{acc_key}' not in mock data. "
                    f"Add it to tests/fixtures/ncbi_mock_data.py\n"
                    f"Available accessions (first 10): {available}..."
                )
                raise MockDataNotFoundError(msg) from None

        return records

    def fetch_taxonomy_record(self, taxid: int) -> NCBITaxonomy | None:
        """Fetch mocked taxonomy record.

        Matches real NCBIClient interface by returning None if not found.

        Args:
            taxid: NCBI Taxonomy ID

        Returns:
            NCBITaxonomy record or None if not mocked
        """
        try:
            return MOCK_TAXONOMY_RECORDS[taxid]
        except KeyError:
            available = ", ".join(
                str(t) for t in sorted(MOCK_TAXONOMY_RECORDS.keys())[:10]
            )
            logger.error(
                "Taxid not in mock data",
                taxid=taxid,
                available=available,
                hint="Add it to tests/fixtures/ncbi_mock_data.py",
            )
            return None
