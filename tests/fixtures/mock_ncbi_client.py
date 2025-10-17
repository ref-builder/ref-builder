"""Mock NCBI client for testing without real API calls or file cache."""

from collections.abc import Collection

from structlog import get_logger

from ref_builder.ncbi.models import NCBIGenbank, NCBITaxonomy
from ref_builder.utils import Accession

logger = get_logger("tests.mock_ncbi_client")


class MockDataNotFoundError(LookupError):
    """Raised when requested data is not available in mock.

    This indicates missing test data, not a test implementation error.
    Add the missing accession/taxid to tests/fixtures/ncbi/ modules.
    """



class MockOTUBuilder:
    """Builder for adding isolates to a mock OTU."""

    def __init__(
        self,
        client: "MockNCBIClient",
        taxid: int,
        refseq: list[str],
    ) -> None:
        """Initialize builder.

        Args:
            client: Parent MockNCBIClient instance
            taxid: NCBI taxonomy ID
            refseq: RefSeq accessions for this OTU

        """
        self._client = client
        self._taxid = taxid
        self._refseq = refseq
        self._isolate_groups: list[list[str]] = []

    def add_isolate(
        self,
        accessions: list[str],
        genbank: list[NCBIGenbank],
    ) -> "MockOTUBuilder":
        """Add an isolate group to this OTU.

        Args:
            accessions: Accessions for this isolate
            genbank: GenBank records for this isolate. Use empty list if none available.

        Returns:
            Self for chaining

        """
        for gb in genbank:
            self._client._genbank_records[gb.accession] = gb

        self._isolate_groups.append(accessions)

        # Update OTU structure
        for i, (tid, refseq, _) in enumerate(self._client._otu_structure):
            if tid == self._taxid:
                self._client._otu_structure[i] = (
                    self._taxid,
                    refseq,
                    self._isolate_groups.copy(),
                )
                break

        return self


class MockNCBIClient:
    """Mock NCBIClient that returns hardcoded test data.

    Implements NCBIClientProtocol for type compatibility with real NCBIClient.
    Raises clear errors when unmocked data is requested to guide test development.

    Also serves as a registry for mock data, populated by per-OTU modules.
    """

    def __init__(self, ignore_cache: bool = False) -> None:
        """Initialize mock client.

        Args:
            ignore_cache: Ignored in mock (for interface compatibility)

        """
        self.ignore_cache = ignore_cache
        self.cache = None
        self._genbank_records: dict[str, NCBIGenbank] = {}
        self._taxonomy_records: dict[int, NCBITaxonomy] = {}
        self._otu_structure: list[tuple[int, list[str], list[list[str]]]] = []

    def add_otu(
        self,
        taxid: int,
        name: str,
        refseq: list[str],
        taxonomy: NCBITaxonomy | None = None,
        genbank: list[NCBIGenbank] | None = None,
    ) -> MockOTUBuilder:
        """Register an OTU and return builder for adding isolates.

        Args:
            taxid: NCBI taxonomy ID
            name: Organism name
            refseq: RefSeq accessions
            taxonomy: Optional taxonomy record. If None, must be filled by refresh.
            genbank: GenBank records for refseq. Use empty list if none available.

        Returns:
            Builder for adding isolate groups

        """
        if taxonomy:
            self._taxonomy_records[taxid] = taxonomy

        for gb in genbank or []:
            self._genbank_records[gb.accession] = gb

        # Initialize OTU structure entry
        self._otu_structure.append((taxid, refseq, []))

        return MockOTUBuilder(self, taxid, refseq)

    def get_otu_structure(self) -> list[tuple[int, list[str], list[list[str]]]]:
        """Get OTU structure for test fixtures.

        Returns:
            List of (taxid, refseq_accessions, isolate_groups) tuples

        """
        return self._otu_structure

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
            MockDataNotFoundError: If any accession is not mocked

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
                records.append(self._genbank_records[acc_key])
            except KeyError:
                available = ", ".join(sorted(self._genbank_records.keys())[:10])
                msg = (
                    f"Accession '{acc_key}' not in mock data. "
                    f"Add it to tests/fixtures/ncbi/ modules\n"
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
            return self._taxonomy_records[taxid]
        except KeyError:
            available = ", ".join(
                str(t) for t in sorted(self._taxonomy_records.keys())[:10]
            )
            logger.error(
                "Taxid not in mock data",
                taxid=taxid,
                available=available,
                hint="Add it to tests/fixtures/ncbi/ modules",
            )
            return None


# Module-level singleton for use by per-OTU modules
mock_ncbi_client = MockNCBIClient()
