"""Mock NCBI client for testing without real API calls or file cache."""

import json
from collections.abc import Collection
from pathlib import Path

from structlog import get_logger

from ref_builder.models.accession import Accession
from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.ncbi.models import (
    NCBIGenbank,
    NCBIRank,
    NCBITaxonomy,
)
from tests.fixtures.ncbi.models import OTURegistry

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

    def add_incompatible_isolate(
        self,
        accessions: list[str],
        genbank: list[NCBIGenbank],
    ) -> "MockOTUBuilder":
        """Add an incompatible isolate that will fail OTU validation.

        These records are registered in mock data (fetchable) but not included
        in OTU structure, allowing tests to verify validation error handling.

        Args:
            accessions: Accessions for this incompatible isolate
            genbank: GenBank records for this isolate. Use empty list if none available.

        Returns:
            Self for chaining

        """
        for gb in genbank:
            self._client._genbank_records[gb.accession] = gb

        if self._taxid not in self._client._incompatible_isolates:
            self._client._incompatible_isolates[self._taxid] = []

        self._client._incompatible_isolates[self._taxid].append(accessions)

        return self


class MockNCBIClient:
    """Mock NCBIClient that returns hardcoded test data.

    Implements NCBIClientProtocol for type compatibility with real NCBIClient.
    Raises clear errors when unmocked data is requested to guide test development.

    Also serves as a registry for mock data, populated by per-OTU modules.
    """

    def __init__(
        self,
        ignore_cache: bool = False,
        manifest: type | None = None,
        data_dir: Path | None = None,
    ) -> None:
        """Initialize mock client.

        Args:
            ignore_cache: Ignored in mock (for interface compatibility)
            manifest: Optional manifest class declaring test OTUs
            data_dir: Optional directory containing JSON data files

        """
        self.ignore_cache = ignore_cache
        self.cache = None
        self._genbank_records: dict[str, NCBIGenbank] = {}
        self._taxonomy_records: dict[int, NCBITaxonomy] = {}
        self._otu_structure: list[tuple[int, list[str], list[list[str]]]] = []
        self._incompatible_isolates: dict[int, list[list[str]]] = {}

        if manifest and data_dir:
            self._load_json_data(manifest, data_dir)
        else:
            self.otus = None

    def _load_json_data(self, manifest: type, data_dir: Path) -> None:
        """Load JSON data files and populate registries.

        Args:
            manifest: Manifest class declaring test OTUs
            data_dir: Directory containing JSON data files

        """
        for json_file in data_dir.glob("*.json"):
            data = json.loads(json_file.read_text())

            taxonomy = NCBITaxonomy.model_validate(data["taxonomy"])
            self._taxonomy_records[taxonomy.id] = taxonomy

            # Mirror OTUService.create() behavior: if taxonomy is not species-level,
            # register the species-level ancestor so it can be fetched
            if taxonomy.rank != NCBIRank.SPECIES:
                species_lineage = taxonomy.species
                # Create a minimal species-level taxonomy record
                species_taxonomy = NCBITaxonomy(
                    id=species_lineage.id,
                    name=species_lineage.name,
                    rank=NCBIRank.SPECIES,
                    lineage=[],
                )
                self._taxonomy_records[species_lineage.id] = species_taxonomy

            for acc, gb_data in data["genbank"].items():
                genbank = NCBIGenbank.model_validate(gb_data)
                self._genbank_records[acc] = genbank

        self.otus = OTURegistry(manifest, data_dir)
        self.otus.validate()

        # Populate _otu_structure from manifest
        for attr_name in dir(manifest):
            attr = getattr(manifest, attr_name)
            if hasattr(attr, "refseq") and hasattr(attr, "isolates"):
                # Get taxid from loaded taxonomy
                handle = getattr(self.otus, attr_name, None)
                if handle:
                    self._otu_structure.append(
                        (handle.taxid, attr.refseq, attr.isolates)
                    )

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

    def get_incompatible_isolates(self, taxid: int) -> list[list[str]]:
        """Get incompatible isolate groups for a specific taxid.

        Args:
            taxid: NCBI taxonomy ID

        Returns:
            List of incompatible isolate groups (each group is a list of accessions)

        """
        return self._incompatible_isolates.get(taxid, [])

    def fetch_genbank_records(
        self,
        accessions: Collection[str | Accession],
    ) -> list[NCBIGenbank]:
        """Fetch mocked GenBank records.

        Accessions starting with 'MISS' or 'NC_MISS' are treated as intentional
        test cases for missing data and silently skipped (matching real NCBI behavior).
        """
        if not accessions:
            return []

        records = []

        for accession in accessions:
            if isinstance(accession, Accession):
                accession_key = accession.key
            else:
                try:
                    accession_key = Accession.from_string(accession).key
                except ValueError:
                    accession_key = accession

            if accession_key.startswith("MISS") or accession_key.startswith("NC_MISS"):
                continue

            try:
                records.append(self._genbank_records[accession_key])
            except KeyError:
                available = ", ".join(sorted(self._genbank_records.keys())[:10])
                msg = (
                    f"Accession '{accession_key}' not in mock data. "
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

    def fetch_lineage(self, taxid: int) -> Lineage:
        """Fetch mocked lineage from species down to target taxon.

        Args:
            taxid: NCBI Taxonomy ID

        Returns:
            Lineage object with complete details for each taxon

        """
        taxonomy = self.fetch_taxonomy_record(taxid)

        if taxonomy is None:
            raise ValueError(f"Could not fetch taxonomy record for taxid {taxid}")

        # If not at species level, fetch the species-level taxonomy
        if taxonomy.rank != NCBIRank.SPECIES:
            species_taxonomy = self.fetch_taxonomy_record(taxonomy.species.id)
            if species_taxonomy is None:
                raise ValueError(
                    f"Could not fetch species taxonomy for taxid {taxonomy.species.id}"
                )
        else:
            species_taxonomy = taxonomy

        # Build list of taxa from target down to species
        taxa_to_fetch = []

        # Add the target taxon if it's not species level
        if taxonomy.rank != NCBIRank.SPECIES:
            taxa_to_fetch.append(taxonomy.id)

        # Always add the species
        taxa_to_fetch.append(species_taxonomy.id)

        # Build Taxon objects
        taxon_objects = []

        for i, tid in enumerate(taxa_to_fetch):
            tax_record = self.fetch_taxonomy_record(tid)

            if tax_record is None:
                logger.warning("Could not fetch taxonomy record", taxid=tid)
                continue

            # Determine parent (None for species, otherwise next in list)
            parent_id = None if i == len(taxa_to_fetch) - 1 else taxa_to_fetch[i + 1]

            taxon = Taxon(
                id=tax_record.id,
                name=tax_record.name,
                parent=parent_id,
                rank=tax_record.rank,
                other_names=TaxonOtherNames(
                    acronym=tax_record.other_names.acronym,
                    synonyms=tax_record.other_names.equivalent_name,
                ),
            )

            taxon_objects.append(taxon)

        if not taxon_objects:
            raise ValueError(f"Could not build lineage for taxid {taxid}")

        return Lineage(taxa=taxon_objects)

    def filter_accessions(self, raw_accessions: Collection[str]) -> set[Accession]:
        """Filter raw accession list and return a set of valid Accession objects.

        Args:
            raw_accessions: Raw accession strings to filter

        Returns:
            Set of valid Accession objects

        """
        from contextlib import suppress

        valid_accessions = set()
        for accession in raw_accessions:
            with suppress(ValueError):
                valid_accessions.add(Accession.from_string(accession))
        return valid_accessions

    def fetch_accessions_by_taxid(
        self,
        taxid: int,
        sequence_min_length: int | None = None,
        sequence_max_length: int | None = None,
        modification_date_start: str | None = None,
        modification_date_end: str | None = None,
        refseq_only: bool = False,
    ) -> list[Accession]:
        """Fetch mock accessions for a given taxid from loaded OTU data.

        Args:
            taxid: NCBI Taxonomy ID
            sequence_min_length: Minimum sequence length filter (ignored in mock)
            sequence_max_length: Maximum sequence length filter (ignored in mock)
            modification_date_start: Start date filter (ignored in mock)
            modification_date_end: End date filter (ignored in mock)
            refseq_only: If True, only return RefSeq accessions

        Returns:
            List of Accession objects from the OTU structure

        """
        accession_keys = []

        for otu_taxid, refseq, isolates in self._otu_structure:
            if otu_taxid == taxid:
                if not refseq_only:
                    # Include RefSeq accessions
                    accession_keys.extend(refseq)
                    # Include all isolate accessions
                    for isolate_group in isolates:
                        accession_keys.extend(isolate_group)
                else:
                    # Only RefSeq accessions
                    accession_keys.extend(refseq)
                break

        # Convert accession keys to Accession objects by looking up in genbank records
        accessions = []
        for key in accession_keys:
            if key in self._genbank_records:
                record = self._genbank_records[key]
                accessions.append(Accession.from_string(record.accession_version))

        return accessions


# Module-level singleton for use by per-OTU modules
mock_ncbi_client = MockNCBIClient()
