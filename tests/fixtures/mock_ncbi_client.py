"""Mock NCBI client for testing without real API calls or file cache."""

import json
from collections.abc import Collection
from contextlib import contextmanager
from pathlib import Path

from structlog import get_logger

from ref_builder.models.accession import Accession
from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.ncbi.models import (
    NCBIGenbank,
    NCBIRank,
    NCBITaxonomy,
)
from tests.fixtures.ncbi.manifest import OTUManifest
from tests.fixtures.ncbi.models import OTURegistry

logger = get_logger("tests.mock_ncbi_client")


class MockDataNotFoundError(LookupError):
    """Raised when requested data is not available in mock.

    This indicates missing test data, not a test implementation error.
    Add the missing accession/taxid to tests/fixtures/ncbi/ modules.
    """


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
            data_dir: Optional directory containing JSON data files

        """
        self.ignore_cache = ignore_cache
        self._genbank_records: dict[str, NCBIGenbank] = {}
        self._taxonomy_records: dict[int, NCBITaxonomy] = {}
        self._otu_structure: list[tuple[int, list[str], list[list[str]]]] = []
        self._incompatible_isolates: dict[int, list[list[str]]] = {}
        self._blocked_accessions: set[str] = set()

        data_dir = Path(__file__).parent / "ncbi" / "otus"

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

        self.otus = OTURegistry(OTUManifest, data_dir)
        self.otus.validate()

        # Populate _otu_structure from manifest
        for attr_name in dir(OTUManifest):
            attr = getattr(OTUManifest, attr_name)
            if hasattr(attr, "refseq") and hasattr(attr, "isolates"):
                # Get taxid from loaded taxonomy
                handle = getattr(self.otus, attr_name, None)
                if handle:
                    self._otu_structure.append(
                        (handle.taxid, attr.refseq, attr.isolates)
                    )

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

    @contextmanager
    def blocking(self, accessions: Collection[str]):
        """Context manager to temporarily block accessions from fetch_accessions_by_taxid.

        Allows testing scenarios where sequences are added to NCBI over time.
        Blocked accessions will not appear in fetch_accessions_by_taxid results,
        but can still be fetched directly via fetch_genbank_records.

        Args:
            accessions: Versioned accessions to block (e.g., ["NC_004452.2", "NC_004452.3"])

        Example:
            with mock_client.blocking(["NC_004452.2", "NC_004452.3"]):
                # Create OTU - only .1 is discoverable
                otu = service.create(["NC_004452.1"])

            # Update OTU - now .2 and .3 are discoverable
            service.update(otu.id)

        """
        self._blocked_accessions.update(accessions)
        try:
            yield
        finally:
            self._blocked_accessions.difference_update(accessions)

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
            # Convert to string, preserving version if present
            if isinstance(accession, Accession):
                accession_str = (
                    f"{accession.key}.{accession.version}"
                    if accession.version
                    else accession.key
                )
            else:
                accession_str = str(accession)

            # Extract unversioned key for MISS check
            try:
                parsed = Accession.from_string(accession_str)
                accession_key = parsed.key
            except ValueError:
                accession_key = accession_str

            if accession_key.startswith("MISS") or accession_key.startswith("NC_MISS"):
                continue

            # Try exact match first (for versioned accessions)
            if accession_str in self._genbank_records:
                records.append(self._genbank_records[accession_str])
            # Fall back to unversioned key for backwards compatibility
            elif accession_key in self._genbank_records:
                records.append(self._genbank_records[accession_key])
            # For unversioned queries, find the highest version
            else:
                # Search for all versioned variants of this accession
                prefix = f"{accession_key}."
                versioned_matches = [
                    (key, record)
                    for key, record in self._genbank_records.items()
                    if key.startswith(prefix)
                ]

                if versioned_matches:
                    # Parse versions and find the highest
                    highest_version = 0
                    highest_record = None

                    for key, record in versioned_matches:
                        try:
                            version_str = key.split(".")[-1]
                            version = int(version_str)
                            if version > highest_version:
                                highest_version = version
                                highest_record = record
                        except (ValueError, IndexError):
                            continue

                    if highest_record:
                        records.append(highest_record)
                    else:
                        available = ", ".join(sorted(self._genbank_records.keys())[:10])
                        msg = (
                            f"Accession '{accession_str}' not in mock data. "
                            f"Add it to tests/fixtures/ncbi/ modules\n"
                            f"Available accessions (first 10): {available}..."
                        )
                        raise MockDataNotFoundError(msg) from None
                else:
                    available = ", ".join(sorted(self._genbank_records.keys())[:10])
                    msg = (
                        f"Accession '{accession_str}' not in mock data. "
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

    def fetch_descendant_taxids(self, species_taxid: int) -> list[int]:
        """Fetch all descendant taxids under a species.

        Returns only subspecific taxids (isolate, "no rank" ranks) from mock data.

        Args:
            species_taxid: A NCBI species-level Taxonomy id

        Returns:
            A list of subspecific taxids under the species

        """
        subspecific_taxids = []

        for taxid, taxonomy in self._taxonomy_records.items():
            if taxid == species_taxid:
                continue

            if taxonomy.rank in (NCBIRank.ISOLATE, NCBIRank.NO_RANK):
                if taxonomy.species.id == species_taxid:
                    subspecific_taxids.append(taxid)

        logger.debug(
            "Found subspecific descendants in mock data",
            species_taxid=species_taxid,
            count=len(subspecific_taxids),
            taxids=subspecific_taxids,
        )

        return subspecific_taxids

    def fetch_lineage(self, taxid: int) -> Lineage:
        """Fetch mocked lineage from species including all descendant taxids.

        Matches the real client by fetching ALL subspecific descendants under
        the species to ensure validation passes for any isolate.

        Args:
            taxid: NCBI Taxonomy ID

        Returns:
            Lineage object with species and all subspecific descendants

        """
        taxonomy = self.fetch_taxonomy_record(taxid)

        if taxonomy is None:
            raise ValueError(f"Could not fetch taxonomy record for taxid {taxid}")

        species_taxonomy = self.fetch_taxonomy_record(taxonomy.species.id)
        if species_taxonomy is None:
            raise ValueError(
                f"Could not fetch species taxonomy for taxid {taxonomy.species.id}"
            )

        species_taxon = Taxon(
            id=species_taxonomy.id,
            name=species_taxonomy.name,
            parent=None,
            rank=species_taxonomy.rank,
            other_names=TaxonOtherNames(
                acronym=species_taxonomy.other_names.acronym,
                synonyms=species_taxonomy.other_names.equivalent_name,
            ),
        )

        descendant_taxids = self.fetch_descendant_taxids(species_taxonomy.id)

        taxon_objects = [species_taxon]

        for descendant_taxid in descendant_taxids:
            tax_record = self.fetch_taxonomy_record(descendant_taxid)
            if tax_record is None:
                logger.warning(
                    "Could not fetch taxonomy record", taxid=descendant_taxid
                )
                continue

            parent_id = species_taxonomy.id
            for lineage_item in tax_record.lineage:
                if lineage_item.rank in ("no rank", "isolate"):
                    if lineage_item.id in descendant_taxids:
                        parent_id = lineage_item.id
                        break

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
        refseq_only: bool = False,
    ) -> list[Accession]:
        """Fetch mock accessions for a given taxid from loaded OTU data.

        Args:
            taxid: NCBI Taxonomy ID
            sequence_min_length: Minimum sequence length filter (ignored in mock)
            sequence_max_length: Maximum sequence length filter (ignored in mock)
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
        # Filter out blocked accessions
        accessions = []
        for key in accession_keys:
            if key in self._genbank_records and key not in self._blocked_accessions:
                record = self._genbank_records[key]
                accessions.append(Accession.from_string(record.accession_version))

        return accessions
