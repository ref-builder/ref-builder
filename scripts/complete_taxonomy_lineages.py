#!/usr/bin/env python3
"""Fetch missing taxonomy records from NCBI to complete lineages.

This script ensures that for every GenBank record in the mock data,
all taxids in its taxonomy lineage are also available in MOCK_TAXONOMY_RECORDS.
"""

import sys
from pathlib import Path

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from ref_builder.ncbi.client import NCBIClient
from tests.fixtures.ncbi_mock_data import MOCK_GENBANK_RECORDS, MOCK_TAXONOMY_RECORDS


def find_missing_lineage_taxids() -> set[int]:
    """Find all taxids referenced in lineages that are missing from mock data.

    Returns:
        Set of taxids that need to be fetched from NCBI
    """
    missing_taxids = set()

    print(f"Checking lineages for {len(MOCK_GENBANK_RECORDS)} GenBank records...")

    for accession, record in MOCK_GENBANK_RECORDS.items():
        taxid = record.source.taxid

        # Check if the main taxid has taxonomy data
        if taxid not in MOCK_TAXONOMY_RECORDS:
            print(f"  Missing taxonomy for {accession} (taxid={taxid})")
            missing_taxids.add(taxid)
            continue

        # Check all lineage taxids
        taxonomy = MOCK_TAXONOMY_RECORDS[taxid]
        for lineage_entry in taxonomy.lineage:
            if lineage_entry.id not in MOCK_TAXONOMY_RECORDS:
                missing_taxids.add(lineage_entry.id)

    return missing_taxids


def fetch_and_display_taxonomy(taxids: set[int]) -> dict:
    """Fetch taxonomy records from NCBI and display them.

    Args:
        taxids: Set of taxids to fetch

    Returns:
        Dict mapping taxid to NCBITaxonomy record
    """
    if not taxids:
        print("\nNo missing taxonomy records found!")
        return {}

    print(f"\nFound {len(taxids)} missing taxonomy records:")
    for taxid in sorted(taxids):
        print(f"  - {taxid}")

    print(f"\nFetching {len(taxids)} taxonomy records from NCBI...")
    client = NCBIClient(ignore_cache=False)

    fetched = {}
    failed = []

    for i, taxid in enumerate(sorted(taxids), 1):
        print(f"  [{i}/{len(taxids)}] Fetching taxid {taxid}...", end=" ")
        try:
            record = client.fetch_taxonomy_record(taxid)
            if record:
                fetched[taxid] = record
                print(f"✓ {record.name}")
            else:
                print("✗ Not found")
                failed.append(taxid)
        except Exception as e:
            print(f"✗ Error: {e}")
            failed.append(taxid)

    if failed:
        print(f"\nFailed to fetch {len(failed)} taxids:")
        for taxid in failed:
            print(f"  - {taxid}")

    return fetched


def main():
    """Find and fetch missing taxonomy records."""
    print("=" * 70)
    print("Complete Taxonomy Lineages for Mock Data")
    print("=" * 70)

    # Find missing taxids
    missing_taxids = find_missing_lineage_taxids()

    if not missing_taxids:
        print("\n✓ All lineages are complete!")
        return

    # Fetch from NCBI
    fetched = fetch_and_display_taxonomy(missing_taxids)

    if not fetched:
        print("\n✗ No records were fetched successfully.")
        return

    print(f"\n{'=' * 70}")
    print(f"Successfully fetched {len(fetched)} taxonomy records.")
    print("Run scripts/extract_ncbi_cache_to_fixtures.py to regenerate mock data.")
    print(f"{'=' * 70}")


if __name__ == "__main__":
    main()
