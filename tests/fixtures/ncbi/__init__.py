"""NCBI mock data organized by OTU.

This package contains per-OTU modules that register mock NCBI data using the
mock_ncbi_client singleton. Each otu_*.py module registers an OTU with its
associated GenBank and taxonomy records.
"""

import importlib
from pathlib import Path

from tests.fixtures.mock_ncbi_client import mock_ncbi_client

# Auto-import all otu_*.py modules to trigger registration
_ncbi_dir = Path(__file__).parent
for module_file in sorted(_ncbi_dir.glob("otu_*.py")):
    importlib.import_module(f"tests.fixtures.ncbi.{module_file.stem}")

__all__ = ["mock_ncbi_client"]
