"""Minimal proof of concept test for manifest-driven mock NCBI architecture."""

from pathlib import Path

from tests.fixtures.ncbi.manifest import OTUManifest
from tests.fixtures.ncbi.models import OTURegistry


def test_otu_registry_loads():
    """Verify OTU registry loads from manifest and JSON data."""
    registry = OTURegistry(OTUManifest, Path("tests/fixtures/ncbi/otus"))

    assert registry.babaco_mosaic_virus.taxid == 3158379
    assert registry.babaco_mosaic_virus.name == "babaco_mosaic_virus"
