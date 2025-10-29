import shutil
from pathlib import Path

import pytest
from polyfactory.factories.pydantic_factory import ModelFactory
from pytest_mock import MockerFixture

from ref_builder.logs import configure_logger
from ref_builder.ncbi.cache import NCBICache
from ref_builder.ncbi.client import NCBIClient, NCBIClientProtocol
from tests.fixtures.factories import (
    NCBIGenbankFactory,
    NCBISourceFactory,
    NCBITaxonomyFactory,
)

configure_logger(True)

# Register fixtures from other modules
pytest_plugins = ["tests.fixtures.repo"]


@pytest.fixture(autouse=True)
def _seed_factories() -> None:
    """Seed the ModelFactory with a fixed seed."""
    ModelFactory.seed_random(1)


@pytest.fixture
def mock_ncbi_client() -> NCBIClientProtocol:
    """A mock NCBI client with hardcoded test data."""
    from tests.fixtures.mock_ncbi_client import MockNCBIClient

    return MockNCBIClient()


@pytest.fixture
def uncached_ncbi_client(scratch_ncbi_cache: NCBICache) -> NCBIClient:
    """An NCBI client that ignores the cache."""
    scratch_ncbi_cache.clear()
    return NCBIClient(ignore_cache=True)


@pytest.fixture
def files_path():
    return Path(__file__).parent / "files"


@pytest.fixture
def scratch_ncbi_cache(mocker: MockerFixture, scratch_user_cache_path: Path):
    """A scratch NCBI cache with preloaded data."""
    mocker.patch(
        "ref_builder.ncbi.cache.user_cache_directory_path",
        scratch_user_cache_path,
    )

    return NCBICache()


@pytest.fixture
def scratch_ncbi_client(
    mocker: MockerFixture,
    scratch_user_cache_path: Path,
) -> NCBIClient:
    """A scratch NCBI client with a preloaded cache."""
    mocker.patch(
        "ref_builder.ncbi.cache.user_cache_directory_path",
        scratch_user_cache_path,
    )

    return NCBIClient(ignore_cache=False)


@pytest.fixture
def scratch_user_cache_path(files_path: Path, tmp_path: Path) -> Path:
    """A path to a user cache that contains preloaded data.

    On Linux, a user cache path would be found at ``~/.cache/ref-builder/ncbi``.

    """
    path = tmp_path / "user_cache"
    path.mkdir()

    shutil.copytree(files_path / "cache_test", path / "ncbi")

    return path


@pytest.fixture
def ncbi_genbank_factory() -> NCBIGenbankFactory:
    """Fixture for a factory that generates NCBIGenbank instances."""
    return NCBIGenbankFactory


@pytest.fixture
def ncbi_source_factory() -> NCBISourceFactory:
    """Fixture for a factory that generates NCBISource instances."""
    return NCBISourceFactory


@pytest.fixture
def ncbi_taxonomy_factory() -> NCBITaxonomyFactory:
    """Fixture for a factory that generates NCBITaxonomy instances."""
    return NCBITaxonomyFactory
