import shutil
from pathlib import Path

import pytest
from polyfactory.factories.pydantic_factory import ModelFactory
from pytest_mock import MockerFixture

from ref_builder.logs import configure_logger
from ref_builder.models.repo import DataType
from ref_builder.ncbi.cache import NCBICache
from ref_builder.ncbi.client import NCBIClient, NCBIClientProtocol
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.repo import Repo
from ref_builder.services.isolate import IsolateService
from ref_builder.services.otu import OTUService
from tests.fixtures.factories import (
    IsolateFactory,
    NCBIGenbankFactory,
    NCBISourceFactory,
    NCBITaxonomyFactory,
    OTUFactory,
    OTUMinimalFactory,
    PlanFactory,
    SequenceFactory,
)
from tests.fixtures.ncbi import OTUManifest

configure_logger(True)


@pytest.fixture(autouse=True)
def _seed_factories() -> None:
    """Seed the ModelFactory with a fixed seed."""
    ModelFactory.seed_random(1)


@pytest.fixture
def mock_ncbi_client() -> NCBIClientProtocol:
    """A mock NCBI client with hardcoded test data."""
    from tests.fixtures.mock_ncbi_client import MockNCBIClient

    return MockNCBIClient(
        manifest=OTUManifest,
        data_dir=Path(__file__).parent / "fixtures" / "ncbi" / "otus",
    )


@pytest.fixture
def uncached_ncbi_client(scratch_ncbi_cache: NCBICache) -> NCBIClient:
    """An NCBI client that ignores the cache."""
    scratch_ncbi_cache.clear()
    return NCBIClient(ignore_cache=True)


@pytest.fixture
def files_path():
    return Path(__file__).parent / "files"


@pytest.fixture
def empty_repo(tmp_path: Path) -> Repo:
    """An empty reference repository."""
    return Repo.new(
        DataType.GENOME,
        "Generic Viruses",
        tmp_path / "test_repo",
        "virus",
    )


@pytest.fixture
def precached_repo(
    mocker: MockerFixture,
    scratch_user_cache_path: Path,
    tmp_path: Path,
) -> Repo:
    """A reference repository with a preloaded NCBI cache."""
    mocker.patch(
        "ref_builder.paths.user_cache_directory_path",
        return_value=scratch_user_cache_path,
    )

    return Repo.new(DataType.GENOME, "Empty", tmp_path / "precached_repo", "virus")


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
def scratch_path(scratch_repo: Repo) -> Path:
    """The path to a scratch reference repository."""
    return scratch_repo.path


@pytest.fixture
def scratch_repo(tmp_path: Path, mock_ncbi_client: NCBIClientProtocol) -> Repo:
    """A prepared scratch repository built from mock data."""
    repo = Repo.new(
        data_type=DataType.GENOME,
        name="test",
        path=tmp_path / "scratch_repo",
        organism="viruses",
    )

    otu_service = OTUService(repo, mock_ncbi_client)
    isolate_service = IsolateService(repo, mock_ncbi_client)

    with repo.lock():
        for (
            taxid,
            plan_accessions,
            isolate_accessions,
        ) in mock_ncbi_client.get_otu_structure():
            otu = otu_service.create(plan_accessions)

            if otu is None:
                raise RuntimeError(
                    f"Failed to create OTU for taxid {taxid} with accessions {plan_accessions}. "
                    f"Check mock data in tests/fixtures/ncbi/ for validation errors."
                )

            for accessions in isolate_accessions:
                isolate_service.create(otu.id, accessions)

    return repo


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
def indexable_otus() -> list[OTUBuilder]:
    """A list of eight OTUs for use in Snapshotter testing."""
    return [
        OTUBuilder.model_validate(OTUFactory.build().model_dump()) for _ in range(8)
    ]


@pytest.fixture
def isolate_factory() -> IsolateFactory:
    """Fixture for a factory that generates IsolateBase instances."""
    return IsolateFactory


@pytest.fixture
def ncbi_genbank_factory() -> NCBIGenbankFactory:
    """Fixture for a factory that generates NCBIGenbank instances."""
    return NCBIGenbankFactory


@pytest.fixture
def ncbi_source_factory() -> NCBISourceFactory:
    """Fixture for a factory that generates NCBISource instances."""
    return NCBISourceFactory


@pytest.fixture
def otu_factory() -> OTUFactory:
    """Fixture for a factory that generates OTUBase instances."""
    return OTUFactory


@pytest.fixture
def otu_minimal_factory() -> OTUMinimalFactory:
    """Fixture for a factory that generates OTUMinimal instances."""
    return OTUMinimalFactory


@pytest.fixture
def plan_factory() -> PlanFactory:
    """Fixture for generating Plan instances."""
    return PlanFactory


@pytest.fixture
def sequence_factory() -> SequenceFactory:
    """Fixture for a factory that generates SequenceBase instances."""
    return SequenceFactory


@pytest.fixture
def ncbi_taxonomy_factory() -> NCBITaxonomyFactory:
    """Fixture for a factory that generates NCBITaxonomy instances."""
    return NCBITaxonomyFactory
