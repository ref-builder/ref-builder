import shutil
from pathlib import Path

import pytest
from polyfactory.factories.pydantic_factory import ModelFactory
from pytest_mock import MockerFixture

from ref_builder.logs import configure_logger
from ref_builder.ncbi.cache import NCBICache
from ref_builder.ncbi.client import NCBIClient, NCBIClientProtocol
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.repo import Repo
from ref_builder.services.cls import Services
from tests.fixtures.factories import (
    IsolateFactory,
    NCBIGenbankFactory,
    NCBISourceFactory,
    NCBITaxonomyFactory,
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
        "Generic Viruses",
        tmp_path / "test_repo",
        "virus",
    )


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


@pytest.fixture(scope="session")
def _session_scratch_repo(tmp_path_factory) -> Repo:
    """A session-scoped scratch repository built once from mock data."""
    from tests.fixtures.mock_ncbi_client import MockNCBIClient

    tmp_path = tmp_path_factory.mktemp("session_scratch")
    repo = Repo.new(
        name="test",
        path=tmp_path / "scratch_repo",
        organism="viruses",
    )

    mock_ncbi_client = MockNCBIClient(
        manifest=OTUManifest,
        data_dir=Path(__file__).parent / "fixtures" / "ncbi" / "otus",
    )

    services = Services(repo, mock_ncbi_client)

    with repo.lock():
        for (
            taxid,
            plan_accessions,
            isolate_accessions,
        ) in mock_ncbi_client.get_otu_structure():
            otu = services.otu.create(plan_accessions)

            if otu is None:
                raise RuntimeError(
                    f"Failed to create OTU for taxid {taxid} with accessions {plan_accessions}. "
                    f"Check mock data in tests/fixtures/ncbi/ for validation errors."
                )

            for accessions in isolate_accessions:
                services.isolate.create(accessions)

    return repo


@pytest.fixture
def scratch_repo(tmp_path: Path, _session_scratch_repo: Repo) -> Repo:
    """A prepared scratch repository built from mock data (copied from session fixture)."""
    repo_path = tmp_path / "scratch_repo"
    shutil.copytree(
        _session_scratch_repo.path,
        repo_path,
        ignore=shutil.ignore_patterns(".git"),
    )
    return Repo(repo_path)


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
def indexable_otus(
    empty_repo: Repo, mock_ncbi_client: NCBIClientProtocol
) -> list[OTUBuilder]:
    """A list of eight OTUs for use in Snapshotter testing."""
    services = Services(empty_repo, mock_ncbi_client)

    with empty_repo.lock():
        otus = [
            services.otu.create(OTUManifest.abaca_bunchy_top_virus.refseq),
            services.otu.create(OTUManifest.babaco_mosaic_virus.refseq),
            services.otu.create(OTUManifest.cabbage_leaf_curl_jamaica_virus.refseq),
            services.otu.create(OTUManifest.dahlia_latent_viroid.refseq),
            services.otu.create(
                OTUManifest.east_african_cassava_mosaic_cameroon_virus.refseq
            ),
            services.otu.create(OTUManifest.oat_blue_dwarf_virus.refseq),
            services.otu.create(OTUManifest.okra_leaf_curl_alphasatellite.refseq),
            services.otu.create(OTUManifest.saccharum_streak_virus.refseq),
        ]

    return [otu for otu in otus if otu is not None]


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
