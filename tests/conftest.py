import shutil
from pathlib import Path

import orjson
import pytest
from polyfactory.factories.pydantic_factory import ModelFactory
from polyfactory.pytest_plugin import register_fixture
from pydantic import BaseModel, TypeAdapter
from pytest_mock import MockerFixture

from ref_builder.legacy.utils import build_legacy_otu
from ref_builder.logs import configure_logger
from ref_builder.ncbi.cache import NCBICache
from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.create import create_otu_with_taxid
from ref_builder.otu.isolate import add_genbank_isolate
from ref_builder.repo import Repo
from ref_builder.utils import DataType
from tests.fixtures.factories import (
    IsolateFactory,
    NCBIGenbankFactory,
    NCBISourceFactory,
    OTUFactory,
    OTUMinimalFactory,
    PlanFactory,
    SequenceFactory,
)

configure_logger(True)


@pytest.fixture(autouse=True)
def _seed_factories() -> None:
    """Seed the ModelFactory with a fixed seed."""
    ModelFactory.seed_random(1)


@pytest.fixture()
def uncached_ncbi_client(scratch_ncbi_cache: NCBICache) -> NCBIClient:
    """An NCBI client that ignores the cache."""
    scratch_ncbi_cache.clear()
    return NCBIClient(ignore_cache=True)


@pytest.fixture()
def files_path():
    return Path(__file__).parent / "files"


@pytest.fixture()
def empty_repo(tmp_path: Path) -> Repo:
    """An empty reference repository."""
    return Repo.new(
        DataType.GENOME,
        "Generic Viruses",
        tmp_path / "test_repo",
        "virus",
    )


@pytest.fixture()
def legacy_otu(legacy_repo_path: Path) -> dict:
    """A legacy OTU."""
    return build_legacy_otu(
        legacy_repo_path / "src" / "a" / "abaca_bunchy_top_virus",
    )


@pytest.fixture()
def legacy_repo_path(
    files_path: Path,
    tmp_path: Path,
) -> Path:
    """The path to a scratch legacy reference."""
    path = tmp_path / "legacy_repo"

    shutil.copytree(files_path / "src_v1", path / "src", dirs_exist_ok=True)

    return path


@pytest.fixture()
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


@pytest.fixture()
def scratch_ncbi_cache(mocker: MockerFixture, scratch_user_cache_path: Path):
    """A scratch NCBI cache with preloaded data."""
    mocker.patch(
        "ref_builder.ncbi.cache.user_cache_directory_path",
        scratch_user_cache_path,
    )

    return NCBICache()


@pytest.fixture()
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


@pytest.fixture()
def scratch_path(scratch_repo: Repo) -> Path:
    """The path to a scratch reference repository."""
    return scratch_repo.path


@pytest.fixture()
def scratch_repo(tmp_path: Path, scratch_event_store_data: dict[str, dict]) -> Repo:
    """A prepared scratch repository."""
    path = tmp_path / "scratch_repo"

    src_path = path / "src"
    src_path.mkdir(parents=True)

    for filename in scratch_event_store_data:
        with open(src_path / filename, "wb") as f:
            f.write(orjson.dumps(scratch_event_store_data[filename]))

    (path / ".cache").mkdir(parents=True)

    return Repo(path)


@pytest.fixture()
def scratch_repo_contents_path(files_path):
    """The path to the scratch repository's table of contents file."""
    return files_path / "src_test_contents.json"


@pytest.fixture()
def scratch_user_cache_path(files_path: Path, tmp_path: Path) -> Path:
    """A path to a user cache that contains preloaded data.

    On Linux, a user cache path would be found at ``~/.cache/ref-builder/ncbi``.

    """
    path = tmp_path / "user_cache"
    path.mkdir()

    shutil.copytree(files_path / "cache_test", path / "ncbi")

    return path


@pytest.fixture()
def scratch_event_store_data(
    pytestconfig,
    scratch_repo_contents_path: Path,
    tmp_path: Path,
) -> dict:
    """Scratch repo events. Cached in .pytest_cache."""
    scratch_src = pytestconfig.cache.get("scratch_src", None)

    if scratch_src:
        return scratch_src

    temp_scratch_repo = Repo.new(
        data_type=DataType.GENOME,
        name="src_test",
        path=tmp_path,
        organism="viruses",
    )

    with open(scratch_repo_contents_path, "rb") as f:
        toc = otu_contents_list_adapter.validate_json(f.read())

    for otu_contents in toc:
        with temp_scratch_repo.lock():
            create_otu_with_taxid(
                temp_scratch_repo,
                otu_contents.taxid,
                otu_contents.plan,
                "",
            )

        otu = temp_scratch_repo.get_otu_by_taxid(otu_contents.taxid)

        for isolate_accessions in otu_contents.contents:
            with temp_scratch_repo.lock():
                add_genbank_isolate(
                    repo=temp_scratch_repo,
                    otu=otu,
                    accessions=isolate_accessions,
                )

    scratch_src = {}

    for event_file_path in (tmp_path / "src").glob("*.json"):
        with open(event_file_path) as f:
            scratch_src[event_file_path.name] = orjson.loads(f.read())

    pytestconfig.cache.set("scratch_src", scratch_src)

    return scratch_src


class OTUContents(BaseModel):
    taxid: int
    """The Taxonomy ID of the OTU."""

    plan: list[str]
    """A list of accessions that determine the OTU's schema."""

    contents: list[list[str]]
    """A list of accessions to be added to the OTU."""


otu_contents_list_adapter = TypeAdapter(list[OTUContents])
"""Aids the serialization of a scratch repo's table of contents."""


@pytest.fixture()
def indexable_otus() -> list[OTUBuilder]:
    """A list of eight OTUs for use in Snapshotter testing."""
    otus = [
        OTUBuilder.model_validate(OTUFactory.build().model_dump()) for _ in range(8)
    ]

    # We want at least one OTU with a `None` legacy ID to test `get_id_by_legacy_id`.
    if all(otu.legacy_id is not None for otu in otus):
        otus[2].legacy_id = None

    # We want at least one OTU with an assigned legacy ID to test `get_id_by_legacy_id`.
    if all(otu.legacy_id is None for otu in otus):
        otus[2].legacy_id = "legacy"

    return otus


isolate_factory = register_fixture(IsolateFactory)
"""Fixture for a factory that generates IsolateBase instances."""

ncbi_genbank_factory = register_fixture(NCBIGenbankFactory)
"""Fixture for a factory that generates NCBIGenbank instances."""

ncbi_source_factory = register_fixture(NCBISourceFactory)
"""Fixture for a factory that generates NCBISource instances."""

otu_factory = register_fixture(OTUFactory)
"""Fixture for a factory that generates OTUBase instances."""

otu_minimal_factory = register_fixture(OTUMinimalFactory)

plan_factory = register_fixture(PlanFactory)
"""Fixture for generating Plan instances."""

sequence_factory = register_fixture(SequenceFactory)
"""Fixture for a factory that generates SequenceBase instances."""
