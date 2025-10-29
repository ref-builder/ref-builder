import shutil
from pathlib import Path

import pytest

from ref_builder.repo import Repo
from ref_builder.services.cls import Services
from tests.fixtures.mock_ncbi_client import MockNCBIClient
from tests.fixtures.ncbi import OTUManifest


@pytest.fixture
def empty_repo(tmp_path: Path):
    """An empty reference repository."""
    return Repo.new(
        "Generic Viruses",
        tmp_path / "test_repo",
        "virus",
    )


@pytest.fixture
def scratch_path(scratch_repo: Repo) -> Path:
    """The path to a scratch reference repository."""
    return scratch_repo.path


@pytest.fixture(scope="session")
def _session_scratch_repo(tmp_path_factory):
    """A session-scoped scratch repository built once from mock data."""
    tmp_path = tmp_path_factory.mktemp("session_scratch")

    repo = Repo.new(
        name="test",
        path=tmp_path / "scratch_repo",
        organism="viruses",
    )

    mock_ncbi_client = MockNCBIClient()

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
def scratch_repo(tmp_path: Path, _session_scratch_repo: Repo):
    """A prepared scratch repository built from mock data."""
    repo_path = tmp_path / "scratch_repo"
    shutil.copytree(
        _session_scratch_repo.path,
        repo_path,
        ignore=shutil.ignore_patterns(".git"),
    )
    return Repo(repo_path)


@pytest.fixture
def indexable_otus(empty_repo: Repo, mock_ncbi_client: MockNCBIClient):
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
