import pytest
from click.testing import CliRunner
from pytest_mock import MockerFixture

from ref_builder.cli.otu import otu as otu_command_group
from ref_builder.repo import Repo

runner = CliRunner()


class TestCreateOTU:
    """Test the behaviour of ``ref-builder otu create``."""

    @pytest.mark.parametrize(
        ("taxid", "accessions"),
        [(3429802, ["NC_020160"]), (3426695, ["DQ178610", "DQ178611"])],
    )
    def test_ok(
        self,
        accessions: list[str],
        empty_repo: Repo,
        taxid: int,
    ):
        """Test that an OTU can be created without the --taxid option."""
        result = runner.invoke(
            otu_command_group,
            ["--path", str(empty_repo.path), "create", *accessions],
        )

        assert result.exit_code == 0

        otus = list(Repo(empty_repo.path).iter_otus())

        assert len(otus) == 1
        assert otus[0].taxid == taxid


class TestCreateOTUWithDuplicateTaxonomy:
    """Test the behaviour of ``ref-builder otu create`` with duplicate taxonomy IDs."""

    def test_without_flag_shows_error(self, empty_repo: Repo):
        """Test that creating OTU with duplicate taxid without flag shows error."""
        # Create first OTU
        result = runner.invoke(
            otu_command_group,
            ["--path", str(empty_repo.path), "create", "NC_001367"],
        )
        assert result.exit_code == 0

        # Try to create another OTU with same taxonomy ID without flag
        result = runner.invoke(
            otu_command_group,
            ["--path", str(empty_repo.path), "create", "V01408"],
        )

        assert result.exit_code == 1
        assert "OTU already exists for taxonomy ID 3432891" in result.output
        assert "Use -i/--create-isolate to create an isolate instead" in result.output

    def test_with_flag_creates_isolate(self, empty_repo: Repo):
        """Test that creating OTU with duplicate taxid and flag creates isolate."""
        # Create first OTU
        result = runner.invoke(
            otu_command_group,
            ["--path", str(empty_repo.path), "create", "NC_001367"],
        )
        assert result.exit_code == 0

        repo = Repo(empty_repo.path)
        otus = list(repo.iter_otus())
        assert len(otus) == 1
        initial_isolate_count = len(otus[0].isolate_ids)

        # Try to create another OTU with same taxonomy ID using -i flag
        result = runner.invoke(
            otu_command_group,
            ["--path", str(empty_repo.path), "create", "-i", "OQ953825"],
        )

        assert result.exit_code == 0
        assert "Created isolate in existing OTU (taxid: 3432891)" in result.output

        # Verify isolate was added to existing OTU
        repo = Repo(empty_repo.path)
        otus = list(repo.iter_otus())
        assert len(otus) == 1  # Still only one OTU
        assert len(otus[0].isolate_ids) == initial_isolate_count + 1  # One more isolate

    def test_with_flag_but_isolate_creation_fails(
        self, empty_repo: Repo, mocker: MockerFixture
    ):
        """Test proper error handling when isolate creation fails even with flag."""
        # Create first OTU
        result = runner.invoke(
            otu_command_group,
            ["--path", str(empty_repo.path), "create", "NC_001367"],
        )
        assert result.exit_code == 0

        # Mock isolate.create to return None (failure)
        mocker.patch(
            "ref_builder.services.isolate.IsolateService.create",
            return_value=None,
        )

        # Try to create with -i flag but isolate creation fails
        result = runner.invoke(
            otu_command_group,
            ["--path", str(empty_repo.path), "create", "-i", "V01408"],
        )

        assert result.exit_code == 1
        assert "Failed to create isolate in existing OTU" in result.output
