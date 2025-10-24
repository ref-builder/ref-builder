from click.testing import CliRunner

from ref_builder.cli.otu import otu as otu_command_group
from ref_builder.repo import Repo

runner = CliRunner()


def test_otu_list(scratch_repo: Repo):
    result = runner.invoke(
        otu_command_group,
        ["--path", str(scratch_repo.path), "list"],
    )

    assert result.exit_code == 0


class TestOTUGet:
    """Test otu console printing commands."""

    def test_plain(self, scratch_repo: Repo):
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path), "get", str(otu_.id)],
            )

            if result.exit_code != 0:
                print(result.output)

            assert result.exit_code == 0

    def test_json(self, scratch_repo: Repo):
        """Test otu console printing commands with JSON output."""
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path), "get", str(otu_.id), "--json"],
            )

            if result.exit_code != 0:
                print(result.output)

            assert result.exit_code == 0

    def test_partial_id(self, scratch_repo: Repo):
        """Test that a partial ID can be used to retrieve OTUs."""
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path), "get", str(otu_.id)[:8]],
            )

            if result.exit_code != 0:
                print(result.output)

            assert result.exit_code == 0

    def test_otu_partial_too_short(self, scratch_repo: Repo):
        """Test that a partial ID that is too short raises an error."""
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path), "get", str(otu_.id)[:5]],
            )

            if result.exit_code != 1:
                print(result.output)

            assert result.exit_code == 1
            assert "OTU not found" in result.output

    def test_empty(self, scratch_repo):
        """Test that an empty ID string raises an error."""
        result = runner.invoke(
            otu_command_group, ["--path", str(scratch_repo.path), "get", ""]
        )

        if result.exit_code != 1:
            print(result.output)

        assert result.exit_code == 1
        assert "OTU not found." in result.output
