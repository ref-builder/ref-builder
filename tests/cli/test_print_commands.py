from uuid import uuid4

from click.testing import CliRunner

from ref_builder.cli.otu import otu as otu_command_group

runner = CliRunner()


class TestPrintCommands:
    """Test otu console printing commands."""

    def test_otu_list(self, scratch_repo):
        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)] + ["list"],
        )

        assert result.exit_code == 0

    def test_otu_get_ok(self, scratch_repo):
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path)] + ["get", str(otu_.id)],
            )

            assert result.exit_code == 0

    def test_otu_get_nonexistent_fail(self, scratch_repo):
        """Test that otu get command on a nonexistent OTU id fails with expected error text."""
        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)] + ["get", str(uuid4())],
        )

        assert result.exit_code == 1

        assert "OTU not found" in result.output

    def test_otu_get_deleted_fail(self, scratch_repo):
        """Test that otu get command on a deleted OTU id fails with expected error text."""
        otu_id = scratch_repo.get_otu_id_by_taxid(438782)

        with scratch_repo.lock(), scratch_repo.use_transaction():
            scratch_repo.delete_otu(
                otu_id,
                rationale="Removing OTU for testing purposes",
                replacement_otu_id=None,
            )

        assert scratch_repo.get_otu(otu_id) is None

        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)] + ["get", str(otu_id)],
        )

        assert result.exit_code == 1

        assert "OTU has been deleted" in result.output

    def test_otu_get_json_ok(self, scratch_repo):
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path)] + ["get", str(otu_.id), "--json"],
            )

            assert result.exit_code == 0

    def test_otu_partial_ok(self, scratch_repo):
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path)] + ["get", str(otu_.id)[:8]],
            )

            assert result.exit_code == 0

    def test_otu_partial_too_short(self, scratch_repo):
        for otu_ in scratch_repo.iter_otus():
            result = runner.invoke(
                otu_command_group,
                ["--path", str(scratch_repo.path)] + ["get", str(otu_.id)[:5]],
            )

            assert result.exit_code == 1

            assert "OTU not found" in result.output

    def test_otu_get_empty(self, scratch_repo):
        result = runner.invoke(
            otu_command_group, ["--path", str(scratch_repo.path)] + ["get", ""]
        )

        assert result.exit_code == 1

        assert "OTU not found." in result.output
