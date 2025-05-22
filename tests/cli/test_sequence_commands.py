from uuid import uuid4

import click.testing

from ref_builder.cli.sequence import sequence as sequence_command_group

runner = click.testing.CliRunner()


class TestSequenceGetCommand:
    """Test `ref-builder sequence get SEQUENCE_ID`` works as expected."""

    def test_ok(self, scratch_repo):
        """Test basic command functionality."""
        otu_id = scratch_repo.get_otu_id_by_taxid(1169032)

        for sequence_id in scratch_repo.get_otu(otu_id).sequence_ids:
            result = runner.invoke(
                sequence_command_group,
                [
                    "--path",
                    str(scratch_repo.path),
                    "get",
                    str(sequence_id),
                ],
            )

            assert result.exit_code == 0

            assert str(sequence_id) in result.stdout

    def test_partial_ok(self, scratch_repo):
        """Test that a partial sequence ID can also be a valid identifier."""
        otu_id = scratch_repo.get_otu_id_by_taxid(1169032)

        for sequence_id in scratch_repo.get_otu(otu_id).sequence_ids:
            result = runner.invoke(
                sequence_command_group,
                [
                    "--path",
                    str(scratch_repo.path),
                    "get",
                    str(sequence_id)[:8],
                ],
            )

            assert result.exit_code == 0

    def test_json_ok(self, scratch_repo):
        otu_id = scratch_repo.get_otu_id_by_taxid(1169032)

        for sequence_id in scratch_repo.get_otu(otu_id).sequence_ids:
            result = runner.invoke(
                sequence_command_group,
                ["--path", str(scratch_repo.path), "get", str(sequence_id), "--json"],
            )

            assert result.exit_code == 0

    def test_empty(self, scratch_repo):
        """Test that an empty sequence identifier string exits with an error."""
        result = runner.invoke(
            sequence_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "get",
                "",
            ],
        )

        assert result.exit_code == 1

        assert "Partial ID segment must be at least 8 characters long" in result.stderr

    def test_nonexistent_id_fail(self, scratch_repo):
        """Test that a nonexistent sequence ID returns nothing in stdout."""
        result = runner.invoke(
            sequence_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "get",
                str(uuid4()),
            ],
        )

        assert result.exit_code == 1

        assert "Sequence ID could not be found" in result.stderr

        assert not result.stdout
