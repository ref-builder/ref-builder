from uuid import UUID, uuid4

import click.testing

from ref_builder.cli.isolate import isolate as isolate_command_group
from ref_builder.cli.main import entry as top_command_group
from ref_builder.cli.otu import otu as otu_command_group
from ref_builder.utils import IsolateName, IsolateNameType

runner = click.testing.CliRunner()


def run_isolate_with_debug_logs(
    repo_path: Path,
    args: list,
    env_mapping: dict[str, str] | None = None,
) -> click.testing.Result:
    """Invoke OTU command with debug logs enabled."""
    env = env_mapping if env_mapping is not None else {}

    return runner.invoke(
        top_command_group,
        [
            "--debug",
            "isolate",
            "--path",
            str(repo_path),
            *args,
        ],
        env=env,
    )


class TestIsolateCreateCommand:
    """Test `ref-builder isolate create`` works as expected."""

    def test_ok(self, precached_repo):
        """Test basic command functionality."""
        path_option_list = ["--path", str(precached_repo.path)]

        taxid = 1169032

        rep_isolate_accessions = ["MF062136", "MF062137", "MF062138"]

        result = runner.invoke(
            otu_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + rep_isolate_accessions,
        )

        assert result.exit_code == 0

        second_isolate_accessions = ["MF062125", "MF062126", "MF062127"]

        result = runner.invoke(
            isolate_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + second_isolate_accessions,
        )

        assert result.exit_code == 0

        otu_after = precached_repo.get_otu_by_taxid(taxid)

        sequence_after = otu_after.get_sequence_by_accession(
            second_isolate_accessions[0]
        )

        second_isolate_id = list(
            otu_after.get_isolate_ids_containing_sequence_id(sequence_after.id)
        )[0]

        assert UUID(result.stdout.strip("\n")) == second_isolate_id

        second_isolate = precached_repo.get_isolate(second_isolate_id)

        assert second_isolate is not None

        assert second_isolate.accessions == set(second_isolate_accessions)

    def test_overwrite_name_option_ok(self, precached_repo):
        """Test that --name option exits smoothly"""
        path_option_list = ["--path", str(precached_repo.path)]

        taxid = 345184

        result = runner.invoke(
            otu_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + ["DQ178610", "DQ178611"],
        )

        assert result.exit_code == 0

        result = runner.invoke(
            isolate_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + ["--name", "isolate", "dummy"]
            + ["DQ178613", "DQ178614"],
        )

        assert result.exit_code == 0

        otu_after = precached_repo.get_otu_by_taxid(taxid)

        sequence_after = otu_after.get_sequence_by_accession("DQ178613")

        second_isolate_id = list(
            otu_after.get_isolate_ids_containing_sequence_id(sequence_after.id)
        )[0]

        assert UUID(result.stdout.strip("\n")) == second_isolate_id

        second_isolate = precached_repo.get_isolate(second_isolate_id)

        assert second_isolate is not None

        assert second_isolate.accessions == {"DQ178613", "DQ178614"}

    def test_unnamed_option_ok(self, precached_repo):
        """Test that --unnamed option exits smoothly."""
        path_option_list = ["--path", str(precached_repo.path)]

        taxid = 345184

        result = runner.invoke(
            otu_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + ["DQ178610", "DQ178611"],
        )

        assert result.exit_code == 0

        result = runner.invoke(
            isolate_command_group,
            path_option_list
            + ["create", "--taxid", str(taxid)]
            + ["--unnamed"]
            + ["DQ178613", "DQ178614"],
        )

        assert result.exit_code == 0

        otu_after = precached_repo.get_otu_by_taxid(taxid)

        sequence_after = otu_after.get_sequence_by_accession("DQ178613")

        second_isolate_id = list(
            otu_after.get_isolate_ids_containing_sequence_id(sequence_after.id)
        )[0]

        assert UUID(result.stdout.strip("\n")) == second_isolate_id

        second_isolate = precached_repo.get_isolate(second_isolate_id)

        assert second_isolate is not None

        assert second_isolate.accessions == {"DQ178613", "DQ178614"}

    def test_duplicate_accessions_error(self, scratch_repo):
        """Test that an error is raised when duplicate accessions are provided."""
        result = runner.invoke(
            isolate_command_group,
            ["--path", str(scratch_repo.path)]
            + ["create", "--taxid", str(345184)]
            + ["DQ178610", "DQ178610"],
        )

        assert result.exit_code == 2

        assert "Duplicate accessions are not allowed." in result.output


class TestIsolateGetCommand:
    """Test `ref-builder isolate get ISOLATE_ID`` works as expected."""

    def test_ok(self, scratch_repo):
        """Test basic command functionality."""
        otu_id = scratch_repo.get_otu_id_by_taxid(1169032)

        for isolate_id in scratch_repo.get_otu(otu_id).isolate_ids:
            result = runner.invoke(
                isolate_command_group,
                [
                    "--path",
                    str(scratch_repo.path),
                    "get",
                    str(isolate_id),
                ],
            )

            assert result.exit_code == 0

    def test_partial_ok(self, scratch_repo):
        """Test that a partial isolate ID can also be a valid identifier."""
        otu_id = scratch_repo.get_otu_id_by_taxid(1169032)

        for isolate_id in scratch_repo.get_otu(otu_id).isolate_ids:
            result = runner.invoke(
                isolate_command_group,
                [
                    "--path",
                    str(scratch_repo.path),
                    "get",
                    str(isolate_id)[:8],
                ],
            )

            assert result.exit_code == 0

    def test_json_ok(self, scratch_repo):
        otu_id = scratch_repo.get_otu_id_by_taxid(1169032)

        for isolate_id in scratch_repo.get_otu(otu_id).isolate_ids:
            result = runner.invoke(
                isolate_command_group,
                ["--path", str(scratch_repo.path), "get", str(isolate_id), "--json"],
            )

            assert result.exit_code == 0

    def test_empty(self, scratch_repo):
        """Test that an empty isolate identifier string exits with an error."""
        result = runner.invoke(
            isolate_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "get",
                "",
            ],
        )

        assert result.exit_code == 1

        assert "Partial ID segment must be at least 8 characters long" in result.stderr

    def test_nonexistent_isolate_id_fail(self, scratch_repo):
        """Test that a nonexistent isolate ID returns nothing in stdout."""
        result = runner.invoke(
            isolate_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "get",
                str(uuid4()),
            ],
        )

        assert result.exit_code == 1

        assert "Isolate ID could not be found" in result.stderr

        assert not result.stdout


class TestIsolateDeleteCommand:
    """Test `ref-builder isolate delete ISOLATE_ID`` works as expected."""

    def test_ok(self, scratch_repo):
        """Test basic command functionality."""
        otu = scratch_repo.get_otu_by_taxid(1169032)

        isolate_id = otu.get_isolate_id_by_name(
            IsolateName(type=IsolateNameType.ISOLATE, value="WMoV-6.3"),
        )

        assert isolate_id in otu.isolate_ids

        result = run_isolate_with_debug_logs(
            scratch_repo,
            ["--path", str(scratch_repo.path), "delete", str(isolate_id)],
        )

        assert result.exit_code == 0

        assert f"Isolate {isolate_id} deleted" in result.stderr

        assert isolate_id not in scratch_repo.get_otu_by_taxid(1169032).isolate_ids

    def test_with_partial_id_ok(self, scratch_repo):
        """Test that a partial isolate ID can also be a valid identifier."""
        otu = scratch_repo.get_otu_by_taxid(1169032)

        isolate_id = otu.get_isolate_id_by_name(
            IsolateName(type=IsolateNameType.ISOLATE, value="WMoV-6.3"),
        )

        assert isolate_id in otu.isolate_ids

        result = runner.invoke(
            isolate_command_group,
            ["--path", str(scratch_repo.path), "delete", str(isolate_id)[:8]],
        )

        assert result.exit_code == 0

        assert f"Isolate {isolate_id} deleted" in result.stderr

        assert isolate_id not in scratch_repo.get_otu_by_taxid(1169032).isolate_ids

    def test_empty(self, scratch_repo):
        """Test that an empty isolate identifier string exits with an error."""
        result = runner.invoke(
            isolate_command_group,
            ["--path", str(scratch_repo.path), "delete", ""],
        )

        assert result.exit_code == 1

    def test_delete_representative_isolate_fail(self, scratch_repo):
        otu = scratch_repo.get_otu_by_taxid(1169032)

        result = runner.invoke(
            isolate_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "delete",
                str(otu.representative_isolate),
            ],
        )

        assert result.exit_code == 1

        assert (
            f"Isolate cannot be deleted, due to being the representative isolate of OTU {otu.id}"
            in result.stderr
        )


class TestIsolateGetOTUCommand:
    """Test that ``ref-builder isolate get-otu-id ISOLATE_ID`` works as expected."""

    def test_ok(self, scratch_repo):
        otu = scratch_repo.get_otu_by_taxid(1169032)

        result = runner.invoke(
            isolate_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "get-otu-id",
                str(otu.representative_isolate),
            ],
        )

        assert result.exit_code == 0

        assert UUID(result.stdout.strip("\n")) == otu.id

    def test_partial_id_ok(self, scratch_repo):
        """Test with a truncated isolate id as identifier."""
        otu = scratch_repo.get_otu_by_taxid(1169032)

        result = runner.invoke(
            isolate_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "get-otu-id",
                str(otu.representative_isolate)[0:8],
            ],
        )

        assert result.exit_code == 0

        assert UUID(result.stdout.strip("\n")) == otu.id

    def test_nonexistent_isolate_id_fail(self, scratch_repo):
        """Test that a nonexistent isolate ID returns nothing in stdout."""
        result = runner.invoke(
            isolate_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "get-otu-id",
                str(uuid4()),
            ],
        )

        assert result.exit_code == 1

        assert "Isolate ID could not be found" in result.stderr

        assert not result.stdout
