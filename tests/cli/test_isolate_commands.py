from click.testing import CliRunner

from ref_builder.cli.isolate import isolate as isolate_command_group
from ref_builder.cli.otu import otu as otu_command_group
from ref_builder.models.isolate import IsolateName, IsolateNameType
from ref_builder.repo import Repo

runner = CliRunner()


class TestIsolateCreateCommand:
    """Test `ref-builder isolate create`` works as expected."""

    def test_ok(self, empty_repo: Repo):
        """Test basic command functionality."""
        first_isolate_accessions = [
            "EF546808",
            "EF546809",
            "EF546810",
            "EF546811",
            "EF546812",
            "EF546813",
        ]

        result = runner.invoke(
            otu_command_group,
            ["--path", str(empty_repo.path), "create", *first_isolate_accessions],
        )

        assert result.exit_code == 0

        second_isolate_accessions = [
            "EF546802",
            "EF546803",
            "EF546804",
            "EF546805",
            "EF546806",
            "EF546807",
        ]

        result = runner.invoke(
            isolate_command_group,
            [
                "--path",
                str(empty_repo.path),
                "create",
                *second_isolate_accessions,
            ],
        )

        assert result.exit_code == 0
        assert "Isolate created" in result.output

    def test_duplicate_accessions(self, scratch_repo: Repo):
        """Test that an error is raised when duplicate accessions are provided."""
        result = runner.invoke(
            isolate_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "create",
                "DQ178610",
                "DQ178610",
            ],
        )

        assert result.exit_code == 2

        assert "Duplicate accessions are not allowed." in result.output


class TestIsolateGetCommand:
    """Test `ref-builder isolate get ISOLATE_ID`` works as expected."""

    def test_ok(self, scratch_repo: Repo, mock_ncbi_client):
        """Test basic command functionality."""
        otu_id = scratch_repo.get_otu_id_by_taxid(
            mock_ncbi_client.otus.wasabi_mottle_virus.taxid
        )

        assert otu_id

        otu = scratch_repo.get_otu(otu_id)

        assert otu

        for isolate_id in otu.isolate_ids:
            print({"isolate_id": isolate_id})
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

    def test_missing_id(self, scratch_repo: Repo):
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

        assert result.exit_code == 2


class TestIsolateDeleteCommand:
    """Test `ref-builder isolate delete ISOLATE_ID`` works as expected."""

    def test_ok(self, scratch_repo: Repo, mock_ncbi_client):
        """Test basic command functionality."""
        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.wasabi_mottle_virus.taxid
        )

        assert otu

        isolate_id = otu.get_isolate_id_by_name(
            IsolateName(type=IsolateNameType.ISOLATE, value="WMoV-6.3"),
        )

        assert isolate_id in otu.isolate_ids

        result = runner.invoke(
            isolate_command_group,
            ["--path", str(scratch_repo.path), "delete", str(isolate_id)],
        )

        assert result.exit_code == 0
        assert "Isolate deleted" in result.output

        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.wasabi_mottle_virus.taxid
        )

        assert otu
        assert isolate_id not in otu.isolate_ids

    def test_missing_id(self, scratch_repo: Repo):
        """Test that an empty isolate identifier string exits with an error."""
        result = runner.invoke(
            isolate_command_group,
            ["--path", str(scratch_repo.path), "delete", ""],
        )

        assert result.exit_code == 2
