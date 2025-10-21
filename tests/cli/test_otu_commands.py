import pytest
from click.testing import CliRunner

from ref_builder.cli.otu import otu as otu_command_group
from ref_builder.repo import Repo

runner = CliRunner()


class TestCreateOTUCommands:
    """Test the behaviour of ``ref-builder otu create``."""

    @pytest.mark.parametrize(
        ("taxid", "accessions"),
        [(3429802, ["NC_020160"]), (3426695, ["DQ178610", "DQ178611"])],
    )
    def test_without_taxid(
        self,
        accessions: list[str],
        precached_repo: Repo,
        taxid: int,
    ):
        """Test that an OTU can be created without the --taxid option."""
        result = runner.invoke(
            otu_command_group,
            ["--path", str(precached_repo.path), "create", *accessions],
        )

        assert result.exit_code == 0

        otus = list(Repo(precached_repo.path).iter_otus())

        assert len(otus) == 1
        assert otus[0].taxid == taxid


class TestUpdateOTUCommand:
    """Test that the ``ref-builder otu update`` command works as planned."""

    def test_ok(self, empty_repo: Repo):
        """Test comprehensive update: promote, upgrade, and add new isolates."""
        path_option = ["--path", str(empty_repo.path)]

        original_accessions = {"MF062125", "MF062126", "MF062127"}

        result = runner.invoke(
            otu_command_group,
            [
                *path_option,
                "create",
                *list(original_accessions),
            ],
        )

        assert result.exit_code == 0

        taxid = int(result.output.split("TAXID")[1].split()[0])

        otu_before = empty_repo.get_otu_by_taxid(taxid)

        assert otu_before
        assert otu_before.accessions == original_accessions

        result = runner.invoke(
            otu_command_group,
            [*path_option, "update", str(taxid)],
        )

        assert result.exit_code == 0

        repo_after = Repo(empty_repo.path)

        otu_after = repo_after.get_otu(otu_before.id)

        assert otu_after

        # Verify promotion occurred: RefSeq accessions present
        refseq_accessions = {"NC_055390", "NC_055391", "NC_055392"}
        assert refseq_accessions.issubset(otu_after.accessions)

        # Verify original accessions were excluded after promotion
        assert original_accessions.issubset(otu_after.excluded_accessions)

        # Verify new isolates were added (more accessions than just RefSeq)
        assert len(otu_after.accessions) > len(refseq_accessions)


class TestExcludeAccessionsCommand:
    """Test that ``ref-builder otu exclude-accessions`` behaves as expected."""

    def test_ok(self, scratch_repo: Repo):
        """Test that command lists out new excluded accessions."""
        result = runner.invoke(
            otu_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "exclude-accessions",
                "3426695",  # Begomovirus brassicajamaicaense (species-level taxid for cabbage leaf curl jamaica virus)
                "DQ178608",
                "DQ178609",
            ],
        )

        assert result.exit_code == 0
        assert "Added accessions to excluded accession list" in result.output
        assert "['DQ178608', 'DQ178609']" in result.output

    def test_redundant(self, scratch_repo: Repo, mock_ncbi_client):
        """Test that the command informs the user there will be not net change."""
        result = runner.invoke(
            otu_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "exclude-accessions",
                str(mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid),
                "DQ178608",
                "DQ178609",
            ],
        )

        assert result.exit_code == 0
        assert "['DQ178608', 'DQ178609']" in result.output

        result = runner.invoke(
            otu_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "exclude-accessions",
                str(mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid),
                "DQ178608",
                "DQ178609",
            ],
        )

        assert result.exit_code == 0
        assert "Excluded accession list already up to date" in result.output


class TestAllowAccessionsCommand:
    """Test that ``ref-builder otu allow-accessions`` behaves as expected."""

    def test_ok(self, scratch_repo: Repo, mock_ncbi_client):
        """Test that command lists out newly allowable accessions"""
        taxid = mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid

        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)]
            + ["exclude-accessions", str(taxid)]
            + ["DQ178608", "DQ178609"],
        )

        assert result.exit_code == 0

        result = runner.invoke(
            otu_command_group,
            ["--path", str(scratch_repo.path)]
            + ["allow-accessions", str(taxid)]
            + ["DQ178608"],
        )

        assert result.exit_code == 0
        assert "Removed accessions from excluded accession list" in result.output
        assert "['DQ178608']" in result.output
        assert "Updated excluded accession list" in result.output
        assert "['DQ178609']" in result.output

    def test_redundant(self, scratch_repo: Repo, mock_ncbi_client):
        """Test that the command informs the user when excluded accessions are already up to date."""
        result = runner.invoke(
            otu_command_group,
            [
                "--path",
                str(scratch_repo.path),
                "allow-accessions",
                str(mock_ncbi_client.otus.cabbage_leaf_curl_jamaica_virus.taxid),
                "DQ178612",
                "DQ178613",
            ],
        )

        assert result.exit_code == 0
        assert "Excluded accession list already up to date" in result.output
