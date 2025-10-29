import pytest
from click.testing import CliRunner

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
