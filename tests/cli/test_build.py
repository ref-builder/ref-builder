"""Tests for the build CLI commands."""

from pathlib import Path

from click.testing import CliRunner

from ref_builder.cli.build import build as build_command_group
from ref_builder.repo import Repo

runner = CliRunner()


class TestBuildVirtool:
    """Test the build virtool command."""

    def test_ok(self, scratch_repo: Repo, tmp_path: Path):
        """Test that build virtool produces a reference.json file."""
        output_path = tmp_path / "reference.json"

        result = runner.invoke(
            build_command_group,
            [
                "virtool",
                "--path",
                str(scratch_repo.path),
                "-o",
                str(output_path),
                "-V",
                "1.0.0",
            ],
        )

        assert result.exit_code == 0
        assert output_path.exists()


class TestBuildFasta:
    """Test the build fasta command."""

    def test_ok(self, scratch_repo: Repo, tmp_path: Path):
        """Test that build fasta produces FASTA and CSV files."""
        output_path = tmp_path / "output.fasta"

        result = runner.invoke(
            build_command_group,
            [
                "fasta",
                "--path",
                str(scratch_repo.path),
                "-o",
                str(output_path),
            ],
        )

        assert result.exit_code == 0
        assert output_path.exists()
        assert (tmp_path / "output.csv").exists()
        assert "FASTA file written to" in result.output
        assert "CSV metadata written to" in result.output

    def test_no_csv(self, scratch_repo: Repo, tmp_path: Path):
        """Test that --no-csv skips CSV generation."""
        output_path = tmp_path / "output.fasta"

        result = runner.invoke(
            build_command_group,
            [
                "fasta",
                "--path",
                str(scratch_repo.path),
                "-o",
                str(output_path),
                "--no-csv",
            ],
        )

        assert result.exit_code == 0
        assert output_path.exists()
        assert not (tmp_path / "output.csv").exists()
        assert "CSV metadata written to" not in result.output
