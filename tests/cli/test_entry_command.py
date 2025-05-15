import shutil

import pytest
from click.testing import CliRunner

from ref_builder.cli.main import entry
from ref_builder.repo import Repo

runner = CliRunner()


@pytest.fixture()
def bad_directory_path(tmp_path):
    bad_directory_path = tmp_path / "bad_directory"

    bad_directory_path.mkdir(exist_ok=True)

    yield bad_directory_path

    shutil.rmtree(bad_directory_path)


@pytest.fixture()
def dummy_reference_output_path(tmp_path):
    reference_output_path = tmp_path / "dummy_reference.json"

    yield reference_output_path

    reference_output_path.unlink(missing_ok=True)


def test_ok():
    """Test that interface loads up as expected"""
    result = runner.invoke(entry, ["--help"])

    assert result.exit_code == 0

    assert (
        "Build and maintain reference sets of pathogen genome sequences."
        in result.output
    )


def test_no_color_ok():
    result = runner.invoke(entry, ["--no-color", "--help"])

    assert result.exit_code == 0

    assert (
        "Build and maintain reference sets of pathogen genome sequences."
        in result.output
    )


class TestBuildCommand:
    """Test that ``ref-builder build`` works as expected."""

    def test_ok(self, scratch_repo: Repo, dummy_reference_output_path):
        result = runner.invoke(
            entry,
            [
                "build",
                "--path",
                str(scratch_repo.path),
                "--version",
                "0.1.0",
                "--target-path",
                str(dummy_reference_output_path)
            ]
        )

        assert result.exit_code == 0

        assert dummy_reference_output_path.exists()

    def test_wrong_directory_fail(self, bad_directory_path, dummy_reference_output_path):
        result = runner.invoke(
            entry,
            [
                "build",
                "--path",
                str(bad_directory_path),
                "--version",
                "0.1.0",
                "--target-path",
                str(dummy_reference_output_path)
            ]
        )

        assert result.exit_code == 1

        assert f'Given path "{bad_directory_path}" is not a reference repository' in result.output


class TestStatusCommand:
    """Test that ``ref-builder status`` works as expected."""

    def test_ok(self, scratch_repo: Repo):
        result = runner.invoke(
            entry,
            ["status", "--path", str(scratch_repo.path)]
        )

        assert result.exit_code == 0

        assert "Data Type" in result.output

    def test_wrong_directory_fail(self, bad_directory_path):
        result = runner.invoke(
            entry,
            ["status", "--path", str(bad_directory_path)]
        )

        assert result.exit_code == 1

        assert f'Given path "{bad_directory_path}" is not a reference repository' in result.output
