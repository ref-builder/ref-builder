from click.testing import CliRunner

from ref_builder.cli.main import entry

runner = CliRunner()


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
