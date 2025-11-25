"""Build commands for ref-builder CLI."""

from pathlib import Path

import click

from ref_builder.build.fasta import build_fasta
from ref_builder.build.virtool import build_json
from ref_builder.cli.options import path_option


@click.group(name="build")
def build() -> None:
    """Build reference files from the repository."""


@build.command(name="virtool")
@click.option(
    "-V",
    "--version",
    default="",
    type=str,
    help="A version string to include in the output file",
)
@click.option(
    "-o",
    "--target-path",
    required=True,
    type=click.Path(exists=False, file_okay=True, path_type=Path),
    help="The path to write the reference.json file to",
)
@path_option
def build_virtool(path: Path, target_path: Path, version: str) -> None:
    """Build a Virtool reference.json file from the reference repository."""
    build_json(target_path, path, version)


@build.command(name="fasta")
@click.option(
    "-o",
    "--target-path",
    required=True,
    type=click.Path(exists=False, file_okay=True, path_type=Path),
    help="The path to write the FASTA file to",
)
@click.option(
    "--no-csv",
    is_flag=True,
    help="Skip generating the CSV metadata file",
)
@path_option
def build_fasta_cmd(path: Path, target_path: Path, no_csv: bool) -> None:
    """Build a FASTA file with versioned accession headers.

    By default, also generates a CSV file with sequence metadata at
    the same location with a .csv extension.
    """
    fasta_path, csv_path = build_fasta(target_path, path, not no_csv)

    click.echo(f"FASTA file written to: {fasta_path}")
    if csv_path:
        click.echo(f"CSV metadata written to: {csv_path}")
