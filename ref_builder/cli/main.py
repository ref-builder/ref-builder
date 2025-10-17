"""Command-line interface for reference builder."""

from pathlib import Path

import click
import structlog
from rich.table import Table

from ref_builder.build import build_json
from ref_builder.cli.dev import dev
from ref_builder.cli.event import event
from ref_builder.cli.isolate import isolate
from ref_builder.cli.options import (
    path_option,
)
from ref_builder.cli.otu import otu
from ref_builder.console import console
from ref_builder.logs import configure_logger
from ref_builder.repo import Repo
from ref_builder.utils import DataType

logger = structlog.get_logger()


@click.group()
@click.option("--debug", is_flag=True, help="Show debug logs")
@click.option("-v", "--verbose", "verbosity", count=True)
@click.option("--no-color", is_flag=True, help="Disable colored output.")
def entry(debug: bool, verbosity: int, no_color: bool) -> None:
    """Build and maintain reference sets of pathogen genome sequences."""
    if debug:
        verbosity = 2

    configure_logger(verbosity, no_color)


@entry.command()
@path_option
def status(path: Path) -> None:
    """Show the status of the reference repository."""
    repo = Repo(path)

    table = Table(show_header=False)

    table.add_column("")
    table.add_column("")

    table.add_row("Name", repo.meta.name)
    table.add_row("Path", str(path.absolute()))
    table.add_row("Data Type", repo.meta.data_type)
    table.add_row("Organism", repo.meta.organism)
    table.add_row("Events", str(repo.last_id))
    table.add_row(
        "Valid",
        "[green]Yes[/green]" if repo.last_id == repo.head_id else "[red]No[/red]",
    )

    console.print(table)


@entry.command()
@click.option(
    "--data-type",
    help="The type of data the reference will contain (eg. genome)",
    required=True,
    type=click.Choice(DataType),
)
@click.option(
    "--name",
    help="A name for the reference",
    required=True,
    type=str,
)
@click.option(
    "--organism",
    default="",
    help="The organism the reference is for (eg. virus)",
    type=str,
)
@click.option(
    "--path",
    default=".",
    help="The path to initialize the repository at",
    type=click.Path(path_type=Path),
)
def init(data_type: DataType, name: str, organism: str, path: Path) -> None:
    """Create a new reference repository.

    If a path is not provided, the repository will be created in the current directory.

    Repositories cannot be created in non-empty directories.
    """
    Repo.new(data_type, name, path, organism)


entry.add_command(otu)
entry.add_command(isolate)
entry.add_command(event)
entry.add_command(dev)


@entry.command()
@click.option(
    "-i",
    "--indent",
    is_flag=True,
    help="Indent the output JSON file",
)
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
def build(indent: bool, path: Path, target_path: Path, version: str) -> None:
    """Build a Virtool reference.json file from the reference repository."""
    build_json(indent, target_path, path, version)
