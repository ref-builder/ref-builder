"""Command-line interface for reference builder."""

from pathlib import Path

import click
import structlog
from rich.table import Table

from ref_builder.cli.build import build
from ref_builder.cli.dev import dev
from ref_builder.cli.event import event
from ref_builder.cli.isolate import isolate
from ref_builder.cli.options import (
    path_option,
)
from ref_builder.cli.otu import otu
from ref_builder.cli.utils import pass_repo
from ref_builder.console import console
from ref_builder.logs import configure_logger
from ref_builder.ncbi.client import NCBIClient
from ref_builder.repo import Repo
from ref_builder.services.cls import Services

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


@entry.command(name="status")
@path_option
def repo_status(path: Path) -> None:
    """Show the status of the reference repository."""
    repo = Repo(path)

    table = Table(show_header=False)

    table.add_column("")
    table.add_column("")

    table.add_row("Name", repo.meta.name)
    table.add_row("Path", str(path.absolute()))
    table.add_row("Organism", repo.meta.organism)
    table.add_row("Events", str(repo.last_id))

    console.print(table)


@entry.command(name="init")
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
def repo_init(name: str, organism: str, path: Path) -> None:
    """Create a new reference repository.

    If a path is not provided, the repository will be created in the current directory.

    Repositories cannot be created in non-empty directories.
    """
    Repo.new(name, path, organism)


@entry.command(name="update")
@pass_repo
def repo_update(repo: Repo) -> None:
    """Update all OTUs with the latest data from NCBI."""
    services = Services(repo, NCBIClient(False))
    services.repo.update()


entry.add_command(build)
entry.add_command(otu)
entry.add_command(isolate)
entry.add_command(event)
entry.add_command(dev)
