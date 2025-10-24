import sys
from pathlib import Path

import click

from ref_builder.cli.options import path_option
from ref_builder.cli.utils import pass_repo
from ref_builder.console import (
    print_event,
    print_event_as_json,
    print_event_list,
    print_otu_event_log,
)
from ref_builder.ncbi.client import NCBIClient
from ref_builder.repo import Repo, locked_repo
from ref_builder.services.cls import Services


@click.group(name="event")
@path_option
@click.pass_context
def event(ctx: click.Context, path: Path) -> None:
    """Read events."""
    ctx.obj = ctx.with_resource(locked_repo(path))


@event.command("list")
@click.option(
    "--otu",
    "otu_identifier_",
    type=str,
    help="An identifier for an OTU, either a UUID or NCBI Taxonomy ID.",
)
@pass_repo
def event_list(repo: Repo, otu_identifier_: str | None) -> None:
    """List the events in this repo."""
    if otu_identifier_ is None:
        print_event_list(repo.iter_event_metadata())

    else:
        services = Services(repo, NCBIClient(False))
        otu = services.otu.get_otu(otu_identifier_)

        if otu is None:
            click.echo("OTU not found.", err=True)
            sys.exit(1)

        print_otu_event_log(list(repo.iter_otu_events(otu.id)))


@event.command(name="get")
@click.argument("IDENTIFIER", type=int)
@click.option(
    "--json",
    "json_",
    is_flag=True,
    help="Output in JSON form",
)
@pass_repo
def event_get(repo: Repo, identifier: int, json_: bool) -> None:
    """Get and view an event."""
    if (event_ := repo.get_event(identifier)) is None:
        click.echo("Event not found.", err=True)
        sys.exit(1)

    if json_:
        print_event_as_json(event_)
    else:
        print_event(event_)
