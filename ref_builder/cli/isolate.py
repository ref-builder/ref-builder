import sys
from pathlib import Path
from uuid import UUID

import click

from ref_builder.cli.options import path_option
from ref_builder.cli.utils import pass_repo
from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.console import print_isolate
from ref_builder.ncbi.client import NCBIClient
from ref_builder.repo import Repo, locked_repo
from ref_builder.services.cls import Services


@click.group(name="isolate")
@click.pass_context
@path_option
def isolate(ctx: click.Context, path: Path) -> None:
    """Manage isolates."""
    ctx.obj = ctx.with_resource(locked_repo(path))


@isolate.command(name="create")
@click.argument(
    "accessions_",
    callback=validate_no_duplicate_accessions,
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@pass_repo
def isolate_create(
    repo: Repo,
    accessions_: list[str],
) -> None:
    """Create a new isolate using the given accessions."""
    services = Services(repo, NCBIClient(False))
    services.isolate.create(accessions_)


@isolate.command(name="delete")
@click.argument("ISOLATE_ID", type=UUID)
@click.option("-m", "--message", required=True, help="Reason for deletion")
@pass_repo
def isolate_delete(repo: Repo, isolate_id: UUID, message: str) -> None:
    """Delete an isolate by UUID.

    ISOLATE_ID is the isolate UUID.
    """
    services = Services(repo, NCBIClient(False))

    if services.isolate.delete(isolate_id, message):
        click.echo("Isolate deleted.")
    else:
        sys.exit(1)


@isolate.command(name="get")
@click.argument("ISOLATE_ID", type=UUID)
@pass_repo
def isolate_get(repo: Repo, isolate_id: UUID) -> None:
    """Get an isolate by UUID.

    ISOLATE_ID is the isolate UUID.
    """
    isolate_ = repo.get_isolate(isolate_id)

    if isolate_ is None:
        click.echo("Isolate could not be found.", err=True)
        sys.exit(1)

    otu_id = repo.get_otu_id_by_isolate_id(isolate_id)
    otu_ = repo.get_otu(otu_id)
    print_isolate(isolate_, otu_.plan)
