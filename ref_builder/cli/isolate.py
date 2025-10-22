import sys
from pathlib import Path

import click

from ref_builder.cli.options import ignore_cache_option, path_option
from ref_builder.cli.utils import get_otu_isolate_ids_from_identifier, pass_repo
from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.console import print_isolate, print_isolate_as_json
from ref_builder.ncbi.client import NCBIClient
from ref_builder.services.isolate import IsolateService
from ref_builder.otu.modify import delete_isolate_from_otu
from ref_builder.repo import Repo, locked_repo


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
@click.option("--taxid", type=int, required=True)
@pass_repo
def isolate_create(
    repo: Repo,
    accessions_: list[str],
    taxid: int,
) -> None:
    """Create a new isolate using the given accessions."""
    otu_ = repo.get_otu_by_taxid(taxid)

    if otu_ is None:
        click.echo(f"OTU {taxid} not found.", err=True)
        sys.exit(1)

    ncbi_client = NCBIClient(False)

    isolate_service = IsolateService(repo, ncbi_client)
    isolate_service.create(otu_.id, accessions_)


@isolate.command(name="delete")
@click.argument("IDENTIFIER", type=str)
@pass_repo
def isolate_delete(repo: Repo, identifier: str) -> None:
    """Delete an isolate with a UUID corresponding to IDENTIFIER.

    IDENTIFIER is an unique isolate ID (>8 characters)
    """
    otu_id, isolate_id = get_otu_isolate_ids_from_identifier(repo, identifier)

    if (otu_id := repo.get_otu_id_by_isolate_id(isolate_id)) is None:
        click.echo("The containing OTU could not be found.", err=True)
        sys.exit(1)

    if (otu_ := repo.get_otu(otu_id)) is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    if delete_isolate_from_otu(repo, otu_, isolate_id):
        click.echo("Isolate deleted.")

    else:
        sys.exit(1)


@isolate.command(name="get")
@click.argument("IDENTIFIER", type=str)
@click.option(
    "--json",
    "json_",
    is_flag=True,
    help="Output in JSON form",
)
@pass_repo
def isolate_get(repo: Repo, identifier: str, json_: bool) -> None:
    """Get an isolate with a UUID corresponding to IDENTIFIER.

    IDENTIFIER is an unique isolate ID (>8 characters)
    """
    otu_id, isolate_id = get_otu_isolate_ids_from_identifier(repo, identifier)

    otu_ = repo.get_otu(otu_id)

    isolate_ = otu_.get_isolate(isolate_id)

    if isolate_ is None:
        click.echo("Isolate could not be found.", err=True)
        sys.exit(1)

    if json_:
        print_isolate_as_json(isolate_)
    else:
        print_isolate(isolate_, otu_.plan)
