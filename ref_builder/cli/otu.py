import sys
from pathlib import Path

import click
import structlog

from ref_builder.cli.options import ignore_cache_option, path_option
from ref_builder.cli.utils import pass_repo
from ref_builder.cli.validate import validate_no_duplicate_accessions
from ref_builder.console import (
    print_otu,
    print_otu_as_json,
    print_otu_event_log,
    print_otu_list,
)
from ref_builder.ncbi.client import NCBIClient
from ref_builder.repo import Repo, locked_repo
from ref_builder.services.cls import Services

logger = structlog.get_logger()


@click.group(name="otu")
@path_option
@click.pass_context
def otu(ctx: click.Context, path: Path) -> None:
    """Manage OTUs."""
    ctx.obj = ctx.with_resource(locked_repo(path))


@otu.command(name="create")
@click.argument(
    "accessions_",
    callback=validate_no_duplicate_accessions,
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@pass_repo
@ignore_cache_option
def otu_create(
    repo: Repo,
    accessions_: list[str],
    ignore_cache: bool,
) -> None:
    """Create a new OTU.

    OTUs are created from a list of accessions and a taxonomy ID. The associated Genbank
    records are fetched and used to create the first isolate and a plan.
    """
    if len(accessions_) != len(set(accessions_)):
        click.echo("Duplicate accessions were provided.", err=True)
        sys.exit(1)

    services = Services(repo, NCBIClient(ignore_cache))

    try:
        otu_ = services.otu.create(accessions_)
    except ValueError as e:
        click.echo(e, err=True)
        sys.exit(1)

    if otu_ is None:
        click.echo("OTU was not created correctly.", err=True)
        sys.exit(1)

    print_otu(otu_)


@otu.command(name="get")
@click.argument("IDENTIFIER", type=str)
@click.option(
    "--as-json",
    "--json",
    metavar="JSON",
    is_flag=True,
    help="Output in JSON form",
)
@pass_repo
def otu_get(repo: Repo, identifier: str, as_json: bool) -> None:
    """Get an OTU by its unique ID or taxonomy ID.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    services = Services(repo, NCBIClient(False))
    otu_ = services.otu.get_otu(identifier)

    if otu_ is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    if as_json:
        print_otu_as_json(otu_)
    else:
        print_otu(otu_)


@otu.command(name="list")
@pass_repo
def otu_list(repo: Repo) -> None:
    """List all OTUs in the repository."""
    print_otu_list(repo.iter_minimal_otus())


@otu.command(name="list-events")
@click.argument("IDENTIFIER", type=str)
@pass_repo
def otu_event_logs(repo: Repo, identifier: str) -> None:
    """Print a log of OTU events to console.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    services = Services(repo, NCBIClient(False))
    otu_ = services.otu.get_otu(identifier)

    if otu_ is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    print_otu_event_log(list(repo.iter_otu_events(otu_.id)))


@otu.command(name="update")
@click.argument("IDENTIFIER", type=str)
@ignore_cache_option
@pass_repo
def otu_update(
    repo: Repo,
    identifier: str,
    ignore_cache: bool,
) -> None:
    """Comprehensively update an OTU with the latest data from NCBI.

    This command performs three operations:
    1. Promotes GenBank accessions to RefSeq equivalents where available
    2. Upgrades outdated sequence versions (e.g., v1 → v2)
    3. Adds new isolates from newly available accessions
    """
    services = Services(repo, NCBIClient(ignore_cache))
    otu_ = services.otu.get_otu(identifier)

    if otu_ is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    services.otu.update(otu_.id)


@otu.command(name="exclude-accessions")  # type: ignore
@click.argument("IDENTIFIER", type=str)
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@pass_repo
def otu_exclude_accessions(
    repo: Repo,
    identifier: str,
    accessions_: list[str],
) -> None:
    """Exclude the given accessions from this OTU.

    Any duplicate accessions will be ignored. Only one exclusion per unique accession
    will be made.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    services = Services(repo, NCBIClient(False))
    otu_ = services.otu.get_otu(identifier)

    if otu_ is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    services.otu.exclude_accessions(otu_.id, set(accessions_))


@otu.command(name="allow-accessions")  # type: ignore
@click.argument("IDENTIFIER", type=str)
@click.argument(
    "accessions_",
    metavar="ACCESSIONS",
    nargs=-1,
    type=str,
    required=True,
)
@pass_repo
def otu_allow_accessions(
    repo: Repo,
    identifier: str,
    accessions_: list[str],
) -> None:
    """Allow the given excluded accessions back into the OTU.

    Any duplicate accessions will be ignored. Only one exclusion per unique accession
    will be made.

    IDENTIFIER is a taxonomy ID or unique OTU ID (>8 characters)
    """
    services = Services(repo, NCBIClient(False))
    otu_ = services.otu.get_otu(identifier)

    if otu_ is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    services.otu.allow_accessions(otu_.id, set(accessions_))
