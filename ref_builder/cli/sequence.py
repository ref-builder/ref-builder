import sys
from pathlib import Path

import click

from ref_builder.cli.options import path_option
from ref_builder.cli.utils import get_otu_sequence_ids_from_identifier, pass_repo
from ref_builder.console import print_sequence, print_sequence_as_json
from ref_builder.repo import Repo, locked_repo


@click.group(name="sequence")
@click.pass_context
@path_option
def sequence(ctx: click.Context, path: Path) -> None:
    """Manage sequences."""
    ctx.obj = ctx.with_resource(locked_repo(path))

@sequence.command(name="get")
@click.argument("IDENTIFIER", type=str)
@click.option("--id-only", is_flag=True, help="Only output sequence ID.")
@click.option(
    "--json",
    "json_",
    is_flag=True,
    help="Output in JSON form.",
)
@pass_repo
def sequence_get(repo: Repo, identifier: str, id_only: bool, json_: bool) -> None:
    """Get an isolate with a UUID corresponding to IDENTIFIER.

    IDENTIFIER can be an accession number or a unique isolate ID (>8 characters)
    """
    otu_id, sequence_id = get_otu_sequence_ids_from_identifier(repo, identifier)

    otu_ = repo.get_otu(otu_id)

    sequence_ = otu_.get_sequence_by_id(sequence_id)

    if sequence_ is None:
        click.echo("Sequence could not be found.", err=True)
        sys.exit(1)

    if id_only:
        click.echo(sequence_.id)

        sys.exit(0)

    if json_:
        print_sequence_as_json(sequence_)

    else:
        print_sequence(sequence_)


@sequence.command(name="get-otu-id")
@click.argument("IDENTIFIER", type=str)
@pass_repo
def sequence_get_otu_id(repo: Repo, identifier: str) -> None:
    """Get the OTU id containing a sequence associated with the given IDENTIFIER.

    IDENTIFIER is an unique isolate ID (>8 characters)
    """
    otu_id, sequence_id = get_otu_sequence_ids_from_identifier(repo, identifier)

    otu_ = repo.get_otu(otu_id)

    if otu_ is None:
        sys.exit(1)

    else:
        click.echo(otu_.id)

@sequence.command(name="get-isolate-ids")
@click.argument("IDENTIFIER", type=str)
@pass_repo
def sequence_get_isolate_ids(repo: Repo, identifier: str) -> None:
    """Get the isolate id/s containing a sequence associated with the given IDENTIFIER.

    IDENTIFIER can be an accession number or a unique isolate ID (>8 characters)
    """
    otu_id, sequence_id = get_otu_sequence_ids_from_identifier(repo, identifier)

    otu_ = repo.get_otu(otu_id)

    if otu_ is None:
        sys.exit(1)

    isolate_ids = otu_.get_isolate_ids_containing_sequence_id(sequence_id)

    if not isolate_ids:
        sys.exit(1)

    else:
        click.echo(", ".join([str(isolate_id) for isolate_id in isolate_ids]))
