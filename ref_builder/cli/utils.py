import sys
from uuid import UUID

import click

from ref_builder.repo import Repo

pass_repo = click.make_pass_decorator(Repo)


def get_otu_isolate_ids_from_identifier(
    repo: Repo, identifier: str
) -> tuple[UUID, UUID]:
    """Return an isolate ID from the repo if identifier matches a single isolate.

    Handles cases where the isolate cannot be found.
    """
    try:
        isolate_id = UUID(identifier)
    except ValueError:
        click.echo("Invalid isolate ID format.", err=True)
        sys.exit(1)

    if (otu_id := repo.get_otu_id_by_isolate_id(isolate_id)) is not None:
        if isolate_id in repo.get_otu(otu_id).isolate_ids:
            return otu_id, isolate_id

    click.echo("Isolate ID could not be found.", err=True)
    sys.exit(1)
