import sys
from uuid import UUID

import click

from ref_builder.errors import InvalidInputError, PartialIDConflictError
from ref_builder.repo import Repo

pass_repo = click.make_pass_decorator(Repo)


def get_otu_isolate_ids_from_identifier(
    repo: Repo, identifier: str
) -> tuple[UUID, UUID]:
    """Return an isolate ID from the repo if identifier matches a single isolate.

    Handles cases where the isolate cannot be found.
    """
    isolate_id = None

    if len(identifier) == 32:
        try:
            isolate_id = UUID(identifier)
        except ValueError:
            isolate_id = None

    else:
        try:
            isolate_id = repo.get_isolate_id_by_partial(identifier)

        except PartialIDConflictError:
            click.echo(
                "Partial ID too short to narrow down results.",
                err=True,
            )

        except InvalidInputError:
            click.echo(
                "Partial ID segment must be at least 8 characters long.", err=True
            )
            sys.exit(1)

    if isolate_id is not None:
        if (otu_id := repo.get_otu_id_by_isolate_id(isolate_id)) is not None:
            if isolate_id in repo.get_otu(otu_id).isolate_ids:
                return otu_id, isolate_id

    click.echo("Isolate ID could not be found.", err=True)
    sys.exit(1)
