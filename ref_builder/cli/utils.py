import sys
import warnings
from uuid import UUID

import click

from ref_builder.errors import InvalidInputError, PartialIDConflictError
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.repo import Repo
from ref_builder.utils import OTUDeletedWarning

pass_repo = click.make_pass_decorator(Repo)

UUID_STRING_LENGTH = 36


def get_otu_from_identifier(repo: Repo, identifier: str) -> OTUBuilder:
    """Return an OTU id from the repo if identifier matches a single OTU.

    The identifier can either be a stringified UUIDv4, a truncated portion
    of a UUID, a NCBI Taxonomy ID or an acronym associated with the OTU.

    Raise a PartialIDConflictError if >1 OTU is found.

    :param repo: the repository to be searched.
    :param identifier: a non-UUID identifier.
        Can be an integer Taxonomy ID, acronym or truncated partial UUID.
    :return: the UUID of the OTU or ``None``
    """
    otu_id = None

    if len(identifier) == UUID_STRING_LENGTH:
        try:
            otu_id = UUID(identifier)
        except ValueError:
            otu_id = None

    elif identifier.isnumeric():
        try:
            taxid = int(identifier)
        except ValueError:
            pass
        else:
            otu_id = repo.get_otu_id_by_taxid(taxid)

    else:
        if (otu_id := repo.get_otu_id_by_acronym(identifier)) is not None:
            return repo.get_otu(otu_id)

        try:
            otu_id = repo.get_otu_id_by_partial(identifier)

        except PartialIDConflictError:
            click.echo(
                "Partial ID too short to narrow down results.",
                err=True,
            )

        except InvalidInputError as e:
            click.echo(
                e,
                err=True,
            )

    if otu_id is None:
        click.echo("OTU not found.", err=True)
        sys.exit(1)

    with warnings.catch_warnings(
        category=OTUDeletedWarning, record=True
    ) as warning_list:
        otu_ = repo.get_otu(otu_id)

    if otu_ is None:
        if warning_list:
            click.echo("OTU has been deleted.", err=True)

        else:
            click.echo("OTU not found.", err=True)

        sys.exit(1)

    return otu_


def get_otu_isolate_ids_from_identifier(repo: Repo, identifier: str) -> (UUID, UUID):
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
