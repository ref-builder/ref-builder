from pydantic import UUID4

from ref_builder.events.base import Event, EventData, RepoQuery
from ref_builder.models.repo import RepoSettings


class CreateRepoData(EventData):
    """The data for a :class:`CreateRepo` event."""

    id: UUID4
    name: str
    organism: str
    settings: RepoSettings


class CreateRepo(Event[CreateRepoData, RepoQuery]):
    """An event that creates a new repository.

    This event is always the first event in a repository's event log.
    """
