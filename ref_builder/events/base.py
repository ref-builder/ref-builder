import datetime

from pydantic import UUID4, BaseModel, computed_field

from ref_builder.otu.builders.otu import OTUBuilder


class EventQuery(BaseModel):
    """A base class for a query that targets an event at a specific resource."""


class EventData(BaseModel):
    """Represents the data for an event."""


class Event(BaseModel):
    """The base event."""

    id: int
    """The unique identifier for the event.

    Event IDs are serially incremented integers.
    """

    data: EventData
    """The data associated with the event."""

    query: EventQuery
    """The query targeting the event at a specific resource."""

    timestamp: datetime.datetime
    """When the event occurred."""

    @computed_field
    @property
    def type(self) -> str:
        """The type of the event as a string."""
        return self.__class__.__name__


class ApplicableEvent(Event):
    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        return otu


class RepoQuery(EventQuery):
    """An event query that targets an event at the repository."""

    repository_id: UUID4


class OTUQuery(EventQuery):
    """An event query that targets an event at an OTU."""

    otu_id: UUID4


class IsolateQuery(OTUQuery):
    """An event query that targets an event at an isolate in a specific OTU."""

    isolate_id: UUID4


class SequenceQuery(OTUQuery):
    """An event query that targets an event at a sequence in a specific isolate and
    OTU.
    """

    sequence_id: UUID4


class LinkSequenceQuery(IsolateQuery):
    """An event query that targets an event at a sequence in a specific isolate and
    OTU.
    """

    sequence_id: UUID4


class EventMetadata(BaseModel):
    """Minimal metadata for an Event."""

    id: int
    """The unique identifier for the event.

    Event IDs are serially incremented integers.
    """

    otu_id: UUID4 | None
    """The unique identifier of an OTU.

    Is None if event scope is repo-wide.
    """

    timestamp: datetime.datetime
    """When the event occurred."""
