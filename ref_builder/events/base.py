import datetime
from typing import Generic, TypeVar

from pydantic import UUID4, BaseModel, computed_field

from ref_builder.models.accession import Accession
from ref_builder.models.otu import OTU

TEventData = TypeVar("TEventData", bound="EventData")
TEventQuery = TypeVar("TEventQuery", bound="EventQuery")


class EventQuery(BaseModel):
    """A base class for a query that targets an event at a specific resource."""


class EventData(BaseModel):
    """Represents the data for an event."""


class Event(BaseModel, Generic[TEventData, TEventQuery]):
    """The base event."""

    id: int
    """The unique identifier for the event.

    Event IDs are serially incremented integers.
    """

    data: TEventData
    """The data associated with the event."""

    query: TEventQuery
    """The query targeting the event at a specific resource."""

    timestamp: datetime.datetime
    """When the event occurred."""

    @computed_field
    @property
    def type(self) -> str:
        """The type of the event as a string."""
        return self.__class__.__name__


class ApplicableEvent(Event[TEventData, TEventQuery]):
    def apply(self, otu: OTU) -> OTU:
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
    """An event query that targets a specific sequence by its accession."""

    accession: Accession


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
