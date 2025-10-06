from pydantic import UUID4

from ref_builder.events.base import (
    ApplicableEvent,
    EventData,
    IsolateQuery,
)
from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.utils import IsolateName


class CreateIsolateData(EventData):
    """The data for a :class:`CreateIsolate` event."""

    id: UUID4
    name: IsolateName | None


class CreateIsolate(ApplicableEvent[CreateIsolateData, IsolateQuery]):
    """An event that creates an isolate for a specific OTU."""

    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        """Add isolate to OTU and return."""
        otu.add_isolate(
            IsolateBuilder(
                id=self.data.id,
                name=self.data.name,
                sequences=[],
            ),
        )

        return otu


class LinkSequenceData(EventData):
    """The data for a :class:`LinkSequence` event."""

    sequence_id: UUID4


class LinkSequence(ApplicableEvent[LinkSequenceData, IsolateQuery]):
    """An event that links an existing sequence to an isolate."""

    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        """Add specified sequence to specified isolate and return."""
        otu.link_sequence(
            isolate_id=self.query.isolate_id,
            sequence_id=self.data.sequence_id,
        )

        return otu


class UnlinkSequenceData(EventData):
    """The data for a :class:`LinkSequence` event."""

    sequence_id: UUID4


class UnlinkSequence(ApplicableEvent[UnlinkSequenceData, IsolateQuery]):
    """An event that unlinks an existing sequence from an isolate."""

    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        """Unlink specified sequence from specified isolate and return OTU."""
        otu.unlink_sequence(
            isolate_id=self.query.isolate_id,
            sequence_id=self.data.sequence_id,
        )

        return otu


class DeleteIsolateData(EventData):
    """The data for a :class:`DeleteIsolate` event."""

    rationale: str


class DeleteIsolate(ApplicableEvent[DeleteIsolateData, IsolateQuery]):
    """An isolate deletion event."""

    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        """Delete the specified isolate and return."""
        otu.delete_isolate(self.query.isolate_id)

        return otu
