from pydantic import UUID4, field_serializer, field_validator

from ref_builder.errors import (
    HydrationIsolateError,
    HydrationSequenceError,
)
from ref_builder.events.base import (
    ApplicableEvent,
    EventData,
    IsolateQuery,
)
from ref_builder.models.accession import Accession
from ref_builder.models.isolate import Isolate, IsolateName
from ref_builder.models.otu import OTU
from ref_builder.models.sequence import Sequence


class CreateIsolateData(EventData):
    """The data for a :class:`CreateIsolate` event."""

    id: UUID4
    name: IsolateName | None
    sequences: list[Sequence]
    taxid: int


class CreateIsolate(ApplicableEvent[CreateIsolateData, IsolateQuery]):
    """An event that creates an isolate for a specific OTU."""

    def apply(self, otu: OTU) -> OTU:
        """Add isolate with sequences to OTU and return."""
        otu.isolates.append(
            Isolate(
                id=self.data.id,
                name=self.data.name,
                taxid=self.data.taxid,
                sequences=self.data.sequences,
            ),
        )

        return otu


class DeleteIsolateData(EventData):
    """The data for a :class:`DeleteIsolate` event."""

    message: str


class DeleteIsolate(ApplicableEvent[DeleteIsolateData, IsolateQuery]):
    """An isolate deletion event."""

    def apply(self, otu: OTU) -> OTU:
        """Delete the specified isolate and return."""
        if otu.get_isolate(self.query.isolate_id) is None:
            raise HydrationIsolateError("Could not find referenced isolate.")

        otu.isolates = [
            isolate for isolate in otu.isolates if isolate.id != self.query.isolate_id
        ]

        return otu.model_validate(otu)


class PromoteIsolateData(EventData):
    """Data for an isolate promotion event.

    Contains a map of old accessions to new RefSeq sequences.
    """

    map: dict[Accession, Sequence]

    @field_serializer("map")
    def serialize_map(self, map_: dict[Accession, Sequence]) -> dict[str, dict]:
        """Serialize map with Accession keys to dict with string keys."""
        return {str(acc): seq.model_dump() for acc, seq in map_.items()}

    @field_validator("map", mode="before")
    @classmethod
    def deserialize_map(cls, value: dict) -> dict[Accession, Sequence]:
        """Deserialize map from dict with string keys to Accession keys."""
        if not value:
            return {}

        # Check if already deserialized (keys are Accession objects)
        first_key = next(iter(value.keys()))
        if isinstance(first_key, Accession):
            return value

        # Deserialize from JSON format
        return {
            Accession.from_string(key): Sequence(**seq_data)
            for key, seq_data in value.items()
        }


class PromoteIsolate(ApplicableEvent[PromoteIsolateData, IsolateQuery]):
    """An event that promotes GenBank sequences to RefSeq in an isolate.

    This event is organized at the isolate-level because all sequences in an isolate
    must be either Refseq or not.
    """

    def apply(self, otu: OTU) -> OTU:
        """Apply promotion by replacing sequences in the isolate.

        For each old accession -> new sequence mapping:

        1. Validate the isolate exists
        2. Validate all old sequences exist in the isolate
        3. Add new RefSeq sequences to the isolate
        4. Remove old GenBank sequences from the isolate
        5. Mark old accessions as excluded
        """
        isolate = otu.get_isolate(self.query.isolate_id)

        if isolate is None:
            raise HydrationIsolateError("Could not find referenced isolate.")

        for old_accession in self.data.map:
            if isolate.get_sequence(old_accession) is None:
                raise HydrationSequenceError(
                    "Could not find referenced sequence accession."
                )

        for old_accession, new_sequence_data in self.data.map.items():
            isolate.sequences.append(
                Sequence(
                    accession=new_sequence_data.accession,
                    definition=new_sequence_data.definition,
                    segment=new_sequence_data.segment,
                    sequence=new_sequence_data.sequence,
                )
            )

            isolate.sequences = [
                seq for seq in isolate.sequences if seq.accession != old_accession
            ]

            otu.promoted_accessions.add(old_accession.key)

        # Trigger Pydantic validation to rebuild OTU lookup dictionaries
        return otu.model_validate(otu)
