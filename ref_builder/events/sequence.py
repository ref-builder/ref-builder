from ref_builder.errors import HydrationIsolateError, HydrationSequenceError
from ref_builder.events.base import (
    ApplicableEvent,
    EventData,
    SequenceQuery,
)
from ref_builder.models.otu import OTU
from ref_builder.models.sequence import Sequence


class UpdateSequenceData(EventData):
    """Data for a sequence version update event."""

    sequence: Sequence
    """The new sequence to replace the old one with."""


class UpdateSequence(ApplicableEvent[UpdateSequenceData, SequenceQuery]):
    """An event that updates a sequence to a newer version.

    The sequence is identified by its accession in the query.
    """

    def apply(self, otu: OTU) -> OTU:
        """Apply update by replacing the targeted sequence with the new one.

        1. Find the isolate containing the sequence (by accession)
        2. Remove the old sequence from the isolate
        3. Add the new sequence to the isolate
        4. Mark the old accession as excluded
        """
        isolate = None

        for iso in otu.isolates:
            if iso.get_sequence(self.query.accession) is not None:
                isolate = iso
                break

        if isolate is None:
            raise HydrationIsolateError("Could not find isolate")

        old_sequence = isolate.get_sequence(self.query.accession)

        if old_sequence is None:
            raise HydrationSequenceError(
                "Could not find referenced sequence accession."
            )

        isolate.sequences = [
            seq for seq in isolate.sequences if seq.accession != self.query.accession
        ]
        isolate.sequences.append(self.data.sequence)

        return otu.model_validate(otu)
