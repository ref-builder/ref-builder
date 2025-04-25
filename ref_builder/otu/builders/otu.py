from uuid import UUID

from pydantic import UUID4

from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.otu.builders.sequence import SequenceBuilder
from ref_builder.otu.models import OTUModel
from ref_builder.utils import IsolateName


class OTUBuilder(OTUModel):
    """Represents an OTU in a Virtool reference repository."""

    isolates: list[IsolateBuilder]
    """Isolates contained in this OTU."""

    representative_isolate: UUID4 | None
    """The UUID of the representative isolate for the OTU."""

    _isolates_by_id: dict[UUID4:IsolateBuilder]
    """A dictionary of isolates indexed by isolate UUID"""

    _sequences_by_id: dict[UUID4:SequenceBuilder]
    """A dictionary of sequences indexed by sequence UUID"""

    def __init__(self, **data) -> None:
        super().__init__(**data)

        self._sequences_by_id = {}
        self._isolates_by_id = {}
        for isolate in self.isolates:
            self._isolates_by_id[isolate.id] = isolate
            for sequence in isolate.sequences:
                self._sequences_by_id[sequence.id] = sequence

    @property
    def accessions(self) -> set[str]:
        """A set of accessions contained in this isolate."""
        return {sequence.accession.key for sequence in self._sequences_by_id.values()}

    @property
    def blocked_accessions(self) -> set[str]:
        """Accessions that should not be considered for addition to the OTU.

        This includes accessions that already exist in the OTU and accessions that have
        been explicitly excluded.
        """
        return self.accessions | self.excluded_accessions

    @property
    def isolate_ids(self) -> set[UUID4]:
        """A set of UUIDs for isolates in the OTU."""
        return set(self._isolates_by_id.keys())

    def add_isolate(self, isolate: IsolateBuilder) -> None:
        """Add an isolate to the OTU."""
        self.isolates.append(isolate)
        self._isolates_by_id[isolate.id] = isolate

        for sequence in isolate.sequences:
            self.add_sequence(sequence)

    def add_sequence(self, sequence: SequenceBuilder) -> None:
        """Add a sequence to a given isolate."""
        self._sequences_by_id[sequence.id] = sequence

    def get_sequence_by_id(self, sequence_id: UUID) -> SequenceBuilder | None:
        return self._sequences_by_id.get(sequence_id)

    def delete_isolate(self, isolate_id: UUID4) -> None:
        """Remove an isolate from the OTU."""
        if self._isolates_by_id.get(isolate_id) is None:
            raise ValueError(f"Isolate {isolate_id} does not exist")

        for isolate in self.isolates:
            if isolate.id == isolate_id:
                for sequence_id in isolate.sequence_ids:
                    if self.get_sequence_by_id(sequence_id) is None:
                        raise ValueError(
                            f"Sequence {sequence_id} not found in the sequence list"
                        )
                    self.delete_sequence(sequence_id)

                self.isolates.remove(isolate)

                self._isolates_by_id.pop(isolate_id)

                break

    def delete_sequence(self, sequence_id: UUID4) -> None:
        """Delete a sequence from a given isolate. Used only during rehydration."""
        self._sequences_by_id.pop(sequence_id)

    def get_isolate(self, isolate_id: UUID4) -> IsolateBuilder | None:
        """Get isolate associated with a given ID.

        Returns None if no such isolate exists.

        :param isolate_id: The UUID of the isolate to retrieve
        :return: the isolate or ``None``
        """
        return self._isolates_by_id.get(isolate_id)

    def get_isolate_id_by_name(self, name: IsolateName) -> UUID4 | None:
        """Get the ID for the isolate with the passed ``name``.

        Returns None if no such isolate exists.

        :param name: The name of the isolate to retrieve
        :return: The isolate ID or ``None``

        """
        for isolate in self.isolates:
            if isolate.name == name:
                return isolate.id

        return None

    def get_sequence_by_accession(
        self,
        accession: str,
    ) -> SequenceBuilder | None:
        """Return a sequence corresponding to given accession
        if it exists in this OTU.
        """
        if accession not in self.accessions:
            return None

        for isolate in self.isolates:
            if (sequence := isolate.get_sequence_by_accession(accession)) is not None:
                return sequence

        raise ValueError(f"Accession {accession} found in index, but not in data")

    def get_isolate_ids_containing_sequence_id(self, sequence_id: UUID4) -> set[UUID4]:
        """Return a set of isolate IDs where the isolate contains the given sequence."""
        containing_isolate_ids = set()

        if sequence_id not in self._sequences_by_id:
            return containing_isolate_ids

        for isolate in self.isolates:
            if sequence_id in isolate.sequence_ids:
                containing_isolate_ids.add(isolate.id)

        if containing_isolate_ids:
            return containing_isolate_ids

        raise ValueError(f"Sequence ID {sequence_id} found in index, but not in data")

    def link_sequence(
        self, isolate_id: UUID4, sequence_id: UUID4
    ) -> SequenceBuilder | None:
        """Link the given sequence to the given isolate."""
        self.get_isolate(isolate_id).add_sequence(self.get_sequence_by_id(sequence_id))

        return self.get_isolate(isolate_id).get_sequence_by_id(sequence_id)

    def unlink_sequence(self, isolate_id: UUID4, sequence_id: UUID4) -> None:
        """Unlink the given sequence from the given isolate. Used only during rehydration."""
        self.get_isolate(isolate_id).delete_sequence(sequence_id)
