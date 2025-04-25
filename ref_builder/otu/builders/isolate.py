from uuid import UUID

from pydantic import field_serializer, field_validator

from ref_builder.otu.builders.sequence import SequenceBuilder
from ref_builder.otu.models import IsolateModel
from ref_builder.utils import IsolateName


class IsolateBuilder(IsolateModel):
    """Represents an isolate in a Virtool reference repository."""

    sequences: list[SequenceBuilder]
    """A list of sequences contained by this isolate."""

    @property
    def accessions(self) -> set[str]:
        """A set of accession numbers for sequences in the isolate."""
        return {sequence.accession.key for sequence in self.sequences}

    @property
    def sequence_ids(self) -> set[UUID]:
        """A set of UUIDs for sequences in the isolate."""
        return {sequence.id for sequence in self.sequences}

    def add_sequence(self, sequence: SequenceBuilder) -> None:
        """Add a sequence to the isolate."""
        self.sequences.append(sequence)

    def delete_sequence(self, sequence_id: UUID) -> None:
        """Delete a sequence from a given isolate."""
        sequence = self.get_sequence_by_id(sequence_id)

        if not sequence:
            raise ValueError("This sequence ID does not exist")

        for sequence in self.sequences:
            if sequence.id == sequence_id:
                self.sequences.remove(sequence)

    def get_sequence_by_accession(
        self,
        accession: str,
    ) -> SequenceBuilder | None:
        """Return a sequence with the given accession if it exists in the isolate,
        else None.
        """
        for sequence in self.sequences:
            if sequence.accession.key == accession:
                return sequence

        return None

    def get_sequence_by_id(self, sequence_id: UUID) -> SequenceBuilder | None:
        """Return a sequence with the given ID if it exists in the isolate,
        else None.
        """
        for sequence in self.sequences:
            if sequence.id == sequence_id:
                return sequence

        return None

    @field_validator("name", mode="before")
    @classmethod
    def convert_name(
        cls: "IsolateBuilder",
        value: dict | IsolateName | None,
    ) -> IsolateName | None:
        """Convert the name to an IsolateName object."""
        if value is None:
            return value

        if isinstance(value, IsolateName):
            return value

        if isinstance(value, dict):
            return IsolateName(**value)

        raise ValueError(f"Invalid type for name: {type(value)}")

    @field_serializer("name")
    def serialize_name(self, name: IsolateName | None) -> dict[str, str] | None:
        """Serialize the isolate name."""
        if name is None:
            return None

        return {
            "type": name.type,
            "value": name.value,
        }
