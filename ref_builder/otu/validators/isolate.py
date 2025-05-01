import warnings
from uuid import UUID

from pydantic import (
    UUID4,
    BaseModel,
    ConfigDict,
    Field,
    field_validator,
)

from ref_builder.otu.models import IsolateModel, LegacyId
from ref_builder.otu.validators.sequence import Sequence, SequenceBase
from ref_builder.otu.validators.utils import IsolateInconsistencyWarning
from ref_builder.utils import Accession, IsolateName, is_refseq


class IsolateBase(IsolateModel):
    """A class representing an isolate with basic validation."""

    id: UUID4
    """The isolate id."""

    legacy_id: LegacyId | None
    """A string based ID carried over from a legacy Virtool reference repository.

    It the isolate was not migrated from a legacy repository, this will be `None`.
    """

    name: IsolateName | None
    """The isolate's name."""

    sequences: list[SequenceBase]
    """The isolates sequences."""

    @property
    def refseq(self) -> bool:
        """Return True if this isolate was sourced from NCBI's RefSeq database."""
        if self.sequences:
            return is_refseq(self.sequences[0].accession.key)

        return False

    def get_sequence_by_accession(
        self,
        accession: Accession,
    ) -> SequenceBase | None:
        """Get a sequence by its accession.

        Return ``None`` if the sequence is not found.

        :param accession: the accession of the sequence to retrieve
        :return: the sequence with the given accession, or None
        """
        for sequence in self.sequences:
            if sequence.accession == accession:
                return sequence

        return None

    def get_sequence_by_id(self, sequence_id: UUID) -> SequenceBase | None:
        """Get a sequence by its ID.

        Return ``None`` if the sequence is not found.

        :param sequence_id: the ID of the sequence to retrieve
        :return: the sequence with the given ID, or None
        """
        for sequence in self.sequences:
            if sequence.id == sequence_id:
                return sequence

        return None

    @field_validator("sequences", mode="before")
    @classmethod
    def convert_sequence_models(
        cls,
        v: list[dict | SequenceBase | Sequence | BaseModel],
    ) -> list[SequenceBase]:
        """Automatically revalidate sequence if not already validated."""
        if not v or isinstance(v[0], dict | SequenceBase):
            return v

        return [SequenceBase.model_validate(sequence.model_dump()) for sequence in v]


class Isolate(IsolateBase):
    """A class representing an isolate with full validation."""

    model_config = ConfigDict(validate_assignment=True)

    sequences: list[Sequence] = Field(min_length=1)
    """The isolates sequences.

    A valid isolate must have at least one sequence.
    """

    @field_validator("sequences", mode="before")
    @classmethod
    def convert_sequence_models(
        cls,
        v: list[dict | SequenceBase | Sequence | BaseModel],
    ) -> list[Sequence]:
        if not v or isinstance(v[0], dict | SequenceBase):
            return v

        return [Sequence.model_validate(sequence.model_dump()) for sequence in v]

    @field_validator("sequences", mode="after")
    @classmethod
    def check_accession_refseq_or_insdc(cls, v: list[Sequence]) -> list[Sequence]:
        """Check if all sequence accessions are all from RefSeq or all from INSDC.
        If not, warn the user.
        """
        if len({s.refseq for s in v}) > 1:
            warnings.warn(
                "Combination of RefSeq and non-RefSeq sequences found in multipartite isolate: "
                + f"{[sequence.accession.key for sequence in v]}",
                IsolateInconsistencyWarning,
                stacklevel=1,
            )

        return v

    @field_validator("name", mode="after")
    @classmethod
    def check_isolate_name_has_value(cls, v: IsolateName | None) -> IsolateName | None:
        """Check if ``name.value`` is a string with at least one character."""
        if v is None:
            return v

        if isinstance(v, IsolateName) and v.value:
            return v

        raise ValueError(f"Invalid isolate name {v}.")
