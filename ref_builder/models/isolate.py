from dataclasses import dataclass
from enum import StrEnum

from pydantic import (
    UUID4,
    BaseModel,
    field_serializer,
    field_validator,
    model_validator,
)

from ref_builder.models.accession import Accession
from ref_builder.models.sequence import Sequence


class IsolateNameType(StrEnum):
    """Possible types for isolate names.

    **Ordered by priority**. Do not reorder attributes.

    Isolate name types were previously called "source types". They are referred to this
    way in Virtool.
    """

    ISOLATE = "isolate"
    STRAIN = "strain"
    CLONE = "clone"
    VARIANT = "variant"
    GENOTYPE = "genotype"
    SEROTYPE = "serotype"


@dataclass(frozen=True)
class IsolateName:
    """A name for an isolate.

    The isolate name consists of a type and a value.

    For example, in the isolate name "Isolate PPSMV2-Badnapur", the type is "isolate"
    and the value is "PPSMV2-Badnapur.

    In the Genbank record for a sequence that belongs to this isolate, the source table
    contains:

    .. code-block:: text

        /isolate="PPSMV2-Badnapur"

    """

    type: IsolateNameType
    """The type of sub-species categorization."""

    value: str
    """The name of this subcategory."""

    def __str__(self) -> str:
        """Return the isolate name as a formatted string."""
        return f"{self.type.capitalize()} {self.value}"


class Isolate(BaseModel):
    """A class representing the fields of an isolate."""

    id: UUID4
    """The isolate id."""

    name: IsolateName | None
    """The isolate's name."""

    sequences: list[Sequence]
    """The isolate's sequences."""

    taxid: int
    """The NCBI Taxonomy ID for this isolate."""

    def get_sequence(
        self,
        accession: Accession,
    ) -> Sequence | None:
        """Return a sequence with the given accession if it exists in the isolate,
        else None.
        """
        for sequence in self.sequences:
            if sequence.accession == accession:
                return sequence

        return None

    @property
    def accessions(self) -> set[str]:
        """A set of accession numbers for sequences in the isolate."""
        return {sequence.accession.key for sequence in self.sequences}

    @property
    def versioned_accessions(self) -> set[Accession]:
        """A set of versioned accessions contained in this isolate."""
        return {sequence.accession for sequence in self.sequences}

    @property
    def is_refseq(self) -> bool:
        """Return True if this isolate was sourced from NCBI's RefSeq database."""
        if self.sequences:
            return self.sequences[0].accession.is_refseq

        return False

    @field_validator("name", mode="before")
    @classmethod
    def convert_name(
        cls,
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

    @model_validator(mode="after")
    def check_accession_consistency(self) -> "Isolate":
        """Validate that all sequences in multipartite isolate use consistent accession provenance."""
        if len(self.sequences) > 1:
            refseq_accessions = []
            genbank_accessions = []

            for sequence in self.sequences:
                if sequence.accession.is_refseq:
                    refseq_accessions.append(sequence.accession.key)
                else:
                    genbank_accessions.append(sequence.accession.key)

            if refseq_accessions and genbank_accessions:
                raise ValueError(
                    f"Combination of RefSeq and non-RefSeq sequences found in multipartite isolate. "
                    f"RefSeq: {', '.join(refseq_accessions)}; Non-RefSeq: {', '.join(genbank_accessions)}"
                )

        return self
