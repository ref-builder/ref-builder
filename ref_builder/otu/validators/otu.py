import warnings
from collections import Counter

from pydantic import (
    UUID4,
    BaseModel,
    ConfigDict,
    Field,
    field_validator,
    model_validator,
)
from pydantic_core import PydanticCustomError

from ref_builder.models import Molecule
from ref_builder.otu.models import LegacyId, OTUModel
from ref_builder.otu.validators.isolate import Isolate, IsolateBase
from ref_builder.otu.validators.sequence import Sequence, SequenceBase
from ref_builder.plan import Plan, PlanWarning


class OTUBase(OTUModel):
    """A class representing an OTU with basic validation."""

    id: UUID4
    """The OTU id."""

    acronym: str
    """The OTU acronym (eg. TMV for Tobacco mosaic virus)."""

    legacy_id: LegacyId | None
    """A string based ID carried over from a legacy Virtool reference repository."""

    molecule: Molecule
    """The type of molecular information contained in this OTU."""

    name: str = Field(min_length=1)
    """The name of the OTU (eg. TMV for Tobacco mosaic virus)"""

    plan: Plan
    """The plan for the OTU."""

    taxid: int = Field(gt=0)
    """The NCBI Taxonomy id for this OTU."""

    excluded_accessions: set[str]
    """A set of accessions that should not be retrieved in future fetch operations."""

    isolates: list[IsolateBase]
    """Isolates contained in this OTU."""

    representative_isolate: UUID4 | None
    """The UUID of the representative isolate of this OTU."""

    @property
    def sequences(self) -> list[SequenceBase]:
        """Sequences contained in this OTU."""
        return [sequence for isolate in self.isolates for sequence in isolate.sequences]

    @field_validator("isolates", mode="before")
    @classmethod
    def convert_isolate_models(
        cls,
        v: list[dict | IsolateBase | Isolate | BaseModel],
    ) -> list[IsolateBase]:
        """Automatically revalidate isolates if not already validated."""
        if not v or isinstance(v[0], dict | IsolateBase):
            return v

        return [IsolateBase.model_validate(isolate.model_dump()) for isolate in v]


class OTU(OTUBase):
    """A class representing an OTU with full validation."""

    model_config = ConfigDict(validate_assignment=True)

    isolates: list[Isolate] = Field(min_length=1)
    """Isolates contained in this OTU.

    A valid OTU must have at least one isolate.
    """

    representative_isolate: UUID4
    """The UUID of the representative isolate of this OTU.

    A valid OTU must have a representative isolate.
    """

    @property
    def sequences(self) -> list[Sequence]:
        """Sequences contained in this OTU."""
        return [sequence for isolate in self.isolates for sequence in isolate.sequences]

    @field_validator("isolates", mode="before")
    @classmethod
    def convert_isolate_models(
        cls,
        v: list[dict | Isolate | IsolateBase | BaseModel],
    ) -> list[Isolate]:
        """Automatically revalidate isolates if not already validated."""
        if not v or isinstance(v[0], dict | Isolate):
            return v

        return [Isolate.model_validate(isolate.model_dump()) for isolate in v]

    @field_validator("plan", mode="after")
    def check_plan_required(cls, plan: Plan) -> Plan:
        """Issue a warning if the plan has no required segments."""
        if not plan.required_segments:
            warnings.warn("Plan has no required segments.", PlanWarning, stacklevel=2)

        return plan

    @model_validator(mode="after")
    def check_excluded_accessions(self) -> "OTU":
        """Ensure that excluded accessions are not in the OTU."""
        if accessions := self.excluded_accessions & {
            sequence.accession.key for sequence in self.sequences
        }:
            raise ValueError(
                f"Excluded accessions found in the OTU: {', '.join(accessions)}"
            )

        return self

    @model_validator(mode="after")
    def check_representative_isolate(self) -> "OTU":
        """Ensure that the default isolate is in the OTU."""
        if self.representative_isolate not in {isolate.id for isolate in self.isolates}:
            raise ValueError("Representative isolate must be in the OTU")

        return self

    @model_validator(mode="after")
    def check_unique_isolate_names(self) -> "OTU":
        """Ensure there are no duplicate isolate names in the OTU."""
        counts = Counter(isolate.name for isolate in self.isolates)

        duplicates = ", ".join(
            [str(name) for name, count in counts.items() if name and count > 1]
        )

        if duplicates:
            raise ValueError(
                f"Isolate names must be unique. Non-unique names: {duplicates}",
            )

        return self

    @model_validator(mode="after")
    def check_isolates_against_plan(self) -> "OTU":
        """Check that all isolates satisfy the OTU's plan."""
        for isolate in self.isolates:
            for sequence in isolate.sequences:
                segment = self.plan.get_segment_by_id(sequence.segment)
                if segment is None:
                    raise PydanticCustomError(
                        "segment_not_found",
                        "Sequence segment {sequence_segment} was not found in "
                        + "the list of segments: {plan_segments}.",
                        {
                            "isolate_id": isolate.id,
                            "sequence_segment": sequence.segment,
                            "plan_segments": list(self.plan.segment_ids),
                        },
                    )

                min_length = int(segment.length * (1.0 - segment.length_tolerance))
                max_length = int(segment.length * (1.0 + segment.length_tolerance))

                if len(sequence.sequence) < min_length:
                    raise PydanticCustomError(
                        "sequence_too_short",
                        "Sequence based on {sequence_accession} does not pass validation "
                        + "against segment {segment_id} "
                        + "({sequence_length} < {min_sequence_length})",
                        {
                            "isolate_id": isolate.id,
                            "sequence_id": sequence.id,
                            "sequence_accession": sequence.accession,
                            "sequence_length": len(sequence.sequence),
                            "segment_id": segment.id,
                            "segment_reference_length": segment.length,
                            "min_sequence_length": min_length,
                        },
                    )

                if len(sequence.sequence) > max_length:
                    raise PydanticCustomError(
                        "sequence_too_long",
                        "Sequence based on {sequence_accession} does not pass validation "
                        + "against segment {segment_id}"
                        + "({sequence_length} > {max_sequence_length})",
                        {
                            "isolate_id": isolate.id,
                            "sequence_id": sequence.id,
                            "sequence_accession": sequence.accession,
                            "sequence_length": len(sequence.sequence),
                            "segment_id": segment.id,
                            "segment_reference_length": segment.length,
                            "max_sequence_length": max_length,
                        },
                    )

        return self
