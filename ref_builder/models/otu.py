import math
import warnings
from uuid import UUID

from pydantic import (
    UUID4,
    BaseModel,
    ConfigDict,
    Field,
    PrivateAttr,
    field_validator,
    model_validator,
)
from pydantic_core import PydanticCustomError

from ref_builder.models.accession import Accession
from ref_builder.models.isolate import Isolate
from ref_builder.models.lineage import Lineage
from ref_builder.models.molecule import Molecule
from ref_builder.models.plan import Plan
from ref_builder.models.sequence import Sequence
from ref_builder.warnings import PlanWarning


class OTUMinimal(BaseModel):
    """A minimal representation of an OTU."""

    acronym: str
    id: UUID
    name: str
    taxid: int


class OTU(BaseModel):
    """A validated OTU."""

    model_config = ConfigDict(validate_assignment=True, revalidate_instances="always")

    _isolates_by_accession: dict[str, Isolate] = PrivateAttr()
    """A dictionary of isolates indexed by accession key"""

    _isolates_by_id: dict[UUID4, Isolate] = PrivateAttr()
    """A dictionary of isolates indexed by isolate UUID"""

    _sequences_by_accession: dict[str, Sequence] = PrivateAttr()
    """A dictionary of sequences indexed by accession key"""

    id: UUID4
    """The OTU id."""

    excluded_accessions: set[str]
    """A set of accessions that should not be retrieved in future fetch operations."""

    promoted_accessions: set[str]
    """GenBank accessions that have been replaced by RefSeq equivalents during promotion."""

    lineage: Lineage
    """The taxonomic lineage from species down to the target taxon."""

    molecule: Molecule
    """The type of molecular information contained in this OTU."""

    isolates: list[Isolate] = Field(min_length=1)
    """Isolate in the OTU."""

    excluded_isolates: list[Isolate] = Field(default_factory=list)
    """Isolates that contain excluded accessions."""

    plan: Plan
    """The plan for the OTU."""

    @model_validator(mode="after")
    def rebuild_lookups(self) -> "OTU":
        """Rebuild lookup dictionaries when isolates change."""
        self._isolates_by_accession = {
            accession: isolate
            for isolate in self.isolates
            for accession in isolate.accessions
        }
        self._isolates_by_id = {isolate.id: isolate for isolate in self.isolates}
        self._sequences_by_accession = {
            sequence.accession.key: sequence
            for isolate in self.isolates
            for sequence in isolate.sequences
        }
        return self

    def get_isolate(self, isolate_id: UUID4) -> Isolate | None:
        """Get isolate associated with a given ID.

        Returns None if no such isolate exists.

        :param isolate_id: The UUID of the isolate to retrieve
        :return: the isolate or ``None``
        """
        return self._isolates_by_id.get(isolate_id)

    def get_isolate_by_accession(self, accession: str) -> Isolate | None:
        """Get the isolate containing a sequence with the given accession.

        Returns None if no such isolate exists.

        :param accession: The accession key to search for
        :return: The isolate containing the accession or ``None``
        """
        return self._isolates_by_accession.get(accession)

    def get_sequence(self, accession: str) -> Sequence | None:
        """Get the sequence with the given accession key.

        Returns None if no such sequence exists.

        :param accession: The accession key to search for
        :return: The sequence or ``None``
        """
        return self._sequences_by_accession.get(accession)

    @property
    def accessions(self) -> set[str]:
        """A set of accessions contained in this isolate."""
        return {accession.key for accession in self.versioned_accessions}

    @property
    def acronym(self) -> str:
        """The OTU acronym computed from the lineage."""
        return self.lineage.acronym

    @property
    def name(self) -> str:
        """The OTU name computed from the species-level taxon in the lineage."""
        return self.lineage.name

    @property
    def blocked_accessions(self) -> set[str]:
        """Accessions that should not be considered for addition to the OTU.

        This includes:
        - Accessions that already exist in the OTU
        - Accessions that have been explicitly excluded
        - Accessions that have been promoted (replaced by RefSeq equivalents)
        """
        return self.accessions | self.excluded_accessions | self.promoted_accessions

    @property
    def isolate_ids(self) -> set[UUID4]:
        """A set of UUIDs for isolates in the OTU."""
        return set(self._isolates_by_id.keys())

    @property
    def sequences(self) -> list[Sequence]:
        """Sequences contained in this OTU."""
        return [sequence for isolate in self.isolates for sequence in isolate.sequences]

    @property
    def synonyms(self) -> set[str]:
        """All possible names derived from the OTU's lineage.

        Returns all taxon names, acronyms, and synonyms from the lineage.
        """
        return self.lineage.synonyms

    @property
    def taxid(self) -> int:
        """The NCBI Taxonomy ID for this OTU (species-level taxon from lineage)."""
        return self.lineage.taxa[0].id

    @property
    def versioned_accessions(self) -> set[Accession]:
        """A set of versioned accessions contained in this OTU."""
        return {
            sequence.accession
            for isolate in self.isolates
            for sequence in isolate.sequences
        }

    @field_validator("plan", mode="after")
    def check_plan_required(cls, plan: Plan) -> Plan:
        """Issue a warning if the plan has no required segments."""
        if not plan.required_segments:
            warnings.warn("Plan has no required segments.", PlanWarning, stacklevel=2)

        return plan

    @model_validator(mode="after")
    def check_promoted_accessions(self) -> "OTU":
        """Ensure that promoted accessions are not in the OTU or excluded."""
        if accessions := self.promoted_accessions & {
            sequence.accession.key for sequence in self.sequences
        }:
            raise ValueError(
                f"Promoted accessions found in the OTU: {', '.join(accessions)}"
            )

        if accessions := self.promoted_accessions & self.excluded_accessions:
            raise ValueError(
                f"Promoted accessions cannot be excluded: {', '.join(accessions)}"
            )

        return self

    @model_validator(mode="after")
    def check_isolate_taxids(self) -> "OTU":
        """Ensure that all isolate taxids are in the OTU's lineage."""
        lineage_taxids = {taxon.id for taxon in self.lineage.taxa}

        for isolate in self.isolates:
            if isolate.taxid not in lineage_taxids:
                raise PydanticCustomError(
                    "isolate_taxid_not_in_lineage",
                    "Isolate taxid {isolate_taxid} is not found in OTU lineage. "
                    "Valid taxids: {lineage_taxids}",
                    {
                        "isolate_id": isolate.id,
                        "isolate_taxid": isolate.taxid,
                        "lineage_taxids": sorted(lineage_taxids),
                    },
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
                        "the list of segments: {plan_segments}.",
                        {
                            "isolate_id": isolate.id,
                            "sequence_segment": sequence.segment,
                            "plan_segments": list(self.plan.segment_ids),
                        },
                    )

                min_length = math.floor(
                    segment.length * (1.0 - segment.length_tolerance)
                )
                max_length = math.ceil(
                    segment.length * (1.0 + segment.length_tolerance)
                )

                if len(sequence.sequence) < min_length:
                    raise PydanticCustomError(
                        "sequence_too_short",
                        "Sequence based on {sequence_accession} does not pass validation "
                        "against segment {segment_id} "
                        "({sequence_length} < {min_sequence_length})",
                        {
                            "isolate_id": isolate.id,
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
                        "against segment {segment_id}"
                        "({sequence_length} > {max_sequence_length})",
                        {
                            "isolate_id": isolate.id,
                            "sequence_accession": sequence.accession,
                            "sequence_length": len(sequence.sequence),
                            "segment_id": segment.id,
                            "segment_reference_length": segment.length,
                            "max_sequence_length": max_length,
                        },
                    )

        return self
