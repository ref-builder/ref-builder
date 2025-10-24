from pydantic import UUID4, ConfigDict, field_serializer

from ref_builder.events.base import ApplicableEvent, Event, EventData, OTUQuery
from ref_builder.models.lineage import Lineage
from ref_builder.models.molecule import Molecule
from ref_builder.models.plan import Plan
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.utils import ExcludedAccessionAction


class CreateOTUData(EventData):
    """The data for a :class:`CreateOTU` event."""

    id: UUID4
    lineage: Lineage
    molecule: Molecule
    taxid: int
    plan: Plan


class CreateOTU(Event[CreateOTUData, OTUQuery]):
    """An event that creates a new OTU."""

    def apply(self) -> OTUBuilder:
        return OTUBuilder(
            id=self.query.otu_id,
            excluded_accessions=set(),
            isolates=[],
            lineage=self.data.lineage,
            molecule=self.data.molecule,
            plan=self.data.plan,
            taxid=self.data.taxid,
        )


class CreatePlanData(EventData):
    """The data for a :class:`CreatePlan` event."""

    plan: Plan


class CreatePlan(ApplicableEvent[CreatePlanData, OTUQuery]):
    """An event that sets the isolate plan for an OTU."""

    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        """Apply changed plan to OTU and return."""
        otu.plan = self.data.plan

        return otu


class UpdateExcludedAccessionsData(EventData):
    """The data for an UpdateAllowedAccessions event."""

    model_config = ConfigDict(use_enum_values=True)

    accessions: set[str]
    action: ExcludedAccessionAction

    @field_serializer("accessions")
    def serialize_accessions(self, accessions: set[str]) -> list[str]:
        return sorted(accessions)


class UpdateExcludedAccessions(ApplicableEvent[UpdateExcludedAccessionsData, OTUQuery]):
    """An event that changes the OTU excluded accessions collection.

    This event is emitted when Genbank accessions are either
    allowed or disallowed from inclusion in the reference.
    """

    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        """Add accession allowance changes to OTU and return."""
        if self.data.action == ExcludedAccessionAction.ALLOW:
            for accession in self.data.accessions:
                otu.excluded_accessions.discard(accession)

        else:
            for accession in self.data.accessions:
                otu.excluded_accessions.add(accession)

        return otu
