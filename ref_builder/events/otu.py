from pydantic import UUID4, ConfigDict, field_serializer

from ref_builder.events.base import (
    ApplicableEvent,
    Event,
    EventData,
    OTUQuery,
)
from ref_builder.events.isolate import CreateIsolateData
from ref_builder.models.isolate import Isolate
from ref_builder.models.lineage import Lineage
from ref_builder.models.molecule import Molecule
from ref_builder.models.otu import OTU
from ref_builder.models.plan import Plan
from ref_builder.utils import ExcludedAccessionAction


class CreateOTUData(EventData):
    """The data for a :class:`CreateOTU` event."""

    id: UUID4
    excluded_accessions: set[str]
    isolate: CreateIsolateData
    lineage: Lineage
    molecule: Molecule
    plan: Plan


class CreateOTU(Event[CreateOTUData, OTUQuery]):
    """An event that creates a new OTU."""

    def apply(self) -> OTU:
        """Apply and OTU creation event.

        Instantiates and returns and OTU builder.
        """
        return OTU(
            id=self.query.otu_id,
            excluded_accessions=self.data.excluded_accessions,
            isolates=[
                Isolate(
                    id=self.data.isolate.id,
                    name=self.data.isolate.name,
                    sequences=self.data.isolate.sequences,
                    taxid=self.data.isolate.taxid,
                )
            ],
            lineage=self.data.lineage,
            molecule=self.data.molecule,
            plan=self.data.plan,
        )


class SetPlanData(EventData):
    """The data for a :class:`SetPlan` event."""

    plan: Plan


class SetPlan(ApplicableEvent[SetPlanData, OTUQuery]):
    """An event that sets the isolate plan for an OTU."""

    def apply(self, otu: OTU) -> OTU:
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

    def apply(self, otu: OTU) -> OTU:
        """Add accession allowance changes to OTU and return."""
        if self.data.action == ExcludedAccessionAction.ALLOW:
            for accession in self.data.accessions:
                otu.excluded_accessions.discard(accession)

        else:
            for accession in self.data.accessions:
                otu.excluded_accessions.add(accession)

        return otu
