from pydantic import UUID4, ConfigDict, field_serializer

from ref_builder.events.base import ApplicableEvent, Event, EventData, OTUQuery
from ref_builder.models import Molecule
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.plan import Plan
from ref_builder.utils import ExcludedAccessionAction


class CreateOTUData(EventData):
    """The data for a :class:`CreateOTU` event."""

    id: UUID4
    acronym: str
    legacy_id: str | None
    molecule: Molecule
    name: str
    taxid: int
    plan: Plan


class CreateOTU(Event):
    """An event that creates a new OTU."""

    data: CreateOTUData
    query: OTUQuery

    def apply(self) -> OTUBuilder:
        return OTUBuilder(
            id=self.query.otu_id,
            acronym=self.data.acronym,
            excluded_accessions=set(),
            isolates=[],
            legacy_id=self.data.legacy_id,
            molecule=self.data.molecule,
            name=self.data.name,
            representative_isolate=None,
            plan=self.data.plan,
            taxid=self.data.taxid,
        )


class CreatePlanData(EventData):
    """The data for a :class:`CreatePlan` event."""

    plan: Plan


class CreatePlan(ApplicableEvent):
    """An event that sets the isolate plan for an OTU."""

    data: CreatePlanData
    query: OTUQuery

    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        """Apply changed plan to OTU and return."""
        otu.plan = self.data.plan

        return otu


class SetRepresentativeIsolateData(EventData):
    """The data for a :class:`SetReprIsolate` event."""

    isolate_id: UUID4


class SetRepresentativeIsolate(ApplicableEvent):
    """An event that sets the representative isolate for an OTU."""

    data: SetRepresentativeIsolateData
    query: OTUQuery

    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        """Update the OTU's representative isolate and return."""
        otu.representative_isolate = self.data.isolate_id

        return otu


class UpdateExcludedAccessionsData(EventData):
    """The data for an UpdateAllowedAccessions event."""

    model_config = ConfigDict(use_enum_values=True)

    accessions: set[str]
    action: ExcludedAccessionAction

    @field_serializer("accessions")
    def serialize_accessions(self, accessions: set[str]) -> list[str]:
        return sorted(accessions)


class UpdateExcludedAccessions(ApplicableEvent):
    """An event that changes the OTU excluded accessions collection.

    This event is emitted when Genbank accessions are either
    allowed or disallowed from inclusion in the reference.
    """

    data: UpdateExcludedAccessionsData
    query: OTUQuery

    def apply(self, otu: OTUBuilder) -> OTUBuilder:
        """Add accession allowance changes to OTU and return."""
        if self.data.action == ExcludedAccessionAction.ALLOW:
            for accession in self.data.accessions:
                otu.excluded_accessions.discard(accession)

        else:
            for accession in self.data.accessions:
                otu.excluded_accessions.add(accession)

        return otu
