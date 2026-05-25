"""Lightweight, validator-bypassing OTU snapshots used by audit tooling.

The OTU model now refuses to construct when it contains duplicate accessions
within an isolate. That guarantees future repos cannot be corrupted in that
way, but it also means an already-corrupt repo can no longer be loaded via
the normal hydration path. The audit machinery here walks events directly via
``Repo.iter_otu_audit_snapshots`` so existing corruption can still be surfaced.
"""

from collections import defaultdict
from dataclasses import dataclass, field
from typing import TYPE_CHECKING
from uuid import UUID

if TYPE_CHECKING:
    from ref_builder.repo import Repo

MIN_OCCURRENCES_FOR_DUPLICATE = 2


@dataclass
class IsolateAuditSnapshot:
    """Per-isolate audit data; ``accession_keys`` is a list so duplicates persist."""

    id: UUID
    name: str | None
    accession_keys: list[str] = field(default_factory=list)


@dataclass
class OTUAuditSnapshot:
    """Per-OTU audit data extracted directly from events."""

    id: UUID
    taxid: int
    name: str
    isolates: list[IsolateAuditSnapshot]


@dataclass
class AccessionPlacement:
    """Where a single accession occurrence lives."""

    otu_id: UUID
    otu_name: str
    otu_taxid: int
    isolate_id: UUID
    isolate_name: str | None


@dataclass
class AccessionAuditReport:
    """Findings from a repo-wide accession scan."""

    cross_otu: dict[str, list[AccessionPlacement]]
    within_otu: dict[str, list[AccessionPlacement]]

    @property
    def has_findings(self) -> bool:
        """True if either kind of duplicate was found."""
        return bool(self.cross_otu or self.within_otu)


def audit_accessions(repo: "Repo") -> AccessionAuditReport:
    """Scan a repo for cross-OTU and within-OTU duplicate accessions."""
    placements: dict[str, list[AccessionPlacement]] = defaultdict(list)

    for otu in repo.iter_otu_audit_snapshots():
        for isolate in otu.isolates:
            for accession_key in isolate.accession_keys:
                placements[accession_key].append(
                    AccessionPlacement(
                        otu_id=otu.id,
                        otu_name=otu.name,
                        otu_taxid=otu.taxid,
                        isolate_id=isolate.id,
                        isolate_name=isolate.name,
                    )
                )

    cross_otu: dict[str, list[AccessionPlacement]] = {}
    within_otu: dict[str, list[AccessionPlacement]] = {}

    for accession_key, occurrences in placements.items():
        if len(occurrences) < MIN_OCCURRENCES_FOR_DUPLICATE:
            continue

        otu_ids = {placement.otu_id for placement in occurrences}
        if len(otu_ids) > 1:
            cross_otu[accession_key] = occurrences
        else:
            within_otu[accession_key] = occurrences

    return AccessionAuditReport(cross_otu=cross_otu, within_otu=within_otu)
