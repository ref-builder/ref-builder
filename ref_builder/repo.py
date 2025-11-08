"""An access layer for a Virtool event-sourced reference repository.

This is a work in progress.


TODO: Check if excluded accessions exist in the repo.
TODO: Check for accessions filed in wrong isolates.
TODO: Check for accession conflicts.

"""

import datetime
import shutil
import uuid
import warnings
from collections import defaultdict
from collections.abc import Collection, Generator, Iterator
from contextlib import contextmanager
from pathlib import Path

import arrow
from structlog import get_logger

from ref_builder.errors import (
    OTUDeletedError,
)
from ref_builder.events.base import (
    ApplicableEvent,
    Event,
    EventData,
    EventMetadata,
    EventQuery,
    IsolateQuery,
    OTUQuery,
    RepoQuery,
    SequenceQuery,
)
from ref_builder.events.isolate import (
    CreateIsolate,
    CreateIsolateData,
    DeleteIsolate,
    DeleteIsolateData,
    PromoteIsolate,
    PromoteIsolateData,
)
from ref_builder.events.otu import (
    CreateOTU,
    CreateOTUData,
    SetPlan,
    SetPlanData,
    UpdateExcludedAccessions,
    UpdateExcludedAccessionsData,
)
from ref_builder.events.repo import (
    CreateRepo,
    CreateRepoData,
)
from ref_builder.events.sequence import UpdateSequence, UpdateSequenceData
from ref_builder.index import Index
from ref_builder.lock import Lock
from ref_builder.models.accession import Accession
from ref_builder.models.isolate import Isolate, IsolateName
from ref_builder.models.lineage import Lineage
from ref_builder.models.molecule import Molecule
from ref_builder.models.otu import OTU, OTUMinimal
from ref_builder.models.plan import Plan
from ref_builder.models.repo import RepoMeta, RepoSettings
from ref_builder.models.sequence import Sequence
from ref_builder.store import EventStore
from ref_builder.utils import (
    ExcludedAccessionAction,
    get_accession_key,
)
from ref_builder.warnings import OTUDeletedWarning

GITIGNORE_CONTENTS = [".cache", "lock"]

logger = get_logger("repo")


class Repo:
    """An event-sourced repository."""

    def __init__(self, path: Path) -> None:
        """Create a new instance of the repository."""
        self.path = path
        """The path to the repo directory."""

        self._event_store = EventStore(self.path)
        """The event store of the event sourced repository."""

        self._index = Index(self.path / ".cache" / "index.db")
        """An index for fast lookups.

        It allows fast lookup of OTUs be key fields, fetching of complete OTU state,
        and the events associated with a given OTU ID.
        """

        self._lock = Lock(self.path)
        """A lock for the repository."""

        # Populate the index if it is empty.
        if not self._index.otu_ids:
            self.rebuild_index()

    @classmethod
    def new(
        cls,
        name: str,
        path: Path,
        organism: str,
        default_segment_length_tolerance: float = 0.03,
    ) -> "Repo":
        """Create a new reference repository."""
        if path.is_file():
            raise ValueError("The target path is a file")

        path.mkdir(parents=True, exist_ok=True)

        if any(path.iterdir()):
            raise ValueError("The target path is not empty")

        with open(path / ".gitignore", "w") as f:
            f.write("\n".join(GITIGNORE_CONTENTS) + "\n")

        shutil.copytree(
            Path(__file__).parent.parent / "assets/github",
            path / ".github",
        )

        (path / ".cache").mkdir()

        repo_id = uuid.uuid4()

        EventStore(path).write_event(
            CreateRepo(
                id=1,
                data=CreateRepoData(
                    id=repo_id,
                    name=name,
                    organism=organism,
                    settings=RepoSettings(
                        default_segment_length_tolerance=default_segment_length_tolerance
                    ),
                ),
                query=RepoQuery(repository_id=repo_id),
                timestamp=arrow.utcnow().naive,
            )
        )

        return Repo(path)

    @property
    def last_id(self) -> int:
        """The id of the most recently added event in the event store."""
        return self._event_store.last_id

    @property
    def meta(self) -> RepoMeta:
        """The metadata for the repository."""
        for event in self._event_store.iter_events():
            if isinstance(event, CreateRepo):
                return RepoMeta(**event.data.model_dump(), created_at=event.timestamp)

        raise ValueError("No repository creation event found")

    @property
    def settings(self) -> RepoSettings:
        """The settings for the repository."""
        for event in self._event_store.iter_events():
            if isinstance(event, CreateRepo):
                return event.data.settings

        raise ValueError("No repository creation event found")

    @contextmanager
    def lock(self) -> Iterator[None]:
        """Lock the repository.

        This prevents read and write access from  other ``ref-builder`` processes.
        """
        self._lock.lock()

        try:
            yield
        finally:
            self._lock.unlock()

    def clear_index(self) -> bool:
        """Delete and replace the repository read index."""
        index_path = self._index.path

        if index_path.exists():
            index_path.unlink()

            (index_path.parent / f"{index_path.stem}.db-shm").unlink()
            (index_path.parent / f"{index_path.stem}.db-wal").unlink()

            return True

        return False

    def rebuild_index(self) -> None:
        """Rebuild the repository read index."""
        logger.info(
            "No repo index found. Rebuilding...",
            path=str(self.path),
        )

        for otu in self.iter_otus_from_events():
            self._index.upsert_otu(otu, self.last_id)

        for event in self._event_store.iter_events():
            try:
                otu_id = event.query.model_dump()["otu_id"]
            except KeyError:
                continue

            self._index.add_event_id(event.id, otu_id, event.timestamp)

    def iter_minimal_otus(self) -> Iterator[OTUMinimal]:
        """Iterate over minimal representations of the OTUs in the repository.

        This is more performant than iterating over full OTUs.
        """
        return self._index.iter_minimal_otus()

    def iter_otus(self) -> Iterator[OTU]:
        """Iterate over the OTUs in the repository."""
        for otu_id in self._index.otu_ids:
            if (otu := self.get_otu(otu_id)) is not None:
                yield otu

    def iter_otus_from_events(self) -> Iterator[OTU]:
        """Iterate over the OTUs, bypassing the index."""
        event_ids_by_otu = defaultdict(list)

        for event in self._event_store.iter_events():
            if hasattr(event.query, "otu_id"):
                event_ids_by_otu[event.query.otu_id].append(event.id)

        for otu_id in event_ids_by_otu:
            try:
                yield self._rehydrate_otu(
                    self._event_store.read_event(event_id)
                    for event_id in event_ids_by_otu[otu_id]
                )
            except OTUDeletedError:
                continue

    def create_otu(
        self,
        isolate: CreateIsolateData,
        lineage: Lineage,
        molecule: Molecule,
        plan: Plan,
        promoted_accessions: set[str],
    ) -> OTU | None:
        """Create an OTU."""
        for taxon in lineage.taxa:
            if (otu_id := self.get_otu_id_by_taxid(taxon.id)) is not None:
                otu = self.get_otu(otu_id)
                msg = f"OTU {otu.id} already contains taxid {taxon.id} ({taxon.name})"
                raise ValueError(msg)

        taxid = lineage.taxa[0].id

        logger.info("Creating new OTU.", taxid=taxid, name=lineage.name)

        otu_id = uuid.uuid4()

        self._write_event(
            CreateOTU,
            CreateOTUData(
                id=otu_id,
                isolate=isolate,
                lineage=lineage,
                molecule=molecule,
                plan=plan,
                promoted_accessions=promoted_accessions,
            ),
            OTUQuery(otu_id=otu_id),
        )

        return self.get_otu(otu_id)

    def create_isolate(
        self,
        otu_id: uuid.UUID,
        isolate_id: uuid.UUID,
        name: IsolateName | None,
        taxid: int,
        sequences: list,
    ) -> Isolate:
        """Create an isolate with sequences atomically.

        :param otu_id: the OTU ID
        :param isolate_id: the isolate ID
        :param name: the isolate name
        :param taxid: the taxid
        :param sequences: list of SequenceData objects
        :return: the created isolate
        """
        event = self._write_event(
            CreateIsolate,
            CreateIsolateData(
                id=isolate_id,
                name=name,
                sequences=sequences,
                taxid=taxid,
            ),
            IsolateQuery(isolate_id=isolate_id, otu_id=otu_id),
        )

        logger.debug(
            "Isolate with sequences written",
            event_id=event.id,
            isolate_id=str(isolate_id),
            name=str(name) if name is not None else None,
            taxid=taxid,
            sequence_count=len(sequences),
        )

        return self.get_otu(otu_id).get_isolate(isolate_id)

    def delete_isolate(
        self,
        otu_id: uuid.UUID,
        isolate_id: uuid.UUID,
        message: str,
    ) -> None:
        """Delete an existing isolate from a given OTU."""
        if (otu_ := self.get_otu(otu_id)) is None:
            raise ValueError(f"OTU does not exist: {otu_id}")

        if isolate_id not in otu_.isolate_ids:
            raise ValueError(f"Isolate does not exist: {isolate_id}")

        self._write_event(
            DeleteIsolate,
            DeleteIsolateData(message=message),
            IsolateQuery(otu_id=otu_id, isolate_id=isolate_id),
        )

    def promote_isolate(
        self,
        otu_id: uuid.UUID,
        isolate_id: uuid.UUID,
        accession_map: dict,
    ) -> None:
        """Promote GenBank sequences to RefSeq in an isolate.

        All sequences in an isolate must be promoted together to maintain consistency.
        If one sequence is RefSeq, all must be RefSeq.

        :param otu_id: The OTU ID
        :param isolate_id: The isolate ID
        :param accession_map: A dict mapping old accessions to new SequenceData
        """
        otu = self.get_otu(otu_id)

        if otu is None:
            raise ValueError(f"OTU does not exist: {otu_id}")

        isolate = otu.get_isolate(isolate_id)

        if isolate is None:
            raise ValueError(f"Isolate does not exist: {isolate_id}")

        self._write_event(
            PromoteIsolate,
            PromoteIsolateData(map=accession_map),
            IsolateQuery(otu_id=otu_id, isolate_id=isolate_id),
        )

    def update_sequence(
        self,
        otu_id: uuid.UUID,
        old_accession: Accession,
        new_sequence: Sequence,
    ) -> None:
        """Update a sequence to a newer version.

        This replaces the old sequence with the new sequence across all isolates
        where the old sequence is linked and deletes the old sequence. Does not
        exclude accessions since the accession key remains the same.
        """
        self._write_event(
            UpdateSequence,
            UpdateSequenceData(sequence=new_sequence),
            SequenceQuery(otu_id=otu_id, accession=old_accession),
        )

    def set_plan(self, otu_id: uuid.UUID, plan: Plan) -> Plan:
        """Set the isolate plan for an OTU."""
        self._write_event(
            SetPlan,
            SetPlanData(plan=plan),
            OTUQuery(otu_id=otu_id),
        )

        return self.get_otu(otu_id).plan

    def exclude_accessions(
        self,
        otu_id: uuid.UUID,
        accessions: Collection[str],
    ) -> set[str]:
        """Add accessions to OTU's excluded accessions."""
        otu = self.get_otu(otu_id)

        try:
            excludable_accessions = {
                get_accession_key(raw_accession) for raw_accession in accessions
            }
        except ValueError as e:
            if "Invalid accession key" in str(e):
                logger.warning(
                    "Invalid accession included in set. "
                    "No changes were made to excluded accessions.",
                    accessions=sorted(accessions),
                )

                return otu.excluded_accessions

            raise

        if unremovable_accessions := excludable_accessions.intersection(otu.accessions):
            logger.warning(
                "Accessions currently in OTU cannot be removed.",
                unremovable_accessions=sorted(unremovable_accessions),
            )
            excludable_accessions -= unremovable_accessions

        if (
            extant_requested_accessions := excludable_accessions
            & otu.excluded_accessions
        ):
            logger.info(
                "Ignoring already excluded accessions",
                requested_exclusions=sorted(extant_requested_accessions),
                old_excluded_accessions=sorted(otu.excluded_accessions),
            )

            excludable_accessions -= otu.excluded_accessions

        if excludable_accessions:
            self._write_event(
                UpdateExcludedAccessions,
                UpdateExcludedAccessionsData(
                    accessions=excludable_accessions,
                    action=ExcludedAccessionAction.EXCLUDE,
                ),
                OTUQuery(otu_id=otu_id),
            )

            logger.info(
                "Added accessions to excluded accession list.",
                taxid=otu.taxid,
                otu_id=str(otu.id),
                new_excluded_accessions=sorted(excludable_accessions),
                old_excluded_accessions=sorted(otu.excluded_accessions),
            )
        else:
            logger.warning("No excludable accessions were given.")

        return self.get_otu(otu_id).excluded_accessions

    def allow_accessions(
        self,
        otu_id: uuid.UUID,
        accessions: Collection[str],
    ) -> set[str]:
        """Remove accessions from OTU's excluded accessions."""
        otu = self.get_otu(otu_id)

        allowable_accessions = set(accessions)

        if redundant_accessions := allowable_accessions - otu.excluded_accessions:
            logger.debug(
                "Ignoring non-excluded accessions",
                non_excluded_accessions=sorted(redundant_accessions),
            )

            allowable_accessions = allowable_accessions - redundant_accessions

        if allowable_accessions:
            self._write_event(
                UpdateExcludedAccessions,
                UpdateExcludedAccessionsData(
                    accessions=set(allowable_accessions),
                    action=ExcludedAccessionAction.ALLOW,
                ),
                OTUQuery(otu_id=otu_id),
            )

            logger.info(
                "Removed accessions from excluded accession list.",
                taxid=otu.taxid,
                otu_id=str(otu.id),
                new_excluded_accessions=sorted(allowable_accessions),
            )

        return self.get_otu(otu_id).excluded_accessions

    def get_otu_id_by_isolate_id(self, isolate_id: uuid.UUID) -> uuid.UUID | None:
        """Get an OTU ID from an isolate ID that belongs to it."""
        return self._index.get_id_by_isolate_id(isolate_id)

    def get_otu(self, otu_id: uuid.UUID) -> OTU | None:
        """Get the OTU with the given ``otu_id``.

        If the OTU does not exist, ``None`` is returned.

        :param otu_id: the id of the OTU
        :return: the OTU or ``None``

        """
        event_index_item = self._index.get_event_ids_by_otu_id(otu_id)

        if event_index_item is None:
            return None

        try:
            events = (
                self._event_store.read_event(event_id)
                for event_id in event_index_item.event_ids
            )

            otu = self._rehydrate_otu(events)

        except OTUDeletedError:
            warnings.warn(
                f"OTU {otu_id} has already been deleted.",
                category=OTUDeletedWarning,
                stacklevel=1,
            )

            return None

        except FileNotFoundError:
            logger.error("Event exists in index, but not in source. Deleting index...")

            self.clear_index()

            raise

        self._index.upsert_otu(otu, self.last_id)

        return otu

    def iter_otu_events(self, otu_id: uuid.UUID) -> Generator[ApplicableEvent]:
        """Iterate through event log."""
        event_index_item = self._index.get_event_ids_by_otu_id(otu_id)

        if event_index_item is not None:
            for event_id in event_index_item.event_ids:
                yield self._event_store.read_event(event_id)

    def iter_event_metadata(self) -> Generator[EventMetadata]:
        """Iterate through the event metadata of all events."""
        yield from self._index.iter_event_metadata()

    def get_otu_by_taxid(self, taxid: int) -> OTU | None:
        """Return the OTU with the given ``taxid``.

        If no OTU is found, return None.

        :param taxid: the taxonomy ID of the OTU
        :return: the matching OTU or ``None``

        """
        if (otu_id := self.get_otu_id_by_taxid(taxid)) is not None:
            return self.get_otu(otu_id)

        return None

    def get_otu_id_by_taxid(self, taxid: int) -> uuid.UUID | None:
        """Return the UUID of the OTU with the given ``taxid``.

        If no OTU is found, return None.

        :param taxid: the taxonomy ID of the OTU
        :return: the UUID of the OTU or ``None``

        """
        return self._index.get_id_by_taxid(taxid)

    def get_isolate(self, isolate_id: uuid.UUID) -> Isolate | None:
        """Return the isolate with the given id if it exists, else None."""
        if otu_id := self.get_otu_id_by_isolate_id(isolate_id):
            return self.get_otu(otu_id).get_isolate(isolate_id)

        return None

    def get_otu_first_created(self, otu_id: uuid.UUID) -> datetime.datetime | None:
        """Get the timestamp of the first event associated with an OTU.
        If no events can be found for this OTU, return None.
        """
        return self._index.get_first_timestamp_by_otu_id(otu_id)

    def get_otu_last_modified(self, otu_id: uuid.UUID) -> datetime.datetime | None:
        """Get the timestamp of the last event associated with an OTU.
        If no events can be found for this OTU, return None.
        """
        return self._index.get_latest_timestamp_by_otu_id(otu_id)

    def get_otu_last_updated(self, otu_id: uuid.UUID) -> datetime.datetime | None:
        """Get the timestamp of the last time this OTU was automatically updated.
        If this OTU has not been updated since this repo was initialized, return None.
        """
        return self._index.get_last_otu_update_timestamp(otu_id)

    def write_otu_update_history_entry(self, otu_id: uuid.UUID) -> int:
        """Add a new entry to the otu update history log and return the primary key
        of the entry.
        """
        if (
            update_id := self._index.add_otu_update_history_entry(
                otu_id,
                arrow.utcnow().naive,
            )
            is None
        ):
            raise SystemError(
                "OTU update history entry could not be retrieved after writing."
            )

        return update_id

    def get_event(self, event_id: int) -> Event | None:
        """Return event from event store."""
        try:
            return self._event_store.read_event(event_id)

        except FileNotFoundError:
            return None

    @staticmethod
    def _rehydrate_otu(events: Iterator[Event]) -> OTU:
        """Rehydrate an OTU from an event iterator."""
        event = next(events)

        with warnings.catch_warnings(record=True) as warning_list:
            if isinstance(event, CreateOTU):
                otu = event.apply()

            else:
                raise TypeError(
                    f"The first event ({event}) for an OTU is not a CreateOTU event",
                )

            for event in events:
                if not isinstance(event, ApplicableEvent):
                    raise TypeError(
                        f"Event {event.id} {event.type} is not an applicable event."
                    )

                otu = event.apply(otu)
                otu = OTU.model_validate(otu)

        for warning_msg in warning_list:
            logger.warning(
                warning_msg.message,
                otu_id=str(otu.id),
                warning_category=warning_msg.category.__name__,
            )

        otu.isolates.sort(
            key=lambda i: f"{i.name.type} {i.name.value}"
            if type(i.name) is IsolateName
            else "",
        )

        for isolate in otu.isolates:
            isolate.sequences.sort(key=lambda s: s.accession)

        # Validate to rebuild OTU lookup dictionaries after all mutations
        return OTU.model_validate(otu)

    def _write_event(
        self, cls: type[Event], data: EventData, query: EventQuery
    ) -> Event:
        """Write an event after validating it would produce a valid state.

        For OTU events, validation is performed by applying the event to an in-memory
        copy of the OTU. If the event would produce an invalid OTU, a ValueError is
        raised and no event is written.

        :param cls: The event class
        :param data: The event data
        :param query: The event query
        :return: The written event
        :raises ValueError: If the event would create invalid state
        """
        # Create event object (not yet written to disk)
        event = cls(
            id=self.last_id + 1,
            data=data,
            query=query,
            timestamp=arrow.utcnow().naive,
        )

        # Validate OTU events by applying in-memory
        if hasattr(event.query, "otu_id"):
            # CreateOTU is always valid (creates new OTU with validated data)
            if not isinstance(event, CreateOTU):
                # Get current OTU state
                otu = self.get_otu(event.query.otu_id)

                if otu is None:
                    msg = f"Cannot apply event to non-existent OTU {event.query.otu_id}"
                    raise ValueError(msg)

                # Apply event to validate it produces valid state
                if isinstance(event, ApplicableEvent):
                    try:
                        otu = event.apply(otu)
                        otu = OTU.model_validate(otu)
                    except Exception as e:
                        msg = f"Event validation failed: {e}"
                        raise ValueError(msg) from e

        # Validation passed - write to disk
        written_event = self._event_store.write_event(event)

        # Update index
        if hasattr(event.query, "otu_id"):
            self._index.add_event_id(
                written_event.id,
                event.query.otu_id,
                written_event.timestamp,
            )

        return written_event


@contextmanager
def locked_repo(path: Path) -> Generator[Repo]:
    """Yield a locked Repo."""
    repo = Repo(path)

    with repo.lock():
        yield repo
