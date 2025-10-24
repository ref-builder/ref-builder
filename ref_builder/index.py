"""An index for improving repository performance."""

import binascii
import datetime
import sqlite3
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path
from uuid import UUID

from ref_builder.events.base import EventMetadata
from ref_builder.models.otu import OTUMinimal
from ref_builder.otu.builders.otu import OTUBuilder


@dataclass
class EventIndexItem:
    """The event IDs associated with an OTU."""

    event_ids: list[int]
    """A list of event IDs"""

    otu_id: UUID
    """The OTU ID this cache index is associated with."""


@dataclass
class Snapshot:
    """A snapshot of an OTU at a specific event."""

    at_event: int
    """The event at which the snapshot was taken."""

    otu: OTUBuilder
    """The OTU that was snapshotted."""


def adapt_uuid(original_uuid: UUID) -> str:
    """Automatically convert UUID to SQLite-compatible string."""
    return str(original_uuid)


def adapt_datetime(original_datetime: datetime.datetime) -> str:
    """Automatically convert datetime to SQLite-compatible string."""
    return original_datetime.isoformat()


def convert_uuid(string_uuid: bytes) -> UUID:
    """Automatically convert stringified UUID back to UUID type."""
    return UUID(string_uuid.decode())


def convert_datetime(string_datetime: bytes) -> datetime.datetime:
    """Automatically convert stringified datetime back to datetime type."""
    return datetime.datetime.fromisoformat(string_datetime.decode())


sqlite3.register_adapter(UUID, adapt_uuid)
sqlite3.register_adapter(datetime.datetime, adapt_datetime)

sqlite3.register_converter("uuid", convert_uuid)
sqlite3.register_converter("datetime", convert_datetime)


class Index:
    """An index for rapid access to repository data."""

    def __init__(self, path: Path) -> None:
        """Initialize the index."""
        self.path = path
        """The path to the index file."""

        self.path.parent.mkdir(exist_ok=True, parents=True)

        self.con = sqlite3.connect(
            path,
            isolation_level=None,
            detect_types=sqlite3.PARSE_DECLTYPES | sqlite3.PARSE_COLNAMES,
        )

        self.con.execute("PRAGMA journal_mode = WAL")
        self.con.execute("PRAGMA synchronous = NORMAL")

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS events (
                event_id INTEGER PRIMARY KEY,
                otu_id TEXT,
                timestamp TEXT
            );
            """,
        )

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS isolates (
                id TEXT PRIMARY KEY,
                name TEXT,
                otu_id TEXT
            );
            """,
        )

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS otus (
                id TEXT PRIMARY KEY,
                acronym TEXT,
                at_event INTEGER,
                name TEXT,
                otu JSONB,
                taxid INTEGER
            );
            """,
        )

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS sequences (
                id TEXT PRIMARY KEY,
                crc TEXT,
                otu_id TEXT,
                sequence TEXT
            );
            """,
        )

        self.con.execute(
            """
            CREATE TABLE IF NOT EXISTS otu_updates (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                otu_id TEXT,
                timestamp_complete TEXT
            );
            """,
        )

        for table, column in [
            ("events", "otu_id"),
            ("isolates", "id"),
            ("otus", "id"),
            ("otus", "name"),
            ("otus", "taxid"),
            ("sequences", "otu_id"),
            ("sequences", "crc"),
            ("otu_updates", "otu_id"),
        ]:
            self.con.execute(
                f"""
                CREATE INDEX IF NOT EXISTS {table}_{column} ON {table} ({column});
                """,
            )

        self.con.commit()

    @property
    def otu_ids(self) -> set[UUID]:
        """A list of all OTUs tracked in the index."""
        return {
            row[0] for row in self.con.execute('SELECT id AS "id [uuid]" FROM otus')
        }

    def add_event_id(
        self,
        event_id: int,
        otu_id: UUID,
        timestamp: datetime.datetime,
    ) -> None:
        """Add a new event ID and associated it with a given OTU ID.

        This allows fast lookup of event IDs for a given OTU ID.

        :param event_id: the event ID to add
        :param otu_id: the OTU ID to add the event ID to
        :param timestamp: the timestamp of the event
        """
        self.con.execute(
            """
            INSERT INTO events (event_id, otu_id, timestamp)
            VALUES (?, ?, ?)
            """,
            (
                event_id,
                otu_id,
                timestamp,
            ),
        )

    def get_event_ids_by_otu_id(self, otu_id: UUID) -> EventIndexItem | None:
        """Get a list of event IDs for an OTU."""
        res = self.con.execute(
            """
            SELECT event_id
            FROM events
            WHERE otu_id = ?
            ORDER BY event_id
            """,
            (otu_id,),
        )

        event_ids = [row[0] for row in res]

        if event_ids:
            return EventIndexItem(
                event_ids,
                otu_id,
            )

        return None

    def get_first_timestamp_by_otu_id(self, otu_id: UUID) -> datetime.datetime | None:
        """Get the timestamp of the earliest event associated with this ID."""
        cursor = self.con.execute(
            """
            SELECT timestamp
            AS "timestamp [datetime]"
            FROM events
            WHERE otu_id = ?
            ORDER BY event_id
            """,
            (otu_id,),
        )

        if result := cursor.fetchone():
            return result[0]

        return None

    def get_latest_timestamp_by_otu_id(self, otu_id: UUID) -> datetime.datetime | None:
        """Get the timestamp of the latest event associated with this ID."""
        cursor = self.con.execute(
            """
            SELECT timestamp
            AS "timestamp [datetime]"
            FROM events
            WHERE otu_id = ?
            ORDER BY event_id DESC
            """,
            (otu_id,),
        )

        if result := cursor.fetchone():
            return result[0]

        return None

    def get_id_by_acronym(self, acronym: str) -> UUID | None:
        """Get an OTU ID by its acronym."""
        if not acronym:
            raise ValueError("Acronym must contain at least one character.")

        cursor = self.con.execute(
            'SELECT id AS "id [uuid]" FROM otus WHERE acronym = ?',
            (acronym,),
        )

        if result := cursor.fetchone():
            return result[0]

        return None

    def get_id_by_name(self, name: str) -> UUID | None:
        """Get an OTU ID by its name."""
        cursor = self.con.execute(
            'SELECT id AS "id [uuid]" FROM otus WHERE name = ?',
            (name,),
        )

        if result := cursor.fetchone():
            return result[0]

        return None

    def get_id_by_taxid(self, taxid: int) -> UUID | None:
        """Get an OTU ID by its taxonomy ID."""
        cursor = self.con.execute(
            'SELECT id AS "id [uuid]" FROM otus WHERE taxid = ?',
            (taxid,),
        )

        if result := cursor.fetchone():
            return result[0]

        return None

    def get_id_by_isolate_id(self, isolate_id: UUID) -> UUID | None:
        """Get an OTU ID from an isolate ID that belongs to it."""
        cursor = self.con.execute(
            'SELECT otu_id AS "otu_id [uuid]" FROM isolates WHERE id = ?',
            (isolate_id,),
        )

        if result := cursor.fetchone():
            return result[0]

        return None

    def delete_otu(self, otu_id: UUID) -> None:
        """Remove an OTU from the index.

        This happens rarely, so we just rewrite the whole index.

        :param otu_id: the ID of the OTU to remove.
        """
        self.con.execute(
            "DELETE FROM isolates WHERE otu_id = ?",
            (otu_id,),
        )

        self.con.execute(
            "DELETE FROM otus WHERE id = ?",
            (otu_id,),
        )

        self.con.execute(
            "DELETE FROM sequences WHERE otu_id = ?",
            (otu_id,),
        )

        self.con.commit()

    def iter_event_metadata(self) -> Iterator[EventMetadata]:
        """Iterate over event metadata."""
        rows = self.con.execute(
            "SELECT event_id, otu_id, timestamp FROM events ORDER BY event_id",
        ).fetchall()

        for row in rows:
            yield EventMetadata(
                id=row[0],
                otu_id=UUID(row[1]) if row[1] else None,
                timestamp=row[2],
            )

    def iter_minimal_otus(self) -> Iterator[OTUMinimal]:
        """Iterate over minimal representations of all OTUs in the index."""
        rows = self.con.execute(
            """
            SELECT acronym, id AS "id [uuid]", name, taxid
            FROM otus ORDER BY name
            """,
        ).fetchall()

        for row in rows:
            yield OTUMinimal(
                acronym=row[0],
                id=row[1],
                name=row[2],
                taxid=row[3],
            )

    def get_last_otu_update_timestamp(self, otu_id: UUID) -> datetime.datetime | None:
        """Get the timestamp of the last event in an update for an OTU."""
        cursor = self.con.execute(
            """
            SELECT timestamp_complete FROM otu_updates
            WHERE otu_id = ? ORDER BY id DESC
            """,
            (otu_id,),
        )

        if result := cursor.fetchone():
            return datetime.datetime.fromisoformat(result[0])

        return None

    def add_otu_update_history_entry(
        self, otu_id: UUID, timestamp: datetime.datetime
    ) -> int | None:
        """Write an entry into the OTU update history."""
        self.con.execute(
            """
            INSERT INTO
            otu_updates(otu_id, timestamp_complete)
            VALUES(?, ?)
            """,
            (
                otu_id,
                timestamp,
            ),
        )

        cursor = self.con.execute(
            "SELECT id FROM otu_updates WHERE otu_id = ? ORDER BY id DESC",
            (otu_id,),
        )

        if result := cursor.fetchone():
            return result[0]

        return None

    def load_snapshot(self, otu_id: UUID) -> Snapshot | None:
        """Load an OTU snapshot."""
        cursor = self.con.execute(
            "SELECT at_event, otu FROM otus WHERE id = ?",
            (otu_id,),
        )

        result = cursor.fetchone()

        if not result:
            return None

        at_event, otu_json = result

        otu = OTUBuilder.model_validate_json(otu_json)

        sequence_ids = [
            sequence.id for isolate in otu.isolates for sequence in isolate.sequences
        ]

        # Fetch all sequences in a single query
        placeholders = ",".join("?" for _ in sequence_ids)

        sequence_map = dict(
            self.con.execute(
                f"SELECT id, sequence FROM sequences WHERE id IN ({placeholders})",
                sequence_ids,
            ).fetchall(),
        )

        # Update the data structure and check for missing sequences
        missing_sequences = []

        for isolate in otu.isolates:
            for sequence in isolate.sequences:
                sequence_id = str(sequence.id)

                if sequence_id in sequence_map:
                    sequence.sequence = sequence_map[sequence_id]
                else:
                    missing_sequences.append(str(sequence_id))

        # Raise an error if any sequences are missing
        if missing_sequences:
            raise ValueError(f"Sequences not found: {', '.join(missing_sequences)}")

        return Snapshot(
            at_event=at_event,
            otu=OTUBuilder.model_validate(otu),
        )

    def upsert_otu(self, otu: OTUBuilder, at_event: int) -> None:
        """Update the index based on a list of OTUs.

        OTUs that aren't already in the index will be added, and those that are will be
        updated.

        Only sequences that have changed will be updated.
        """
        sequence_ids = {seq.id for isolate in otu.isolates for seq in isolate.sequences}

        placeholders = ",".join("?" for _ in sequence_ids)

        # Delete any sequences that are no longer in the OTU.
        self.con.execute(
            f"""
            DELETE FROM sequences WHERE otu_id = ? AND NOT id IN ({placeholders});
            """,
            (otu.id, *sequence_ids),
        )

        batch = []

        # Insert or update the isolates.
        for isolate in otu.isolates:
            self.con.execute(
                "INSERT OR REPLACE INTO isolates (id, name, otu_id) VALUES (?, ?, ?)",
                (
                    isolate.id,
                    str(isolate.name),
                    otu.id,
                ),
            )

            for sequence in isolate.sequences:
                crc = _calculate_crc32(sequence.sequence)
                batch.append(
                    (
                        sequence.id,
                        crc,
                        otu.id,
                        sequence.sequence,
                    ),
                )

        # Perform batch insert
        self.con.executemany(
            """
            INSERT OR REPLACE INTO sequences (id, crc, otu_id, sequence)
            VALUES (?, ?, ?, ?)
            """,
            batch,
        )

        # Update only if CRC has changed
        self.con.executemany(
            """
            UPDATE sequences
            SET crc = ?, otu_id = ?, sequence = ?
            WHERE id = ? AND crc != ?
            """,
            [(crc, otu_id, seq, seq_id, crc) for seq_id, crc, otu_id, seq in batch],
        )

        self.con.execute(
            """
            INSERT OR REPLACE INTO otus (id, acronym, at_event, name, otu, taxid)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                otu.id,
                otu.acronym,
                at_event,
                otu.name,
                otu.model_dump_json(),
                otu.taxid,
            ),
        )

        self.con.commit()

    def prune(self, event_id: int) -> None:
        """Rollback the index to a previous event.

        * Deletes all events after the given event ID.
        * Deletes all OTU snapshots that have ``at_event`` later than ``event_id``.

        Removed snapshots will be automatically updated next time an OTU is fetched from
        the repository.

        :param event_id: the event ID to rollback to
        """
        self.con.execute(
            """
            DELETE FROM events WHERE event_id > ?
            """,
            (event_id,),
        )

        self.con.execute(
            """
            DELETE FROM otus WHERE at_event > ?
            """,
            (event_id,),
        )

        self.con.commit()


def _calculate_crc32(sequence: str) -> str:
    """Calculate the CRC32 checksum of a sequence.

    :param sequence: the sequence to calculate the checksum of.
    :return: the CRC32 checksum as a hexadecimal string.

    """
    crc = binascii.crc32(sequence.encode("utf-8"))

    # Convert CRC as a hexadecimal string.
    return hex(crc & 0xFFFFFFFF)[2:].zfill(8)
