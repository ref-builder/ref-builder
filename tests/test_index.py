"""Tests for :class:`.EventIndex`."""

import uuid
from pathlib import Path

import arrow
import pytest

from ref_builder.index import EventIndexItem, Index
from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.models.molecule import Molecule
from ref_builder.models.otu import OTUMinimal
from ref_builder.models.plan import Plan
from ref_builder.ncbi.client import NCBIClientProtocol
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.repo import Repo
from ref_builder.services.cls import Services

SNAPSHOT_AT_EVENT = (
    31,
    22,
    21,
    24,
    11,
    32,
    32,
    30,
)
"""Hardcoded at_event values for the snapshotter_otus fixture."""


@pytest.fixture
def index(indexable_otus: list[OTUBuilder], tmp_path: Path) -> Index:
    """An index with eight OTUs already cached."""
    _index = Index(tmp_path / "index.db")

    for otu, at_event in zip(indexable_otus, SNAPSHOT_AT_EVENT, strict=True):
        _index.upsert_otu(otu, at_event)

    return _index


class TestDeleteOTU:
    """Test the `delete_otu` method of the Snapshotter class."""

    def test_ok(self, index: Index, indexable_otus: list[OTUBuilder]):
        """Test that an OTU can be deleted from the index."""
        otu = indexable_otus[2]

        index.delete_otu(otu.id)

        assert index.get_event_ids_by_otu_id(otu.id) is None
        assert index.get_id_by_taxid(otu.taxid) is None

    def test_not_found(self, index: Index, indexable_otus: list[OTUBuilder]):
        """Test that nothing happens when the OTU ID is not found."""
        index.delete_otu(uuid.uuid4())
        assert index.otu_ids == {otu.id for otu in indexable_otus}


def test_iter_otus(index: Index, indexable_otus: list[OTUBuilder]):
    """Test that the index iterates over all OTUs ordered by name."""
    assert list(index.iter_minimal_otus()) == sorted(
        [
            OTUMinimal(
                acronym=otu.acronym,
                id=otu.id,
                name=otu.name,
                taxid=otu.taxid,
            )
            for otu in indexable_otus
        ],
        key=lambda otu: otu.name,
    )


class TestLoadSnapshot:
    """Test the `load_snapshot` method of the Snapshotter class."""

    def test_ok(self, index: Index, indexable_otus: list[OTUBuilder]):
        """Test that we can load a snapshot from the index."""
        for at_event, otu in zip(SNAPSHOT_AT_EVENT, indexable_otus, strict=False):
            snapshot = index.load_snapshot(otu.id)

            assert snapshot.at_event == at_event
            assert snapshot.otu == otu

    def test_not_found(self, index: Index):
        """Test that `None` is returned when the OTU ID is not found."""
        assert index.load_snapshot(uuid.uuid4()) is None


def test_otu_ids(index: Index, indexable_otus: list[OTUBuilder]):
    """Test that stored OTU IDs are correct."""
    assert index.otu_ids == {otu.id for otu in indexable_otus}


class TestEvents:
    """Test the event index functionality of the repository Index."""

    def test_ok(self, index: Index, indexable_otus: list[OTUBuilder]):
        """Test that we can set and get events IDs for an OTU."""
        otu = indexable_otus[1]

        index.add_event_id(100, otu.id, arrow.utcnow().naive)
        index.add_event_id(101, otu.id, arrow.utcnow().naive)
        index.add_event_id(104, otu.id, arrow.utcnow().naive)

        assert index.get_event_ids_by_otu_id(otu.id) == EventIndexItem(
            event_ids=[100, 101, 104],
            otu_id=otu.id,
        )

    def test_otu_id_not_found(self, index: Index):
        """Test that we get ``None`` when an OTU ID is not found."""
        assert index.get_event_ids_by_otu_id(uuid.uuid4()) is None

    def test_get_first_timestamp_ok(
        self, index: Index, indexable_otus: list[OTUBuilder]
    ):
        """Test ``.get_first_timestamp_by_otu_id()`` retrieves the first timestamp."""
        otu = indexable_otus[1]

        first_timestamp = arrow.utcnow().naive

        index.add_event_id(100, otu.id, first_timestamp)

        assert index.get_first_timestamp_by_otu_id(otu.id) == first_timestamp

        second_timestamp = arrow.utcnow().naive

        index.add_event_id(101, otu.id, second_timestamp)

        assert index.get_first_timestamp_by_otu_id(otu.id) == first_timestamp

    def test_get_latest_timestamp_ok(
        self, index: Index, indexable_otus: list[OTUBuilder]
    ):
        """Test ``.get_latest_timestamp_by_otu_id()`` retrieves the latest timestamp."""
        otu = indexable_otus[1]

        index.add_event_id(100, otu.id, arrow.utcnow().naive)

        first_timestamp = arrow.utcnow().naive

        index.add_event_id(101, otu.id, first_timestamp)

        assert index.get_latest_timestamp_by_otu_id(otu.id) == first_timestamp

        second_timestamp = arrow.utcnow().naive

        index.add_event_id(104, otu.id, second_timestamp)

        assert second_timestamp > first_timestamp

        assert index.get_latest_timestamp_by_otu_id(otu.id) == second_timestamp


class TestGetIDByTaxid:
    """Test `Index.get_id_by_taxid`."""

    def test_ok(self, index: Index, indexable_otus: list[OTUBuilder]):
        """Test that the correct OTU ID is retrieved by taxid."""
        for otu in indexable_otus:
            assert index.get_id_by_taxid(otu.taxid) == otu.id

    def test_any_lineage_taxid(self, index: Index, indexable_otus: list[OTUBuilder]):
        """Test that OTU can be found by ANY taxid in its lineage."""
        for otu in indexable_otus:
            for taxon in otu.lineage.taxa:
                assert index.get_id_by_taxid(taxon.id) == otu.id

    def test_not_found(self, index: Index):
        """Test that `None` is returned when the taxid is not found."""
        assert index.get_id_by_taxid(999999999999999) is None

    def test_subspecies_taxid_lookup(
        self, empty_repo: Repo, mock_ncbi_client: NCBIClientProtocol
    ):
        """Test that OTU can be found by sub-species taxid.

        Wasabi mottle virus has:
        - Target taxon (sub-species): taxid=1169032, rank="no rank"
        - Species: taxid=3432896, rank="species"

        Both taxids should find the same OTU.
        """
        from tests.fixtures.ncbi import OTUManifest

        services = Services(empty_repo, mock_ncbi_client)

        with empty_repo.lock():
            otu = services.otu.create(OTUManifest.wasabi_mottle_virus.refseq)

        assert otu is not None

        index = Index(empty_repo.path / "index.db")
        index.upsert_otu(otu, 1)

        subspecies_taxid = 1169032
        species_taxid = 3432896

        otu_id_subspecies = index.get_id_by_taxid(subspecies_taxid)
        otu_id_species = index.get_id_by_taxid(species_taxid)

        assert otu_id_subspecies is not None
        assert otu_id_species is not None
        assert otu_id_subspecies == otu_id_species == otu.id


class TestGetIDByIsolateID:
    def test_ok(self, index: Index, indexable_otus: list[OTUBuilder]):
        """Test the `get_id_by_isolate_id` method of the Index class."""
        for otu in indexable_otus:
            first_isolate = next(iter(otu.isolate_ids))

            assert index.get_id_by_isolate_id(first_isolate) == otu.id
