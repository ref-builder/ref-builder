import datetime

import arrow
from structlog.testing import capture_logs

from ref_builder.models.accession import Accession
from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.otu.promote import (
    promote_otu_accessions_from_records,
    replace_otu_sequence_from_record,
    upgrade_outdated_sequences_in_otu,
)
from ref_builder.repo import Repo
from ref_builder.services.cls import Services
from tests.fixtures.factories import IsolateFactory


def test_replace_sequence(empty_repo: Repo):
    """Test OTU sequence replacement."""
    services = Services(empty_repo, NCBIClient(False))

    with empty_repo.lock():
        otu_before = services.otu.create(["MF062125", "MF062126", "MF062127"])

    assert otu_before
    assert otu_before.accessions == {"MF062125", "MF062126", "MF062127"}

    isolate = otu_before.isolates[0]

    assert isolate

    sequence_ids_before = [sequence.id for sequence in isolate.sequences]

    refseq_records = NCBIClient(ignore_cache=False).fetch_genbank_records(
        ["NC_055390", "NC_055391", "NC_055392"]
    )

    record_by_sequence_id = {
        sequence_ids_before[i]: refseq_records[i] for i in (0, 1, 2)
    }

    with empty_repo.lock():
        for sequence_id, record in record_by_sequence_id.items():
            assert replace_otu_sequence_from_record(
                empty_repo,
                otu_before,
                sequence_id,
                replacement_record=record,
                exclude_accession=True,
            )

    otu_after = empty_repo.get_otu(otu_before.id)

    assert otu_after
    assert otu_after.accessions == {
        "NC_055390",
        "NC_055391",
        "NC_055392",
    }


def test_multi_linked_promotion(empty_repo: Repo):
    """Test the promotion of a sequence that is linked to more than one isolate."""
    services = Services(empty_repo, NCBIClient(False))

    with empty_repo.lock():
        otu_before = services.otu.create(["MF062125", "MF062126", "MF062127"])

    assert otu_before
    assert otu_before.accessions == {"MF062125", "MF062126", "MF062127"}

    segment_before = otu_before.plan.get_segment_by_name_key("L")
    sequence_before = otu_before.get_sequence_by_accession("MF062125")

    assert segment_before
    assert sequence_before

    mock_isolate = IsolateBuilder.model_validate(
        IsolateFactory.build_on_plan(otu_before.plan).model_dump()
    )

    for i in range(len(mock_isolate.sequences)):
        mock_sequence = mock_isolate.sequences[i]

        if mock_sequence.segment == segment_before.id:
            mock_isolate.sequences[i] = sequence_before.model_copy()

    with empty_repo.lock(), empty_repo.use_transaction():
        isolate_before = empty_repo.create_isolate(
            otu_before.id, name=mock_isolate.name
        )

        assert isolate_before

        empty_repo.link_sequence(
            otu_before.id,
            isolate_before.id,
            sequence_before.id,
        )

        accession_counter = 1

        for mock_sequence in mock_isolate.sequences:
            # Skip segment L.
            if mock_sequence.segment == segment_before.id:
                continue

            assert mock_sequence.segment

            sequence_before = empty_repo.create_sequence(
                otu_before.id,
                accession=f"FA00000{accession_counter}.1",
                definition=mock_sequence.definition,
                segment=mock_sequence.segment,
                sequence=mock_sequence.sequence,
            )

            assert sequence_before

            empty_repo.link_sequence(
                otu_before.id,
                isolate_before.id,
                sequence_before.id,
            )

            accession_counter += 1

    otu_after = empty_repo.get_otu(otu_before.id)

    assert otu_after
    assert otu_after.accessions == {
        "MF062125",
        "MF062126",
        "MF062127",
        "FA000001",
        "FA000002",
    }
    assert isolate_before.id in otu_after.isolate_ids

    refseq_records = NCBIClient(ignore_cache=False).fetch_genbank_records(
        ["NC_055390", "NC_055391", "NC_055392"],
    )

    with empty_repo.lock():
        assert promote_otu_accessions_from_records(
            empty_repo, otu_after, refseq_records
        ) == {"NC_055390", "NC_055391", "NC_055392"}

    otu_after_promote = empty_repo.get_otu(otu_before.id)

    assert otu_after_promote
    assert otu_after_promote.accessions == {
        "NC_055390",
        "NC_055391",
        "NC_055392",
        "FA000001",
        "FA000002",
    }

    first_isolate = otu_after_promote.isolates[0]

    assert first_isolate.accessions == {"NC_055390", "NC_055391", "NC_055392"}

    assert otu_after_promote.get_isolate(isolate_before.id).accessions == (
        {"NC_055390", "FA000001", "FA000002"}
    )


class TestUpgradeSequencesInOTU:
    """Test OTU-wide outdated sequence version upgrade."""

    def test_ok(self, scratch_repo: Repo):
        """Test a simple fetch and replace upgrade."""
        services = Services(scratch_repo, NCBIClient(False))

        with scratch_repo.lock():
            otu_before = services.otu.create(["NC_004452.1"])

        assert otu_before
        assert "NC_004452" in otu_before.accessions

        sequence_before = otu_before.get_sequence_by_accession("NC_004452")

        assert sequence_before
        assert sequence_before.accession.version == 1

        containing_isolate_id = list(
            otu_before.get_isolate_ids_containing_sequence_id(sequence_before.id)
        )[0]

        isolate_before = otu_before.get_isolate(containing_isolate_id)

        assert isolate_before

        with scratch_repo.lock():
            upgraded_sequence_ids = upgrade_outdated_sequences_in_otu(
                scratch_repo, otu_before
            )

        updated_sequence_id = list(upgraded_sequence_ids)[0]

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert otu_after
        assert isolate_before.id in otu_after.isolate_ids

        isolate_after = otu_after.get_isolate(containing_isolate_id)

        assert isolate_after
        assert (
            isolate_after.sequence_ids
            == {updated_sequence_id}
            != isolate_before.sequence_ids
        )

        sequence_after = otu_after.get_sequence_by_id(updated_sequence_id)

        assert sequence_after
        assert sequence_after.accession.key == sequence_before.accession.key
        assert sequence_after.accession.version > sequence_before.accession.version

    def test_with_future_date_limit(self, scratch_repo: Repo):
        """Test that setting modification_date_start to a future date does returns no new sequences."""
        services = Services(scratch_repo, NCBIClient(False))

        with scratch_repo.lock():
            otu_before = services.otu.create(["NC_004452.1"])

        assert otu_before
        assert "NC_004452" in otu_before.accessions

        sequence = otu_before.get_sequence_by_accession("NC_004452")

        assert sequence
        assert sequence.accession == (Accession(key="NC_004452", version=1))

        with scratch_repo.lock(), capture_logs() as captured_logs:
            upgrade_outdated_sequences_in_otu(
                scratch_repo,
                otu_before,
                modification_date_start=arrow.utcnow().naive
                + datetime.timedelta(days=1),
            )

        assert any(
            log.get("event") == "All sequences are up to date." for log in captured_logs
        )
