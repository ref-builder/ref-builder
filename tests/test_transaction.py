from structlog.testing import capture_logs

from ref_builder.repo import Repo
from tests.fixtures.factories import OTUFactory


def test_commit(empty_repo: Repo, otu_factory: OTUFactory):
    """Test a successful transaction."""
    fake_otu = otu_factory.build()

    with empty_repo.lock(), empty_repo.use_transaction():
        otu = empty_repo.create_otu(
            fake_otu.acronym,
            molecule=fake_otu.molecule,
            name=fake_otu.name,
            plan=fake_otu.plan,
            taxid=fake_otu.taxid,
        )

        assert otu

        isolate = empty_repo.create_isolate(
            otu.id,
            name=fake_otu.isolates[0].name,
        )

        for fake_sequence in fake_otu.isolates[0].sequences:
            sequence = empty_repo.create_sequence(
                otu.id,
                accession=str(fake_sequence.accession),
                definition=fake_sequence.definition,
                segment=fake_sequence.segment,
                sequence=fake_sequence.sequence,
            )

            assert isolate
            assert sequence

            empty_repo.link_sequence(otu.id, isolate.id, sequence.id)

    assert empty_repo.last_id == 5
    assert len(list(empty_repo.iter_otus())) == 1


def test_fail(empty_repo: Repo, otu_factory: OTUFactory):
    """Test auto-validation behaviour. If the transaction results in an invalid OTU,
    the repo should roll back all events.
    """
    fake_otu = otu_factory.build()

    assert empty_repo.last_id == 1

    with capture_logs() as cap_logs, empty_repo.lock(), empty_repo.use_transaction():
        empty_repo.create_otu(
            fake_otu.acronym,
            molecule=fake_otu.molecule,
            name=fake_otu.name,
            plan=fake_otu.plan,
            taxid=fake_otu.taxid,
        )

    assert any(
        "OTU does not pass validation" in cap_log["event"] for cap_log in cap_logs
    )

    assert empty_repo.last_id == 1
    assert len(list(empty_repo.iter_otus())) == 0


def test_abort(empty_repo: Repo, otu_factory: OTUFactory):
    """Test manual transaction abort. The repo should roll back all events."""
    otu = otu_factory.build()

    with empty_repo.lock(), empty_repo.use_transaction() as transaction:
        empty_repo.create_otu(
            otu.acronym,
            molecule=otu.molecule,
            name=otu.name,
            plan=otu.plan,
            taxid=otu.taxid,
        )

        assert empty_repo.last_id == 2
        assert len(list(empty_repo.iter_otus())) == 1

        transaction.abort()

    assert empty_repo.last_id == 1
    assert len(list(empty_repo.iter_otus())) == 0
