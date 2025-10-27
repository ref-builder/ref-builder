from structlog.testing import capture_logs

from ref_builder.models.isolate import IsolateName, IsolateNameType
from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.models.molecule import Molecule, MoleculeType, Strandedness, Topology
from ref_builder.models.plan import Plan, Segment
from ref_builder.ncbi.models import NCBIRank
from ref_builder.repo import Repo

TMV_LINEAGE = Lineage(
    taxa=[
        Taxon(
            id=12242,
            name="Tobacco mosaic virus",
            parent=3432891,
            rank=NCBIRank.NO_RANK,
            other_names=TaxonOtherNames(acronym=["TMV"], synonyms=[]),
        ),
        Taxon(
            id=3432891,
            name="Tobamovirus tabaci",
            parent=None,
            rank=NCBIRank.SPECIES,
            other_names=TaxonOtherNames(acronym=[], synonyms=[]),
        ),
    ]
)


def test_commit(empty_repo: Repo):
    """Test a successful transaction."""
    plan = Plan.new(
        [
            Segment.new(
                length=15,
                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                name=None,
            )
        ]
    )

    with empty_repo.lock(), empty_repo.use_transaction():
        otu = empty_repo.create_otu(
            lineage=TMV_LINEAGE,
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MoleculeType.RNA,
                topology=Topology.LINEAR,
            ),
            plan=plan,
        )

        assert otu

        isolate = empty_repo.create_isolate(
            otu.id,
            IsolateName(IsolateNameType.ISOLATE, "Test"),
            taxid=12242,
        )

        sequence = empty_repo.create_sequence(
            otu.id,
            accession="NC_001367.1",
            definition="Tobacco mosaic virus",
            segment=plan.segments[0].id,
            sequence="ACGTACGTACGTACG",
        )

        assert isolate
        assert sequence

        empty_repo.link_sequence(otu.id, isolate.id, sequence.id)

    assert empty_repo.last_id == 5
    assert len(list(empty_repo.iter_otus())) == 1


def test_fail(empty_repo: Repo):
    """Test auto-validation behaviour. If the transaction results in an invalid OTU,
    the repo should roll back all events.
    """
    plan = Plan.new(
        [
            Segment.new(
                length=15,
                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                name=None,
            )
        ]
    )

    assert empty_repo.last_id == 1

    with capture_logs() as cap_logs, empty_repo.lock(), empty_repo.use_transaction():
        empty_repo.create_otu(
            lineage=TMV_LINEAGE,
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MoleculeType.RNA,
                topology=Topology.LINEAR,
            ),
            plan=plan,
        )

    assert any(
        "OTU does not pass validation" in cap_log["event"] for cap_log in cap_logs
    )

    assert empty_repo.last_id == 1
    assert len(list(empty_repo.iter_otus())) == 0


def test_abort(empty_repo: Repo):
    """Test manual transaction abort. The repo should roll back all events."""
    plan = Plan.new(
        [
            Segment.new(
                length=15,
                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                name=None,
            )
        ]
    )

    with empty_repo.lock(), empty_repo.use_transaction() as transaction:
        empty_repo.create_otu(
            lineage=TMV_LINEAGE,
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MoleculeType.RNA,
                topology=Topology.LINEAR,
            ),
            plan=plan,
        )

        assert empty_repo.last_id == 2
        assert len(list(empty_repo.iter_otus())) == 1

        transaction.abort()

    assert empty_repo.last_id == 1
    assert len(list(empty_repo.iter_otus())) == 0
