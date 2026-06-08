"""Tests for ``ref-builder validate accessions``."""

from pathlib import Path
from uuid import uuid4

import arrow
import pytest
from click.testing import CliRunner

from ref_builder.cli.validate_cmd import validate as validate_command_group
from ref_builder.events.base import IsolateQuery
from ref_builder.events.isolate import CreateIsolate, CreateIsolateData
from ref_builder.models.accession import Accession
from ref_builder.models.isolate import IsolateName, IsolateNameType
from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.models.molecule import Molecule, MoleculeType, Strandedness, Topology
from ref_builder.models.plan import Plan, Segment
from ref_builder.models.sequence import Sequence
from ref_builder.ncbi.models import NCBIRank
from ref_builder.repo import Repo

runner = CliRunner()

SEGMENT_LENGTH = 15

TMV_LINEAGE = Lineage(
    taxa=[
        Taxon(
            id=3432891,
            name="Tobamovirus tabaci",
            parent=None,
            rank=NCBIRank.SPECIES,
            other_names=TaxonOtherNames(acronym=[], synonyms=[]),
        )
    ]
)

CMV_LINEAGE = Lineage(
    taxa=[
        Taxon(
            id=12305,
            name="Cucumovirus mosaiccucumeris",
            parent=None,
            rank=NCBIRank.SPECIES,
            other_names=TaxonOtherNames(acronym=[], synonyms=[]),
        )
    ]
)

MOLECULE = Molecule(
    strandedness=Strandedness.SINGLE,
    type=MoleculeType.RNA,
    topology=Topology.LINEAR,
)


def _plan(repo: Repo) -> Plan:
    return Plan.new(
        [
            Segment.new(
                length=SEGMENT_LENGTH,
                length_tolerance=repo.settings.default_segment_length_tolerance,
                name=None,
            )
        ]
    )


def _seed_isolate(plan: Plan, accession_key: str, taxid: int) -> CreateIsolateData:
    return CreateIsolateData(
        id=uuid4(),
        name=IsolateName(IsolateNameType.ISOLATE, "Seed"),
        taxid=taxid,
        sequences=[
            Sequence(
                accession=Accession(key=accession_key, version=1),
                definition="seed",
                segment=plan.segments[0].id,
                sequence="ACGTACGTACGTACG",
            )
        ],
    )


@pytest.fixture
def clean_repo(empty_repo: Repo) -> Repo:
    """A two-OTU repo with no duplicate accessions."""
    plan = _plan(empty_repo)

    with empty_repo.lock():
        empty_repo.create_otu(
            isolate=_seed_isolate(plan, "TM000001", taxid=3432891),
            lineage=TMV_LINEAGE,
            molecule=MOLECULE,
            plan=plan,
            promoted_accessions=set(),
        )
        empty_repo.create_otu(
            isolate=_seed_isolate(plan, "AB000001", taxid=12305),
            lineage=CMV_LINEAGE,
            molecule=MOLECULE,
            plan=plan,
            promoted_accessions=set(),
        )

    return empty_repo


def _raw_write_create_isolate(
    repo: Repo,
    otu_id,
    accession_key: str,
    plan: Plan,
    taxid: int,
    isolate_label: str,
) -> None:
    """Write a CreateIsolate event directly, bypassing all validation guards."""
    event = CreateIsolate(
        id=repo.last_id + 1,
        data=CreateIsolateData(
            id=uuid4(),
            name=IsolateName(IsolateNameType.ISOLATE, isolate_label),
            taxid=taxid,
            sequences=[
                Sequence(
                    accession=Accession(key=accession_key, version=1),
                    definition=isolate_label,
                    segment=plan.segments[0].id,
                    sequence="ACGTACGTACGTACG",
                )
            ],
        ),
        query=IsolateQuery(isolate_id=uuid4(), otu_id=otu_id),
        timestamp=arrow.utcnow().naive,
    )
    written = repo._event_store.write_event(event)
    repo._index.add_event_id(written.id, otu_id, written.timestamp)


def _corrupt_repo_with_dupes(repo: Repo) -> None:
    """Build a repo that contains both cross-OTU and within-OTU duplicate accessions.

    Writes the conflicting events directly to the event store so the corruption
    mimics what an already-broken on-disk reference (e.g. plant-viruses) looks like.
    """
    plan = _plan(repo)

    with repo.lock():
        otu_tmv = repo.create_otu(
            isolate=_seed_isolate(plan, "TM000001", taxid=3432891),
            lineage=TMV_LINEAGE,
            molecule=MOLECULE,
            plan=plan,
            promoted_accessions=set(),
        )
        otu_cmv = repo.create_otu(
            isolate=_seed_isolate(plan, "AB000001", taxid=12305),
            lineage=CMV_LINEAGE,
            molecule=MOLECULE,
            plan=plan,
            promoted_accessions=set(),
        )

        _raw_write_create_isolate(
            repo, otu_cmv.id, "TM000001", plan, 12305, "DupeAcross"
        )
        # Distinct accession reused within TMV so the within-OTU report fires.
        _raw_write_create_isolate(
            repo, otu_tmv.id, "TN000001", plan, 3432891, "WithinA"
        )
        _raw_write_create_isolate(
            repo, otu_tmv.id, "TN000001", plan, 3432891, "WithinB"
        )


class TestValidateAccessionsCommand:
    def test_clean_repo_exits_zero(self, clean_repo: Repo):
        """A repo with no duplicates exits 0 and reports cleanly."""
        result = runner.invoke(
            validate_command_group, ["accessions", "--path", str(clean_repo.path)]
        )

        assert result.exit_code == 0
        assert "No duplicate accessions found" in result.output

    def test_corrupt_repo_reports_and_exits_nonzero(self, tmp_path: Path):
        """A repo containing both kinds of duplicate exits non-zero and reports both."""
        repo = Repo.new("Generic Viruses", tmp_path / "corrupt_repo", "virus")
        _corrupt_repo_with_dupes(repo)

        result = runner.invoke(
            validate_command_group, ["accessions", "--path", str(repo.path)]
        )

        assert result.exit_code == 1
        assert "Cross-OTU duplicate accessions" in result.output
        assert "Within-OTU duplicate accessions" in result.output
        assert "TM000001" in result.output
