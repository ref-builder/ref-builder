from pathlib import Path
from uuid import UUID, uuid4

import orjson
import pytest
from structlog.testing import capture_logs

from ref_builder.events.isolate import CreateIsolateData
from ref_builder.models.accession import Accession
from ref_builder.models.isolate import Isolate, IsolateName, IsolateNameType
from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.models.molecule import Molecule, MoleculeType, Strandedness, Topology
from ref_builder.models.otu import OTU
from ref_builder.models.plan import Plan, Segment, SegmentRule
from ref_builder.models.sequence import Sequence
from ref_builder.ncbi.models import NCBIRank
from ref_builder.repo import GITIGNORE_CONTENTS, Repo

SEGMENT_LENGTH = 15

TMV_LINEAGE = Lineage(
    taxa=[
        Taxon(
            id=3432891,
            name="Tobamovirus tabaci",
            parent=None,
            rank=NCBIRank.SPECIES,
            other_names=TaxonOtherNames(acronym=[], synonyms=[]),
        ),
        Taxon(
            id=12242,
            name="Tobacco mosaic virus",
            parent=3432891,
            rank=NCBIRank.NO_RANK,
            other_names=TaxonOtherNames(acronym=["TMV"], synonyms=[]),
        ),
        Taxon(
            id=12227,
            name="Tobacco mosaic virus",
            parent=3432891,
            rank=NCBIRank.NO_RANK,
            other_names=TaxonOtherNames(acronym=[], synonyms=[]),
        ),
    ]
)


@pytest.fixture
def initialized_repo(tmp_path: Path):
    """Return a pre-initialized mock Repo."""
    repo = Repo.new(
        "Generic Viruses",
        tmp_path / "initialized_repo",
        "virus",
    )

    with repo.lock():
        plan = Plan.new(
            [
                Segment.new(
                    length=SEGMENT_LENGTH,
                    length_tolerance=repo.settings.default_segment_length_tolerance,
                    name=None,
                )
            ]
        )

        isolate_data = CreateIsolateData(
            id=uuid4(),
            name=IsolateName(IsolateNameType.ISOLATE, "A"),
            taxid=12227,
            sequences=[
                Sequence(
                    accession=Accession(key="TM000001", version=1),
                    definition="TMV",
                    segment=plan.segments[0].id,
                    sequence="ACGTACGTACGTACG",
                )
            ],
        )

        otu = repo.create_otu(
            isolate=isolate_data,
            lineage=TMV_LINEAGE,
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MoleculeType.RNA,
                topology=Topology.LINEAR,
            ),
            plan=plan,
            promoted_accessions=set(),
        )

        assert otu is not None

    return repo


def init_otu(repo: Repo) -> OTU:
    """Create an OTU with one isolate."""
    plan = Plan.new(
        segments=[
            Segment.new(
                length=SEGMENT_LENGTH,
                length_tolerance=repo.settings.default_segment_length_tolerance,
                name=None,
            )
        ]
    )

    isolate_data = CreateIsolateData(
        id=uuid4(),
        name=None,
        taxid=12227,
        sequences=[],
    )

    result = repo.create_otu(
        isolate=isolate_data,
        lineage=TMV_LINEAGE,
        molecule=Molecule(
            strandedness=Strandedness.SINGLE,
            type=MoleculeType.RNA,
            topology=Topology.LINEAR,
        ),
        plan=plan,
        promoted_accessions=set(),
    )
    assert result is not None
    return result


class TestNew:
    def test_ok(self, empty_repo: Repo, tmp_path: Path):
        """Test that creating a new ``Repo`` object returns the expected object and
        creates the expected directory structure.
        """
        assert empty_repo.path == tmp_path / "test_repo"
        assert empty_repo.last_id == 1

        assert empty_repo.meta.name == "Generic Viruses"
        assert empty_repo.meta.organism == "virus"

        assert empty_repo.settings.default_segment_length_tolerance == 0.03

        assert (empty_repo.path / ".gitignore").exists()

        with open(empty_repo.path / ".gitignore") as f:
            assert f.read() == "\n".join(GITIGNORE_CONTENTS) + "\n"

    def test_alternate_settings(self, tmp_path: Path):
        """Test retrieval of non-default settings."""
        repo = Repo.new(
            "Generic Viruses",
            tmp_path / "alt_setting_repo",
            "virus",
            default_segment_length_tolerance=0.05,
        )

        assert repo.settings.default_segment_length_tolerance == 0.05


class TestCreateOTU:
    def test_empty(self, empty_repo: Repo):
        """Test that creating an OTU.

        The method should create the correct even and return an OTU.
        """
        plan = Plan.new(
            segments=[
                Segment.new(
                    length=SEGMENT_LENGTH,
                    length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                    name=None,
                )
            ]
        )

        isolate_data = CreateIsolateData(
            id=uuid4(),
            name=None,
            taxid=12227,
            sequences=[],
        )

        with empty_repo.lock():
            otu = empty_repo.create_otu(
                isolate=isolate_data,
                lineage=TMV_LINEAGE,
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MoleculeType.RNA,
                    topology=Topology.LINEAR,
                ),
                plan=plan,
                promoted_accessions=set(),
            )

            assert otu
            assert otu == OTU.model_validate(
                {
                    "id": otu.id,
                    "excluded_accessions": set(),
                    "promoted_accessions": set(),
                    "lineage": TMV_LINEAGE,
                    "molecule": Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MoleculeType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    "plan": Plan(
                        id=plan.id,
                        segments=[
                            Segment(
                                id=plan.segments[0].id,
                                length=SEGMENT_LENGTH,
                                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                                name=None,
                                rule=SegmentRule.REQUIRED,
                            )
                        ],
                    ),
                    "isolates": [
                        Isolate(
                            id=isolate_data.id,
                            name=None,
                            taxid=12227,
                            sequences=[],
                        )
                    ],
                }
            )

            with open(empty_repo.path.joinpath("src", "00000002.json")) as f:
                event = orjson.loads(f.read())

            del event["timestamp"]

            assert event == {
                "data": {
                    "id": str(otu.id),
                    "isolate": {
                        "id": str(isolate_data.id),
                        "name": None,
                        "taxid": 12227,
                        "sequences": [],
                    },
                    "lineage": {
                        "taxa": [
                            {
                                "id": 3432891,
                                "name": "Tobamovirus tabaci",
                                "parent": None,
                                "rank": "species",
                                "other_names": {"acronym": [], "synonyms": []},
                            },
                            {
                                "id": 12242,
                                "name": "Tobacco mosaic virus",
                                "parent": 3432891,
                                "rank": "no rank",
                                "other_names": {"acronym": ["TMV"], "synonyms": []},
                            },
                            {
                                "id": 12227,
                                "name": "Tobacco mosaic virus",
                                "parent": 3432891,
                                "rank": "no rank",
                                "other_names": {"acronym": [], "synonyms": []},
                            },
                        ]
                    },
                    "molecule": {
                        "strandedness": "single",
                        "type": "RNA",
                        "topology": "linear",
                    },
                    "plan": {
                        "id": str(plan.id),
                        "segments": [
                            {
                                "id": str(plan.segments[0].id),
                                "length": SEGMENT_LENGTH,
                                "length_tolerance": empty_repo.settings.default_segment_length_tolerance,
                                "name": None,
                                "rule": "required",
                            }
                        ],
                    },
                    "promoted_accessions": [],
                },
                "id": 2,
                "query": {
                    "otu_id": str(otu.id),
                },
                "type": "CreateOTU",
            }

            assert empty_repo.last_id == 2

    def test_duplicate_taxid(self, empty_repo: Repo):
        """Test that creating an OTU with an existing taxid fails."""
        with empty_repo.lock():
            isolate_data = CreateIsolateData(
                id=uuid4(),
                name=None,
                taxid=12227,
                sequences=[],
            )

            empty_repo.create_otu(
                isolate=isolate_data,
                lineage=TMV_LINEAGE,
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MoleculeType.RNA,
                    topology=Topology.LINEAR,
                ),
                plan=Plan.new(
                    segments=[
                        Segment.new(
                            length=SEGMENT_LENGTH,
                            length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                            name=None,
                        )
                    ]
                ),
                promoted_accessions=set(),
            )

            isolate_data_2 = CreateIsolateData(
                id=uuid4(),
                name=None,
                taxid=12227,
                sequences=[],
            )

            with pytest.raises(
                ValueError,
                match="already contains taxid",
            ):
                empty_repo.create_otu(
                    isolate=isolate_data_2,
                    lineage=TMV_LINEAGE,
                    molecule=Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MoleculeType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    plan=Plan.new(
                        segments=[
                            Segment.new(
                                length=SEGMENT_LENGTH,
                                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                                name=None,
                            )
                        ]
                    ),
                    promoted_accessions=set(),
                )

    def test_plan_required_segment_warning(self, empty_repo: Repo):
        """Test that missing required segments raises a warning."""
        plan = Plan.new(
            segments=[
                Segment(
                    id=uuid4(),
                    length=SEGMENT_LENGTH,
                    length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                    name=None,
                    rule=SegmentRule.RECOMMENDED,
                )
            ]
        )

        isolate_data = CreateIsolateData(
            id=uuid4(),
            name=None,
            taxid=12227,
            sequences=[],
        )

        with (
            capture_logs() as captured_logs,
            empty_repo.lock(),
        ):
            otu = empty_repo.create_otu(
                isolate=isolate_data,
                lineage=TMV_LINEAGE,
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MoleculeType.RNA,
                    topology=Topology.LINEAR,
                ),
                plan=plan,
                promoted_accessions=set(),
            )

            assert otu

            isolate = empty_repo.create_isolate(
                otu.id,
                isolate_id=uuid4(),
                name=IsolateName(IsolateNameType.ISOLATE, "A"),
                taxid=12227,
                sequences=[
                    Sequence(
                        accession=Accession(key="TM000001", version=1),
                        definition="TMV",
                        segment=otu.plan.segments[0].id,
                        sequence="ACGTACGTACGTACG",
                    )
                ],
            )

            assert isolate

        assert any(
            log.get("warning_category") == "PlanWarning" for log in captured_logs
        )

    def test_ok(self, initialized_repo: Repo):
        """Test that a complete OTU can be created, validated and retrieved."""
        assert next(iter(initialized_repo.iter_otus())) is not None


class TestCreateIsolate:
    """Test the creation and addition of new isolates in Repo."""

    def test_ok(self, empty_repo: Repo):
        """Test creating an isolate.

        The method should return the expected ``Isolate`` create an event.
        """
        with empty_repo.lock():
            otu = init_otu(empty_repo)

            isolate_id = uuid4()
            isolate = empty_repo.create_isolate(
                otu.id,
                isolate_id=isolate_id,
                name=IsolateName(IsolateNameType.ISOLATE, "A"),
                taxid=12227,
                sequences=[],
            )

            assert isolate
            assert isinstance(isolate.id, UUID)
            assert isolate.sequences == []
            assert isolate.name
            assert isolate.name.value == "A"
            assert isolate.name.type == "isolate"

            with open(empty_repo.path.joinpath("src", "00000003.json")) as f:
                event = orjson.loads(f.read())

            del event["timestamp"]

            assert event == {
                "data": {
                    "id": str(isolate.id),
                    "name": {"type": "isolate", "value": "A"},
                    "taxid": 12227,
                    "sequences": [],
                },
                "id": 3,
                "query": {
                    "otu_id": str(otu.id),
                    "isolate_id": str(isolate.id),
                },
                "type": "CreateIsolate",
            }

            assert empty_repo.last_id == 3

    def test_create_unnamed(self, empty_repo: Repo):
        """Test that creating an isolate returns the expected ``Isolate`` object and
        creates the expected event file.
        """
        with empty_repo.lock():
            otu = init_otu(empty_repo)

            isolate_id = uuid4()
            isolate = empty_repo.create_isolate(
                otu.id,
                isolate_id=isolate_id,
                name=None,
                taxid=12227,
                sequences=[],
            )

            assert isolate
            assert isinstance(isolate.id, UUID)
            assert isolate.sequences == []
            assert isolate.name is None


class TestGetOTU:
    """Test the retrieval of OTU data."""

    def test_ok(self, empty_repo: Repo):
        """Test that getting an OTU returns the expected ``OTU`` object including
        two isolates with one sequence each.
        """
        monopartite_plan = Plan.new(
            segments=[
                Segment.new(
                    length=SEGMENT_LENGTH,
                    length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                    name=None,
                )
            ]
        )

        isolate_data = CreateIsolateData(
            id=uuid4(),
            name=None,
            taxid=12227,
            sequences=[],
        )

        with empty_repo.lock():
            otu = empty_repo.create_otu(
                isolate=isolate_data,
                lineage=TMV_LINEAGE,
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MoleculeType.RNA,
                    topology=Topology.LINEAR,
                ),
                plan=monopartite_plan,
                promoted_accessions=set(),
            )

            assert otu

            segment_id = otu.plan.segments[0].id

            isolate_a_id = uuid4()
            isolate_a = empty_repo.create_isolate(
                otu.id,
                isolate_id=isolate_a_id,
                name=IsolateName(IsolateNameType.ISOLATE, "A"),
                taxid=12227,
                sequences=[
                    Sequence(
                        accession=Accession(key="TM000001", version=1),
                        definition="TMV",
                        segment=segment_id,
                        sequence="ACGTACGTACGTACG",
                    )
                ],
            )

            isolate_b_id = uuid4()
            isolate_b = empty_repo.create_isolate(
                otu.id,
                isolate_id=isolate_b_id,
                name=IsolateName(IsolateNameType.ISOLATE, "B"),
                taxid=12227,
                sequences=[
                    Sequence(
                        accession=Accession(key="TN000001", version=1),
                        definition="TMV",
                        segment=segment_id,
                        sequence="TTACGTGGAGAGACC",
                    )
                ],
            )

            otu = empty_repo.get_otu(otu.id)

        otu_contents = [
            Isolate(
                id=isolate_data.id,
                name=None,
                taxid=12227,
                sequences=[],
            ),
            Isolate(
                id=isolate_a.id,
                name=IsolateName(type=IsolateNameType.ISOLATE, value="A"),
                taxid=12227,
                sequences=[
                    Sequence(
                        accession=Accession(key="TM000001", version=1),
                        definition="TMV",
                        segment=segment_id,
                        sequence="ACGTACGTACGTACG",
                    ),
                ],
            ),
            Isolate(
                id=isolate_b.id,
                name=IsolateName(type=IsolateNameType.ISOLATE, value="B"),
                taxid=12227,
                sequences=[
                    Sequence(
                        accession=Accession(key="TN000001", version=1),
                        definition="TMV",
                        segment=segment_id,
                        sequence="TTACGTGGAGAGACC",
                    ),
                ],
            ),
        ]

        assert otu
        assert (
            otu.model_dump()
            == OTU.model_validate(
                {
                    "id": otu.id,
                    "excluded_accessions": set(),
                    "promoted_accessions": set(),
                    "lineage": TMV_LINEAGE,
                    "molecule": Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MoleculeType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    "plan": Plan(
                        id=monopartite_plan.id,
                        segments=[
                            Segment(
                                id=segment_id,
                                length=SEGMENT_LENGTH,
                                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                                name=None,
                                rule=SegmentRule.REQUIRED,
                            )
                        ],
                    ),
                    "isolates": otu_contents,
                }
            ).model_dump()
        )
        assert empty_repo.last_id == 4

    def test_retrieve_nonexistent_otu(self, initialized_repo: Repo):
        """Test that getting an OTU that does not exist returns ``None``."""
        assert initialized_repo.get_otu(uuid4()) is None

    def test_accessions(self, initialized_repo: Repo):
        """Test that the `accessions` property returns the expected accessions."""
        assert next(initialized_repo.iter_otus()).accessions == {"TM000001"}

    def test_blocked_accessions(self, initialized_repo: Repo):
        """Test that the `blocked_accessions` property returns the expected set of
        accessions.
        """
        otu = initialized_repo.get_otu_by_taxid(3432891)

        assert otu

        excludable_accessions = {"GR33333", "TL44322"}

        with initialized_repo.lock():
            isolate = initialized_repo.create_isolate(
                otu.id,
                isolate_id=uuid4(),
                name=IsolateName(type=IsolateNameType.ISOLATE, value="B"),
                taxid=12227,
                sequences=[
                    Sequence(
                        accession=Accession(key="TN000001", version=1),
                        definition="TMV",
                        segment=otu.plan.segments[0].id,
                        sequence="TTACGTGGAGAGACC",
                    )
                ],
            )

            assert isolate

        otu_before = initialized_repo.get_otu(otu.id)

        assert otu_before
        assert otu_before.blocked_accessions == {
            "TM000001",
            "TN000001",
        }

        with initialized_repo.lock():
            initialized_repo.exclude_accessions(otu.id, excludable_accessions)

        otu_after = initialized_repo.get_otu(otu.id)

        assert otu_after
        assert otu_after.blocked_accessions == {
            "TM000001",
            "TN000001",
            "GR33333",
            "TL44322",
        }


def test_get_otu_id_from_isolate_id(initialized_repo: Repo):
    """Test that the OTU id can be retrieved from a isolate ID contained within."""
    otu = next(initialized_repo.iter_otus())

    isolate = otu.isolates[0]

    assert initialized_repo.get_otu_id_by_isolate_id(isolate.id) == otu.id


class TestGetIsolate:
    def test_ok(self, initialized_repo: Repo):
        """Test that getting an isolate returns the expected ``Isolate``."""
        otu = next(initialized_repo.iter_otus())

        for isolate in otu.isolates:
            assert otu.get_isolate(isolate.id) in otu.isolates


class TestCreateIsolateValidation:
    """Test the validation of new added isolates."""

    @pytest.mark.parametrize("bad_sequence", ["ACGTACGTA", "ACGTACGTACGTACGTTTTTT"])
    def test_bad_segment_length_fail(self, initialized_repo: Repo, bad_sequence: str):
        """Test that a new isolate with segments that are too long or short raises a validation error."""
        otu_before = next(iter(initialized_repo.iter_otus()))

        with initialized_repo.lock():
            with pytest.raises(ValueError, match="Event validation failed"):
                initialized_repo.create_isolate(
                    otu_before.id,
                    isolate_id=uuid4(),
                    name=IsolateName(IsolateNameType.ISOLATE, "B"),
                    taxid=12227,
                    sequences=[
                        Sequence(
                            accession=Accession(key="MK000001", version=1),
                            definition="TMV",
                            segment=otu_before.plan.segments[0].id,
                            sequence=bad_sequence,
                        )
                    ],
                )


class TestExcludeAccessions:
    """Test that excluded accessions are reflected in events.

    Excluded accessions are exempt from future queries to the NCBI Nucleotide.
    """

    def test_ok(self, initialized_repo: Repo):
        """Test that Repo.exclude_accessions() writes the correct event."""
        otu = initialized_repo.get_otu_by_taxid(3432891)
        id_at_creation = initialized_repo.last_id

        assert otu
        assert otu.excluded_accessions == set()
        assert id_at_creation == 2

        accessions = {"TM100021", "TM100022", "TM100023"}

        with initialized_repo.lock():
            initialized_repo.exclude_accessions(otu.id, accessions)

        otu_before = initialized_repo.get_otu(otu.id)

        assert otu_before
        assert otu_before.excluded_accessions == accessions
        assert initialized_repo.last_id == id_at_creation + 1 == 3

        with open(
            initialized_repo.path / "src" / f"{initialized_repo.last_id:08}.json"
        ) as f:
            event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "accessions": ["TM100021", "TM100022", "TM100023"],
                "action": "exclude",
            },
            "id": 3,
            "query": {
                "otu_id": str(otu.id),
            },
            "type": "UpdateExcludedAccessions",
        }

        with initialized_repo.lock():
            initialized_repo.exclude_accessions(otu.id, {"TM100024"})

        assert initialized_repo.last_id == id_at_creation + 2 == 4

        otu_after_second = initialized_repo.get_otu(otu.id)
        assert otu_after_second
        assert otu_after_second.excluded_accessions == accessions | {"TM100024"}

    def test_existing_accession(self, initialized_repo: Repo):
        """Test that excluding an accession moves its isolate to excluded_isolates."""
        otu = initialized_repo.get_otu_by_taxid(3432891)
        assert otu
        assert otu.excluded_accessions == set()
        assert len(otu.isolates) == 1
        assert len(otu.excluded_isolates) == 0

        # Add a second isolate so we can exclude one without violating min_length=1
        with initialized_repo.lock():
            initialized_repo.create_isolate(
                otu_id=otu.id,
                isolate_id=uuid4(),
                name=IsolateName(IsolateNameType.ISOLATE, "B"),
                taxid=3432891,
                sequences=[
                    Sequence(
                        accession=Accession(key="TN000001", version=1),
                        definition="Second isolate",
                        segment=otu.plan.segments[0].id,
                        sequence="ACGTACGTACGTACG",
                    )
                ],
            )

        otu = initialized_repo.get_otu_by_taxid(3432891)
        assert len(otu.isolates) == 2
        assert len(otu.excluded_isolates) == 0

        # Now exclude the first accession
        accession_to_exclude = "TM000001"
        id_before_exclude = initialized_repo.last_id

        with initialized_repo.lock():
            initialized_repo.exclude_accessions(otu.id, {accession_to_exclude})

        # Event should be written
        assert initialized_repo.last_id == id_before_exclude + 1

        otu = initialized_repo.get_otu_by_taxid(3432891)
        assert otu

        # Verify the accession is excluded
        assert accession_to_exclude in otu.excluded_accessions

        # Verify one isolate is active, one is excluded
        assert len(otu.isolates) == 1
        assert len(otu.excluded_isolates) == 1

        # Verify the excluded isolate contains the excluded accession
        assert accession_to_exclude in otu.excluded_isolates[0].accessions

        # Verify the active isolate does not
        assert accession_to_exclude not in otu.isolates[0].accessions

    def test_already_excluded(self, initialized_repo: Repo):
        """Test that an attempted redundant exclusion does not create a new event."""
        repo = initialized_repo

        otu_before = initialized_repo.get_otu_by_taxid(3432891)

        assert otu_before
        assert otu_before == repo.get_otu(otu_before.id)
        assert otu_before.excluded_accessions == set()

        event_id_after_creation = repo.last_id

        accessions = {"TM100021", "TM100022", "TM100023"}

        with repo.lock():
            repo.exclude_accessions(otu_before.id, accessions)

        # Event ID should increment for the first successful exclusion.
        assert repo.last_id == event_id_after_creation + 1

        with repo.lock():
            repo.exclude_accessions(otu_before.id, {"TM100021"})

        otu_after = repo.get_otu(otu_before.id)

        assert otu_after
        assert otu_after.excluded_accessions == accessions

        # Event ID should not increment.
        assert repo.last_id == event_id_after_creation + 1

    def test_partially_already_excluded(self, initialized_repo: Repo):
        """Test that a partially redundant list of exclusions creates a new event with
        the redundant accession omitted.
        """
        otu_before = initialized_repo.get_otu_by_taxid(3432891)

        assert otu_before

        event_id_after_creation = initialized_repo.last_id

        accessions = {"TM100021", "TM100022", "TM100023"}

        with initialized_repo.lock():
            initialized_repo.exclude_accessions(otu_before.id, accessions)

        otu_after_first = initialized_repo.get_otu(otu_before.id)

        assert otu_after_first
        assert otu_after_first.excluded_accessions == accessions
        assert initialized_repo.last_id == event_id_after_creation + 1

        with initialized_repo.lock():
            initialized_repo.exclude_accessions(otu_before.id, {"TM100023", "TM100024"})

        otu_after_second = initialized_repo.get_otu(otu_before.id)

        assert otu_after_second
        assert otu_after_second.excluded_accessions == (
            otu_after_first.excluded_accessions | {"TM100024"}
        )
        assert initialized_repo.last_id == event_id_after_creation + 2

        with open(
            initialized_repo.path.joinpath("src", f"{initialized_repo.last_id:08}.json")
        ) as f:
            event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "accessions": ["TM100024"],
                "action": "exclude",
            },
            "id": initialized_repo.last_id,
            "query": {
                "otu_id": str(otu_before.id),
            },
            "type": "UpdateExcludedAccessions",
        }


class TestAllowAccessions:
    """Test that accessions allowed back into the OTU are no longer contained
    in the excluded accessions set.
    """

    def test_ok(self, initialized_repo: Repo):
        """Test that Repo.allow_accessions() produces the correct event and
        creates the expected OTU.excluded_accessions set.
        """
        target_repo = initialized_repo

        otu = target_repo.get_otu_by_taxid(3432891)

        assert otu
        assert (id_at_creation := target_repo.last_id) == 2

        accessions = {"TM100021", "TM100022", "TM100023"}

        with target_repo.lock():
            target_repo.exclude_accessions(otu.id, accessions)

            assert target_repo.last_id == id_at_creation + 1 == 3

            otu = target_repo.get_otu(otu.id)

            assert otu
            assert otu.excluded_accessions == accessions

            target_repo.allow_accessions(otu.id, ["TM100021", "TM100022"])

            assert target_repo.last_id == id_at_creation + 2 == 4

        otu_after = target_repo.get_otu(otu.id)

        assert otu_after

        with open(
            target_repo.path.joinpath("src", f"0000000{target_repo.last_id}.json")
        ) as f:
            event = orjson.loads(f.read())

            del event["timestamp"]

            assert event == {
                "data": {
                    "accessions": ["TM100021", "TM100022"],
                    "action": "allow",
                },
                "id": 4,
                "query": {
                    "otu_id": str(otu_after.id),
                },
                "type": "UpdateExcludedAccessions",
            }

        assert otu_after.excluded_accessions == {"TM100023"}

    def test_skip_redundant_accessions(self, initialized_repo: Repo):
        """Test that an event only gets written if the accession exists
        in the exclusion list, avoiding the creation of a redundant event.
        """
        otu = initialized_repo.get_otu_by_taxid(3432891)

        assert otu

        accessions = {"TM100021", "TM100022", "TM100023"}

        with initialized_repo.lock():
            initialized_repo.exclude_accessions(otu.id, accessions)

        assert (id_after_first_exclusion := initialized_repo.last_id) == 3

        otu = initialized_repo.get_otu(otu.id)

        assert otu
        assert otu.excluded_accessions == accessions

        with initialized_repo.lock():
            initialized_repo.allow_accessions(otu.id, ["TM100024"])

        assert initialized_repo.last_id == id_after_first_exclusion == 3
        assert not initialized_repo.path.joinpath("src", "00000004.json").exists()

        otu = initialized_repo.get_otu(otu.id)

        assert otu
        assert otu.excluded_accessions == accessions

        with initialized_repo.lock():
            initialized_repo.allow_accessions(otu.id, accessions)

        otu = initialized_repo.get_otu(otu.id)

        assert otu
        assert otu.excluded_accessions == set()
        assert initialized_repo.last_id == id_after_first_exclusion + 1 == 4
        assert initialized_repo.path.joinpath("src", "00000004.json").exists()

    def test_restore_isolates(self, initialized_repo: Repo):
        """Test that allowing an excluded accession restores its isolate."""
        otu = initialized_repo.get_otu_by_taxid(3432891)
        assert otu
        assert len(otu.isolates) == 1
        assert len(otu.excluded_isolates) == 0

        # Add a second isolate
        with initialized_repo.lock():
            initialized_repo.create_isolate(
                otu_id=otu.id,
                isolate_id=uuid4(),
                name=IsolateName(IsolateNameType.ISOLATE, "B"),
                taxid=3432891,
                sequences=[
                    Sequence(
                        accession=Accession(key="TN000001", version=1),
                        definition="Second isolate",
                        segment=otu.plan.segments[0].id,
                        sequence="ACGTACGTACGTACG",
                    )
                ],
            )

        otu = initialized_repo.get_otu_by_taxid(3432891)
        assert len(otu.isolates) == 2
        assert len(otu.excluded_isolates) == 0

        # Exclude one isolate
        accession_to_exclude = "TM000001"
        with initialized_repo.lock():
            initialized_repo.exclude_accessions(otu.id, {accession_to_exclude})

        otu = initialized_repo.get_otu_by_taxid(3432891)
        assert len(otu.isolates) == 1
        assert len(otu.excluded_isolates) == 1
        assert accession_to_exclude in otu.excluded_accessions
        assert accession_to_exclude in otu.excluded_isolates[0].accessions

        # Allow the accession back
        with initialized_repo.lock():
            initialized_repo.allow_accessions(otu.id, [accession_to_exclude])

        otu = initialized_repo.get_otu_by_taxid(3432891)

        # Verify the isolate is restored
        assert len(otu.isolates) == 2
        assert len(otu.excluded_isolates) == 0
        assert accession_to_exclude not in otu.excluded_accessions

        # Verify both isolates are active
        all_accessions = {acc for isolate in otu.isolates for acc in isolate.accessions}
        assert accession_to_exclude in all_accessions


class TestDeleteIsolate:
    """Test that an isolate can be deleted from an OTU."""

    def test_ok(self, initialized_repo: Repo):
        """Test basic functionality."""
        otu_before = initialized_repo.get_otu_by_taxid(3432891)

        assert otu_before

        with initialized_repo.lock():
            isolate_b_id = uuid4()
            isolate_b = initialized_repo.create_isolate(
                otu_before.id,
                isolate_id=isolate_b_id,
                name=IsolateName(IsolateNameType.ISOLATE, "B"),
                taxid=12227,
                sequences=[
                    Sequence(
                        accession=Accession(key="TN000001", version=1),
                        definition="TMV",
                        segment=otu_before.plan.segments[0].id,
                        sequence="TTACGTGGAGAGACC",
                    )
                ],
            )

            assert isolate_b
            assert "TN000001" in isolate_b.accessions

        event_id_before_delete = initialized_repo.last_id

        otu_before = initialized_repo.get_otu(otu_before.id)

        assert otu_before
        assert isolate_b.id in otu_before.isolate_ids
        assert "TN000001" in otu_before.accessions
        assert event_id_before_delete == initialized_repo.last_id == 3

        with initialized_repo.lock():
            initialized_repo.delete_isolate(
                otu_before.id,
                isolate_b.id,
                message="Testing deletion",
            )

        otu_after = initialized_repo.get_otu(otu_before.id)

        assert otu_after
        assert otu_after.get_isolate(isolate_b.id) is None
        assert otu_before != otu_after
        assert len(otu_after.isolates) == len(otu_before.isolates) - 1
        assert isolate_b.id not in otu_after.isolate_ids
        assert isolate_b.accessions not in otu_after.accessions

        assert initialized_repo.last_id == event_id_before_delete + 1


class TestMalformedEvent:
    """Test that malformed events cannot be rehydrated."""

    def test_bad_event_typing(self, initialized_repo: Repo):
        """Test that an event with an invalid event type discriminator does not attempt
        to rehydrate.
        """
        file_path = initialized_repo.path.joinpath("src", "00000002.json")

        with open(file_path, "rb") as f:
            event = orjson.loads(f.read())

        otu = initialized_repo.get_otu_by_taxid(3432891)

        assert type(otu) is OTU

        event["type"] = "MalformedEvent"

        with open(file_path, "wb") as f:
            f.write(orjson.dumps(event))

        with pytest.raises(ValueError, match="Unknown event type: MalformedEvent"):
            initialized_repo.get_otu_by_taxid(3432891)

    def test_bad_event_data(self, initialized_repo: Repo):
        """Test that an event with bad data cannot be rehydrated."""
        path = initialized_repo.path.joinpath("src", "00000002.json")

        with open(path, "rb") as f:
            event = orjson.loads(f.read())

        with initialized_repo.lock():
            otu = initialized_repo.get_otu_by_taxid(3432891)

        assert type(otu) is OTU

        # Corrupt the lineage data
        event["data"]["lineage"]["taxa"] = "popcorn"

        with open(path, "wb") as f:
            f.write(orjson.dumps(event))

        with (
            initialized_repo.lock(),
            pytest.raises(
                ValueError,
                match="Input should be a valid list",
            ),
        ):
            initialized_repo.get_otu_by_taxid(3432891)
