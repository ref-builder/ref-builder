import warnings
from pathlib import Path
from uuid import UUID, uuid4

import orjson
import pytest
from structlog.testing import capture_logs

from ref_builder.errors import InvalidInputError
from ref_builder.models import Molecule, MolType, Strandedness, Topology
from ref_builder.otu.builders.isolate import IsolateBuilder
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.builders.sequence import SequenceBuilder
from ref_builder.plan import Plan, Segment, SegmentRule
from ref_builder.repo import GITIGNORE_CONTENTS, Repo
from ref_builder.utils import (
    Accession,
    DataType,
    IsolateName,
    IsolateNameType,
    OTUDeletedWarning,
)

SEGMENT_LENGTH = 15


def init_otu_with_contents(repo: Repo, otu: OTUBuilder):
    """Add an unvalidated OTU to a given repo."""
    isolate_to_copy = otu.isolates[0]

    assert isolate_to_copy

    otu_before = repo.create_otu(
        acronym=otu.acronym,
        molecule=otu.molecule,
        name=otu.name,
        plan=otu.plan,
        taxid=otu.taxid,
    )
    assert otu_before is not None

    isolate_before = repo.create_isolate(
        otu_before.id,
        name=isolate_to_copy.name,
    )
    assert isolate_before is not None

    for sequence in isolate_to_copy.sequences:
        if len(otu_before.plan.segments):
            matching_segment = otu_before.plan.segments[0].id
        else:
            segment_by_name = otu_before.plan.get_segment_by_name_key(
                str(sequence.segment)
            )
            assert segment_by_name is not None
            matching_segment = segment_by_name.id

        sequence_init = repo.create_sequence(
            otu_before.id,
            accession=sequence.accession,
            definition=sequence.definition,
            segment=matching_segment,
            sequence=sequence.sequence,
        )
        assert sequence_init is not None

        repo.link_sequence(otu_before.id, isolate_before.id, sequence_init.id)

    result = repo.get_otu(otu_before.id)

    assert result is not None

    return result


@pytest.fixture
def initialized_repo(tmp_path: Path):
    """Return a pre-initialized mock Repo."""
    repo = Repo.new(
        DataType.GENOME,
        "Generic Viruses",
        tmp_path / "initialized_repo",
        "virus",
    )

    with repo.lock(), repo.use_transaction():
        otu = repo.create_otu(
            acronym="TMV",
            molecule=Molecule(
                strandedness=Strandedness.SINGLE,
                type=MolType.RNA,
                topology=Topology.LINEAR,
            ),
            name="Tobacco mosaic virus",
            plan=Plan.new(
                [
                    Segment.new(
                        length=SEGMENT_LENGTH,
                        length_tolerance=repo.settings.default_segment_length_tolerance,
                        name=None,
                    )
                ]
            ),
            taxid=12242,
        )
        assert otu is not None

        sequence_1 = repo.create_sequence(
            otu.id,
            "TM000001.1",
            "TMV",
            otu.plan.segments[0].id,
            "ACGTACGTACGTACG",
        )
        assert sequence_1 is not None

        isolate_a = repo.create_isolate(
            otu.id,
            IsolateName(IsolateNameType.ISOLATE, "A"),
        )

        assert isolate_a is not None

        repo.link_sequence(otu.id, isolate_a.id, sequence_1.id)

    return repo


def init_otu(repo: Repo) -> OTUBuilder:
    """Create an empty OTU."""
    result = repo.create_otu(
        acronym="TMV",
        molecule=Molecule(
            strandedness=Strandedness.SINGLE,
            type=MolType.RNA,
            topology=Topology.LINEAR,
        ),
        name="Tobacco mosaic virus",
        plan=Plan.new(
            segments=[
                Segment.new(
                    length=SEGMENT_LENGTH,
                    length_tolerance=repo.settings.default_segment_length_tolerance,
                    name=None,
                )
            ]
        ),
        taxid=12242,
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

        assert empty_repo.meta.data_type == DataType.GENOME
        assert empty_repo.meta.name == "Generic Viruses"
        assert empty_repo.meta.organism == "virus"

        assert empty_repo.settings.default_segment_length_tolerance == 0.03

        assert (empty_repo.path / ".gitignore").exists()

        with open(empty_repo.path / ".gitignore") as f:
            assert f.read() == "\n".join(GITIGNORE_CONTENTS) + "\n"

    def test_alternate_settings(self, tmp_path: Path):
        """Test retrieval of non-default settings."""
        repo = Repo.new(
            DataType.GENOME,
            "Generic Viruses",
            tmp_path / "alt_setting_repo",
            "virus",
            default_segment_length_tolerance=0.05,
        )

        assert repo.settings.default_segment_length_tolerance == 0.05


class TestCreateOTU:
    def test_empty_ok(self, empty_repo: Repo):
        """Test that creating an OTU returns the expected ``OTUBuilder`` object and creates
        the expected event file.
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

        with empty_repo.lock(), empty_repo.use_transaction():
            otu = empty_repo.create_otu(
                acronym="TMV",
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                name="Tobacco mosaic virus",
                plan=plan,
                taxid=12242,
            )
            assert otu is not None

            assert otu == OTUBuilder(
                id=otu.id,
                acronym="TMV",
                excluded_accessions=set(),
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                name="Tobacco mosaic virus",
                plan=Plan(
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
                taxid=12242,
                isolates=[],
            )

            with open(empty_repo.path.joinpath("src", "00000002.json")) as f:
                event = orjson.loads(f.read())

            del event["timestamp"]

            assert event == {
                "data": {
                    "id": str(otu.id),
                    "acronym": "TMV",
                    "molecule": {
                        "strandedness": "single",
                        "type": "RNA",
                        "topology": "linear",
                    },
                    "name": "Tobacco mosaic virus",
                    "plan": {
                        "id": str(plan.id),
                        "segments": [
                            {
                                "id": str(plan.segments[0].id),
                                "length": SEGMENT_LENGTH,
                                "length_tolerance": empty_repo.settings.default_segment_length_tolerance,
                                "name": None,
                                "rule": SegmentRule.REQUIRED,
                            }
                        ],
                    },
                    "taxid": 12242,
                },
                "id": 2,
                "query": {
                    "otu_id": str(otu.id),
                },
                "type": "CreateOTU",
            }

            assert empty_repo.last_id == 2

    def test_duplicate_name_fail(self, empty_repo: Repo):
        """Test that creating an OTU with a name that already exists raises a
        ``ValueError``.
        """
        with empty_repo.lock(), empty_repo.use_transaction():
            empty_repo.create_otu(
                acronym="TMV",
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                name="Tobacco mosaic virus",
                plan=Plan.new(
                    segments=[
                        Segment.new(
                            length=SEGMENT_LENGTH,
                            length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                            name=None,
                        )
                    ]
                ),
                taxid=12242,
            )

            with pytest.raises(
                ValueError,
                match="An OTU with the name 'Tobacco mosaic virus' already exists",
            ):
                empty_repo.create_otu(
                    acronym="TMV",
                    molecule=Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MolType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    name="Tobacco mosaic virus",
                    plan=Plan.new(
                        segments=[
                            Segment.new(
                                length=SEGMENT_LENGTH,
                                length_tolerance=empty_repo.settings.default_segment_length_tolerance,
                                name=None,
                            )
                        ]
                    ),
                    taxid=438782,
                )

    def test_plan_required_segment_warning(self, empty_repo: Repo):
        """Test that creating an OTU without required segments raises a warning
        once the transaction exits.
        """
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
        with (
            capture_logs() as captured_logs,
            empty_repo.lock(),
            empty_repo.use_transaction(),
        ):
            otu = empty_repo.create_otu(
                acronym="TMV",
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                name="Tobacco mosaic virus",
                plan=plan,
                taxid=12242,
            )

            assert otu

            sequence = empty_repo.create_sequence(
                otu.id,
                "TM000001.1",
                "TMV",
                otu.plan.segments[0].id,
                "ACGTACGTACGTACG",
            )

            assert sequence

            isolate = empty_repo.create_isolate(
                otu.id,
                IsolateName(IsolateNameType.ISOLATE, "A"),
            )

            assert isolate

            empty_repo.link_sequence(otu.id, isolate.id, sequence.id)

        assert any(
            log.get("warning_category") == "PlanWarning" for log in captured_logs
        )


class TestCreateIsolate:
    """Test the creation and addition of new isolates in Repo."""

    def test_ok(self, empty_repo: Repo):
        """Test that creating an isolate returns the expected ``IsolateBuilder`` object and
        creates the expected event file.
        """
        with empty_repo.lock(), empty_repo.use_transaction():
            otu = init_otu(empty_repo)

            isolate = empty_repo.create_isolate(
                otu.id,
                IsolateName(IsolateNameType.ISOLATE, "A"),
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
        """Test that creating an isolate returns the expected ``IsolateBuilder`` object and
        creates the expected event file.
        """
        with empty_repo.lock(), empty_repo.use_transaction():
            otu = init_otu(empty_repo)

            isolate = empty_repo.create_isolate(
                otu.id,
                None,
            )

            assert isolate
            assert isinstance(isolate.id, UUID)
            assert isolate.sequences == []
            assert isolate.name is None


def test_create_sequence(empty_repo: Repo):
    """Test that creating a sequence returns the expected ``SequenceBuilder`` object and
    creates the expected event file.
    """
    with empty_repo.lock(), empty_repo.use_transaction():
        otu = init_otu(empty_repo)

        sequence = empty_repo.create_sequence(
            otu.id,
            "TM000001.1",
            "TMV",
            otu.plan.segments[0].id,
            "ACGTACGTACGTACG",
        )

        assert sequence is not None

        assert sequence == SequenceBuilder(
            id=sequence.id,
            accession=Accession(key="TM000001", version=1),
            definition="TMV",
            segment=otu.plan.segments[0].id,
            sequence="ACGTACGTACGTACG",
        )

        with open(empty_repo.path.joinpath("src", "00000003.json")) as f:
            event = orjson.loads(f.read())

        del event["timestamp"]

        assert event == {
            "data": {
                "id": str(sequence.id),
                "accession": {"key": "TM000001", "version": 1},
                "definition": "TMV",
                "segment": str(otu.plan.segments[0].id),
                "sequence": "ACGTACGTACGTACG",
            },
            "id": 3,
            "query": {
                "otu_id": str(otu.id),
                "sequence_id": str(sequence.id),
            },
            "type": "CreateSequence",
        }

        assert empty_repo.last_id == 3


class TestGetOTU:
    """Test the retrieval of OTU data."""

    def test_ok(self, empty_repo: Repo):
        """Test that getting an OTU returns the expected ``OTUBuilder`` object including
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

        with empty_repo.lock():
            with empty_repo.use_transaction():
                otu = empty_repo.create_otu(
                    acronym="TMV",
                    molecule=Molecule(
                        strandedness=Strandedness.SINGLE,
                        type=MolType.RNA,
                        topology=Topology.LINEAR,
                    ),
                    name="Tobacco mosaic virus",
                    taxid=12242,
                    plan=monopartite_plan,
                )

                assert otu

                segment_id = otu.plan.segments[0].id

                sequence_1 = empty_repo.create_sequence(
                    otu.id,
                    "TM000001.1",
                    "TMV",
                    segment_id,
                    "ACGTACGTACGTACG",
                )

                isolate_a = empty_repo.create_isolate(
                    otu.id,
                    IsolateName(IsolateNameType.ISOLATE, "A"),
                )

                empty_repo.link_sequence(
                    otu.id, isolate_a.id, sequence_id=sequence_1.id
                )

                sequence_2 = empty_repo.create_sequence(
                    otu.id,
                    "TN000001.1",
                    "TMV",
                    segment_id,
                    "TTACGTGGAGAGACC",
                )

                isolate_b = empty_repo.create_isolate(
                    otu.id,
                    IsolateName(IsolateNameType.ISOLATE, "B"),
                )

                empty_repo.link_sequence(
                    otu.id, isolate_b.id, sequence_id=sequence_2.id
                )

            otu = empty_repo.get_otu(otu.id)

        otu_contents = [
            IsolateBuilder(
                id=isolate_a.id,
                name=IsolateName(type=IsolateNameType.ISOLATE, value="A"),
                sequences=[
                    SequenceBuilder(
                        id=otu.isolates[0].sequences[0].id,
                        accession=Accession(key="TM000001", version=1),
                        definition="TMV",
                        segment=segment_id,
                        sequence="ACGTACGTACGTACG",
                    ),
                ],
            ),
            IsolateBuilder(
                id=isolate_b.id,
                name=IsolateName(type=IsolateNameType.ISOLATE, value="B"),
                sequences=[
                    SequenceBuilder(
                        id=otu.isolates[1].sequences[0].id,
                        accession=Accession(key="TN000001", version=1),
                        definition="TMV",
                        segment=segment_id,
                        sequence="TTACGTGGAGAGACC",
                    ),
                ],
            ),
        ]

        assert (
            otu.model_dump()
            == OTUBuilder(
                id=otu.id,
                acronym="TMV",
                excluded_accessions=set(),
                molecule=Molecule(
                    strandedness=Strandedness.SINGLE,
                    type=MolType.RNA,
                    topology=Topology.LINEAR,
                ),
                name="Tobacco mosaic virus",
                plan=Plan(
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
                taxid=12242,
                isolates=otu_contents,
            ).model_dump()
        )

        assert empty_repo.last_id == 8

    def test_acronym_ok(self, initialized_repo: Repo):
        """Test that getting an OTU ID from the exact acronym of the OTU returns the expected ID."""
        otu = next(initialized_repo.iter_otus())

        assert initialized_repo.get_otu_id_by_acronym("TMV") == otu.id

    def test_acronym_fail(self, initialized_repo: Repo):
        """Test that a non-matching acronym cannot retrieve the OTU ID"""
        assert initialized_repo.get_otu_id_by_acronym("TM") is None

        assert initialized_repo.get_otu_id_by_acronym("TMVV") is None

    def test_acronym_empty_fail(self, initialized_repo: Repo):
        """Test that an attempt to search indexed OTUs using an empty acronym
        logs an error message and returns None.
        """
        with capture_logs() as logs:
            assert initialized_repo.get_otu_id_by_acronym("") is None

        assert any(["Bad input" in log["event"] for log in logs])

    def test_partial_ok(self, initialized_repo: Repo):
        """Test that getting an OTU ID starting with a truncated 8-character portion
        returns an ID.
        """
        otu = next(initialized_repo.iter_otus())

        partial_id = str(otu.id)[:8]

        assert initialized_repo.get_otu_id_by_partial(partial_id) == otu.id

    def test_partial_too_short(self, initialized_repo: Repo):
        """Test that getting an OTU ID starting with a truncated 7-character portion
        does not return an ID.
        """
        otu = next(initialized_repo.iter_otus())

        partial_id = str(otu.id)[:7]

        with pytest.raises(
            InvalidInputError,
            match="Partial ID segment must be at least 8 characters long",
        ):
            initialized_repo.get_otu_id_by_partial(partial_id)

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
        otu = initialized_repo.get_otu_by_taxid(12242)

        assert otu

        excludable_accessions = {"GR33333", "TL44322"}

        with initialized_repo.lock(), initialized_repo.use_transaction():
            sequence = initialized_repo.create_sequence(
                otu.id,
                "TN000001.1",
                "TMV",
                otu.plan.segments[0].id,
                "TTACGTGGAGAGACC",
            )

            isolate = initialized_repo.create_isolate(
                otu.id,
                IsolateName(type=IsolateNameType.ISOLATE, value="B"),
            )

            assert isolate

            initialized_repo.link_sequence(otu.id, isolate.id, sequence_id=sequence.id)

        assert initialized_repo.get_otu(otu.id).blocked_accessions == {
            "TM000001",
            "TN000001",
        }

        with initialized_repo.lock(), initialized_repo.use_transaction():
            initialized_repo.exclude_accessions(otu.id, excludable_accessions)

        assert initialized_repo.get_otu(otu.id).blocked_accessions == {
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


class TestDeleteOTU:
    """Test OTU deletion."""

    def test_ok(self, initialized_repo: Repo):
        """Test that .delete_otu() returns True if OTU can be deleted."""
        otu_before = next(initialized_repo.iter_otus())

        assert isinstance(otu_before, OTUBuilder)

        with initialized_repo.lock(), initialized_repo.use_transaction():
            assert initialized_repo.delete_otu(
                otu_before.id,
                rationale="Testing OTU deletion",
                replacement_otu_id=None,
            )

        assert initialized_repo.get_otu_id_by_taxid(otu_before.taxid) is None

        with warnings.catch_warnings(record=True) as warning_list:
            assert initialized_repo.get_otu(otu_before.id) is None

            for warning_msg in warning_list:
                assert warning_msg.category.__name__ == "OTUDeletedWarning"

        assert list(initialized_repo.iter_otus()) == []

    def test_not_found(self, initialized_repo: Repo):
        """Test that .delete_otu() returns False if OTU didn't exist to begin with."""
        dummy_otu_id = uuid4()

        with initialized_repo.lock(), initialized_repo.use_transaction():
            assert not initialized_repo.delete_otu(
                dummy_otu_id,
                rationale="Testing OTU deletion",
                replacement_otu_id=None,
            )

        with warnings.catch_warnings(record=True) as warning_list:
            assert initialized_repo.get_otu(dummy_otu_id) is None

        assert not warning_list

    def test_already_deleted(self, initialized_repo: Repo):
        """Test that .delete_otu() returns False and warns if an already-deleted OTU
        is marked for deletion a second time.
        """
        otu_before = next(initialized_repo.iter_otus())

        assert isinstance(otu_before, OTUBuilder)

        with initialized_repo.lock(), initialized_repo.use_transaction():
            assert initialized_repo.delete_otu(
                otu_before.id,
                rationale="Testing OTU deletion of nonexistent OTU",
                replacement_otu_id=None,
            )

        assert list(initialized_repo.iter_otus()) == []

        with initialized_repo.lock(), initialized_repo.use_transaction():
            result = initialized_repo.delete_otu(
                otu_before.id,
                rationale="Testing OTU deletion of already-deleted OTU",
                replacement_otu_id=None,
            )

            assert result is False

        with warnings.catch_warnings(
            category=OTUDeletedWarning, record=True
        ) as warning_messages:
            assert initialized_repo.get_otu(otu_before.id) is None

        assert len(warning_messages) == 1
        warning_message = warning_messages[0]

        assert (
            str(warning_message.message)
            == f"OTU {otu_before.id} has already been deleted."
        )
        assert warning_message.category == OTUDeletedWarning


class TestGetIsolate:
    def test_by_id(self, initialized_repo: Repo):
        """Test that getting an isolate returns the expected ``IsolateBuilder`` object."""
        otu = next(initialized_repo.iter_otus())

        for isolate in otu.isolates:
            assert otu.get_isolate(isolate.id) in otu.isolates

    def test_by_name(self, initialized_repo: Repo):
        """Test that getting an isolate ID by name returns the expected ID."""
        otu = next(initialized_repo.iter_otus())

        isolate_ids = {isolate.id for isolate in otu.isolates}

        assert (
            otu.get_isolate_id_by_name(
                IsolateName(type=IsolateNameType.ISOLATE, value="A"),
            )
            in isolate_ids
        )

    def test_get_with_unnamed_isolate(self, initialized_repo: Repo):
        """Test that getting an OTU with an unnamed isolate ID behaves as expected."""
        otu = next(initialized_repo.iter_otus())

        isolate_before = otu.isolates[0]

        with initialized_repo.lock():
            with initialized_repo.use_transaction():
                sequence = initialized_repo.create_sequence(
                    otu.id,
                    accession="NP000001.1",
                    definition="TMV B",
                    segment=otu.plan.segments[0].id,
                    sequence="TTGACCACGTGGAGA",
                )

                assert sequence

                isolate_unnamed = initialized_repo.create_isolate(
                    otu.id,
                    None,
                )

                assert isolate_unnamed

                initialized_repo.link_sequence(
                    otu.id, isolate_unnamed.id, sequence_id=sequence.id
                )

            otu_after = next(initialized_repo.iter_otus())

        assert len(otu_after.isolate_ids) == len(otu.isolate_ids) + 1
        assert otu_after.accessions == {"TM000001", "NP000001"}
        assert otu_after.isolate_ids == {isolate_before.id, isolate_unnamed.id}

        isolate_unnamed_after = otu_after.get_isolate(isolate_unnamed.id)

        assert isolate_unnamed_after
        assert isolate_unnamed_after.name is None
        assert isolate_unnamed_after.accessions == {"NP000001"}


def test_get_isolate_id_from_partial(initialized_repo: Repo):
    """Test that an isolate id can be retrieved from a truncated ``partial`` string."""
    otu = next(initialized_repo.iter_otus())

    isolate = otu.isolates[0]

    assert initialized_repo.get_isolate_id_by_partial(str(isolate.id)[:8]) == isolate.id


class TestCreateOTUWithValidation:
    def test_ok(self, initialized_repo: Repo):
        """Test that a complete OTU can be created, validated and retrieved."""
        assert next(iter(initialized_repo.iter_otus())) is not None

    def test_bad_otu_data_fail(self, empty_repo: Repo, initialized_repo: Repo):
        """Test that a bad OTU with bad data fails at validation stage."""
        base_otu = next(iter(initialized_repo.iter_otus()))

        assert isinstance(base_otu, OTUBuilder)

        base_otu.taxid = -3000

        with capture_logs() as logs:
            with empty_repo.lock(), empty_repo.use_transaction():
                init_otu_with_contents(
                    empty_repo,
                    base_otu,
                )

        assert any(
            [
                log["event"] == "ValidationError" and log["log_level"] == "warning"
                for log in logs
            ]
        )


class TestCreateIsolateValidation:
    """Test the validation of new added isolates."""

    def test_name_exists_fail(self, initialized_repo: Repo):
        """Test that an isolate is not created if the name is already taken."""
        otu_init = next(iter(initialized_repo.iter_otus()))

        last_id_before_test = initialized_repo.last_id

        with initialized_repo.lock(), initialized_repo.use_transaction():
            sequence_1 = initialized_repo.create_sequence(
                otu_init.id,
                "MK000001.1",
                "TMV",
                otu_init.plan.segments[0].id,
                "ACGTACGTACGTACG",
            )

            assert sequence_1

            isolate_b = initialized_repo.create_isolate(
                otu_init.id,
                IsolateName(IsolateNameType.ISOLATE, "A"),
            )

            assert isolate_b

            initialized_repo.link_sequence(otu_init.id, isolate_b.id, sequence_1.id)

            assert initialized_repo.last_id > last_id_before_test

        otu_after = initialized_repo.get_otu(otu_init.id)

        assert otu_after
        assert isolate_b.id not in otu_after.isolate_ids

    @pytest.mark.parametrize("bad_sequence", ["ACGTACGTA", "ACGTACGTACGTACGTTTTTT"])
    def test_bad_segment_length_fail(self, initialized_repo: Repo, bad_sequence: str):
        """Test that a new isolate with segments that are too long or short does not pass validation."""
        otu_init = next(iter(initialized_repo.iter_otus()))

        last_id_before_test = initialized_repo.last_id

        with initialized_repo.lock(), initialized_repo.use_transaction():
            sequence_1 = initialized_repo.create_sequence(
                otu_init.id,
                "MK000001.1",
                "TMV",
                otu_init.plan.segments[0].id,
                bad_sequence,
            )

            assert sequence_1

            isolate_b = initialized_repo.create_isolate(
                otu_init.id,
                IsolateName(IsolateNameType.ISOLATE, "B"),
            )

            assert isolate_b

            initialized_repo.link_sequence(otu_init.id, isolate_b.id, sequence_1.id)

            assert initialized_repo.last_id > last_id_before_test

        otu_after = initialized_repo.get_otu(otu_init.id)

        assert otu_after
        assert isolate_b.id not in otu_after.isolate_ids


class TestExcludeAccessions:
    """Test that excluded accessions are reflected in events.

    Excluded accessions are exempt from future queries to the NCBI Nucleotide.
    """

    def test_ok(self, initialized_repo: Repo):
        """Test that Repo.exclude_accessions() writes the correct event."""
        otu = initialized_repo.get_otu_by_taxid(12242)
        id_at_creation = initialized_repo.last_id

        assert otu
        assert otu.excluded_accessions == set()
        assert id_at_creation == 5

        accessions = {"TM100021", "TM100022", "TM100023"}

        with initialized_repo.lock(), initialized_repo.use_transaction():
            initialized_repo.exclude_accessions(otu.id, accessions)

        assert initialized_repo.get_otu(otu.id).excluded_accessions == accessions
        assert initialized_repo.last_id == id_at_creation + 1 == 6

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
            "id": 6,
            "query": {
                "otu_id": str(otu.id),
            },
            "type": "UpdateExcludedAccessions",
        }

        with initialized_repo.lock(), initialized_repo.use_transaction():
            initialized_repo.exclude_accessions(otu.id, {"TM100024"})

        assert initialized_repo.last_id == id_at_creation + 2 == 7
        assert initialized_repo.get_otu(otu.id).excluded_accessions == accessions | {
            "TM100024"
        }

    def test_existing_accession(self, initialized_repo: Repo):
        """Test that the accession of a sequence already contained in the OTU cannot
        be added to the excluded accessions.
        """
        assert (id_before_exclude := initialized_repo.last_id) == 5

        otu = initialized_repo.get_otu_by_taxid(12242)

        assert otu.excluded_accessions == set()

        accession = next(iter(otu.accessions))

        with initialized_repo.lock(), initialized_repo.use_transaction():
            initialized_repo.exclude_accessions(otu.id, {f"{accession}.1"})

        assert initialized_repo.last_id == id_before_exclude == 5

        assert initialized_repo.get_otu_by_taxid(12242).excluded_accessions == set()

        with open(
            initialized_repo.path / "src" / f"{initialized_repo.last_id:08}.json"
        ) as f:
            event = orjson.loads(f.read())

        assert event["type"] != "UpdateExcludedAccessions"

    def test_already_excluded(self, initialized_repo: Repo):
        """Test that an attempted redundant exclusion does not create a new event."""
        repo = initialized_repo

        otu_before = initialized_repo.get_otu_by_taxid(12242)

        assert otu_before
        assert otu_before == repo.get_otu(otu_before.id)
        assert otu_before.excluded_accessions == set()

        event_id_after_creation = repo.last_id

        accessions = {"TM100021", "TM100022", "TM100023"}

        with repo.lock(), repo.use_transaction():
            repo.exclude_accessions(otu_before.id, accessions)

        # Event ID should increment for the first successful exclusion.
        assert repo.last_id == event_id_after_creation + 1

        with repo.lock(), repo.use_transaction():
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
        otu_before = initialized_repo.get_otu_by_taxid(12242)

        assert otu_before

        event_id_after_creation = initialized_repo.last_id

        accessions = {"TM100021", "TM100022", "TM100023"}

        with initialized_repo.lock(), initialized_repo.use_transaction():
            initialized_repo.exclude_accessions(otu_before.id, accessions)

        otu_after_first = initialized_repo.get_otu(otu_before.id)

        assert otu_after_first
        assert otu_after_first.excluded_accessions == accessions
        assert initialized_repo.last_id == event_id_after_creation + 1

        with initialized_repo.lock(), initialized_repo.use_transaction():
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
        creates the expected OTUBuilder.excluded_accessions set.
        """
        target_repo = initialized_repo

        otu = target_repo.get_otu_by_taxid(12242)

        assert otu
        assert (id_at_creation := target_repo.last_id) == 5

        accessions = {"TM100021", "TM100022", "TM100023"}

        with target_repo.lock():
            with target_repo.use_transaction():
                target_repo.exclude_accessions(otu.id, accessions)

            assert target_repo.last_id == id_at_creation + 1 == 6

            assert target_repo.get_otu(otu.id).excluded_accessions == accessions

            with target_repo.use_transaction():
                target_repo.allow_accessions(otu.id, ["TM100021", "TM100022"])

            assert target_repo.last_id == id_at_creation + 2 == 7

        otu_after = target_repo.get_otu(otu.id)

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
                "id": 7,
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
        target_repo = initialized_repo

        otu = target_repo.get_otu_by_taxid(12242)

        accessions = {"TM100021", "TM100022", "TM100023"}

        with target_repo.lock(), target_repo.use_transaction():
            target_repo.exclude_accessions(otu.id, accessions)

        assert (id_after_first_exclusion := target_repo.last_id) == 6

        assert target_repo.get_otu(otu.id).excluded_accessions == accessions

        # Attempt to allow an accession not on the exclusion list
        with target_repo.lock(), target_repo.use_transaction():
            target_repo.allow_accessions(otu.id, ["TM100024"])

        assert target_repo.last_id == id_after_first_exclusion == 6
        assert not target_repo.path.joinpath("src", "00000007.json").exists()

        assert target_repo.get_otu(otu.id).excluded_accessions == accessions

        # Clear the excluded accessions list
        with target_repo.lock(), target_repo.use_transaction():
            target_repo.allow_accessions(otu.id, accessions)

        assert target_repo.get_otu(otu.id).excluded_accessions == set()

        assert target_repo.last_id == id_after_first_exclusion + 1 == 7
        assert target_repo.path.joinpath("src", "00000007.json").exists()


class TestDeleteIsolate:
    """Test that an isolate can be redacted from an OTU."""

    def test_ok(self, initialized_repo: Repo):
        """Test basic functionality."""
        otu_before = initialized_repo.get_otu_by_taxid(12242)

        assert otu_before

        with initialized_repo.lock(), initialized_repo.use_transaction():
            sequence_2 = initialized_repo.create_sequence(
                otu_before.id,
                "TN000001.1",
                "TMV",
                otu_before.plan.segments[0].id,
                "TTACGTGGAGAGACC",
            )

            isolate_b = initialized_repo.create_isolate(
                otu_before.id,
                IsolateName(IsolateNameType.ISOLATE, "B"),
            )

            assert isolate_b
            assert sequence_2

            initialized_repo.link_sequence(
                otu_id=otu_before.id,
                isolate_id=isolate_b.id,
                sequence_id=sequence_2.id,
            )

        event_id_before_delete = initialized_repo.last_id

        otu_before = initialized_repo.get_otu(otu_before.id)

        assert otu_before
        assert isolate_b.id in otu_before.isolate_ids
        assert "TN000001" in otu_before.accessions
        assert "TN000001" in otu_before.get_isolate(isolate_b.id).accessions
        assert event_id_before_delete == initialized_repo.last_id == 8

        with initialized_repo.lock(), initialized_repo.use_transaction():
            initialized_repo.delete_isolate(
                otu_before.id,
                isolate_b.id,
                rationale="Testing redaction",
            )

        otu_after = initialized_repo.get_otu(otu_before.id)

        assert otu_after
        assert otu_after.get_isolate(isolate_b.id) is None
        assert otu_before != otu_after
        assert len(otu_after.isolates) == len(otu_before.isolates) - 1
        assert isolate_b.id not in otu_after.isolate_ids
        assert isolate_b.accessions not in otu_after.accessions

        assert initialized_repo.last_id == event_id_before_delete + 1


def test_replace_sequence(initialized_repo: Repo):
    """Test the replacement of an existing sequence using a new Genbank accession and
    record.
    """
    otu_init = initialized_repo.get_otu_by_taxid(12242)

    replaceable_sequence = otu_init.get_sequence_by_accession("TM000001")

    isolate_id = next(
        iter(otu_init.get_isolate_ids_containing_sequence_id(replaceable_sequence.id))
    )

    with initialized_repo.lock():
        with initialized_repo.use_transaction():
            new_sequence = initialized_repo.create_sequence(
                otu_init.id,
                "RP000001.1",
                "TMV edit",
                otu_init.plan.segments[0].id,
                "TACGTGGAGAGACCA",
            )

        with initialized_repo.use_transaction():
            initialized_repo.replace_sequence(
                otu_id=otu_init.id,
                isolate_id=isolate_id,
                sequence_id=new_sequence.id,
                replaced_sequence_id=replaceable_sequence.id,
                rationale="Testing sequence redaction",
            )

        otu_after = initialized_repo.get_otu(otu_init.id)

    assert otu_after.get_sequence_by_id(new_sequence.id) is not None
    assert otu_after.get_sequence_by_id(replaceable_sequence.id) is None


class TestMalformedEvent:
    """Test that malformed events cannot be rehydrated."""

    def test_bad_event_typing(self, initialized_repo: Repo):
        """Test that an event with an invalid event type discriminator does not attempt
        to rehydrate.
        """
        filepath = initialized_repo.path.joinpath("src", "00000002.json")

        with open(filepath, "rb") as f:
            event = orjson.loads(f.read())

        otu = initialized_repo.get_otu_by_taxid(12242)

        assert type(otu) is OTUBuilder

        event["type"] = "MalformedEvent"

        with open(filepath, "wb") as f:
            f.write(orjson.dumps(event))

        with pytest.raises(ValueError, match="Unknown event type: MalformedEvent"):
            initialized_repo.get_otu_by_taxid(12242)

    def test_bad_event_data(self, initialized_repo: Repo):
        """Test that an event with bad data cannot be rehydrated."""
        path = initialized_repo.path.joinpath("src", "00000002.json")

        with open(path, "rb") as f:
            event = orjson.loads(f.read())

        with initialized_repo.lock():
            otu = initialized_repo.get_otu_by_taxid(12242)

        assert type(otu) is OTUBuilder

        event["data"]["taxid"] = "popcorn"

        with open(path, "wb") as f:
            f.write(orjson.dumps(event))

        with (
            initialized_repo.lock(),
            pytest.raises(
                ValueError,
                match="Input should be a valid integer, unable to parse string",
            ),
        ):
            initialized_repo.get_otu_by_taxid(12242)


def test_prune_on_load(initialized_repo: Repo):
    """Test that the Repo properly prunes the events and the index even before read-only actions,
    such as .get_otu().
    """
    assert initialized_repo.get_otu_by_taxid(12242)

    assert (
        initialized_repo.path / "src" / f"0000000{initialized_repo.head_id}.json"
    ).exists()

    assert initialized_repo.last_id == initialized_repo.head_id

    head_path = initialized_repo.path / "head"

    # Set head file to before OTU addition and reload Repo.
    head_path.unlink()
    with open(head_path, "w") as f:
        f.write("1")

    reinitialized_repo = Repo(initialized_repo.path)

    assert reinitialized_repo.head_id < initialized_repo.head_id

    assert not (
        reinitialized_repo.path / "src" / f"0000000{initialized_repo.head_id}.json"
    ).exists()

    assert reinitialized_repo.last_id == reinitialized_repo.head_id

    assert reinitialized_repo.get_otu_by_taxid(12242) is None
