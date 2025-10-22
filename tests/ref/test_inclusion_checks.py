import pytest
from faker import Faker
from pydantic import ValidationError

from ref_builder.ncbi.client import NCBIClient
from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.utils import IsolateName, IsolateNameType
from ref_builder.repo import Repo
from ref_builder.services.isolate import IsolateService
from tests.fixtures.factories import NCBIGenbankFactory, NCBISourceFactory
from tests.fixtures.providers import AccessionProvider, SequenceProvider

faker = Faker()
faker.add_provider(AccessionProvider)
faker.add_provider(SequenceProvider)


class TestAddMultipartiteIsolate:
    """Test that new multipartite isolates follow plan length limits."""

    @pytest.fixture(autouse=True)
    def _setup(
        self,
        ncbi_genbank_factory: NCBIGenbankFactory,
        ncbi_source_factory: NCBISourceFactory,
    ):
        def func(
            otu: OTUBuilder,
            sequence_length_multiplier: float = 1.0,
        ) -> list[NCBIGenbank]:
            """Generate a collection of mock Genbank records capable of passing the isolate inclusion checks
            on a given multipartite OTU, as long as sequence_length_multiplier is 1.0.
            Builds one mock record per segment.

            :param otu: A monopartite OTU that the generated record must match.
            :param sequence_length_multiplier: A float multiplier for the generated mock sequence length.
            """
            records = []

            accession_starter = 100000

            for i in range(len(otu.plan.required_segments)):
                segment = otu.plan.required_segments[i]

                source = ncbi_source_factory.build(
                    taxid=otu.taxid,
                    organism=otu.name,
                    segment=str(segment.name),
                )

                sequence_length = int(segment.length * sequence_length_multiplier)

                records.append(
                    ncbi_genbank_factory.build(
                        accession=f"AB{accession_starter + i}",
                        sequence=faker.sequence(
                            min=sequence_length, max=sequence_length
                        ),
                        source=source,
                    )
                )

            return records

        self.create_mock_isolate_records = func

    @pytest.mark.parametrize("sequence_length_multiplier", [1.0, 1.03, 0.98])
    def test_ok(self, scratch_repo: Repo, sequence_length_multiplier: float):
        """Test that sequences within recommended length tolerance
        are added without issue.
        """
        otu_before = scratch_repo.get_otu_by_taxid(438782)

        records = self.create_mock_isolate_records(
            otu_before,
            sequence_length_multiplier,
        )

        with scratch_repo.lock():
            isolate_service = IsolateService(scratch_repo, NCBIClient(False))
            isolate = isolate_service.create_from_records(
                otu_before.id,
                IsolateName(type=IsolateNameType.ISOLATE, value="mock"),
                records,
            )

        otu_after = scratch_repo.get_otu_by_taxid(438782)

        assert isolate.id in otu_after.isolate_ids
        assert isolate.accessions.issubset(otu_after.accessions)

    @pytest.mark.parametrize("sequence_length_multiplier", [0.5, 1.035, 20.0])
    def test_fail(self, scratch_repo: Repo, sequence_length_multiplier: float):
        """Test that sequences that exceed recommended length tolerance are
        automatically rejected.
        """
        otu_before = scratch_repo.get_otu_by_taxid(438782)

        records = self.create_mock_isolate_records(
            otu_before,
            sequence_length_multiplier,
        )

        with scratch_repo.lock():
            isolate_service = IsolateService(scratch_repo, NCBIClient(False))
            try:
                isolate_service.create_from_records(
                    otu_before.id,
                    IsolateName(type=IsolateNameType.ISOLATE, value="mock"),
                    records,
                )
            except ValidationError as exc:
                for error in exc.errors():
                    assert error["type"] in ("sequence_too_short", "sequence_too_long")

        otu_after = scratch_repo.get_otu(otu_before.id)

        assert otu_after.isolate_ids == otu_before.isolate_ids
        assert otu_after.accessions == otu_before.accessions
