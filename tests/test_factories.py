"""Tests for data factories used in testing."""

import uuid
from uuid import uuid4

from syrupy import SnapshotAssertion

from ref_builder.models.accession import Accession
from ref_builder.ncbi.models import NCBISource
from ref_builder.otu.builders.sequence import SequenceBuilder
from ref_builder.otu.validators.isolate import Isolate
from ref_builder.otu.validators.otu import OTU
from ref_builder.otu.validators.sequence import Sequence
from tests.fixtures.factories import (
    IsolateFactory,
    NCBIGenbankFactory,
    NCBISourceFactory,
    OTUFactory,
    SequenceFactory,
)


def test_ncbi_source_factory(ncbi_source_factory: NCBISourceFactory):
    """Test that NCBISourceFactory creates valid mock NCBISource objects."""
    assert all(
        isinstance(dummy_source, NCBISource)
        for dummy_source in ncbi_source_factory.coverage()
    )


def test_ncbi_genbank_factory(
    ncbi_genbank_factory: NCBIGenbankFactory, snapshot: SnapshotAssertion
):
    """Test that NCBIGenbankFactory creates valid fake Genbank records."""
    records = list(ncbi_genbank_factory.coverage())

    assert [
        SequenceBuilder(
            id=uuid4(),
            accession=Accession.from_string(record.accession_version),
            definition=record.definition,
            segment=uuid.uuid4(),
            sequence=record.sequence,
        )
        for record in records
    ]

    assert [record.model_dump() for record in records] == snapshot


class TestNCBIGenbankFactoryBuildIsolate:
    """Tests for NCBIGenbankFactory.build_isolate()."""

    def test_sequential_accessions(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that build_isolate() generates sequential accessions."""
        records = ncbi_genbank_factory.build_isolate(segment_count=3, refseq=True)

        accessions = [record.accession for record in records]

        # Extract numeric parts
        numbers = [int(acc.split("_")[1]) for acc in accessions]

        # Verify they are sequential
        assert numbers == [numbers[0], numbers[0] + 1, numbers[0] + 2]

    def test_consistent_version(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that all records have the same version."""
        records = ncbi_genbank_factory.build_isolate(segment_count=3, version=2)

        versions = [record.accession_version.split(".")[1] for record in records]

        assert all(v == "2" for v in versions)

    def test_letter_segment_keys(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that letter segment keys are generated correctly."""
        records = ncbi_genbank_factory.build_isolate(
            segment_count=3,
            segment_key_style="letter",
            segment_prefix="RNA",
            segment_delimiter="-",
        )

        segment_names = [record.source.segment for record in records]

        assert segment_names == ["RNA-A", "RNA-B", "RNA-C"]

    def test_number_segment_keys(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that number segment keys are generated correctly."""
        records = ncbi_genbank_factory.build_isolate(
            segment_count=3,
            segment_key_style="number",
            segment_prefix="DNA",
            segment_delimiter=" ",
        )

        segment_names = [record.source.segment for record in records]

        assert segment_names == ["DNA 1", "DNA 2", "DNA 3"]

    def test_shared_metadata(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that all records share the same organism, taxid, etc."""
        records = ncbi_genbank_factory.build_isolate(segment_count=3)

        # Check organism is consistent
        organisms = [record.organism for record in records]
        assert len(set(organisms)) == 1

        # Check taxid is consistent
        taxids = [record.source.taxid for record in records]
        assert len(set(taxids)) == 1

        # Check mol_type is consistent
        moltypes = [record.moltype for record in records]
        assert len(set(moltypes)) == 1

        # Check host is consistent
        hosts = [record.source.host for record in records]
        assert len(set(hosts)) == 1

    def test_refseq_accessions(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that RefSeq accessions are generated when refseq=True."""
        records = ncbi_genbank_factory.build_isolate(segment_count=3, refseq=True)

        assert all(record.accession.startswith("NC_") for record in records)
        assert all(record.refseq for record in records)

    def test_genbank_accessions(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that GenBank accessions are generated when refseq=False."""
        records = ncbi_genbank_factory.build_isolate(segment_count=3, refseq=False)

        assert all(not record.accession.startswith("NC_") for record in records)
        assert all(not record.refseq for record in records)

    def test_segment_count_validation(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that ValueError is raised when segment_count > 26 for letter style."""
        import pytest

        with pytest.raises(ValueError, match="segment_count cannot exceed 26"):
            ncbi_genbank_factory.build_isolate(
                segment_count=27, segment_key_style="letter"
            )

    def test_number_style_large_count(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that number style supports large segment counts."""
        records = ncbi_genbank_factory.build_isolate(
            segment_count=30,
            segment_key_style="number",
            segment_delimiter=" ",
        )

        assert len(records) == 30

        # Verify keys go up to 30
        segment_keys = [record.source.segment.split()[-1] for record in records]
        assert segment_keys[-1] == "30"
        assert segment_keys[0] == "1"

    def test_auto_prefix_from_moltype(self, ncbi_genbank_factory: NCBIGenbankFactory):
        """Test that segment prefix is auto-detected from mol_type."""
        from ref_builder.ncbi.models import NCBISourceMolType

        # Test DNA
        base_source = NCBISourceFactory.build(mol_type=NCBISourceMolType.GENOMIC_DNA)
        records = ncbi_genbank_factory.build_isolate(
            segment_count=2,
            base_source=base_source,
            segment_delimiter="-",
        )

        assert all(record.source.segment.startswith("DNA-") for record in records)

        # Test RNA
        base_source = NCBISourceFactory.build(mol_type=NCBISourceMolType.GENOMIC_RNA)
        records = ncbi_genbank_factory.build_isolate(
            segment_count=2,
            base_source=base_source,
            segment_delimiter="-",
        )

        assert all(record.source.segment.startswith("RNA-") for record in records)

    def test_snapshot(
        self,
        ncbi_genbank_factory: NCBIGenbankFactory,
        snapshot: SnapshotAssertion,
    ):
        """Test that build_isolate() generates expected data structure."""
        records = ncbi_genbank_factory.build_isolate(
            segment_count=3,
            refseq=True,
            segment_prefix="RNA",
            segment_key_style="letter",
            segment_delimiter="-",
            version=1,
        )

        assert [record.model_dump() for record in records] == snapshot


def test_sequence_factory(sequence_factory: SequenceFactory):
    """Test that SequenceFactory creates valid mock sequence data."""
    assert all(
        Sequence.model_validate(sequence.model_dump())
        for sequence in sequence_factory.coverage()
    )


def test_isolate_factory():
    """Test that IsolateFactory creates valid mock isolate data."""
    assert all(
        Isolate.model_validate(isolate.model_dump())
        for isolate in IsolateFactory.coverage()
    )


def test_otu_factory(otu_factory: OTUFactory):
    """Test that OTUFactory creates valid mock OTU data."""
    mock_otus = otu_factory.coverage()
    assert all(OTU.model_validate(otu.model_dump()) for otu in mock_otus)
