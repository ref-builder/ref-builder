"""Tests for FASTA and CSV build functionality."""

import csv
from pathlib import Path

from ref_builder.build.fasta import build_fasta
from ref_builder.repo import Repo


class TestBuildFasta:
    """Test the build_fasta function."""

    def test_ok(self, scratch_repo: Repo, tmp_path: Path):
        """Test that build_fasta produces valid FASTA and CSV files."""
        output_path = tmp_path / "output.fasta"

        fasta_path, csv_path = build_fasta(output_path, scratch_repo.path, True)

        assert fasta_path.exists()
        assert csv_path is not None
        assert csv_path.exists()
        assert csv_path == tmp_path / "output.csv"

        fasta_content = fasta_path.read_text()
        lines = fasta_content.strip().split("\n")

        headers = [line for line in lines if line.startswith(">")]
        sequences = [line for line in lines if not line.startswith(">")]

        assert len(headers) > 0
        assert len(headers) == len(sequences)

        for header in headers:
            accession = header[1:]
            assert "." in accession

        with open(csv_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)

        assert len(rows) == len(headers)

        expected_columns = {
            "accession",
            "otu_name",
            "otu_taxid",
            "isolate_name",
            "segment_name",
            "definition",
            "otu_id",
            "isolate_id",
        }
        assert set(reader.fieldnames) == expected_columns

    def test_no_csv(self, scratch_repo: Repo, tmp_path: Path):
        """Test that build_fasta can skip CSV generation."""
        output_path = tmp_path / "output.fasta"

        fasta_path, csv_path = build_fasta(output_path, scratch_repo.path, False)

        assert fasta_path.exists()
        assert csv_path is None
        assert not (tmp_path / "output.csv").exists()

    def test_accession_matches_repo(self, scratch_repo: Repo, tmp_path: Path):
        """Test that FASTA headers match accessions from repo."""
        output_path = tmp_path / "output.fasta"

        fasta_path, _ = build_fasta(output_path, scratch_repo.path, True)

        fasta_content = fasta_path.read_text()

        repo = Repo(scratch_repo.path)
        expected_accessions = set()

        for otu in repo.iter_otus():
            for isolate in otu.isolates:
                for sequence in isolate.sequences:
                    expected_accessions.add(str(sequence.accession))

        for accession in expected_accessions:
            assert f">{accession}" in fasta_content
