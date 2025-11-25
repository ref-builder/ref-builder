"""Build FASTA and CSV files from a reference repository."""

import csv
from pathlib import Path

from ref_builder.repo import Repo


def build_fasta(
    output_path: Path,
    path: Path,
    include_csv: bool,
) -> tuple[Path, Path | None]:
    """Build a FASTA file from a reference repository.

    :param output_path: Path for the output FASTA file
    :param path: Path to the reference repository
    :param include_csv: Whether to generate accompanying CSV file
    :return: Tuple of (fasta_path, csv_path or None)
    """
    repo = Repo(path)

    fasta_lines: list[str] = []
    csv_rows: list[dict[str, str]] = []

    for otu in repo.iter_otus():
        for isolate in otu.isolates:
            for sequence in isolate.sequences:
                fasta_lines.append(f">{sequence.accession}")
                fasta_lines.append(sequence.sequence)

                if include_csv:
                    segment = otu.plan.get_segment_by_id(sequence.segment)
                    segment_name = (
                        str(segment.name) if segment and segment.name else "Unnamed"
                    )

                    csv_rows.append(
                        {
                            "accession": str(sequence.accession),
                            "otu_name": otu.name,
                            "otu_taxid": str(otu.taxid),
                            "isolate_name": (str(isolate.name) if isolate.name else ""),
                            "segment_name": segment_name,
                            "definition": sequence.definition,
                            "otu_id": str(otu.id),
                            "isolate_id": str(isolate.id),
                        }
                    )

    with open(output_path, "w") as f:
        f.write("\n".join(fasta_lines))
        if fasta_lines:
            f.write("\n")

    csv_path = None

    if include_csv and csv_rows:
        csv_path = output_path.with_suffix(".csv")
        fieldnames = [
            "accession",
            "otu_name",
            "otu_taxid",
            "isolate_name",
            "segment_name",
            "definition",
            "otu_id",
            "isolate_id",
        ]

        with open(csv_path, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(csv_rows)

    return output_path, csv_path
