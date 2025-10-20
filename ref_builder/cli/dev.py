"""Development commands for ref-builder maintainers."""

import json
from pathlib import Path

import click
from rich.table import Table

from ref_builder.console import console
from ref_builder.ncbi.client import NCBIClient
from tests.fixtures.ncbi import OTUManifest, OTUSpec


@click.group(name="dev")
def dev() -> None:
    """Development commands for ref-builder maintainers."""


@dev.command(name="list")
def list_otus() -> None:
    """List all mock OTUs with their taxids and names."""
    table = Table(title="Mock NCBI Data")

    table.add_column("Taxid", style="cyan", justify="right")
    table.add_column("Name", style="green")
    table.add_column("# Segments", justify="right")
    table.add_column("# Isolates", justify="right")

    data_dir = Path("tests/fixtures/ncbi/otus")
    otu_data = []

    for attr_name in dir(OTUManifest):
        if attr_name.startswith("_"):
            continue
        attr = getattr(OTUManifest, attr_name)
        if not isinstance(attr, OTUSpec):
            continue

        json_path = data_dir / f"{attr_name}.json"
        if not json_path.exists():
            continue

        data = json.loads(json_path.read_text())
        taxid = data["taxonomy"]["id"]
        name = data["taxonomy"]["name"]
        segment_count = len(attr.refseq)
        isolate_count = len(attr.isolates) + 1

        otu_data.append((taxid, name, segment_count, isolate_count))

    for taxid, name, segment_count, isolate_count in sorted(otu_data):
        table.add_row(
            str(taxid),
            name,
            str(segment_count),
            str(isolate_count),
        )

    console.print(table)
    console.print(f"\n[bold]Total:[/bold] {len(otu_data)} OTUs")


@dev.command(name="refresh")
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default="tests/fixtures/ncbi/otus",
    help="Output directory for generated JSON/FASTA files",
)
def refresh(output_dir: Path) -> None:
    """Regenerate mock NCBI data from manifest by fetching from NCBI.

    Reads OTUManifest and fetches real data from NCBI for each OTU.
    Generates JSON (taxonomy + genbank) and FASTA files.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Initialize NCBI client (ignore cache to get fresh data)
    client = NCBIClient(ignore_cache=True)

    generated_count = 0

    # Iterate over all OTUSpec attributes in the manifest
    for attr_name in dir(OTUManifest):
        if attr_name.startswith("_"):
            continue

        attr = getattr(OTUManifest, attr_name)
        if not isinstance(attr, OTUSpec):
            continue

        click.echo(f"Fetching data for {attr_name}...")

        # Fetch all GenBank records
        all_accessions = attr.all_accessions
        genbank_records = client.fetch_genbank_records(all_accessions)

        if not genbank_records:
            click.echo(f"  Warning: No GenBank records found for {attr_name}", err=True)
            continue

        # Build genbank dict keyed by accession
        genbank_dict = {
            record.accession: record.model_dump() for record in genbank_records
        }

        # Fetch taxonomy from first refseq accession
        first_refseq_record = next(
            (r for r in genbank_records if r.accession in attr.refseq), None
        )

        if not first_refseq_record:
            click.echo(f"  Warning: No refseq record found for {attr_name}", err=True)
            continue

        taxid = first_refseq_record.source.taxid
        taxonomy = client.fetch_taxonomy_record(taxid)

        if not taxonomy:
            click.echo(f"  Warning: No taxonomy found for taxid {taxid}", err=True)
            continue

        # Write JSON file
        json_data = {
            "taxonomy": taxonomy.model_dump(),
            "genbank": genbank_dict,
        }

        json_path = output_path / f"{attr_name}.json"
        json_path.write_text(json.dumps(json_data, indent=2))
        click.echo(f"  Wrote {json_path}")

        generated_count += 1

    # Generate type stub
    stub_path = Path("tests/fixtures/ncbi/otus.pyi")
    stub_lines = [
        "from tests.fixtures.ncbi.models import OTUHandle",
        "",
        "",
        "class OTURegistry:",
    ]

    for attr_name in dir(OTUManifest):
        if attr_name.startswith("_"):
            continue
        attr = getattr(OTUManifest, attr_name)
        if isinstance(attr, OTUSpec):
            stub_lines.append("    @property")
            stub_lines.append(f"    def {attr_name}(self) -> OTUHandle: ...")
            stub_lines.append("")

    stub_path.write_text("\n".join(stub_lines))
    click.echo(f"  Wrote {stub_path}")

    click.echo(f"\n✓ Generated {generated_count} OTUs in {output_path}")
