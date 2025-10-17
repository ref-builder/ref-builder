"""Development commands for ref-builder maintainers."""

from pathlib import Path

import click
from jinja2 import Environment, FileSystemLoader

from ref_builder.dev.format_model import format_model, format_value_list
from tests.fixtures.ncbi import mock_ncbi_client


@click.group(name="dev")
def dev() -> None:
    """Development commands for ref-builder maintainers."""


@dev.command(name="refresh")
@click.option(
    "--output-dir",
    type=click.Path(path_type=Path),
    default="tests/fixtures/ncbi",
    help="Output directory for generated modules",
)
def refresh(output_dir: Path) -> None:
    """Regenerate mock NCBI data modules from current registry.

    Reads the current MockNCBIClient registry and regenerates all per-OTU
    module files in tests/fixtures/ncbi/.

    This is useful when:
    - Adding new mock data to existing OTUs
    - Updating mock data structure
    - Fixing formatting issues in generated modules
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Load Jinja2 template
    template_dir = Path(__file__).parent.parent / "dev"
    env = Environment(loader=FileSystemLoader(template_dir))
    template = env.get_template("otu_module.py.j2")

    generated_count = 0

    for taxid, refseq, isolate_groups in mock_ncbi_client.get_otu_structure():
        # Get taxonomy record
        taxonomy = mock_ncbi_client._taxonomy_records.get(taxid)

        if taxonomy is None:
            click.echo(
                f"Warning: No taxonomy record for taxid {taxid}, skipping",
                err=True,
            )
            continue

        name = taxonomy.name

        # Create safe variable name from organism name
        var_name = name.lower().replace(" ", "_").replace("-", "_")

        # Get GenBank records for RefSeq accessions
        refseq_genbank = []
        for acc in refseq:
            if acc in mock_ncbi_client._genbank_records:
                refseq_genbank.append(mock_ncbi_client._genbank_records[acc])

        # Build isolate data
        isolates = []
        for accessions in isolate_groups:
            isolate_genbank = []
            for acc in accessions:
                if acc in mock_ncbi_client._genbank_records:
                    isolate_genbank.append(mock_ncbi_client._genbank_records[acc])

            isolates.append(
                {
                    "accessions": repr(accessions),
                    "genbank": format_value_list(isolate_genbank, indent=3)
                    if isolate_genbank
                    else "[]",
                }
            )

        # Render template
        content = template.render(
            taxid=taxid,
            name=name,
            var_name=var_name,
            refseq=repr(refseq),
            taxonomy=format_model(taxonomy, indent=2) if taxonomy else None,
            genbank=format_value_list(refseq_genbank, indent=2)
            if refseq_genbank
            else "[]",
            isolates=isolates,
        )

        # Write module file
        module_path = output_path / f"otu_{var_name}.py"
        module_path.write_text(content)

        generated_count += 1
        click.echo(f"Generated {module_path}")

    click.echo(f"\n✓ Generated {generated_count} OTU modules in {output_path}")
