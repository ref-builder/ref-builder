"""``ref-builder validate`` — retroactive integrity checks for a repository."""

import sys
from pathlib import Path

import click
from rich.table import Table

from ref_builder.audit import audit_accessions
from ref_builder.cli.options import path_option
from ref_builder.console import console
from ref_builder.repo import Repo


@click.group(name="validate")
def validate() -> None:
    """Run integrity checks against an existing repository."""


@validate.command(name="accessions")
@path_option
def validate_accessions(path: Path) -> None:
    """Report duplicate accessions across and within OTUs.

    Exits non-zero when duplicates are found so the command can gate CI.
    """
    repo = Repo(path)

    report = audit_accessions(repo)

    if not report.has_findings:
        console.print("[green]No duplicate accessions found.[/green]")
        return

    if report.cross_otu:
        table = Table(
            title="Cross-OTU duplicate accessions",
            title_style="bold red",
            show_lines=True,
        )
        table.add_column("Accession key")
        table.add_column("OTU (taxid)")
        table.add_column("Isolate")

        for accession_key in sorted(report.cross_otu):
            placements = report.cross_otu[accession_key]
            for index, placement in enumerate(placements):
                table.add_row(
                    accession_key if index == 0 else "",
                    f"{placement.otu_name} ({placement.otu_taxid})",
                    f"{placement.isolate_name or 'Unnamed'} [{placement.isolate_id}]",
                )

        console.print(table)

    if report.within_otu:
        table = Table(
            title="Within-OTU duplicate accessions",
            title_style="bold red",
            show_lines=True,
        )
        table.add_column("Accession key")
        table.add_column("OTU (taxid)")
        table.add_column("Isolate")

        for accession_key in sorted(report.within_otu):
            placements = report.within_otu[accession_key]
            for index, placement in enumerate(placements):
                table.add_row(
                    accession_key if index == 0 else "",
                    f"{placement.otu_name} ({placement.otu_taxid})",
                    f"{placement.isolate_name or 'Unnamed'} [{placement.isolate_id}]",
                )

        console.print(table)

    sys.exit(1)
