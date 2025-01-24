from collections.abc import Iterator

import rich.console
from rich.table import Table

from ref_builder.models import OTUMinimal
from ref_builder.plan import SegmentRule
from ref_builder.resources import RepoOTU


def _render_taxonomy_id_link(taxid: int) -> str:
    return f"[link=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxid}]{taxid}[/link]"


def _render_nucleotide_link(accession: str) -> str:
    return f"[link=https://www.ncbi.nlm.nih.gov/nuccore/{accession}]{accession}[/link]"


def print_otu(otu: RepoOTU) -> None:
    """Print the details for an OTU to the console.

    :param otu: The OTU to print.

    """
    console.print(f"[bold][underline]{otu.name}[/bold][/underline]")
    console.line()

    table = Table(
        box=None,
        show_header=False,
    )

    table.add_row("[bold]ACRONYM[/bold]", otu.acronym)
    table.add_row("[bold]ID[/bold]", str(otu.id))

    if otu.legacy_id:
        table.add_row("[bold]LEGACY ID[/bold]", otu.legacy_id)

    table.add_row("[bold]TAXID[/bold]", _render_taxonomy_id_link(otu.taxid))

    max_accession_length = max(
        len(str(sequence.accession))
        for isolate in otu.isolates
        for sequence in isolate.sequences
    )

    max_segment_name_length = max(
        len(str(segment.name)) for segment in otu.plan.segments
    )

    console.print(table)

    console.line()
    console.print("[bold]PLAN[/bold]")
    console.line()

    plan_table = Table(box=None)

    plan_table.add_column("NAME")
    plan_table.add_column("REQUIRED")
    plan_table.add_column("LENGTH")
    plan_table.add_column("TOLERANCE")
    plan_table.add_column("ID")

    if not otu.plan.monopartite:
        for segment in otu.plan.segments:
            plan_table.add_row(
                str(segment.name),
                "[red]Yes[/red]"
                if segment.required == SegmentRule.REQUIRED
                else "[grey]No[/grey]",
                str(segment.length),
                str(segment.length_tolerance),
                str(segment.id),
            )
    else:
        monopartite_segment = otu.plan.segments[0]

        segment_name = (
            monopartite_segment.name
            if monopartite_segment.name is not None
            else "Unnamed"
        )

        plan_table.add_row(
            segment_name,
            "[red]Yes[/red]",
            str(monopartite_segment.length),
            str(monopartite_segment.length_tolerance),
        )

    console.print(plan_table)

    console.line()
    console.print("[bold]ISOLATES[/bold]")

    for isolate in otu.isolates:
        console.line()
        console.print(str(isolate.name)) if isolate.name is not None else console.print(
            "[UNNAMED]",
        )
        console.line()

        isolate_table = Table(
            box=None,
        )

        isolate_table.add_column("ACCESSION", width=max_accession_length)
        isolate_table.add_column("LENGTH")
        isolate_table.add_column("SEGMENT", min_width=max_segment_name_length)
        isolate_table.add_column("DEFINITION")

        for sequence in sorted(isolate.sequences, key=lambda s: s.accession):
            isolate_table.add_row(
                _render_nucleotide_link(str(sequence.accession)),
                str(len(sequence.sequence)),
                str(otu.plan.get_segment_by_id(sequence.segment).name),
                sequence.definition,
            )

        console.print(isolate_table)


def print_otu_list(otus: Iterator[OTUMinimal]) -> None:
    """Print a list of OTUs to the console.

    :param otus: The OTUs to print.

    """
    table = Table(
        box=None,
    )

    table.add_column("NAME")
    table.add_column("ACRONYM")
    table.add_column("TAXID")
    table.add_column("ID")

    added_rows = False

    for otu in otus:
        table.add_row(
            otu.name,
            otu.acronym,
            str(otu.taxid),
            str(otu.id),
        )

        added_rows = True

    if added_rows:
        console.print(table)
    else:
        console.print("No OTUs found")


console = rich.console.Console()
