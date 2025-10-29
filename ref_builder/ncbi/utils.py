import re
from collections import defaultdict

from ref_builder.models.accession import Accession
from ref_builder.models.isolate import IsolateName, IsolateNameType
from ref_builder.ncbi.models import NCBIGenbank


def parse_refseq_comment(comment: str) -> tuple[str, str]:
    """Parse a standard RefSeq comment."""
    if not comment:
        raise ValueError("Empty comment")

    pattern = re.compile(r"^(\w+ REFSEQ): [\w ]+. [\w ]+ (\w+).")

    if match := pattern.search(comment):
        return match.group(1), match.group(2)

    raise ValueError("Invalid RefSeq comment")


def group_genbank_records_by_isolate(
    records: list[NCBIGenbank],
) -> dict[IsolateName, dict[Accession, NCBIGenbank]]:
    """Indexes Genbank records by isolate name."""
    isolates = defaultdict(dict)

    for record in records:
        isolate_name = _extract_isolate_name_from_record(record)
        if isolate_name is None:
            # Assume this is a monopartite OTU and do not group.
            continue

        versioned_accession = Accession.from_string(record.accession_version)
        isolates[isolate_name][versioned_accession] = record

    return isolates


def _extract_isolate_name_from_record(record: NCBIGenbank) -> IsolateName | None:
    """Get the isolate name from a Genbank record."""
    for source_type in IsolateNameType:
        if value := getattr(record.source, source_type, None):
            return IsolateName(type=source_type, value=value)

    return None
