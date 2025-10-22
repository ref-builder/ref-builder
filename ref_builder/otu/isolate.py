from uuid import UUID

from ref_builder.ncbi.models import NCBIGenbank
from ref_builder.otu.builders.otu import OTUBuilder
from ref_builder.otu.builders.sequence import SequenceBuilder
from ref_builder.repo import Repo


def create_sequence_from_record(
    repo: Repo,
    otu: OTUBuilder,
    record: NCBIGenbank,
    segment_id: UUID,
) -> SequenceBuilder | None:
    """Take a NCBI Nucleotide record and create a new sequence."""
    return repo.create_sequence(
        otu.id,
        accession=record.accession_version,
        definition=record.definition,
        segment=segment_id,
        sequence=record.sequence,
    )
