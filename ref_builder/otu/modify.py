from collections.abc import Collection
from uuid import UUID

from pydantic import ValidationError
from structlog import get_logger

from ref_builder.ncbi.client import NCBIClient, TaxonLevelError
from ref_builder.ncbi.models import NCBITaxonomy
from ref_builder.otu.utils import (
    DeleteRationale,
    assign_segment_id_to_record,
    create_segments_from_records,
)
from ref_builder.plan import (
    Plan,
    SegmentName,
    SegmentRule,
)
from ref_builder.repo import Repo
from ref_builder.resources import RepoOTU, RepoSequence

logger = get_logger("otu.modify")


def exclude_accessions_from_otu(
    repo: Repo,
    otu: RepoOTU,
    accessions: Collection[str],
) -> None:
    """Exclude accessions from future addition to an OTU."""
    original_excluded_accessions = otu.excluded_accessions.copy()

    with repo.use_transaction():
        excluded_accessions = repo.exclude_accessions(
            otu_id=otu.id, accessions=accessions
        )

    if excluded_accessions != original_excluded_accessions:
        logger.info(
            "Updated excluded accession list.",
            otu_id=str(otu.id),
            excluded_accessions=sorted(excluded_accessions),
        )

    if excluded_accessions == original_excluded_accessions:
        logger.info(
            "Excluded accession list already up to date.",
            excluded_accessions=excluded_accessions,
        )


def allow_accessions_into_otu(
    repo: Repo,
    otu: RepoOTU,
    accessions: Collection[str],
) -> None:
    """Allow accessions for future addition to an OTU.

    This reverses the effect of exclude_accessions_from_otu.
    """
    original_excluded_accessions = otu.excluded_accessions.copy()

    with repo.use_transaction():
        excluded_accessions = repo.allow_accessions(
            otu_id=otu.id, accessions=accessions
        )

    if excluded_accessions != original_excluded_accessions:
        logger.info(
            "Updated excluded accession list.",
            otu_id=str(otu.id),
            excluded_accessions=sorted(excluded_accessions),
        )

    if excluded_accessions == original_excluded_accessions:
        logger.info(
            "Excluded accession list already up to date.",
            excluded_accessions=excluded_accessions,
        )


def update_otu_identifiers(
    repo: Repo,
    otu: RepoOTU,
    taxid: int,
    ignore_cache: bool = False,
) -> RepoOTU | None:
    """Update OTU with new Taxonomy ID and associated name."""
    client = NCBIClient(ignore_cache)

    otu_logger = logger.bind(
        id=str(otu.id), original_taxid=otu.taxid, original_name=otu.name
    )

    if repo.get_otu_id_by_taxid(taxid):
        raise ValueError(
            f"Taxonomy ID {taxid} has already been added to this reference.",
        )

    taxonomy = client.fetch_taxonomy_record(taxid)

    if taxonomy is None:
        otu_logger.fatal(f"Could not retrieve {taxid} from NCBI Taxonomy")
        return None

    otu_logger.info("Updating OTU with new Taxonomy data...")

    return replace_otu_taxonomy_from_record(
        repo,
        otu,
        taxonomy,
    )


def delete_isolate_from_otu(repo: Repo, otu: RepoOTU, isolate_id: UUID) -> bool:
    """Remove an isolate from a specified OTU."""
    otu_logger = logger.bind(otu_id=str(otu.id), taxid=otu.taxid)

    if isolate_id == otu.representative_isolate:
        otu_logger.error(
            "The representative isolate cannot be deleted from the OTU.",
            isolate_id=str(isolate_id),
        )
        return False

    isolate = otu.get_isolate(isolate_id)

    if not isolate:
        otu_logger.error("Isolate not found.", isolate_id=str(isolate_id))
        return False

    with repo.use_transaction():
        repo.delete_isolate(otu.id, isolate.id, rationale=DeleteRationale.USER)

    otu_logger.info(
        "Isolate removed.",
        name=isolate.name,
        removed_isolate_id=str(isolate_id),
        current_isolate_ids=[str(isolate_id_) for isolate_id_ in otu.isolate_ids],
    )

    repo.get_otu(otu.id)

    return True


def set_plan(
    repo: Repo,
    otu: RepoOTU,
    plan: Plan,
) -> Plan | None:
    """Set an OTU's plan to the passed ``plan``."""
    log = logger.bind(name=otu.name, taxid=otu.taxid, plan=plan.model_dump())

    try:
        with repo.use_transaction():
            repo.set_plan(otu.id, plan)
    except ValueError:
        log.exception()
        return None

    return repo.get_otu(otu.id).plan


def set_plan_length_tolerances(
    repo: Repo,
    otu: RepoOTU,
    tolerance: float,
) -> Plan | None:
    """Sets a plan's length tolerances to a new float value."""
    try:
        new_plan = otu.plan.model_copy()
        for segment in new_plan.segments:
            segment.length_tolerance = tolerance
    except ValidationError as exc:
        for error in exc.errors():
            logger.error(
                "Length tolerance must be between 0.0 and 1.0.",
                error_type=error["type"],
                name=otu.name,
                taxid=otu.taxid,
                requested_tolerance=tolerance,
            )
        return None

    with repo.use_transaction():
        repo.set_plan(
            otu.id,
            plan=new_plan,
        )

    return repo.get_otu(otu.id).plan


def add_segments_to_plan(
    repo: Repo,
    otu: RepoOTU,
    rule: SegmentRule,
    accessions: Collection[str],
    ignore_cache: bool = False,
) -> set[UUID]:
    """Add new segments to a multipartite plan."""
    log = logger.bind(name=otu.name, taxid=otu.taxid, rule=str(rule))

    sequence_length_tolerance = repo.settings.default_segment_length_tolerance

    if otu.plan.monopartite and otu.plan.segments[0].name is None:
        log.error(
            "Cannot add new segments unless all existing segments are named.",
            segment_id=str(otu.plan.segments[0].id),
            segment_name=otu.plan.segments[0].name,
        )

        return set()

    client = NCBIClient(ignore_cache)

    if not (records := client.fetch_genbank_records(accessions)):
        log.error(
            "Could not fetch records associated with requested accessions.",
            accessions=sorted(accessions),
        )
        return set()

    if len(records) < len(accessions):
        log.error(
            "Not all requested accessions could be fetched.",
            requested_accessions=sorted(accessions),
            accessible_record_accessions=sorted(
                record.accession_version for record in records
            ),
        )
        return set()

    new_segments = create_segments_from_records(
        records,
        rule,
        length_tolerance=sequence_length_tolerance,
    )
    if not new_segments:
        log.warning("No segments can be added.")
        return set()

    if len(new_segments) < len(accessions):
        log.error("Could not create all new segments.")
        return set()

    new_plan = otu.plan.model_copy()
    new_plan.segments.extend(new_segments)

    set_plan(repo, otu, new_plan)

    log.info("Added new segments", ids=[str(segment.id) for segment in new_segments])

    return {segment.id for segment in new_segments}


def rename_plan_segment(
    repo: Repo,
    otu: RepoOTU,
    segment_id: UUID,
    segment_name: SegmentName,
) -> Plan | None:
    """Set a new name for the segment matching the given ID if possible, else return None."""
    log = logger.bind(
        name=otu.name,
        taxid=otu.taxid,
        segment_id=str(segment_id),
        segment_name=str(segment_name),
    )

    modified_plan = otu.plan.model_copy()

    for segment in modified_plan.segments:
        if segment.id == segment_id:
            segment.name = segment_name

            return set_plan(repo, otu, modified_plan)

    log.error("Segment with requested ID not found.")

    return None


def replace_sequence_in_otu(
    repo: Repo,
    otu: RepoOTU,
    new_accession: str,
    replaced_accession: str,
    ignore_cache: bool = False,
) -> RepoSequence | None:
    """Replace a sequence in an OTU."""
    ncbi = NCBIClient(ignore_cache)

    replaced_sequence = otu.get_sequence_by_accession(replaced_accession)

    if replaced_sequence is None:
        logger.error(
            "This accession does not exist in this OTU.", accession=replaced_accession
        )
        return None

    affected_isolate_ids = otu.get_isolate_ids_containing_sequence_id(
        replaced_sequence.id
    )
    if not affected_isolate_ids:
        logger.warning(
            "This sequence is not linked to any isolates.",
            accession=str(replaced_sequence.accession),
            sequence_id=replaced_sequence.id,
        )
        return None

    record = ncbi.fetch_genbank_records([new_accession])[0]

    segment_id = assign_segment_id_to_record(record, otu.plan)
    if segment_id is None:
        logger.error("This segment does not match the plan.")

    with repo.use_transaction():
        new_sequence = repo.create_sequence(
            otu.id,
            accession=record.accession_version,
            definition=record.definition,
            legacy_id=None,
            segment=segment_id,
            sequence=record.sequence,
        )

        for isolate_id in affected_isolate_ids:
            repo.replace_sequence(
                otu.id,
                isolate_id=isolate_id,
                sequence_id=new_sequence.id,
                replaced_sequence_id=replaced_sequence.id,
                rationale="Requested by user",
            )

    if new_sequence is not None:
        logger.info(
            f"{replaced_accession} replaced by {new_sequence.accession}.",
            new_sequence_id=new_sequence.id,
        )
        return new_sequence

    logger.error(f"{replaced_accession} could not be replaced.")


def set_representative_isolate(
    repo: Repo,
    otu: RepoOTU,
    isolate_id: UUID,
) -> UUID | None:
    """Set an OTU's representative isolate to a given existing isolate ID.

    Returns the isolate ID if successful, else None.
    """
    otu_logger = logger.bind(name=otu.name, taxid=otu.taxid)

    new_representative_isolate = otu.get_isolate(isolate_id)
    if new_representative_isolate is None:
        otu_logger.error(
            "Isolate not found. Consider adding a new isolate.",
            requested_isolate=str(isolate_id),
        )
        return None

    if otu.representative_isolate is not None:
        if otu.representative_isolate == new_representative_isolate.id:
            otu_logger.warning(
                "Redundant replacement attempt.",
                isolate_id=str(otu.representative_isolate),
            )
            return otu.representative_isolate

        otu_logger.warning(
            "Replacing representative isolate.",
            old_isolate_id=str(otu.representative_isolate),
            new_isolate_id=str(new_representative_isolate.id),
        )

    with repo.use_transaction():
        repo.set_representative_isolate(otu.id, new_representative_isolate.id)

    otu_logger.info(
        "New representative isolate set.",
        representative_isolate_id=str(new_representative_isolate.id),
    )

    return new_representative_isolate.id


def replace_otu_taxonomy_from_record(
    repo: Repo, otu: RepoOTU, taxon_record: NCBITaxonomy
) -> RepoOTU | None:
    """Replace an OTU taxonomy data using a given NCBI taxonomy record."""
    with repo.use_transaction():
        repo.update_otu_identifiers(
            otu_id=otu.id, taxid=taxon_record.id, name=taxon_record.name
        )

    return repo.get_otu(otu.id)


def fetch_taxonomy_based_on_otu_contents(
    otu: RepoOTU, ignore_cache: bool = False
) -> NCBITaxonomy:
    """Given an OTU, fetch a taxonomy record based on its contents."""
    ncbi = NCBIClient(ignore_cache)
    rep_isolate = otu.get_isolate(otu.representative_isolate)

    records = ncbi.fetch_genbank_records(rep_isolate.accessions)

    record_taxid = records[0].source.taxid

    return ncbi.fetch_taxonomy_record(record_taxid)


def correct_otu_ranks(repo: Repo, ignore_cache: bool = False) -> None:
    """Correct OTU metadata if the Taxonomy rank is above or below species-level."""
    ncbi = NCBIClient(ignore_cache)

    sub_species_taxids = set()
    over_ranked_taxids = set()

    updated_otu_ids = set()

    for otu in repo.iter_otus():
        try:
            taxon_record = ncbi.fetch_taxonomy_record(otu.taxid)
        except TaxonLevelError:
            logger.warning(
                "OTU is too high level to work with.",
                msg=str(TaxonLevelError),
            )

            over_ranked_taxids.add(otu.taxid)

            taxon_record = fetch_taxonomy_based_on_otu_contents(otu, ignore_cache)

        if taxon_record.rank != "species":
            logger.debug(
                "OTU has a sub-species level taxonomy listing.",
                rank=taxon_record.rank,
                species=taxon_record.species,
            )
            sub_species_taxids.add(otu.taxid)

            logger.debug(
                "Fetching species-level taxonomy...",
                species_taxid=taxon_record.species.id,
                species_name=taxon_record.species.name,
            )

            species_taxon_record = ncbi.fetch_taxonomy_record(taxon_record.species.id)

            logger.debug(
                "Updating OTU with new Taxonomy data...",
                otu_id=otu.id,
                current_name=otu.name,
                current_taxid=otu.taxid,
                new_name=species_taxon_record.name,
                new_taxid=species_taxon_record.id,
            )

            updated_otu = replace_otu_taxonomy_from_record(
                repo,
                otu=otu,
                taxon_record=species_taxon_record,
            )

            if updated_otu is not None:
                updated_otu_ids.add(updated_otu.id)

    if updated_otu_ids:
        logger.info(
            "OTU identifiers updated",
            otu_ids=updated_otu_ids,
            over_ranked_count=len(over_ranked_taxids),
            sub_species_count=len(sub_species_taxids),
        )
