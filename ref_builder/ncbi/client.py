"""A client for fetching data from NCBI databases."""

import datetime
import os
from collections.abc import Collection
from contextlib import contextmanager, suppress
from enum import StrEnum
from http import HTTPStatus
from typing import Protocol
from urllib.error import HTTPError

from Bio import Entrez
from pydantic import ValidationError
from structlog import get_logger

from ref_builder.models.accession import Accession
from ref_builder.models.lineage import Lineage, Taxon, TaxonOtherNames
from ref_builder.ncbi.cache import NCBICache
from ref_builder.ncbi.models import (
    NCBIDatabase,
    NCBIGenbank,
    NCBIRank,
    NCBITaxonomy,
)

if email := os.environ.get("NCBI_EMAIL"):
    Entrez.email = email

if api_key := os.environ.get("NCBI_API_KEY"):
    Entrez.api_key = os.environ.get("NCBI_API_KEY")

logger = get_logger("ncbi")

ESEARCH_PAGE_SIZE = 1000
"""The number of results to fetch per page in an Entrez esearch query."""

DATE_TEMPLATE = "%Y/%m/%d"
"""The standard date format used by NCBI Entrez."""


class NCBIClientProtocol(Protocol):
    """Protocol defining the interface for NCBI client implementations."""

    ignore_cache: bool

    def fetch_genbank_records(
        self,
        accessions: Collection[str | Accession],
    ) -> list[NCBIGenbank]: ...

    def fetch_taxonomy_record(self, taxid: int) -> NCBITaxonomy | None: ...

    def fetch_descendant_taxids(self, species_taxid: int) -> list[int]: ...

    @staticmethod
    def fetch_accessions_by_taxid(
        taxid: int,
        sequence_min_length: int = 0,
        sequence_max_length: int = 0,
        refseq_only: bool = False,
    ) -> list[Accession]: ...

    @staticmethod
    def filter_accessions(raw_accessions: Collection[str]) -> set[Accession]: ...


class TaxonLevelError(ValueError):
    """Raised when a fetched taxonomy record is above species level."""


class GenbankRecordKey(StrEnum):
    """Genbank record keys."""

    ACCESSION_VERSION = "GBSeq_accession-version"
    COMMENT = "GBSeq_comment"
    DEFINITION = "GBSeq_definition"
    FEATURE_TABLE = "GBSeq_feature-table"
    LENGTH = "GBSeq_length"
    PRIMARY_ACCESSION = "GBSeq_primary-accession"
    SEQUENCE = "GBSeq_sequence"


class NCBIClient:
    """A client for fetching, caching, and validating data from NCBI databases."""

    def __init__(self, ignore_cache: bool) -> None:
        """Initialize the NCBI client with a cache path and an ignore_cache flag.

        :param ignore_cache: If True, does not allow the return of cached data
        """
        self.cache = NCBICache()
        self.ignore_cache = ignore_cache

    def fetch_genbank_records(
        self,
        accessions: Collection[str | Accession],
    ) -> list[NCBIGenbank]:
        """Fetch or load NCBI Nucleotide records corresponding to a list of accessions.

        Cache fetched records if found. Returns validated records.

        :param accessions: A list of accessions to be fetched
        :return: A list of validated NCBIGenbank records
        """
        if not accessions:
            return []

        records = []
        cached_keys = set()

        log = logger.bind(accessions=accessions)

        if not self.ignore_cache:
            uncached_accessions = []

            for accession in accessions:
                if isinstance(accession, Accession):
                    accession_key = accession.key
                    version = accession.version
                else:
                    try:
                        versioned_accession = Accession.from_string(accession)
                        accession_key = versioned_accession.key
                        version = versioned_accession.version
                    except ValueError:
                        accession_key = accession
                        version = "*"

                record = self.cache.load_genbank_record(accession_key, version)

                if record is not None:
                    records.append(record)
                    cached_keys.add(accession_key)

            if records:
                log.debug(
                    "Loaded cached records",
                    record_count=len(records),
                    cached_records=[
                        record.get(GenbankRecordKey.PRIMARY_ACCESSION)
                        for record in records
                    ],
                )
            if uncached_accessions:
                log.debug(
                    "Uncached accessions found",
                    uncached_accessions=uncached_accessions,
                )

        # Build list of accessions to fetch by excluding those found in cache
        accessions_to_fetch = []
        for accession in accessions:
            if isinstance(accession, Accession):
                key = accession.key
            else:
                try:
                    key = Accession.from_string(accession).key
                except ValueError:
                    key = accession

            if key not in cached_keys:
                accessions_to_fetch.append(accession)

        if accessions_to_fetch:
            log.debug("Fetching accessions...", fetch_list=accessions_to_fetch)
            new_records = NCBIClient.fetch_unvalidated_genbank_records(
                accessions_to_fetch
            )

            for record in new_records:
                versioned_accession = Accession.from_string(
                    record[GenbankRecordKey.ACCESSION_VERSION],
                )
                self.cache.cache_genbank_record(
                    record,
                    versioned_accession.key,
                    versioned_accession.version,
                )

            if new_records:
                records.extend(new_records)

        if records:
            return sorted(
                NCBIClient._validate_genbank_records(records),
                key=lambda r: r.accession,
            )

        return []

    @staticmethod
    def fetch_unvalidated_genbank_records(accessions: Collection[str]) -> list[dict]:
        """Fetch a list of Genbank records given a list of accessions.

        :param accessions: List of accession numbers to fetch from GenBank
        :return: A list of deserialized XML records from NCBI Nucleotide
        """
        log = logger.bind(accessions=accessions)

        try:
            with log_http_error():
                try:
                    handle = Entrez.efetch(
                        db=NCBIDatabase.NUCCORE,
                        id=list(accessions),
                        rettype="gb",
                        retmode="xml",
                    )
                except RuntimeError as e:
                    log.warning("Bad ID.", exception=e)
                    return []

        except HTTPError as e:
            if e.code == HTTPStatus.BAD_REQUEST:
                log.exception("Accessions not found")
            else:
                log.exception("HTTPError")

            return []

        try:
            records = Entrez.read(handle)
        except RuntimeError:
            log.exception("NCBI returned unparseable data")
            return []

        if records:
            # Handle cases where not all accessions can be fetched
            if len(records) != len(accessions):
                log.debug("Partial results fetched. Returning results...")

            return records

        return []

    @staticmethod
    def fetch_accessions_by_taxid(
        taxid: int,
        sequence_min_length: int = 0,
        sequence_max_length: int = 0,
        refseq_only: bool = False,
    ) -> list[Accession]:
        """Fetch all accessions associated with the given ``taxid``.

        :param taxid: A Taxonomy ID
        :param sequence_min_length: The minimum length of a fetched sequence.
        :param sequence_max_length: The maximum length of a fetched sequence.
        :param refseq_only: Only fetch accessions from NCBI RefSeq database.:
        :return: A list of Accession objects
        """
        log = logger.bind(taxid=taxid)

        page = 1
        accessions = []

        search_terms = [f"txid{taxid}[orgn]"]

        if sequence_min_length > 0 and sequence_max_length > 0:
            search_terms.append(
                NCBIClient.generate_sequence_length_filter_string(
                    sequence_min_length,
                    sequence_max_length,
                )
            )

        if refseq_only:
            search_terms.append("refseq[filter]")

        search_term_string = " AND ".join(search_terms)

        log.debug(
            "Fetching NCBI Nucleotide accessions associated with NCBI Taxonomy ID...",
            search_string=search_term_string,
        )

        while True:
            retstart = (page - 1) * ESEARCH_PAGE_SIZE

            with log_http_error():
                handle = Entrez.esearch(
                    db=NCBIDatabase.NUCCORE,
                    term=search_term_string,
                    idtype="acc",
                    retstart=retstart,
                    retmax=ESEARCH_PAGE_SIZE,
                )

            result = Entrez.read(handle)

            if not result["IdList"]:
                break

            result_count = int(result["Count"])

            accessions += result["IdList"]

            if result_count - retstart <= ESEARCH_PAGE_SIZE:
                break

            log.debug(
                "Large fetch. May take longer than expected...",
                result_count=result_count,
                page=page,
                page_size=ESEARCH_PAGE_SIZE,
            )

            page += 1

        return list(NCBIClient.filter_accessions(accessions))

    @staticmethod
    def _validate_genbank_records(records: list[dict]) -> list[NCBIGenbank]:
        """Process a list of raw Genbank dicts into validated NCBIGenbank records.
        Logs an error if there is an issue with validation or parsing,
        but does not fail out.

        :param records: A list of unvalidated NCBI Genbank records
        :return: A list of validated records as NCBIGenbank
        """
        clean_records = []

        for record in records:
            accession = record.get(GenbankRecordKey.PRIMARY_ACCESSION, "?")

            try:
                clean_records.append(NCBIGenbank.model_validate(record))

            except ValidationError as exc:
                logger.debug(
                    "Encountered validation errors",
                    accession=accession,
                    count=exc.error_count(),
                    errors=exc.errors(),
                )
                continue

        return clean_records

    def fetch_taxonomy_record(self, taxid: int) -> NCBITaxonomy | None:
        """Fetch and validate a taxonomy record from NCBI Taxonomy.

        If the record rank has an invalid rank (e.g. "no data"), makes an additional
        docsum fetch and attempts to extract the rank data.

        :param taxid: A NCBI Taxonomy id
        :return: A validated NCBI Taxonomy record NCBITaxonomy if possible,
            else None
        """
        log = logger.bind(taxid=taxid)

        record = None if self.ignore_cache else self.cache.load_taxonomy(taxid)

        if record is None:
            with log_http_error():
                try:
                    handle = Entrez.efetch(
                        db=NCBIDatabase.TAXONOMY,
                        id=taxid,
                        rettype="null",
                    )
                except HTTPError:
                    return None

            try:
                records = Entrez.read(handle)
            except RuntimeError:
                log.exception("NCBI returned unparseable data.")
                return None

            if records:
                record = records[0]
                self.cache.cache_taxonomy_record(record, taxid)
            else:
                return None

        try:
            return NCBITaxonomy.model_validate(record)

        except ValidationError as e:
            for error in e.errors():
                log.warning(
                    "ValidationError",
                    msg=error["msg"],
                    loc=error["loc"],
                    type=error["type"],
                )

                if error["type"] == "taxon_rank_too_high":
                    raise TaxonLevelError(error["msg"])

        return None

    def fetch_descendant_taxids(self, species_taxid: int) -> list[int]:
        """Fetch all descendant taxids under a species.

        Uses NCBI taxonomy subtree search to find all taxa under the given species.
        Returns only subspecific taxids (isolate, "no rank" ranks).

        :param species_taxid: A NCBI species-level Taxonomy id
        :return: A list of subspecific taxids under the species
        """
        log = logger.bind(species_taxid=species_taxid)

        log.debug("Fetching descendant taxids...")

        with log_http_error():
            handle = Entrez.esearch(
                db=NCBIDatabase.TAXONOMY,
                term=f"txid{species_taxid}[Subtree]",
                retmax=10000,
            )

        result = Entrez.read(handle)

        if not result["IdList"]:
            return []

        descendant_taxids = [int(tid) for tid in result["IdList"]]

        subspecific_taxids = []
        for taxid in descendant_taxids:
            if taxid == species_taxid:
                continue

            tax_record = self.fetch_taxonomy_record(taxid)
            if tax_record and tax_record.rank in (NCBIRank.ISOLATE, NCBIRank.NO_RANK):
                subspecific_taxids.append(taxid)

        log.debug(
            "Found subspecific descendants",
            count=len(subspecific_taxids),
            taxids=subspecific_taxids,
        )

        return subspecific_taxids

    def fetch_lineage(self, taxid: int) -> Lineage:
        """Fetch a complete lineage from species including all descendant taxids.

        Fetches the taxonomy record for the given taxid and builds a lineage
        from species level including ALL subspecific descendants. This ensures
        that any isolate under the species can be added to an OTU without
        validation errors.

        :param taxid: A NCBI Taxonomy id
        :return: A Lineage object with species and all subspecific descendants
        """
        log = logger.bind(taxid=taxid)

        taxonomy = self.fetch_taxonomy_record(taxid)

        if taxonomy is None:
            raise ValueError(f"Could not fetch taxonomy record for taxid {taxid}")

        species_taxonomy = self.fetch_taxonomy_record(taxonomy.species.id)
        if species_taxonomy is None:
            raise ValueError(
                f"Could not fetch species taxonomy for taxid {taxonomy.species.id}"
            )

        species_taxon = Taxon(
            id=species_taxonomy.id,
            name=species_taxonomy.name,
            parent=None,
            rank=species_taxonomy.rank,
            other_names=TaxonOtherNames(
                acronym=species_taxonomy.other_names.acronym,
                synonyms=species_taxonomy.other_names.equivalent_name,
            ),
        )

        descendant_taxids = self.fetch_descendant_taxids(species_taxonomy.id)

        taxon_objects = [species_taxon]

        for descendant_taxid in descendant_taxids:
            tax_record = self.fetch_taxonomy_record(descendant_taxid)
            if tax_record is None:
                log.warning("Could not fetch taxonomy record", taxid=descendant_taxid)
                continue

            parent_id = species_taxonomy.id
            for lineage_item in tax_record.lineage:
                if lineage_item.rank in ("no rank", "isolate"):
                    if lineage_item.id in descendant_taxids:
                        parent_id = lineage_item.id
                        break

            taxon = Taxon(
                id=tax_record.id,
                name=tax_record.name,
                parent=parent_id,
                rank=tax_record.rank,
                other_names=TaxonOtherNames(
                    acronym=tax_record.other_names.acronym,
                    synonyms=tax_record.other_names.equivalent_name,
                ),
            )

            taxon_objects.append(taxon)

        if not taxon_objects:
            raise ValueError(f"Could not build lineage for taxid {taxid}")

        return Lineage(taxa=taxon_objects)

    @staticmethod
    def filter_accessions(raw_accessions: Collection[str]) -> set[Accession]:
        """Filter raw eSearch accession list and return a set of compliant Nucleotide accessions."""
        valid_accessions = set()
        for accession in raw_accessions:
            with suppress(ValueError):
                valid_accessions.add(Accession.from_string(accession))

        return valid_accessions

    @staticmethod
    def generate_sequence_length_filter_string(
        min_length: int = 0, max_length: int = 0
    ) -> str:
        """Return a term filter string delimiting a given length range.

        Returns an empty string if not given a min or max length parameter.
        """
        if min_length > 0:
            if min_length > 0:
                return f'"{min_length}"[SLEN] : "{max_length}"[SLEN]'

            return f'"{min_length}"[SLEN]'

        if max_length > 0:
            return f'"{max_length}"[SLEN]'

        return ""

    @staticmethod
    def generate_date_filter_string(
        filter_type: str,
        start_date: datetime.date | None = None,
        end_date: datetime.date | None = None,
    ) -> str:
        """Return a term filter string delimiting a given time range.

        Returns an empty string if not given a start or end date parameter.

        :param filter_type: The search term filter type. Can be "MDAT" or "PDAT".:
        :param start_date: The start date of the search range. If None, assume no lower bound.:
        :param end_date: The end date of the search range. If None, assume no upper bound.:
        :return: A formatted date filter string or an empty string.
        """
        if filter_type not in {"MDAT", "PDAT"}:
            raise ValueError(
                "Invalid filter type. Only ``MDAT``, ``PDAT`` are supported."
            )

        if start_date is None and end_date is None:
            return ""

        start_date_string = "0001/01/01"
        end_date_string = "3000/12/31"

        if start_date:
            start_date_string = start_date.strftime(DATE_TEMPLATE)

            if end_date:
                end_date_string = end_date.strftime(DATE_TEMPLATE)

        return " ".join(
            [
                f'"{start_date_string}"[{filter_type}]',
                ":",
                f'"{end_date_string}"[{filter_type}]',
            ]
        )


@contextmanager
def log_http_error() -> None:
    """Log detailed HTTPError info for debugging before throwing the HTTPError."""
    try:
        yield
    except HTTPError as e:
        logger.exception(
            "HTTPError raised",
            body=e.read(),
            code=e.code,
            reason=e.reason,
        )
        raise
