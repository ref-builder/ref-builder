"""Functions for working with NCBI taxonomy dumps."""

import tarfile
import tempfile
from collections import defaultdict
from datetime import UTC, datetime, timedelta
from pathlib import Path
from urllib.request import urlretrieve

from ref_builder.paths import user_cache_directory_path

TAXDUMP_URL = "https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz"
VIRUSES_TAXID = 10239
CACHE_MAX_AGE = timedelta(days=7)

_taxdump_cache_path = user_cache_directory_path / "ncbi" / "taxdump"


def fetch_plant_virus_taxids() -> set[int]:
    """Fetch taxids for plant-infecting virus species.

    Downloads NCBI taxdump if cache is older than 1 week, then:
    1. Parses nodes.dmp to find virus taxids (descendants of 10239) at species rank
    2. Parses host.dmp for taxids with plant hosts
    3. Returns intersection

    :return: Set of taxids for plant virus species
    """
    taxdump_path = _get_current_taxdump()

    if taxdump_path is None:
        taxdump_path = _download_and_extract_taxdump()

    virus_species_taxids = _parse_virus_species(taxdump_path / "nodes.dmp")
    plant_host_taxids = _parse_plant_hosts(taxdump_path / "host.dmp")

    return virus_species_taxids & plant_host_taxids


def _get_current_taxdump() -> Path | None:
    """Get the path to a valid (non-stale) taxdump cache directory.

    Looks for directories matching `taxdump_YYYY-MM-DDTHH-MM-SS` pattern
    and returns the most recent one if it's not stale.

    :return: Path to valid taxdump directory, or None if stale/missing
    """
    if not _taxdump_cache_path.exists():
        return None

    taxdump_dirs = sorted(
        _taxdump_cache_path.glob("taxdump_*"),
        reverse=True,
    )

    if not taxdump_dirs:
        return None

    latest = taxdump_dirs[0]
    timestamp_str = latest.name.removeprefix("taxdump_")

    try:
        timestamp = datetime.strptime(timestamp_str, "%Y-%m-%dT%H-%M-%S").replace(
            tzinfo=UTC,
        )
    except ValueError:
        return None

    if datetime.now(UTC) - timestamp > CACHE_MAX_AGE:
        return None

    return latest


def _download_and_extract_taxdump() -> Path:
    """Download and extract the NCBI taxdump to the cache directory.

    :return: Path to the extracted taxdump directory
    """
    _taxdump_cache_path.mkdir(exist_ok=True, parents=True)

    timestamp = datetime.now(UTC).strftime("%Y-%m-%dT%H-%M-%S")
    taxdump_path = _taxdump_cache_path / f"taxdump_{timestamp}"

    with tempfile.TemporaryDirectory() as tmpdir:
        tar_path = Path(tmpdir) / "new_taxdump.tar.gz"
        urlretrieve(TAXDUMP_URL, tar_path)

        with tarfile.open(tar_path, "r:gz") as tar:
            tar.extractall(taxdump_path)

    return taxdump_path


def _parse_virus_species(nodes_path: Path) -> set[int]:
    """Parse nodes.dmp to find all virus species taxids.

    Finds all descendants of VIRUSES_TAXID (10239) that have rank 'species'.

    nodes.dmp format: taxid | parent_taxid | rank | ...

    :param nodes_path: Path to nodes.dmp
    :return: Set of virus species taxids
    """
    children: dict[int, list[int]] = defaultdict(list)
    ranks: dict[int, str] = {}

    with open(nodes_path, encoding="utf-8") as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = int(parts[0])
            parent_taxid = int(parts[1])
            rank = parts[2]

            children[parent_taxid].append(taxid)
            ranks[taxid] = rank

    virus_species: set[int] = set()
    stack = [VIRUSES_TAXID]

    while stack:
        taxid = stack.pop()

        if ranks.get(taxid) == "species":
            virus_species.add(taxid)

        stack.extend(children[taxid])

    return virus_species


def _parse_plant_hosts(host_path: Path) -> set[int]:
    """Parse host.dmp to find taxids with plant hosts.

    host.dmp format: taxid | host_name | ...

    :param host_path: Path to host.dmp
    :return: Set of taxids that have plant hosts
    """
    plant_taxids: set[int] = set()

    with open(host_path, encoding="utf-8") as f:
        for line in f:
            if "plants" in line.lower():
                parts = line.split("\t|\t")
                taxid = int(parts[0])
                plant_taxids.add(taxid)

    return plant_taxids
