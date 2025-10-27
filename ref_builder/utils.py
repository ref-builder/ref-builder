import re
from enum import StrEnum

from ref_builder.models.accession import Accession

ZERO_PADDING_MAX = 99999999
"""The maximum number that can be padded with zeroes in event IDs and filenames."""

GENBANK_ACCESSION_PATTERN = re.compile(pattern=r"^[A-Z]{1,2}[0-9]{5,6}$")

REFSEQ_ACCESSION_PATTERN = re.compile(pattern=r"^NC_[0-9A-Z]+$")
"""RefSeq accession pattern for viral complete genomic molecules (NC_ prefix).

The identifier following the underscore can be alphanumeric and of variable length.

Examples: NC_003619, NC_010314, NC_ABC123
"""


class ExcludedAccessionAction(StrEnum):
    """Possible actions that can be taken on the excluded/allowed status of an accession."""

    ALLOW = "allow"
    EXCLUDE = "exclude"


def is_accession_key_valid(accession_key: str) -> bool:
    """Return True if the given accession is a valid Genbank or RefSeq accession."""
    return (
        GENBANK_ACCESSION_PATTERN.match(accession_key) is not None
        or REFSEQ_ACCESSION_PATTERN.match(accession_key) is not None
    )


def get_accession_key(raw: str) -> str:
    """Parse a string to check if it follows Genbank or RefSeq accession parameters and
    return the key part only.
    """
    if is_accession_key_valid(raw):
        return raw

    try:
        versioned_accession = Accession.from_string(raw)
    except ValueError:
        raise ValueError("Invalid accession key")

    if GENBANK_ACCESSION_PATTERN.match(
        versioned_accession.key
    ) or REFSEQ_ACCESSION_PATTERN.match(versioned_accession.key):
        return versioned_accession.key

    raise ValueError("Invalid accession key")


def generate_natural_sort_key(string: str) -> list[int | str]:
    """Generate a natural order sorting key for a string.

    This list: ["1", "10", "2"] will be sorted as ["1", "2", "10"], as opposed to
    lexographical order, which would result in ["1", "10", "2"].

    :param string: the string to convert to a sorting key
    :return: the sorting key
    """

    def _convert(string_: str) -> int | str:
        return int(string_) if string_.isdigit() else string_.lower()

    return [_convert(c) for c in re.split("([0-9]+)", string)]


def filter_accessions(
    accessions: list[str],
    blocked: set[str],
) -> list[str]:
    """Filter a list of accessions by removing blocked accessions.

    :param accessions: list of accession strings to filter
    :param blocked: set of blocked accession strings
    :return: filtered list of accessions
    """
    return [acc for acc in accessions if acc not in blocked]


def pad_zeroes(number: int) -> str:
    """Pad a number with zeroes to make it 8 characters long.

    :param number: the number to pad
    :return: the padded number
    """
    if number > ZERO_PADDING_MAX:
        raise ValueError("Number is too large to pad")

    return str(number).zfill(8)
