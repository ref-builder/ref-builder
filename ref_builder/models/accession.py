import re
from dataclasses import dataclass

REFSEQ_ACCESSION_PATTERN = re.compile(pattern=r"^NC_[0-9A-Z]+$")
"""RefSeq accession pattern for viral complete genomic molecules (NC_ prefix).

The identifier following the underscore can be alphanumeric and of variable length.

Examples: NC_003619, NC_010314, NC_ABC123
"""


@dataclass(frozen=True, order=True)
class Accession:
    """A Genbank accession number."""

    key: str
    """The accession key.

    In the accession "MN908947.3", the key is "MN908947".
    """

    version: int
    """The version number.

    In the accession "MN908947.3", the version is 3.
    """

    @classmethod
    def from_string(cls, string: str) -> "Accession":
        """Create an Accession from a versioned accession string,
        e.g. "MN908947.3".
        """
        if not string or not string.strip():
            raise ValueError("Accession string cannot be empty or whitespace.")

        parts = string.split(".")

        if len(parts) < 2:
            raise ValueError(
                f'Accession string "{string}" does not contain two parts '
                "delimited by a period. Expected format: KEY.VERSION (e.g., NC_123456.1)"
            )

        if len(parts) > 2:
            raise ValueError(
                f'Accession string "{string}" contains multiple periods. '
                "Expected format: KEY.VERSION (e.g., NC_123456.1)"
            )

        key, string_version = parts

        if not string_version.isdigit():
            raise ValueError(
                f"Accession version ({string_version}) is not an integer. "
                "Expected format: KEY.VERSION (e.g., NC_123456.1)"
            )

        version = int(string_version)

        return Accession(key=key, version=version)

    def __str__(self) -> str:
        """Return the accession as a string."""
        return f"{self.key}.{self.version}"

    def __repr__(self) -> str:
        """Return the accession as a repr string."""
        return f"Accession({self.key}, {self.version})"

    @property
    def is_refseq(self) -> bool:
        """Return True if this accession is from NCBI's RefSeq database."""
        return REFSEQ_ACCESSION_PATTERN.match(self.key) is not None
