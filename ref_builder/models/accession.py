from dataclasses import dataclass


@dataclass(frozen=True)
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
        try:
            key, string_version = string.split(".")
        except ValueError as e:
            if "not enough values to unpack" in str(e):
                raise ValueError(
                    f'Given accession string "{string}" does not contain two parts'
                    "delimited by a period."
                ) from e

            raise

        if string_version.isdigit():
            version = int(string_version)
        else:
            raise ValueError(f"Accession version ({string_version})is not an integer.")

        return Accession(key=key, version=version)

    def __eq__(self, other: "Accession") -> bool:
        if isinstance(other, Accession):
            return self.key == other.key and self.version == other.version

        raise ValueError(
            f"Invalid comparison against invalid value {other} (type {type(other)})"
        )

    def __lt__(self, other: "Accession") -> bool:
        if isinstance(other, Accession):
            return self.key < other.key or self.version < other.version

        raise ValueError(
            f"Invalid comparison against invalid value {other} (type {type(other)})"
        )

    def __gt__(self, other: "Accession") -> bool:
        if isinstance(other, Accession):
            return self.key > other.key or self.version > other.version
        raise ValueError(
            f"Invalid comparison against invalid value {other} (type {type(other)})"
        )

    def __str__(self) -> str:
        """Return the accession as a string."""
        return f"{self.key}.{self.version}"
