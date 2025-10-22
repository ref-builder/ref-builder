from dataclasses import dataclass

from pydantic import UUID4, BaseModel

from ref_builder.isolate import IsolateNameType


@dataclass(frozen=True)
class IsolateName:
    """A name for an isolate.

    The isolate name consists of a type and a value.

    For example, in the isolate name "Isolate PPSMV2-Badnapur", the type is "isolate"
    and the value is "PPSMV2-Badnapur.

    In the Genbank record for a sequence that belongs to this isolate, the source table
    contains:

    .. code-block:: text

        /isolate="PPSMV2-Badnapur"

    """

    type: IsolateNameType
    """The type of sub-species categorization."""

    value: str
    """The name of this subcategory."""

    def __str__(self) -> str:
        """Return the isolate name as a formatted string."""
        return f"{self.type.capitalize()} {self.value}"


class IsolateModel(BaseModel):
    """A class representing the fields of an isolate."""

    id: UUID4
    """The isolate id."""

    name: IsolateName | None
    """The isolate's name."""
