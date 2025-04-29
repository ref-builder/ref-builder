import datetime

from pydantic import UUID4, BaseModel, Field

from ref_builder.utils import DataType


class RepoMeta(BaseModel):
    """Represents the metadata for a Virtool reference repository."""

    id: UUID4
    """The repository id."""

    created_at: datetime.datetime
    """The date and time the repository was created."""

    data_type: DataType
    """The repository data type."""

    name: str
    """The repository name."""

    organism: str
    """The repository organism."""


class RepoSettings(BaseModel):
    """The default settings of a Virtool reference repository."""

    default_segment_length_tolerance: float = Field(0.03, ge=0.0, le=1.0)
    """The deviation a sequence is allowed from its plan segment's length before it
    fails validation.
    """
