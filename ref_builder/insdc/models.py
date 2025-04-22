"""Models for data following INSDC standrds."""

import contextlib
from collections.abc import Iterator
from enum import StrEnum
from typing import Annotated, Literal

from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
    field_validator,
)

from ref_builder.models import Molecule, MolType


class FTGenbankLexicon(StrEnum):
    """Table keys used in Genbank."""

    FEATURE_KEY = "GBFeature_key"
    FEATURE_INTERVALS = "GBFeature_intervals"
    FEATURE_QUALIFIERS = "GBFeature_quals"
    INTERVAL_FROM = "GBInterval_from"
    INTERVAL_TO = "GBInterval_to"
    QUALIFIER_NAME = "GBQualifier_name"
    QUALIFIER_VALUE = "GBQualifier_value"


class INSDCDatabaseCrossReference(BaseModel):
    """A cross reference qualifier.

    Based on the INSDC guidelines for the /db_xref qualifier

    Reference:
    https://www.insdc.org/submitting-standards/dbxref-qualifier-vocabulary/
    """

    database: str

    identifier: int | str

    @classmethod
    def from_string(cls, db_xref: str) -> "INSDCDatabaseCrossReference":
        """Return labelled data from a standard ``db_xref`` string."""
        database, identifier = db_xref.split(":")[:2]

        return INSDCDatabaseCrossReference.model_validate(
            {"database": database, "identifier": identifier}
        )

    @field_validator("identifier", mode="after")
    @classmethod
    def convert_integer_identifier(cls, v: str | int) -> str | int:
        """Check if identifier can be converted to an integer and convert if so."""
        with contextlib.suppress(ValueError):
            return int(v)

        return v


class INSDCMolType(StrEnum):
    """In vivo molecule types.

    Based on the INSDC controlled vocabulary list for the /mol_type qualifier

    Reference:
    https://www.insdc.org/submitting-standards/controlled-vocabulary-moltype-qualifier/
    """

    GENOMIC_DNA = "genomic DNA"
    OTHER_DNA = "other DNA"
    UNASSIGNED_DNA = "unassigned DNA"
    GENOMIC_RNA = "genomic RNA"
    MRNA = "mRNA"
    TRNA = "tRNA"
    TRANSCRIBED_RNA = "transcribed RNA"
    VIRAL_CRNA = "viral cRNA"
    OTHER_RNA = "other RNA"

    @classmethod
    def from_molecule(cls: type, molecule: Molecule) -> "INSDCMolType":
        """Return an INSDC moltype value that matches the passed ref-builder molecule."""

        match molecule.type:
            case MolType.DNA:
                return INSDCMolType.GENOMIC_DNA
            case MolType.RNA:
                return INSDCMolType.GENOMIC_RNA
            case MolType.CRNA:
                return INSDCMolType.VIRAL_CRNA
            case MolType.MRNA:
                return INSDCMolType.MRNA
            case MolType.TRNA:
                return INSDCMolType.TRANSCRIBED_RNA
            case _:
                raise ValueError(
                    f"Molecule cannot be matched to INSDCMolType: {molecule.type}"
                )

    def to_molecule(self) -> MolType:
        """Return ref-builder molecule type."""
        if "DNA" in self.value:
            return MolType.DNA

        match self:
            case INSDCMolType.GENOMIC_RNA:
                return MolType.RNA
            case INSDCMolType.VIRAL_CRNA:
                return MolType.CRNA
            case INSDCMolType.MRNA:
                return MolType.MRNA
            case INSDCMolType.TRANSCRIBED_RNA:
                return MolType.TRNA
            case _:
                raise ValueError(
                    f"INSDC molecule type cannot be matched to MolType: {self}"
                )


class INSDCFeatureQualifier(BaseModel):
    """A feature qualifier in a GenBank source table."""

    model_config = ConfigDict(populate_by_name=True, use_enum_values=True)

    name: Annotated[str, Field(validation_alias=FTGenbankLexicon.QUALIFIER_NAME)]
    """The name of the qualifier."""

    value: Annotated[
        str | None, Field(validation_alias=FTGenbankLexicon.QUALIFIER_VALUE)
    ] = None
    """The value of the qualifier. 
    If this field is not present, the qualifier is a flag.
    """


class INSDCFeatureInterval(BaseModel):
    """An interval in a GenBank source table."""

    model_config = ConfigDict(populate_by_name=True, use_enum_values=True)

    start: Annotated[
        int,
        Field(validation_alias=FTGenbankLexicon.INTERVAL_FROM),
    ]
    """The start of the interval."""

    end: Annotated[
        int,
        Field(validation_alias=FTGenbankLexicon.INTERVAL_TO),
    ]
    """The end of the interval."""


class INSDCFeature(BaseModel):
    """An item in a feature table."""

    model_config = ConfigDict(populate_by_name=True, use_enum_values=True)

    key: Annotated[
        str,
        Field(validation_alias=FTGenbankLexicon.FEATURE_KEY),
    ]
    """The name of the feature."""

    intervals: Annotated[
        list[INSDCFeatureInterval],
        Field(validation_alias=FTGenbankLexicon.FEATURE_INTERVALS),
    ]
    """The intervals of the feature."""

    qualifiers: Annotated[
        list[INSDCFeatureQualifier] | None,
        Field(validation_alias=FTGenbankLexicon.FEATURE_QUALIFIERS),
    ] = None
    """The qualifiers of the feature."""


class INSDCFeatureSource(INSDCFeature):
    """Normalized data from a source feature."""

    key: Annotated[
        Literal["source"],
        Field(validation_alias=FTGenbankLexicon.FEATURE_KEY),
    ]

    qualifiers: Annotated[
        list[INSDCFeatureQualifier],
        Field(validation_alias=FTGenbankLexicon.FEATURE_QUALIFIERS),
    ]


class INSDCFeatureTable(RootModel):
    """An unsorted INSDC feature table."""

    root: list[INSDCFeature]

    def __iter__(self) -> Iterator[INSDCFeature]:
        return iter(self.root)
