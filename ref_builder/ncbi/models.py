"""Models for NCBI Genbank and Taxonomy data."""

from enum import StrEnum
from typing import Annotated, Any

from pydantic import (
    AliasChoices,
    BaseModel,
    ConfigDict,
    Field,
    computed_field,
    field_validator,
    model_validator,
)
from pydantic_core import PydanticCustomError

from ref_builder.insdc.models import (
    INSDCDatabaseCrossReference,
    INSDCFeature,
    INSDCFeatureSource,
    INSDCFeatureTable,
)
from ref_builder.models import Molecule, MolType, Strandedness, Topology


class NCBIDatabase(StrEnum):
    """NCBI databases used in ref-builder."""

    NUCCORE = "nuccore"
    TAXONOMY = "taxonomy"


class NCBIRank(StrEnum):
    """NCBI ranks used in ref-builder."""

    FAMILY = "family"
    ORDER = "order"
    GENUS = "genus"
    SPECIES = "species"
    ISOLATE = "isolate"
    NO_RANK = "no rank"


class NCBISourceMolType(StrEnum):
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
    def from_molecule(cls: type, molecule: Molecule) -> "NCBISourceMolType":
        """Return a NCBI moltype value that matches the passed ref-builder molecule."""
        match molecule.type:
            case MolType.DNA:
                return NCBISourceMolType.GENOMIC_DNA
            case MolType.RNA:
                return NCBISourceMolType.GENOMIC_RNA
            case MolType.CRNA:
                return NCBISourceMolType.VIRAL_CRNA
            case MolType.MRNA:
                return NCBISourceMolType.MRNA
            case MolType.TRNA:
                return NCBISourceMolType.TRANSCRIBED_RNA
            case _:
                raise ValueError(
                    f"Molecule cannot be matched to NCBISourceMolType: {molecule.type}"
                )


class NCBISource(BaseModel):
    """An NCBI source table."""

    model_config = ConfigDict(populate_by_name=True)

    taxid: int = Field(validation_alias=AliasChoices("taxon", "db_xref"))
    organism: str
    mol_type: NCBISourceMolType
    isolate: str = ""
    host: str = ""
    segment: str = ""
    strain: str = ""
    clone: str = ""
    proviral: bool = False
    macronuclear: bool = False
    focus: bool = False
    transgenic: bool = False

    @classmethod
    def from_feature(cls, feature: INSDCFeatureSource | INSDCFeature) -> "NCBISource":
        """Return a NCBISource from a source feature."""
        qualifier_data = {}

        for raw_qualifier in feature.qualifiers:
            if raw_qualifier.name == "db_xref":
                crossref = INSDCDatabaseCrossReference.from_string(raw_qualifier.value)
                qualifier_data[crossref.database] = crossref.identifier

            else:
                qualifier_data[raw_qualifier.name] = (
                    True if raw_qualifier.value is None
                    else raw_qualifier.value
                )

        return NCBISource.model_validate(qualifier_data)

    @field_validator("taxid", mode="before")
    def convert_taxid(cls, v: int | str):
        """Extract TaxId from db_xref string as an integer."""
        if type(v) == int:
            return v

        crossref = INSDCDatabaseCrossReference.from_string(v.value)

        return crossref.identifier

class GenbankRecordKey(StrEnum):
    """Keys used in Genbank XML records."""

    ACCESSION_KEY = "GBSeq_primary-accession"
    ACCESSION = "GBSeq_accession-version"
    COMMENT = "GBSeq_comment"
    DEFINITION = "GBSeq_definition"
    LENGTH = "GBSeq_length"
    FEATURE_TABLE = "GBSeq_feature-table"
    MOLTYPE = "GBSeq_moltype"
    ORGANISM = "GBSeq_organism"
    SEQUENCE = "GBSeq_sequence"
    STRANDEDNESS = "GBSeq_strandedness"
    TOPOLOGY = "GBSeq_topology"


class NCBIGenbank(BaseModel):
    """An NCBI Genbank record."""

    model_config = ConfigDict(populate_by_name=True, use_enum_values=True)

    accession: Annotated[str, Field(validation_alias=GenbankRecordKey.ACCESSION_KEY)]
    accession_version: Annotated[str, Field(validation_alias=GenbankRecordKey.ACCESSION)]
    strandedness: Annotated[
        Strandedness,
        Field(validation_alias=GenbankRecordKey.STRANDEDNESS),
    ]
    moltype: Annotated[MolType, Field(validation_alias=GenbankRecordKey.MOLTYPE)]
    topology: Annotated[Topology, Field(validation_alias=GenbankRecordKey.TOPOLOGY)]
    definition: Annotated[str, Field(validation_alias=GenbankRecordKey.DEFINITION)]
    organism: Annotated[str, Field(validation_alias=GenbankRecordKey.ORGANISM)]
    sequence: Annotated[
        str,
        Field(
            validation_alias=GenbankRecordKey.SEQUENCE,
            pattern=r"^[ATCGRYKMSWBDHVNatcgrykmswbdhvn]+$",
        ),
    ]
    source: Annotated[
        NCBISource,
        Field(validation_alias=GenbankRecordKey.FEATURE_TABLE)
    ]
    comment: Annotated[str, Field(validation_alias=GenbankRecordKey.COMMENT)] = ""
    features: Annotated[
        INSDCFeatureTable,
        Field(
            default_factory=list,
            exclude=True,
            validation_alias=GenbankRecordKey.FEATURE_TABLE,
        )
    ]


    @computed_field()
    def refseq(self) -> bool:
        """Whether this is a RefSeq record.

        RefSeq records have accessions that start with "NC_".
        """
        return self.accession.startswith("NC_")

    @field_validator("sequence", mode="after")
    @classmethod
    def uppercase_sequence(cls, raw: str) -> str:
        """Force the sequence field to uppercase."""
        return raw.upper()

    @field_validator("source", mode="before")
    @classmethod
    def extract_source_from_feature_table(
        cls,
        raw: NCBISource | list[dict[str:Any]],
    ) -> NCBISource:
        """If the source field isn't a ``NCBISource`` object, extract the data from
        the feature table and convert.
        """
        if isinstance(raw, NCBISource):
            return raw

        for raw_feature in raw:
            feature = INSDCFeature.model_validate(raw_feature)

            if feature.key == "source":
                return NCBISource.from_feature(feature)

        raise ValueError("Feature table contains no ``source`` table.")

    @model_validator(mode="after")
    def check_source(self) -> "NCBIGenbank":
        """Check that the source organism matches the record organism."""
        if self.source.organism == self.organism:
            return self

        raise ValueError("Non-matching organism fields on record and source")


class NCBILineage(BaseModel):
    """An NCBI lineage record."""

    model_config = ConfigDict(populate_by_name=True)

    id: Annotated[int, Field(validation_alias="TaxId")]
    name: Annotated[str, Field(validation_alias="ScientificName")]
    rank: Annotated[str, Field(validation_alias="Rank")]


class NCBITaxonomyOtherNames(BaseModel):
    """An NCBI taxonomy record's other names."""

    acronym: Annotated[list[str], Field(validation_alias="Acronym")] = []
    genbank_acronym: Annotated[list[str], Field(validation_alias="GenbankAcronym")] = []
    equivalent_name: Annotated[list[str], Field(validation_alias="EquivalentName")] = []
    synonym: Annotated[list[str], Field(validation_alias="Synonym")] = []
    includes: Annotated[list[str], Field(validation_alias="Includes")] = []


class NCBITaxonomy(BaseModel):
    """An NCBI taxonomy record."""

    model_config = ConfigDict(populate_by_name=True)

    id: Annotated[int, Field(validation_alias="TaxId")]
    name: Annotated[str, Field(validation_alias="ScientificName")]
    other_names: Annotated[
        NCBITaxonomyOtherNames,
        Field(validation_alias="OtherNames"),
    ] = NCBITaxonomyOtherNames()
    lineage: Annotated[list[NCBILineage], Field(validation_alias="LineageEx")]
    rank: Annotated[NCBIRank, Field(validation_alias=AliasChoices("rank", "Rank"))]

    @field_validator("rank", mode="before")
    @classmethod
    def check_rank_level(cls, v: str | NCBIRank) -> str | NCBIRank:
        if v not in (NCBIRank.SPECIES, NCBIRank.NO_RANK, NCBIRank.ISOLATE):
            raise PydanticCustomError(
                "taxon_rank_too_high",
                "Taxon rank ({rank}) is too high.",
                {"rank": v},
            )

        return v

    @computed_field
    def species(self) -> NCBILineage:
        """Return the species level taxon in the lineage."""

        if self.rank is NCBIRank.SPECIES:
            return NCBILineage(id=self.id, name=self.name, rank=self.rank)

        for item in self.lineage:
            if item.rank == "species":
                return item

        raise PydanticCustomError(
            "irrelevant_taxon_lineage",
            "No species level taxon found in lineage."
            + "Taxonomy Id {taxid} level ({rank}) is too high.",
            {"taxid": self.id, "rank": self.rank},
        )
