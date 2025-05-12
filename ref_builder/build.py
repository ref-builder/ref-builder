import datetime
from pathlib import Path
from typing import Annotated

import arrow
import orjson
from pydantic import (
    BaseModel,
    ConfigDict,
    Field,
    RootModel,
    ValidationError,
    field_serializer,
    field_validator,
)

from ref_builder.models import Molecule, Strandedness
from ref_builder.otu.validate import get_validated_otu
from ref_builder.otu.validators.isolate import Isolate
from ref_builder.otu.validators.otu import OTU
from ref_builder.otu.validators.sequence import Sequence
from ref_builder.plan import Plan, Segment, SegmentName, SegmentRule
from ref_builder.repo import Repo


def _get_molecule_string(molecule: Molecule) -> str:
    """Return a string representation of the molecule."""
    string = ""

    if molecule.strandedness == Strandedness.SINGLE:
        string += "ss"
    else:
        string += "ds"

    if "DNA" in molecule.type:
        string += "DNA"
    else:
        string += "RNA"

    return string


class ProductionSegment(BaseModel):
    """A segment belonging to a production OTU schema."""

    molecule: str
    name: str
    required: bool

    @field_validator("molecule", mode="before")
    @classmethod
    def convert_from_molecule_object(cls, v: str | Molecule) -> str:
        """Convert molecule field from object to string."""
        if isinstance(v, Molecule):
            return _get_molecule_string(v)

        if type(v) is str:
            return v

        raise ValueError(f"Invalid molecule field input: {v}, {type(v)}")

    @field_validator("name", mode="before")
    @classmethod
    def convert_from_segment_name(cls, v: str | SegmentName | None) -> str:
        """Convert segment name from either a SegmentName or a null."""
        if v is None:
            return "Unnamed"

        if isinstance(v, SegmentName):
            return str(v)

        if isinstance(v, str):
            return v

        raise ValueError(f"Invalid segment name input: {v}, {type(v)}")

    @field_validator("required", mode="before")
    @classmethod
    def convert_required_field_to_boolean(cls, v: bool | SegmentRule) -> bool:
        """Convert required field from SegmentRule to boolean."""
        if isinstance(v, SegmentRule | str):
            return v == SegmentRule.REQUIRED

        if isinstance(v, bool):
            return v

        raise ValueError(f"Invalid required field input: {v}, {type(v)}")

    @classmethod
    def build_from_segment_and_molecule(
        cls, segment: Segment, molecule: Molecule
    ) -> "ProductionSegment":
        """Return a production segment from Segment and Molecule data."""
        return ProductionSegment.model_validate(
            {
                "molecule": molecule,
                "name": segment.name,
                "required": segment.rule,
            }
        )


class ProductionSchema(RootModel):
    """A ``schema`` plan belonging to a production OTU."""

    root: list[ProductionSegment]

    @classmethod
    def build_from_plan_and_molecule(cls, plan: Plan, molecule: Molecule):
        molecule_string = ""

        if molecule.strandedness == Strandedness.SINGLE:
            molecule_string += "ss"
        else:
            molecule_string += "ds"

        if "DNA" in molecule.type:
            molecule_string += "DNA"
        else:
            molecule_string += "RNA"

        return ProductionSchema(
            [
                ProductionSegment(
                    molecule=molecule_string,
                    name=str(segment.name) if segment.name else "Unnamed",
                    required=segment.rule == SegmentRule.REQUIRED,
                )
                for segment in plan.segments
            ]
        )


class ProductionResource(BaseModel):
    """Parent class for production-ready resources."""

    model_config = ConfigDict(serialize_by_alias=True)

    id: str


class ProductionSequence(ProductionResource):
    """A production sequence belonging to a production isolate."""

    id: Annotated[str, Field(serialization_alias="_id")]
    accession: str
    definition: str
    sequence: str
    segment: str = "Unnamed"
    host: str = ""

    @classmethod
    def build_from_validated_sequence_and_plan(
        cls,
        validated_sequence: Sequence,
        otu_plan: Plan,
    ) -> "ProductionSequence":
        segment = otu_plan.get_segment_by_id(validated_sequence.segment)

        return ProductionSequence(
            id=str(validated_sequence.id),
            accession=str(validated_sequence.accession),
            definition=validated_sequence.definition,
            host="",
            segment="Unnamed" if segment.name is None else str(segment.name),
            sequence=validated_sequence.sequence,
        )


class ProductionIsolate(ProductionResource):
    """A production isolate belonging to a production OTU."""

    id: Annotated[str, Field(serialization_alias="id")]
    default: bool
    sequences: list[ProductionSequence]
    source_name: str = "unknown"
    source_type: str = "unknown"

    @classmethod
    def build_from_validated_isolate(
        cls, validated_isolate: Isolate, plan: Plan, is_representative: bool
    ) -> "ProductionIsolate":
        return ProductionIsolate(
            id=str(validated_isolate.id),
            default=is_representative,
            sequences=[
                ProductionSequence.build_from_validated_sequence_and_plan(
                    validated_sequence, plan
                )
                for validated_sequence in validated_isolate.sequences
            ],
            source_name=(
                validated_isolate.name.value
                if validated_isolate.name is not None
                else "unknown"
            ),
            source_type=(
                str(validated_isolate.name.type)
                if validated_isolate.name is not None
                else "unknown"
            ),
        )


class ProductionOTU(ProductionResource):
    """A production OTU belonging to a production reference file."""

    id: Annotated[str, Field(serialization_alias="_id")]
    abbreviation: str
    isolates: list[ProductionIsolate]
    name: str
    plan: Annotated[ProductionSchema, Field(serialization_alias="schema")]
    taxid: int

    @classmethod
    def build_from_validated_otu(cls, validated_otu: OTU):
        return ProductionOTU(
            id=str(validated_otu.id),
            abbreviation=validated_otu.acronym,
            isolates=[
                ProductionIsolate.build_from_validated_isolate(
                    validated_isolate,
                    is_representative=(
                        validated_isolate.id == validated_otu.representative_isolate
                    ),
                    plan=validated_otu.plan,
                )
                for validated_isolate in validated_otu.isolates
            ],
            name=validated_otu.name,
            plan=ProductionSchema.build_from_plan_and_molecule(
                validated_otu.plan,
                molecule=validated_otu.molecule,
            ),
            taxid=validated_otu.taxid,
        )

    @field_validator("isolates", mode="after")
    @classmethod
    def sort_isolates(cls, v: list[ProductionIsolate]) -> list[ProductionIsolate]:
        """Sort isolates by named alphabetical order."""
        v.sort(key=lambda x: x.source_type + x.source_name)

        return v


class ProductionReference(BaseModel):
    """A production-ready reference."""

    created_at: datetime.datetime
    data_type: str
    name: str
    organism: str
    otus: list[ProductionOTU]

    @field_serializer("created_at")
    def convert_datetime_to_iso(self, v: datetime.datetime) -> str:
        """Convert ``created_at`` datetime format to an ISO 8601 formatted string."""
        return v.isoformat()

    @field_validator("otus", mode="after")
    @classmethod
    def sort_otus(cls, v: list[ProductionOTU]) -> list[ProductionOTU]:
        """Sort OTUs by named alphabetical order."""
        v.sort(key=lambda x: x.name)

        return v


def build_json(indent: bool, output_path: Path, path: Path, version: str) -> Path:
    """Build a Virtool reference JSON file from a data directory.

    :param indent: whether to indent the JSON output
    :param output_path: The path to write the output JSON file to
    :param path: The path to a reference repository
    :param version: the version string to include in the reference.json file
    """
    repo = Repo(path)

    otus = []

    for unvalidated_otu in repo.iter_otus():
        try:
            validated_otu = get_validated_otu(unvalidated_otu)
        except ValidationError:
            continue

        otus.append(ProductionOTU.build_from_validated_otu(validated_otu))

    production_reference = ProductionReference(
        created_at=arrow.utcnow().datetime,
        data_type=repo.meta.data_type,
        name=version,
        organism=repo.meta.organism,
        otus=otus,
    )

    with open(output_path, "wb") as f:
        f.write(
            orjson.dumps(
                production_reference.model_dump(mode="json"),
                option=orjson.OPT_INDENT_2 if indent else 0,
            ),
        )

    if not output_path.exists():
        raise FileNotFoundError(f"Built reference not found at {output_path}.")

    return output_path
