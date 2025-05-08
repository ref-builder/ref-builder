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
    field_validator,
)

from ref_builder.models import Molecule, Strandedness
from ref_builder.otu.validate import get_validated_otu
from ref_builder.otu.validators.isolate import Isolate
from ref_builder.otu.validators.otu import OTU
from ref_builder.otu.validators.sequence import Sequence
from ref_builder.plan import Plan, SegmentRule
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

        return ProductionSchema([
            ProductionSegment(
                molecule=molecule_string,
                name=str(segment.name) if segment.name else "Unnamed",
                required=segment.rule == SegmentRule.REQUIRED,
            ) for segment in plan.segments
        ])


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
        cls, validated_sequence: Sequence, otu_plan: Plan,
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
                ProductionSequence.build_from_validated_sequence_and_plan(validated_sequence, plan)
                for validated_sequence in validated_isolate.sequences
            ],
            source_name=(
                validated_isolate.name.value if validated_isolate.name is not None else "unknown"
            ),
            source_type=(
                str(validated_isolate.name.type)
                if validated_isolate.name is not None
                else "unknown"
            )
        )


class ProductionOTU(ProductionResource):
    """A production OTU belonging to a production reference file."""

    id: Annotated[str, Field(serialization_alias="_id")]
    abbreviation: str
    isolates: list[ProductionIsolate]
    name: str
    plan: Annotated[
        ProductionSchema,
        Field(serialization_alias="schema")
    ]
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


def build_json(indent: bool, output_path: Path, path: Path, version: str) -> None:
    """Build a Virtool reference JSON file from a data directory.

    :param indent: whether to indent the JSON output
    :param output_path: The path to write the output JSON file to
    :param path: The path to a reference repository
    :param version: the version string to include in the reference.json file
    """
    repo = Repo(path)

    otus = []

    for otu in repo.iter_otus():
        isolates = []

        for isolate in otu.isolates:
            sequences = []

            for sequence in isolate.sequences:
                segment_name = otu.plan.get_segment_by_id(sequence.segment).name

                sequences.append(
                    {
                        "_id": sequence.legacy_id or str(sequence.id),
                        "accession": str(sequence.accession),
                        "definition": sequence.definition,
                        "host": "",
                        "segment": "Unnamed"
                        if segment_name is None
                        else str(segment_name),
                        "sequence": sequence.sequence,
                    }
                )

            isolates.append(
                {
                    "id": isolate.legacy_id or str(isolate.id),
                    "default": isolate.id == otu.representative_isolate,
                    "sequences": sequences,
                    "source_name": isolate.name.value
                    if isolate.name is not None
                    else "unknown",
                    "source_type": str(isolate.name.type)
                    if isolate.name is not None
                    else "unknown",
                },
            )

        isolates.sort(key=lambda x: x["source_type"] + x["source_name"])

        molecule = _get_molecule_string(otu.molecule)

        schema = [
            {
                "molecule": molecule,
                "name": str(segment.name) if segment.name else "Unnamed",
                "required": segment.rule == SegmentRule.REQUIRED,
            }
            for segment in otu.plan.segments
        ]

        otus.append(
            {
                "_id": otu.legacy_id or otu.id,
                "abbreviation": otu.acronym,
                "isolates": isolates,
                "name": otu.name,
                "schema": schema,
                "taxid": otu.taxid,
            },
        )

    otus.sort(key=lambda x: x["name"])

    with open(output_path, "wb") as f:
        f.write(
            orjson.dumps(
                {
                    "created_at": arrow.utcnow().isoformat(),
                    "data_type": repo.meta.data_type,
                    "name": version,
                    "organism": repo.meta.organism,
                    "otus": otus,
                },
                option=orjson.OPT_INDENT_2 if indent else 0,
            ),
        )
