"""Factories for generating quasi-realistic NCBISource and NCBIGenbank data."""

from typing import Literal

from faker.providers import lorem
from polyfactory import Use
from polyfactory.decorators import post_generated
from polyfactory.factories.pydantic_factory import ModelFactory

from ref_builder.models.molecule import Molecule, MoleculeType
from ref_builder.models.plan import Plan
from ref_builder.ncbi.models import (
    NCBIGenbank,
    NCBILineage,
    NCBIRank,
    NCBISource,
    NCBISourceMolType,
    NCBITaxonomy,
    NCBITaxonomyOtherNames,
)
from ref_builder.plan import get_segments_max_length, get_segments_min_length
from tests.fixtures.providers import (
    AccessionProvider,
    OrganismProvider,
    SegmentProvider,
    SequenceProvider,
)

DNA_MOLTYPES = {
    NCBISourceMolType.GENOMIC_DNA,
    NCBISourceMolType.OTHER_DNA,
    NCBISourceMolType.UNASSIGNED_DNA,
}
"""NCBISourceMolTypes that map to MolType.DNA"""


ModelFactory.__faker__.add_provider(AccessionProvider)
ModelFactory.__faker__.add_provider(OrganismProvider)
ModelFactory.__faker__.add_provider(SegmentProvider)
ModelFactory.__faker__.add_provider(SequenceProvider)
ModelFactory.__faker__.add_provider(lorem)


class NCBISourceFactory(ModelFactory[NCBISource]):
    """NCBISource Factory with quasi-realistic data."""

    host = Use(ModelFactory.__faker__.host)
    """A fake host name for the virus."""

    @classmethod
    def isolate(cls) -> str:
        """Raw isolate name faker."""
        if cls.__faker__.boolean(80):
            delimiter = cls.__faker__.random_element(["", "-", "_"])
            components = cls.__faker__.random_elements(
                [
                    cls.__faker__.country().replace(" ", ""),
                    cls.__faker__.last_name(),
                    cls.__faker__.country_code(),
                    str(cls.__faker__.random_int(0, 9)),
                    str(cls.__faker__.random_int(10, 99)),
                    str(cls.__faker__.random_int(100, 9999)),
                ],
                2,
                unique=True,
            )

            return delimiter.join(components)

        return ""

    organism = Use(ModelFactory.__faker__.organism)
    """A fake name for the virus."""

    taxid = Use(ModelFactory.__faker__.random_int, min=1000, max=999999)

    @classmethod
    def segment(cls) -> str:
        """Raw segment name faker."""
        if cls.__faker__.boolean(80):
            return (
                f"{cls.__faker__.segment_prefix()}"
                f"{cls.__faker__.segment_delimiter()}"
                f"{cls.__faker__.segment_key()}"
            )
        return cls.__faker__.segment_key()

    @classmethod
    def clone(cls) -> str:
        """Raw clone name faker."""
        if cls.__faker__.boolean(10):
            delimiter = cls.__faker__.random_element(["-", "_", " ", "/"])
            return delimiter.join(cls.__faker__.words(2))

        return ""

    @classmethod
    def strain(cls) -> str:
        """Raw strain name faker."""
        if cls.__faker__.boolean(10):
            delimiter = cls.__faker__.random_element(["-", "_", " ", "/"])
            return delimiter.join(cls.__faker__.words(2))

        return ""

    proviral = Use(ModelFactory.__faker__.boolean, chance_of_getting_true=5)
    """Pseudorandom proviral flag for DNA records only."""

    @post_generated
    @classmethod
    def macronuclear(cls, mol_type: NCBISourceMolType) -> bool:
        """Pseudorandom macronuclear flag for DNA records only."""
        if mol_type in DNA_MOLTYPES:
            return cls.__faker__.boolean(5)

        return False

    @post_generated
    @classmethod
    def focus(cls, mol_type: NCBISourceMolType) -> bool:
        """Pseudorandom focus flag for DNA records only.
        Mutually exclusive with transgenic.
        """
        if mol_type in DNA_MOLTYPES:
            return cls.__faker__.boolean(5)

        return False

    @post_generated
    @classmethod
    def transgenic(cls, focus: bool) -> bool:
        """Transgenic flag set to False if focus is True."""
        if focus and cls.__faker__.boolean(5):
            return not focus

        return False


class NCBIGenbankFactory(ModelFactory[NCBIGenbank]):
    """NCBIGenbank Factory with quasi-realistic data."""

    source = Use(NCBISourceFactory.build)

    accession = Use(ModelFactory.__faker__.accession)
    """A fake accession that conforms to NCBI standards."""

    @post_generated
    @classmethod
    def accession_version(cls, accession: str) -> str:
        """Raw accession_version faker."""
        return f"{accession}.{cls.__faker__.random_int(1, 3)}"

    @post_generated
    @classmethod
    def moltype(cls, source: NCBISource) -> MoleculeType:
        """Map moltype field to source.moltype equivalent."""
        if source.mol_type in DNA_MOLTYPES:
            return MoleculeType.DNA

        try:
            return {
                NCBISourceMolType.GENOMIC_RNA: MoleculeType.RNA,
                NCBISourceMolType.MRNA: MoleculeType.MRNA,
                NCBISourceMolType.TRANSCRIBED_RNA: MoleculeType.RNA,
                NCBISourceMolType.VIRAL_CRNA: MoleculeType.CRNA,
                NCBISourceMolType.TRNA: MoleculeType.TRNA,
                NCBISourceMolType.OTHER_RNA: MoleculeType.RNA,
                NCBISourceMolType.UNASSIGNED_RNA: MoleculeType.RNA,
            }[source.mol_type]
        except KeyError as err:
            raise ValueError(
                f"Source moltype {source.mol_type} cannot be matched to MolType",
            ) from err

    @post_generated
    @classmethod
    def organism(cls, source: NCBISource) -> str:
        """Organism faker."""
        return source.organism

    @classmethod
    def sequence(cls) -> str:
        """Sequence faker."""
        return cls.__faker__.sequence()

    @post_generated
    @classmethod
    def taxid(cls, source: NCBISource) -> int:
        """Match taxid field to source.taxid."""
        return source.taxid

    @classmethod
    def build_isolate(
        cls,
        segment_count: int,
        refseq: bool = True,
        segment_prefix: str | None = None,
        segment_key_style: Literal["letter", "number"] = "letter",
        segment_delimiter: str | None = None,
        version: int | None = None,
        base_source: NCBISource | None = None,
        **kwargs,
    ) -> list[NCBIGenbank]:
        """Build a list of NCBIGenbank records representing one virus isolate.

        All records will share the same organism, taxid, and other metadata,
        with sequential accessions and consistent segment naming.

        Args:
            segment_count: Number of segments to generate
            refseq: Use RefSeq accessions (NC_*) if True, GenBank otherwise
            segment_prefix: Prefix for segment names (e.g., "RNA", "DNA")
                          Auto-detected from mol_type if None
            segment_key_style: "letter" for A,B,C... or "number" for 1,2,3...
            segment_delimiter: Delimiter between prefix and key
                             Random choice if None
            version: Version number for all records (random 1-3 if None)
            base_source: Pre-built NCBISource to use for all records
            **kwargs: Additional overrides for NCBIGenbank.build()

        Returns:
            List of NCBIGenbank records representing one isolate

        Raises:
            ValueError: If segment_count > 26 for letter style

        """
        if segment_key_style == "letter" and segment_count > 26:
            raise ValueError("segment_count cannot exceed 26 for letter style")

        if refseq:
            accessions = cls.__faker__.refseq_accessions(segment_count)
        else:
            accessions = cls.__faker__.genbank_accessions(segment_count)

        if version is None:
            version = cls.__faker__.random_int(1, 3)

        if base_source is None:
            base_source = NCBISourceFactory.build()

        if segment_prefix is None:
            segment_prefix = "DNA" if base_source.mol_type in DNA_MOLTYPES else "RNA"

        if segment_delimiter is None:
            segment_delimiter = cls.__faker__.random_element([" ", "-", "_"])

        if segment_key_style == "letter":
            segment_keys = [chr(ord("A") + i) for i in range(segment_count)]
        else:
            segment_keys = [str(i + 1) for i in range(segment_count)]

        records = []

        for i in range(segment_count):
            segment_name = f"{segment_prefix}{segment_delimiter}{segment_keys[i]}"

            source = NCBISourceFactory.build(
                organism=base_source.organism,
                taxid=base_source.taxid,
                mol_type=base_source.mol_type,
                host=base_source.host,
                isolate=base_source.isolate,
                strain=base_source.strain,
                clone=base_source.clone,
                proviral=base_source.proviral,
                macronuclear=base_source.macronuclear,
                focus=base_source.focus,
                transgenic=base_source.transgenic,
                segment=segment_name,
            )

            record = cls.build(
                accession=accessions[i],
                accession_version=f"{accessions[i]}.{version}",
                source=source,
                **kwargs,
            )

            records.append(record)

        return records

    @classmethod
    def build_from_plan(
        cls,
        plan: Plan,
        refseq: bool = True,
        version: int | None = None,
        base_source: NCBISource | None = None,
        organism: str | None = None,
        taxid: int | None = None,
        molecule: Molecule | None = None,
        **kwargs,
    ) -> list[NCBIGenbank]:
        """Build NCBIGenbank records satisfying an OTU Plan.

        All records will share the same organism, taxid, and other metadata,
        with sequential accessions. Segment names and sequence lengths match the plan.

        Args:
            plan: The Plan object defining segment requirements
            refseq: Use RefSeq accessions (NC_*) if True, GenBank otherwise
            version: Version number for all records (random 1-3 if None)
            base_source: Pre-built NCBISource to use for all records
            organism: Override organism name (uses base_source if not specified)
            taxid: Override taxid (uses base_source if not specified)
            molecule: Molecule to derive mol_type from (uses base_source if not specified)
            **kwargs: Additional overrides for NCBIGenbank.build()

        Returns:
            List of NCBIGenbank records matching the plan

        """
        if base_source is None:
            base_source = NCBISourceFactory.build()

        if organism is not None:
            base_source = NCBISourceFactory.build(
                organism=organism,
                taxid=base_source.taxid,
                mol_type=base_source.mol_type,
                host=base_source.host,
                isolate=base_source.isolate,
                strain=base_source.strain,
                clone=base_source.clone,
                proviral=base_source.proviral,
                macronuclear=base_source.macronuclear,
                focus=base_source.focus,
                transgenic=base_source.transgenic,
                segment=base_source.segment,
            )

        if taxid is not None:
            base_source = NCBISourceFactory.build(
                organism=base_source.organism,
                taxid=taxid,
                mol_type=base_source.mol_type,
                host=base_source.host,
                isolate=base_source.isolate,
                strain=base_source.strain,
                clone=base_source.clone,
                proviral=base_source.proviral,
                macronuclear=base_source.macronuclear,
                focus=base_source.focus,
                transgenic=base_source.transgenic,
                segment=base_source.segment,
            )

        if molecule is not None:
            base_source = NCBISourceFactory.build(
                organism=base_source.organism,
                taxid=base_source.taxid,
                mol_type=NCBISourceMolType.from_molecule(molecule),
                host=base_source.host,
                isolate=base_source.isolate,
                strain=base_source.strain,
                clone=base_source.clone,
                proviral=base_source.proviral,
                macronuclear=base_source.macronuclear,
                focus=base_source.focus,
                transgenic=base_source.transgenic,
                segment=base_source.segment,
            )

        if version is None:
            version = cls.__faker__.random_int(1, 3)

        if refseq:
            accessions = cls.__faker__.refseq_accessions(len(plan.segments))
        else:
            accessions = cls.__faker__.genbank_accessions(len(plan.segments))

        records = []
        for i, segment in enumerate(plan.segments):
            min_length = get_segments_min_length([segment])
            max_length = get_segments_max_length([segment])
            sequence = cls.__faker__.sequence(min=min_length, max=max_length)

            segment_name = str(segment.name) if segment.name is not None else ""

            source = NCBISourceFactory.build(
                organism=base_source.organism,
                taxid=base_source.taxid,
                mol_type=base_source.mol_type,
                host=base_source.host,
                isolate=base_source.isolate,
                strain=base_source.strain,
                clone=base_source.clone,
                proviral=base_source.proviral,
                macronuclear=base_source.macronuclear,
                focus=base_source.focus,
                transgenic=base_source.transgenic,
                segment=segment_name,
            )

            record = cls.build(
                accession=accessions[i],
                accession_version=f"{accessions[i]}.{version}",
                source=source,
                sequence=sequence,
                **kwargs,
            )

            records.append(record)

        return records


class NCBILineageFactory(ModelFactory[NCBILineage]):
    """NCBILineage Factory with quasi-realistic data."""

    id = Use(ModelFactory.__faker__.random_int, min=1000, max=999999)
    """A realistic taxonomy ID."""

    name = Use(ModelFactory.__faker__.organism)
    """Generate a realistic organism name."""

    rank = Use(
        ModelFactory.__faker__.random_element,
        elements=[NCBIRank.FAMILY, NCBIRank.ORDER, NCBIRank.GENUS, NCBIRank.SPECIES],
    )
    """A realistic NCBI rank."""


class NCBITaxonomyOtherNamesFactory(ModelFactory[NCBITaxonomyOtherNames]):
    """NCBITaxonomyOtherNames Factory with quasi-realistic data."""

    @classmethod
    def acronym(cls) -> list[str]:
        """Generate a realistic acronym list."""
        if cls.__faker__.boolean(30):
            return [
                "".join(
                    [word[0].upper() for word in cls.__faker__.organism().split(" ")]
                )
            ]
        return []

    genbank_acronym = Use(list)
    """Empty genbank acronym list."""

    equivalent_name = Use(list)
    """Empty equivalent name list."""

    synonym = Use(list)
    """Empty synonym list."""

    includes = Use(list)
    """Empty includes list."""


class NCBITaxonomyFactory(ModelFactory[NCBITaxonomy]):
    """NCBITaxonomy Factory with quasi-realistic data."""

    id = Use(ModelFactory.__faker__.random_int, min=1000, max=999999)
    """A realistic taxonomy ID."""

    name = Use(ModelFactory.__faker__.organism)
    """Generate a realistic organism name."""

    other_names = Use(NCBITaxonomyOtherNamesFactory.build)
    """Generate realistic other names."""

    rank = Use(
        ModelFactory.__faker__.random_element,
        elements=[NCBIRank.SPECIES, NCBIRank.ISOLATE, NCBIRank.NO_RANK],
    )
    """A realistic NCBI rank (limited to species or below)."""

    @post_generated
    @classmethod
    def lineage(cls, id: int, name: str, rank: NCBIRank) -> list[NCBILineage]:
        """Generate a realistic lineage with species-level taxon."""
        lineage = [
            NCBILineageFactory.build(rank=NCBIRank.FAMILY),
            NCBILineageFactory.build(rank=NCBIRank.ORDER),
            NCBILineageFactory.build(rank=NCBIRank.GENUS),
        ]

        if rank == NCBIRank.SPECIES:
            return lineage
        species = NCBILineageFactory.build(rank=NCBIRank.SPECIES)
        lineage.append(species)
        return lineage
