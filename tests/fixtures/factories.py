"""Factories for generating quasi-realistic NCBISource and NCBIGenbank data."""

from typing import Literal

from faker.providers import lorem
from polyfactory import PostGenerated, Use
from polyfactory.decorators import post_generated
from polyfactory.factories.pydantic_factory import ModelFactory

from ref_builder.models import Molecule, MolType, OTUMinimal
from ref_builder.ncbi.models import NCBIGenbank, NCBISource, NCBISourceMolType
from ref_builder.otu.utils import get_segments_max_length, get_segments_min_length
from ref_builder.otu.validators.isolate import IsolateBase
from ref_builder.otu.validators.otu import OTUBase
from ref_builder.otu.validators.sequence import SequenceBase
from ref_builder.plan import (
    Plan,
    Segment,
    SegmentName,
    SegmentRule,
)
from ref_builder.utils import Accession, IsolateName, IsolateNameType
from tests.fixtures.providers import (
    AccessionProvider,
    BusinessProvider,
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
ModelFactory.__faker__.add_provider(BusinessProvider)
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
    def moltype(cls, source: NCBISource) -> MolType:
        """Map moltype field to source.moltype equivalent."""
        if source.mol_type in DNA_MOLTYPES:
            return MolType.DNA

        try:
            return {
                NCBISourceMolType.GENOMIC_RNA: MolType.RNA,
                NCBISourceMolType.MRNA: MolType.MRNA,
                NCBISourceMolType.TRANSCRIBED_RNA: MolType.RNA,
                NCBISourceMolType.VIRAL_CRNA: MolType.CRNA,
                NCBISourceMolType.TRNA: MolType.TRNA,
                NCBISourceMolType.OTHER_RNA: MolType.RNA,
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


def derive_acronym(_: str, values: dict[str, str]) -> str:
    """Derive an acronym from an OTU name."""
    name = values["name"]
    return "".join([part[0].upper() for part in name.split(" ")])


class SegmentFactory(ModelFactory[Segment]):
    """Segment Factory with quasi-realistic data."""

    ModelFactory.__faker__.add_provider(SequenceProvider)

    length = Use(ModelFactory.__faker__.sequence_length)
    """Generate a quasi-realistic length for a sequence."""

    @classmethod
    def name(cls) -> SegmentName | None:
        """Generate a quasi-realistic segment name or null."""
        if cls.__faker__.random_int(0, 10) > 5:
            return SegmentName(
                prefix=cls.__faker__.random_element(["DNA", "RNA"]),
                key=cls.__faker__.segment_key(),
            )

        return None

    @classmethod
    def length_tolerance(cls) -> float:
        """Generate a realistic length tolerance."""
        return cls.__faker__.random_int(0, 5) * 0.01

    @staticmethod
    def build_series(n: int) -> list[Segment]:
        """Generate a matching series of segments"""
        segment_name_keys = "ABCDEF"

        if n > len(segment_name_keys):
            raise ValueError("Generation of over 6 segments is unsupported.")

        return [
            SegmentFactory.build(
                name=SegmentName(prefix="DNA", key=segment_name_keys[i]),
                required=SegmentRule.REQUIRED,
            )
            for i in range(n)
        ]


class PlanFactory(ModelFactory[Plan]):
    """A Polyfactory that generates valid instances of :class:`Plan`.

    The factory generates a random number of segments, with a 75% chance of generating a
    monopartite plan.
    """

    @classmethod
    def segments(cls) -> list[Segment]:
        """Return a set of quasi-realistic segments."""
        # The segment represents a monopartite OTU 75% of the time.
        if cls.__faker__.random_int(0, 3):
            return [SegmentFactory.build(name=None, rule=SegmentRule.REQUIRED)]

        return SegmentFactory.build_series(cls.__faker__.random_int(2, 5))


class SequenceFactory(ModelFactory[SequenceBase]):
    """Sequence factory with quasi-realistic data."""

    id = Use(ModelFactory.__faker__.uuid4)
    """Generate a UUID."""

    definition = Use(ModelFactory.__faker__.sentence)
    """Generate a mock sentence to serve as the definition field."""

    legacy_id = Use(ModelFactory.__faker__.legacy_id)
    """Generate an 8-character unique identifier as used in virtool-cli."""

    segment = Use(ModelFactory.__faker__.uuid4)
    """Generate a quasi-realistic mock segment string."""

    sequence = Use(ModelFactory.__faker__.sequence)
    """Generate a quasi-realistic mock genomic sequence."""

    @classmethod
    def accession(cls) -> Accession:
        """Generate a quasi-realistic accession."""
        ModelFactory.__faker__.add_provider(AccessionProvider)
        return Accession(key=ModelFactory.__faker__.accession(), version=1)

    @classmethod
    def build_on_segment(
        cls, segment: Segment, accession: Accession | None = None
    ) -> SequenceBase:
        """Build a sequence based on a given segment. Takes an optional accession."""
        min_length = get_segments_min_length([segment])
        max_length = get_segments_max_length([segment])

        if accession is None:
            accession = Accession(key=ModelFactory.__faker__.accession(), version=1)

        return SequenceFactory.build(
            accession=accession,
            sequence=cls.__faker__.sequence(min=min_length, max=max_length),
            segment=segment.id,
        )


class IsolateFactory(ModelFactory[IsolateBase]):
    """Isolate factory with quasi-realistic data."""

    ModelFactory.__faker__.add_provider(BusinessProvider)

    id = Use(ModelFactory.__faker__.uuid4, cast_to=None)
    """Generate a UUID."""

    legacy_id = Use(ModelFactory.__faker__.legacy_id)
    """Generate an 8-character unique identifier as used in virtool-cli."""

    @classmethod
    def name(cls) -> IsolateName:
        """Generate a quasi-realistic isolate name."""
        return IsolateName(
            type=IsolateNameType.ISOLATE,
            value=cls.__faker__.word(part_of_speech="noun").capitalize(),
        )

    @classmethod
    def sequences(cls) -> list[SequenceBase]:
        """Generate between 1 and 6 sequences with numerically sequential accessions."""
        sequence_count = cls.__faker__.random_int(1, 6)

        return [
            SequenceFactory.build(accession=Accession(key=accession, version=1))
            for accession in cls.__faker__.accessions(sequence_count)
        ]

    @classmethod
    def build_on_plan(cls, plan: Plan, refseq: bool = False):
        """Take a plan and return a matching isolate."""
        if refseq:
            sequential_accessions = cls.__faker__.refseq_accessions(len(plan.segments))
        else:
            sequential_accessions = cls.__faker__.genbank_accessions(len(plan.segments))

        return IsolateFactory.build(
            sequences=[
                SequenceFactory.build_on_segment(
                    segment=plan.segments[counter],
                    accession=Accession(sequential_accessions[counter], 1),
                )
                for counter in range(len(plan.segments))
            ]
        )


class OTUFactory(ModelFactory[OTUBase]):
    """OTU Factory with quasi-realistic data."""

    ModelFactory.__faker__.add_provider(BusinessProvider)
    ModelFactory.__faker__.add_provider(OrganismProvider)

    acronym = PostGenerated(derive_acronym)
    """Generate an acronym for the OTU derived from its name."""

    @classmethod
    def excluded_accessions(cls) -> set[str]:
        """Generate a set of excluded accessions."""
        return set()

    @post_generated
    @classmethod
    def isolates(cls, plan: Plan) -> list[IsolateBase]:
        """Derive a list of isolates from a plan."""
        isolates = []

        for _ in range(cls.__faker__.random_int(2, 5)):
            isolates.append(IsolateFactory.build_on_plan(plan=plan))

        return isolates

    id = Use(ModelFactory.__faker__.uuid4, cast_to=None)
    """Generate a UUID."""

    legacy_id = Use(ModelFactory.__faker__.legacy_id)
    """Generate an 8-character unique identifier as used in virtool-cli."""

    name = Use(ModelFactory.__faker__.organism)
    """Generate a realistic name for a plant virus."""

    plan = Use(PlanFactory.build)
    """Generate a quasi-realistic plan for the OTU."""

    @post_generated
    @classmethod
    def sequences(cls, isolates: list[IsolateBase]) -> list[SequenceBase]:
        """Derive a list of sequences from a list of isolates."""
        return [sequence for isolate in isolates for sequence in isolate.sequences]

    taxid = Use(ModelFactory.__faker__.random_int, min=1000, max=999999)
    """A realistic taxonomy ID."""


class OTUMinimalFactory(ModelFactory[OTUMinimal]):
    """OTUMinimal Factory with quasi-realistic data."""

    ModelFactory.__faker__.add_provider(BusinessProvider)
    ModelFactory.__faker__.add_provider(OrganismProvider)

    acronym = PostGenerated(derive_acronym)
    """An acronym for the OTU derived from its name."""

    legacy_id = Use(ModelFactory.__faker__.legacy_id)
    """Generate a realistic 8-character ``legacy_id`` for the OTU."""

    name = Use(ModelFactory.__faker__.organism)
    """Generate a realistic name for the OTU."""

    taxid = Use(ModelFactory.__faker__.random_int, min=1000, max=999999)
    """A realistic taxonomy ID."""
