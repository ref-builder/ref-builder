"""Faba bean necrotic stunt alphasatellite 1 (taxid=1441799)."""

from ref_builder.ncbi.models import (
    MolType,
    NCBIGenbank,
    NCBILineage,
    NCBIRank,
    NCBISource,
    NCBISourceMolType,
    NCBITaxonomy,
    NCBITaxonomyOtherNames,
    Strandedness,
    Topology,
)
from tests.fixtures.ncbi import mock_ncbi_client

faba_bean_necrotic_stunt_alphasatellite_1 = mock_ncbi_client.add_otu(
    taxid=1441799,
    name="Faba bean necrotic stunt alphasatellite 1",
    refseq=["NC_023881"],
    taxonomy=NCBITaxonomy(
        id=1441799,
        name="Faba bean necrotic stunt alphasatellite 1",
        other_names=NCBITaxonomyOtherNames(
            acronym=[], genbank_acronym=[], equivalent_name=[], synonym=[], includes=[]
        ),
        lineage=[
            NCBILineage(id=10239, name="Viruses", rank="superkingdom"),
            NCBILineage(id=1458186, name="Alphasatellitidae", rank="family"),
            NCBILineage(id=1441744, name="Nanoalphasatellitinae", rank="subfamily"),
            NCBILineage(
                id=2683341, name="unclassified Nanoalphasatellitinae", rank="no rank"
            ),
        ],
        rank=NCBIRank.SPECIES,
        species=NCBILineage(
            id=1441799, name="Faba bean necrotic stunt alphasatellite 1", rank="species"
        ),
    ),
    genbank=[
        NCBIGenbank(
            accession="NC_023881",
            accession_version="NC_023881.1",
            strandedness=Strandedness.SINGLE,
            moltype=MolType.DNA,
            topology=Topology.CIRCULAR,
            definition="Faba bean necrotic stunt alphasatellite 1 isolate Peshtatuek_12b, complete sequence",
            organism="Faba bean necrotic stunt alphasatellite 1",
            sequence="ACCCCGCCTTGGAACACCTCCTTGGAACGGGTATAAATAGGATTTTAATTTCTAGAAAATTAAAATGGCCTGTGCGAATTGGGTTTTCACTCGCAATTTTCAAGGAGCTCTCCCTTCTCTCTCGTTCGACGAGAGAGTTCAATACGCTGTCTGGCAACACGAAAGAGGAACTCATGACCATATCCAGGGAGTAATTCAATTGAAGAAAAAAGCTCGATTTTCGACTGTTAAGGAGATAATTGGGGGAAATCCTCATGTGGAGAAAATGAAAGGTACAATTGAAGAAGCTTCAGCTTATGTCCAGAAAGAAGAAACAAGAGTTGCAGGTCCCTGGAGTTATGGTGATTTATTGAAGAGAGGATCTCACAGGAGGAAGACGATGGAGAGATATTTAGAAGACCCAGAAGAAATGAAATTGAAGGACCCAGATGTTGCTCTTCGCTGTAATGCGAAGAGACTGAAGGAAGATTATTGTTCTTGTTTTTCTTCTTTTAAACTTCGTCCTTGGCAAATTGAGCTTCATCGGGTTTTAATGAATGAACCAGATGATCGTTCTATCATCTGGGTTTATGGTCCGGACGGAGGAGAAGGAAAGAGTACATTTGCTAAGGAATTAATTAAATATGGATGGTTTTATACTGCTGGAGGGAAGACGCAAGACGTTCTGTATATGTATGCTCAAGATCCAGAGAGAAATATTGCTTTTGATGTACCCAGATGTTCTTCGGAAATGATGAACTATCAAGCTATGGAGATGATGAAGAATAGATGCTTTGCAAGTACGAAATATAGATCGATAGATCTTTGTGTTAGGAAAAATGTATTTTTAGTTGTTTTTGCAAACGTGGAACCAGACCCCACAAAAATAAGTGGGGACAGAATTGTCATTATCAACTGTTGAATTTGAATTTCAAACTAAGCGAAGCGGGAAATTTCCCGCTCTTTTGTTCATTTAGCAAAACGCTGTCGTTTTTACCTTGGACCAAGGCGGGTTTAGTATT",
            source=NCBISource(
                taxid=1441799,
                organism="Faba bean necrotic stunt alphasatellite 1",
                mol_type=NCBISourceMolType.GENOMIC_DNA,
                isolate="Peshtatuek_12b",
                host="Lens culinaris",
                segment="",
                strain="",
                clone="",
                proviral=False,
                macronuclear=False,
                focus=False,
                transgenic=False,
            ),
            comment="PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review. The reference sequence is identical to KC978991.; GenBank Accession Numbers KC978982-KC978991 represent sequences from the 10 segments of Faba bean necrotic stunt virus.; ##Assembly-Data-START## ; Sequencing Technology :: Sanger dideoxy sequencing ; ##Assembly-Data-END##; COMPLETENESS: full length.",
            refseq=True,
        )
    ],
)
