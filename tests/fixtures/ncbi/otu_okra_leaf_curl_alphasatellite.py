"""Okra leaf curl alphasatellite (taxid=518829)."""

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

okra_leaf_curl_alphasatellite = mock_ncbi_client.add_otu(
    taxid=518829,
    name="Okra leaf curl alphasatellite",
    refseq=["NC_005954"],
    taxonomy=NCBITaxonomy(
        id=518829,
        name="Okra leaf curl alphasatellite",
        other_names=NCBITaxonomyOtherNames(
            acronym=[], genbank_acronym=[], equivalent_name=[], synonym=[], includes=[]
        ),
        lineage=[
            NCBILineage(id=10239, name="Viruses", rank="acellular root"),
            NCBILineage(id=3403061, name="Viruses incertae sedis", rank="no rank"),
            NCBILineage(id=1458186, name="Alphasatellitidae", rank="family"),
            NCBILineage(id=283494, name="Geminialphasatellitinae", rank="subfamily"),
            NCBILineage(
                id=1231296,
                name="unclassified Begomovirus-associated alphasatellites",
                rank="no rank",
            ),
        ],
        rank=NCBIRank.SPECIES,
        species=NCBILineage(
            id=518829, name="Okra leaf curl alphasatellite", rank="species"
        ),
    ),
    genbank=[
        NCBIGenbank(
            accession="NC_005954",
            accession_version="NC_005954.1",
            strandedness=Strandedness.DOUBLE,
            moltype=MolType.DNA,
            topology=Topology.CIRCULAR,
            definition="Okra leaf curl alphasatellite rep gene for replication associated protein, clone NOB-3",
            organism="Okra leaf curl alphasatellite",
            sequence="ACCCCGCCTCGAGACTTTAGCCCGAGACCCTCCGTACAGTTTTTTGTGCTCTCGCTTATATTTCTGTCTTCTGCGATAGATGGCTGCAATTAAGTCACATTGGTGGTGTTTCACGGTGTTCTTCCTCTCTGCTACTGCACCTGACTTGGTGCCTCTGTTCGAAAACACTCACGTGAGTTATGCTTGTTGGCAAGAAGAGGAGTCTCCGACGACGAAACGTCGCCACCTTCAAGGATACCTGCAATTGAAGGGTCAGAGGACCCTGAACCAGGTGAAGGCCATATTTGGGGATTTGAAGCCCCATCTTGAGAAACAGCGAGCTCGTAAGACAGACGATGCTCGCGATTACTGTATGAAACCCGAAACTAGGGTTTCTGGCCCCTTTGAATTTGGGGAATACTGTCCTGCTGGTTCTCACAAACGACGACAAAGGGAACTCGTAATTCGATCTCCGGTGAGAATGGCAGAGGAAAATCCGTCCGTCTTCCGACGAGTAAAGGCAAAGATTGCTGAGGAAGAATTCCAGAAGAGCGCTCATGAGATTCAAATTTCAAATTTGAAATCTTGGCAATCGCGCCTAAAGACGCTCCTGGAGAGGGACCCAGACGACCGCACTATTCTCTGGGTTTACGGACCTAATGGTGGGGAAGGGAAATCCACTTTCGCCAGAGACCTATACAGAAGTGGGTCCTGGTTCTATACACGTGGAGGATCTGCAGATAATGTAGCTTACCAGTACATAGGATGTTTAGGCAATAATATTGTATTTGATATTCCTCGTGATAAGAAGGAGTATCTTCAGTATAGTTTAATAGAGATGTTTAAGGATAGATTAATAGTTAGTAATAAGTACGAGCCTCTTATGGCCCCTTTGCTTAATTGTATTCATGTTGTAGTTATGTCTAATTTTCTCCCAGACTTTGAGAAGATTAGTCAGGATAGAGTCCATGTAATCCCATGTATACCATGTGGTGTTTGTCTTAAACACCATAATATTAATGATAAATGTGATGATTATGTGGATTAATTATTTTTGTTTCTTGAATACAAGAAAGAATGAAATGAAAAAAAAACAAAGAAAAAACATGAGCCTTCTATTATTTTAAGATAACAGTTGCCGCGCAGCGGCATTCAAAAAAAAAATGGAAAAGAAAAATAAAATGAAAAAGAAAAATAAAATGTCTTATTAAAACGACGACGTATTGAATAAATTGATCTGGGAGTAAAATGTAAATATTGGAAACTTGCGGGCTTGCGCTATAAATATTTAAAAGTTACGGTTAAGCGGAGTAAATGTAGAAGTCTTCAGGGGTATGTGTTTGCAATTATTTTGGTCTCTCAGCTATAAATAACACGTCACGAGGCGGGTGTAGTATT",
            source=NCBISource(
                taxid=518829,
                organism="Okra leaf curl alphasatellite",
                mol_type=NCBISourceMolType.GENOMIC_DNA,
                isolate="",
                host="Okra",
                segment="",
                strain="",
                clone="NOB-3",
                proviral=False,
                macronuclear=False,
                focus=False,
                transgenic=False,
            ),
            comment="PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review. The reference sequence is identical to AJ512954.; COMPLETENESS: full length.",
            refseq=True,
        )
    ],
)
