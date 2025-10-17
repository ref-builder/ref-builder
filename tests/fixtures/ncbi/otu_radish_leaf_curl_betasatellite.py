"""Radish leaf curl betasatellite (taxid=662596)."""

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

radish_leaf_curl_betasatellite = mock_ncbi_client.add_otu(
    taxid=662596,
    name="Radish leaf curl betasatellite",
    refseq=["NC_010239"],
    taxonomy=NCBITaxonomy(
        id=662596,
        name="Radish leaf curl betasatellite",
        other_names=NCBITaxonomyOtherNames(
            acronym=[], genbank_acronym=[], equivalent_name=[], synonym=[], includes=[]
        ),
        lineage=[
            NCBILineage(id=10239, name="Viruses", rank="superkingdom"),
            NCBILineage(id=1993640, name="Tolecusatellitidae", rank="family"),
            NCBILineage(id=190729, name="Betasatellite", rank="genus"),
        ],
        rank=NCBIRank.SPECIES,
        species=NCBILineage(
            id=662596, name="Radish leaf curl betasatellite", rank="species"
        ),
    ),
    genbank=[
        NCBIGenbank(
            accession="NC_010239",
            accession_version="NC_010239.1",
            strandedness=Strandedness.SINGLE,
            moltype=MolType.DNA,
            topology=Topology.CIRCULAR,
            definition="Radish leaf curl virus satellite DNA beta, complete genome",
            organism="Radish leaf curl betasatellite",
            sequence="ACCGGTGGCGAGCGGATTTTTTTGCGTCGTTGTGGGACCCACTTGTTAGATCTGCAGAATAAATGGGCCTTTTTATGCTATTGGGCCTATATTTATTGGGCTTTGTATAATATATTTAAATTGGGCTTTGTAGTCATGAGAATAAGATAATAAAACAGTTTCATTGATTATCTGTATAAATACGCATTTATACAGATGAACGCGTATACACATCGTATTCATCCATTACATTAATGTCTATTACAGGAGCCTCTTCCATCATCAGCATATCAATTCCTTCTATCATGTCTTCTTGCTTGAATTCTCCTATAGTAGACTCTTTGTACATGATTGCTAACAGGTTGTGTATTCTTTCTTCCAGAGTGTTGAAGTTGAATGGTGGTATGATTCCACTATGCCCGTATGGAATCATGAATCTCCTTGTTGCTAGCCCCGGAGATCGTGTTGAGAATATTTGTATTTGCACTAGAATGGAGTCATGTTCTTTCAGGCGTACGTCTATGATGAATTCCATCCCTTTCTGGTTTTTGTATTTGATCGTCATTGTGTTGTGTGAGTGTTCATATATGATGAACACTCATTAAATAGAACATATATTTGTGATGTATGGTTATTATTGATATACGTGGTTTGTTTATGTTCATATATTATGTAGTGGAGATATTTGTTATGGATGTGGACGTTGATATTATTTGGTTTCGTGAGAGTGAGTTTTAATGAATCCTAATTATGATATAATTAGGAAGAAAGAAAGAAAAAGAAAAGAATGAAAAAATAGAAAATGAAAAAAATAGAAAACTACATAAACTATAAATTCTATTCTATGAAAATGGTCGCGCAGCGAAACAGAAAACCGAAACTAACAAAAGAAATAAAAAAAATGAAAAAAGAAATAAACAACAGATGAAACGTGAGAAGAGAAAATAAAAACAAGAAACATCAGTGGTCCCCACTGAATAAATAAAAATAAAAAAAGAAAAAGAAGAAAATGATTTAATCCTTATACTTGAGTTACTTACTGTTTTACCGAACGGTAAAATAGTAATTGATATTTACACAAGGGTAAATATCTGAGTCCCCGATAGGTAAATGTGTCCCCAATATATCGGGGACTCAATTGGTGACCCATATTTTGCTTTCCTAAAATACCCTTGTCTCTGTGTCTGGTAGGCGCGTGGGAGTGGACTGAAAAAGTAGAATTTCTCTCTCCTAAAACTCGTCGGAACTCCGATTAAGGCACTTCCGGTCACCAATTTGCGACACGCGCGGCGGTGCGTACCCCTGGGAGGGTAGGGTACCACTACGCTACGCAGCAGCCTTAGCTACGCCGGAGCTTAGCTCGCCACCGTTCTAATATT",
            source=NCBISource(
                taxid=662596,
                organism="Radish leaf curl betasatellite",
                mol_type=NCBISourceMolType.GENOMIC_DNA,
                isolate="Varanasi",
                host="radish cultivar Japanese white",
                segment="",
                strain="",
                clone="",
                proviral=False,
                macronuclear=False,
                focus=False,
                transgenic=False,
            ),
            comment="PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review. The reference sequence was derived from EF175734.; COMPLETENESS: full length.",
            refseq=True,
        )
    ],
)
