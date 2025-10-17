"""Vallota mosaic virus (taxid=430059)."""

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

vallota_mosaic_virus = mock_ncbi_client.add_otu(
    taxid=430059,
    name="Vallota mosaic virus",
    refseq=["NC_043170"],
    taxonomy=NCBITaxonomy(
        id=430059,
        name="Vallota mosaic virus",
        other_names=NCBITaxonomyOtherNames(
            acronym=[], genbank_acronym=[], equivalent_name=[], synonym=[], includes=[]
        ),
        lineage=[
            NCBILineage(id=10239, name="Viruses", rank="superkingdom"),
            NCBILineage(id=2559587, name="Riboviria", rank="clade"),
            NCBILineage(id=2732396, name="Orthornavirae", rank="kingdom"),
            NCBILineage(id=2732408, name="Pisuviricota", rank="phylum"),
            NCBILineage(id=2732507, name="Stelpaviricetes", rank="class"),
            NCBILineage(id=2732550, name="Patatavirales", rank="order"),
            NCBILineage(id=39729, name="Potyviridae", rank="family"),
            NCBILineage(id=12195, name="Potyvirus", rank="genus"),
        ],
        rank=NCBIRank.SPECIES,
        species=NCBILineage(id=430059, name="Vallota mosaic virus", rank="species"),
    ),
    genbank=[
        NCBIGenbank(
            accession="NC_043170",
            accession_version="NC_043170.1",
            strandedness=Strandedness.SINGLE,
            moltype=MolType.RNA,
            topology=Topology.LINEAR,
            definition="Vallota mosaic virus polyprotein gene, partial cds",
            organism="Vallota mosaic virus",
            sequence="AGCCATCTACGGTAGTTGATAATACGCTGATGGTTATCCTAACCATGCAATATGCACTAGCCAAGCAGGACATCGGTTTCCAGGAACAAGAGGATATAATACGATATTTTGCCAATGGAGATGACTTACTGATTGCAGTGAATGGAGATAAGGGGATCGCCTTATTAAACACGCTGCAAGAATCATTCAGTGAGATGGGGTTGAACTATGACTTTAACGATCGTACTCACAATAAAAGCGAACTTAGTTTTATGTCTCACCAAGCGCTCGAATATGATGGTATGTATATACCAAAAATTAAGAAGGAGAGAATAGTCTCAATTTTGGAATGGGATAGGAGTGTGGAGCCAGAGCACAGAATGGAAGCTATTTGTGCCGCAATGGTTGAGGCATGGGGATATCCAGAACTGTTACATGAGATCCGTAAGTTCTATGCATTCATGTTGGATCAGGAGCCTTTCTCAGAATTGAATGCACAGGGTAGAGCGCATTACATATCAGAACAAGCTTTGAAGACGCTGTACATGGACGGCAAAGTAACTTTACTGGACATTGAACCATATCTTCAAGAGATAGCACATCTGAGCTTGGTAGATTTGGATGAGATGGTCTACCATCAGGCGGACAAAACGATAGATGCTGGCACCTCTACTTCTCAATCAGATCGAGCGCCACAGGTGGATCGAGACATAAACGCAGGTACTTTCGTCATTCCTCGCATTAAAGCACTAGGTGGGAAGATGGCGCTCCCAAAGGTGCGAGGGAAAAGTGTGATGAATCTACAACATCTCCTCACTTATTCCCCAGAACAGACTGACATTTCGAACACGCGTGCAACCCACAAACAATTTGCAACGTGGTATGATCGTGTCATGGAGAGTTATGGAGTGACTGATGCACAAATGGAAATTATTTTAAATGGGCTTATGGTTTGGTGTATTGAAAATGGAACTTCACCGAACTTGAGCGGAATGTGGACGATGATGGATAAAGATGAACAAGTTGAGTATCCATTGAAACCGATCTTAGAAAATGCGCAACCCACATTTAGACAAATTATGGCGCATTTTAGTAACGCAGCCGAAGCTTACATCGAAAAGAGAAATTCGGAAAGAAGGTATATGCCAAGGTTTGGCAGCCAACGTAATCTAACAGACTACAGCTTAGCACGCTATGCTTTCGATTTCTATGAAATAACTTCGCACACACCAGTGCGAGCAAGGGAAGCACACATCCAAATGAAAGCGGCGGCTCTTCGAAACACCAAAACGAGGATGTTTGGATTAGATGGAAAAGTTGGGACCGAGGAAGAGGACACGGAACGGCATGTGACGACAGATGTCAATCGGAACATGCATTCACTGTTAGGTGTTAACATGTGATTTAGTGAAGCAGGTATTTTACTTAATATCTATATTTGCCCGTCAGCATCCTTTTTCTTTTGCTTTAATGTATTAGCGTTTCAGTAAAGGTGATCCAGTTATTTCGTATGATTAACTGGAGTGACTTTGTTGTGAACACACATATGTTTTAGTCGAACTGGAATGAGTTGTGTATGTGATTGAACTTGAGCTACATAACACATTGCAGCAAGAGAAAAAAAAAAAAAAAA",
            source=NCBISource(
                taxid=430059,
                organism="Vallota mosaic virus",
                mol_type=NCBISourceMolType.GENOMIC_RNA,
                isolate="",
                host="Nerine sp.",
                segment="",
                strain="",
                clone="",
                proviral=False,
                macronuclear=False,
                focus=False,
                transgenic=False,
            ),
            comment="PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review. The reference sequence is identical to FJ618540.; COMPLETENESS: full length.",
            refseq=True,
        )
    ],
)
