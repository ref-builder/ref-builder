"""Dahlia latent viroid (taxid=1278205)."""

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

dahlia_latent_viroid = mock_ncbi_client.add_otu(
    taxid=1278205,
    name="Dahlia latent viroid",
    refseq=["NC_020160"],
    taxonomy=NCBITaxonomy(
        id=1278205,
        name="Dahlia latent viroid",
        other_names=NCBITaxonomyOtherNames(
            acronym=[], genbank_acronym=[], equivalent_name=[], synonym=[], includes=[]
        ),
        lineage=[
            NCBILineage(id=10239, name="Viruses", rank="superkingdom"),
            NCBILineage(id=185751, name="Pospiviroidae", rank="family"),
            NCBILineage(id=147262, name="Hostuviroid", rank="genus"),
        ],
        rank=NCBIRank.SPECIES,
        species=NCBILineage(id=1278205, name="Dahlia latent viroid", rank="species"),
    ),
    genbank=[
        NCBIGenbank(
            accession="NC_020160",
            accession_version="NC_020160.1",
            strandedness=Strandedness.SINGLE,
            moltype=MolType.RNA,
            topology=Topology.LINEAR,
            definition="Dahlia latent viroid, complete sequence",
            organism="Dahlia latent viroid",
            sequence="CAGGTCTTCTAAGGGTTCCTGTGGTGCCTCCCTCCAAGGCCGCGTAGGGAAAGAAAAAGTGAAAGAAAGTGTACCTGAAGAGAGAGAGACTCTCTGAGGAGCCCCGGGGCAACTCCGAGAGTGCTGCGGCAGGAGGAGCCGGGGGCGGAGGTTGTCTCGCTCTTGAGAGCTCCGCTCCTTGTAGCTTTGAGACTACCGCCCTTTTGCTTCCTTCTCGCTGGCTGACTCGAGGACGCGACCGGTGGTACACCAAGAGGTCTCCACCTCCTAGGGTACTTTTTTCTAACACCGATTCCGTACCGAGCGCCGGAGAGTGAAGCGCCCAGGGCTCCAAAGAAGCCC",
            source=NCBISource(
                taxid=1278205,
                organism="Dahlia latent viroid",
                mol_type=NCBISourceMolType.GENOMIC_RNA,
                isolate="4706174",
                host="Dahlia pinnata Cav.",
                segment="",
                strain="",
                clone="",
                proviral=False,
                macronuclear=False,
                focus=False,
                transgenic=False,
            ),
            comment="PROVISIONAL REFSEQ: This record has not yet been subject to final NCBI review. The reference sequence is identical to JX263426.; COMPLETENESS: full length.",
            refseq=True,
        )
    ],
)
