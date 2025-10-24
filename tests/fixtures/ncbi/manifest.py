"""Registry of mock OTU data for testing.

Each attribute declares an OTU with its required accessions.
Run `ref-builder dev refresh` to fetch data from NCBI.
"""

from tests.fixtures.ncbi.models import OTUSpec


class OTUManifest:
    """Manifest declaring all test OTUs."""

    abaca_bunchy_top_virus = OTUSpec(
        refseq=[
            "NC_010314",
            "NC_010315",
            "NC_010316",
            "NC_010317",
            "NC_010318",
            "NC_010319",
        ],
        isolates=[
            ["EF546808", "EF546809", "EF546810", "EF546811", "EF546812", "EF546813"],
            ["EF546802", "EF546803", "EF546804", "EF546805", "EF546806", "EF546807"],
        ],
    )

    babaco_mosaic_virus = OTUSpec(
        refseq=["NC_036587"],
        isolates=[["MT240490"]],
    )

    cabbage_leaf_curl_jamaica_virus = OTUSpec(
        refseq=["DQ178610", "DQ178611"],
        isolates=[["DQ178608", "DQ178609"]],
    )

    dahlia_latent_viroid = OTUSpec(
        refseq=["NC_020160"],
    )

    east_african_cassava_mosaic_cameroon_virus = OTUSpec(
        refseq=["NC_004625", "NC_004630"],
    )

    oat_blue_dwarf_virus = OTUSpec(
        refseq=["NC_001793"],
    )

    okra_leaf_curl_alphasatellite = OTUSpec(
        refseq=["NC_005954"],
    )

    saccharum_streak_virus = OTUSpec(
        refseq=["NC_013464"],
    )

    tobacco_mosaic_virus = OTUSpec(
        refseq=["NC_001367"],
        isolates=[["OQ953825"], ["HE818414"], ["AJ011933"], ["AF395128"], ["V01408"]],
    )

    wasabi_mottle_virus = OTUSpec(
        refseq=["KJ207375"],
        isolates=[
            ["MK431779"],
            ["NC_003355"],
            ["AB017504"],
            ["MH200607"],
        ],
    )
