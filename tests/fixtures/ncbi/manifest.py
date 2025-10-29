"""Registry of mock OTU data for testing.

Each attribute declares an OTU with its required accessions.
Run `ref-builder dev refresh` to fetch data from NCBI.
"""

from tests.fixtures.ncbi.models import OTUSpec


class OTUManifest:
    """Manifest declaring all test OTUs."""

    abaca_bunchy_top_virus = OTUSpec(
        refseq=[
            "NC_010314.1",
            "NC_010315.1",
            "NC_010316.1",
            "NC_010317.1",
            "NC_010318.1",
            "NC_010319.1",
        ],
        isolates=[
            [
                "EF546808.1",
                "EF546809.1",
                "EF546810.1",
                "EF546811.1",
                "EF546812.1",
                "EF546813.1",
            ],
            [
                "EF546802.1",
                "EF546803.1",
                "EF546804.1",
                "EF546805.1",
                "EF546806.1",
                "EF546807.1",
            ],
        ],
    )

    babaco_mosaic_virus = OTUSpec(
        refseq=["NC_036587.1"],
        isolates=[],
    )

    beet_black_scorch_virus = OTUSpec(
        refseq=[
            "NC_004452.1",
        ],
        isolates=[
            ["FN565520.1"],
            ["JN635327.1"],
            ["MH399846.1"],
            ["EF153268.1"],
            ["NC_004452.2"],
            ["NC_004452.3"],
        ],
    )

    cabbage_leaf_curl_jamaica_virus = OTUSpec(
        refseq=["DQ178610.1", "DQ178611.1"],
        isolates=[["DQ178608.1", "DQ178609.1"]],
    )

    dahlia_latent_viroid = OTUSpec(
        refseq=["NC_020160.1"],
    )

    east_african_cassava_mosaic_cameroon_virus = OTUSpec(
        refseq=["NC_004625.1", "NC_004630.1"],
    )

    oat_blue_dwarf_virus = OTUSpec(
        refseq=["NC_001793.1"],
    )

    okra_leaf_curl_alphasatellite = OTUSpec(
        refseq=["NC_005954.1"],
    )

    saccharum_streak_virus = OTUSpec(
        refseq=["NC_013464.1"],
    )

    tobacco_mosaic_virus = OTUSpec(
        refseq=["NC_001367.1"],
        isolates=[
            ["OQ953825.1"],
            ["HE818414.1"],
            ["AJ011933.1"],
            ["AF395128.1"],
            ["V01408.1"],
        ],
    )

    wasabi_mottle_virus = OTUSpec(
        refseq=["KJ207375.1"],
        isolates=[
            ["MK431779.1"],
            ["NC_003355.1"],
            ["AB017504.1"],
            ["MH200607.1"],
        ],
    )
