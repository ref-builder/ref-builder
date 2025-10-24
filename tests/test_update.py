
from ref_builder.ncbi.client import NCBIClient
from ref_builder.otu.promote import promote_otu_accessions
from ref_builder.repo import Repo
from ref_builder.services.cls import Services


class TestPromoteOTU:
    """Test OTU accession promotion from Genbank to RefSeq."""

    def test_ok(self, empty_repo: Repo):
        """Test that RefSeq accessions can be promoted automatically."""
        ncbi_client = NCBIClient(ignore_cache=False)
        services = Services(empty_repo, ncbi_client)

        with empty_repo.lock():
            otu = services.otu.create(["MF062136", "MF062137", "MF062138"])

            assert otu

            isolate = services.isolate.create(
                otu.id, ["MF062125", "MF062126", "MF062127"]
            )

            assert isolate

        otu_before = empty_repo.get_otu(otu.id)

        assert otu_before
        assert otu_before.accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
            "MF062136",
            "MF062137",
            "MF062138",
        }

        isolate_before = otu_before.get_isolate(isolate.id)

        assert isolate_before
        assert isolate_before.accessions == {
            "MF062125",
            "MF062126",
            "MF062127",
        }

        with empty_repo.lock():
            promoted_accessions = promote_otu_accessions(empty_repo, otu_before)

        assert promoted_accessions == {"NC_055390", "NC_055391", "NC_055392"}

        otu_after = empty_repo.get_otu(otu.id)

        assert otu_after
        assert otu_after.isolate_ids == otu_before.isolate_ids
        assert otu_after.accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
            "MF062136",
            "MF062137",
            "MF062138",
        }
        assert otu_after.excluded_accessions == {"MF062125", "MF062126", "MF062127"}

        isolate_after = otu_after.get_isolate(isolate_before.id)

        assert isolate_after
        assert isolate_after.accessions == {
            "NC_055390",
            "NC_055391",
            "NC_055392",
        }
