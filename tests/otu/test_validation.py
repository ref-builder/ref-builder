from ref_builder.otu.validate import check_otu_is_valid
from ref_builder.repo import Repo


class TestValidateOTU:
    """Test OTU validation."""

    def test_ok(self, scratch_repo: Repo):
        """Validate scratch_repo contents."""
        for minimal_otu in scratch_repo.iter_minimal_otus():
            assert minimal_otu

            otu = scratch_repo.get_otu(minimal_otu.id)

            assert otu
            assert check_otu_is_valid(otu)

    def test_no_isolates(self, scratch_repo: Repo, mock_ncbi_client):
        """Test that validating an OTU with no isolates returns False."""
        otu = scratch_repo.get_otu_by_taxid(
            mock_ncbi_client.otus.saccharum_streak_virus.taxid
        )

        assert otu
        assert check_otu_is_valid(otu)
        assert len(otu.isolates) == 1

        isolate = otu.isolates[0]

        assert isolate

        otu_invalid = otu.model_copy()
        otu_invalid.delete_isolate(isolate.id)

        assert otu_invalid.isolates == []

        assert not check_otu_is_valid(otu_invalid)
