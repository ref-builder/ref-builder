from ref_builder.otu.validate import check_otu_is_valid
from ref_builder.repo import Repo


class TestValidateOTU:
    """Test OTU validation."""

    def test_ok(self, scratch_repo: Repo):
        """Validate scratch_repo contents"""
        for otu_metadata in scratch_repo.iter_minimal_otus():
            assert check_otu_is_valid(scratch_repo.get_otu_by_taxid(otu_metadata.taxid))

    def test_no_isolates(self, scratch_repo: Repo):
        """Test that validating an OTU with no isolates returns False."""
        otu = scratch_repo.get_otu_by_taxid(223262)

        assert otu
        assert check_otu_is_valid(otu)

        isolate = otu.isolates[0]

        assert isolate

        otu_invalid = otu.model_copy()
        otu_invalid.delete_isolate(isolate.id)

        assert otu_invalid.isolates == []

        assert not check_otu_is_valid(otu_invalid)
