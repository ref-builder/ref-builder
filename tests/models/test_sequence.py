"""Tests for Sequence model."""

from uuid import uuid4

from pydantic import ValidationError

from ref_builder.models.accession import Accession
from ref_builder.models.sequence import Sequence


class TestSequence:
    """Test the ``Sequence`` model which is used for complete validation of sequences."""

    def test_ok(self):
        """Test that a valid sequence passes validation."""
        assert Sequence.model_validate(
            {
                "id": uuid4(),
                "accession": Accession("NC_001234", 1),
                "definition": "Test virus, complete genome",
                "segment": uuid4(),
                "sequence": "ATCGATCGATCG",
            }
        )

    def test_invalid_accession(self):
        """Test that an invalid accession fails validation."""
        bad_sequence_data = {
            "id": uuid4(),
            "accession": Accession("BADSEQUENCE", 1),
            "definition": "Test virus, complete genome",
            "segment": uuid4(),
            "sequence": "ATCGATCGATCG",
        }

        try:
            Sequence.model_validate(bad_sequence_data)
        except ValidationError as e:
            for error in e.errors():
                assert "accession" in error["loc"]
                assert (
                    "Accession BADSEQUENCE.1 does not match a valid accession pattern"
                    in error["msg"]
                )
