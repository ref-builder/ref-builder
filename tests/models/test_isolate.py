"""Tests for Isolate model."""

from uuid import uuid4

import pytest
from pydantic import ValidationError

from ref_builder.models.accession import Accession
from ref_builder.models.isolate import Isolate


class TestIsolate:
    """Test the ``Isolate`` model which is used for complete validation of isolates."""

    def test_accession_consistency_warning(self):
        """Test that a validation error is raised if accession provenances are mixed."""
        segment_id = uuid4()

        # Create isolate with mixed RefSeq and GenBank accessions
        mixed_isolate_data = {
            "id": uuid4(),
            "name": None,
            "taxid": 12345,
            "sequences": [
                {
                    "id": uuid4(),
                    "accession": Accession("NC_000001", 1),  # RefSeq
                    "definition": "Test virus segment A",
                    "segment": segment_id,
                    "sequence": "ATCGATCGATCG",
                },
                {
                    "id": uuid4(),
                    "accession": Accession("BD000001", 1),  # GenBank
                    "definition": "Test virus segment B",
                    "segment": segment_id,
                    "sequence": "GCTAGCTAGCTA",
                },
            ],
        }

        with pytest.raises(
            ValidationError,
            match="Combination of RefSeq and non-RefSeq sequences found in multipartite isolate",
        ):
            Isolate.model_validate(mixed_isolate_data)
