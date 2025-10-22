from pathlib import Path

import arrow
import orjson

from ref_builder.build import build_json
from ref_builder.repo import Repo


def test_ok(scratch_repo: Repo, tmp_path: Path):
    """Test that build_json produces deterministic output with correct business logic."""
    # Build twice to test determinism
    output_1 = tmp_path / "ref1.json"
    output_2 = tmp_path / "ref2.json"

    build_json(output_1, scratch_repo.path, "2.1.0")
    build_json(output_2, scratch_repo.path, "2.1.0")

    ref_1 = orjson.loads(output_1.read_bytes())
    ref_2 = orjson.loads(output_2.read_bytes())

    # Test determinism (same input produces same output)
    assert {**ref_1, "created_at": ""} == {**ref_2, "created_at": ""}

    # Test timestamp is current
    assert (arrow.utcnow() - arrow.get(ref_1["created_at"])).seconds == 0

    # Test OTUs are sorted alphabetically
    otu_names = [otu["name"] for otu in ref_1["otus"]]
    assert otu_names == sorted(otu_names)

    # Test each OTU has exactly one default isolate
    for otu in ref_1["otus"]:
        default_count = sum(1 for iso in otu["isolates"] if iso["default"])
        assert default_count == 1
