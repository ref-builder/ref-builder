"""Core models for manifest-driven mock NCBI data."""

import json
from dataclasses import dataclass
from pathlib import Path


class OTUSpec:
    """Declares the GenBank accessions needed for a mock OTU."""

    def __init__(
        self,
        refseq: list[str],
        isolates: list[list[str]] | None = None,
        incompatible: list[str] | None = None,
    ):
        self.refseq = refseq
        self.isolates = isolates or []
        self.incompatible = incompatible or []
        self.name: str = ""

    def __set_name__(self, owner, name):
        """Capture attribute name for file mapping."""
        self.name = name

    @property
    def all_accessions(self) -> list[str]:
        """All accessions that should be fetchable."""
        result = list(self.refseq)
        for isolate in self.isolates:
            result.extend(isolate)
        result.extend(self.incompatible)
        return result


@dataclass
class OTUHandle:
    """Handle to a mock OTU's metadata."""

    name: str
    taxid: int


class OTURegistry:
    """Dynamic registry populated from manifest and JSON files."""

    def __init__(self, manifest: type, data_dir: Path):
        """Load all available OTUs without validation."""
        self._manifest = manifest
        self._data_dir = data_dir

        for attr_name in dir(manifest):
            if attr_name.startswith("_"):
                continue
            attr = getattr(manifest, attr_name)
            if isinstance(attr, OTUSpec):
                handle = self._load_otu(attr_name, attr, data_dir)
                if handle:
                    setattr(self, attr_name, handle)

    def _load_otu(self, name: str, spec: OTUSpec, data_dir: Path) -> OTUHandle | None:
        """Load OTU if JSON exists, return None otherwise."""
        json_path = data_dir / f"{name}.json"
        if not json_path.exists():
            return None

        data = json.loads(json_path.read_text())

        # Mirror OTUService.create() behavior: return species-level taxid
        # If taxonomy rank is not 'species', find species-level ancestor in lineage
        taxonomy = data["taxonomy"]
        if taxonomy["rank"] != "species":
            for lineage_entry in taxonomy["lineage"]:
                if lineage_entry["rank"] == "species":
                    return OTUHandle(name=name, taxid=lineage_entry["id"])

        return OTUHandle(name=name, taxid=taxonomy["id"])

    def validate(self) -> None:
        """Validate that all manifest entries have data and match."""
        for attr_name in dir(self._manifest):
            if attr_name.startswith("_"):
                continue
            attr = getattr(self._manifest, attr_name)
            if not isinstance(attr, OTUSpec):
                continue

            if not hasattr(self, attr_name):
                raise RuntimeError(
                    f"Missing data for OTU '{attr_name}'. Run: ref-builder dev refresh"
                )

            json_path = self._data_dir / f"{attr_name}.json"
            data = json.loads(json_path.read_text())
            genbank_accessions = list(data["genbank"].keys())

            # Match manifest entries (versioned or unversioned) to data
            refseq_in_data = []
            for manifest_acc in attr.refseq:
                # Exact match (for versioned entries like "NC_004452.1")
                if manifest_acc in genbank_accessions:
                    refseq_in_data.append(manifest_acc)
                # Unversioned match (for entries like "NC_001367" matching "NC_001367.1")
                else:
                    for data_acc in genbank_accessions:
                        if data_acc.startswith(f"{manifest_acc}."):
                            refseq_in_data.append(data_acc)
                            break

            if sorted(attr.refseq) != sorted(refseq_in_data):
                raise RuntimeError(
                    f"Manifest/data mismatch for '{attr_name}'. Run: ref-builder dev refresh"
                )
