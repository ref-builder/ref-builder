# Mock NCBI Client - Remaining Work

## Current Status

New manifest-based system operational with 9 test OTUs. Tests: **334 passing, 19 failing**.

## Phase 3: Cleanup ✅ Complete

### Remove Old System Artifacts ✅

1. ✅ ~~Delete `ref_builder/dev/otu_module.py.j2` Jinja template~~
2. ✅ ~~Remove old `dev refresh` logic that used template generation~~ (already migrated to JSON generation)

### Fix Test Failures - Progress: 11/30 Complete

**Completed fixes:**
- ✅ `tests/test_resources.py` - All 3 failures fixed (TestSequence, TestIsolate, TestOTU)
- ✅ `tests/cli/test_otu_commands.py` - 4 failures fixed (TestExcludeAccessionsCommand, TestAllowAccessionsCommand)
- ✅ `tests/test_add.py` - 1 failure fixed (test_plan_mismatch)
- ✅ `tests/test_modify.py` - All 11 failures fixed by using `mock_ncbi_client` fixture
- ✅ `tests/services/test_isolate.py` - 2 failures fixed (test_multipartite, test_refseq_exclusion)
- ✅ Fixed root cause: Changed `NCBISource` fields from `str = ""` to `str | None = None`
  - Empty NCBI fields no longer create invalid IsolateNames
  - `_extract_isolate_name_from_record` now correctly returns None when no isolate name exists
- ✅ Added missing EF accessions to `abaca_bunchy_top_virus` manifest

**Remaining failures (19 tests):**
- `tests/services/test_isolate.py` (~3 remaining)
- `tests/test_mock_ncbi_minimal.py` (~1)
- `tests/otu/test_utils.py` (~2)
- `tests/otu/test_validation.py` (~1)
- Others (~12)

**Solution pattern established:**
- Use `mock_ncbi_client` fixture parameter in test methods
- Access taxids via `mock_ncbi_client.otus.<otu_name>.taxid`
- Add missing accessions to manifest when needed
- Run `dev refresh` to regenerate mock data after manifest changes

## Documentation

Update relevant docs:
- How to add new test OTUs to manifest
- How to use `mock_ncbi_client.otus` in tests
- `dev refresh` command usage

## Benefits Achieved

- ✅ LSP autocomplete: `mock_ncbi_client.otus.wasabi_mottle_virus.taxid`
- ✅ Single source of truth in `manifest.py`
- ✅ Real NCBI data via `dev refresh`
- ✅ No circular dependencies
- ✅ No Jinja templates for Python code
- ✅ Species-level taxid handling matches production behavior
