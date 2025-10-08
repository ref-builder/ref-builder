# TODO: Complete NCBI Mock Client Migration

## Context
We're replacing the fragile file-based NCBI cache system with hardcoded Pydantic models for testing. This eliminates issues with:
- NCBI data changing over time (taxid merges, accession updates)
- pytest cache fragility
- Silent failures when cache data is missing
- Difficulty adding new test data

## Completed ✓
1. Created extraction script: `scripts/extract_ncbi_cache_to_fixtures.py`
2. Generated `tests/fixtures/ncbi_mock_data.py` with 75 GenBank + 26 Taxonomy records as Pydantic models
3. Created `tests/fixtures/mock_ncbi_client.py`:
   - `MockNCBIClient` class matching `NCBIClient` interface
   - `MockDataNotFoundError` exception for missing mock data
   - Clear error messages guiding developers to add missing data

## Next Steps

### 1. Update test fixtures in `tests/conftest.py`
- [ ] Create `mock_ncbi_client` fixture that returns `MockNCBIClient()`
- [ ] Update `scratch_event_store_data` fixture:
  - Replace `NCBIClient(False)` with `MockNCBIClient()`
  - Remove pytest caching logic (`pytestconfig.cache.get/set`)
  - Remove file I/O for event store persistence
  - Generate scratch repo data in-memory each time
- [ ] Remove `scratch_user_cache_path` fixture (no longer needed)
- [ ] Remove `precached_repo` fixture or update to use `MockNCBIClient`
- [ ] Remove `scratch_ncbi_cache` fixture (no longer needed)
- [ ] Remove `scratch_ncbi_client` fixture or update to use `MockNCBIClient`
- [ ] Update `uncached_ncbi_client` fixture to use `MockNCBIClient(ignore_cache=True)`

### 2. Validate test dataset
- [ ] Create validation script to check test data integrity:
  - Detect duplicates between plan and contents arrays
  - Identify taxid mismatches (e.g., merged taxids)
  - Verify all accessions in test data exist in mock data
  - Check for orphaned or unused mock data entries
  - Validate data structure consistency
- [ ] Run validation and document all issues found

### 3. Fix test data issues
- [ ] Fix `tests/files/src_test_contents.json`:
  - Remove `["NC_003355"]` from taxid 1169032 contents (duplicate of plan)
  - Verify no other plan/contents duplicates exist
- [ ] Handle taxid mismatches in mock data:
  - Current: taxid 1169032 in test data
  - Reality: NCBI merged it to taxid 3432896
  - Either: Update test data to use 3432896, or add both taxids to mock

### 4. Update test imports
- [ ] Search for `from ref_builder.ncbi.client import NCBIClient` in tests
- [ ] Replace with `from tests.fixtures.mock_ncbi_client import MockNCBIClient` where appropriate
- [ ] Keep real `NCBIClient` for integration tests marked with `@pytest.mark.ncbi`

### 5. Run and fix tests
- [ ] Run `pytest tests/cli/ -xvs` and fix failures
- [ ] Add missing accessions/taxids to `ncbi_mock_data.py` as `MockDataNotFoundError` is raised
- [ ] Verify all 37 CLI tests pass

### 6. Cleanup
- [ ] Remove `tests/files/cache_test/` directory (no longer needed)
- [ ] Remove `.pytest_cache` references from `.gitignore` if appropriate
- [ ] Update test documentation about using mocked NCBI data

### 7. Staged changes
- [ ] Apply staged changes to `tests/conftest.py`:
  ```diff
  -            otu_service.create(otu_contents.plan)
  -
  -        otu = repo.get_otu_by_taxid(otu_contents.taxid)
  +            otu = otu_service.create(otu_contents.plan)

           assert otu

           for isolate_accessions in otu_contents.contents:
               with repo.lock():
  -                assert add_genbank_isolate(
  +                add_genbank_isolate(
  ```

## Benefits After Completion
- ✅ No external dependencies during tests
- ✅ Test data is version-controlled Python code
- ✅ Easy to see what data exists (grep the source)
- ✅ Loud failures guide developers to add missing data
- ✅ No cache invalidation issues
- ✅ Faster tests (no file I/O)
- ✅ NCBI data updates don't break tests

## Current Blockers
- Test failures due to taxid mismatch (1169032 vs 3432896)
- Duplicate accessions in test data (plan vs contents)
- Need to remove pytest caching from fixtures
