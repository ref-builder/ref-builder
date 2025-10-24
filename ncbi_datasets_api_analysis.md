# NCBI Datasets API v2 Integration Analysis for ref-builder

## Executive Summary

This report analyzes the potential integration of NCBI Datasets API v2 into the ref-builder project's NCBIClient for handling plant virus data. The analysis compares Datasets API capabilities with the current EUtils/BioPython implementation and proposes a dual-operation strategy for resilience.

## Current NCBIClient Implementation Overview

### Architecture
- **Primary Interface**: BioPython's Entrez module wrapping EUtils API
- **Databases Used**: `nuccore` (GenBank/RefSeq) and `taxonomy`
- **Caching Layer**: Local JSON-based cache for GenBank records and taxonomy data
- **Key Features**:
  - Fetches GenBank records by accession (ref_builder/ncbi/client.py:99-178)
  - Retrieves taxonomy records and builds lineages (ref_builder/ncbi/client.py:342-462)
  - Searches for accessions by taxid with filtering (ref_builder/ncbi/client.py:225-312)
  - Handles pagination for large result sets
  - Robust error handling with HTTP error logging

### Data Models
- **NCBIGenbank**: Validated Pydantic models for GenBank records
- **NCBITaxonomy**: Taxonomy records with rank validation
- **Lineage**: Hierarchical structure from species to isolate/strain
- **Focus**: Plant viruses with subspecific lineage tracking

## NCBI Datasets API v2 Capabilities

### Strengths for Virus Data
1. **Bulk Download Support**: Efficient retrieval of multiple genomes in single requests
2. **Modern Response Formats**: JSON, JSON Lines (streaming), TSV, Excel
3. **Dehydrated Archives**: Support for large downloads with delayed file retrieval
4. **Rich Metadata**: Comprehensive genome and annotation reports
5. **Higher Rate Limits**: With API key authentication (same as EUtils)

### Limitations for Plant Viruses
1. **Protein Dataset Restriction**: "Protein datasets are available for SARS-CoV-2 only"
2. **Virus Bias**: Special features focus on coronaviruses, not plant viruses
3. **Limited Taxonomy Integration**: Less granular than EUtils for lineage building
4. **No Direct Sequence Search**: Cannot search by sequence length or modification date like EUtils
5. **Incomplete API Documentation**: Virus-specific endpoints not fully documented in OpenAPI spec

## Comparative Analysis: Datasets API vs EUtils

| Feature | EUtils (Current) | Datasets API v2 | Winner for Plant Viruses |
|---------|-----------------|-----------------|-------------------------|
| **Accession Fetching** | efetch with XML parsing | Direct JSON response | Datasets (cleaner) |
| **Bulk Operations** | Sequential or batched | Native bulk support | Datasets |
| **Taxonomy Lineage** | Complete hierarchy access | Limited lineage data | EUtils |
| **Search Flexibility** | Complex queries with filters | Basic filtering | EUtils |
| **Response Parsing** | XML (requires BioPython) | JSON (native Python) | Datasets |
| **Error Handling** | HTTP errors with details | Structured error responses | Tie |
| **Caching** | Manual implementation | Manual implementation | Tie |
| **Plant Virus Support** | Full support | Generic virus support | EUtils |

## Recommended Integration Strategy

### 1. Hybrid Approach - Use Datasets API Where Superior

```python
class NCBIClient:
    def __init__(self, ignore_cache: bool, use_datasets_api: bool = True):
        self.use_datasets_api = use_datasets_api
        # ... existing initialization

    def fetch_genbank_records(self, accessions):
        if self.use_datasets_api and self._is_datasets_available():
            try:
                return self._fetch_via_datasets_api(accessions)
            except (HTTPError, ConnectionError) as e:
                logger.warning("Datasets API failed, falling back to EUtils", error=e)
                return self._fetch_via_eutils(accessions)
        return self._fetch_via_eutils(accessions)
```

### 2. Areas to Maintain EUtils

**Must Keep EUtils For:**
- Taxonomy lineage building (fetch_lineage method)
- Complex search queries with date/length filters (fetch_accessions_by_taxid)
- Plant virus-specific metadata that may not be in Datasets API

**Can Migrate to Datasets API:**
- Bulk GenBank record fetching (when available)
- Simple accession lookups
- Genome sequence retrieval

### 3. Implementation Priorities

#### Phase 1: Add Datasets API as Optional Path
- Implement `_fetch_via_datasets_api()` for GenBank records
- Add feature flag to enable/disable Datasets API
- Maintain all existing EUtils functionality

#### Phase 2: Optimize for Performance
- Use Datasets API for bulk operations (>10 accessions)
- Use EUtils for single accessions or complex queries
- Implement smart routing based on query type

#### Phase 3: Enhanced Error Handling
- Implement automatic failover between APIs
- Add monitoring for API health/performance
- Cache API availability status with TTL

## Dual Operation Strategy for Resilience

### 1. Health Check System

```python
class APIHealthChecker:
    def __init__(self):
        self.datasets_healthy = True
        self.eutils_healthy = True
        self.last_check = {}
        self.check_interval = 300  # 5 minutes

    def is_healthy(self, api: str) -> bool:
        if time.time() - self.last_check.get(api, 0) > self.check_interval:
            self._check_health(api)
        return self.datasets_healthy if api == "datasets" else self.eutils_healthy
```

### 2. Intelligent Routing

```python
def select_api(self, operation: str, params: dict) -> str:
    # Prefer Datasets for bulk, EUtils for complex queries
    if operation == "bulk_fetch" and params["count"] > 10:
        if self.health_checker.is_healthy("datasets"):
            return "datasets"

    # EUtils for taxonomy and complex searches
    if operation in ["taxonomy", "complex_search"]:
        return "eutils"

    # Default with failover
    return "datasets" if self.health_checker.is_healthy("datasets") else "eutils"
```

### 3. Graceful Degradation

- Primary: Datasets API for supported operations
- Fallback: EUtils when Datasets fails or lacks features
- Cache: Serve from cache when both APIs unavailable
- Error: Clear error messages with suggested retry

## Special Considerations for Plant Viruses

### 1. Data Availability
- Plant viruses are underrepresented in specialized Datasets API features
- No special metadata fields like those for coronaviruses
- Standard genome/annotation data should be sufficient

### 2. Taxonomy Challenges
- Plant virus taxonomy often includes subspecific ranks (strain, isolate)
- EUtils provides better granularity for lineage construction
- Maintain EUtils for taxonomy operations

### 3. Reference Building Requirements
- Focus on sequence quality and completeness
- Both APIs provide adequate sequence data
- Annotation quality varies regardless of API

## Implementation Recommendations

### Immediate Actions
1. **Keep EUtils as Primary**: Current implementation is stable and comprehensive
2. **Add Datasets API Experimentally**: Implement for bulk GenBank fetches only
3. **Monitor Performance**: Track success rates and response times for both APIs

### Future Enhancements
1. **Implement Smart Caching**: Cache API availability status
2. **Add Retry Logic**: Exponential backoff for both APIs
3. **Create API Abstraction Layer**: Unified interface hiding API details
4. **Performance Metrics**: Track and log API performance for optimization

### Code Organization
```
ref_builder/ncbi/
├── client.py          # Main client interface
├── cache.py           # Existing cache implementation
├── models.py          # Data models (unchanged)
├── eutils.py          # Extract EUtils-specific code
├── datasets.py        # New Datasets API implementation
└── health.py          # API health checking and routing
```

## Conclusion

While the NCBI Datasets API v2 offers modern features and better performance for bulk operations, it currently lacks the granularity needed for comprehensive plant virus data retrieval. The recommended approach is to:

1. **Maintain EUtils as the primary API** for its proven reliability with plant viruses
2. **Selectively integrate Datasets API** for bulk GenBank record fetching where it excels
3. **Implement intelligent failover** to ensure resilience when either API is unavailable
4. **Monitor both APIs** to make data-driven decisions about future migration

This hybrid approach maximizes reliability while positioning the codebase to take advantage of Datasets API improvements as they become available for plant virus research.

## Appendix: Key Differences Summary

| Aspect | Impact on Plant Virus Work |
|--------|---------------------------|
| **Coronavirus bias in Datasets** | Limited impact - basic genome data sufficient |
| **No protein datasets for plant viruses** | Not critical - nucleotide data is primary need |
| **Better bulk operations in Datasets** | High value - improves large-scale fetching |
| **Limited taxonomy in Datasets** | High impact - must retain EUtils for lineage |
| **Modern JSON responses** | Medium value - easier parsing but not critical |
| **Streaming with JSON Lines** | Low value - current batch sizes manageable |