# Pyright Type Checking Fixes

## Overview
This document organizes the 841 pyright errors found in the codebase and provides a plan to fix them systematically, starting with the most impactful and easiest fixes first.

## Error Categories

### 1. Pydantic Field Override Issues (🔥 High Impact - ~200+ errors)
**Root Cause:** Event classes inherit from `Event` but override `data` and `query` fields with incompatible types, violating Pydantic's invariant typing rules.

**Files Affected:**
- `ref_builder/events/isolate.py` (12 errors)
- `ref_builder/events/otu.py`
- `ref_builder/events/sequence.py`
- `ref_builder/events/repo.py`

**Pattern:**
```python
# Base class
class Event(BaseModel):
    data: EventData
    query: EventQuery

# Child class (PROBLEM)
class CreateIsolate(ApplicableEvent):
    data: CreateIsolateData  # ❌ Incompatible override
    query: IsolateQuery      # ❌ Incompatible override
```

**Solution:** Use Generic base classes or field validators instead of direct field overrides.

### 2. Optional Member Access (🚨 Critical - ~400+ errors)
**Root Cause:** Accessing attributes on values that could be `None` without proper null checking.

**Common Patterns:**
- `.name` attribute access on Optional objects
- `.isolates`, `.accessions`, `.plan` access without None checks
- UUID field access on Optional[UUID] types

**Files Most Affected:**
- Test files (`tests/test_*.py`)
- CLI modules (`ref_builder/cli/*.py`)
- Console module (`ref_builder/console.py`)

**Example Errors:**
```python
otu.name  # ❌ if otu could be None
# Should be:
if otu is not None:
    otu.name
```

### 3. Type Compatibility Issues (⚡ Quick Fixes - 5 errors)

#### Date/DateTime Mismatches
- `ref_builder/cli/otu.py:263` - passing `date` where `datetime` expected

#### Tuple Type Syntax
- `ref_builder/cli/utils.py:86` - using old tuple syntax instead of `tuple[T1, T2]`

#### UUID Type Inconsistencies
- Mixing `UUID4` and `UUID` types in function signatures

### 4. Return Type Mismatches (🔧 Medium Priority - ~10 errors)
**Root Cause:** Functions declared to return non-Optional types but actually returning Optional values.

**Example:**
```python
def get_otu() -> OTUBuilder:  # ❌ Says non-Optional
    return builder_or_none    # ❌ Actually Optional
```

## Implementation Plan

### Phase 1: Low Hanging Fruit (Quick Wins)
1. **Fix tuple type syntax** (1 error - 2 minutes)
2. **Fix date/datetime conversion** (2 errors - 5 minutes)
3. **Standardize UUID types** (5-10 errors - 10 minutes)

### Phase 2: Systematic Fixes (High Impact)
1. **Refactor Event class hierarchy** to fix Pydantic field overrides (~200 errors)
2. **Add systematic None checks** using search/replace patterns (~400 errors)

### Phase 3: Cleanup (Final Polish)
1. **Fix remaining return type inconsistencies**
2. **Verify all fixes with incremental pyright runs**

## Quick Reference

### Most Common Error Patterns
1. `"X" is not a known attribute of "None"` → Add None check
2. `Variable is mutable so its type is invariant` → Pydantic field override issue
3. `Type "X | None" is not assignable to type "X"` → Handle Optional properly
4. `Tuple expression not allowed` → Use `tuple[T1, T2]` syntax

### Files with Most Errors
- Test files (highest count due to test data access)
- Event definition files (Pydantic issues)
- CLI modules (Optional handling)

### Search Patterns for Bulk Fixes
- `\.name\b` - Find potential None attribute access
- `\.isolates\b` - Find isolates attribute access
- `\.accessions\b` - Find accessions attribute access
- `data: \w+Data$` - Find Pydantic field overrides