# Implementation Notes: Replication API Naming Unification

## Objective

Unify the naming scheme for "host" vs "node" replication across WorldContainer and derived classes, add string-to-enum converters, and maintain backward compatibility.

## Design Decisions

### 1. Naming Convention Choice

**Selected:** `node` over `host` terminology

**Rationale:**
- The enum is already `NodeReplicated` (not `HostReplicated`)
- Consistency with existing enum naming
- "Node" is the more common term in HPC literature for compute nodes
- "Host" can be ambiguous (host vs client, host system, etc.)

**Canonical Method Names:**
- `rank_replication()` - Query method for rank-level replication
- `node_replication()` - Query method for node-level replication
- `replicate_on_ranks()` - Action method for rank-level replication
- `replicate_on_nodes()` - Action method for node-level replication

### 2. String Converter Design

**Pattern Used:** Class-specific converter names to avoid overload ambiguity

**Exchange::Algorithm:**
- Primary: `Exchange<T,NDIM>::from_string_algorithm(std::string)`
- Deprecated: `Exchange<T,NDIM>::from_string(std::string)`
- Free function: `algorithm_from_string<T,NDIM>(std::string)`
- Wrapper: `AlgorithmFromString<T,NDIM>`
- User literal: `operator"" _alg` (for double,3D)

**DistributionType:**
- Primary: `distribution_type_from_string(std::string)`
- Wrapper: `DistributionTypeFromString`

**Normalization Strategy:**
All converters apply consistent normalization:
1. Convert to lowercase
2. Replace spaces with underscores
3. Replace dashes with underscores
4. Match against known variants

This allows inputs like:
- "Small-Memory" → `small_memory`
- "rank replicated" → `RankReplicated`
- "NodeReplicated" → `NodeReplicated`

### 3. Backward Compatibility

**Approach:** Deprecated aliases with `[[deprecated]]` attribute

**Advantages:**
- Existing code continues to work
- Compiler warnings guide users to new names
- Clean migration path
- Can be removed in future major version

**Implementation:**
```cpp
// New canonical method
bool rank_replication() const { ... }

// Deprecated alias
[[deprecated("Use rank_replication() instead")]]
bool is_replicated() const { return rank_replication(); }
```

### 4. Internal Implementation

**WorldContainerImpl Methods:** Kept unchanged

**Rationale:**
- Private implementation details
- No API breakage
- Updated comments for consistency
- Public API provides the unified interface

## Testing Strategy

### Unit Tests Created

1. **test_naming_unification.cc**
   - Tests DistributionType converter with all variants
   - Tests WorldContainer query methods (old and new)
   - Tests WorldContainer replication methods
   - Validates deprecated aliases work correctly

2. **test_algorithm_converter.cc**
   - Tests Exchange::Algorithm converter with all variants
   - Tests deprecated from_string wrapper
   - Tests free function and wrapper classes
   - Tests user-defined literal
   - Validates exception on invalid input

### Syntax Validation

Created standalone test programs:
- `/tmp/test_syntax_worlddc.cc` - Validates worlddc.h changes
- `/tmp/test_syntax_exchange.cc` - Validates SCFOperators.h changes

Both compiled successfully with g++ -std=c++17.

## Documentation

### Files Created/Updated

1. **REPLICATION_API_MIGRATION.md** - Comprehensive user guide
   - Migration strategy
   - Code examples
   - Terminology clarification
   - Troubleshooting

2. **parallel_runtime.dox** - Developer documentation
   - Added note about API changes
   - References migration guide

3. **Test files** - Inline documentation
   - Usage examples
   - Edge case testing

## Known Usages of Deprecated APIs

Found in existing codebase:
- `src/madness/mra/test6.cc` - Uses `replicate()`
- `src/madness/mra/funcimpl.h` - Uses `replicate()`, `replicate_on_hosts()`, `is_replicated()`, `is_host_replicated()`
- `src/madness/mra/macrotaskq.h` - Uses `replicate()`
- `src/madness/world/cloud.h` - Uses `replicate_on_hosts()`

**Impact:** These will generate deprecation warnings but continue to work.

## Security Analysis

**CodeQL Scan:** No issues detected

**Considerations:**
- String converters throw `std::invalid_argument` on invalid input
- No buffer overflows (using std::string)
- No injection vulnerabilities
- Input validation via enum matching

## Future Work

### Phase 1 (Current Release)
- [x] Add canonical APIs
- [x] Add deprecated aliases
- [x] Add string converters
- [x] Add tests
- [x] Add documentation

### Phase 2 (Next Release)
- [ ] Migrate internal usage to canonical names
- [ ] Update examples and tutorials
- [ ] Consider adding similar converters for other enums

### Phase 3 (Future Major Release)
- [ ] Remove deprecated aliases
- [ ] Update ABI version
- [ ] Final cleanup

## Lessons Learned

1. **Class-specific converter names** prevent overload ambiguity issues
2. **Comprehensive normalization** improves user experience
3. **[[deprecated]] attribute** provides clean migration path
4. **Wrapper classes** enable implicit conversion where needed
5. **User-defined literals** offer convenience for common cases

## References

- Problem statement: GitHub issue #XXX
- WorldContainer API: `src/madness/world/worlddc.h`
- Exchange operator: `src/madness/chem/SCFOperators.h`
- Migration guide: `doc/REPLICATION_API_MIGRATION.md`
