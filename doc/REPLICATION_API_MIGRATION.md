# Replication API Naming Migration Guide

## Overview

The replication-related APIs in MADNESS have been unified to use consistent "node" and "rank" terminology instead of the previously mixed "host"/"node" naming. This document provides guidance for migrating existing code.

## Changes Summary

### WorldContainer API Changes

#### New Canonical Method Names

**Query Methods:**
- `rank_replication()` - Check if container is replicated across all ranks (new canonical name)
- `node_replication()` - Check if container is replicated across all nodes/hosts (new canonical name)

**Replication Methods:**
- `replicate_on_ranks()` - Replicate container to all ranks (new canonical name)
- `replicate_on_nodes()` - Replicate container to all nodes/hosts (new canonical name)

#### Deprecated Method Names (Still Supported)

The following methods are deprecated but remain functional with backward compatibility:

- `is_replicated()` → Use `rank_replication()` instead
- `is_host_replicated()` → Use `node_replication()` instead
- `replicate()` → Use `replicate_on_ranks()` instead
- `replicate_on_hosts()` → Use `replicate_on_nodes()` instead

### String-to-Enum Converters

#### DistributionType Converter

New function: `distribution_type_from_string(std::string s)`

Supported string variants (case-insensitive, accepts spaces/dashes/underscores):
- **Distributed**: "distributed"
- **RankReplicated**: "rank_replicated", "rankreplicated", "rank"
- **NodeReplicated**: "node_replicated", "nodereplicated", "node", "host_replicated", "hostreplicated", "host"

Example usage:
```cpp
auto dist_type = distribution_type_from_string("rank-replicated");
// Or using the wrapper for implicit conversion:
DistributionTypeFromString wrapper("node");
DistributionType dt = wrapper;
```

#### Exchange::Algorithm Converter

New method: `Exchange<T,NDIM>::from_string_algorithm(std::string s)`

Supported string variants (case-insensitive, accepts spaces/dashes/underscores):
- **small_memory**: "small_memory", "smallmemory", "small"
- **large_memory**: "large_memory", "largememory", "large"
- **multiworld_efficient**: "multiworld_efficient", "multiworldefficient", "multiworld"
- **multiworld_efficient_row**: "multiworld_efficient_row", "multiworld_efficientrow"
- **fetch_compute**: "fetch_compute", "fetchcompute", "fetch"

Example usage:
```cpp
using AlgType = Exchange<double, 3>::Algorithm;
auto alg = Exchange<double, 3>::from_string_algorithm("small-memory");

// Or using the deprecated wrapper:
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"
auto alg2 = Exchange<double, 3>::from_string("large");

// Free function:
auto alg3 = algorithm_from_string<double, 3>("multiworld");

// User-defined literal (for double,3D):
auto alg4 = "small"_alg;

// Implicit conversion wrapper:
AlgorithmFromString<double, 3> wrapper("fetch");
AlgType alg5 = wrapper;
```

## Migration Strategy

### Phase 1: Update to New Names (Recommended)

1. Replace deprecated method calls with canonical names:
   ```cpp
   // Old code
   if (container.is_replicated()) { ... }
   container.replicate();
   
   // New code
   if (container.rank_replication()) { ... }
   container.replicate_on_ranks();
   ```

2. Update string-based configuration parsing to use new converters:
   ```cpp
   // Old code
   std::string mode = get_config("replication_mode");
   if (mode == "host") { /* replicate on hosts */ }
   
   // New code
   auto dist_type = distribution_type_from_string(get_config("replication_mode"));
   if (dist_type == NodeReplicated) { /* replicate on nodes */ }
   ```

### Phase 2: Continue Using Deprecated APIs (Temporary)

If you need time to migrate, you can suppress deprecation warnings:

```cpp
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wdeprecated-declarations"

// Use deprecated APIs here
container.replicate();
if (container.is_replicated()) { ... }

#pragma GCC diagnostic pop
```

### Phase 3: Clean Up (Future)

In a future major release, the deprecated methods may be removed. Migrate all code to use canonical names before then.

## Terminology Clarification

- **Rank**: A single MPI process (ProcessID)
- **Node/Host**: A physical compute node that may contain multiple ranks
- **Rank Replication**: Data is replicated to every rank in the world
- **Node Replication**: Data is replicated once per node (to the lowest rank on each node)

## Testing

Two new test files have been added to validate the changes:
- `src/madness/world/test_naming_unification.cc` - Tests WorldContainer naming and DistributionType converter
- `src/madness/chem/test_algorithm_converter.cc` - Tests Exchange::Algorithm converter

Run these tests to verify correct behavior:
```bash
cd build
make test_naming_unification test_algorithm_converter
./src/madness/world/test_naming_unification
./src/madness/chem/test_algorithm_converter
```

## Questions or Issues?

If you encounter any problems during migration, please file an issue on the MADNESS GitHub repository with:
1. The deprecated API you're using
2. The context in which it's used
3. Any error messages or unexpected behavior

## See Also

- WorldContainer API documentation: `src/madness/world/worlddc.h`
- Exchange operator documentation: `src/madness/chem/SCFOperators.h`
- DistributionType enum: `src/madness/world/worlddc.h` (lines 81-86)
