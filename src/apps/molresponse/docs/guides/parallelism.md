# Parallelism & subworlds

How molresponse scales a response calculation across MPI ranks and nodes. This is
a design/architecture reference; the recipes for launching runs are in
[`madqc/RESPONSE_PROPERTIES.md`](../../../madqc/RESPONSE_PROPERTIES.md).

## The shape of the problem

A response calculation is **embarrassingly parallel across independent solves**:
different excited-state roots, different perturbation directions, and different
frequencies are almost entirely independent of one another. Each *individual*
solve, by contrast, is a sequence of MRA operations (BSH applications, the
two-electron kernel) that MADNESS already threads within a World.

The key measured finding that sets the architecture: **spreading one solve across
more ranks walls quickly.** Taking a single excited-state bundle from 8 → 32 ranks
made the per-iteration cost *worse* (≈18% parallel efficiency) — the intra-solve
work does not strong-scale past a node. So the design does the opposite: **one
solve per World, many Worlds**, and parallelizes over the independent work items.

## The two-axis decomposition

1. **Within a World** — the existing MADNESS thread pool + rank set runs a single
   solve (one root, or one perturbation@frequency). This is where BSH/exchange
   threading lives; it is left as-is.
2. **Across Worlds** — the independent solves are distributed over **subworlds**
   via the MacroTask queue (`MacroTaskQ`). Each subworld pulls a work item
   (a root, or an FD point), solves it to convergence, and returns the result;
   the scheduler hands out the next item as subworlds free up.

This is the same MacroTask machinery the exchange operator already uses, applied
one level up — at the granularity of whole response solves.

## Subworlds and the ground state

Subworlds are created **node-aligned**: `MPI_Comm_split_type(MPI_COMM_TYPE_SHARED)`
puts every rank on a physical node into the same subworld, so a subworld is a
node (or a sub-node group). Two MADNESS realities shape this:

- **There is no shared-memory Function storage.** A subworld cannot point at
  another World's coefficient data; MADNESS has intra-node *communicators* but not
  shared coefficient *memory*.
- **A subworld must be given its data by copy-in.** The **Cloud**
  (`cloud.h`) serializes Functions across the World boundary
  (`store`/`load`), which is how the ground state reaches each subworld.

So the ground state is **Cloud-copied into each node-subworld once** and left
distributed there for every work item that subworld processes — the practical
substitute for the (impossible) single shared copy, and far cheaper than
replicating it per task.

## Status

| capability | status |
|---|---|
| MacroTask exchange (within-solve two-electron fan-out) | ✅ in use |
| node-aligned subworld creation (`Split_type` SHARED) | ✅ |
| Cloud-shared ground state into subworlds | ✅ |
| FD state-parallel (perturbation × frequency across subworlds) | 🚧 landing incrementally |
| ES root-parallel across subworlds | 🚧 design complete, incremental |

Single-World runs (one rank set, threaded) are always available and are the
default for small systems; the subworld fan-out is the path for many-root /
many-frequency / large-system runs across nodes.

## Notes on launching

On this deployment, multi-node MPI uses `srun --mpi=pmix` with ≤2 ranks per node;
`--mpi=pmi2` silently produces singletons. A single rank with threads (no
launcher) is the simplest mode and exercises the within-World path only. See the
cluster workflow notes in
[`madqc/SEAWULF_INTERACTIVE_WORKFLOW.md`](../../../madqc/SEAWULF_INTERACTIVE_WORKFLOW.md).
