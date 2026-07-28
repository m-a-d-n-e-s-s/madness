# StepContext — passing data between workflow steps

How one madqc task hands results to the next, why the mechanism exists, and what
it takes for a method (CC2, OEP, TDHF) to participate.

## The problem it solves

A madqc workflow is a list of tasks: an SCF, then a response calculation, and in
principle more. Each runs in its own `task_<i>/` directory and contributes a
summary to the run's JSON record. But the tasks are not independent — the response
step needs the ground state the SCF produced.

Before StepContext, that data moved through a **build-time side channel**: the
workflow builder captured a live `shared_ptr<SCF>` from the reference task and
handed it to the downstream driver when the pipeline was *constructed*. That works,
but it has three problems:

- It is invisible in the `Driver` interface — nothing in a driver's signature says
  it depends on an upstream step, so the coupling is implicit.
- It only carries what the builder happened to capture (a live engine pointer), so
  a step cannot publish anything else it computed — an optimized geometry, an
  archive path, a derived quantity.
- It does not chain: a third stage cannot consume the second stage's outputs,
  because the wiring was fixed at construction for one specific pair.

## The mechanism

A single `StepContext` is created once per workflow and **threaded through every
task in order** (`Drivers.hpp`: one context, passed to each driver's `execute`).
It is a small struct of named, optional artifacts
(`src/madness/chem/Applications.hpp`):

```cpp
struct StepContext {
  std::optional<Molecule>                            molecule;   // possibly optimized/displaced
  std::shared_ptr<SCF>                               reference;  // live ground-state engine
  std::map<std::string, std::filesystem::path>       archives;   // named paths, absolute
  nlohmann::json                                     blob;       // not-yet-first-class artifacts
};
```

Every field is optional, and **unset means "no upstream step provided this"** — so
a step that runs first, or runs alone, simply sees an empty context and falls back
to its own inputs. There is no ordering requirement to satisfy and no null-checking
protocol to get wrong.

Applications interact through two hooks, both defaulting to no-ops:

```cpp
virtual void consume_context(const StepContext& ctx) {}   // read what upstream published
virtual void publish_to_context(StepContext& ctx)   {}    // publish this step's outputs
```

The driver calls `consume_context` before running the step and
`publish_to_context` after (`Drivers.hpp:66,69`). Because both default to nothing,
**an application that does not participate needs no changes at all** — which is
what made this adoptable incrementally rather than as a flag-day refactor.

## What it enables

- **Explicit, typed dependencies.** A step declares what it consumes; the coupling
  is in the code rather than in builder wiring.
- **Chains longer than two.** Any step can read what any earlier step published, so
  a third and fourth stage become possible — geometry optimization followed by a
  property calculation *at the optimized geometry*, or a displacement fan-out
  feeding finite-difference properties.
- **Publishing more than a pointer.** A step can hand forward an optimized
  molecule, an archive path, or arbitrary JSON — not just a live engine.
- **Restart-friendliness.** Publishing an *archive path* rather than an in-memory
  object means the downstream step can reload from disk, which is what makes a
  chained workflow survive being restarted.

## Current state

Phase 1 is implemented and in use: the context is threaded through the workflow,
and the response path participates — the ground-state step publishes its archive
path and the response step consumes the reference it needs. The older build-time
capture remains as a compatibility path so nothing broke during adoption.

## Adopting it in CC2 / OEP / TDHF

These methods still use the build-time `shared_ptr` capture. Migrating one is
small and mechanical, and needs no new mechanism:

1. Override `consume_context` to take `reference` (and `molecule`, if the method
   should honour an upstream optimized geometry) from the context instead of from
   a builder-captured pointer.
2. Override `publish_to_context` to publish whatever a downstream step could want
   — at minimum the archive path of the converged state, so the next step can
   reload rather than depend on live memory.
3. Leave the build-time capture in place until the migration is verified, then
   remove it.

The response path is the worked example to copy: it consumes a reference and
publishes an archive path, and it reloads from that archive rather than relying on
the in-memory handoff — which is also why the response path is unaffected by the
state of the others.

The roadmap items that build on this (a first-class `optimize` workflow, properties
at every optimization step, finite-difference displacement fan-out) are described
in [`ARCHITECTURE_ROADMAP.md`](ARCHITECTURE_ROADMAP.md); each assumes a step can
publish to and consume from the context, so this is the enabling piece.
