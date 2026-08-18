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
  std::optional<Molecule>                            molecule;        // possibly optimized/displaced
  std::shared_ptr<SCF>                               reference;       // live ground-state engine
  std::shared_ptr<Nemo>                              nemo_reference;  // live nemo ground state
  std::map<std::string, std::filesystem::path>       archives;        // named paths, absolute
  nlohmann::json                                     blob;            // not-yet-first-class artifacts
};
```

`reference` and `nemo_reference` are separate fields rather than one field of a
common type because `SCF` and `Nemo` share no base class. A producer publishes
whichever it has: `SCFApplication<moldft_lib>` publishes `reference`,
`SCFApplication<nemo_lib>` publishes `nemo_reference`. Aliasing a Nemo's inner
SCF into `reference` would silently change what the response step consumes.

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

The context is threaded through the workflow and **all four consumers read it**.
The response step adopts the reference and geometry it is handed; CC2, TDHF and OEP
verify the context against what they were constructed with but cannot adopt it
(next section). The older build-time capture remains as the mechanism those three
actually use.

## CC2 / OEP / TDHF: verify-only, and why

All four consumers now read the context, but these three can only **verify** it.
They cannot adopt a reference or a geometry published by an upstream step, and no
amount of `consume_context` code changes that, because each builds its engine from
a reference captured when the workflow is **assembled**, not when it runs:

- each inherits from both `Application` and its engine, and the engine base is
  constructed in the member-init list from that captured `shared_ptr<const Nemo>`;
- `CC2`'s constructor calls `set_protocol` on the reference and builds its own
  `TDHF` and `CCPotentials` from it — three separate captures;
- `TDHF`'s constructor freezes its derived parameters off the reference;
- `OEP`'s own `Nemo` base is constructed *from* the reference's parameters and
  molecule.

`TDHF::set_reference`, `OEP::set_reference` and `CCPotentials::reset_nemo` look
like the migration and are not: each swaps only part of that state, leaving a
silently inconsistent object. **A half-swap is worse than a loud refusal.**

So each of the three overrides `consume_context` with a single call to
`check_context_matches_reference` (`Applications.hpp`), which throws if the
threaded context disagrees with what the step was built on — a different reference
object, or an upstream geometry that differs by more than round-trip noise — and
also refuses a reference that has no orbitals. All three conditions are
unreachable in the workflows the builders can currently assemble; they exist so
that a future `optimize → cc2` chain fails loudly instead of writing out a
correlation energy computed at the wrong geometry.

**The remaining work** is to defer engine construction from the builder to
`run()`. Then these three can consume a reference like the response step does, and
change 2's `optimize → cc2` acceptance criterion becomes reachable.

The response path is the worked example: it consumes a reference and reloads the
ground state from the on-disk archive rather than relying on the in-memory
handoff — which is also why the response path is unaffected by the state of the
others.

The roadmap items that build on this (a first-class `optimize` workflow, properties
at every optimization step, finite-difference displacement fan-out) are described
in [`ARCHITECTURE_ROADMAP.md`](ARCHITECTURE_ROADMAP.md); each assumes a step can
publish to and consume from the context, so this is the enabling piece.
