# 16 — molresponse_v3 Master Architecture (release target)

Status: **design, drafted 2026-06-08.** The single document that defines the
complete v3 architecture for the MADNESS Release. Supersedes nothing — it
*assembles* the agreed pieces (docs 05/06/07/12/13/14/15) into one layer cake,
marks what already exists vs what is new, and defines the build sequence and the
public contracts everything fills. Read this after `00_status.md`; read the
referenced numbered docs for the detail of each layer.

## Purpose & scope

molresponse_v3 is a **unified response solver** for **mixed linear and
non-linear response properties** (α, β, vibrational/resonance Raman, 2PA, …),
built to **scale to large molecular systems** and to emit **research-quality,
exportable data** — for visualization and for **multiresolution ML models
targeting excited states**. It will be the response engine behind
`madqc --wf=response`. See `project_v3_vision` (memory) for the standing goals.

Near-term deliverables that this architecture must serve:
1. Benchmark **Raman, 2PA, and resonant-Raman** data.
2. **MRA-fidelity ES data export** for ML.
3. **Large-system scaling** (the C6H6 / naphthalene OOM wall in the root
   `CLAUDE.md`).
4. Clean **madqc integration** (replace the v2 `MolresponseLib.hpp` path).

## Design stance — re-layer + complete in place

The *solver core* already exists and is tested. The "ground-up rebuild" is
**(a)** formalizing that core into clean, contract-driven layers, and **(b)**
building the missing top layers properly, designed up front. **No throwaway.**
v2 (`src/madness/chem/MolresponseLib.hpp`) stays the production path and the
parity reference until v3 replaces it.

---

## The layer cake

```
 L5  madqc integration + benchmark/validation harness        [NEW + cm.sh exists]
 L4  Export / Viz / ML  (Tier-B sibling: post-processing)    [NEW]
 L3  Observability  (timing · diagnostics)                    [NEW; per-point metrics exist]
 L2  State-parallel  (memory-driven subworld sizing)          [NEW — doc 15 · 15c/15d]
 L1  Public contract + Orchestrator  (run_response seam)      [NEW — doc 05]
 ─────────────────────────────────────────────────────────────────────────────
 L0  Foundations (EXIST, tested)
     types/kernels · solvers · persistence · scheduler · properties · metrics
```

The spine is **L1's contract**: a single `ResponseWorkflowOutput` whose slots
(`properties`, `metadata`, `timing`, `diagnostics`, `exports`) are *filled* by
L2–L4. Designing that schema first means every later layer has a defined place
to write, and madqc/test/Python all consume one object.

### L0 — Foundations (exist; this is what we re-layer, not rewrite)

| Concern | Files | Doc |
|---|---|---|
| Typed solvers (no runtime `calc_type` branch) | `kernels/*`, `solvers/{fd_solver,es_solver,iterate_protocol,convergence_policy}.hpp` | 12 |
| Persistence + restart (just hardened, 2b) | `solvers/{response_metadata,fd_save_load,es_save_load,vbc_save_load}.hpp` | 13 |
| Scheduler DAG / reconcile / two-phase | `calc/{calc_manager,calc_executor}.hpp` | 15 (15a/15b) |
| Properties α/β/Raman (Tier-A) | `kernels/{beta,vbc,assembly}.hpp`, `ResponsePropertyPlanner.hpp` | 14 |
| Per-state metrics (coeffs/bytes/rss/iters) | `solvers/state_metrics.hpp` | 15 |

Re-layering work here is **consolidation, not redesign**: remove the
top-level/`solvers/` duplicate solver headers (`./FDSolver.hpp` vs
`solvers/fd_solver.hpp` etc.), settle one home per concern, and make `main.cpp`
stop being an ad-hoc orchestrator (it becomes an L1 Input builder).

---

## L1 — Public contract + Orchestrator (doc 05)

The one entry point. Pure of madqc/SCF knowledge: callers build the Input;
the orchestrator returns a self-contained Output.

```cpp
struct ResponseWorkflowInput {
  Molecule            molecule;
  std::string         archive_file;     // SCF checkpoint (ground state)
  std::vector<double> protocols;        // truncation-threshold ladder
  ResponseParameters  response_params;  // perturbations, freqs, properties,
                                        // ES knobs, parallel knobs, export knobs
};

struct ResponseWorkflowOutput {
  nlohmann::json properties;   // Tier-A Cartesian tensors + provenance (per protocol_key)
  nlohmann::json metadata;     // response_metadata.json mirror (state status/restart)
  nlohmann::json timing;       // L3: 3-level structured timing
  nlohmann::json diagnostics;  // L3: convergence trajectories, divergence, mem-HWM, schedule
  nlohmann::json exports;      // L4: manifest of viz/ML artifacts written (paths + descriptors)
  nlohmann::json debug_log;    // optional iteration-level
};

ResponseWorkflowOutput run_response(World&, const ResponseWorkflowInput&);
```

The orchestrator is **thin** (doc 06 Stage A shape) and owns the Output +
stage timers:

```
run_response(world, input):
  T_total, T_stage = timers
  ground   = GroundState::from_archive(...)                 # stage: load
  plan     = merge_plans(requests_from(response_params))    # stage: plan      (Tier-A planner, doc 14)
  mgr      = CalcManager(plan, calc_dir);  mgr.build(mol)
  exec     = make_executor(world, ground, input)            # serial or parallel (L2)
  mgr.run(world, exec)                                      # stage: solve     (existing loop)
  out.properties   = assemble_properties(...)               # stage: properties (Tier-A assembly)
  out.exports      = ExportManager(...).run(...)            # stage: export    (L4, off critical path)
  out.timing       = collect_timing();  out.diagnostics = collect_diagnostics()  # L3
  return out
```

`main.cpp` and the test runner become **Input builders** that call
`run_response`; so does madqc (L5). This is doc-07 Increment 10's groundwork and
the single seam that makes the engine testable, scriptable, and Python-bindable.

**Contract:** the orchestrator contains *no* scheduling, ownership, or
file-format logic — those live in CalcManager (L2) and the persistence layer.
It only sequences stages and assembles the Output.

---

## L2 — State-parallel (doc 15 · 15c/15d; doc 06 Stage B)

The large-system scaling layer. **Memory-driven, protocol-dependent** subworld
sizing — not a fixed group count.

- `groups(P) = clamp( node_mem_budget / mem_per_task(state, P), 1, n_natural )`.
- `mem_per_task(state, P_high)` is **predicted from a measured low-protocol
  solve** (15d): the seed-state solve gives the *actual* adaptive leaf structure
  for this molecule; scale by `k³ × leaf-growth(P_low→P_high)` (uses
  `state_metrics`). → a clean **pre-flight abort** instead of exit-137.
- Mechanism: a `ParallelExecutor` implementing `ICalcExecutor`. CalcManager.run
  already runs **one wave per pass**; the parallel layer fans *that one wave*
  across subworlds (`with_subworld`), each subworld runs the existing serial
  `CalcExecutor` locally, then rank 0 merges metadata shards. **Serial executor =
  `groups(P)=1`** — same code path, no special-casing.
- **Parity contract (doc 06):** serial vs parallel results must match within
  tolerance; any drift is a parallel-layer bug, not the solver.

Increment split (this is **R5, the last layer**): **15d (memory estimate +
abort) first** — fed by the R4 diagnostic study's measured per-state footprints,
and de-risks every large run — then **15c (subworld sizing + ParallelExecutor)**.

Open mechanism decisions (deferred, below): v2 `MacroTaskQ` subgroups vs raw
`with_subworld`; node-local **shared ground-orbital replica** (MPI shared-memory
windows) to cut the `n_occ` replication that is the measured OOM driver.

---

## L3 — Observability (timing · diagnostics)

What makes the scaling work *measurable* and the code *research-quality*. Must
exist **before/with** L2 (it feeds the memory model).

### Timing — three levels (doc 05)

| Level | Where recorded | Status |
|---|---|---|
| **Stage** (load/plan/solve[+per-protocol]/properties/export/total, wall+cpu) | orchestrator timers → `Output.timing` | NEW |
| **Point** (per `(state, protocol, freq)` wall_s, + coeffs/rss/iters) | CalcManager owns the `(node,protocol)` loop → record `wall_s` alongside `state_metrics` | partial (metrics exist; add wall_s) |
| **Operation** (density/coulomb/exchange/bsh/kain per iter) | `iterate_protocol` `post_step` / a `TimedScope`; **debug-gated** | NEW (optional) |

### Diagnostics — one structured record (→ `Output.diagnostics` + metadata)

- **Convergence trajectory** per point (bsh/density residual history) — solver
  `State` already carries `last_*_residual`; persist a compact history.
- **Divergence / stall** flags — exist (`diverged`, `accept_at_maxiter`).
- **Memory HWM** per `(rank, protocol)` — `state_metrics.rss_gb`; emit
  `MEMORY_HWM rank=N protocol=K rss_gb=X` (root `CLAUDE.md` ask).
- **Scheduler trace** — reconcile action per `(node,protocol)`, wave structure,
  gate holds (the run already prints these ad hoc — formalize as JSON).
- **Machine-readable protocol lines** — `PROTOCOL_START index=N thresh=X k=N` /
  `PROTOCOL_DONE index=N iters=N converged=true` (root `CLAUDE.md` ask;
  formalize the existing `PROTOCOL_SET`).

All collective (norms/sizes/reduces) stay on **all ranks**; only printing is
rank-0 gated (standing contract).

---

## L4 — Export / Viz / ML

A **Tier-B sibling** (doc 14): post-processing that **consumes already-converged
states from persistence and NEVER triggers a solve**. Runs off the critical
path, after assembly. One `ExportManager` with pluggable writers, driven by
`response_params` export knobs, keyed by the same **state identity / protocol_key**
as persistence (provenance), recording every artifact in `Output.exports`.

### What objects get exported (the physics)

- **Linear response:** response density ρ⁽¹⁾ and response orbitals x (and y) per
  `(perturbation, ω, protocol)`.
- **Excited states (closed-shell focus):** ES eigenvectors X (and Y), **transition
  densities ρ₀ₙ**, optionally NTOs — the **ML priority** (multiresolution ES data).
- **Ground:** occupied orbitals (ML reference).

### The four output targets (all on MADNESS primitives)

| Target | Mechanism | Consumer | Fidelity |
|---|---|---|---|
| **A. Native MRA coefficients** | `Function::save` (parallel archive) + descriptor manifest | multiresolution ES-ML (the differentiator) | full |
| **B. Sampled grids → numpy/json** | regular-grid sampler / `plot_line` → `.npy`/`.json` | gecko + Python ML pipelines | lossy |
| **C. Volumetric cube/DX/VTK** | `plotdx` / `plotvtk_*` (+ cube) | VMD / ParaView | viz |
| **D. Structured property/state JSON** | generalize `es_analysis.hpp` (ω, transition moments, osc. strengths, density descriptors: multipoles, ⟨r²⟩) | analysis + ML features | scalar/tensor |

**Priority for the near-term ML goal:** Target **A** on ES transition
densities/vectors, plus Target **D** (generalize the already-built—but
disabled—`es_analysis.hpp`). C and B follow for viz/baseline pipelines.

**Contract:** an export descriptor records `{object, state_id, omega,
protocol_key, target, path, grid_spec?}` so every artifact is reproducible and
joinable to the metadata. Exports are *additive* — turning them off changes
nothing in the solve.

---

## L5 — madqc integration + harness

- **madqc:** an analog of v2's `ResponseApplication<molresponse_lib>`
  (`WorkflowBuilders.hpp`) that builds `ResponseWorkflowInput` from
  `CalculationParameters` + `SCF` and serializes `ResponseWorkflowOutput` into
  `<prefix>.calc_info.json`. v3's Output is a strict superset of v2's `Results`
  (adds `timing`/`diagnostics`/`exports` as first-class). doc-07 Increment 10.
- **Harness (exists):** `cm.sh` + the results ledger. Extend with a
  **timing-comparison** and **serial/parallel parity** path (doc 05's
  benchmark-skill sketch) once L2/L3 land.

---

## Build sequence (re-layer + complete)

Dependency-ordered; each increment is independently testable via `cm.sh` on the
allocation and recorded in the ledger. **State-parallel (L2) is deliberately
last** — built on a serial pipeline that is already validated, instrumented (L3),
and driven through madqc, with a real per-state memory model from the diagnostic
study (R4) to size subworlds against.

- **R0a — Contract + thin orchestrator. ✅ DONE (2026-06-08).**
  `orchestrator/response_workflow.hpp`: `ResponseWorkflowInput` /
  `ResponseWorkflowOutput` + `run_response()` wrapping today's flow (ground load →
  CalcManager → assemble → Output); `ExecutorContext` split into plain
  `ExecutorSettings` + `ExecutorContext : ExecutorSettings` (World-bound). New
  clean driver `tests/test_run_response.cpp`. **Validated:** h2o α bit-identical
  to the cm_resume baseline (7.903656 / 9.191244 / 8.532639). Two scoped
  deviations from the sketch above: (1) `Input` carries a pre-built `ResponsePlan`
  (not `requests`) — requests-lowering moves inside `run_response` in R0b/R3 once
  the TDA/Full choice is a request field; (2) a *new* driver exercises the seam
  rather than rewriting the entangled `test_calc_manager_run` (which has dead
  analyze branches + double gs-loads) — that migration is R0b.
- **R0b — main.cpp onto the seam + L0 consolidation. ✅ DONE (2026-06-09).**
  `main.cpp` rewritten as a thin `run_response` app (the installed `molresponse_v3`
  binary now drives the seam; `--archive` CLI, legacy input-file/`fd_solve` path
  dropped — nothing invoked it). Retired the legacy `test_v3_solver`/`test_es_solver`
  binaries and removed the now-dead top-level `FDSolver.hpp`/`ESSolver.hpp`. Kept the
  **live** ES-guess chain (`ESSolverGuess.hpp` → `ResponseFunctions`/`ResponseKernel`,
  used by `build_response_ground_state` + the skeletons). Correction: `PrintLevel`
  was already single-homed in `kernels/tags.hpp` (no duplication). **Not migrated:**
  `test_calc_manager_run` stays as the low-level CalcManager test (seam already
  proven by `test_run_response`); migrating the entangled driver wasn't worth the
  risk. Validated: full v3 build + h2 app smoke (α_zz=6.449).
- **R1 — Observability.** Stage + point(`wall_s`) timing → `Output.timing`;
  diagnostics record (convergence/divergence/mem-HWM/schedule) → `Output.diagnostics`;
  formalize `PROTOCOL_START/DONE`. Makes everything below measurable.
- **R2 — Export/Viz/ML.** `ExportManager` + writers; **A (ES transition densities)
  + D (generalize `es_analysis`)** first, then C, then B. Off critical path.
  (Re-enabling `es_analysis` requires the ES-bundle load-path heap-OOB fix —
  tracked in `project_v3_es_analysis_parked`.)
- **R3 — madqc integration.** `ResponseApplication<v3>`; build
  `ResponseWorkflowInput` from `CalculationParameters`+`SCF`; serialize Output →
  `<prefix>.calc_info.json`; serial-vs-v2 parity. **Then route `cm.sh` through
  madqc** so the harness exercises the production entry point uniformly.
- **R4 — Diagnostic study.** Using madqc-driven runs + L3 observability,
  characterize timing / memory(`rss`,`coeffs`) / convergence / size-scaling across
  the molecule set (h2o → c2h4 → ch3oh → c6h6 …). Output: the empirical per-state
  `mem_per_task(state,P)` model and timing breakdowns — i.e. the **input to L2's
  subworld sizing and pre-flight abort**.
- **R5 — State-parallel (L2, last).** **15d** (memory estimate + pre-flight abort,
  fed by R4) → **15c** (subworld sizing + `ParallelExecutor` fanning one
  CalcManager wave across subworlds); serial/parallel parity test. Optionally the
  node-local shared ground-orbital replica (deferred decision below).
- **Ongoing — properties.** Finish **2PA** (`TPA_DESIGN.md`) and **resonant
  Raman** on the β/VBC core (the benchmark deliverables); they slot into Tier-A
  assembly + the L4 export.

R0 and R1 are the "base design" to settle first; R2 (export) can proceed in
parallel with R3/R4 since all only consume the contract. The chain
**R1 → R3 → R4 → R5** is the spine of the scaling effort: instrument, ship through
madqc, measure, then parallelize against measured data.

---

## Cross-cutting contracts (do not regress)

- **Collective ops on all ranks**; printing only under `if (world.rank()==0)`.
- **Closed-shell ES only**; open-shell ES + Full-RPA open-shell out of scope.
  ES/diffuse references use **d-aug-cc-pVQZ**.
- **Never write `response_metadata.json` directly** — go through the metadata
  layer.
- **Exports never trigger a solve** and are keyed by state identity / protocol_key.
- **Propose a diff + approval** before non-trivial solver/runtime edits.
- Commit with `git -c core.hooksPath=/dev/null`; branch `molresponse-feature-next`.

## Open design decisions (to settle as each layer is reached)

1. **Subworld mechanism (L2):** v2 `MacroTaskQ` subgroups vs raw `with_subworld`;
   and whether to add a **node-local shared ground-orbital replica** (MPI
   shared-memory windows) to cut the `n_occ` replication OOM driver.
2. **Export defaults (L4):** which states export by default, and the default
   grid spec (bbox/npts) for Targets B/C.
3. **Output.metadata (L1):** embed the full `response_metadata.json` or reference
   it by path (size vs self-containment).
4. **Complex/damped response:** scope for damped 2PA (resonance-divergent
   undamped vs complex damped) — `TPA_DESIGN.md` open question.
5. **Operation-level timing overhead (L3):** always-on lightweight vs debug-gated.
