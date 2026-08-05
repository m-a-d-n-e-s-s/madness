# madqc — architecture roadmap: the comprehensive QC driver

Distilled from the 2026-07 pre-release code review of the workflow/chaining and
I/O layers. This is the DESIGN for where madqc goes after the release; the
review's correctness fixes (rank-0/atomic `calc_info.json`, failure records,
v3-shaped `.out` summary, per-row property accuracy, the `response.hdf5` deck
knob + io provenance) already landed.

## Where we are (as-built)

- `qcapp::Workflow` is a **linear list of single points**: each `Driver` runs in
  a positional `task_<i>/` dir; its `summary()` JSON is appended to `tasks[]`.
  Data *can* now flow through the Driver interface: change 1 landed, and
  `Workflow::run` threads one `StepContext` task-to-task (`Drivers.hpp:124`).
- Migration of the consumers onto that context is **partial — one of four**.
  `ResponseApplication::consume_context` (`Applications.hpp:417`) takes the
  upstream reference and geometry from the context; CC2, TDHF and OEP still take
  their ground state from the **build-time side-channel**, where the builder
  captures `reference->calc()` (a live `shared_ptr<SCF/Nemo>`) before adding the
  downstream driver. See the deferred finding at the bottom.
- Geometry optimization WORKS today but only *inside* one SCF driver
  (`dft gopt=1` → `moldft_lib` → `MolOpt::optimize_app`; writes `_opt.xyz`,
  emits `optimization_results`). The `optimize` WORKFLOW is a stub (enum +
  disabled dispatch + fully-commented `OptimizeDriver` + a dead `SCFTarget`
  adapter whose result schema doesn't match `SCFApplication`).
- Restart machinery is strong and reusable: SCF checkpoint validation with a
  geometry-match guard; response restart via `response_metadata.json` +
  best-usable-seed selection; per-task `PathManager` dirs (absolute paths).

## The five changes (ordered by enablement)

### 1. StepContext dataflow (M) — the enabling change
**Status: mechanism LANDED; consumer migration partial (response only).** The
struct is `Applications.hpp:31`, the threading is `Drivers.hpp:124-129`, and
`SinglePointDriver::execute` runs consume → run → publish (`Drivers.hpp:66-69`).
What is left is moving CC2/TDHF/OEP off the build-time capture — see the deferred
finding at the bottom. Nothing downstream in this roadmap is blocked on that.

Extend `Driver::execute(workdir)` to `execute(workdir, StepContext&)` where
`StepContext` carries named artifacts: `molecule` (possibly optimized),
`reference` (`shared_ptr<SCF/Nemo>`), `archives` (name → path), and free-form
JSON. `Workflow::run` threads it task-to-task. Replaces the hidden shared_ptr
capture (which stays as a compatibility path until builders migrate).
*Acceptance:* SCF→response runs unchanged, but the response driver obtains the
reference from the context; a third stage can consume stage 2's outputs.

### 2. First-class `optimize` workflow (M)
Unblocked: the one bullet here that needs change 1 (publishing the molecule) has
the mechanism it needs, and the other two are independent of it.

Revive `OptimizeDriver` on top of `MolOpt` (the working optimizer):
- fix `SCFTarget` to read the real result schema (`SCFTargetAdapter.hpp:32-33`
  reads `res.at("energy")` / `res.at("gradient")`; neither key exists at that
  level). The right keys are `results["properties"]["energy"]` and
  `results["properties"]["gradient"]` — **not** `scf_total_energy`, which the
  plain `scf` path leaves at 0.0 and only the nemo/mp2/cc2 engines fill. The
  gradient key is present only when `derivatives` is on.
- enable `WorkflowKind::Optimize` dispatch + add `"optimize"` to
  `runnable_workflows` (its `std::array` size must go 7 → 8); consume the
  (currently unread) `OptimizationParameters`. That group is already in the
  `Params` tuple, but carries two latent bugs: `geometry_tolerence` is declared
  `bool` yet the commented `OptimizeDriver` passes it as MolOpt's `double gtol`
  (`Drivers.hpp:238`), and `get_method()` reads a `"method"` key that is never
  initialized. It is also the natural home for the geometry tolerances, which
  today are split: `MolOpt` gets values derived from `protocol().back()`
  (`Applications.hpp:890-899`) while `valid()` tests the stored optimization
  against the deck's own `gtol` (`:777`) — one knob, two meanings.
- publish the optimized `Molecule` into the StepContext. Already true for the
  in-SCF `gopt` path: `moldft_lib::run` puts the optimized geometry in
  `scf_molecule`, and `publish_to_context` forwards it as `ctx.molecule`.
*Acceptance:* `--wf=optimize` optimizes and a following stage (SCF or response)
runs AT the optimized geometry — the geopt→Raman pipeline's first half. Assert
against an SCF or response follow-on, not cc2/cis/oep: those builders call
`reference->calc()` at assembly time, which constructs the engine and freezes the
geometry into `mad.in` before any task runs, so they cannot honour an upstream
molecule until construction is deferred.

### 3. Properties at every optimization step (M)
Add a post-accepted-step hook to `MolOpt::optimize_app` (callback carrying the
current geometry + SCF state). `OptimizeDriver` routes each hook invocation to
an optional list of per-step property drivers, writing under
`task_<i>/step_<k>/`. *Acceptance:* geometry optimization emitting α at every
accepted geometry, all recorded in `calc_info.json` with step indices.

### 4. Finite-difference displacement driver (L)
`DisplacementDriver`: given a base geometry + displacement spec (atoms, axes,
±h), generate displaced `Molecule`s, run a sub-workflow (SCF+response) per
displacement in `task_<i>/disp_<k>/` (reuse the checkpoint geometry guard for
restart), gather results, and run an FD-combine step (e.g. dα/dR → full-tensor
Raman — the validation target is the archived 81-component v2 tensor in
`madness_studies/refs/raman_v2_h2o_full_tensor.json`). Needs per-step
geometry/parameter overrides in the StepContext (change 1) — `Params` itself
stays shared. *Acceptance:* full-tensor h2o Raman from displacements matching
the archived v2 reference.

### 5. Parallel HDF5 (L, design decision first)
Today: native parallelism is nio-routed multi-file POSIX; the HDF5 restart
backend is a rank-0 gather at nio=1, covering closed-shell response states
only (no SCF path). A real parallel story needs:
- MPI-IO file access (`H5Pset_fapl_mpio`) + hyperslab-partitioned datasets;
- coverage: open-shell + ES/VBC states + (decision) the SCF `.restartdata`;
- one deck-level I/O block (`io { backend hdf5; nio N; }`) governing all apps,
  stamped into `calc_info.json` (the response side already stamps `io`).

## Shared primitives to reuse (don't rebuild)
- `MolOpt` / `molecular_optimizer` — working optimizers.
- `PathManager` (absolute per-task dirs) + `ScopedCWD`.
- `checkpoint_geometry_matches` — the geometry guard; also the natural basis
  for ground-state-archive fingerprinting in response restart metadata.
- `best_usable_fd_source_key` — single source of truth for restart selection.
- The subworld fan-out layer (doc 32) — `DisplacementDriver`'s displaced solves
  are exactly the independent-states pattern it distributes.

## CC2 / OEP / TDHF chaining — the deferred review finding

The pre-release review flagged that CC2/TDHF/OEP **chaining** (feeding one method's
result into the next) still rides the build-time `shared_ptr` side-channel rather
than a typed interface. This is deliberately **out of the response path** and is
**design work, not a response-release bug**:

- **The response path is immune.** The response driver reloads its ground state
  from the on-disk archive via the adapter, so it never depends on the fragile
  in-memory handoff; a response run is unaffected regardless of chaining.
- **The fix is StepContext (change 1 above), applied to CC2/OEP.** CC2, OEP, and
  TDHF would migrate off the hidden `reference->calc()` capture exactly as the
  response driver does: obtain the `reference`/`molecule` from the `StepContext`
  the workflow threads task-to-task, and publish their own outputs into it for the
  next stage. No new mechanism is needed — StepContext already carries the named
  artifacts (`molecule`, `reference`, `archives`) these methods hand off, and is
  live in the tree; only these three consumers are unmigrated.
- **It is more than adding `consume_context` overrides.** Their builders take the
  reference by calling `reference->calc()` while assembling the workflow
  (`WorkflowBuilders.hpp:99,108,126`), and that call constructs the engine, which
  freezes the parameters by writing and re-reading `mad.in`. A context value read
  later cannot affect an engine already built, so the migration has to defer
  engine construction to `run()`.
- **Scope.** Implementing this belongs to the madqc chaining workstream (this
  roadmap), not the molresponse response release. It is recorded here as the
  design that closes the finding; the code lands when the chaining workstream is
  taken up.
