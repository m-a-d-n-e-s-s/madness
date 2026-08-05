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
  `Workflow::run` threads one `StepContext` task-to-task.
- Migration of the consumers onto that context is **partial — one of four**.
  `ResponseApplication::consume_context` takes the upstream reference and
  geometry from the context; CC2, TDHF and OEP still take their ground state from
  the **build-time side-channel**, where the builder captures `reference->calc()`
  (a live `shared_ptr<SCF/Nemo>`) before adding the downstream driver. See the
  deferred finding at the bottom.
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
Extend `Driver::execute(workdir)` to `execute(workdir, StepContext&)` where
`StepContext` carries named artifacts: `molecule` (possibly optimized),
`reference` (`shared_ptr<SCF/Nemo>`), `archives` (name → path), and free-form
JSON. `Workflow::run` threads it task-to-task. Replaces the hidden shared_ptr
capture (which stays as a compatibility path until builders migrate).
*Acceptance:* SCF→response runs unchanged, but the response driver obtains the
reference from the context; a third stage can consume stage 2's outputs.

**Status: mechanism LANDED; consumer migration partial (response only).** The
struct is in `Applications.hpp`, the threading in `Drivers.hpp`, and
`SinglePointDriver::execute` runs consume → run → publish. What is left is moving
CC2/TDHF/OEP off the build-time capture — see the deferred finding at the bottom.
Nothing downstream in this roadmap is blocked on that.

### 2. First-class `optimize` workflow (M)
**Status: LANDED.** `--optimize --wf=<scf|nemo>` is one task,
`qcapp::OptimizeDriver<Library>` (`Drivers.hpp`), driving `MolOpt` over the reference
engine the Library policy supplies — the same code optimizes on moldft
(`Calc = SCF`) and nemo (`Calc = Nemo`). `--wf` names the reference method and
`--optimize` asks for its geometry to be optimized; naming the engine in the
`optimization` group as well was removed, and `--wf=optimize` now answers with a
migration message. It reads the
`optimization` group (which nothing read before), derives any threshold the deck
leaves unset from `protocol().back()`, writes `<prefix>_opt.xyz`, and publishes the
optimized geometry into the StepContext. Covered by
`src/examples/qc/{scf,nemo}_lih_optimize`; the moldft case reproduces the
in-SCF `dft gopt` geometry to the last digit (r = 3.035071 bohr), and the two
engines agree on the minimum to 1.3e-3 bohr.

Three deviations from the original plan, each deliberate:

- It is a **Driver**, not an Application. Numerical gradients are to be computed
  from *displaced sub-runs* — one sub-run per ±h Cartesian displacement, each in its
  own directory with its own `calc_info`, so displacements are restartable and can
  later be spread across subworlds — and owning sub-runs is a Driver's job. The
  gradient source is chosen at one point, `GeometryTarget` in
  `SCFTargetAdapter.hpp`: `AnalyticTarget` today, a `DisplacedEnergyTarget` later,
  which will run one `SCFApplication<Library>` per displacement under
  `task_<i>/step_<k>/disp_<...>/` and so inherit that class's checkpoint and
  `checkpoint_geometry_matches` guard for free. Its energies must come from
  `results["properties"]["energy"]` — **not** `scf_total_energy`, which the plain
  `scf` path leaves at 0.0. This is the same seam changes 3 and 4 need.
- `SCFTarget` was **replaced rather than repaired**, but its idea returns for the
  numerical path above: an Application per displaced geometry is exactly right when
  the point is a restartable sub-run, and exactly wrong for walking one engine
  along a path. So the analytic target drives the engine directly — both engines
  reach `MolOpt` through `OptimizationTargetInterface`, which `MolecularEnergy`
  (for an SCF) and `Nemo` (itself) already implement — and the per-geometry
  Application returns only where sub-runs are wanted.
- A **bug had to be fixed for nemo to work**: nothing dropped
  `converged_for_thresh` when the nuclei moved, so `Nemo::value`'s `skip_solve`
  short-circuited the SCF at every new geometry and the optimizer walked a frozen
  wavefunction — energy constant to every digit, gradients from stale orbitals. That
  affected the standalone `nemo` app's optimizer (`MolecularOptimizer`) too, not
  just this workflow.

Still open: `initial_hessian` is declined with an explicit error rather than
silently ignored; the optimizer keeps no checkpoint of its own, so an interrupted
optimization restarts from the input geometry; and the acceptance criterion's second
half — *a following stage runs AT the optimized geometry* — needs a consumer, i.e.
`SCFApplication::consume_context` adopting `ctx.molecule` before it builds its
engine. The geometry is published; nothing reads it yet.

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
  artifacts (`molecule`, `reference`, `archives`) these methods hand off.
- **Scope.** Implementing this belongs to the madqc chaining workstream (this
  roadmap), not the molresponse response release. It is recorded here as the
  design that closes the finding; the code lands when the chaining workstream is
  taken up.
