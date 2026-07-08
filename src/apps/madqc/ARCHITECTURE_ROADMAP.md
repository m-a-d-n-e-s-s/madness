# madqc — architecture roadmap: the comprehensive QC driver

Distilled from the 2026-07 pre-release code review of the workflow/chaining and
I/O layers. This is the DESIGN for where madqc goes after the release; the
review's correctness fixes (rank-0/atomic `calc_info.json`, failure records,
v3-shaped `.out` summary, per-row property accuracy, the `response.hdf5` deck
knob + io provenance) already landed.

## Where we are (as-built)

- `qcapp::Workflow` is a **linear list of independent single points**: each
  `Driver` runs in a positional `task_<i>/` dir; its `summary()` JSON is
  appended to `tasks[]`. No data flows through the Driver interface.
- The working SCF→{response, cc2, cis, oep} pipelines pass data via a
  **build-time side-channel**: the builder captures `reference->calc()`
  (a live `shared_ptr<SCF/Nemo>`) before adding the downstream driver.
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

### 2. First-class `optimize` workflow (M)
Revive `OptimizeDriver` on top of `MolOpt` (the working optimizer):
- fix `SCFTarget` to read the real result schema
  (`results["scf"]["scf_total_energy"]` / gradient — not `res.at("energy")`);
- enable `WorkflowKind::Optimize` dispatch + add `"optimize"` to
  `runnable_workflows`; consume the (currently unread) `OptimizationParameters`;
- publish the optimized `Molecule` into the StepContext.
*Acceptance:* `--wf=optimize` optimizes and a following stage (SCF or response)
runs AT the optimized geometry — the geopt→Raman pipeline's first half.

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
