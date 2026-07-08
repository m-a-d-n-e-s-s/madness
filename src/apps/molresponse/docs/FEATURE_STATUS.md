# molresponse / madqc — feature status

The living status of the modern MADNESS response stack (`molresponse` solver +
`madqc` driver). This is the **user-facing** source of truth for what is ready,
what is opt-in/experimental, and what is on the roadmap. Update it as features
graduate — bump the badge, flip the default, add the doc link.

**Legend**
- ✅ **Stable** — validated, on by default (or a stable opt-in flag); safe for production.
- 🟡 **Experimental** — works, opt-in and OFF by default; interface/behavior may change.
- 🔬 **In development** — not in this release; tracked for a follow-up.

> **The v2 engine is removed.** `molresponse_v3` is *the* response engine
> (`engine v3` is the default; `engine v2` in an old deck fails with a
> migration hint — see `MIGRATION_FROM_V2.md`). Pre-removal reference data
> (v2↔v3 parity ≤2e-5 rel on α; the full 81-component v2 Raman tensor,
> re-verified to 1.6e-08) is archived in `madness_studies/refs/`.

## Features

| Feature | Status | Enable | Why you'd use it | Guide |
|---|---|---|---|---|
| Response solver (α, β) | ✅ Stable | default (`engine v3`) | Modern unified solver; α/β validated (H₂O α_zz≈8.53, β_zzz≈7.76 vs refs; v2 parity ≤2e-5), restart-safe, structured diagnostics | Getting started with molresponse |
| madqc orchestration | ✅ Stable | default | One SCF→response driver, one input deck; `.out` human summary + `calc_info.json` (rank-0, atomic, failure records) | `madqc/README.md`, `MIGRATION_FROM_V2.md` |
| Protocol ladder (honest-climb) | ✅ Stable | automatic | Coarse rungs seed finer ones even when unconverged (budget-exhausted rungs advance; `converged` stays honest); properties carry per-row accuracy of the states that built them | Protocol & convergence guide |
| State-parallel FD | ✅ Stable | `--fd-subworlds=P` (per node; 0 = single-world) | 2–2.56× faster + ~8× less per-rank memory; subworld α matches canonical to 3.5e-7 | Scaling FD across nodes |
| FD exchange tensor | ✅ Stable (opt-in) | `--fd-tensor` (FD), `--es-tensor` (ES γ) | Shared convolution tensors + per-protocol g0 cache + tiling (bounds n² memory; required at n≳34) | Exchange tensor path, `operator_contracts.md` |
| Single-component Raman | ✅ Stable (validated) | `nuclear` block, one (atom, axis) | β(dipole; dipole, nuclear) — validated vs the archived v2 value (0.13%) | RESPONSE_PROPERTIES.md |
| HDF5 restart / I/O | 🟡 Experimental (opt-in, closed-shell) | deck `response.hdf5 true` / `--hdf5` / env `MADRESPONSE_IO_HDF5` (build `-DMADNESS_ENABLE_HDF5`) | 15–40% faster restart reads; auto-detecting readers; effective backend stamped into metadata (`io` block) | Accelerating restart with HDF5 |
| MRA visualization export | ✅ Stable (NP=1) | build `-DMADNESS_ENABLE_VTK`; `dump_mra_trees --htg/--amr/--coeffs` | Inspect trees + response orbitals/ρ⁽¹⁾ in ParaView; native `.mad.h5` coefficient archive | Visualizing MRA trees & response orbitals |
| Performance cost model | 🟡 Experimental | build `-DENABLE_WORLD_PROFILE`; `perf_model_fit.py` | Per-phase profile (incl. exchange meters) → predict wall-time / find bottlenecks | Profiling & the cost model |
| State-parallel auto-select | 🔬 In development | — (manual `--fd-subworlds` is stable) | Pick subworld count from problem shape + nodes (needs weak-scale calibration) | (roadmap) |

## Response properties

| Property | Status | Notes |
|---|---|---|
| Polarizability α (static + dynamic) | ✅ Stable | per-row accuracy (`row_accuracy`) recorded with every tensor |
| First hyperpolarizability β | ✅ Stable | 2n+1 / VBC; rows carry `row_accuracy` (FD inputs) + `vbc_accuracy` |
| Single-component Raman | ✅ Stable | validated vs archived v2 reference (0.13%); rows carry accuracy like β |
| Excited states (TDA / Full RPA) | 🟡 Experimental | eigenpairs compute; transition analysis (oscillator strengths) always emitted for converged bundles; protocol climb stabilizing |
| Full-tensor per-atom Raman | 🔬 In development | v2-only capability at removal; the 81-component v2 reference tensor is archived (`refs/raman_v2_h2o_full_tensor.json`) as the validation target |
| Two-photon absorption (2PA) | 🔬 In development | kernel designed; gated on ES convergence |
| Resonance Raman | 🔬 In development | gated on ES convergence |

## Known limits (release notes)

- **HDF5**: closed-shell response states only; the structured (interop) writer is
  NP=1; the restart backend gathers at nio=1 (no parallel HDF5 yet — see roadmap).
- **Large-molecule dynamics on small subworlds**: 20+-occupied Full solves at k8+
  can exceed the hang-watchdog on 2-rank subworlds (`MAD_WAIT_TIMEOUT` raises it);
  prefer more ranks per subworld or the single-world path for the finest rung.
- **Visualization export** is single-rank (NP=1); production restart I/O is not.

## Roadmap — madqc as the comprehensive QC driver

From the 2026-07 architecture review (`madqc/ARCHITECTURE_ROADMAP.md`):

1. **StepContext dataflow** — drivers pass geometry/reference/archives through
   the workflow interface instead of build-time capture; enables 3+-stage pipelines.
2. **First-class geometry optimization** (`optimize` workflow) — revive
   `OptimizeDriver` on top of the working `MolOpt`/`gopt` machinery.
3. **Properties at every optimization step** — per-accepted-geometry property hook.
4. **Finite-difference higher-order properties** — displacement fan-out driver +
   per-step geometry overrides + FD combine.
5. **Parallel HDF5** — MPI-IO file access + partitioned datasets; one I/O story
   for SCF + response state.

Also tracked: state-parallel auto-selector + multi-node-per-state; ES
state-parallel; full-tensor Raman via the displacement/state-parallel machinery;
ground-state archive fingerprinting in restart metadata (prevents mixed-GS dirs).
