# molresponse / madqc — feature status

The living status of the modern MADNESS response stack (`molresponse` solver +
`madqc` driver). This is the **user-facing** source of truth for what is ready,
what is opt-in/experimental, and what is on the roadmap. Update it as features
graduate — bump the badge, flip the default, add the guide link.

This matrix states capability and how to enable each feature; it does not report
cross-code validation numbers (those are documented separately). Per-property
detail is in the [feature guides](guides/).

**Legend**
- ✅ **Stable** — validated, on by default (or a stable opt-in flag); safe for production.
- 🟡 **Experimental / preliminary** — works, opt-in and/or interface may still change.
- 🔬 **In development** — not in this release; tracked for a follow-up.

> **The v2 engine is removed.** `molresponse` is *the* response engine (`engine v3`
> is the default; `engine v2` in an old deck fails with a migration hint — see
> `MIGRATION_FROM_V2.md`). Pre-removal v2↔v3 parity reference data is archived
> outside the release tree.

## Features

| Feature | Status | Enable | Why you'd use it | Guide |
|---|---|---|---|---|
| Response solver (α, β) | ✅ Stable | default (`engine v3`) | Modern unified solver; restart-safe, structured diagnostics | [polarizability](guides/polarizability.md), [hyperpolarizability](guides/hyperpolarizability.md) |
| madqc orchestration | ✅ Stable | default | One SCF→response driver, one input deck; `.out` human summary + `calc_info.json` (rank-0, atomic, failure records) | [`madqc/README.md`](../../madqc/README.md), `MIGRATION_FROM_V2.md` |
| Protocol ladder (honest-climb) | ✅ Stable | automatic | Coarse rungs seed finer ones even when unconverged; properties carry per-row accuracy of the states that built them | [excited_states](guides/excited_states.md) |
| State-parallel FD | ✅ Stable | `--fd-subworlds=P` (per node; 0 = single-world) | Faster + much less per-rank memory via node-aligned subworlds | [parallelism](guides/parallelism.md) |
| FD exchange tensor | ✅ Stable (opt-in) | `--fd-tensor` (FD), `--es-tensor` (ES γ) | Shared convolution tensors + per-protocol g0 cache + tiling (bounds n² memory; required at n≳34) | [parallelism](guides/parallelism.md), [formalism](guides/formalism.md) |
| DALTON warm-start | ✅ Stable | `seed_from_dalton` / `seed_moldft_from_dalton` tools | Seed the excited-state solve (and SCF) from a DALTON calculation — fewer iterations and correct root selection; can skip the coarse protocol rung | [DALTON warm-start](guides/dalton_warm_start.md) |
| HDF5 restart / I/O | 🟡 Experimental (opt-in, closed-shell) | deck `response.hdf5 true` or `io { backend hdf5 }` (build `-DMADNESS_ENABLE_HDF5=ON`) | Interop/inspection with HDF5 tooling; auto-detecting readers; backend stamped per entry in metadata | [`madqc/HDF5_IO.md`](../../madqc/HDF5_IO.md) |
| MRA visualization export | ✅ Stable (NP=1) | build `-DMADNESS_ENABLE_VTK`; `dump_mra_trees --htg/--amr/--coeffs` | Inspect trees + response orbitals/ρ⁽¹⁾ in ParaView; native `.mad.h5` coefficient archive | — |

## Response properties

| Property | Status | Notes |
|---|---|---|
| Polarizability α (static + dynamic) | ✅ Stable | per-row accuracy (`row_accuracy`) recorded with every tensor; benchmarked in the published literature |
| First hyperpolarizability β (static + SHG) | ✅ Stable | 2n+1 / VBC; rows carry accuracy; benchmarked in the published literature |
| Single-component Raman | ✅ Stable | β(dipole; dipole, nuclear) for one (atom, axis); rows carry accuracy |
| Excited states (TDA / Full RPA) | ✅ Stable (closed shell) | eigenpairs + oscillator strengths; protocol climb + optional DALTON seeding |
| Two-photon absorption (2PA) | 🟡 Preliminary (closed-shell singlets) | single-residue tensor + Monson–McClain observables; working equations being finalized ([guide](guides/two_photon_absorption.md)) |
| Full-tensor per-atom Raman | 🔬 In development | per-atom/axis tensor; gated on the state-parallel layer |
| Resonance Raman | 🔬 In development | gated on the excited-state/response path |
| Open-shell response | 🔬 Out of scope | closed-shell only in this release |

## Known limits (release notes)

- **HDF5**: closed-shell response states only; the restart backend gathers at
  nio=1 (no parallel HDF5 yet — see roadmap).
- **Large-molecule dynamics on small subworlds**: 20+-occupied Full solves at k8+
  can exceed the hang-watchdog on 2-rank subworlds (`MAD_WAIT_TIMEOUT` raises it);
  prefer more ranks per subworld or the single-world path for the finest rung.
- **Visualization export** is single-rank (NP=1); production restart I/O is not.

## Roadmap — madqc as the comprehensive QC driver

From the architecture review ([`madqc/ARCHITECTURE_ROADMAP.md`](../../madqc/ARCHITECTURE_ROADMAP.md)):
StepContext dataflow (typed task-to-task handoff — also the path to close CC2/OEP
chaining), a first-class `optimize` workflow, properties at every optimization
step, finite-difference higher-order properties, and parallel HDF5. Also tracked:
state-parallel auto-selection, ES state-parallel, and full-tensor Raman.
