# molresponse_v3

A unified MADNESS solver for **molecular response properties** — linear and
non-linear — built on multiresolution analysis (MRA). It is the response engine
behind `madqc --wf=response` (opt in with `engine v3`), and is designed to scale
to larger systems and to export research-quality, MRA-fidelity data.

> **Status snapshot (2026-06): active development.** The solver core, the
> orchestration seam, observability, and madqc integration are in place and
> validated; several properties and the large-system scaling layer are still in
> flight. For the live, per-workstream picture see
> [`docs/00_status.md`](docs/00_status.md); for the architecture see
> [`docs/16_architecture.md`](docs/16_architecture.md).

## What works today

The capability matrix below is the honest current state — what is **validated
end-to-end through madqc** vs. what is **experimental / in development**. (Run
recipes and examples: [`../madqc/RESPONSE_PROPERTIES.md`](../madqc/RESPONSE_PROPERTIES.md).)

| Property | Status | Notes |
|----------|--------|-------|
| **Polarizability α** (static + frequency-dependent) | ✅ **supported** | validated vs v2 / reference (e.g. H₂O α_zz=8.53) |
| **First hyperpolarizability β** (static + dynamic SHG/OR) | ✅ **supported** | 2n+1 / VBC path (e.g. H₂O β_zzz≈7.76); validate vs published benchmark β. Dynamic runs need all VBC contributions to reach the top protocol |
| **Vibrational Raman** — single component | 🟡 **implemented; validation pending** | β(dipole; dipole, nuclear) for one (atom, axis); runs end-to-end but no reference yet — confirm the printed values |
| **Excited states** (TDA / full RPA) | 🟡 **experimental** | eigenpairs compute; convergence/climb across protocols is being stabilized |
| **2-photon absorption (2PA)** | 🔬 **in development** | kernel designed (`docs/TPA_DESIGN.md`); gated on ES convergence |
| **Resonance Raman** | 🔬 **in development** | gated on ES convergence |
| **Full-tensor Raman** | 🔬 **in development** | per-atom/axis tensor; gated on the state-parallel layer |
| **Large-system scaling** (state-parallel) | 🔬 **in development** | memory-driven subworld sizing; addresses the C₆H₆/naphthalene OOM wall |

**Scope:** closed-shell only. Open-shell excited states and open-shell full-RPA
are out of scope for v3. Excited-state / diffuse-state references use
**d-aug-cc-pVQZ** (single augmentation manufactures a phantom ~3% error).

## Using v3

**Through madqc (recommended).** Add `engine v3` to the `response` block of a
madqc input deck and run the response workflow:

```text
response
    engine              v3
    dipole              true
    dipole.directions   z
    dipole.frequencies  [0.0, 0.04]
    quadratic           true        # also compute beta
end
```
```bash
madqc --wf=response h2o.in
```

Outputs land in `<prefix>.calc_info.json` (machine-readable) and `<prefix>.out`
(human-readable). See [`../madqc/RESPONSE_PROPERTIES.md`](../madqc/RESPONSE_PROPERTIES.md)
for per-property knobs, example inputs, and the response result schema, and
[`../madqc/README.md`](../madqc/README.md) for the general madqc driver.

**Standalone binary** (development / low-level testing) — drives the
`run_response` seam directly on a ground-state archive:

```bash
molresponse_v3 --archive=<path_to_moldft.restartdata>
mpirun -n <N> molresponse_v3 --archive=<path>
```

## Relationship to the other response codes

- **`molresponse_v2`** — the current **production** path (`MolresponseLib.hpp`)
  and the **parity reference**: v3 is validated against v2 (same input →
  matching `calc_info.json` values) at every increment. v2 remains the default;
  v3 is opt-in via `engine v3`.
- **`molresponse_legacy`** — frozen reference for legacy excited-state behavior.
- v3 still reuses a few working v2 utilities during development; those
  dependencies shrink as v3 matures.

v3's result object (`ResponseWorkflowOutput`) is a strict **superset** of v2's
`Results` — it adds first-class `timing`, `diagnostics`, and `exports`.

Switching an existing v2 deck to v3 is a one-line `engine v3` opt-in — see
[`MIGRATION_FROM_V2.md`](MIGRATION_FROM_V2.md) for what changes (and what to keep
on v2 for now).

## Architecture & design docs

- [`docs/16_architecture.md`](docs/16_architecture.md) — the master architecture:
  the L0→L5 layer cake and the R0–R5 build sequence (read after `00_status.md`).
- [`docs/00_status.md`](docs/00_status.md) — live status: active workstreams,
  hot-file conflict map, standing contracts. **Read this first** as a developer.
- [`CLAUDE.md`](CLAUDE.md) — the module's working contracts and build/run loop.
- `docs/` (numbered) — per-layer design detail: type system, persistence,
  scheduler/calc-manager, properties, excited-state and parallelism designs.

## Build / test

The canonical build/run/validate loop is the `cm.sh` harness
(`/gpfs/scratch/ahurtado/madness_es_bench/README.md`); MRA solves run on a
compute node. To build just the response binary from a MADNESS build directory:

```bash
cmake --build . --target molresponse_v3
```
