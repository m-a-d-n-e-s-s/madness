# Migrating from the v2 response engine to v3

Short version: **your input deck does not change.** v3 is selected by a single
line in the `response` block (`engine v3`); the `dft`, `molecule`, and the rest
of the `response` knobs are identical, and v3's output is a **superset** of v2's.
Today v3 supports a subset of properties (α, static β, single-component Raman) —
for everything else, stay on v2. This guide says exactly what changes and how to
check the two engines agree.

> v2 (`MolresponseLib.hpp`) is the **production default and the parity
> reference**. v3 is opt-in and validated against v2 at every increment. Nothing
> below removes or changes the v2 path.

---

## The only change: select the engine

`madqc --wf=response` picks the engine from the `engine` knob in the `response`
block. The default is `v2`. To use v3:

```text
response
    engine              v3          # <- the one change; default is v2
    dipole              true
    dipole.directions   z
    dipole.frequencies  [0.0]
end
```

Everything else — the `dft` block (xc, protocol, dconv, …), the `molecule`
block, and the response perturbation/property knobs (`dipole.*`, `quadratic`,
`beta.shg`/`beta.or`, `excited.*`, `requested_properties`) — is **the same**.
They are parsed by the same `ResponseParameters`, so a deck that runs under v2
runs under v3 unchanged; only the solver behind it differs.

---

## Should you migrate? (per property)

| Property | Recommendation |
|----------|----------------|
| **Polarizability α** (static + dynamic) | ✅ v3 ready — validated vs v2 |
| **First hyperpolarizability β** (static + dynamic SHG/OR) | ✅ v3 ready — validate vs published benchmark β; for dynamic β confirm all VBC contributions climb |
| **Single-component Raman** | 🟡 v3 runs it, but unvalidated — no reference yet; confirm the printed values |
| **Full-tensor Raman** | ⚠️ stay on v2 — v3 full tensor is in development |
| **Excited states, 2PA, resonant Raman** | ⚠️ stay on v2 — v3 ES is experimental, 2PA/resonant in development |
| **Large systems (state-parallel)** | ⚠️ stay on v2 — v3 scaling layer is in development |

See [`README.md`](README.md) for the full capability matrix and
[`../madqc/RESPONSE_PROPERTIES.md`](../madqc/RESPONSE_PROPERTIES.md) for the
per-property recipes.

---

## What you get from v3 (output differences)

v3 writes the **same** response keys v2 does into `<prefix>.calc_info.json`
(`properties.response_properties[]`, `raman_spectra`, `vibrational_analysis`) —
so existing parsers keep working unchanged. v3 **adds** fields under the response
task's `metadata` (additive — they don't displace anything v2 wrote):

- `engine` = `"molresponse_v3"` (so you can tell which engine produced a record),
- `timing` — structured stage/point wall+cpu timing,
- `diagnostics` — convergence trajectories, divergence/stall flags, memory HWM,
  scheduler trace,
- (forthcoming) `exports` — a manifest of MRA/grid/viz artifacts.

If you parse `calc_info.json`, **no change is required**; the new fields are
opt-in extras you can ignore or consume.

---

## Numerical agreement & a parity check

v3 is validated to reproduce v2 within the protocol's target precision (e.g.
H₂O α_zz matches to the converged digits). If you want to confirm it for your
system, run the **same deck through both engines** and diff the values:

```bash
# engine v2
madqc --wf=response --response="engine=v2" mymol.in
# engine v3
madqc --wf=response --response="engine=v3" mymol.in
```

Compare the `properties.response_properties[].value` entries (and `raman_spectra`)
between the two `calc_info.json` files. The `cm.sh` harness automates exactly
this: `cm_madqc_parity <input.in>` runs both engines on one deck and prints the
α values side by side (`/gpfs/scratch/ahurtado/madness_es_bench/README.md`).

A mismatch beyond the protocol tolerance is a v3 bug — report it against the
parity reference, don't silently prefer one engine.

---

## See also

- [`README.md`](README.md) — v3 status + capability matrix.
- [`../madqc/RESPONSE_PROPERTIES.md`](../madqc/RESPONSE_PROPERTIES.md) —
  per-property knobs, examples, and the response result schema.
- [`../madqc/README.md`](../madqc/README.md) — the general madqc driver.
