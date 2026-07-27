# Migrating from the v2 response engine to v3

**The v2 response engine (`MolresponseLib.hpp`) has been removed.** v3 is now the
sole response engine behind `madqc --wf=response`, and it is the default — there
is nothing to opt into. If you have an old deck that still selects v2, the
migration is a single edit: **delete the `engine v2` line.** The rest of the deck
(`dft`, `molecule`, and the `response` knobs) is unchanged, and v3's output is a
superset of what v2 wrote.

> During v3's development v2 was the parity reference — v3 was validated to
> reproduce v2's `calc_info.json` values (e.g. H₂O α_zz) at every increment.
> That parity work is done and v2 is gone; this guide is now just the deck-level
> "what changed" note that `madqc` points old decks at.

---

## The only change: drop the engine line

`madqc --wf=response` no longer takes an engine choice — it always runs v3. A
deck that still contains `engine v2` fails loudly:

```text
ERROR: response engine 'v2' (MolresponseLib) was removed.
       Delete the `engine v2` line (v3 is the default) — the
       input deck is otherwise unchanged; see
       molresponse_v3/MIGRATION_FROM_V2.md.
```

So just remove it (an explicit `engine v3` is still accepted but redundant):

```text
response
    dipole              true
    dipole.directions   z
    dipole.frequencies  [0.0]
end
```

Everything else — the `dft` block (xc, protocol, dconv, …), the `molecule`
block, and the response perturbation/property knobs (`dipole.*`, `quadratic`,
`beta.shg`/`beta.or`, `excited.*`, `requested_properties`) — is parsed by the
same `ResponseParameters`, so a deck that ran under v2 runs under v3 unchanged;
only the solver behind it differs.

---

## Property status

v3 does not yet cover every property v2 did; the honest per-property state
(supported / experimental / in development) lives in the capability matrix in
[`README.md`](README.md). Since v2 is gone there is no fallback engine — for a
property still marked experimental or in-development, treat v3's output
accordingly and confirm the printed values.

---

## What you get from v3 (output differences)

v3 writes the **same** response keys v2 did into `<prefix>.calc_info.json`
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

## See also

- [`README.md`](README.md) — v3 status + capability matrix.
- [`../madqc/RESPONSE_PROPERTIES.md`](../madqc/RESPONSE_PROPERTIES.md) —
  per-property knobs, examples, and the response result schema.
- [`../madqc/README.md`](../madqc/README.md) — the general madqc driver.
