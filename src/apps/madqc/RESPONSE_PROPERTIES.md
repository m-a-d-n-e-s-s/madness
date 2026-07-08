# Computing response properties with madqc

How to compute molecular **response properties** — polarizability, hyper­polar­iz­ability,
Raman, and excited states — through `madqc --wf=response`. This is the per-property
companion to the general driver guide ([`README.md`](README.md)); for the engine itself
and its current capability matrix see
[`../molresponse/README.md`](../molresponse/README.md).

> **Engines.** `madqc --wf=response` has two response engines. `engine v2`
> (default) is the production path; `engine v3` is the next-generation unified
> solver (opt-in). The recipes below target **v3**; set `engine v2` to use the
> production engine. Validate new results against v2 (same input → matching
> `calc_info.json`).

---

## How property selection works

The `response` block selects properties two ways (you can mix them):

1. **High-level flags** — `dipole`, `quadratic`, etc. madqc derives the
   `requested_properties` list from them.
2. **Explicit list** — `requested_properties [polarizability, hyperpolarizability, raman]`.

A response run always needs the **perturbation directions and frequencies**
(`dipole.directions`, `dipole.frequencies`) and a **protocol ladder** (the
truncation-threshold schedule, set in the `dft` block as `protocol`). The run
solves the ground-state SCF first, then the response.

```text
dft
    xc        hf
    protocol  [1e-4, 1e-6]      # threshold ladder (k climbs 6 → 8)
    dconv     1e-4
end
response
    engine              v3
    dipole              true
    dipole.directions   xyz
    dipole.frequencies  [0.0, 0.04]
end
molecule
    units angstrom
    O   0.0  0.0     0.1173
    H   0.0  0.7572 -0.4692
    H   0.0 -0.7572 -0.4692
end
```

---

## Polarizability α  ✅ supported (static + dynamic)

The linear dipole–dipole response α(ω). Static (ω=0) and frequency-dependent.

| Knob | Meaning |
|------|---------|
| `dipole true` | request linear dipole response (→ `polarizability`) |
| `dipole.directions xyz` | Cartesian components to solve (`x`, `y`, `z`, or any subset) |
| `dipole.frequencies [0.0, 0.04]` | driving frequencies ω (a.u.); `[0.0]` = static only |

```text
response
    engine              v3
    dipole              true
    dipole.directions   z          # α_zz only (fast); use xyz for the full tensor
    dipole.frequencies  [0.0]
end
```

**Result** (`<prefix>.calc_info.json` → response task):
`properties.response_properties[]` entries with `property="polarizability"`,
`component=["z","z"]`, `freqB=<ω>`, `value=<α>`. Example (H₂O/HF): α_zz(0)=8.53.

> Requesting a single direction (`z`) solves one perturbation channel — much
> faster, and the natural choice for a quick check.

---

## First hyperpolarizability β  ✅ supported (static + dynamic)

The quadratic response β, static and dynamic (SHG/OR) — a benchmark-quality
capability; validate against published reference β values. **Operational
caveat:** for dynamic β, make sure all VBC contributions reach the top protocol
(check the convergence diagnostics) — otherwise the assembled β is incomplete.

| Knob | Meaning |
|------|---------|
| `quadratic true` | compute β from the defined perturbations (adds `hyperpolarizability`) |
| `beta.shg true` | second-harmonic generation β(−2ω; ω, ω) (default) |
| `beta.or true` | optical rectification β(0; ω, −ω) |
| `dipole.frequencies` | the ω used to form the β process |

```text
response
    engine              v3
    dipole              true
    dipole.directions   z
    dipole.frequencies  [0.0]      # static β
    quadratic           true
end
```

**Result:** `response_properties[]` with `property="hyperpolarizability"`,
`component=["z","z","z"]`, `freqB`/`freqC` = the driving frequencies, `value=<β>`.
Example (H₂O/HF): β_zzz(0,0)≈7.76.

---

## Vibrational Raman  🟡 single component (validation pending)

Raman intensity from β(dipole; dipole, nuclear) — the response to a nuclear
displacement. The method is benchmark-quality in principle, but in v3 it is
**not yet validated**: there is no reference value on file, and the correctness
of the printed Raman output has not been confirmed — treat results as
provisional until checked. **Currently single-component** (one `(atom, axis)`);
the full per-atom/axis tensor is in development (it depends on the state-parallel
layer).

| Knob | Meaning |
|------|---------|
| `requested_properties [raman]` | request the Raman (polarizability-gradient) property |
| `dipole.directions` / `dipole.frequencies` | the dipole perturbation |
| `nuclear.directions` / `nuclear.frequencies` | the nuclear-displacement perturbation |

```text
response
    engine              v3
    requested_properties [raman]
    dipole.directions   z
    dipole.frequencies  [0.0]
    nuclear.directions  z
    nuclear.frequencies [0.0]
end
```

The nuclear-displacement FD source is cusped and **resolution-limited**: it
floors at coarse protocols and only sharpens by climbing (e.g.
`protocol [1e-4, 1e-6]`). **Result:** `properties.raman_spectra`.

---

## Excited states (TDA / RPA)  🟡 experimental

Excitation energies and (with analysis) transition properties. Eigenpairs
compute; **convergence and climbing across protocols are being stabilized**, so
treat ES results as experimental and check residuals.

| Knob | Meaning |
|------|---------|
| `excited.enable true` | run the excited-state bundle |
| `excited.num_states N` | number of roots to target |
| `excited.tda true` | Tamm-Dancoff approximation; `false` = full RPA (X,Y) |
| `excited.maxiter` / `excited.maxsub` | iteration cap / subspace size |

```text
response
    engine          v3
    excited.enable  true
    excited.num_states 3
    excited.tda     true
end
```

**Scope:** closed-shell only. Open-shell ES / open-shell full-RPA are out of
scope. ES / diffuse-state references use **d-aug-cc-pVQZ**.

---

## Convergence & protocol knobs

| Where | Knob | Effect |
|-------|------|--------|
| `dft` | `protocol [1e-4, 1e-6]` | truncation-threshold ladder; each rung raises k and reprojects |
| `dft` | `dconv` | SCF density convergence |
| `response` | `maxiter` | response iteration cap |
| `response` | `dconv` | response convergence target |
| `response` | `print_level` | 0–3, response verbosity |

A coarse single rung (`protocol [1e-4]`) is fastest for a smoke test; climb
(`[1e-4, 1e-6]`) for converged values (k 6 → 8).

---

## The response result (`calc_info.json`)

The response task in `<prefix>.calc_info.json`:

| Key | Meaning |
|-----|---------|
| `type` | `"response"` |
| `properties.response_properties[]` | computed tensor components |
| `…[].property` | `"polarizability"` or `"hyperpolarizability"` |
| `…[].component` | axis labels, e.g. `["z","z"]` (α) or `["z","z","z"]` (β) |
| `…[].freqB` (+ `freqC` for β) | driving frequencies (a.u.) |
| `…[].value` | the tensor component |
| `properties.raman_spectra` | per-mode Raman intensities (when computed) |
| `properties.vibrational_analysis` | normal-mode frequencies (when computed) |
| `metadata` | per-state / per-protocol convergence + (v3) `timing`, `diagnostics` |

The human-readable `<prefix>.out` renders the same numbers (see the example in
[`README.md`](README.md)). For the full top-level schema (SCF/CIS/correlated
tasks, tensor-json encoding) see the **Output** section of [`README.md`](README.md).

---

## See also

- [`README.md`](README.md) — the general madqc driver: workflows, CLI, input
  format, full `calc_info.json` schema, visualization.
- [`../molresponse/README.md`](../molresponse/README.md) — the v3 response
  engine: capability matrix, status, v2 relationship.
- [`../molresponse/MIGRATION_FROM_V2.md`](../molresponse/MIGRATION_FROM_V2.md)
  — switching a v2 deck to `engine v3` (and what to keep on v2 for now).
- `cm.sh` (`/gpfs/scratch/ahurtado/madness_es_bench/README.md`) — the build/run/
  validate harness; `cm_mq <mol>` drives exactly this response path.
