# Property Architecture: Response Properties vs Chemical Properties

Status: **design, agreed 2026-05-29.** Captures the two-tier property split
and the user-facing specification model. Builds on the response-property
planner (`ResponsePropertyPlanner.hpp`) and the unified persistence layer
(`13_unified_persistence_schema.md`). Supersedes the property framing
sketched in `04_property_design.md` (now a pointer to this doc).

---

## The core idea

molresponse's fundamental job is narrow and well-defined: **solve excited
states and response states, and contract them into the Cartesian tensor
components of response properties.** Everything else a chemist cares about is
a *transformation* of those tensors — often combined with non-electronic data
(normal modes, atomic masses, geometry) — into observable quantities with
physical units.

This motivates a hard split into two layers:

- **Tier A — Response Properties (the core driver).** Cartesian tensors that
  come *directly* from contracting response/excited states. Specifying which
  response properties you want, at which frequencies / directions / accuracy,
  IS the specification of a molresponse calculation — it determines exactly
  which states must be solved.

- **Tier B — Chemical Properties (post-processing + fancy printing).** Take
  Tier-A tensors (plus auxiliary data) and transform them into chemist-facing
  observables with units. This is reporting, not solving.

The split is not cosmetic. It draws a contract boundary that keeps the engine
simple and the reporting flexible, and it matches how mature QC codes
(Dalton, NWChem, Q-Chem) separate "compute polarizability derivatives" from
"generate the Raman spectrum."

---

## The invariant that defines the boundary

> **A Tier-B chemical property never triggers a new electronic solve.**

If producing a quantity requires solving another response or excited state,
that solve is a Tier-A response property and must be expressed as one (so the
planner can schedule it and the persistence layer can cache/restart it).
Tier B only consumes already-computed Tier-A tensors and auxiliary,
non-electronic ingredients.

Auxiliary ingredients Tier B may pull in (none of which are electronic-
response solves):
- normal modes / Hessian eigenvectors (themselves assembled from Tier-A
  nuclear-displacement response, but the *transformation* to modes is
  post-processing),
- atomic masses, molar mass,
- geometry, point-group / symmetry data,
- empirical broadening widths, temperature, laser frequency for prefactors.

---

## The canonical example: Raman

Raman is the clearest illustration of the boundary, and the reason the
planner was recently corrected.

| | Tier A (response property) | Tier B (chemical property) |
|---|---|---|
| **Vibrational Raman** | ∂α/∂Q — polarizability gradient w.r.t. nuclear Cartesian displacement. `GradientMode::Nuclear`. | Transform Cartesian ∂α/∂Q → normal-mode activities via the modes; depolarization ratios; broadened cm⁻¹ spectrum. |
| **Resonance Raman** | transition polarizability gradient via excited states. `GradientMode::Resonant`. | Resonance-Raman excitation profile; mode-resolved enhancement. |

Both Ramans are the **same** Tier-A property — the polarizability gradient
∂α/∂Q. They differ only in *how* the derivative is obtained (nuclear-
displacement response vs excited-state transition densities). That is exactly
what `ResponsePropertyKind::PolarizabilityGradient` + `GradientMode`
encodes. The spectroscopic observable (the thing a chemist plots) is Tier B
in both cases.

The clean rule the example yields: **Tier A produces Cartesian tensors;
Tier B transforms to spectroscopic observables, applies units, broadens, and
prints.**

---

## Tier A: Response Properties (what the planner plans)

Implemented in `ResponsePropertyPlanner.hpp`. `plan_one(req)` turns one
`ResponsePropertyRequest` into a deduplicated `ResponsePlan`; `merge_plans`
combines many requests (deduping shared states, unioning protocol ramps).

| `ResponsePropertyKind` | Tensor | States planned |
|---|---|---|
| `Polarizability` | α(−ω;ω), rank 2 | dipole FD at each ω, per axis |
| `Hyperpolarizability` | β(−ω_σ;ω₁,ω₂), rank 3 | dipole FD at the `BetaProcess` triplet freqs |
| `PolarizabilityGradient` (`Nuclear`) | ∂α/∂Q, rank 2 / mode | dipole FD at optical ω + nuclear-disp FD at 0 (all atoms) |
| `PolarizabilityGradient` (`Resonant`) | transition-α gradient | ES bundle + dipole FD at excitation energies |

Future Tier-A additions (not yet planned): second hyperpolarizability γ,
optical-rotation tensor G′ (needs magnetic-perturbation FD), magnetizability.

`ResponsePlan` carries four state lists:
- `fd` — concrete dipole FD points.
- `es` — excited-state bundles.
- `derived_fd` — FD whose frequency comes from an ES root (`es_root_id == "*"`
  expands to one per converged root post-ES).
- `nuclear_fd` — nuclear-displacement FD; `pert.atom < 0` is the all-atoms
  sentinel, expanded to one per atom against the concrete molecule.

The two sentinels (`"*"` for ES roots, `nuc_*` for atoms) keep the planner
pure and molecule-independent; the calc manager resolves them at execution
time.

---

## Tier B: Chemical Properties (post-processing)

Not yet implemented — a future `ChemicalProperty` module. It reads Tier-A
tensors from `response_metadata.json`'s `properties/` subtree (doc 13),
combines them with auxiliary data, and emits a separate human-facing artifact
(report tables, spectrum files). It must NOT write back into the
state-driving metadata.

Representative chemical properties and their Tier-A inputs:

| Chemical property | Tier-A input | + Auxiliary | Output units |
|---|---|---|---|
| Mean polarizability ᾱ, anisotropy Δα | α(ω) | — | a.u. / Å³ |
| Vibrational Raman spectrum | ∂α/∂Q (nuclear) | normal modes, masses | activities, depol. ratios, cm⁻¹ |
| Resonance Raman profile | ∂α/∂Q (resonant) | ES energies, modes | excitation profile |
| First hyperpolarizability figures (β_HRS, β_∥) | β | — | a.u. / esu |
| Specific optical rotation [α]_D | G′(ω) | molar mass | deg·cm³·g⁻¹·dm⁻¹ |
| Two-photon cross-section σ | 2-photon tensor | line-shape | GM |

Simplest first report when Tier B begins: ᾱ / Δα from α — no auxiliary data,
exercises the read-tensors-and-print path end to end.

---

## User-facing specification: three tiers of control

"Reasonable defaults, but full control" → a layered input model. A thin
input-binding layer (YAML / madqc input parsing) translates the higher tiers
down to the explicit `ResponsePropertyRequest` the planner consumes. The
binding lives OUTSIDE the planner — it is a UX concern, not a planning one.

### UX-1 — name only (matches input-deck QC convention)

```yaml
property: polarizability             # α at default freqs, all axes, default ramp
property: hyperpolarizability_shg    # β-SHG at default freqs
property: vibrational_raman          # ∂α/∂Q nuclear; expands to per-atom FD
property: resonant_raman             # ∂α/∂Q via ES; n_roots from response.excited
```

Defaults (globally configurable):
- axes = {x, y, z}
- frequencies = `[0.0]` for static-implied; a small dynamic set (e.g. the
  sodium-D line ≈ 0.0773 a.u. and/or a standard vis set) for dynamic
- protocol ramp = the same ramp the SCF/response run used
- β process = SHG

### UX-2 — structured options (keep UX-1 defaults, override pieces)

```yaml
property: hyperpolarizability
  process:   shg | or | eope | static | all_triplets
  freqs:     [0.057, 0.0656]
  axes:      [x, z]
  protocols: [1e-4, 1e-6, 1e-8]
```

### UX-3 — explicit ResponsePropertyRequest (programmatic escape hatch)

The C++ struct directly, for bespoke requests or driver code. UX-1 / UX-2 are
sugar that lower to a list of UX-3 requests, which `merge_plans` then dedupes.

This mirrors how Dalton/NWChem/Q-Chem expose properties (a high-level keyword
plus an override block), and the dedupe + protocol-keyed provenance from the
persistence layer is strictly better than the typical code that re-solves the
same FD point at a finer threshold without reusing the coarse one.

---

## End-to-end data flow

```
                 ┌──────────────── Tier A: Response Properties ────────────────┐
 user input  →   ResponsePropertyRequest(s)
   (UX-1/2/3)        │  plan_one / merge_plans
                     ▼
                 ResponsePlan { fd, es, derived_fd, nuclear_fd }
                     │  calc manager: reconcile vs response_metadata.json,
                     │  expand "*" / nuc_* sentinels, dispatch solves
                     ▼
                 FD/ES solves  →  save_fd_state / save_es_roots
                     │            (doc 13 persistence + restart)
                     ▼
                 properties/<name>/<protocol_key>/  (Cartesian tensors)
                 └──────────────────────────────────────────────────────────┘
                     │
                 ┌───▼──────────── Tier B: Chemical Properties ───────────────┐
                 ChemicalProperty reporter: read tensors + aux (modes, masses)
                     │  transform → units → broaden
                     ▼
                 report.json / printed tables / spectra   (NO new solves)
                 └──────────────────────────────────────────────────────────┘
```

The **calc manager** sits between the planner and persistence: it is the
executor that turns a `ResponsePlan` into actual solves, honoring the
restart-precedence and aggregate-metadata machinery already built. Its design
is the next increment.

---

## Why this split matters

1. **The engine stays small.** molresponse only ever reasons about response
   properties and the states they need. New chemical observables are added in
   Tier B without touching the solver or the planner.
2. **Caching/restart is automatic for everything physical.** Because every
   electronic solve is a Tier-A response property keyed by protocol, the
   persistence layer caches and restarts it — including states shared between
   several chemical properties.
3. **Reporting is cheap to iterate.** Units, broadening, spectral conventions,
   and output formatting change in Tier B alone, re-reading cached tensors —
   no recomputation.
4. **The boundary is testable.** Tier A is pure planning + numerical assembly;
   Tier B is pure transformation of stored tensors. Both unit-test without a
   live solve.
