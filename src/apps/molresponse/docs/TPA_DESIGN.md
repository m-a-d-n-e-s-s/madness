# Two-Photon Absorption (TPA) — design

## Context
TPA is "β with different ABC": it is the **single residue of the quadratic
response** `⟨⟨μ_α; μ_β, μ_γ⟩⟩`. Where β (hyperpolarizability) contracts three FD
response states, TPA replaces **one index with an excited-state transition
vector** `X_f` (the residue of the response at the pole ω → ω_f) and evaluates at
the **two-photon resonance** ω₁ = ω₂ = ω_f/2. So TPA reuses the same contraction
core (`beta_abc` / VBC machinery) plus the converged ES bundle.

## Physics / target quantity
The two-photon transition-moment tensor between ground and excited state f:

    S_αβ^{0→f}  =  Σ_p [ ⟨0|μ_α|p⟩⟨p|μ_β|f⟩ + ⟨0|μ_β|p⟩⟨p|μ_α|f⟩ ] / (ω_p − ω_f/2)

In response theory this is the **single residue** of `⟨⟨μ_α; μ_β, μ⟩⟩_{ω₁,ω₂}` as
ω₂ → ω_f, with ω₁ = −ω_f/2. The rotationally-invariant TPA observables (what to
validate) are:

    δ_TPA = Σ_αβ ( F·S_αα S_ββ* + G·S_αβ S_αβ* + H·S_αβ S_βα* )

with (F,G,H) = (2,2,2) for parallel linear polarization (the common δ^∥). The raw
S_αβ tensor per state is the primary computed object.

## Reuse map (what already exists)
- **ES bundle**: converged `X_f`, ω_f persist via `solvers/es_save_load.hpp`
  (`es__<pkey>/roots.json` + per-root archives); load with `load_es_roots`.
- **Two-photon resonance FD**: the calc manager already builds **derived dipole
  FD at ω_f/2** (`DerivedFD`, `es_freq_factor = 0.5`, used today for resonance
  Raman / off-pole seeding). ω_f/2 IS the two-photon resonance frequency — exactly
  what S_αβ needs. So the linear responses are already produced by the resonant
  pipeline.
- **Contraction core**: `apply_gamma_raw` (shared two-electron action) + the
  `beta_abc` term structure. TPA is the same contraction with the C-FD response
  swapped for `X_f` and the A/B dipole responses taken at ω_f/2.

## Proposed implementation
1. **Property kind**: add `ResponsePropertyKind::TwoPhotonAbsorption` (+ request
   fields: `n_roots`, `axes`). Plan = ES bundle (TDA or Full) + derived dipole FD
   at ω_f/2 for each axis — i.e. the existing `PolarizabilityGradient/Resonant`
   plan, re-tagged. (Likely a thin variant of that branch in
   `ResponsePropertyPlanner.hpp`.)
2. **Kernel** `kernels/tpa.hpp`:
   ```cpp
   // Per excited state f: the 3x3 two-photon transition-moment tensor S_αβ.
   madness::Tensor<double>
   tpa_moment(World&, const ResponseGroundState& g0,
              const ResponseStateX<Shell>& Xf,          // ES eigenvector (root f)
              const std::array<ResponseStateXY<Shell>,3>& mu_response_half_omega, // FD μ_a at ω_f/2
              const std::array<real_function_3d,3>& mu_op);  // raw dipole ops
   ```
   The body is the β residue: contract `Xf` against the dipole-dressed linear
   responses (and the ground-orbital transition-dipole terms), symmetrized in α↔β.
   It reuses `apply_gamma_raw` for any two-electron pieces and `common_ops`
   inner products — NOT a fresh electronic solve.
3. **Executor + assembly**: after the ES bundle + ω_f/2 FD converge, a Tier-B
   `assemble_tpa(ctx, plan)` loads `X_f` + the FD responses, calls `tpa_moment`
   per root, records `properties/tpa` (S tensor + δ^∥) in `response_metadata.json`.
   Pure post-processing (off the critical solve path), mirroring `assemble_beta`.

## Open questions (settle before coding the kernel)
- **Exact residue form / normalization**: the precise set of terms (which
  relaxation/`zeta`-type terms survive in the single residue vs the full β) and
  the ES-vector normalization convention (⟨X_f|X_f⟩ metric for TDA vs RPA).
  Derive against Olsen–Jørgensen quadratic-response residues; pin with a small
  case before trusting numbers.
- **TDA vs Full (RPA)**: TDA gives a simpler residue (no Y); the RPA residue uses
  the symplectic metric. Start TDA (matches the validated ES path).
- **Damping**: undamped (resonance-divergent) vs damped (complex) TPA. Start
  undamped, off-resonance states.

## Validation
- No TPA reference in the repo yet. Options: Dalton quadratic-response single
  residue (`.QUADRA`/two-photon module) for a small molecule (water, LiH) at HF,
  same geometry as the v3 fixtures; or a literature δ_TPA. The α-isolation method
  applies: if α + ES + β all PASS but TPA fails, the bug is in `tpa_moment`.
- Effort: kernel + theory deep-dive is the bulk; the planning/execution/assembly
  scaffolding is a thin variant of the resonance-Raman / β paths.

## Relation to Raman
Raman (β with a nuclear C) is being implemented first — it validates the
"different-ABC" generalization with a non-dipole **operator** but still all-FD
states. TPA is the next step up: a different **state** (ES eigenvector) in place of
an FD response. Doing Raman first de-risks the operator/contraction plumbing.
