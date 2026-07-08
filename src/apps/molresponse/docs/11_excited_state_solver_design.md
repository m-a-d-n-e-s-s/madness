# Excited-State Solver Design (molresponse_v3)

Status: **draft for review** — no code yet.

Goal: an ES solver that handles **closed-shell (RHF) and open-shell (UHF)**
ground states, in both **TDA** and **Full TDDFT/TDHF** flavors, using the
same building blocks as the FD solver and starting from the same
`GroundState` data.

This doc focuses on what's new for ES on top of the FD machinery already
in place. It is meant to be read alongside `02_type_system_design.md`
(types, building blocks, FD↔ES sharing) and `07_increment_plan.md`
(increment ordering, test plan).

---

## What's already in place (for free)

The FD work has produced most of the infrastructure ES needs:

| Piece | File | Used by ES as-is? |
|---|---|---|
| `RealResponseState` (x_alpha, y_alpha, x_beta, y_beta) | `ResponseFunctions.hpp` | yes — one per root |
| `compute_response_density` (per-state, both spins) | `ResponseKernel.hpp` | yes |
| `compute_gamma` (per-spin, with TDA/Full/Static branches) | `ResponseKernel.hpp` | yes |
| `make_bsh_operators` / `_y` (per-state, omega-dependent) | `ResponseKernel.hpp` | yes — rebuilt each iter as ω changes |
| `fd_iteration_step` (BSH apply per spin) | `ResponseKernel.hpp` | yes — reused for theta→x_new |
| `compute_lambda` (V0·x − ε·x + γ) | `ResponseKernel.hpp` | yes |
| `compute_subspace_matrix` / `compute_overlap_matrix` | `ResponseKernel.hpp` | yes |
| `diagonalize_subspace` (sygvp) | `ResponseKernel.hpp` | yes |
| `GroundState` (orbitals, Q-projectors, Fock no-diag, V_local) | `GroundState.hpp` | yes — same instance feeds FD and ES |
| `QProjector` per spin | `GroundState.hpp` | yes |

What this means: there is **no new operator-layer code** required.
The work is at the bundle/iteration-loop level.

---

## What ES adds

1. **Bundle storage** — N excited states, each a `RealResponseState`,
   plus per-root metadata (current ω estimate, residual norm, root_id).
2. **Eigenvalue solve in subspace** — diagonalize Λ in the current
   bundle each iteration to get current ω estimates and the rotation.
3. **Per-root BSH** — BSH operators are ω-dependent and ω is a moving
   target per root, so BSH operators are rebuilt per root per iteration.
4. **Initial guess generation** — for ES, no perturbation drives the
   RHS; we need to seed the bundle. Two viable approaches:
   - **(a) Single-orbital primitives**: x_i = phi_i × g(r) for each
     pair (i, basis function) — produces O(n_occ × n_guess) candidates,
     orthonormalize and keep the lowest N.
   - **(b) HOMO-LUMO style**: x_i = phi_HOMO × phi_a where phi_a is
     a virtual approximation. Legacy uses this with localized Gaussian
     guess primitives. Simpler in the absence of virtuals: guess from
     localized Gaussians multiplied by occupied orbitals.
5. **Bundle re-orthogonalization** — after each BSH update, the bundle
   can drift to linear dependence; Gram-Schmidt across roots before
   the next Λ build.
6. **Convergence per-root + collective** — track per-root residual
   and energy change; converged when *all* roots meet the threshold.
7. **Root identity tracking** — across iterations roots can swap due
   to rotation; a stable `root_id` is assigned at first iteration and
   tracked by maximum-overlap matching with the previous iteration's
   bundle.
8. **Step restriction** — ES is more prone to oscillation than FD,
   especially when nearby roots cross. A simple step-restriction on
   the BSH update norm per root mirrors the legacy `maxrotn` knob.

---

## Proposed API

Mirrors `fd_solve`. Lives in `ESSolver.hpp`.

```cpp
namespace molresponse_v3 {

struct ESSolveResult {
    std::vector<RealResponseState> roots;   // converged X bundle
    Tensor<double> omegas;                  // excitation energies
    Tensor<double> residuals;               // per-root final residual
    std::vector<int> root_ids;              // stable root labels
    int iterations;
    bool converged;
};

ESSolveResult es_solve(
    World& world,
    ResponseType type,        // TDA or Full (Static is FD-only)
    int num_states,           // bundle size
    GroundState& gs,
    int maxiter = 25,
    double dconv = 1e-4,
    double maxrotn = 0.5,
    int maxsub = 10,
    PrintLevel print_level = PrintLevel::Normal,
    const std::vector<RealResponseState>* initial_guess = nullptr);

} // namespace
```

Closed-shell vs open-shell is **not** an API parameter — it follows
from `gs.is_spin_restricted()`, exactly as in `fd_solve`.

---

## Bundle storage

Use `std::vector<RealResponseState>` directly. No new wrapper type
in the first cut.

Reasoning: `RealResponseState` already carries x_alpha/y_alpha/x_beta/y_beta
with the right closed/open-shell branching. A wrapper buys us nothing
at the kernel level — root metadata (ω, residual, root_id) lives in
parallel arrays in `ESSolveResult` and in the iteration loop itself.

If/when persistence and naming arrive (Increment 4 / 9), a `Bundle`
type with archive I/O may be worth introducing. Not now.

---

## Iteration loop (pseudocode)

```
es_solve(type, N, gs):
    X = initial_guess(N, gs, type)            # vector<ResponseState>, size N
    omega = initial_omegas(N)                  # rough spacings, refined iter 1
    root_id = [0, 1, ..., N-1]                 # initial labels

    for iter in 0..maxiter:

        # 1. Build Λ for each root using its CURRENT ω
        #    Λ_s = compute_lambda(X_s, gs, type)   per spin like fd_iteration
        Lambda_X = [compute_lambda_full(world, type, X_s, gs) for s in 0..N]

        # 2. Build A and S
        #    For Full: A and S use the symplectic metric (x·x' − y·y')
        A = compute_subspace_matrix(X.x_flat, Lambda_X)   # plus y contribution for Full
        S = compute_overlap_matrix(type, X.x, X.y)

        # 3. Diagonalize → new ω, rotation U
        omega_new, U = diagonalize_subspace(A, S)

        # 4. Track roots: maximum-overlap match between U and previous U
        root_id = match_roots(U, prev_U, root_id)

        # 5. Rotate the bundle:  X_s = sum_t U_ts * X_t
        X = rotate(X, U)
        Lambda_X = rotate(Lambda_X, U)

        # 6. Per-root BSH update:
        #    theta_s = -2 * (Lambda_X_s − ω_s · X_s)
        #    BSH operators rebuilt with current ω_s
        X_new = []
        for s in 0..N:
            bsh_alpha_x = make_bsh_operators(eps_alpha, omega_new[s], lo)
            bsh_alpha_y = (Full) ? make_bsh_operators_y(eps_alpha, omega_new[s], lo) : []
            (similarly for beta)
            theta_s.x_alpha = -2 * (Lambda_X_s.x_alpha − omega_new[s] * X_s.x_alpha)
            (similarly for y_alpha, x_beta, y_beta)
            X_new_s = apply_bsh_step(theta_s, bsh_*)    # Q-project inside

        # 7. Step-restrict per root
        for s: if ||X_new_s − X_s|| > maxrotn: damp toward X_s

        # 8. Re-orthogonalize the bundle (Gram-Schmidt with Full metric)
        X_new = orthonormalize_bundle(X_new, type)

        # 9. Convergence
        residual_s = ||X_new_s − X_s||,  domega_s = |omega_new[s] − omega[s]|
        if max(residual) < res_tol AND max(domega) < omega_tol: converged

        X = X_new;  omega = omega_new;  prev_U = U

    return ESSolveResult{X, omega, residual, root_id, iter, converged}
```

This is structurally close to legacy `iterate_excited.cc` but uses
v3 building blocks. KAIN per-root is deferred until after the basic
loop converges (legacy uses a `kain_x_space[b]` per state — easy to add
once the no-KAIN path is correct).

---

## Open-shell handling

For UHF, each `RealResponseState` already carries x_alpha and x_beta
(and y_alpha/y_beta if Full). The kernel pieces handle this:

- `compute_response_density` sums alpha + beta channels with
  spin_factor=1 for unrestricted.
- `compute_gamma` is called per-spin (same pattern as `fd_iteration`).
- Λ for the bundle is computed per-spin and the spins are stitched into
  one flat vector for the subspace matrix dot product:
  `inner(X_s.flat(), Lambda_X_t.flat())`.
- BSH operators are built per-spin from `gs.energies_alpha()` /
  `gs.energies_beta()`.

The only UHF-specific question is **the initial guess**: alpha and
beta channels need correlated guesses to avoid spin contamination
artefacts in the early iterations. Simplest first cut: independent
guesses per spin, let convergence sort it out. Refine later if it
causes problems on Li / OH / NO.

---

## Initial guess strategy (concrete proposal)

**TDA RHF, first cut**: for N requested states and n_occ occupied,
seed each root s ∈ [0, N-1] as

```
X_s.x_alpha[i] = δ_{i, s mod n_occ} × g_s(r) × phi_i(r)
```

where `g_s(r)` is a localized Gaussian primitive at the molecule's
center of charge, with width chosen by atomic number. Then orthonormalize.

This is crude but enough to get the solver running. Refinements:

- Use multiple Gaussian widths to span the low-energy spectrum.
- Use the dipole-perturbation states (already producible via
  `dipole_perturbation`) as a guess basis for the first 3 roots.
- For UHF, generate alpha and beta guesses independently from the
  respective orbital sets.

Concretely: borrow the legacy `LocalizedGaussianGuess` functor (in
v2 `ExcitedStateBundleSolver.cpp`) and adapt it to v3 types.

---

## Convergence criteria

A root is converged when:

```
||X_new_s − X_s|| < residual_target
AND
|ω_new[s] − ω[s]|  < omega_target
```

with

```
residual_target = max(thresh × 10, dconv) × max(5, n_atom)   # same as FD
omega_target    = dconv                                       # tighter
```

Bundle is converged when **all** roots are individually converged.
Stalled roots (residual oscillating without improvement for K iterations)
get a step-restriction tightening, mirroring FD's stall handling
in v2.

---

## Increment ordering

I'd build this in four small sub-increments rather than one big "ES
solver" milestone. Each is independently testable.

| Sub-inc | What | Test |
|---|---|---|
| 7a | TDA RHF only, simple guess, no KAIN | H2 lowest TDA root → match legacy ω_1 |
| 7b | TDA UHF | Li atom lowest few roots vs Dalton/legacy |
| 7c | Full RHF (TDDFT) | H2 lowest TDDFT root vs legacy |
| 7d | Full UHF | OH or Li dynamic — vs legacy |

Within each sub-increment, the loop is:
1. Wire up the iteration on a known case
2. Validate against legacy or analytic value
3. Iterate on numerical issues (step restriction, guess quality)
4. Add KAIN per-root once basic convergence is solid

---

## Files added / changed

| File | Change |
|---|---|
| `ESSolver.hpp` | new — top-level `es_solve` and helpers |
| `ESSolverGuess.hpp` | new — `LocalizedGaussianGuess` adapted from v2 |
| `ResponseKernel.hpp` | small additions: bundle Λ wrapper, root-matcher |
| `test_es_solver.cpp` | new — H2 TDA / Li UHF TDA validation |
| `CMakeLists.txt` | add `test_es_solver` target |
| `docs/11_excited_state_solver_design.md` | this file |

No changes to `ResponseFunctions.hpp`, `GroundState.hpp`,
`Perturbations.hpp`, or `FDSolver.hpp` are anticipated.

---

## Test plan (first pass)

Mirroring the FD test setup:

| System | Type | What we check |
|---|---|---|
| H2 (RHF, k=6, dconv=1e-4) | TDA | ω_1 vs legacy (~0.46 hartree); root identity stable |
| H2 (RHF, k=6) | Full TDDFT | ω_1 vs legacy; check that Full < TDA (correlation) |
| Li (UHF, 1 atom, k=6) | TDA | first 3 roots vs Dalton cc-pV5Z; alpha/beta channel population |
| Li (UHF) | Full TDDFT | spin contamination check; ω_1 vs legacy |

Each test reuses the `mad.restartdata` checkpoints already produced
for the FD validation runs. No new fixture inputs needed.

---

## Open questions for review

1. **Bundle wrapper now or later?** I propose `vector<RealResponseState>`
   for the first cut and revisit with persistence. Acceptable?
2. **KAIN per-root in 7a, or deferred?** Legacy uses one KAIN per root;
   I'd defer it to keep 7a small. Risk: poorer convergence on Li-like
   close-spacing spectra. Acceptable to defer?
3. **Initial guess source.** Adapt v2's `LocalizedGaussianGuess`
   (faster) vs port legacy directly (cleaner reference). Adapt v2 unless
   that locks us into v2-isms.
4. **Symplectic metric for Full.** `compute_overlap_matrix` already
   subtracts ⟨y|y⟩ for Full. Does the BSH update path also need to
   account for this when computing residuals? I think yes — needs
   verification against legacy `iterate_excited.cc` lines around the
   `bsh_residualsX` / `bsh_residualsY` split.
5. **Davidson-style subspace expansion?** Legacy expands the subspace
   beyond N roots and re-collapses. Not needed for first cut (we work
   with exactly N roots) but worth noting as a future improvement.

Sign-off on these five lets me start on 7a.
