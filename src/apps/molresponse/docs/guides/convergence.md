# Convergence: gates, floors, the plateau detector, and what a residual buys

How a molresponse solve decides it is done, what the reported residual means for
the property, and the knobs that control it. Companion to
[`polarizability.md`](polarizability.md) and [`excited_states.md`](excited_states.md);
the derivations are in the release report (`madness-workspace/reports/2026-09-09_release_report`,
section "Numerical analysis").

## Gates

`ConvergencePolicy` (`solvers/convergence_policy.hpp`) resolves, at each protocol
rung with wavelet threshold `thresh`,

```
dconv          = max(thresh, dconv_user)
bsh_target     = 5 * dconv        # FD: ||x_new - x_old|| per channel
density_target = 5 * dconv        # FD and ES: ||rho_new - rho_old||
omega_target   = 5 * dconv        # ES: |omega_new - omega_old| per root
```

FD converges when every channel's BSH residual and density change are below
target. ES converges when every active root's density change **and** eigenvalue
step are below target — the amplitude residual is printed but not gated, because
it jitters in near-degenerate directions long after ω (a second-order quantity)
has settled.

`dconv_user` comes from the deck's `response { dconv }`. When the deck does not
set it, `ParameterManager::set_derived_values` derives `100 * thresh` of the finest
rung (gate `500 * thresh`), and `maxiter` follows the dft block.

## The residual floor

Within a rung the BSH residual does not go to zero. Measured on the seeded
closeout legs (H2O, LiH, C2H4; static and ω = 0.1 dipole perturbations):

| rung        | BSH residual floor | density residual |
|-------------|--------------------|------------------|
| 1e-6 / k8   | 2–6e-5             | keeps falling    |
| 1e-8 / k10  | 4–8e-6             | to ~1e-7         |

A `dconv` whose gate lies below the floor is never met (the first 1e-8 closeout
attempt with `dconv = thresh`, gate 5e-8, ran every leg to the 60-iteration cap
with a flat residual). Rule of thumb: keep `5 * dconv` at least 3x above the floor
of the rung; the derived default satisfies this at both measured rungs.

The floor is 20–60x thresh at 1e-6/k8 and 400–800x at 1e-8/k10; the two-point
exponent is 0.42 (≈ "√thresh"), a fit through two points that also differ in k,
not a derived law. The experiment that isolates the mechanism (thresh sweep at
fixed k, k sweep at fixed thresh, tightened operator/product tolerances) is
described in the report.

## Plateau (stall) detector

`stall.window` (default 6) and `stall.ratio` (default 0.1), deck block `response`.
Each iteration appends the normalised gate distance

```
FD: g = max_c max(bsh_c / bsh_target, drho_c / density_target)
ES: g = max_s max(drho_s / density_target, |dw_s| / omega_target)   (active roots)
```

to a per-rung history (reset when thresh changes). When `g > 1` and
`g_now > (1 - ratio) * g_{now - window}` the solve is **stalled**: `step()` sets
`State::stalled`, `converged()` returns true so the loop exits, the log prints
`[STALL] iter N: ...` and `Stopped at iter N (stalled ...)`, and the metadata entry
carries `"stalled": true`. The verdict is unchanged with respect to running to the
cap: FD applies the same best-effort acceptance as maxiter (`--accept-at-maxiter`,
metadata `accepted: true`, log `ACCEPTED best-effort @ stall`); ES reports the
bundle as not converged. `stall.window 0` disables the detector. An iteration
without a measurable |dw| records +inf, which never counts as a stall, so a cold
ES start cannot trip it.

## What the residual means for a property

With `r_B` the reported BSH residual and `theta = v - W x` the BSH source (the
quantity applied to the Green's function each iteration):

- Plain estimate `alpha = c <v|x>`: the error is **linear**,
  `|delta alpha| <= |c| ||theta|| ||r_B||`, with `||theta||` a few atomic units for
  the closeout legs. At the 1e-6 floor this is ≲ 2e-4 au on H2O's alpha ≈ 9, an
  order of magnitude below the 0.04 % agreement with the DALTON reference.
- Every property assembled from first-order vectors without a stationary
  functional (β and Raman via the 2n+1 contraction, 2PA via the single residue)
  inherits the same linear law per leg. The 0.6–5 % MADNESS–DALTON differences of
  the 2PA tables are therefore not iteration error.
- The stationary functional `J[x] = 2<v|x> - <x|(A - ω)|x>` has a **quadratic**
  error and can be evaluated one iteration late for free:
  `J[x_{n+1}] = <v|x_{n+1}> + <theta_{n+1} - theta_n | x_{n+1}>` (raw BSH output).
  This is the planned property-based stopping criterion.

The FD solver prints the source overlap `<v|x>` of each channel per iteration
(Verbose banner, `print_final`, and the convergence CSV columns
`property, gate, stalled`) so the linear law can be checked against a run's log.
