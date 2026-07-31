# Two-Photon Absorption (2PA)

The strength with which a molecule jumps to an excited state by absorbing two
photons at once, each carrying half the transition energy — a genuinely nonlinear
spectroscopy, and the newest response property in the stack.

## What it computes

For each excited state *f*, a symmetric two-photon transition-moment tensor
between the ground state and *f*: the amplitude for absorbing two photons of given
polarizations, each at half the transition energy. The experimentally relevant
outputs are its rotational invariants and the linear/circular two-photon
strengths and cross-sections (Monson–McClain averages).

Conceptually, 2PA is obtained as a **residue of the quadratic (hyperpolarizability)
response** taken at an excitation, rather than as a separate property built from
scratch — so it reuses the excited-state and frequency-response machinery. The
precise working equations are still being finalized and will be documented in
[`formalism.md`](formalism.md); this guide intentionally does not pin them down.

## Status

| configuration | status |
|---|---|
| closed-shell, singlet final states | ✅ supported |
| two-photon tensor + Monson–McClain observables | ✅ |
| open-shell final states | 🚧 out of scope (future) |

## Preliminary results

*Preliminary — the authoritative validation campaign lives in a separate results
record, not this code release.* Cross-checked against DALTON's quadratic-response
two-photon module. The comparison is on the **rotationally-invariant** observables,
which are the physical quantities — within a degenerate manifold the raw Cartesian
tensor components rotate freely between codes while the invariants are conserved.

| system | basis | preliminary agreement with DALTON |
|---|---|---|
| H₂O | d-aug-cc-pVQZ | dominant tensor elements within **0.3–1.2%** |
| LiH (ionic; includes a degenerate manifold) | aug-cc-pVQZ | Monson–McClain invariants (D<sub>f</sub>, D<sub>g</sub>, δ, R) within **≤0.4%** |
| C₂H₄ (π) | aug-cc-pVQZ | bright roots within **3–4%** |

The C₂H₄ residual is a basis effect, not a method one: single augmentation is
short of the complete-basis limit for diffuse states (see
[`excited_states.md`](excited_states.md)), which MADNESS reaches natively.

### DALTON warm-start

2PA consumes converged excited states, so it benefits directly from seeding the
excited-state solve from a DALTON calculation (see [`excited_states.md`](excited_states.md)):
the DALTON eigenvectors both cut the excited-state iteration count and steer the
solver onto the intended root manifold — which matters here because the 2PA tensor
is taken at a specific excited state.

## Run it

2PA is requested as an excited-state residue property in the `response` block; see
[`madqc/RESPONSE_PROPERTIES.md`](../../../madqc/RESPONSE_PROPERTIES.md). It needs
the excited state *f* plus the dipole responses at half the transition energy.
Results land in `response_metadata.json` under `properties/tpa` — the two-photon
tensor and the Monson–McClain observables per state.

## References

- Quadratic-response / two-photon formalism; Monson & McClain for the rotational
  averages. *(full citations to add once the formalism section is finalized)*
- Companion: [excited states](excited_states.md), [hyperpolarizability](hyperpolarizability.md).
