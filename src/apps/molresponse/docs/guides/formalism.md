# Response formalism

The shared machinery every property in this engine specializes. The individual
feature guides ([polarizability](polarizability.md),
[hyperpolarizability](hyperpolarizability.md), [Raman](raman.md),
[two-photon absorption](two_photon_absorption.md)) point back here rather than
re-deriving it. This is a working summary; the full derivation and benchmarks are
in the accompanying paper(s).

## The perturbation expansion

Response properties are derivatives of the energy (equivalently, of the density
matrix) with respect to one or more perturbing fields. Expanding the density
matrix γ in the field strengths λ<sub>C</sub>,

    γ = γ⁰ + Σ_C λ_C γ^C + Σ_BC λ_B λ_C γ^BC + …

gives the ground state γ⁰, the **first-order** responses γ^C (linear), the
**second-order** responses γ^BC (quadratic), and so on. Each order satisfies a
response equation obtained by collecting terms at that order in the
Fock/commutator relation.

## First-order (linear) response

For a perturbation v^C at frequency ω<sub>C</sub>, the first-order response is
carried by a pair of orbital-space functions (x<sub>i</sub><sup>C</sup>,
y<sub>i</sub><sup>C</sup>) — the excitation and de-excitation parts — one pair per
occupied orbital *i*. They solve the coupled response equations, cast in MADNESS
as a bound-state Helmholtz (BSH / Green's-function) fixed-point iteration on the
resolution ladder:

    x_p^C = −2 Ĝ(k_p) * [ V⁰ x_p^C − Σ ε_ip x_i^C + g′[γ^C] φ_p + V_p^C ]

and the conjugate equation for y. Here g′[γ^C] is the first derivative of the
electron-interaction operator evaluated with the response density (the coupling
that makes the equations self-consistent). This linear solve is the **FD solver**;
it is the single most reused component in the engine.

## Property contraction — the (A, B, C) form

A property is a trace of a perturbation operator against a response density. In
the generalized notation used throughout the engine,

    P_ABC(−ω_A; ω_B, ω_C) = Tr[ v^A γ^BC ]

with the frequency sum rule ω_A = −(ω_B + ω_C). Choosing the operators and orders
selects the property. The quadratic response density γ^BC contains a lower-order
part γ_L (products of first-order responses, plus an occupied-space relaxation
term ζ) and a genuinely second-order part γ_Q built from the **VBC source** — the
quadratic right-hand side assembled from the two first-order responses and their
operators.

## Specializations

**Polarizability α** — the linear special case: A and C are dipole operators, no
quadratic source. α<sub>ij</sub>(ω) = −2(⟨x<sub>i</sub>(ω)|μ<sub>j</sub>⟩ +
⟨y<sub>i</sub>(ω)|μ<sub>j</sub>⟩). See [polarizability](polarizability.md).

**Hyperpolarizability β** — the quadratic case with A, B, C all dipole operators;
β<sub>ABC</sub>(−ω<sub>A</sub>; ω<sub>B</sub>, ω<sub>C</sub>) contracts the dipole
against the VBC-sourced γ^BC. See [hyperpolarizability](hyperpolarizability.md).

**Raman** — the same quadratic machinery with a **nuclear-displacement** operator
substituted for one leg (the polarizability gradient), rather than a third dipole.
See [Raman](raman.md).

## Two-photon absorption — being finalized

<!-- PLACEHOLDER. The two-photon transition moment is obtained as the single
     residue of the quadratic response at an excitation energy. The precise
     working composition — in particular the two-electron E[3] residue term that
     distinguishes it from an ordinary β evaluation — is still being finalized and
     will be written out here once settled. Until then the two-photon guide
     describes the property and its preliminary results without committing to the
     equations. -->

Two-photon absorption is obtained as the **single residue of the quadratic
response** (β) taken at an excited state, rather than as a separate property. It
therefore reuses the excited-state and frequency-response machinery above. The
precise working equations — including the two-electron residue correction that
makes it distinct from an ordinary β evaluation — are **being finalized** and will
be added to this section; see [two-photon absorption](two_photon_absorption.md)
for the current descriptive treatment and preliminary results.

## References

- The MADNESS response / correlation-consistent basis benchmark papers (α, β).
- Quadratic-response and two-photon literature (to be cited when the 2PA section
  is finalized).
