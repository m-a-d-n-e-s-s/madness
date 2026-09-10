# Raman

Raman scattering intensities from polarizability derivatives with respect to
nuclear motion — the same quadratic-response machinery as β, with a
nuclear-displacement perturbation in place of one optical field.

## What it computes

The polarizability gradient ∂α/∂Q along the nuclear displacement Q, from which
Raman activities follow. Formally it is the quadratic response
P<sub>ABC</sub> with A, B dipole operators and **C a nuclear-displacement
operator** — the polarizability differentiated with respect to that displacement
— so it slots directly into the (A, B, C) contraction the engine already runs for
β.

## Status

| configuration | status |
|---|---|
| static Raman (polarizability gradients) | ✅ supported |
| resonant Raman | 🚧 future (scoping done) |

## Validation

Raman reuses the β quadratic-response contraction with a different operator on one
leg; it inherits the shared kernel's validation. A dedicated Raman benchmark (and
resonant Raman) is future work — no results are reproduced in this release guide.

## Run it

Request Raman in the `response` block; because it needs the nuclear-displacement
perturbation as well as the dipole responses, see the recipe and operator setup in
[`madqc/RESPONSE_PROPERTIES.md`](../../../madqc/RESPONSE_PROPERTIES.md). Results
land in `response_metadata.json` under `properties/raman`.

## Under the hood

The polarizability-gradient contraction reuses the **quadratic source** and
the same first-order (FD) responses as [β](hyperpolarizability.md); only one
perturbation operator changes (dipole → nuclear displacement), and the operator's
sign/normalization convention for the mixed case is handled in assembly. Shared
(A, B, C) machinery: [`formalism.md`](formalism.md).

The quadratic source is the one second-order right-hand side of the theory
(Hurtado eq. 19 = Parker 2018 eq. 28 = Sałek 2002 eq. 68), built by
`tpa::quadratic_source` (`kernels/tpa_source_spec.hpp`); `vbc::compute_vbc`
(`kernels/vbc.hpp`) is the same source written term by term and is asserted
equal to it by `test_vbc_spec_equivalence`. Until 2026-09-09 the Raman and β
paths used a builder whose response densities were transposed, which is invisible
in the static limit (x = y) but wrong at finite frequency; see
[hyperpolarizability](hyperpolarizability.md#validation-at-finite-frequency).
Validated at finite frequency on 2026-09-10 against a DALTON HF/d-aug-cc-pVQZ
finite-difference reference (α_zz(ω) at O displaced ±0.005/±0.010 bohr,
Richardson-extrapolated) at the HF/aug-cc-pVQZ optimized H₂O geometry:
∂α_zz/∂z_O = 4.28802 vs 4.28783 (ω = 0, +0.004 %) and 4.45217 vs 4.45349
(ω = 0.1, −0.03 %), protocol 1e-6/k8. The Raman correctness tolerance in
`madness-workspace/refs/madness_results.json` is 2 %.

The response exchange operator must use moldft's `multiworld_row` algorithm
(`kernels/common_ops.hpp`, override `MADRESPONSE_HFEXALG`): the tiled
`multiworld` algorithm gives α 0.9 % low and this derivative 18.5 % low on the
asymmetric bra/ket of response densities.

## References

- Quadratic-response Raman / polarizability-gradient literature. *(full citation to add)*
- Companion guides: [hyperpolarizability](hyperpolarizability.md),
  [polarizability](polarizability.md).
