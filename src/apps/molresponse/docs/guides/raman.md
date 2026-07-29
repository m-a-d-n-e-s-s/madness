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

The polarizability-gradient contraction reuses the **VBC quadratic source** and
the same first-order (FD) responses as [β](hyperpolarizability.md); only one
perturbation operator changes (dipole → nuclear displacement), and the operator's
sign/normalization convention for the mixed case is handled in assembly. Shared
(A, B, C) machinery: [`formalism.md`](formalism.md).

## References

- Quadratic-response Raman / polarizability-gradient literature. *(full citation to add)*
- Companion guides: [hyperpolarizability](hyperpolarizability.md),
  [polarizability](polarizability.md).
