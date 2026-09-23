# Hyperpolarizability β

The first nonlinear response of the dipole — the leading term behind
second-harmonic generation and the electro-optic effect, and the quadratic-response
reference the two-photon property builds on.

## What it computes

The first hyperpolarizability tensor β<sub>ABC</sub>(−ω<sub>A</sub>;
ω<sub>B</sub>, ω<sub>C</sub>): the quadratic response of the dipole to two applied
fields, with ω<sub>A</sub> = −(ω<sub>B</sub> + ω<sub>C</sub>). The common
frequency triplets are static (0; 0, 0), second-harmonic generation
(−2ω; ω, ω), electro-optic (−ω; ω, 0), and optical rectification (0; ω, −ω).

## Status

| configuration | status |
|---|---|
| static (0; 0, 0) | ✅ supported |
| second-harmonic generation (−2ω; ω, ω) | ✅ supported |

MADNESS is basis-set-free: β is converged in resolution on the protocol ladder
rather than extrapolated in a Gaussian basis.

## Validation

The MADNESS hyperpolarizability has been benchmarked against reference Gaussian
calculations across a large molecule set and the correlation-consistent basis
families; see the author's polarizability/hyperpolarizability benchmark
papers [1,2]. Results are not reproduced here — this guide documents the feature
and how to run it.

## Run it

```text
dft
    xc        hf
    protocol  [1e-4, 1e-6]
end
response
    quadratic          true                 # request hyperpolarizability
    dipole.directions  xyz
    dipole.frequencies [0.0, 0.04]           # SHG uses ω and 2ω internally
end
```

`madqc --wf=response <deck>`; the tensor lands in `response_metadata.json` under
`properties/beta`, keyed by (A, B, C, freqs). Full recipe:
[`madqc/RESPONSE_PROPERTIES.md`](../../../madqc/RESPONSE_PROPERTIES.md).

## Under the hood

β is the quadratic (2n+1) contraction: from the two first-order dipole responses
the engine assembles the **quadratic source**, forms the second-order response
density, and contracts it against the third dipole — no explicit second-order
solve. It reuses the linear (FD) responses that [α](polarizability.md) produces.
Shared (A, B, C) machinery: [`formalism.md`](formalism.md). Parallel/subworld
execution: [`parallelism.md`](parallelism.md).

### One source, two builders

The theory has one second-order right-hand side (Hurtado eq. 19 = Parker 2018
eq. 28 = Sałek 2002 eq. 68); β, Raman and two-photon absorption all contract it.
Two builders exist in the code and are asserted equal (up to Q-projection and
truncation order) by `test_vbc_spec_equivalence` and gate 5 of
`test_tpa_pq_vs_vbc`:

- `tpa::quadratic_source` (`kernels/tpa_source_spec.hpp`) — the (P,Q) spec; one
  call sums both photon orderings. **Default for β and Raman** since 2026-09-09.
- `vbc::compute_vbc` (`kernels/vbc.hpp`) — the same source written term by term
  in the equation's leg orientation (`gamma_legs` / `gamma_dagger_legs` /
  `ketbra` from `kernels/source_spec.hpp`).

The test driver keeps `--beta-vbc-source` (select the other builder) and
`--beta-compare-sources` (build both, print both β values) as diagnostics.

### Validation at finite frequency

Before 2026-09-09 `vbc.hpp` (inherited from molresponse_v2) built every response
density transposed: |φ><x| where the equation has |x><φ|. Static β is blind to
this (x = y, Kleinman-exact), which is why it was never caught; SHG at ω = 0.1 au
for H₂O against DALTON HF/d-aug-cc-pVQZ at the same geometry was off by +12.9 %
(β_zxx) and −4.0 % (β_xzx). After the fix (recontraction of the same saved
first-order legs, protocol 1e-6/k8):

| component (ω = 0.1) | MADNESS | DALTON | dev |
|---|---|---|---|
| β(z;x,x) | 1.3595 | 1.3542 | +0.39 % |
| β(x;x,z) | 4.3700 | 4.3553 | +0.34 % |
| β(z;y,y) | 12.674 | 12.663 | +0.08 % |
| β(y;y,z) | 12.492 | 12.482 | +0.08 % |
| β(z;z,z) | 11.442 | 11.437 | +0.04 % |

These DALTON values are the `beta_shg` references at ω = 0.1 in
`madness-workspace/refs/madness_results.json`. Finite-frequency MRA SHG tables
produced with the v2 builder should be re-examined before being relied on.

## References

1. Hurtado *et al.*, correlation-consistent basis-set benchmark of MADNESS
   hyperpolarizabilities. *(full citation to add)*
2. Companion polarizability benchmark. *(full citation to add)*
- Companion guides: [polarizability](polarizability.md), [Raman](raman.md),
  [two-photon absorption](two_photon_absorption.md).
