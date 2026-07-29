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
the engine assembles the **VBC quadratic source**, forms the second-order response
density, and contracts it against the third dipole — no explicit second-order
solve. It reuses the linear (FD) responses that [α](polarizability.md) produces.
Shared (A, B, C) machinery: [`formalism.md`](formalism.md). Parallel/subworld
execution: [`parallelism.md`](parallelism.md).

## References

1. Hurtado *et al.*, correlation-consistent basis-set benchmark of MADNESS
   hyperpolarizabilities. *(full citation to add)*
2. Companion polarizability benchmark. *(full citation to add)*
- Companion guides: [polarizability](polarizability.md), [Raman](raman.md),
  [two-photon absorption](two_photon_absorption.md).
