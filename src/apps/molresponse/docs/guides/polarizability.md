# Polarizability α

The linear response of the electronic dipole to an applied electric field — how
readily the molecule's charge cloud distorts, and the foundation every other
response property builds on.

## What it computes

The dipole–dipole response tensor α<sub>ij</sub>(ω): the *i*-component of the
induced dipole per unit field applied along *j*, at frequency ω,

    α_ij(ω) = −2 ( ⟨x_i(ω) | μ_j⟩ + ⟨y_i(ω) | μ_j⟩ )   (closed shell),

the paired-metric contraction of the first-order dipole response (x, y) against
the dipole operator. Static (ω = 0) gives the ordinary polarizability; finite ω
gives the dynamic polarizability, whose poles are the excitation energies.

## Status

| configuration | status |
|---|---|
| static (ω = 0) | ✅ supported |
| dynamic α(ω) | ✅ supported |

MADNESS is basis-set-free: α is converged in resolution on the protocol ladder
rather than extrapolated in a Gaussian basis.

## Validation

The MADNESS polarizability has been benchmarked against reference Gaussian
calculations across a large molecule set and the correlation-consistent basis
families; see the author's polarizability/hyperpolarizability benchmark
papers [1,2]. Results are not reproduced here — this guide documents the feature
and how to run it.

## Run it

```text
dft
    xc        hf
    protocol  [1e-4, 1e-6]        # resolution ladder (k climbs 6 → 8)
end
response
    dipole             true       # linear dipole response → polarizability
    dipole.directions  xyz
    dipole.frequencies [0.0, 0.04]  # a.u.; [0.0] = static only
end
```

`madqc --wf=response <deck>`; the tensor lands in `response_metadata.json` under
`properties/alpha`, keyed by (ω, directions). Full recipe:
[`madqc/RESPONSE_PROPERTIES.md`](../../../madqc/RESPONSE_PROPERTIES.md).

## Under the hood

The pure linear-response step: for each (direction, ω) the FD solver converges
the first-order response (x, y) on the ladder, then `assemble_alpha` contracts it
against the dipole operator (closed-shell factor −2). No quadratic source, no
excited-state input — which is why it is the reference path the higher properties
(β, Raman, 2PA) build on. Shared response formalism: [`formalism.md`](formalism.md).
Parallel/subworld execution: [`parallelism.md`](parallelism.md).

## References

1. Hurtado *et al.*, correlation-consistent basis-set benchmark of MADNESS
   polarizabilities. *(full citation to add)*
2. Companion hyperpolarizability benchmark. *(full citation to add)*
