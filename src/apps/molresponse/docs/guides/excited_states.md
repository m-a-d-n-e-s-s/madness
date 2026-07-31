# Excited States

Vertical excitation energies and transition properties from linear-response
theory — the spectrum a molecule absorbs, and the input the two-photon residue
is built on.

## What it computes

The excitation energies ω<sub>f</sub> and their response eigenvectors, as the
poles of the linear (dipole) response — i.e. the linear-response / random-phase
eigenvalue problem. From the converged roots the engine reports transition
dipoles and oscillator strengths per state. Two closed-shell flavors:

- **TDA** — the Tamm–Dancoff (excitation-only) approximation.
- **Full / RPA** — the paired (X, Y) response; the flavor the two-photon residue
  and dynamic polarizability poles use.

## Status

| configuration | status |
|---|---|
| TDA, closed shell | ✅ supported |
| Full / RPA, closed shell | ✅ supported |
| open shell | 🚧 out of scope (future) |

Roots are solved on the resolution ladder (k climbs with the protocol) with a
warmup + KAIN acceleration; MADNESS is basis-set-free, so the energies converge
in resolution toward the complete-basis limit rather than in a Gaussian basis.

## Preliminary results

*Preliminary — the authoritative campaign lives in a separate results record.*

**Basis-free excitation energies.** Because MADNESS refines toward the
complete-basis limit, it is the reference the finite-basis DALTON numbers
converge onto. Water, lowest excitations (eV), MADNESS at three resolutions vs
DALTON's Gaussian basis ladder:

| water root | 1 | 2 | 3 | 4 (diffuse) |
|---|---|---|---|---|
| DALTON aug-cc-pVDZ | 8.53 | 10.21 | 10.86 | 12.06 |
| DALTON aug-cc-pVQZ | 8.54 | 10.21 | 10.81 | 11.54 |
| DALTON aug-cc-pV5Z | 8.54 | 10.21 | 10.80 | 11.35 |
| DALTON d-aug-cc-pVQZ | 8.54 | 10.20 | 10.77 | 11.12 |
| **MADNESS (k6 → k10)** | **8.54** | **10.20** | **10.77** | **11.11** |

MADNESS is converged (its k6/k8/k10 agree to ~0.003 eV) and matches DALTON's best
basis; the diffuse root 4 shows that single augmentation is short of the limit by
up to ~0.9 eV even at 5Z, while double augmentation reaches it.

### DALTON warm-start (seeding)

The excited-state solve can be warm-started from a converged DALTON calculation
(its CIS/RPA eigenvectors projected onto the MRA grid). Preliminary, water:

| start | iterations to converge |
|---|---|
| standard `solid_harmonics` guess (TDA) | 18 |
| DALTON seed (TDA) | **6** |
| standard guess, Full/RPA | 28 |
| DALTON seed, Full/RPA | **8** |

Two benefits, both observed: fewer iterations (≈3×), and **correct root
selection** — the seed steers the solver onto the states DALTON found, recovering
a root the bare guess skipped. A near-converged seed also lets the solve **start
directly at the production resolution** and skip the coarse warmup rung entirely,
which is the fastest route overall. This is the same warm-start the two-photon
property relies on.

## Run it

```text
dft
    xc        hf
    protocol  [1e-4, 1e-6]
end
response
    excited_states  true
    es.n_states     4          # number of roots
    es.full         true       # Full/RPA (omit for TDA)
end
```

`madqc --wf=response <deck>`; roots and transition properties land in
`response_metadata.json` under the excited-state records. Full recipe:
[`madqc/RESPONSE_PROPERTIES.md`](../../../madqc/RESPONSE_PROPERTIES.md). Seeding
from DALTON uses the `seed_from_dalton` tool.

## Under the hood

The ES solver runs an oversampled TDA warmup to seed the roots, then iterates the
response eigenvalue problem with KAIN on the resolution ladder, locking converged
roots. The initial guess is `solid_harmonics` by default, or a projected DALTON
eigenvector bundle when seeding. Parallel/subworld execution across roots:
[`parallelism.md`](parallelism.md).

## References

- MADNESS linear-response / excited-state literature. *(full citations to add)*
- Companion: [two-photon absorption](two_photon_absorption.md),
  [polarizability](polarizability.md).
