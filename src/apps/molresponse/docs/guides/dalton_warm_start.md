# DALTON warm-start

Seed a MADNESS calculation from a converged DALTON calculation: project DALTON's
orbitals or response vectors onto the multiresolution grid and start the MADNESS
solve from there instead of from a cold guess.

## Why

MADNESS converges to its own basis-set-free solution regardless of where it
starts, so a seed **never changes the answer**. It buys two other things:

1. **Cost** — far fewer iterations, because a converged Gaussian-basis solution is
   already close to the MRA solution.
2. **Correct state selection** — for excited states this matters more than the
   speed. A cold guess can converge onto a *different set of roots* than intended
   (skipping a state and picking up a higher one); seeding from DALTON's
   eigenvectors steers the solver onto the intended manifold.

The second point is the reason to care even when runtime is not the constraint: if
you need *the same states* DALTON found — to compare, or to compute a residue
property like two-photon absorption at a specific state — the seed is what makes
those runs comparable.

## What can be seeded

| target | tool | DALTON input | status |
|---|---|---|---|
| Ground state (SCF orbitals) | `seed_moldft_from_dalton` | `molden.inp` | ✅ |
| Excited states, TDA | `seed_from_dalton` | `RSPVEC` (excitation eigenvectors) | ✅ |
| Excited states, Full/RPA | `seed_from_dalton --full` | `RSPVEC` (X and Y blocks) | ✅ |
| Frequency-dependent response | — | `RSPVEC` (linear-response N(ω)) | 🔬 future |

DALTON writes two files we consume, both from an ordinary run:

- **`molden.inp`** — MO coefficients, basis, and geometry, in portable text. The
  source for the ground-state seed and for the AO→MO projection used by all seeds.
- **`RSPVEC`** — converged response vectors: the excitation eigenvectors (with the
  de-excitation `Y` block for RPA) and the linear-response vectors at whatever
  frequencies were requested.

DALTON's own binary restart (`SIRIUS.RST`) is *not* needed — `molden.inp` carries
the orbitals portably.

## Ground state

```bash
seed_moldft_from_dalton --molden=<dalton>/molden.inp --n-occ=5 --out=mad \
                        --L=200 --thresh=1e-4
moldft --input=gs.in          # deck contains `restart true`, matching prefix/L
```

This writes a `mad.restartdata` archive in moldft's own format, so moldft resumes
from it with **no change to the SCF code**. The MADNESS SCF then re-converges to
its own MRA minimum — the DALTON orbitals are a starting guess, not a constraint.

Two couplings to respect:
- **`--L` must match the deck's box size.** moldft's loader hard-errors on a box
  mismatch (it tolerates and re-projects a `k`/threshold mismatch).
- `--thresh` sets the projection resolution; seeding at the rung you intend to
  start from avoids an immediate re-projection.

## Excited states

```bash
seed_from_dalton --rspvec=<dalton>/RSPVEC --molden=<dalton>/molden.inp \
                 --n-occ=5 --roots=0,1,2,3 --omegas=<w0>,<w1>,<w2>,<w3> \
                 --full --calc-dir=<calc> --thresh=1e-4
```

Writes one N-root excited-state bundle into the calc directory; the next response
run resumes from it. Notes:

- **Match the bundle type to the solve.** `--full` writes a Full/RPA (X,Y) bundle
  for an `--es-full` solve; without it the bundle is TDA. A TDA bundle handed to a
  Full solve is rejected and the run falls back to the cold guess, so this must
  line up.
- **`--omegas` are the DALTON excitation energies** for the roots you list — they
  seed the eigenvalues alongside the vectors.
- If the DALTON run was CIS/TDA (no `Y` block), `--full` promotes with `Y = 0`.

The tool reports the RPA metric `‖X‖²−‖Y‖²` per root as a sanity check; it should
be ~0.5 in the spatial-orbital normalization used here.

### A seeded excited-state solve tracks the seeded states

**A seeded ES solve refines the states you hand it — it does not search for the
N lowest.** If the seeding basis cannot describe a low-lying state (e.g. a diffuse
state seeded from a non-augmented basis), its eigenvector for that root is a
*different, higher* state, and the seeded MRA solve will faithfully converge that
state — fast, `converged=true`, and silently in place of the true N-th lowest
(observed on H2O/cc-pVDZ: seeded root 3 → 0.4626 au in 2 iterations, while the
true 4th-lowest is 0.4096 au). This is a feature when you target a *specific*
state (2PA residues) and a trap when you mean "give me the N lowest". Every
seeded solve therefore prints a `[SEED-GUARD]` block — per-root seed overlap,
seed ω, and ω shift, recorded under `excited_states/<key>/seed_guard` in
`response_metadata.json` — with a loud warning on basin escape (overlap < 0.5
with every seed root) and a note when the solve was pure tracking. When a
reference ladder from a better basis is available (validation campaigns),
compare the reported energies against it by hand — the guard records them in
the metadata for exactly that purpose.

## Choosing where to start the protocol ladder

MADNESS normally climbs a resolution ladder — a cheap coarse rung to reach the
right basin, then the production rung. **A seeded run does not need the coarse
rung**, and skipping it is the fastest route: the seed already provides what the
coarse rung exists to produce.

This also removes a failure mode rather than just saving time: the coarse rung is
where a cold guess can latch onto the wrong set of roots, so a seeded run that
enters directly at the production resolution both runs faster and selects states
more reliably. A cold run, by contrast, *needs* the coarse rung — forced to skip
it, it spends every iteration at the expensive resolution and is the slowest
configuration of all.

Practical recipe: seed at the resolution you intend to run, and start there.

## Verifying a seed did what you think

Seeds fail quietly — a rejected bundle or a mis-scaled projection produces a run
that looks normal but is really a cold start (or worse). Three cheap checks, worth
doing before trusting any new seeded workflow:

1. **Confirm it loaded.** The log states whether the state was resumed and whether
   the protocol key matched exactly. A silent fall-back to the cold guess is the
   most common failure.
2. **Check an invariant.** Electron count / kinetic energy for a ground state (a
   doubled or halved density shows up immediately); the RPA metric for an
   excited-state bundle.
3. **Compare the converged result to a cold run.** They must agree — the seed
   changes the path, not the destination. If they differ, the seed is wrong, not
   better.

## Scope

Closed-shell, matching the response engine's scope. The projection tools run on a
single rank (`NP=1`); the seeded solve itself has no such restriction.

See also: [excited states](excited_states.md),
[two-photon absorption](two_photon_absorption.md).
