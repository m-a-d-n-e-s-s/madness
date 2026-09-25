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
| Frequency-dependent response | in-solver, `dalton.dir` | `RSPVEC` (linear-response N(ω), exact frequency; `seed.freq_tol` for nearest) | ✅ |
| Excited states, in-solver | `dalton.dir` | `RSPVEC` (EXCITLAB records) | ✅ |
| Ground state, in-solver | `dalton.dir` (madqc pre-run hook) | `molden.inp` | ✅ |

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

## In-solver seeding from the deck (`dalton.dir`)

Since 2026-09-09 every stage can be seeded from ONE DALTON directory named in the
deck (`io.dalton.dir`), without running the projection tools by hand:

```
dft
  protocol [1e-6]        # one fine rung: the seed is basis-set quality; no econv/dconv (the
end                      # protocol is the knob — SCF::solve floors dconv at the threshold)
io
  backend    hdf5                          # run-wide archive backend (GS and response)
  dalton.dir /path/to/dalton/seed          # run-wide seed: loose RSPVEC + molden.inp, or a unique *.tar.gz
end
response
  protocol        [1e-6]
  seed.start_rung fine                      # moot with a single rung; kept for ladders
  seed.freq_tol   0.01                      # nearest DALTON N(ω) for legs not in the RSPVEC
  seed.es_y       zero                      # zero | dalton: write the DALTON Y block into the ES seed?
  seed.es_warmup  false                     # false: the Full solve starts from the seed directly
end
```

`dalton.dir` and the backend are deck-level (`io`) because they serve every
stage — the SCF, the frequency-dependent legs and the excited states; the
response-block `dalton.dir` / `hdf5` remain as aliases and lose when both are
given.

- **Ground state.** `madqc` installs a pre-run hook on the SCF application: when
  `dalton.dir` is set and no `<prefix>.restartdata` exists, the molden orbitals are
  projected (Löwdin-orthonormalized) and written as the SCF restart, so moldft
  starts from the DALTON orbitals. A preserved copy `<prefix>.gs_seed.restartdata`
  (plus `.h5` when `hdf5` is on) and a `<prefix>.gs_seed.json` record what was
  seeded; the converged GS is mirrored to `<archive>.h5` when `hdf5` is on.
- **Frequency-dependent response.** A dipole FD leg whose frequency matches an
  `XDIPLEN` record to 1e-9 au starts from that vector. Derived legs (the two-photon
  legs at ω_f/2) never match; `seed.freq_tol > 0` takes the closest record within
  the tolerance as the initial guess (`seed_kind = dalton_nearest` in the metadata).
- **Excited states.** `EXCITLAB` records (sorted by energy) become the initial
  `es__<key>` bundle at the active protocol, gauge-rotated onto the MRA orbitals,
  Q-projected, marked unconverged (iteration 0); an existing bundle wins. The Full
  solver does not start from that bundle directly: an iteration-0 bundle is routed
  seed → KAIN-free TDA warm-up (`excited.tda_warmup_iters` steps from the seed's X
  block) → promotion to Full with y = 0, the same path the cold guess takes. A
  seed handed straight to the Full solver has a residual of order 0.1 against the
  MRA operator and ran to the −ε_core ghost within three iterations on H₂O, LiH
  and C₂H₄ (2026-09-10), with or without the DALTON de-excitation block; the
  TDA-converged start the cold path provides is what the Full solver needs.
  `seed.es_y` (`zero` default, `dalton`) selects whether the DALTON Y block is
  written into the bundle at all; it is discarded by the warm-up either way.

  The standalone tool `seed_from_dalton --roots=... --omegas=...` indexes
  `EXCITLAB` records in **file** order, which is not energy order for
  `**PROPERTIES .EXCITA` outputs; pass the omegas in that order or prefer the
  in-solver `dalton.dir` path, which sorts by energy.

DALTON side, learned the hard way:

- `*LINEAR` and `*QUADRA` cannot share one `**RESPONSE` block ("RSPSYM: BOTH LR AND
  QR CALC SPECIFIED"). Use `**PROPERTIES .EXCITA` for the excitation vectors (it
  writes `EXCITLAB` records) together with `**RESPONSE *LINEAR .FREQUENCIES` for the
  N(ω) grid, or split into two runs and merge.
- The scratch tarball is `<dal>_<mol>.tar.gz`, or `<dal>.tar.gz` when the two
  stems coincide.

The gecko `SeededCampaign` (madness-workspace/workflow) writes these decks and the
SLURM chain (DALTON seed → madqc) for the closeout molecules.

## Scope

Closed-shell, matching the response engine's scope. The projection tools run on a
single rank (`NP=1`); the seeded solve itself has no such restriction.

See also: [excited states](excited_states.md),
[two-photon absorption](two_photon_absorption.md).
