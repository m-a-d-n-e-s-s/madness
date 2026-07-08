# Excited-state benchmark suite (legacy vs v1)

Times each codebase on a ladder of systems, all running the same 4-root
TDA + RPA at protocol [1e-6], starting from a freshly-converged
moldft archive. The point is to measure the **rotate-vs-recompute gap**:
v1 carries the gamma matrix through subspace rotations, legacy recomputes
it. The gap should grow with `n_occ × n_states`, so we ramp `n_occ`.

## Layout

```
benchmarks/
├── README.md              (this file)
├── h2o/                   (n_occ=5,  4 nodes × 8 tasks/node × 10 threads, 4 hr,  hbm-short-96core)
│   ├── benchmark.slurm
│   ├── response_tda.in
│   └── response_rpa.in
├── c2h4/                  (n_occ=8,  4 nodes × 8 tasks/node × 10 threads, 8 hr,  hbm-short-96core)
│   ├── benchmark.slurm
│   ├── response_tda.in
│   └── response_rpa.in
└── c6h6/                  (n_occ=21, 8 nodes × 4 tasks/node × 20 threads, 24 hr, needs a long partition)
    ├── benchmark.slurm
    ├── response_tda.in
    └── response_rpa.in
```

The ground-state geometries and `moldft.in` files come from the existing
fixtures in `../systems/{h2o_hf, c2h4_hf, c6h6_hf}/`.

## What each `benchmark.slurm` does

1. Stages `moldft.in` into a per-job work dir under `$WORK_BASE`.
2. Runs **moldft** once to converge the ground state at protocol [1e-4, 1e-6].
3. For each `(code, mode) ∈ {legacy, v1} × {tda, rpa}`:
   - Symlinks the moldft archive into `run_<code>_<mode>/`.
   - Copies `response_<mode>.in` (named `input` for legacy, `response.in`
     for v1).
   - Runs the binary, captures wall-clock + iteration count + final omegas
     into `timing.csv`.
4. Prints a side-by-side summary at the end of the SLURM stdout.

## Before submitting

Defaults are set for Seawulf Xeon Max HBM nodes (account `ahurtado`,
partition `hbm-short-96core`). Likely edits per host:

- **C6H6 needs a longer partition** — `hbm-short-96core` caps walltime
  before 24 hr. Switch its `--partition` line to your long-running queue
  (e.g. `hbm-long-96core` or `extended-96core`) before submitting.
- The launcher defaults to `mpirun --map-by numa numactl --preferred-many=8-15`
  (HBM tier preferred). To fall back to plain `srun`:
  `LAUNCHER=srun sbatch benchmark.slurm`.
- `source ~/load_xeonmax.sh` must exist and set up your build's runtime
  environment. If you use modules instead, replace that line.

Path overrides (all optional, take env vars):
- `$REPO=$HOME/Projects/madness`
- `$BUILD=$REPO/build`
- `$WORK_BASE=/gpfs/scratch/$USER/madness_es_bench`

## Build prerequisites

On Seawulf, before submitting:

```bash
cd $REPO/build
make moldft -j8
make molresponse_legacy -j8
make molresponse        -j8           # this is v1
```

## Submit

All three at once:

```bash
cd src/apps/molresponse/tests/fixtures/benchmarks
./submit_all.sh                # all three: h2o c2h4 c6h6
./submit_all.sh h2o            # just one
./submit_all.sh h2o c2h4       # subset, in order
```

Or by hand:

```bash
cd src/apps/molresponse/tests/fixtures/benchmarks/h2o   && sbatch benchmark.slurm
cd ../c2h4                                                 && sbatch benchmark.slurm
cd ../c6h6                                                 && sbatch benchmark.slurm
```

Each job is independent. C6H6 will take the longest; consider running
H2O first as a smoke test before committing C6H6's 24-hour reservation.
The script writes `jobs.txt` next to itself with one `<sys> <jobid>` line
per submission.

## After they finish

Each job leaves a `timing.csv` in its work dir with this schema:

```
system,code,mode,wall_seconds,iters,omega1,omega2,omega3,omega4,exit
```

A simple comparison across systems:

```bash
WORK_BASE=/gpfs/scratch/$USER/madness_es_bench
( head -1 $WORK_BASE/h2o_*/timing.csv | head -1
  for sys in h2o c2h4 c6h6; do
    tail -n+2 $WORK_BASE/${sys}_*/timing.csv 2>/dev/null
  done ) | column -t -s,
```

You should see `wall_seconds` for `code=v1` consistently below `code=legacy`
within the same `(system, mode)` pair. The ratio is the rotate-vs-recompute
speedup — the headline number we want for the design discussion.

## What we're looking for

The rotation savings live inside `compute_gamma_full` (the J/K/XC build).
With protocol [1e-6] and 4 roots, expected costs scale roughly as:

| n_occ | gamma cost  | rotate-savings expectation |
|-------|-------------|-----------------------------|
| 5  (H2O)  | small        | ~5-15% (small but measurable) |
| 8  (C2H4) | medium       | ~15-30%                       |
| 21 (C6H6) | dominates    | ~40-60% (this is the target)  |

If the gap on C6H6 is < 30% something is off and we re-think before
porting the rotate-everything pattern to v3. If it's > 50%, port it.
