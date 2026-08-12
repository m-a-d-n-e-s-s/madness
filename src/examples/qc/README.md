# Quantum-chemistry example decks (`qctest`)

Every directory here is a **complete, CI-verified `madqc` calculation**: an input
deck, the one-liner that runs it, the results it must produce, and a reference
run. They serve three purposes at once, which is the point — a deck that is
documentation but not a test rots, and a test that is not a deck teaches nobody.

1. **Copy-paste starting points.** Each `*.in` is a working deck for one
   workflow. Nothing here is templated or generated.
2. **Integration tests.** Registered with ctest under the `qctest` label, so a
   deck cannot silently stop working.
3. **Ground truth for agents.** When asked "what does a CIS input look like",
   the answer is a file in this directory that provably ran, not prose that may
   describe knobs which no longer exist.

## Running a case

```bash
cd nemo_he_hf
./run.sh                  # uses `madqc` from PATH
MADQC=/path/to/build/src/apps/madqc/madqc ./run.sh    # or point it at a build tree
```

Via ctest, from a configured build directory:

```bash
ctest -L qctest                        # every case
ctest -L qctest -LE verylong           # everything but the expensive ones
ctest -R "madness/test/qc/nemo_he_hf"  # one case
cmake --build . --target check-qctest-madness
```

**Every run starts from an empty working directory.** `madqc` writes its per-task
scratch to `<case>/task_0/<calc>/` and restarts from `mad.restartdata.*` if it finds
it, so a second run in the same directory converges onto the first run's orbitals
rather than solving the deck. `run_qctest.py` therefore erases the working directory
before staging the case. It only does so where it can tell the directory is a scratch
one -- empty, carrying the `.qctest_workdir` marker, or holding this case's own output
-- and otherwise refuses rather than deleting; `--clean` overrides.

Cases tagged `short` or `medium` are also picked up by `check-short-madness`
(`ctest -L "short|medium"`), so they gate CI along with the unit tests.

## Cases

Wall times are measured on node26 (96 cores, `MAD_NUM_THREADS=20`).

| Case | `--wf=` | System | Demonstrates | Time | Tier |
|------|---------|--------|--------------|------|------|
| `scf_he_hf` | `scf` | He | the minimal deck — start here | 5 s | short |
| `scf_he_hf_mpi` | `scf` | He | same deck on 2 MPI ranks; thread budget and `--bind-to none` | 7 s | short |
| `nemo_he_hf` | `nemo` | He | regularized (nuclear-cusp-free) orbitals | 12 s | short |
| `nemo_he_varyk` | `nemo` | He | `k` unpinned, so it varies 6 → 8 across the ladder; must match `nemo_he_hf` | 13 s | medium |
| `oep_be_oaep` | `oep` | Be | optimized effective potential, OAEP model; virial diagnostics | 28 s | medium |
| `cis_he_singlets` | `cis` | He | CIS excited states; the `tdhf` group | 9 s | medium |
| `scf_h2o_hf` | `scf` | H₂O | `protocol` ladder 1e-4 → 1e-6 | 38 s | long |
| `scf_lih_gopt` | `scf` | LiH | geometry optimization inside an SCF task (`dft gopt`) | 38 s | long |
| `scf_lih_optimize` | `scf` + `--optimize` | LiH | geometry optimization as its own task, moldft reference | 37 s | long |
| `nemo_lih_optimize` | `nemo` + `--optimize` | LiH | the same, on a nemo reference | 77 s | long |
| `nemo_h2o_canon` | `nemo` | H₂O | `localize canon` | 75 s | long |
| `nemo_h2o_boys` | `nemo` | H₂O | `localize boys` | 83 s | long |
| `nemo_h2o_new` | `nemo` | H₂O | `localize new` | 79 s | long |
| `mp2_he_corr` | `mp2` | He | MP2 correlation energy; fixed `k` | 790 s | verylong |
| `lrcc2_he_excited` | `cc2` | He | CC2 ground state + LR-CC2 excitation | 6744 s | verylong |

The three `nemo_h2o_*` cases are one test split three ways: localization is a
unitary rotation among occupied orbitals, so all three must give the same total
energy. They do, to 2.5e-6 (−76.0681689 / −76.0681714 / −76.0681707), and each
case pins its own value rather than the invariance being asserted in a comment.

`scf_he_hf_mpi` reproduces `scf_he_hf`'s energy to 2.2e-15 on 2 ranks. It gates on
`requires: {"mpi": true, "threads": 16}` and **skips** below that: 2 ranks x (3
workers + 1 comm) is 8 MADNESS threads before MPI's own progress threads, so on an
8-core machine it oversubscribes and never finishes — 7 s on a 96-core node was a
timeout on an 8-core laptop. Its `run.sh` therefore sets the per-rank worker count
from `MPI_WORKERS` rather than inheriting `MAD_NUM_THREADS`, which by convention
means *per job*, not per rank.

The three optimization cases cover the two forms a geometry optimization takes.
`scf_lih_gopt` uses `dft gopt`, which optimizes *inside* one SCF task;
`scf_lih_optimize` and `nemo_lih_optimize` use `--optimize --wf=<scf|nemo>`, the
composable form that is its own task and publishes the optimized geometry for a
later step. The moldft pair
agree to the last digit of the geometry (r = 3.034046 bohr), which is the check
that the first-class optimizer really drives the same MolOpt as the in-SCF path;
`nemo_lih_optimize` lands 1.4e-4 bohr away, the difference between a regularized
and a plain SCF reference at these thresholds. The nemo case is also the
regression test for a wavefunction that is *not* re-solved as the nuclei move —
`Nemo::value` skips the SCF when the orbitals are already converged to the
requested threshold, and if that flag survives a geometry change the optimizer
walks a frozen wavefunction. The symptom is `delta-e` of exactly 0.0 across
iterations, which is why `final_energy` is checked to 1e-6 there rather than the
usual 1e-5: the frozen-wavefunction answer is off by only 9e-6.

`scf_lih_gopt` was the first case to optimize a geometry, and the two keys it
leans on are behavioral rather than energetic: `optimization_results.nsteps` and
`optimization_results.max_gradient`. MolOpt projects translations and rotations
out of the gradient before testing convergence; drop that projection and the raw
net force (8.4e-5 here, an identical dE/dz on both atoms) never reaches
`gtol`, so the run burns all 10 `gmaxiter` cycles at a geometry that stopped
moving after two — `nsteps` 2 → 10 and `max_gradient` 4.8e-7 → 8.4e-5, while the
energy (3.8e-8) and the bond length (3.6e-6 bohr) stay well inside any sane
tolerance. Asserting the energy alone would not have noticed. The case pays for
the 1e-6 protocol rung for the same reason: `gtol` is derived from
`protocol().back()`, and a single 1e-4 rung puts it at 2e-4, *above* the
spurious net force — a criterion that cannot tell a projected gradient from an
unprojected one. `max_gradient` carries `"allow_zero": true` because a converged
gradient is legitimately near zero.

There is deliberately **no plain-CC2 case**: a `calc_type cc2` ground state on
helium measured 1833 s, and `lrcc2_he_excited` already solves that ground state on
its way to the excitation.

## Anatomy of a case

```
scf_h2o_hf/
  scf_h2o_hf.in                       the deck
  run.sh                              one line: exec ${MADQC:-madqc} --wf=scf scf_h2o_hf.in
  check.json                          which results are checked, and how tightly
  reference/scf_h2o_hf.calc_info.json compared numerically against a fresh run
  reference/scf_h2o_hf.out            what a correct run looks like -- never compared
```

**The deck is named after the case.** `madqc` derives its output prefix from the
deck's filename, so `scf_h2o_hf.in` produces `scf_h2o_hf.calc_info.json` and
`scf_h2o_hf.out`. A deck named literally `input` is special-cased to the default
prefix `mad` (see `ParameterManager.hpp`) — hence the convention.

**`run.sh` is one line** and honours `${MADQC}` so the same file works for a user
with an installed binary and for ctest against a build tree. MPI cases use
`${MPIEXEC}` / `${NP}` the same way. A case may `exit 77` to skip itself.

**`reference/*.out` is checked in but never diffed.** It carries wall times,
paths and host details; comparing it would be permanently red. It is there to be
read.

**`check.json` is the assertion list**, which also documents which
`calc_info.json` keys carry the meaningful results:

```json
{
  "description": "Hartree-Fock single point on water",
  "checks": [
    {"key": ["tasks", 0, "scf_total_energy"], "tol": 1e-4},
    {"key": ["tasks", 0, "scf_eigenvalues_a", "vals", 0], "tol": 1e-4},
    {"key": ["tasks", 0, "model"], "tol": 0}
  ]
}
```

`key` is a path of keys and list indices into `<prefix>.calc_info.json`. `tol` is
an absolute tolerance; `0` means "must match exactly", and ints and strings always
compare exactly. A key absent from either file is a failure, not a skip. Optional
`requires` gates a case on resources — `{"threads": 20}`, `{"mpi": true}`,
`{"env": ["MAD_ROOT_DIR"]}` — and turns it into a ctest skip rather than a
failure.

**A reference value no larger than its own tolerance is rejected.** Such a check
passes for every conceivable result, so it survives any regression. Exactly `0.0`
is the common case: not every field of `calc_info.json` is populated by every
workflow — `tasks[0].scf_total_energy` is filled by the nemo/mp2/cc2 paths but
left at `0.0` by the plain `scf` and `response` paths, and several of the older
scripted tests compare precisely that. A dark state's `1e-26` oscillator strength
checked to `1e-3` is just as empty. Compare a key that carries a value, or set
`"allow_zero": true` where the near-zero is the physics (a symmetry-vanishing
dipole component, a gradient at a stationary point).

## Adding a case

No code, four data files:

1. `mkdir <wf>_<system>_<method>/` and write the deck, named after the directory.
2. Write the one-line `run.sh` and `chmod +x`.
3. Write `check.json` — prefer a handful of physically meaningful quantities
   (total energy, an orbital eigenvalue, an excitation, a tensor component) over
   dumping every key.
4. Generate the reference from a run you have inspected:

   ```bash
   mkdir -p /tmp/qc && cd /tmp/qc
   MADQC=/path/to/build/src/apps/madqc/madqc \
     python3 /path/to/madness/bin/run_qctest.py --case /path/to/src/examples/qc/<case> --update
   ```

5. Register it in `CMakeLists.txt` with a tier you have actually measured:
   `add_qctest(<case> madqc "qctest;short;applications")`.

Tiers, by measured wall time:

| Tier | Wall time | In `check-short-madness`? |
|------|-----------|---------------------------|
| `short` | < 10 s | yes |
| `medium` | < 30 s | yes |
| `long` | < 2 min | no |
| `verylong` | everything slower | no |

`verylong` cases typically also carry `requires.threads`, so small machines skip
them instead of timing out.

These are deliberately tight. `check-short-madness` is `ctest -L "short|medium"`,
so every second in those two tiers is a second on every CI run — a QC calculation
that takes minutes belongs in `long` or `verylong` no matter how interesting it
is. Measure before tagging: the same case can differ several-fold between a
laptop and a cluster node (`cc2_he_corr` is ~7 min on a 96-core node), so tag
from the slowest machine you expect to gate on, not the fastest.

## What does *not* belong here

A qctest case is "run one deck, compare declared result keys". Tests that
exercise *behavior* stay as hand-written scripted tests in `src/apps/madqc/`
(registered with `add_scripted_tests`), because expressing them here would mean
growing `check.json` into a programming language:

| Test | Why it stays Python |
|------|--------------------|
| `test_cis_symmetry_h2o.py` | two chained runs — the second restarts from the first with `no_compute=1` and a different `irrep` |
| `test_molresponse_h2o_alpha_beta_z.py`, `test_molresponse_lih_alpha_raman_beta_xyz.py` | flattens the nested v3 `response_properties` object (property → protocol → rows, tensors expanded per element) and rejects a v2-shaped reference outright, so a stale reference fails loudly instead of passing |
| `test_moldft_restart_e2e.py` (currently disabled) | three sequential runs, including recovery from a corrupted checkpoint |

The line to hold: if the assertion is "these numbers came out right", it is a
qctest case. If it is "this sequence of operations behaves correctly", it is a
scripted test.

**Regenerating a reference is a numerical claim.** `--update` overwrites the
checked-in values; a commit that does so should say why the numbers moved and
that the new ones are right. It is not a way to make a red test green.

## Deck gotchas worth knowing before writing one

- **`protocol` ladder vs. fixed `k`.** `scf` and `response` want
  `protocol [1e-4, 1e-6]`. The correlated (`mp2`, `cc2`) and `oep` paths take
  the reference at a *single* resolution: hand them a reference that climbed the
  ladder and the run mixes resolutions and aborts in the tensor layer. Those
  decks set `k` explicitly.
- **A workflow generally needs its own group block** in the deck (even empty,
  as in `nemo_he_hf`), not just `--wf=`.
- **`molecule` and `geometry` are interchangeable** block names
  (`Molecule::GeometryParameters::input_tag`).
- **Every knob in a group is discoverable**: `madqc --print_parameters=dft`,
  `--print_parameters=response`, and so on. Prefer that over guessing.
- **Under MPI, disable the launcher's default binding.** OpenMPI binds each rank
  to one core, so the MADNESS worker thread and the comm thread that every rank
  spawns when `nranks > 1` end up on the same hwthread. The symptom is not
  slowness but a stall: an SCF wedged in `SCF::loadbal` →
  `WorldDCPmapInterface::redistribute` → `WorldGopInterface::fence`, burning CPU
  and making no progress. `mpiexec --bind-to none` (see
  `scf_he_hf_mpi/run.sh`) fixes it; other launchers spell it differently
  (`--cpu-bind` under SLURM, `PARSEC_MCA_bind_threads=0` for PaRSEC).

Full driver documentation: [`src/apps/madqc/README.md`](../../apps/madqc/README.md).
