# MADQC — Multiresolution Quantum Chemistry Driver

`madqc` is the unified command-line driver for MADNESS quantum-chemistry
calculations. It parses one input deck, assembles a **workflow** (a chain of
calculation steps), runs them, and writes both a machine-readable result file
(`<prefix>.calc_info.json`) and a human-readable report (`<prefix>.out`).

This is the **user guide**. To *extend* madqc with a new workflow, see
[`WORKFLOWS.md`](WORKFLOWS.md). For the Seawulf build/run loop, see

For **runnable** decks, prefer
[`src/examples/qc/`](../../examples/qc/README.md) over the examples quoted in
this file: each case there is a complete deck plus its reference results, run in
CI under the `qctest` label, so it cannot drift out of date.

---

## Quick start

```bash
# Hartree-Fock / DFT single point (the default workflow)
madqc --wf=scf h2o.in

# show usage / per-workflow help
madqc --help              # general
madqc --help=response     # workflow-specific

# list every knob in a parameter group, then exit
madqc --print_parameters=dft
madqc --print_parameters=response
```

Two files are produced (here `prefix = h2o`, taken from the input file name):

| File | Audience | Contents |
|------|----------|----------|
| `h2o.out` | humans | formatted energies, dipole, excitations, response tensors |
| `h2o.calc_info.json` | machines | full structured results (the source of truth) |

A `h2o/` directory is also created holding per-step work dirs
(`h2o/task_0/moldft/`, `h2o/task_1/molresponse/`, …) with restart archives and
any density/orbital plots.

---

## Workflows

Select with `--wf=<name>` (aliases: `--workflow`, `-w`). Default: `scf`.

| `--wf=` | Steps | Description |
|---------|-------|-------------|
| `scf`   | SCF | Hartree-Fock / DFT single point (`moldft` engine) |
| `nemo`  | SCF | HF/DFT with regularized (nuclear-cusp-free) orbitals |
| `response` | SCF → response | Frequency-dependent polarizability / hyperpolarizability / Raman / excited states — see [`RESPONSE_PROPERTIES.md`](RESPONSE_PROPERTIES.md) |
| `mp2` / `cc2` | SCF → correlated | MP2 / CC2 correlation energies and excitations |
| `cis`   | SCF → CIS | Configuration-interaction-singles excited states (`tdhf` engine) |
| `oep`   | SCF → OEP | Optimized effective potential |

Every multi-step workflow runs the ground-state SCF first, then hands the
reference wavefunction to the downstream method.

Geometry optimization is not a workflow of its own but a flag on one: `--optimize`
optimizes the geometry of whatever reference `--wf` names, so
`madqc --optimize --wf=nemo deck.in` optimizes on a nemo reference and
`--optimize --wf=scf` on a moldft one. It runs as its own task, honours the
`optimization` parameter group, writes `<prefix>_opt.xyz`, and publishes the
optimized geometry for a later step. Only `scf` and `nemo` are optimizable today.

Alternatively `gopt 1` in the `dft` group optimizes *inside* an SCF task; that is
the older form and remains supported. Both drive the same optimizer and give the
same geometry. Working decks: `src/examples/qc/{scf,nemo}_lih_optimize` and
`src/examples/qc/scf_lih_gopt`.

---

## Command-line options

| Option | Effect |
|--------|--------|
| `--wf=<name>` | choose the workflow (default `scf`) |
| `--optimize` | optimize the geometry of that workflow's reference (`scf`, `nemo`); takes no value |
| `--help[=<workflow>]` | usage; optionally workflow-specific examples |
| `--print_parameters=<group>` | print all knobs in a group and exit |
| `--geometry=<name>` | molecule from a built-in name or `.xyz` file |
| `--<group>="k1=v1; k2=v2"` | override knobs on the command line, e.g. `--dft="k=8; econv=1e-6"` |
| `[input_file]` | path to an input deck (plain-text or JSON) |

`<group>` is one of: `dft`, `nemo`, `response`, `cc2`, `cis`, `oep`, `optimization`,
`geometry`.

---

## Input format

An input deck is a set of `group … end` blocks. Block names are the
**parameter-group tags** (which differ slightly from the `--wf`/print-group
names — see the table below). Anything not in the deck takes its default; use
`--print_parameters=<group>` to see all defaults.

```text
dft
    xc          hf
    protocol    [1e-4, 1e-6]
    dconv       1e-4
    localize    new
    dipole      true
end

response
    dipole              true
    dipole.directions   z
    dipole.frequencies  [0.0, 0.02]
    quadratic           true        # also compute hyperpolarizability (beta)
    engine              v3          # optional: select the molresponse_v3 engine
end

molecule
    units     angstrom
    eprec     1e-6
    O   0.0000000   0.0000000   0.1172999
    H   0.0000000   0.7572000  -0.4692000
    H   0.0000000  -0.7572000  -0.4692000
end
```

A JSON input deck with the same group keys is also accepted (the driver
auto-detects JSON vs. plain text).

### Example decks for the other workflows

Each deck below has been run as written. Note the **fixed `k`**: the correlated
and OEP paths take the reference at a single resolution, so these decks set `k`
explicitly rather than giving `protocol` a ladder to climb. (Handing a downstream
correlated method a reference that climbed the ladder mixes resolutions and
aborts; use the ladder for `scf`/`response`, a fixed `k` here.)

**`--wf=nemo`** — HF/DFT with regularized, nuclear-cusp-free orbitals:

```text
dft
    xc       hf
    k        7
    maxiter  10
    econv    1e-5
    dconv    1e-3
end
nemo
end
molecule
    units atomic
    He 0.0 0.0 0.0
end
```

**`--wf=cis`** — configuration-interaction-singles excited states (`tdhf` group):

```text
dft
    xc       hf
    k        6
    maxiter  10
    econv    1e-4
    dconv    1e-3
end
tdhf
    thresh                     1e-3
    maxiter                    10
    nexcitations               2
    guess_excitations          2
    guess_excitation_operators dipole+
end
molecule
    units atomic
    He 0.0 0.0 0.0
end
```

**`--wf=cc2`** (and `--wf=mp2`, same `cc2` group) — correlated energies:

```text
dft
    xc       hf
    k        5
    maxiter  10
    econv    1e-4
    dconv    1e-3
end
cc2
    freeze 1
end
molecule
    units atomic
    He 0.0 0.0 0.0
end
```

**`--wf=oep`** — optimized effective potential:

```text
dft
    xc       hf
    k        7
    maxiter  10
    econv    1e-4
    dconv    1e-3
end
oep
    model       oaep
    oep_maxiter 3
end
molecule
    units atomic
    Be 0.0 0.0 0.0
end
```

Every knob in a group is discoverable with `madqc --print_parameters=<group>`
(`dft`, `nemo`, `tdhf`, `cc2`, `oep`, `response`). The same settings can be given
on the command line instead of in a deck, e.g.
`madqc --molecule=he --wf=cc2 --dft="k=5; econv=1.e-4" --cc2="freeze 1"`, which is
the form the scripted tests in this directory use.

### Input block ↔ parameter group mapping

| Input block (`… end`) | `--print_parameters` group | Used by |
|-----------------------|----------------------------|---------|
| `dft`      | `dft` (alias `scf`) | scf, and the reference of every workflow |
| `molecule` | `geometry` | all |
| `response` | `response` | response |
| `cc2`      | `cc2` | mp2, cc2 |
| `tdhf`     | `cis` | cis |
| `nemo`     | `nemo` | nemo |
| `oep`      | `oep` | oep |

---

## Output: `<prefix>.out` (human-readable)

Rendered at the end of the run from the aggregated results. Example
(SCF + CIS + response):

```text
======================================================================
  MADQC RESULTS SUMMARY
======================================================================

  Task 0 : SCF  (model = scf)
  ----------------------------------------------------------------------
    Molecule         : H2O
    Total energy     :      -76.067441563 Ha       -2069.9005 eV
    Dipole (a.u.)    :  x -0.000010  y +0.000000  z -0.389700
                     : |mu| = 0.389700 a.u. = 0.990519 Debye
    Orbital eps (a)  : -20.5500 -0.9700 -0.9700 -0.7100 -0.5100  (Ha)
    Converged        :  thresh = 1e-06  dconv = 0.0001
    Wall time        :  5.0 s   (1 MPI x 7 threads)

  Task 1 : EXCITED STATES  (model = cis)
  ----------------------------------------------------------------------
    Root  Irrep        E (Ha)      E (eV)   lambda(nm)    f(length)
       1      a      0.851800     23.1787       53.49     0.185600

  Task 2 : RESPONSE
  ----------------------------------------------------------------------
    alpha_zz(w=0.000)  =       8.532200
    beta_zzz(wB=0.000,wC=0.020) =       7.754200
======================================================================
```

The report only prints sections that are present, so it adapts to whatever the
workflow produced.

---

## Output: `<prefix>.calc_info.json` (schema reference)

Top level is a list of task results in execution order:

```json
{ "tasks": [ <task_0>, <task_1>, ... ] }
```

Numeric arrays are stored as **tensor-json**:
`{"dims":[...], "ndim":N, "size":M, "vals":[...]}`. An absent/empty tensor is
`{"ndim":-1, "size":0}` (no `vals`).

### SCF task (`scf`, `nemo`)

| Key | Meaning |
|-----|---------|
| `model` | `"scf"` |
| `properties.energy` | **total SCF energy (Ha)** — use this one |
| `scf_total_energy` | 0 for `moldft`; populated by `nemo`. Prefer `properties.energy` |
| `properties.dipole` | dipole vector `[x,y,z]` (a.u.), tensor-json |
| `properties.gradient` | nuclear gradient (Ha/bohr), tensor-json |
| `scf_eigenvalues_a` / `_b` | orbital energies (Ha), tensor-json |
| `scf_fock_a` / `_b` | Fock matrix, tensor-json |
| `convergence_info.converged_for_thresh` / `_dconv` | the protocol / density thresholds actually met |
| `molecule` | `{symbols, geometry, parameters}` |
| `optimization_results` | geometry-opt summary (`nsteps`, `final_energy`, `max_gradient`, …) |
| `metadata` | `elapsed_time`, `finished_at`, `git_hash`, `host`, `mpi_size`, `nthreads` |

### Excited-state task (`cis`)

| Key | Meaning |
|-----|---------|
| `model` | `"cis"` |
| `excitations[]` | one entry per root |
| `excitations[].omega` | excitation energy (Ha) |
| `excitations[].irrep` | irreducible representation |
| `excitations[].oscillator_strength_length` / `_velocity` | oscillator strengths |
| `excitations[].current_error` | residual at convergence |
| `nfreeze` | frozen-core count |

### Correlated task (`mp2`, `cc2`)

| Key | Meaning |
|-----|---------|
| `model` | `"mp2"` or `"cc2"` |
| `correlation_energy` | correlation energy (Ha) |
| `<model>_total_energy` | SCF + correlation (Ha) |
| `excitations[]` | LR-CC2 excited states (when requested) |

### Response task

| Key | Meaning |
|-----|---------|
| `type` | `"response"` |
| `properties.response_properties[]` | computed tensor components |
| `…[].property` | `"polarizability"` or `"hyperpolarizability"` |
| `…[].component` | axis labels, e.g. `["z","z"]` (α) or `["z","z","z"]` (β) |
| `…[].freqB` (and `freqC` for β) | driving frequencies (a.u.) |
| `…[].value` | the tensor component |
| `properties.raman_spectra` | per-mode Raman intensities (when computed) |
| `properties.vibrational_analysis` | normal-mode frequencies (when computed) |
| `metadata` | per-state / per-protocol convergence detail |

---

## Visualization

`madqc` writes volumetric plots when the corresponding knobs are set in the
`dft` group, and indexes everything it produced in a manifest.

Plot knobs (`dft` group):

| Knob | Effect |
|------|--------|
| `plotdens` | write the converged density (`total_density`, `spin_density`) |
| `plotcoul` | write the total Coulomb potential (`coulomb`) |
| `plotcube` | also emit Gaussian `.cube` (for Avogadro/VMD) alongside `.dx` |
| `plotlo` / `plothi` | range of molecular orbitals to plot (`amo-NNNNN`, `bmo-NNNNN`) |
| `npt_plot` | grid points per dimension (default 101) |
| `plot_cell` | `lo hi` per dimension (default: whole box) |

Files land in the per-task work dir, e.g. `h2o/task_0/moldft/total_density.cube`.

> **Fixed (was: plotdx segfault) — runtime-confirmed:** `SCF::do_plots` used to
> pass an empty plot cell to the OpenDX writer (`plotdx`, `mraimpl.h`) whenever
> the `plot_cell` knob was unset — the common case — because
> `CalculationParameters::plot_cell()` returns a copy, so the defaulting
> assignment hit a discarded temporary. plotdx then dereferenced an empty cell and
> crashed. `do_plots` now builds a local cell defaulted to the simulation cell, so
> `plotdens`/`plotcoul` work without an explicit `plot_cell` (confirmed on h2o:
> `total_density.dx` + `coulomb.dx` written, no crash). `plotdx` itself now also
> `MADNESS_CHECK_THROW`s on an empty/ill-shaped cell, so a future caller passing
> one gets a clear error instead of a segfault.

After the run, the driver scans the output tree and writes
**`<prefix>.viz_manifest.json`** — a single discovery index:

```json
{
  "prefix": "h2o",
  "molecule": { "formula": "H2O", "symbols": [...], "geometry": [...] },
  "results": { "calc_info": "h2o.calc_info.json", "summary": "h2o.out" },
  "artifacts": [
    { "task": 0, "method": "moldft", "kind": "density",
      "format": "cube", "path": "h2o/task_0/moldft/total_density.cube" },
    { "task": 0, "method": "moldft", "kind": "orbital_a", "index": 2,
      "format": "dx", "path": "h2o/task_0/moldft/amo-00002.dx" }
  ]
}
```

`kind` is one of `density`, `potential`, `orbital_a`/`orbital_b`, `field`;
`format` is `cube`, `dx`, or `line`. If no plot knobs were enabled, no manifest
is written (and the driver says so).

Viewing the artifacts:
- **ParaView / VMD / Avogadro** open `.cube` directly; `.dx` via ParaView or OpenDX.
- **gecko** (Python) discovers a run from its `.calc_info.json` for analysis/plots.

## Worked example: water SCF

`h2o.in`:
```text
dft
    xc        hf
    protocol  [1e-4, 1e-6]
    dipole    true
end
molecule
    units angstrom
    O   0.0  0.0  0.1173
    H   0.0  0.7572 -0.4692
    H   0.0 -0.7572 -0.4692
end
```

Run and read the result:
```bash
madqc --wf=scf h2o.in
cat h2o.out                                   # human report
```

Extract a value from the JSON (Python):
```python
import json
d = json.load(open("h2o.calc_info.json"))
scf = d["tasks"][0]
print("E(SCF) =", scf["properties"]["energy"], "Ha")
print("dipole =", scf["properties"]["dipole"]["vals"])
```

---

## See also

- [`RESPONSE_PROPERTIES.md`](RESPONSE_PROPERTIES.md) — how to compute each response property (α, β, Raman, excited states)
- [`../molresponse/README.md`](../molresponse/README.md) — the v3 response engine: status + capability matrix
- [`WORKFLOWS.md`](WORKFLOWS.md) — how to add a new workflow
- `madqc --print_parameters=<group>` — authoritative, always-current knob list
- `doc/quantum.md` (repo root) — broader MADNESS quantum-chemistry overview
