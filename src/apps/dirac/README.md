# Dirac Hartree-Fock Solver (`DFdriver`)

`DFdriver` is a 4-component molecular Dirac Hartree-Fock (DHF) solver implemented in the MADNESS (Multiresolution Adaptive Numerical Environment for Scientific Simulation) framework. It computes relativistic electronic structures using multiwavelet adaptive representations.

---

## Overview

`DFdriver` solves the 4-component Dirac-Fock equations:
$$\begin{pmatrix} V + J - K & c\, \boldsymbol{\sigma} \cdot \mathbf{p} \\ c\, \boldsymbol{\sigma} \cdot \mathbf{p} & V + J - K - 2mc^2 \end{pmatrix} \begin{pmatrix} \psi^L \\ \psi^S \end{pmatrix} = \varepsilon \begin{pmatrix} \psi^L \\ \psi^S \end{pmatrix}$$

Key features:
- **Initial Guess from `moldft`**: Initializes 4-component spinors from a converged non-relativistic ground-state calculation in `moldft` via the kinetic balance condition:
  $$\psi^L = \phi_{\text{NR}}, \quad \psi^S \approx \frac{-i}{2c} (\boldsymbol{\sigma} \cdot \boldsymbol{\nabla}) \phi_{\text{NR}}$$
- **BSH Integral Operator Formulation**: Employs the bound-state Helmholtz (BSH) integral representation to solve for orbital updates without finite-difference grid artifacts.
- **KAIN Acceleration**: Uses the Krylov Accelerated Inexact Newton (KAIN) method for robust convergence.
- **Flexible Nuclear Models**: Supports Gaussian, Fermi, and point-charge nuclear distributions.
- **Restart Capability**: Can restart from `moldft` restart files or resume previous `DFdriver` runs.

---

## Compilation

From your build directory:

```bash
cmake ../madness -GNinja -DENABLE_MPI=OFF \
  -DLAPACK_LIBRARIES="-L/home/rjh/Devel/mTxm/lib -lmTxmq -L/opt/arm/armpl_23.10_gcc-12.2/lib/ -larmpl_lp64 -lamath -lastring"

ninja src/apps/dirac/DFdriver
```

---

## Basic Workflow

Running a Dirac-Fock calculation is a two-step process:

### Step 1: Run `moldft` for the Non-Relativistic Ground State

Run `moldft` to generate the molecular orbital restart archive (default prefix: `mad.restartdata.00000`):

```bash
# Example: Helium atom
MAD_NUM_THREADS=20 ./src/apps/moldft/moldft --geometry=he
```

Or for a custom geometry using an input deck (e.g. Beryllium):
```
dft
  xc hf
end

geometry
  be 0.0 0.0 0.0
end
```
```bash
MAD_NUM_THREADS=20 ./src/apps/moldft/moldft
```

### Step 2: Run `DFdriver`

Create an `input` file in the working directory:

```
DiracFock
  archive mad.restartdata
  small 1e-5
  thresh 1e-6
  k 8
  kain
  maxsub 5
  no_save
end
```

Run `DFdriver`:

```bash
MAD_NUM_THREADS=20 ./src/apps/dirac/DFdriver
```

By default, `DFdriver` reads `input` in the current working directory. You can specify a different input file via:
```bash
MAD_NUM_THREADS=20 ./src/apps/dirac/DFdriver -input my_input
```

---

## Input Options

All parameters are specified within the `DiracFock ... end` block. Lines starting with `#` are treated as comments.

| Parameter | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `archive <file>` | `string` | *(required)* | Path or name of the input restart archive. Accepts `mad.restartdata`, `mad.restartdata.00000`, or `mad`. |
| `job <n>` | `int` | `0` | Calculation type. `0`: DF on occupied orbitals only. |
| `thresh <val>` | `double` | `1e-6` | Wavelet truncation/refinement threshold. Automatically updates default `dconv`. |
| `dconv <val>` | `double` | `1e-6` | Density convergence threshold ($\lVert \Delta \rho \rVert$). |
| `k <order>` | `int` | `8` | Multiwavelet polynomial order. If different from archive, orbitals are automatically projected. |
| `small <val>` | `double` | `1e-5` | Smallest length scale to be resolved. |
| `max_iter <n>` | `int` | `20` | Maximum number of SCF iterations. |
| `min_iter <n>` | `int` | `2` | Minimum number of SCF iterations. |
| `kain` | flag | `false` | Enable KAIN nonlinear accelerator. |
| `maxsub <n>` | `int` | `10` | Maximum subspace size for KAIN. |
| `maxrotn <val>` | `double` | `0.25` | Maximum orbital rotation step allowed by KAIN. |
| `nucleus <n>` | `int` | `0` | Nuclear charge distribution: `0` = Gaussian, `1` = Fermi, `2` = Point. |
| `restart` | flag | `false` | Resume from a previous `DFdriver` run rather than from `moldft`. |
| `savefile <file>`| `string` | `DFrestartdata`| Name of the archive to write restart data to. |
| `no_save` | flag | `false` | Disable saving restart data at the end of each iteration. |
| `lb_iter <n>` | `int` | `20` | Interval (in iterations) between adaptive load balancing passes. |
| `speed_of_light`| `double` | `137.035999...`| Speed of light in atomic units (CODATA 2022). |
| `bohr_rad <val>` | `double` | `52917.72...`  | Bohr radius in fm (CODATA 2022). |
| `lineplot` | flag | `false` | Generate 1D lineplots of large and small spinor components along the x-axis. |
| `no_compute` | flag | `false` | Skip SCF iterations and exit after setup. |

---

## Restarting Calculations

### From `moldft`
Provide the base name or chunk filename of the `moldft` calculation:
```
DiracFock
  archive mad.restartdata
  ...
end
```
Both `archive mad.restartdata` and `archive mad.restartdata.00000` are accepted.

### From a Previous `DFdriver` Calculation
To continue an interrupted or unconverged `DFdriver` calculation:
```
DiracFock
  archive DFrestartdata
  restart
  kain
  max_iter 20
end
```
Ensure the previous run did not specify `no_save`.

---

## Example Inputs

### 1. Standard Helium Atom (`DFinput_sample`)
```
DiracFock
  archive mad.restartdata
  small 1e-5
  thresh 1e-6
  k 8
  kain
  maxsub 5
  no_save
end
```

### 2. High-Precision Calculation with Fermi Nuclear Model
```
DiracFock
  archive /path/to/calc.restartdata
  thresh 1e-7
  dconv 1e-7
  k 10
  nucleus 1
  kain
  maxsub 8
  savefile my_df_checkpoint
end
```

