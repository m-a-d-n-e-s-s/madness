# Running Pseudopotential Calculations with `madqc`

`madqc` supports **Goedecker–Teter–Hutter (GTH)** dual-space separable norm-conserving pseudopotentials (including Hartwigsen–Goedecker–Hutter, HGH, relativistic extensions) as well as **Model Core Potentials (MCP)**. Pseudopotentials replace the singular $-Z/r$ nuclear attraction and tightly bound core electrons with a smooth effective potential and short-ranged nonlocal projectors, substantially reducing the number of orbitals and the multiresolution grid refinement needed near heavy nuclei.

---

## 1. Quick Start: Specifying Pseudopotentials in `madqc`

There are two ways to enable GTH pseudopotentials in a `madqc` input file:

### Option A: Mixed All-Electron / Pseudopotential Calculation (Per-Atom `ps` Prefix)
To apply pseudopotentials only to specific atoms (for example, heavy atoms while keeping hydrogen all-electron), prefix the atomic symbol in the `geometry` block with `ps` (case-insensitive, e.g., `psO`, `psC`, `psMg`):

```text
dft
  xc pbe
  L 30.0
  k 8
  protocol [1e-4, 1e-6]
  dconv 1e-6
end

geometry
  units bohr
  no_orient true
  psO  0.00000000000000  0.00000000000000  0.71002389091700
  H   -1.42650808202133  0.00000000000000 -0.35501194545850
  H    1.42650808202133  0.00000000000000 -0.35501194545850
end
```

In this mode:
- Atoms prefixed with `ps` (`psO`) have their nuclear charge replaced by the valence charge $Z_{\text{eff}}$ from the GTH parameter file ($Z_{\text{eff}} = 6$ for oxygen, removing the $1s^2$ core) and experience the GTH local + nonlocal pseudopotential.
- Standard atomic symbols (`H`) are treated as **all-electron** atoms with their smoothed point-nuclear Coulomb attraction and core electrons intact.
- Nuclear-nuclear repulsion, nuclear dipole moments, and nuclear gradients automatically use $Z_{\text{eff}}$ on pseudo-atoms and $Z$ on all-electron atoms.

---

### Option B: All-Atom Pseudopotential Calculation (`psp_calc true`)
To apply GTH pseudopotentials to **all** atoms in the molecule without modifying the atomic symbols in `geometry`, set `psp_calc true` inside the `dft` block (or inside the `geometry` block):

```text
dft
  xc lda
  psp_calc true
  L 30.0
  k 8
  protocol [1e-4, 1e-6]
  dconv 1e-6
end

geometry
  units bohr
  no_orient true
  Be 0.0 0.0 0.0
end
```

---

## 2. Supported DFT Exchange-Correlation Functionals and Parameter Files

GTH pseudopotentials in MADNESS are parameterized for specific exchange-correlation functionals:

| `xc` Setting in `dft` | Default GTH XML File | Parameter Source | Elements Covered |
| :--- | :--- | :--- | :--- |
| `lda` | `gth.xml` | Goedecker–Teter–Hutter (1996) / Hartwigsen–Goedecker–Hutter (1998) LDA (Padé/VWN5) parameters | $Z = 1\text{ (H)}$ to $Z = 86\text{ (Rn)}$ |
| `pbe` | `gth_pbe.xml` | Krack (2005) / CP2K GTH-PBE parameters | $Z = 1\text{ (H)}$ to $Z = 86\text{ (Rn)}$ |

### Automatic Selection & Custom Override (`psp_file`)
- By default (`psp_file none` in `geometry`), `madqc` automatically selects `gth_pbe.xml` when `xc` contains `pbe` (case-insensitive) and `gth.xml` otherwise.
- You can explicitly override the GTH XML file in the `geometry` block using `psp_file`:
  ```text
  geometry
    psp_file /path/to/custom_gth.xml
    ...
  end
  ```
- **Data Directory Lookup**: `madqc` searches for `gth.xml` and `gth_pbe.xml` in:
  1. The current working directory (or an explicit path if `psp_file` contains `/`),
  2. The directory specified by the `MRA_CHEMDATA_DIR` environment variable,
  3. The compiled-in `MRA_CHEMDATA_DIR` installation directory (`src/madness/chem`), and
  4. `src/apps/moldft/` in the source tree.

---

## 3. Model Core Potentials (`core_type mcp`)

In addition to GTH pseudopotentials, `madqc` also supports **Model Core Potentials (MCP)** (stored in `src/madness/chem/mcp` and `mcp_guess`) for all-atom core-potential calculations:

```text
dft
  xc lda
  core_type mcp
  ...
end
```

---

## 4. Directory Contents and Benchmark Suite

This directory (`src/apps/madqc/PS`) contains benchmark inputs, outputs, and verification scripts comparing `madqc` against **NWChem `nwpw` (`pspw`)** (aperiodic plane-wave DFT with identical GTH pseudopotentials on heavy atoms and Coulomb potential on all-electron H):

- `README.md`: This user guide.
- `DETAILS.md`: Technical documentation of the pseudopotential implementation, bugs fixed, capabilities, and complete numerical comparison tables.
- `madqc/`: Input (`.in`) and reference output (`.out`) files for `madqc` across 8 test configurations (`be`, `ne`, `mg`, `water`, `water2`, `methane`, `ethylene` with AE H, and `ethylene_allpsp` with all-atom GTH pseudopotentials) for both `lda` and `pbe` (16 calculations total).
- `nwchem/`: Input (`.nw`) and reference output (`.out`) files for NWChem `pspw`.
- `run_and_compare.py`: Script to run `madqc` on the benchmark suite (or compare existing outputs) and verify numerical agreement against the reference `madqc` and NWChem results.

### Threading & Avoiding Core Oversubscription (`MAD_NUM_THREADS`)

MADNESS uses spin-locks in its task-queue thread pool, so oversubscribing CPU cores (or letting BLAS/OpenMP spawn nested threads) can cause severe slowdowns. When running `madqc` directly, always restrict the MADNESS thread count via `MAD_NUM_THREADS` and set `OMP_NUM_THREADS=1` and `MKL_NUM_THREADS=1`:
```bash
export MAD_NUM_THREADS=4
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
```
`run_and_compare.py` sets these automatically (controllable via `--threads`, default `4`).

### Running the Verification Script

To compare the existing outputs against the reference values without re-running `madqc`:
```bash
python3 run_and_compare.py --compare-only
```

To run `madqc` on all benchmark inputs (or selected tests) and compare against the reference results:
```bash
python3 run_and_compare.py --madqc /path/to/madqc --threads 4
# Or run specific systems:
python3 run_and_compare.py --madqc /path/to/madqc --threads 4 --tests be_lda ne_lda ethylene_allpsp_lda ethylene_allpsp_pbe
```
