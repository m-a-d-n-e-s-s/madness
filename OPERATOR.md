# Operator Application Stencil Cutoff Optimization & Coefficient Norm Screening

## Overview

In MADNESS, convolution with Coulomb and Yukawa/BSH kernels represented as sums of separated Gaussians is performed using the non-standard (NS) form of the input function. In 1D, the NS form of an operator has the block representation
```math
\begin{pmatrix} T & C \\ B & A \end{pmatrix}
```
where $T$ maps scaling coefficients to scaling coefficients, $B$ and $C$ couple scaling and wavelet coefficients, and $A$ maps wavelet coefficients to wavelet coefficients. In $D > 1$ dimensions, rather than decomposing the operator into $2^D$ sub-blocks, a single $(2k)^D$ block $R$ (which includes $T$) is applied to the NS-form coefficients at level $n$, and the application of $T$ alone to the $k^D$ block of scaling coefficients is subtracted out (`at.t_term = true`), except at the coarsest level where the telescoping sum terminates (`at.t_term = false`).

Previously, the NS telescoping summation was carried all the way up the tree to level $n = 0$ for every Gaussian term in the separated expansion, subtracting $T$ at all levels $n > 0$ and retaining $T$ only at $n = 0$. However, for a Gaussian with exponent $\alpha_\mu$, once the refinement level $n$ is coarse enough that the operator stencil has narrowed to $\pm 1$ boxes in a given dimension, coarsening further to levels $n < n_{\text{cut}}$ cannot narrow the stencil any further (it remains $\pm 1$ up to level 0). Thus, terminating the NS sum for term $\mu$ at $n_{\text{cut}}(\mu)$, applying the $T$ block at that level (`at.t_term = false`), and skipping term $\mu$ for $n < n_{\text{cut}}(\mu)$ avoids redundant coarse-level applications.

---

## Summary of Changes

### 1. Operator Stencil Cutoff & $T$-Block Application
- **Per-Gaussian Cutoff Level (`src/madness/mra/convolution1d.h`, `src/madness/mra/operator.h`)**:
  - In `GaussianConvolution1D::cutoff_level()`, the 1D cutoff level is computed from the effective cell-scaled exponent $\beta = \alpha_\mu L_d^2$ against the $\pm 1$ box width threshold ($\beta r^2 = 49.0$ at $r = 2^{-n}$):
    ```math
    n_{\text{cut}, d}(\mu) = \begin{cases}
    \lfloor 0.5 \log_2(\alpha_\mu L_d^2 / 49.0) \rfloor & \text{if } \alpha_\mu L_d^2 \ge 49.0 \\
    0 & \text{otherwise}
    \end{cases}
    ```
  - For a $D$-dimensional separated Gaussian term (`ConvolutionND::cutoff_level()`), $n_{\text{cut}}(\mu) = \min_d n_{\text{cut}, d}(\mu)$, representing the lowest level at which the stencil narrows to $\pm 1$ across the dimensions.
  - Periodic remainder operators (`LatticeImages::exclude_home` / `images_only_`) return $n_{\text{cut}} = 0$ so lattice-summed background terms telescope to level 0 as required.
- **Adaptive-Tree Leaf Depth Guard (`src/madness/mra/funcimpl.h`, `src/madness/mra/operator.h`)**:
  - On spatially adaptive trees, diffuse outer branches may terminate at leaf levels $n_{\text{leaf}} < n_{\text{cut}}(\mu)$. Unconditionally skipping levels $n < n_{\text{cut}}(\mu)$ would omit scaling potential contributions on branches that never reach $n_{\text{cut}}(\mu)$.
  - `FunctionImpl::apply` computes the global minimum leaf level `tree_min_level` across all MPI ranks (`world.gop.min`) and passes it statelessly by value through `do_apply` into `SeparatedConvolution::apply`, `apply2`, and `apply2_lowdim`.
  - For each Gaussian term $\mu$, the effective cutoff level is:
    ```math
    n_{\text{eff}}(\mu) = \min\big(n_{\text{cut}}(\mu),\; \text{tree\_min\_level}\big)
    ```
  - When $n_{\text{eff}}(\mu) > 0$:
    - For `source.level() < eff_ncut`: the Gaussian term $\mu$ is skipped (`continue`).
    - At `source.level() == eff_ncut`: `at_mu.t_term = false` (the full $R$ block including $T$ is applied without subtracting $T$).
    - For `source.level() > eff_ncut`: `at_mu.t_term = true` ($T$ is subtracted as usual).
- **Thread Safety & Const-Correctness**:
  - `SeparatedConvolution` stores `cutoff_levels` immutably upon construction and takes `tree_min_level` by value in `apply`, `apply2`, and `apply2_lowdim`, maintaining strict `const`-correctness and thread safety across multi-threaded orbital evaluations.

### 2. Input Function Coefficient Norm Screening
- **Global & Per-Level Coefficient Norm Estimation (`src/madness/mra/funcimpl.h`, `src/madness/mra/mra.h`)**:
  - Added `max_coeff_norm_local()`, `max_coeff_norm()`, and `max_coeff_norm_level()` to `FunctionImpl` and `Function`.
  - `max_coeff_norm()` computes the global supremum Frobenius norm $c_{\max} = \max_{\text{nodes}} \|c\|_F$ over all local tree nodes and MPI ranks via `world.gop.max`.
  - In `apply(op, f)` (`src/madness/mra/mra.h`), `cnorm_max` is evaluated on the redundant/NS tree prior to `apply_only` and passed into `FunctionImpl::do_apply`, which screens separated term applications using `tol / fac / std::max(cmax, 1e-15)`.

### 3. Modified Files
- `src/madness/mra/convolution1d.h`
- `src/madness/mra/operator.h`
- `src/madness/mra/funcimpl.h`
- `src/madness/mra/mra.h`
- `src/madness/mra/testbsh.cc`

---

## Verification & Benchmark Results

### 1. MRA Unit & Regression Test Suites (`MAD_NUM_THREADS=16`)

| Test Target | Scope | Result |
| :--- | :--- | :---: |
| **`testbsh`** | BSH operator convolution across all 9 polynomial orders ($k = 4 \dots 12$) + `max_coeff_norm` checks | **100% PASS** |
| **`testopdir`** | Directional derivative & separated convolutions across all 27 permutations ($m \in [0, 2]^3$) on non-cubic cells and combined operators ($f_{12}, fg_{12}, f_{12}^2, f_{12}^2 g$) | **100% PASS** |
| **`test_images_operator`** | Periodic boundary condition & lattice image convolutions across 12 periodic scenarios | **100% PASS** |
| **`testsuite`** | Full 1D–4D multiwavelet regression test suite ($3.08 \times 10^6$ tasks) | **100% PASS** |

### 2. Single-Rank Water Benchmark (`madqc --geometry=water`, 1 MPI $\times$ 16 threads)
- **Baseline (`master`) Energy**: `-76.067934303751 Ha`
- **Modified (`feature/operator-stencil-cutoff`) Energy**: `-76.067934304252 Ha`
- **Absolute Energy Difference**: $|\Delta E| = 5.02 \times 10^{-10}\text{ Ha}$ (converged in 10 SCF iterations, 37.5s wall time).

---

### 3. Multi-Rank `madqc` Benchmarks (4 MPI Processes $\times$ 8 Threads Each)

All runs below were executed from scratch in clean isolated directories using:
```bash
MAD_NUM_THREADS=8 mpirun -np 4 madqc --geometry=<geom>
```
Comparing baseline (`master`, commit `73ad2cf89`) against `feature/operator-stencil-cutoff` (commit `86949712a`).
For Strontium (`sr`, $Z = 38$), `--dft="aobasis=sto-3g"` is required because the default `6-31g` atomic guess basis only covers elements up to $Z = 36$ (Kr); results for `sr` are reported for both canonical orbitals (`--dft="aobasis=sto-3g; localize canon"`) and default localized orbitals (`--dft="aobasis=sto-3g"`).

#### Table 3A: Converged Total Energies (Before vs. After)

| Geometry | Electrons ($N_e$) | SCF Iters (`1e-4` + `1e-6`) Before / After | Energy Before (`master`) [Ha] | Energy After (`cutoff`) [Ha] | $\Delta E = E_{\text{after}} - E_{\text{before}}$ [Ha] | Relative $\|\Delta E / E\|$ |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`mg`** | 12 | 7 (5+2) / 7 (5+2) | `-199.614535094458` | `-199.614534942790` | $+1.52 \times 10^{-7}$ | $7.6 \times 10^{-10}$ |
| **`ca`** | 20 | 8 (5+3) / 8 (5+3) | `-676.758075526346` | `-676.758076556433` | $-1.03 \times 10^{-6}$ | $1.5 \times 10^{-9}$ |
| **`sr`** (`localize canon`) | 38 | 15 (10+5) / 15 (10+5) | `-3131.545562138397` | `-3131.545562039156` | $+9.92 \times 10^{-8}$ | $3.2 \times 10^{-11}$ |
| **`sr`** (`localize new`) | 38 | 17 (12+5) / 22 (16+6) | `-3131.545552123907` | `-3131.545561561119` | $-9.44 \times 10^{-6}$ | $3.0 \times 10^{-9}$ |
| **`c2h4`** | 16 | 10 (6+4) / 10 (6+4) | `-78.069525275792` | `-78.069525275788` | $+3.31 \times 10^{-12}$ | $4.2 \times 10^{-14}$ |
| **`water2`** | 20 | 10 (6+4) / 10 (6+4) | `-152.141501177296` | `-152.141501177274` | $+2.27 \times 10^{-11}$ | $1.5 \times 10^{-13}$ |
| **`water5`** | 50 | 10 (6+4) / 10 (6+4) | `-380.380408551116` | `-380.380408551152` | $-3.51 \times 10^{-11}$ | $9.2 \times 10^{-14}$ |

#### Table 3B: Wall-Clock & Operator Timings (Before vs. After)

| Geometry | Total Wall Before [s] | Total Wall After [s] | Wall $\Delta t$ (%) | `Apply BSH` Total Before / After [s] | `Apply BSH` Total $\Delta t$ (%) | `Apply BSH` Per-Iter (`1e-4`) Before / After [s] | `Apply BSH` (`1e-4`) $\Delta t$ (%) | `HF exchange` Before / After [s] |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`mg`** | 29.4 | 28.7 | -0.7s (**-2.4%**) | 2.99 / 2.69 | **-10.0%** | 0.258 / 0.220 | **-14.7%** | 15.52 / 15.54 |
| **`ca`** | 83.4 | 81.8 | -1.6s (**-1.9%**) | 9.94 / 8.79 | **-11.6%** | 0.728 / 0.650 | **-10.7%** | 55.72 / 55.58 |
| **`sr`** (`localize canon`) | 470.7 | 470.9 | +0.2s (+0.04%) | 50.84 / 48.58 | **-4.4%** | 1.963 / 1.763 | **-10.2%** | 372.46 / 373.20 |
| **`sr`** (`localize new`)* | 496.9 | 643.3* | +146.4s (+29.5%)* | 60.29 / 83.33* | +38.2%* | 1.762 / 1.556 | **-11.7%** | 386.94 / 496.16* |
| **`c2h4`** | 44.1 | 43.5 | -0.6s (**-1.4%**) | 5.02 / 4.27 | **-14.9%** | 0.422 / 0.340 | **-19.4%** | 23.49 / 23.94 |
| **`water2`** | 59.3 | 61.0 | +1.7s (+2.9%) | 7.88 / 7.18 | **-8.9%** | 0.567 / 0.502 | **-11.5%** | 31.12 / 31.82 |
| **`water5`** | 195.1 | 195.8 | +0.7s (+0.4%) | 19.49 / 19.14 | **-1.8%** | 1.797 / 1.662 | **-7.5%** | 124.87 / 121.98 |

*\*Note on `sr` with default orbital localization (`localize new`):* For spherical `Sr` ($Z=38$) with degenerate $2p, 3p, 3d, 4p$ shells, orbital localization forms symmetry-broken localized hybrids on the Cartesian mesh whose orientation is sensitive to $\sim 10^{-6}$ perturbations, causing the KAIN solver to take 4 extra iterations in the `1e-4` protocol and 1 extra iteration in the `1e-6` protocol (`22` vs `17` total iterations), while per-iteration `Apply BSH` time at `1e-4` is **11.7% faster** (`1.556s` vs `1.762s`). With canonical orbitals (`localize canon`), both before and after runs converge in **15 iterations (10 + 5)** to within $+9.92 \times 10^{-8}\text{ Ha}$ (relative difference $3.2 \times 10^{-11}$), with total `Apply BSH` time reduced from `50.84s` to `48.58s` (**-4.4%** overall, **-10.2%** in the `1e-4` protocol).
