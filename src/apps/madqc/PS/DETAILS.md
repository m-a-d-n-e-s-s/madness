# Technical Details: Pseudopotential Capability in `madqc`

This document describes the investigation, bug fixes, capabilities, benchmark methodology, and numerical verification results for pseudopotential calculations in `madqc` (PR [#793](https://github.com/m-a-d-n-e-s-s/madness/pull/793)).

---

## 1. Capabilities of the MADNESS Pseudopotential Code

### 1.1 Types of Pseudopotentials Supported
1. **Goedecker–Teter–Hutter (GTH) / Hartwigsen–Goedecker–Hutter (HGH) Pseudopotentials (`gth_pseudopotential.h`, `gth_pseudopotential.cc`)**:
   - Separable dual-space norm-conserving pseudopotentials consisting of an analytical local potential $V_{\text{loc}}(r)$ and short-ranged Gaussian-polynomial nonlocal projectors $p_i^{lm}(\mathbf{r})$ up to $l_{\max} = 3$ ($s, p, d, f$ channels, up to 3 projectors per $(l,m)$ channel):
     $$V_{\text{loc}}(r) = -\frac{Z_{\text{eff}}}{r}\operatorname{erf}\!\left(\frac{r}{\sqrt{2}\,r_{\text{loc}}}\right) + \exp\!\left(-\frac{1}{2}\left(\frac{r}{r_{\text{loc}}}\right)^2\right)\sum_{k=1}^{4} C_k \left(\frac{r}{r_{\text{loc}}}\right)^{2k-2}$$
     $$\hat{V}_{\text{nl}} = \sum_{l=0}^{l_{\max}} \sum_{m=-l}^{l} \sum_{i,j=1}^{3} |p_i^{lm}\rangle\, h_{ij}^{(l)}\, \langle p_j^{lm}|$$
   - Supports **purely local** pseudopotentials ($l_{\max} = -1$, e.g., `H`, `He` in `gth.xml` and `H`, `He`, `Li`, `Be` in `gth_pbe.xml`), diagonal nonlocal projectors, and off-diagonal coupled projectors ($h_{01}^{(l)}, h_{02}^{(l)}, h_{12}^{(l)} \neq 0$, e.g., `Ne`, `Mg`, `Ar`, and transition metals).
   - Supports **all-atom** pseudopotential calculations (`psp_calc true` in `dft` or `geometry`) as well as **mixed all-electron / pseudopotential** calculations (`ps<Symbol>` prefix on selected atoms in `geometry`, such as `psO`, `psC`, `psMg`, while keeping `H` all-electron).

2. **Model Core Potentials (`corepotential.h`, `corepotential.cc`)**:
   - Invoked via `core_type mcp` in the `dft` block using parameter files `mcp` and `mcp_guess` in `src/madness/chem/`.

### 1.2 Supported DFT Exchange-Correlation Functionals and Parameter Libraries
- **`lda`** (Vosko–Wilk–Nusair `LDA_X` + `LDA_C_VWN`): Uses `gth.xml` (Goedecker–Teter–Hutter 1996 / Hartwigsen–Goedecker–Hutter 1998 parameters, $Z=1\text{–}86$).
- **`pbe`** (`GGA_X_PBE` + `GGA_C_PBE` via Libxc): Uses `gth_pbe.xml` (Krack 2005 / CP2K GTH-PBE parameters, $Z=1\text{–}86$).
- `CalculationParameters::select_gth_file(xc_str, user_psp_file)` automatically selects `gth_pbe.xml` when `xc` contains `"pbe"` (case-insensitive) and `gth.xml` otherwise, unless overridden by `psp_file <filename>` in the `geometry` block.

### 1.3 Initial Guess, Valence Electron Count, and Potential Assembly
- During parameter initialization, `Molecule::update_pseudo_atom_charges(gth_file)` updates `Atom::q` to the valence charge $Z_{\text{eff}}$ for every atom with `pseudo_atom == true` and sets `Molecule::pure_ae = false` (and `Molecule::psp_calc = true` if all atoms are pseudo-atoms).
- `CalculationParameters::set_derived_values(molecule)` sums `molecule.total_nuclear_charge()` using the updated $Z_{\text{eff}}$ values so `nelec`, `nalpha`, and `nbeta` reflect only the active valence + all-electron electrons.
- In `SCF::initial_guess`, the minimal-basis AO superposition guess (`atomic_guess_density` / `load_esp_rom`) skips core shells for pseudo-atoms (retaining only the top $Z_{\text{eff}}/2$ valence orbitals per pseudo-atom) while retaining all shells on all-electron atoms.
- In `SCF::make_nuclear_potential`:
  - If all atoms are pseudo-atoms (`psp_calc == true`), `vnuc` is solely the GTH local potential `gthpseudopotential->vlocalpot()`.
  - If some atoms are pseudo-atoms and some are all-electron (`pure_ae == false && psp_calc == false`), `PotentialManager::make_nuclear_potential` builds the smoothed $-Z/r$ nuclear potential on all-electron atoms only (`if (molecule.get_pseudo_atom(i)) continue;`) and adds `gthpseudopotential->vlocalpot()` on the pseudo-atoms.
- During each SCF iteration (`SCF::apply_potential`), `gthpseudopotential->apply_potential` evaluates the nonlocal projector action $\hat{V}_{\text{nl}}|\psi_n\rangle$ and accumulates the nonlocal pseudopotential energy $E_{\text{nl}} = \sum_n f_n \langle \psi_n | \hat{V}_{\text{nl}} | \psi_n \rangle$.

---

## 2. Bugs Discovered and Fixed

Seven bugs and regressions were identified and fixed to restore full pseudopotential functionality in `madqc`:

1. **`madqc` Workflow Regression in `ParameterManager` / `SCF` (`ParameterManager.hpp`, `CalculationParameters.h`, `SCF.cc`, `Applications.hpp`)**:
   - **Issue**: When `madqc` was refactored to construct `ParameterManager` before `SCF`, `cparam.set_derived_values(molecule)` computed `nelec`, `nalpha`, and `nbeta` *before* `SCF::SCF` propagated `psp_calc` to `molecule` or updated `Atom::q` to $Z_{\text{eff}}$. As a result, `psp_calc true` left `nelec` equal to the all-electron count, and `ps<Symbol>` atoms had `Atom::q = 0` when `nelec` was computed! Furthermore, `SCF::SCF` only updated `molecule` on MPI rank 0 without broadcasting `this->molecule` to other ranks, and hardcoded `"gth.xml"` instead of selecting `"gth_pbe.xml"` for PBE.
   - **Fix**: Added `CalculationParameters::select_gth_file` and `Molecule::update_pseudo_atom_charges(gth_file)`, invoked them in `ParameterManager::set_derived_values` (before `cparam.set_derived_values(molecule)`) and in `SCF::SCF`, broadcast `this->molecule` across all MPI ranks, and passed `gth_file` to `make_pseudo_potential(world, gth_file)`.

2. **JSON Serialization & `from_json` Stripping of Pseudo-Atom State (`molecule.cc`, `molecule.h`)**:
   - **Issue**: `Molecule::print` and `Molecule::insert_symbols_and_geometry` emitted `"O"` instead of `"psO"` for pseudo-atoms, and `Molecule::from_json` re-parsed the geometry strings without preserving `pseudo_atom`, `pure_ae`, `psp_calc`, `psp_file`, or $Z_{\text{eff}}$. Because `ParameterManager` round-trips `Molecule` through JSON during setup, any `psO`/`psC` prefix was lost before `SCF` started.
   - **Fix**: Updated `Molecule::print` and `Molecule::insert_symbols_and_geometry` to prefix pseudo-atom symbols with `"ps"`, and updated `Molecule::from_json` and `SCFApplication::checkpoint_hamiltonian` to preserve and restore `psp_file`, `pure_ae`, `psp_calc`, and $Z_{\text{eff}}$.

3. **Crash on Semicore (`-sc`) Entries in `gth_pbe.xml` (`gth_pseudopotential.cc`)**:
   - **Issue**: `gth_pbe.xml` contains semicore entries such as `<atom symbol="Mg-sc">`, `<atom symbol="Au-sc">`, and `<atom symbol="Ga-sc">`. `get_charge_from_file` and `load_pseudo_from_file` passed every XML `symbol` directly to `symbol_to_atomic_number`, which threw a fatal exception on any `-sc` symbol.
   - **Fix**: Added `if (symbol && std::string(symbol).find('-') != std::string::npos) continue;` in `get_charge_from_file` and `load_pseudo_from_file`.

4. **Crash on Purely Local Pseudopotentials (`lmax = -1`, e.g., `Be` in `gth_pbe.xml`) (`gth_pseudopotential.h`)**:
   - **Issue**: In `gth_pbe.xml`, `Be` ($Z_{\text{eff}}=4$), `Li` ($Z_{\text{eff}}=3$), `He` ($Z_{\text{eff}}=2$), and `H` ($Z_{\text{eff}}=1$) are purely local GTH pseudopotentials with `lmax = -1` and zero nonlocal projectors. When `atoms_with_projectors` was empty (`natoms == 0`), `GTHPseudopotential::apply_potential` divided by `localproj.size() == 0` and indexed into an empty `nlmat(0,0)`, causing a segmentation fault.
   - **Fix**: Initialized `enl = 0.0` and returned `vpsi` immediately when `natoms == 0` in both `apply_potential` and `apply_potential_simple`.

5. **Off-Diagonal Nonlocal Projector Coupling Coefficients $h_{01}^{(l)}, h_{02}^{(l)}, h_{12}^{(l)}$ (`gth_pseudopotential.cc`)**:
   - **Issue**: In the original GTH/HGH formalism (Hartwigsen, Goedecker, and Hutter, *Phys. Rev. B* **58**, 3641 (1998)), `gth.xml` and `gth_pbe.xml` store the diagonal elements $h_{00}^{(l)}, h_{11}^{(l)}, h_{22}^{(l)}$, while off-diagonal elements are defined analytically when not explicitly given in XML:
     - $l = 0$: $h_{01}^{(0)} = -\frac{1}{2}\sqrt{\frac{3}{5}}\,h_{11}^{(0)}$, $h_{02}^{(0)} = \frac{1}{2}\sqrt{\frac{5}{21}}\,h_{22}^{(0)}$, $h_{12}^{(0)} = -\frac{1}{2}\sqrt{\frac{100}{63}}\,h_{22}^{(0)}$
     - $l = 1$: $h_{01}^{(1)} = -\frac{1}{2}\sqrt{\frac{5}{7}}\,h_{11}^{(1)}$, $h_{02}^{(1)} = \frac{1}{2}\sqrt{\frac{35}{11}}\,h_{22}^{(1)}$, $h_{12}^{(1)} = -\frac{1}{2}\frac{14}{\sqrt{11}}\,h_{22}^{(1)}$
     - $l = 2$: $h_{01}^{(2)} = -\frac{1}{2}\sqrt{\frac{7}{9}}\,h_{11}^{(2)}$, $h_{02}^{(2)} = \frac{1}{2}\sqrt{\frac{63}{143}}\,h_{22}^{(2)}$, $h_{12}^{(2)} = -\frac{1}{2}\frac{18}{\sqrt{143}}\,h_{22}^{(2)}$
     Previously, `load_pseudo_from_file` hardcoded $h_{01} = h_{02} = h_{12} = 0$ for all $l$, causing errors of $\sim 0.14\text{ Ha}$ on `Ne` (`h11 = -1.07624528` for $l=0$) and $\sim 0.01\text{ Ha}$ on `Mg` (`h11 = 1.32994091` for $l=0$ in LDA, `0.92534825` in PBE).
   - **Fix**: Updated `load_pseudo_from_file` to read `h01`, `h02`, `h12` from XML if present and otherwise compute them from the analytical HGH relations for $l = 0, 1, 2$.

6. **Real Spherical Harmonic Normalization Bug for $l=2, m=+1$ (`itmp2 == 7`) (`gth_pseudopotential.h`)**:
   - **Issue**: In `ProjRLMFunctor::operator()`, the $l=2, m=+1$ ($d_{xz}$) component (`itmp2 == 7`) was missing the factor of $\sqrt{2}$ present in $d_{yz}$ (`itmp2 == 5`) and $d_{xy}$ (`itmp2 == 4`), breaking rotational invariance for $d$-projector atoms.
   - **Fix**: Multiplied `itmp2 == 7` by `std::sqrt(2.0)`.

7. **GTH XML Installation & Path Resolution (`CMakeLists.txt`, `gth_pseudopotential.cc`)**:
   - **Issue**: `gth.xml` and `gth_pbe.xml` resided in `src/apps/moldft/` and were neither installed into `MRA_CHEMDATA_DIR` nor located when running `madqc` outside `src/apps/moldft/`. Also, `nuclear_repulsion_derivative` and `nuclear_repulsion_second_derivative` in `molecule.cc` used `atomic_number` instead of `q` for pseudo-atoms.
   - **Fix**: Added `gth.xml` and `gth_pbe.xml` to `install(FILES ...)` in `src/madness/chem/CMakeLists.txt`, added `resolve_gth_filename()` in `gth_pseudopotential.cc`, and fixed `nuclear_repulsion_derivative` and `nuclear_repulsion_second_derivative` in `molecule.cc` to use `atoms[i].q` when `atoms[i].pseudo_atom` is true.

---

## 3. Benchmark Methodology

### 3.1 `madqc` Settings
- **Atoms (`Be`, `Ne`)**: Cubic box half-width `L = 30.0` Bohr, wavelet order `k = 8`, adaptive refinement protocol `[1e-4, 1e-6, 1e-8]`, density convergence `dconv = 1e-7`.
- **Mixed AE/PSP Molecular Suite (`mg`, `water`, `water2`, `methane`, `ethylene`)**: Exact Cartesian geometries (in Bohr) from `src/madness/chem/structure_library` with `no_orient true`, pseudopotentials on heavy atoms (`psMg`, `psO`, `psC`) and all-electron hydrogen (`H`), `L = 30.0` Bohr, `k = 8`, `protocol [1e-4, 1e-6]`, `dconv = 1e-6`.
- **All-Atom PSP Molecular Test (`ethylene_allpsp`)**: Same `ethylene` geometry from `structure_library` with `psp_calc true` so that **both** `C` ($Z_{\text{eff}} = 4, l_{\max} = 1$) and `H` ($Z_{\text{eff}} = 1, l_{\max} = -1$) use their GTH pseudopotentials (`L = 30.0` Bohr, `k = 8`, `protocol [1e-4, 1e-6]`, `dconv = 1e-6`).

### 3.2 NWChem `nwpw` (`pspw`) Settings
- Executed with `/home/rjh/Devel/nwchem-7.3.1/bin/LINUX64/nwchem` using the `pspw` module (`task pspw energy`).
- **Boundary Conditions**: Free-space cluster / aperiodic convolution (`boundary_conditions aperiodic`) so electrostatic and local pseudopotential terms are evaluated without periodic images, matching MADNESS's free-space boundary conditions.
- **Grid & Supercell**:
  - `Be` and `Mg`: $128 \times 128 \times 128$ FFT grid with `SC = 30.0` Bohr ($E_{\text{cut}} = 89.8\text{ Ha}$; `SC = 25.0` Bohr for `be_pbe`).
  - `Ne` (`ne_lda`, `ne_pbe`): Two-stage calculation on `SC = 16.0` Bohr starting from $96 \times 96 \times 96$ ($E_{\text{cut}} = 177.65\text{ Ha}$) followed by $192 \times 192 \times 192$ ($E_{\text{cut}} = 710.61\text{ Ha}$, grid spacing $\Delta x = 0.0833\text{ Bohr}$, `time_step 1.0`, `tolerances 1.0d-9 1.0d-9`). See Section 3.3 below for the plane-wave convergence analysis on `Ne`.
  - `water`, `methane`, `ethylene`, `ethylene_allpsp`: $128 \times 128 \times 128$ FFT grid with `SC = 20.0` Bohr ($E_{\text{cut}} = 202.1\text{ Ha}$, $\Delta x = 0.15625\text{ Bohr}$).
  - `water2`: $128 \times 128 \times 128$ FFT grid with `SC = 22.0` Bohr ($E_{\text{cut}} = 167.0\text{ Ha}$).
- **Pseudopotentials**:
  - Pseudopotential atoms (`Be`, `Ne`, `Mg`, `O`, `C`, and `H` in `ethylene_allpsp`) use HGH/GTH format (`psp_type = 1`) `.psp` files constructed with the **exact parameters** from MADNESS's `gth.xml` (for `xc lda`) and `gth_pbe.xml` (for `xc pbe96`).
  - All-electron `H` in `water`, `water2`, `methane`, and `ethylene` uses an HGH-format Coulomb potential (`H_ae.psp`: $Z_{\text{ion}} = 1, r_{\text{loc}} = 0.00275\text{ Bohr}, C_i = 0, h_{ij} = 0$), which NWChem evaluates analytically in Fourier space as $-4\pi/G^2 \exp(-\frac{1}{2}G^2 r_{\text{loc}}^2)$, matching MADNESS's smoothed all-electron nuclear potential on `H`.

### 3.3 Analysis of `Ne` Plane-Wave Cutoff Convergence and Exact 1D Radial Verification

In initial NWChem `pspw` runs using `SC = 20.0` Bohr and `ngrid 128 128 128` ($\Delta x = 0.15625\text{ Bohr}, E_{\text{cut}} = 202.13\text{ Ha}$), `ne_lda` and `ne_pbe` exhibited discrepancies of $3.16\text{ mHa}$ and $3.95\text{ mHa}$ relative to `madqc`.

To determine whether the error resided in `madqc` or in the uniform plane-wave basis of NWChem `pspw`:
1. **Exact 1D Radial Kohn–Sham Chebyshev Spectral Benchmark**:
   We solved the spherical radial Kohn–Sham equations for `Ne` with the exact GTH-LDA and GTH-PBE Hamiltonians using a parity-symmetric Chebyshev spectral collocation solver ($N = 701$ points on $R_{\max} = 10.0\text{ Bohr}$). Comparing `madqc` (`k = 8, thresh = 1e-8`) against the exact 1D radial solution:

   | Quantity (Ha) | `Ne` GTH-LDA (Exact 1D) | `Ne` GTH-LDA (`madqc`) | `Ne` GTH-PBE (Exact 1D) | `Ne` GTH-PBE (`madqc`) |
   | :--- | ---: | ---: | ---: | ---: |
   | Total Energy ($E_{\text{total}}$) | `-34.85667256` | `-34.85667201` | `-35.05245104` | `-35.05245086` |
   | Kinetic ($T$) | `28.38126155` | `28.38125889` | `28.18645919` | `28.18644288` |
   | Local PSP ($V_{\text{loc}}$) | `-90.62689960` | `-90.62689642` | `-90.43412396` | `-90.43409204` |
   | Nonlocal PSP ($V_{\text{nl}}$) | `2.15542813` | `2.15542743` | `2.03282670` | `2.03282929` |
   | Coulomb ($J$) | `31.03994009` | `31.03994018` | `31.14730825` | `31.14728688` |
   | Exchange-Correlation ($E_{\text{xc}}$) | `-5.80640272` | `-5.80640209` | `-5.98492122` | `-5.98491786` |
   | $\epsilon_{2s}$ | `-1.32730072` | `-1.32730122` | `-1.39513326` | `-1.39513643` |
   | $\epsilon_{2p}$ | `-0.49748515` | `-0.49748503` | `-0.48289489` | `-0.48289730` |

   `madqc` agrees with the exact 1D radial solution to **$+0.55\ \mu\text{Ha}$ for LDA** and **$+0.18\ \mu\text{Ha}$ for PBE**!

2. **Origin of the Plane-Wave Basis Truncation Error on `Ne` in NWChem `pspw`**:
   Among all elements tested (`Be`, `Mg`, `C`, `O`, `Ne`), `Ne` sits at the far right of the second period ($Z_{\text{eff}} = 8$) and has by far the **hardest GTH pseudopotential**, with compact Gaussian radii $r_{\text{loc}} = 0.1900\text{ Bohr}$, $r_0 = 0.1795\text{ Bohr}$ ($0.1761\text{ Bohr}$ in PBE), and a strong repulsive $s$-channel projector ($h_{00}^{(0)} = +28.51\text{ Ha}$). At $\Delta x = 0.15625\text{ Bohr}$ ($E_{\text{cut}} = 202.1\text{ Ha}$), the grid spacing is nearly equal to $r_0$, leaving a $+3.16\text{ mHa}$ (LDA) and $+3.95\text{ mHa}$ (PBE) plane-wave basis incompleteness error:
   - $E_{\text{cut}} = 50.53\text{ Ha}$: $E_{\text{LDA}} = -34.66230314\text{ Ha}$ ($\Delta E = +194.37\text{ mHa}$)
   - $E_{\text{cut}} = 101.07\text{ Ha}$: $E_{\text{LDA}} = -34.81172127\text{ Ha}$ ($\Delta E = +44.95\text{ mHa}$)
   - $E_{\text{cut}} = 151.60\text{ Ha}$: $E_{\text{LDA}} = -34.84661503\text{ Ha}$ ($\Delta E = +10.06\text{ mHa}$)
   - $E_{\text{cut}} = 202.13\text{ Ha}$: $E_{\text{LDA}} = -34.85350735\text{ Ha}$ ($\Delta E = +3.16\text{ mHa}$)
   - $E_{\text{cut}} = 710.61\text{ Ha}$ (`SC = 16.0`, `ngrid 192 192 192`, $\Delta x = 0.0833\text{ Bohr}$): $E_{\text{LDA}} = -34.85668876\text{ Ha}$ ($\Delta E = -0.017\text{ mHa}$) and $E_{\text{PBE}} = -35.05250345\text{ Ha}$ ($\Delta E = -0.053\text{ mHa}$).

3. **All-Electron `H` vs. GTH Pseudopotential `H` in Molecules (`ethylene` vs. `ethylene_allpsp`)**:
   Similarly, in `water`, `water2`, `methane`, and `ethylene` (where `H` is treated as all-electron with a $-1/r$ Coulomb potential), NWChem's uniform plane-wave grid at $E_{\text{cut}} = 202.1\text{ Ha}$ underbinds each all-electron `H` atom by $\sim 0.3\text{–}0.4\text{ mHa}$ due to the sharp nuclear cusp on `H`, producing a total discrepancy of $-1.19\text{ mHa}$ (LDA) and $-1.30\text{ mHa}$ (PBE) on `ethylene`. When `H` also uses its smooth GTH pseudopotential (`ethylene_allpsp_lda` and `ethylene_allpsp_pbe`), `madqc` and NWChem `pspw` agree to **$+6.46\ \mu\text{Ha}$ (LDA)** and **$+2.25\ \mu\text{Ha}$ (PBE)**, and `madqc`'s wall time drops by nearly $3\times$ (from `107.1s` to `35.8s` in LDA) with half as many wavelet coefficients (`2.45e6` vs. `4.70e6`).

> [!NOTE]
> **Note on `Be` GTH-PBE (`be_pbe`)**: In `gth_pbe.xml`, `Be` is a 4-electron purely local pseudopotential ($Z_{\text{eff}} = 4, l_{\max} = -1, r_{\text{loc}} = 0.325\text{ Bohr}$) with a large repulsive barrier ($C_2 = +17.279$) at $r \sim r_{\text{loc}}$ that pushes the $2s$ orbital energy ($\epsilon_{2s} = -0.2061\text{ Ha}$) very close to the $2p$ orbital energy, breaking $2s$–$2p$ degeneracy in the opposite direction of a bare Coulomb potential. Starting from `madqc`'s spherical $1s^2 2s^2$ AO guess preserves the physical $(1s)^2(2s)^2\ (^1S)$ state ($E = -14.535436\text{ Ha}$), which agrees with an exact 1D spherical radial Kohn-Sham solver ($E = -14.536212\text{ Ha}$) to $0.78\text{ mHa}$. In contrast, NWChem `pspw` starts from an unsymmetrized 3D random plane-wave guess without angular momentum constraints and spontaneously breaks spherical symmetry by mixing $2s$ and $2p^2$ character ($E = -14.623508\text{ Ha}$).

---

## 4. Complete Numerical Results

### 4.1 Total Energy and Nonlocal Pseudopotential Energy Comparison (`madqc` vs. NWChem `pspw`)

All energies are in Hartrees ($\text{Ha}$).

| Test Case | XC | PSP / AE Specification | $E_{\text{total}}$ (`madqc`) | $E_{\text{total}}$ (`NWChem`) | $\Delta E_{\text{total}}$ (`madqc - NWChem`) | $E_{\text{nl}}$ (`madqc`) | $E_{\text{nl}}$ (`NWChem`) |
| :--- | :--- | :--- | ---: | ---: | ---: | ---: | ---: |
| `be_lda` | LDA | `psBe` ($Z_{\text{eff}}=2$) | `-0.99124457` | `-0.99124509` | `+0.00000052` | `0.23908933` | `0.23908705` |
| `be_pbe` | PBE | `psBe` ($Z_{\text{eff}}=4, l_{\max}=-1$) | `-14.53543649` | `-14.53621201`* | `+0.00077552`* | `0.00000000` | `0.00000000` |
| `ne_lda` | LDA | `psNe` ($Z_{\text{eff}}=8$) | `-34.85667201` | `-34.85668876` | `+0.00001675` | `2.15542743` | `2.15541540` |
| `ne_pbe` | PBE | `psNe` ($Z_{\text{eff}}=8$) | `-35.05245086` | `-35.05250345` | `+0.00005259` | `2.03282929` | `2.03285566` |
| `mg_lda` | LDA | `psMg` ($Z_{\text{eff}}=2$) | `-0.83302448` | `-0.83302537` | `+0.00000089` | `0.19798418` | `0.19798145` |
| `mg_pbe` | PBE | `psMg` ($Z_{\text{eff}}=2$) | `-0.78554523` | `-0.78554835` | `+0.00000312` | `0.11239781` | `0.11245962` |
| `water_lda` | LDA | `psO` ($Z_{\text{eff}}=6$), `2` AE `H` | `-17.18646296` | `-17.18585569` | `-0.00060727` | `1.14500062` | `1.14484248` |
| `water_pbe` | PBE | `psO` ($Z_{\text{eff}}=6$), `2` AE `H` | `-17.22421418` | `-17.22340527` | `-0.00080891` | `1.14844619` | `1.14836213` |
| `water2_lda` | LDA | `2 psO` ($Z_{\text{eff}}=6$), `4` AE `H` | `-34.38495077` | `-34.38315683` | `-0.00179394` | `2.28610942` | `2.28606754` |
| `water2_pbe` | PBE | `2 psO` ($Z_{\text{eff}}=6$), `4` AE `H` | `-34.45656437` | `-34.45421688` | `-0.00234749` | `2.29245510` | `2.29261075` |
| `methane_lda` | LDA | `psC` ($Z_{\text{eff}}=4$), `4` AE `H` | `-8.03824351` | `-8.03696971` | `-0.00127380` | `0.43639087` | `0.43643367` |
| `methane_pbe` | PBE | `psC` ($Z_{\text{eff}}=4$), `4` AE `H` | `-8.07840368` | `-8.07701913` | `-0.00138455` | `0.43232153` | `0.43237216` |
| `ethylene_lda` | LDA | `2 psC` ($Z_{\text{eff}}=4$), `4` AE `H` | `-13.69641777` | `-13.69522395` | `-0.00119382` | `0.91657666` | `0.91667292` |
| `ethylene_pbe` | PBE | `2 psC` ($Z_{\text{eff}}=4$), `4` AE `H` | `-13.72929315` | `-13.72799118` | `-0.00130197` | `0.90865018` | `0.90865911` |
| `ethylene_allpsp_lda` | LDA | `2 psC` ($Z_{\text{eff}}=4$), `4 psH` ($Z_{\text{eff}}=1$) | `-13.69574202` | `-13.69574848` | `+0.00000646` | `0.91657412` | `0.91654668` |
| `ethylene_allpsp_pbe` | PBE | `2 psC` ($Z_{\text{eff}}=4$), `4 psH` ($Z_{\text{eff}}=1$) | `-13.72852252` | `-13.72852477` | `+0.00000225` | `0.90865311` | `0.90867941` |

*\*For `be_pbe`, the $(1s)^2(2s)^2\ (^1S)$ reference energy from the 1D spherical radial Kohn-Sham solver is `-14.53621201 Ha`; NWChem `pspw` without angular symmetry constraints converges to a symmetry-broken $2s/2p$ mixed state at `-14.62350849 Ha`.*

---

### 4.2 Detailed Energy Components (`madqc` vs. NWChem `pspw`)

All values in Hartrees ($\text{Ha}$). In each cell `A / B`, `A` is `madqc` and `B` is NWChem `pspw`.

| Test Case | Kinetic ($T$) | Local PSP + Nuclear ($V_{\text{loc}}$) | Nonlocal PSP ($V_{\text{nl}}$) | Coulomb ($J$) | Exchange-Correlation ($E_{\text{xc}}$) | Nuclear Repulsion ($E_{\text{nuc}}$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| `be_lda` | `0.36065597` / `0.36065241` | `-1.91276905` / `-1.91276073` | `0.23908933` / `0.23908705` | `0.69128142` / `0.69127670` | `-0.36950225` / `-0.36950052` | `0.00000000` / `0.00000000` |
| `be_pbe` | `12.85774087` / `13.04152245` | `-31.80415516` / `-31.98633274` | `0.00000000` / `0.00000000` | `7.06125052` / `7.07484575` | `-2.65027272` / `-2.75354395` | `0.00000000` / `0.00000000` |
| `ne_lda` | `28.38125889` / `28.38148370` | `-90.62689642` / `-90.62729630` | `2.15542743` / `2.15541540` | `31.03994018` / `31.04016163` | `-5.80640209` / `-5.80645319` | `0.00000000` / `0.00000000` |
| `ne_pbe` | `28.18644288` / `28.18654457` | `-90.43409204` / `-90.43375928` | `2.03282929` / `2.03285566` | `31.14728688` / `31.14674308` | `-5.98491786` / `-5.98488747` | `0.00000000` / `0.00000000` |
| `mg_lda` | `0.26491062` / `0.26491308` | `-1.55692757` / `-1.55692611` | `0.19798418` / `0.19798145` | `0.57672941` / `0.57672726` | `-0.31572114` / `-0.31572105` | `0.00000000` / `0.00000000` |
| `mg_pbe` | `0.27250860` / `0.27245637` | `-1.40399767` / `-1.40401004` | `0.11239781` / `0.11245962` | `0.54246856` / `0.54247198` | `-0.30892252` / `-0.30892628` | `0.00000000` / `0.00000000` |
| `water_lda` | `13.83670233` / `13.83520731` | `-56.58949930` / `-56.58709505` | `1.14500062` / `1.14484248` | `21.47323595` / `21.47297224` | `-4.14310251` / `-4.14298269` | `7.09119995` / `7.09120001` |
| `water_pbe` | `13.79532817` / `13.79140406` | `-56.50804394` / `-56.50196437` | `1.14844619` / `1.14836213` | `21.49494953` / `21.49317922` | `-4.24609407` / `-4.24558632` | `7.09119995` / `7.09120001` |
| `water2_lda` | `27.67662965` / `27.67265846` | `-136.04268827` / `-136.03672650` | `2.28610942` / `2.28606754` | `54.37499837` / `54.37462937` | `-8.29390766` / `-8.29369359` | `25.61390772` / `25.61390787` |
| `water2_pbe` | `27.60372933` / `27.59553530` | `-135.89584164` / `-135.88239490` | `2.29245510` / `2.29261075` | `54.42780408` / `54.42349254` | `-8.49861897` / `-8.49736839` | `25.61390772` / `25.61390787` |
| `methane_lda` | `6.68579932` / `6.68184411` | `-36.79103286` / `-36.78515489` | `0.43639087` / `0.43643367` | `15.19161138` / `15.19039944` | `-3.09162748` / `-3.09110744` | `9.53061526` / `9.53061540` |
| `methane_pbe` | `6.70787785` / `6.70379154` | `-36.78380966` / `-36.77837485` | `0.43232153` / `0.43237216` | `15.21306444` / `15.21256239` | `-3.17847309` / `-3.17798578` | `9.53061526` / `9.53061540` |
| `ethylene_lda` | `10.89461636` / `10.89093610` | `-69.65364698` / `-69.64818331` | `0.91657666` / `0.91667292` | `29.45543925` / `29.45427705` | `-4.75330584` / `-4.75282975` | `19.44390278` / `19.44390305` |
| `ethylene_pbe` | `10.89136829` / `10.88806718` | `-69.58472277` / `-69.58077704` | `0.90865018` / `0.90865911` | `29.47606452` / `29.47639395` | `-4.86455614` / `-4.86423742` | `19.44390278` / `19.44390305` |
| `ethylene_allpsp_lda` | `10.88470978` / `10.88448944` | `-69.64285043` / `-69.64240135` | `0.91657412` / `0.91654668` | `29.45513644` / `29.45485763` | `-4.75321472` / `-4.75314394` | `19.44390278` / `19.44390305` |
| `ethylene_allpsp_pbe` | `10.88051459` / `10.88064641` | `-69.57298546` / `-69.57347607` | `0.90865311` / `0.90867941` | `29.47571653` / `29.47608013` | `-4.86432407` / `-4.86435770` | `19.44390278` / `19.44390305` |

---

### 4.3 Canonical Kohn-Sham Orbital Eigenvalues (`madqc`)

| Test Case | Occupied Orbital Eigenvalues $\epsilon_i$ (Ha) |
| :--- | :--- |
| `be_lda` | `-0.20568755` |
| `be_pbe` | `-3.90296780`, `-0.20613816` |
| `ne_lda` | `-1.32730122`, `-0.49748503` ($\times 3$) |
| `ne_pbe` | `-1.39513643`, `-0.48289730` ($\times 3$) |
| `mg_lda` | `-0.17558933` |
| `mg_pbe` | `-0.16607189` |
| `water_lda` | `-0.93372820`, `-0.49649386`, `-0.34215903`, `-0.27335106` |
| `water_pbe` | `-0.93951696`, `-0.49395024`, `-0.33919915`, `-0.26787302` |
| `water2_lda` | `-0.95889247`, `-0.90762371`, `-0.51933176`, `-0.47181980`, `-0.37172963`, `-0.32133683`, `-0.29816199`, `-0.25123964` |
| `water2_pbe` | `-0.96332380`, `-0.91345116`, `-0.51553840`, `-0.46886344`, `-0.36720597`, `-0.31785183`, `-0.29136816`, `-0.24562712` |
| `methane_lda` | `-0.62249778`, `-0.34720523` ($\times 3$) |
| `methane_pbe` | `-0.62575369`, `-0.34655670` ($\times 3$) |
| `ethylene_lda` | `-0.68750405`, `-0.52022468`, `-0.42273132`, `-0.37760080`, `-0.31128787`, `-0.25461475` |
| `ethylene_pbe` | `-0.68998977`, `-0.52243633`, `-0.42066109`, `-0.37455429`, `-0.31243454`, `-0.24807320` |
| `ethylene_allpsp_lda` | `-0.68748396`, `-0.52015969`, `-0.42267853`, `-0.37758028`, `-0.31123654`, `-0.25462323` |
| `ethylene_allpsp_pbe` | `-0.68997239`, `-0.52237251`, `-0.42060966`, `-0.37453759`, `-0.31238509`, `-0.24808590` |
