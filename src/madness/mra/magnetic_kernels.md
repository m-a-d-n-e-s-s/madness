# Magnetic Green's-function kernels

This page documents the code in `src/madness/mra/magnetic_*.h`. It is the integral-equation
machinery for an electron in a uniform magnetic field: four kernels for (T − E)⁻¹, one interface, and
two front ends (znemo, pymadness).

## Context

### The problem

The free bound-state Helmholtz kernel (T − E)⁻¹ does not describe an electron in a field B. With the
free kernel, znemo's default (`kernel bsh`) has to put the paramagnetic term (B/2)L_z and the
diamagnetic term B²ρ²/8 on the right-hand side, and both grow with the field. The kernels here invert
the magnetic kinetic operator itself:

```
T_B = ½ (p + A)² = −½∇² + (B/2) L_z + (B²/8) ρ²,     A = (B/2)(−ỹ, x̃, 0)
```

### Conventions

- atomic units; B along +z; π = p + A, so the lowest Landau level has m ≤ 0;
- gauge origin O = (o_x, o_y), with x̃ = x − o_x and ỹ = y − o_y (znemo and the FT kernel use O = 0);
- inf spec T_B = |B|/2, and every kernel requires E < |B|/2.

### How all four kernels are built

Each kernel is a Laplace transform of a heat kernel that factorizes exactly:

```
G(E) = (T − E)⁻¹ = ∫₀^∞ dτ e^{Eτ} e^{−τT}      (E < inf spec T)
```

- The resolvent and the Laplace transform are explained in `notes/resolvent_kernel_laplace.md` of the
  article repository (see References).
- The τ integral becomes a quadrature in ln τ.
- Each node is a short product of Gaussian convolutions, phase multiplications and masks.

## Files

| File | Content |
|---|---|
| `magnetic_primitives.h` | Shared pieces: Gaussian convolutions on one or two axes of a 3D function (low-dimensional apply), one-term anisotropic Gaussian operators (2D or 3D), the bilinear phase e^{iκ(x₀−o₀)(x₁−o₁)} with two refinement policies, Mehler and splitting coefficients, and the τ quadratures. |
| `magnetic_kernel.h` | `MagneticKernelParameters`, the abstract `MagneticKernel` (`apply`, `bottom`, `absorbs_paramagnetic`, `absorbs_diamagnetic`, `gauge`, `vector_potential`, gauge transformations). |
| `magnetic_kernel_chain.h` | `MagneticKernelChain`: real-space splitting chain, symmetric or Landau gauge. |
| `magnetic_kernel_ft.h` | `MagneticKernelFT`: Fourier approach. |
| `magnetic_kernel_mehler.h` | `MagneticKernelMehler`: diamagnetic (harmonic-oscillator) kernel. |
| `magnetic_kernels.h` | Includes all of the above; `make_magnetic_kernel(world, params)`. |
| `test_magnetic_kernels.cc` | Unit test (`ctest -R mra_magnetic`) and a driver for comparisons: `mode=make / apply / compare`. |
| `chem/znemo.{h,cc}` | Front end: `kernel` in the `complex` block. |
| `pymadness/.../bind_magnetic.cpp`, `pymadness/examples/magnetic_kernels.py` | Python front end. |

Library changes the kernels rely on, both small and in `funcimpl.h` / `operator.h` / `mra.h`:

- **Low-dimensional apply on full tensors** for any operator dimension below the function dimension:
  `SeparatedConvolution::norm_lowdim` and `apply2_lowdim_full`, with the displacement loop continued
  outward while contributions are significant. Used by the chain kernels.
- **`hartree_product` / `project_out` with mixed dimensions** (2D × 1D ↔ 3D). Used by the FT kernel.

## The kernels

| `kind` | Contains | Gauge | Paramagnetic term | τ quadrature | Origin |
|---|---|---|---|---|---|
| `chain_symmetric` | T_B | symmetric | in the kernel | Gauss–Legendre in ln τ | `MagneticGreensFunction3D` in `hydrogen_magnetic_pymadness.py` (article repo, pymadness branch `magnetic`) |
| `chain_landau` | T_B | Landau, A_L = (−Bỹ, 0, 0) | in the kernel | Gauss–Legendre in ln τ | `Gop` of `h2p_landau.py --gauge landau` |
| `ft` | T_B | symmetric, O = 0 | in the kernel | trapezoid in ln τ | `MagneticGreen` in `magnetic_green.h`, branch `magnetic-fourier` (A. Armen) |
| `mehler` | −½∇² + B²ρ̃²/8 | symmetric | **right-hand side** | Gauss–Legendre in ln τ | new in 3D; 2D `MehlerGreensFunction2D` in `landau_mehler_pymadness.py` |

### Real-space chain (`chain_symmetric`, `chain_landau`)

The transverse heat kernel factorizes exactly at every τ
(`notes/magnetic_kernel_realspace_splitting.md`, `notes/mra_magnetic_scheme.md`):

```
e^{−τT⊥} = e^{−a π_y²} e^{−s π_x²} e^{−a π_y²},    a = tanh(Bτ/2)/(2B),  s = sinh(Bτ)/(2B)
```

Each factor e^{−cπ²} is a normalized Gaussian of variance 2c, conjugated by a phase. The free z motion
is merged into the x convolution:

```
e^{−τT} f = e^{−iBx̃ỹ/2} G_y(2a) [ e^{+iBx̃ỹ} G_xz(2s, τ) [ e^{−iBx̃ỹ} G_y(2a) [ e^{+iBx̃ỹ/2} f ]]]
```

**Chain.**
- The bracket without its outer pair e^{∓iBx̃ỹ/2} is the Landau-gauge propagator, with χ = (B/2)x̃ỹ and
  φ_L = e^{+iχ}φ_sym. So `chain_landau` is the same code without the two outer phase multiplications.
  This removes the London plane wave from off-axis orbitals; see `notes/mra_magnetic_scheme.md` §4 and
  the He₂ notes.
- The chain is propagated incrementally from node to node: S_k = e^{−(τ_k−τ_{k−1})T} S_{k−1}.
- Each step is split into sub-steps with |B| Δτ ≤ `max_btau` (default 2), and adjacent y half-steps are
  merged.
- S carries the Landau zero point e^{|B|τ/2}.
- Beyond τ = `nlll`/|B| (default 28), only the free z motion is propagated: the excited Landau levels
  are damped by e^{−nlll}.

**Axis order.** The code runs internally in the axis order (y, x, z), exactly as the Python original
did (`mapdim` at entry and exit). Then the y convolution is a 1D low-dimensional apply on axis 0 and
the x+z convolution a 2D one on axes (1, 2).

### Fourier approach (`ft`)

T_B = T⊥ − ½∂_z², and T⊥ separates in the y wavenumber q, by the Mehler factorization with frequency
B in x for each q:

```
U⊥ = e^{iBxy/2} (1/2π) ∫dq e^{iq(y−y')} g(x−X_q) N(x−x') g(x'−X_q) e^{−iBx'y'/2}
X_q = −q/B,  g = e^{−m u²},  N = √(c/π) e^{−c d²},  m = (B/2) tanh(Bτ/2),  c = B/(2 sinh Bτ)
```

**Steps.**
1. Multiply by e^{−iBxy/2}. The phase is kept only on the cylinder ρ < `ft_rcut` (erfc taper of width
   `ft_rtaper`).
2. Fourier project in y, giving F_q(x, z) (`project_out`).
3. Per q and node: x mask, then a 2D (x, z) Gaussian, then x mask.
4. Assemble Σ_q e^{iqy} Ψ_q (`hartree_product`), then multiply by e^{+iBxy/2}.

**Settings and costs.**
- y is treated as periodic with period P, and the modes are processed in batches (`ft_qbatch`).
- The cost scales with the number of q modes, about (|B| r_ext/2 + `ft_qextra`) P/π.
- The default `ft_rcut` assumes the nuclei lie near the z axis. Pass `ft_rcut` and `ft_rtaper` for
  off-axis nuclei at strong field, as Armen's B = 425 jobs did (rcut 1.0, rtaper 0.1).
- Armen's `magnetic/README.md` and `magnetic/docs/parameters_physics.pdf` on `magnetic-fourier`
  explain every parameter.

### Mehler kernel (`mehler`)

T_HO = −½∇² + (B²/8)ρ̃² is a 2D oscillator of frequency w = |B|/2, plus the free z motion. Its heat
kernel factorizes at rank 1 (`notes/mehler_kernel_factorization.md`):

```
⟨r|e^{−τT_HO}|r'⟩ = (c/π)(2πτ)^{−1/2} g(ρ̃) e^{−c|ρ−ρ'|² − (z−z')²/2τ} g(ρ̃'),
g = e^{−m ρ̃²},  m = (w/2) tanh(wτ/2),  c = w/(2 sinh wτ)
```

so each node is one 3D anisotropic Gaussian convolution between two masks.

**Properties.**
- inf spec T_HO = |B|/2 = inf spec T_B.
- On functions with L_z = 0 about O, the Mehler kernel equals G_B.
- The paramagnetic term is not in the kernel. znemo adds 0.5 B L_z φ to the right-hand side.
- The notes report slow convergence at strong field for this route (150–245 iterations);
  `kernel_nbroaden` and a level shift are the knobs.

### Shared primitives and the original policies

The kernels share one implementation of each building block. Where the two original codes differed,
both policies are kept as options, and each kernel defaults to its own. This keeps the original numbers
bitwise (see Verification):

| Primitive | Policies / parameters | Used by |
|---|---|---|
| `mul_bilinear_phase` | `poly_error`: refine where the best degree-(k−1) fit of the phase misses `tolfac`·truncate_tol relative to the box norm (one recursive sweep); `box_phase`: refine boxes within ρ < r_ext whose phase varies by more than `dphi`, with an optional erfc taper | chain / ft |
| `gauss_conv_1d`, `gauss_conv_2d` | normalized Gaussians on outer axes, `nbroaden` passes | chain |
| `gauss_product_operator<NDIM>` | unnormalized anisotropic one-term operator, exponent floor 1e−10 | ft (2D), mehler (3D) |
| `gauss_legendre_log` | `nodes`, `tmin`, tmax = −ln(`eps`)/(|B|/2 − E) | chain, mehler |
| `trapezoid_log` | `ft_h`, `ft_tmin`, tmax = `ft_tmax_factor`/(|B|/2 − E) | ft |
| `mehler_m`, `mehler_log_c`, `splitting_coefficients` | stable for Bτ → 0 and Bτ → ∞ | all |

## Using the kernels

### C++

```cpp
#include <madness/mra/magnetic_kernels.h>
FunctionDefaults<3>::set_cubic_cell(-L, L);      // cubic cell; k and thresh as usual
MagneticKernelParameters p;
p.kind = "chain_landau";
p.B = 10.0;
auto K = make_magnetic_kernel(world, p);         // sets FunctionDefaults<1>, <2> to match
complex_function_3d g = K->apply(f, E);          // (T_k - E)^{-1} f, E < K->bottom()
```

### znemo: the `complex` block

`kernel bsh`, the default, leaves znemo unchanged. With a magnetic kernel:

- **Kinetic terms.** The kinetic, paramagnetic and diamagnetic terms are in the Green's function (for
  `mehler`, the paramagnetic term stays on the right-hand side).
- **Energy and Fock matrix.** They use ½⟨π_dφ_i|π_dφ_j⟩, with π = −i∇ + A in the kernel's gauge.
- **Spin Zeeman.** The spin Zeeman term ±B/2 moves into the energy argument of G.
- **No NCF, no box.** No nuclear correlation factor and no confining box are used.
- **Kernel energies.** As in the original drivers, the kernel energies E_i are started from the Fock
  matrix and updated by E − Re⟨Vψ|ψ − ψ̃⟩/‖ψ̃‖². They are kept below |B|/2 − `kernel_gap`
  (`kernel_gap_start` in the first iteration).
  - Guess orbitals often have Fock energies above |B|/2.
  - Near the threshold the range of G(E) along B grows as 1/√(2(|B|/2 − E)), and so does the memory.
    With E = B/2 − 1e−3 at B = 10, a single application exceeded 45 GB.
- **Output.** `calc_info.json` gets `kernel` and `kernel_energy`, which is Σ_i E_i − (J − K) + E_nuc,
  the energy the original drivers reported.

| Parameter | Default | Meaning |
|---|---|---|
| `kernel` | `bsh` | `bsh`, `ft`, `chain_symmetric`, `chain_landau`, `mehler` |
| `thresh` | −1 (= econv) | truncation threshold, decoupled from the energy convergence |
| `nuclear_smoothing` | −1 (from eprec) | smoothing radius c of the nuclear potential (MADNESS `smoothed_potential`) |
| `spin_zeeman` | true | include ±B/2 |
| `guess_type`, `guess2` | `default`, `sigma` | `magnetic_lcao`: magnetically translated, Landau-squeezed 1s LCAO, second orbital σ_u or π = (x − iy)σ_g |
| `kernel_gap`, `kernel_gap_start` | 1e−3, 0.05 | E ≤ B/2 − gap |
| `kernel_nodes`, `kernel_tmin`, `kernel_nlll`, `kernel_nbroaden`, `phase_tolfac` | 40, 1e−4, 28, 0, 0.1 | chain and mehler |
| `ft_h`, `ft_tmin`, `ft_tmax_factor`, `ft_qextra`, `ft_rcut`, `ft_rtaper`, `ft_dphi` | 0.25, 1e−8, 45, 10, derived, derived, 3 | ft |

With a magnetic kernel znemo also uses truncate mode 0 and initial level 2, and a Coulomb operator
precision of 0.1·thresh, as the original codes did.

Example, triplet H₂ at B = 10 with the bond perpendicular to B:

```
dft
  k 10
  l 20
  econv 3e-7
  dconv 5.5e-3
  maxiter 45
  nopen 2
  lo 1e-4
end
complex
  physical_B 10
  kernel chain_landau
  thresh 1e-4
  nuclear_smoothing 1e-2
  spin_zeeman false
  guess_type magnetic_lcao
end
geometry
  units atomic
  H  0.4 0.0 0.0
  H -0.4 0.0 0.0
end
```

Map from the original command-line options:

| Original | znemo |
|---|---|
| `h2p_madness.py --nodes --tmin --nlll --phase-tolfac --c --gap` | `kernel_nodes`, `kernel_tmin`, `kernel_nlll`, `phase_tolfac`, `nuclear_smoothing`, `kernel_gap` |
| `h2p_madness.py --guess2 sigma/pi` | `guess_type magnetic_lcao`, `guess2` |
| `magneticatom h= tmin= tmax_factor= qextra= rcut= rtaper= dphi=` | `ft_h`, `ft_tmin`, `ft_tmax_factor`, `ft_qextra`, `ft_rcut`, `ft_rtaper`, `ft_dphi` |
| `magneticatom eprec=` | `nuclear_smoothing` = (eprec/1.3)^{1/3} for Z = 1 |

### Python

```python
import pymadness as mad
K = mad.make_magnetic_kernel(world, mad.MagneticKernelParameters(kind="ft", B=10.0))
g = K.apply(f, E)                       # or K.apply([f1, f2], [E1, E2])
```

See `src/pymadness/examples/magnetic_kernels.py`.

## Verification

**1. Against the original implementations, one G application.**
- Set-up: B = 10, E = 4, an off-axis Gaussian with its London phase; k = 8, thresh 1e−5, L = 20.
- Method: `test_magnetic_kernels mode=make/apply/compare` against Armen's `MagneticGreen` (built on
  `magnetic-fourier`) and against the pure-Python chain (`chain_ref.py` in
  `/huge/aiuser/flbiai/test/magnetic_kernel/`).
- Results:
  - ft: 6e−16 (bitwise);
  - chain_symmetric: 1.1e−12;
  - chain_landau: 1.2e−12.

**2. Kernels against each other.** `ctest -R mra_magnetic`, B = 2, thresh 1e−4, nbroaden 6.
- Test functions: an axial Gaussian, an m = −1 Gaussian, an off-axis Gaussian.
- ft vs chain_symmetric, Landau (conjugated) vs symmetric, and Mehler vs symmetric (m = 0): all agree
  to 1.7–6.2e−4 relative, below 10·thresh.
- The weak resolvent identity ⟨f|(T_B − E)Gf⟩ = ⟨f|f⟩ holds to 1e−6 to 3e−6.
- With `nbroaden` = 0 (the original setting), the pairwise differences level off at about 3e−5 as
  thresh decreases. That is the edge-term error of narrow kernels between boxes of different level;
  `nbroaden` = 6 removes it (Landau vs symmetric: 3e−5 → 4e−6 at thresh 1e−6).

**3. SCF.** The triplet-H₂ test is in `/huge/aiuser/flbiai/test/magnetic_kernel/`, with references from
the B ≤ 100 scan (`h2_h2plus_field_scan.md`) and new runs of the original codes.
- Agreement at thresh 1e−4 is at the thresh level (1e−5 to 1e−4 Eh), not bitwise. J, K, V and the SCF
  algorithm (KAIN and canonical orbitals vs. a plain fixed point) differ between znemo and the drivers,
  although the kernels are identical.
- Example: B = 10, 90°, R = 0.8. znemo chain_symmetric gives a kernel-consistent E_T of 6.506523; the
  scan gives 6.506437.

## References

**Article repository** (github.com/fbischoff/2026-magnetic-kernel):
- `notes/magnetic_greens_function.md`: the magnetic Green's function, the gauge phase, the
  decomposition of T_B.
- `notes/magnetic_kernel_realspace_splitting.md`: the exact three-factor product (chain).
- `notes/mra_magnetic_scheme.md`: the MRA scheme of the chain (sub-steps, LLL shortcut, gauges);
  `notes/decoupling_magnetic_phase.md`.
- `notes/mehler_kernel_factorization.md`: the rank-1 Mehler kernel and its use as an iteration kernel.
- `notes/greens_function_approaches.md`, `notes/greens_function_general_molecules.md`: comparison of
  the approaches, molecules off the field axis.
- `notes/resolvent_kernel_laplace.md`: resolvent, Green's function, heat kernel, Laplace transform.
- `h2_h2plus_field_scan.md`: H₂⁺ and H₂ field scan with the chain (the SCF references).

**Branch `magnetic-fourier`** (github.com/armen-chemist/madness):
- `magnetic/README.md`, `magnetic/docs/parameters_physics.pdf`, `src/madness/mra/magnetic_green.h`.
