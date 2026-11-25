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
| `chain_symmetric` | T_B | symmetric | in the kernel | trapezoid in ln τ | `MagneticGreensFunction3D` in `hydrogen_magnetic_pymadness.py` (article repo, pymadness branch `magnetic`) |
| `chain_landau` | T_B | Landau, A_L = (−Bỹ, 0, 0) | in the kernel | trapezoid in ln τ | `Gop` of `h2p_landau.py --gauge landau` |
| `ft` | T_B | symmetric, O = 0 | in the kernel | trapezoid in ln τ | `MagneticGreen` in `magnetic_green.h`, branch `magnetic-fourier` (A. Armen) |
| `mehler` | −½∇² + B²ρ̃²/8 | symmetric | **right-hand side** | trapezoid in ln τ | new in 3D; 2D `MehlerGreensFunction2D` in `landau_mehler_pymadness.py` |

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
- Beyond τ = −ln(eps)/|B|, only the free z motion is propagated: the excited Landau levels are damped
  below eps.

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
- The cost scales with the number of q modes, about qmax P/π. `ft_qmax` = max(|B| L/2 + `ft_qextra`, Q + |B x_A|/2),
  L = min(cell half width, r_ext): Q from the nuclear cusps (`ft_qextra` < 0, the default). Cutting |q| > Q drops the cusp momenta of Vψ and raises the
  energy by (4/3) S/Q³, S = Σ_A n_A Z_A² ρ_A(0) (ρ_A(0) = Z_A³ g(|B|/Z_A²)/π, n_A = min(2, Z_A)); Q is chosen for a
  third of the ft precision floor 79 µEh S^0.75 (thresh/1e−4)^0.69 (fit, H atom, B = 1–100), plus |B x_A|/2.
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
| `trapezoid_log` | step h from `eps`, t_min = `lo`²/4 (`MagneticKernelParameters::set_derived_values`, derived as the Coulomb and BSH fits: relative error (4π/√h) e^{−π²/h} = eps for every spectral component; t_min·f shifts the energy by ≈ −10 % of the nuclear smoothing error of radius lo); tmax = −ln(eps)/(|B|/2 − E). znemo: lo = smallest smoothing radius, eps = 0.001 thresh | all |
| `mehler_m`, `mehler_log_c`, `splitting_coefficients` | stable for Bτ → 0 and Bτ → ∞ | all |

### Short proper times

None of the kernels has a dedicated short-time limit. The only short-time treatment is the cut of the
Laplace integral at t_min:

```
∫₀^{t_min} dτ e^{Eτ} e^{−τT} f  ≈  t_min · f        (heat kernel → δ, e^{Eτ} → 1)
```

| | chain | ft | mehler |
|---|---|---|---|
| t_min | lo²/4 | lo²/4 | lo²/4 |
| t_min · f in the code | `quad.tmin * a`, added unphased after the outer phase | `Psi = s0 * F` in every q mode | `quad.tmin * f` |
| first node | at t_min, weight ½ (trapezoid) | at t_min, weight ½ | at t_min, weight ½ |
| how the first node is reached | one sub-step from τ = 0, with the full phase pair e^{∓iBx~y~} | independent node | independent node |
| short-time errors | inherited by all later nodes (incremental chain) | local to the node | local to the node |

- **Error of t_min · f.** A component with T − E = κ gets the weight t_min instead of (1 − e^{−κ t_min})/κ,
  a relative error of about κ t_min/2. It matters only for κ ≳ 1/t_min, i.e. the cusp of Vψ; the bias
  scales with Z²ρ(R_A) (`magnetic_gf_consolidated.md` §7).
- **No τ → 0 branch.** The node formulas are evaluated in closed form at every τ (tanh, sinh, log-form
  normalizations, no clamp on narrow Gaussians). The special cases in `splitting_coefficients`
  (|Bτ| < 1e−8), `mehler_log_c` (wτ = 0) and `tanh_over_x` (x = 0) are B → 0 guards: with τ ≥ t_min they
  never fire at nonzero field. The `log_w < −60` skips and the 1e−10 exponent floor act on long times only.
- **Not the free kernel.** At short τ the symmetric-gauge kernel is the free Gaussian times the phase
  e^{i(B/2)(xy′−yx′)}, of size (B/2)ρ√τ (about 2 rad at B = 425, ρ = 1, τ = 1e−4). A short-time piece
  must keep this phase (Peierls form) or be local to a nucleus.
- **Cost.** In the chain the narrow short-time kernels act on phase-multiplied functions, which inflates
  the trees over the whole support. This dominates the chain cost (86 % of a G at τ < 1e−3 for a
  displaced H atom at B = 10, thresh 1e−5). The Python prototypes had cheaper short-time variants
  (Strang and exact Mehler steps for definite m, local Mehler, free short piece `--short free`); none is
  in the C++ code (`magnetic_gf_consolidated.md` §5, §10).
- **Choosing t_min.** The notes use t_min = 1e−4 with 40 nodes at thresh 1e−4 and t_min = 1e−5 with
  50 nodes at thresh 1e−5 (`notes/mra_magnetic_scheme.md` §2.3); the C++ defaults do not adapt to the
  threshold. Lowering t_min alone spreads the same `kernel_nodes` over a wider ln τ range, so raise the
  node count with it. A much smaller t_min makes the chain trees and memory explode
  (`notes/hydrogen_solvers_handoff.md`).

### Regularized right-hand side (nemo block, `ncf`)

With a nuclear correlation factor R (the `ncf` of the nemo block, as in nemo: default `(slater, 2.0)`, `(none, 1.0)` for the
bare orbitals) the kernels act on the nemos F = ψ/R. Since R is real,

```
R⁻¹ (π²/2) R F = (π²/2) F + U₁·∇F + i U₁·A F + (U₂ − V) F,   U₁ = −∇R/R
```

so the kernel G_B(E) is unchanged and its right-hand side is U₂F + U₁·∇F + iU₁·A F + (J − K)F (`Nuclear(ncf)`
and `zeeman_R_comm`, the latter in the gauge of the kernel), with the metric R² in the Fock matrix and the energy. The kinetic matrix is symmetric, π(RF) = R(π + iU₁)F:
T_ij = ½ Σ_d ⟨R² w_id|w_jd⟩, w_d = (π_d + iU₁_d)F (`kinetic_pi_matrix`); the Fock matrix and the energy then take the
bare smoothed potential V and no R commutator.

The bare right-hand side Vψ carries the nuclear cusp, |ĝ(p)|² ∝ Z²ρ(0)/p⁴; it sets tmin = lo²/4 (lo the smoothing
radius), the ft window, (4/3)S/Q³, and the ft precision floor. With R the cusp is gone: the leading non-smooth term is
U₁·∇F ≈ Z r̂·d (d = ∇F at the nucleus; d = 0 for s orbitals), |ĝ|² ∝ Z²|d|²/p⁶, so

- ft window: the cut costs (32/45) Z²|d|²/Q⁵; with |d|² ≤ 2n_Aρ_A(0), Q = [64 S/(45·0.1 thresh)]^{1/5};
- tmin: ⟨g|T − E|g⟩ = C is finite and the energy shift is −(tmin²/2) C: tmin = √(2·0.001 thresh/C), with C from the
  right-hand sides of every iteration (`adapt_kernel_tmin`: tmin grows in the first, shrinks later when C grows),
  no longer tied to the smoothing radius. (0.1 and 0.01 thresh left the kernel-consistent energy up to 25× further
  off than predicted; with 0.001 thresh H at B = 100 matches bare accuracy, E_k +4.6 vs +2.2 µEh, at 1/5 the time.)

Single applications ⟨g|G g⟩ on model right-hand sides of the H atom (thresh 1e−5; `test_magnetic_kernels mode=cusp`
in the experimental worktree): ft cut and floor −1150 → −25 (p) / +0.3 (s) µEh at qextra 0, B = 1; floor 9 / 21 / 100
µEh (bare) → below 0.3 µEh (p, s at B ≤ 10); chain at tmin 1e−3: +82 (bare) → +1.2 (p) / +3.5 (s) µEh at B = 1.

## Using the kernels

### C++

```cpp
#include <madness/mra/magnetic_kernels.h>
FunctionDefaults<3>::set_cubic_cell(-L, L);      // cubic cell; k and thresh as usual
MagneticKernelParameters p;                      // or p(world, parser): the magnetic_kernel block
p.set_derived_values(10.0, FunctionDefaults<3>::get_thresh());   // B; + nuclei, cusp_free
auto K = make_magnetic_kernel(world, "chain_landau", 10.0, p);   // sets FunctionDefaults<1>, <2> to match
complex_function_3d g = K->apply(f, E);          // (T_k - E)^{-1} f, E < K->bottom()
```

### znemo: the `complex` block

`kernel bsh`, the default, leaves znemo unchanged. With a magnetic kernel:

- **Kinetic terms.** The kinetic, paramagnetic and diamagnetic terms are in the Green's function (for
  `mehler`, the paramagnetic term stays on the right-hand side).
- **Energy and Fock matrix.** They use ½⟨π_dφ_i|π_dφ_j⟩, with π = −i∇ + A in the kernel's gauge.
- **Spin Zeeman.** The spin Zeeman term ±B/2 moves into the energy argument of G.
- **NCF, no box.** The nuclear correlation factor is the `ncf` of the nemo block (see "Regularized right-hand
  side"); no confining box is used.
- **Kernel energies.** As in the original drivers, the kernel energies E_i are started from the Fock
  matrix and updated by E − Re⟨Vψ|ψ − ψ̃⟩/‖ψ̃‖². An energy bound by at least 10⁻³
  is used as it is; one closer to |B|/2 is replaced by |B|/2 − 0.05 (internal; the output reports it).
  - Guess orbitals often have Fock energies above |B|/2.
  - Near the threshold the range of G(E) along B grows as 1/√(2(|B|/2 − E)), and so does the memory.
    With E = B/2 − 1e−3 at B = 10, a single application exceeded 45 GB.
- **Output.** `calc_info.json` gets `kernel` and `kernel_energy`, which is Σ_i E_i − (J − K) + E_nuc,
  the energy the original drivers reported.

| Parameter | Default | Meaning |
|---|---|---|
| `kernel` | `bsh` | `bsh`, `ft`, `chain_symmetric`, `chain_landau`, `mehler` |
| `guess_lz` | `[]` | initial orbitals: the restart sequence of moldft and nemo (dft `restart`: `auto` reads `<prefix>.restartdata` of a previous znemo calculation at the same geometry, eprec and B, and does not iterate if it is converged; `iterate`, `read_only`, `none` as in moldft), otherwise the hcore guess: diagonalize the one-electron Hamiltonian of the solver (kernel gauge, cut potentials) in magnetically translated, Landau-squeezed AOs of `aobasis` plus, in a field, the s and p₋ = (x − i sgn B · y) Slater functions e^{−√2 Z r} of each nucleus; `guess_lz` picks the occupied orbitals by L_z (B along a linear molecule; e.g. `[0, -1]` for σ_g π_u), empty: aufbau |

### The `magnetic_kernel` block

`MagneticKernelParameters`; the kind and the field are those of the `complex` block (`kernel`, `physical_B`),
passed to `make_magnetic_kernel` and `set_derived_values`. Derived parameters are set in
`MagneticKernelParameters::set_derived_values` (formulas there) and may be overridden like any other parameter.

| Parameter | Default | Meaning |
|---|---|---|
| `origin` | `[0, 0, 0]` | gauge origin (x, y; ft: 0) |
| `lo` | derived: min_A c_A | smallest length scale |
| `eps` | derived: 10⁻³ thresh | precision of (T_k − E)⁻¹ |
| `tau_step` | derived: (4π/√h) e^{−π²/h} = eps | step h in ln τ |
| `tmin` | derived: lo²/4; nemo: adapted, √(2 eps/C) | proper-time cut-off |
| `t_lll` | derived: −ln(eps)/\|B\| | chain: lowest Landau level only beyond |
| `max_btau`, `phase_tolfac`, `phase_recursive`, `nbroaden` | 2, 0.1, true, 0 | chain sub-steps, phase refinement, broadening |
| `ft_qextra` | −1 | ≥ 0: \|q\| ≤ \|B\| L/2 + qextra; < 0: from the nuclear cusps |
| `ft_qmax` | derived (see "Fourier approach") | y-wavenumber window |
| `ft_rcut`, `ft_rtaper` | derived: 1.5 √(8 ln(1/thresh)/\|B\|), min(0.5, rcut/4) | phase cylinder and taper |
| `ft_period` | derived: min(2 × cell width, 4 r_ext) | y-period |
| `ft_dphi`, `ft_qbatch` | 3, 64 | phase refinement, batch size |
| `ft_fq_tol` | derived: 0.01 thresh | dropped modes |

The smoothing radius c of the nuclear potential comes from `eprec` in the geometry block and the field:
c = [min(1e−3, eprec/2)/(0.65 Z⁵ g(|B|/Z²))]^{1/3}, g(β) = (1 + 0.2128β)(1 + 0.874 ln(1 + 0.3394β)).
The truncation threshold is the last value of the dft block's `protocol`, independent of `econv`.
The spin Zeeman term ±B/2 per α/β electron is not part of the orbital energies or the iteration; it is
added to the total energy only (printed as `spin zeeman term`). The print level is the dft block's `print_level` (3 adds timings and the kernel's own output).

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
  protocol [1e-4]
end
complex
  physical_B 10
  kernel chain_landau
end
geometry
  units atomic
  eprec 1e-6
  H  0.4 0.0 0.0
  H -0.4 0.0 0.0
end
```

Map from the original command-line options:

| Original | znemo |
|---|---|
| `h2p_madness.py --nodes --tmin --nlll --phase-tolfac --c --gap` | `tau_step`, `tmin`, `t_lll` (derived from thresh and eprec), `phase_tolfac`, `eprec` (geometry), — (internal) |
| `h2p_madness.py --guess2 sigma/pi` | `guess_lz` (π: `[0, -1]`) |
| `magneticatom h= tmin= tmax_factor= qextra= rcut= rtaper= dphi=` | `tau_step`, `tmin` (derived), — (derived), `ft_qextra`, `ft_rcut`, `ft_rtaper`, `ft_dphi` (magnetic_kernel block) |
| `magneticatom eprec=` | `eprec` (geometry); znemo adds the field dependence of c |

### Python

```python
import pymadness as mad
K = mad.make_magnetic_kernel(world, "ft", 10.0)   # optional: params=mad.MagneticKernelParameters(lo=...), nuclei=
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
- `notes/magnetic_kernels_common_origin.md`: ft, chain_symmetric and chain_landau derived from one
  Mehler factorization; why the z Gaussian is merged into the x convolution.
- `notes/mra_magnetic_scheme.md`: the MRA scheme of the chain (sub-steps, LLL shortcut, gauges);
  `notes/decoupling_magnetic_phase.md`.
- `notes/mehler_kernel_factorization.md`: the rank-1 Mehler kernel and its use as an iteration kernel.
- `notes/greens_function_approaches.md`, `notes/greens_function_general_molecules.md`: comparison of
  the approaches, molecules off the field axis.
- `notes/resolvent_kernel_laplace.md`: resolvent, Green's function, heat kernel, Laplace transform.
- `h2_h2plus_field_scan.md`: H₂⁺ and H₂ field scan with the chain (the SCF references).

**Branch `magnetic-fourier`** (github.com/armen-chemist/madness):
- `magnetic/README.md`, `magnetic/docs/parameters_physics.pdf`, `src/madness/mra/magnetic_green.h`.
