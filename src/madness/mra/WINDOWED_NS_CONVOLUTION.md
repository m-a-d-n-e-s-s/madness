# Tail-Completed, Two-Threshold Windowed Non-Standard 3D Convolution in MADNESS

This document describes the mathematical formulation, data structures, parallel/distributed execution pipeline, `mTxmq` kernel lowering, verification suite, and benchmark performance of the **Tail-Completed, Two-Threshold Windowed Non-Standard (NS) 3D Convolution** implemented in `src/madness/mra/`.

---

## 1. Motivation & Design Overview

Applying isotropic 3D integral operators represented as sums of separated Gaussians—most notably the free-space Coulomb kernel $1/r$, the bound-state Helmholtz (BSH) kernel $e^{-\mu r}/(4\pi r)$, and the Yukawa kernel $e^{-\mu r}/r$—on deeply refined, unbalanced octrees presents two classical bottlenecks in standard multiwavelet NS form:

1. **Gibbs (`erf`) Interface Oscillations Across Coarse-to-Fine Jumps**: When a coarse leaf box abuts a finer subtree across a face, edge, or corner without 2:1 tree grading, extrapolating high-order Legendre polynomials *outside* the coarse leaf's support creates artificial boundary discontinuities. Convolving those discontinuities with sharp Gaussians induces spurious error-function (`erf`) oscillations that force unnecessary adaptive refinement in Molecular Orbitals (MOs), inflating leaf counts around nuclear cores and valence-core transitions.
2. **Deep-Tree Cusp Quadrature Cost**: Resolving an all-electron nuclear cusp of charge $Z$ to accuracy $\epsilon$ requires octree depths $n_{\max} \approx 20\text{--}26$ and $O(10^2)$ Gaussian exponents spanning $15+$ orders of magnitude. Sweeping every Gaussian across every tree level scales poorly with cusp depth and wastes work on levels where a Gaussian is either ultra-sharp ("spiky", acting as a local differential moment operator) or ultra-smooth ("diffuse", already resolved on coarser scales).

The Windowed NS algorithm resolves both bottlenecks simultaneously by enforcing **five invariants**:

- **Invariant 1 — Zero Exterior Polynomial Extrapolation**: Every virtual halo box $V$ at level $j$ is strictly a sub-box of an existing coarser real leaf $B \in \mathcal{L}$ ($\text{level}(B) < j$), prolonged *solely inside* $B$'s support one dyadic step at a time using exact two-scale Legendre filters $(h_0, h_1)$.
- **Invariant 2 — Two-Threshold Scale Separation & Constant Window Width $W$**: At tree level $j$, a Gaussian of simulation-cell exponent $\alpha_{\text{sim},p} = e^{2t_p}$ has dimensionless stiffness $\gamma_{p,j} = \alpha_{\text{sim},p} 4^{-j}$. It is active in 3D NS sum-factorization **only** inside a narrow scale window $j \in [j_0(p), j_1(p)]$ of constant width $W \approx 4\text{--}6$ levels bounded by an exact Plancherel diffuse threshold $\gamma_{\text{diff}}(k,\epsilon)$ and a closed-form spiky threshold $\gamma_{\text{spiky}}(k) = k^2(k+1)^2/30$.
- **Invariant 3 — Overflow-Free Gauss–Christoffel Spiky Compression & Discrete Tail Completion**: All Gaussians spiky at leaf level $n$ ($n_{\text{spiky}}(p) \ge n$, including the analytic $p > p_{\text{tail}}$ cusp tail) are compressed via log-rescaled Stieltjes–Golub–Welsch quadrature into $Q \in \{4, 5\}$ local Legendre moment operators $\bm{M}(\tilde{\gamma}_{n,q})^{\otimes 3}$ per leaf. At the opposite end ($\mu = 0$), the omitted $p < p_{\min}$ infrared tail is completed by matching exact discrete geometric moments $(M_0, M_1, M_2)$ with a single synthetic Gaussian plus a root $(0,0,0)$ constant shift $w_{\text{const}}$, eliminating Euler–Maclaurin tail bumps.
- **Invariant 4 — Exact Boundary-Halo Commutation (Theorem 3.1)**: When a Gaussian enters its NS window at base level $j_0(p)$, any coarser leaf $B$ ($\text{level}(B) < j_0(p)$) has *already* accumulated $p$'s local moment contribution inside Stage 2's Gauss–Christoffel leaf rule. Subtracting $W_p \bm{M}(\gamma_{p,j_0})^{\otimes 3} s_{\text{virt}}[V]$ on primary virtual halo boxes $V \in \partial\mathcal{H}_{j_0}^{(1,T)}$ before adjoint restriction cancels the interior polynomial contribution to machine precision while preserving exact smooth cross-interface flux.
- **Invariant 5 — Dyadic Scale-Locking ($\beta = 4^{1/M}$)**: Choosing quadrature spacing $\Delta t = \frac{\ln 2}{M}$ locks $\gamma_{p,j}$ to a universal table of $W \times M$ (`~15–35`) reference stiffnesses $\gamma_{\text{ref}}(w_{\text{eff}}, m)$ independent of tree depth $n_{\max}$, screening parameter $\mu$, and cubic cell width $L$.

---

## 2. Source Layout & Public Entry Points

### 2.1 File Map

| File | Role in Implementation |
| :--- | :--- |
| `windowed_quadrature.h` | Two-threshold solver (`compute_gamma_spiky`, `compute_gamma_diff`), parity-sparse `LegendreMomentOp` ($\bm{M}(\gamma)$ / $\bm{M}(\gamma)^{\mathsf{T}}$), dyadic-locked `WindowedQuadraturePlan`, discrete low/high tail completion, and log-scale Stieltjes–Golub–Welsch `GCLevelRule` compression |
| `windowed_ns_stencil.h` | 64-bit SplitMix64 `WindowedKeyHash` (`Level <= 62`), universal process-wide `(k, W, M, eps)` memoized `WindowedStencilCache1D`, `StencilEntry1D` transition blocks (`T`, `R`, `R_2k_to_k`, `R_k_to_2k`), contiguous `U_by_rank` 1D SVDs (`StencilSVD1D`), and exact 3D NS sub-block/shell/constant-wavelet norms |
| `windowed_ns_halo.h` | Idle-thread-aware `parallel_for_chunks`, 64-bit packed spatial key indexing (`j <= 20`, exact `Key<3>` fallback for `j > 20`), `WindowedNSHaloWorkspace` (serial & 3-fence MPI build), single-pass interior `virt_cache` prolongation, bottom-up adjoint restriction, and dynamic subtree promotion |
| `windowed_ns_sumfact.h` | Pruned dual-direction 3D sweep graphs (`build_sumfact_graph_3d`), Tier-1 `mTxmq` (`libxsmm` JIT) + L1 `fast_transpose` + Z-weight-folded `execute_sumfact_3d`, self-box Tucker contractions, 3-tier source-norm bucketing, `candidate_exceeds_tol` constant-subtracted screening, and 5-stage `apply_windowed_ns_3d<T>` |
| `windowed_convolution.h`, `windowed_convolution.cc` | `WindowedSeparatedConvolution<double,3>`, factory functions (`WindowedCoulombOperator[Ptr]`, `WindowedBSHOperator[Ptr][3D]`, `WindowedYukawaOperator[Ptr]`), thread-safe plan/stencil cache, LPT-scheduled batched vector `apply`, and opt-in `SeparatedConvolution` bridge |
| `operator.h`, `mra.h`, `vmra.h` | Non-intrusive opt-in bridge (`SeparatedConvolution::enable_windowed_ns(bool)` and `MADNESS_WINDOWED_NS=1`) for scalar and vector `apply(...)`, plus `OT_BSH` metadata population in `BSHOperatorPtr3D` |
| `mraimpl.h` | Standalone `project_refine_op` fix (`else if` in `newspecialpts` filtering) eliminating $O(2^{\text{special\_level}})$ special-point list duplication on deep nuclear cusps |
| `test_windowed_convolution.cc` | 4-Gate unit, unbalanced interface jump, 22–25 level nuclear cusp, heterogeneous batched `apply`, and multi-rank MPI parity suite (`--small`, `--full`, `--gate=1..4`, `--write-parity`, `--check-parity`) |

### 2.2 Usage & Runtime Opt-In

1. **Direct C++ API (`windowed_convolution.h`)**:
   ```cpp
   auto pop = WindowedCoulombOperatorPtr(world, eps);          // 1 / r
   auto bop = WindowedBSHOperatorPtr(world, mu, eps);          // exp(-mu*r) / (4*pi*r)
   auto yop = WindowedYukawaOperatorPtr(world, mu, eps);       // exp(-mu*r) / r

   real_function_3d V = apply(*pop, rho);                      // scalar real or complex
   vecfuncT Bpsi = apply(world, bop_ptrs, psi);                // LPT-batched vector apply
   ```
2. **Per-Operator Bridge on Legacy `SeparatedConvolution<double,3>`**:
   ```cpp
   SeparatedConvolution<double, 3> op = CoulombOperator(world, lo, eps);
   op.enable_windowed_ns(true);                                // delegates scalar & vector apply
   ```
3. **Global Environment Opt-In (Zero Application Code Changes)**:
   ```bash
   MADNESS_WINDOWED_NS=1 ./madqc --geometry=water
   ```
   Eligible 3D free-space cubic-cell `OT_G12` (Coulomb) and `OT_BSH` (BSH) operators automatically delegate both scalar `apply(op, f)` (`mra.h`) and batched `apply(world, op, vec_f)` / `apply(world, vec_ops, vec_f)` (`vmra.h`) to `WindowedSeparatedConvolution<double, 3>`, while periodic, range-restricted, modified, or non-3D operators remain on the legacy code path.

---

## 3. Detailed Algorithm as Implemented

### 3.1 Two-Threshold Scale Separation (`windowed_quadrature.h`)

Let $L > 0$ be the cubic simulation cell width (`FunctionDefaults<3>::get_cell_width()`), $\mu_{\text{sim}} = \mu_{\text{user}} L$, and $k$ be the Legendre scaling order on each octree box. In unit simulation coordinates $\bm{x} \in [0, 1]^3$, a 3D Gaussian convolution term $c_p e^{-\alpha_{\text{sim},p} \|\bm{x}-\bm{y}\|^2}$ with $\alpha_{\text{sim},p} = e^{2t_p}$ has physical 3D operator weight
$$W_p^{\text{scaled}} = \text{op\_scale} \cdot c_p \left(\frac{\pi}{\alpha_{\text{sim},p}}\right)^{3/2}, \qquad \text{op\_scale} = \begin{cases} L^2 / (4\pi), & \text{BSH}, \\ L^2, & \text{Coulomb / YukawaRaw}. \end{cases}$$

At octree level $j$ (box width $2^{-j}$), its dimensionless stiffness is $\gamma_{p,j} = \alpha_{\text{sim},p} 4^{-j} = e^{2t_p - 2j\ln 2}$. Two exact thresholds partition all scales:

1. **Spiky Threshold $\gamma_{\text{spiky}}(k)$ (`WindowedQuadraturePlan::compute_gamma_spiky`)**:
   $$\gamma_{\text{spiky}}(k) = \frac{k^2 (k + 1)^2}{30}.$$
   Whenever $\gamma_{p,n} \ge \gamma_{\text{spiky}}(k)$ on a leaf box at level $n$, the Gaussian's support is strictly narrower than the degree-$(k-1)$ Legendre grid resolution, and its action on smooth polynomial scaling coefficients reduces (up to $\lesssim \epsilon$) to the truncated asymptotic differential moment series $\bm{M}(\gamma_{p,n})^{\otimes 3}$.

2. **Diffuse Threshold $\gamma_{\text{diff}}(k, \epsilon)$ (`WindowedQuadraturePlan::compute_gamma_diff`)**:
   Defined as the unique root of the exact 3D spherical-Bessel / Plancherel wavelet projection bound:
   $$F_k(\gamma) = 2 \sup_{z > 0} \left[ e^{-z^2 / \gamma} \sqrt{\sum_{n=k}^{\infty} (2n + 1) j_n(z)^2} \right] = \epsilon.$$
   - The tail ratio $S_k(z) = \sum_{n=k}^{\infty} (2n+1) j_n(z)^2 / \sum_{n=0}^{\infty} (2n+1) j_n(z)^2$ (where the denominator equals $1$ by Parseval's identity) is evaluated in `spherical_bessel_tail_sum(z, k)` via **Miller's backward recurrence** starting at $n_{\text{top}} = \max(k+30, \lceil z + 30 + 5 z^{1/3}\rceil)$ with $10^{-100}$ dynamic overflow rescaling, avoiding catastrophic cancellation as $z \to 0^+$.
   - `plancherel_bound_value(gamma, k)` brackets the global maximum around $z_{\text{peak}} \approx \sqrt{k\gamma/2}$ on a 40-interval log-grid $[ \ln z_{\text{peak}} - 3, \ln z_{\text{peak}} + 2.5 ]$ followed by 60 golden-section refinement steps.
   - `compute_gamma_diff(k, eps)` seeds from the closed-form double-factorial asymptotic root and runs 64 bisection steps in $\ln \gamma$, memoizing `(k, eps) -> gamma_diff` in a static thread-safe cache.

For each quadrature term $p$, the transition levels are:
$$n_{\text{spiky}}(p) = \left\lfloor \frac{2t_p - \ln\gamma_{\text{spiky}}(k)}{2\ln 2} \right\rfloor, \qquad n_{\text{diff}}(p) = \left\lfloor \frac{2t_p - \ln\gamma_{\text{diff}}(k,\epsilon)}{2\ln 2} \right\rfloor,$$
giving the active NS window $j \in [j_0(p), j_1(p)]$ with $j_0(p) = \max(0, n_{\text{spiky}}(p) + 1)$ and $j_1(p) = \min(n_{\max} - 1, n_{\text{diff}}(p))$. For $j > n_{\text{diff}}(p)$, wavelet outputs vanish below $\epsilon$ and scaling outputs follow automatically from two-scale reconstruction of coarser levels.

---

### 3.2 Dyadic-Locked Quadrature, Tail Completion & GC Compression (`WindowedQuadraturePlan`)

1. **Dyadic Scale-Locking (`select_oversampling_M`)**:
   Using the exponential substitution $t = \ln s$ in $\frac{e^{-\mu r}}{r} = \frac{2}{\sqrt{\pi}} \int_{-\infty}^{\infty} \exp\!\big(-e^{2t} r^2 - \frac{\mu^2}{4} e^{-2t} + t\big)\,dt$ with step size
   $$\Delta t = \frac{\ln 2}{M}, \qquad \beta = e^{2\Delta t} = 4^{1/M},$$
   the smallest integer $M \in [2, 8]$ is chosen so Proposition 4.1's trapezoidal discretization error $2/\sqrt{\cosh(\pi^2 / \Delta t)} \le 0.1\epsilon$ (typically $M = 3\text{--}5$). Decomposing $p = q_{\text{oct}} M + m$ with Euclidean residue $m \in \{0, \dots, M-1\}$ (`divmod_euclidean`), we precompute exact integer base offsets:
   $$s_{\text{base,spiky}}[m] = \left\lfloor \frac{m}{M} - \frac{\ln\gamma_{\text{spiky}}}{2\ln 2} \right\rfloor, \qquad s_{\text{base,diff}}[m] = \left\lfloor \frac{m}{M} - \frac{\ln\gamma_{\text{diff}}}{2\ln 2} \right\rfloor.$$
   Then $n_{\text{spiky}}(p) = q_{\text{oct}} + s_{\text{base,spiky}}[m]$ and $n_{\text{diff}}(p) = q_{\text{oct}} + s_{\text{base,diff}}[m]$ hold with **zero floating-point floor drift** across all octaves $q_{\text{oct}} \in \mathbb{Z}$, and at every active level $j \in [j_0(p), j_1(p)]$:
   $$\gamma_{p,j} = \gamma_{\text{ref}}(w_{\text{eff}}, m) = \exp\!\Big(2\ln 2 \big(\tfrac{m}{M} - (s_{\text{base,spiky}}[m] + 1 + w_{\text{eff}})\big)\Big), \qquad w_{\text{eff}} = j - (n_{\text{spiky}}(p) + 1) \in [0, W - 1],$$
   where $W = \max_{0 \le m < M}(s_{\text{base,diff}}[m] - s_{\text{base,spiky}}[m])$.

2. **Discrete Low-Tail Completion at $\mu = 0$**:
   On $[0, 1]^3$ ($r_{\max} = \sqrt{3}$), truncating the $t_p = p\Delta t$ sum cleanly at $p_{\min} = \lfloor \frac{1}{7\Delta t} \ln(\frac{525\sqrt{\pi}\epsilon_{\text{quad}}}{4 r_{\max}^7}) \rfloor$ ($\epsilon_{\text{quad}} = 0.1\epsilon$) leaves an infinite geometric series over $p < p_{\min}$ whose exact discrete moments in $r^2$ are:
   $$M_r = \frac{2\Delta t}{\sqrt{\pi}} \frac{e^{(2r+1) t_{p_{\min}}}}{\operatorname{expm1}((2r+1)\Delta t)}, \qquad r \in \{0, 1, 2\}.$$
   Matching $M_0, M_1, M_2$ exactly via a single off-grid Gaussian $\alpha_{\text{low}} = M_2 / M_1$, $c_{\text{low}} = M_1^2 / M_2$ plus a constant kernel shift $w_{\text{const,sim}} = M_0 - c_{\text{low}}$ (`w_const_scaled = op_scale * w_const_sim`) reproduces $1/r$ across the entire domain $[10^{-8}, \sqrt{3}]$ to $O(\epsilon_{\text{quad}})$ without any continuous-vs-discrete Euler–Maclaurin mismatch. For $\mu_{\text{sim}} > 0$, super-exponential decay $\exp(-\frac{1}{4}\mu_{\text{sim}}^2 e^{-2t_p})$ allows direct truncation at $p_{\min}$ with $w_{\text{const}} = 0$.

3. **High-$\alpha$ Cusp Completion & Log-Scale Stieltjes–Golub–Welsch Compression (`build_gc_rule_for_level`)**:
   Explicit dyadic terms are generated up to $p_{\text{tail}} = p_{\max} + 10M$, where $p_{\max} + 1$ is the first grid index spiky at $n_{\max}$. At each tree level $n \ge 0$, all spiky terms ($n_{\text{spiky}}(p) \ge n$) have inverse stiffness $u_{p,n} = \gamma_{p,n}^{-1} = 4^n e^{-2t_p} \le \gamma_{\text{spiky}}(k)^{-1}$ and level-normalized positive weights $\widehat{W}_{p,n} = 4^n W_p^{\text{scaled}}$ (which remain $O(u_{p,n})$ and completely free of overflow/underflow at arbitrary $n$). When $\mu_{\text{sim}} = 0$ and $n_{\text{spiky}}(p_{\text{tail}}) \ge n$, the exact infinite geometric tail $p > p_{\text{tail}}$ is appended as a single remainder spike:
   $$\widehat{W}_{\text{rem},n} = \frac{\widehat{W}_{p_{\text{tail}},n}}{\beta - 1}, \qquad u_{\text{rem},n} = \frac{u_{p_{\text{tail}},n}}{\beta + 1}.$$
   Rescaling $s_i = u_i / u_{\max} \in (0, 1]$ and $w_i = \widehat{W}_i / \sum_l \widehat{W}_l$, discrete Stieltjes orthogonalization builds the $Q_{\text{eff}} \times Q_{\text{eff}}$ symmetric tridiagonal Jacobi matrix $\bm{J}$ ($Q = 4$ for $k \le 6$, $Q = 5$ for $7 \le k \le 14$), diagonalized via `madness::syev`. Eigenvalues $s_q$ and top eigenvector components $V_{0,q}^2$ yield $Q_{\text{eff}} \le 5$ compressed nodes $\tilde{\gamma}_{n,q} = (u_{\max} s_q)^{-1}$ and weights $\widetilde{W}_{n,q} = 4^{-n} (\sum_l \widehat{W}_l) V_{0,q}^2$, each storing its precomputed $k \times k$ lower-triangular matrix $\bm{M}(\tilde{\gamma}_{n,q})^{\mathsf{T}}$.

4. **Parity-Sparse 1D Local Moment Operator (`LegendreMomentOp`)**:
   On orthonormalized Legendre polynomials $\tilde{\phi}_i(x) = \sqrt{2i+1}\,P_i(2x-1)$ on $[0, 1]$, the exact second-derivative operator has entries
   $$D_2(i, q) = \frac{1}{2}\sqrt{(2i+1)(2q+1)}\,(q - i)(i + q + 1) \quad \text{for } q > i,\; (q - i) \text{ even},$$
   and zero otherwise. `LegendreMomentOp` precomputes $\frac{1}{r!}\bm{D}_2^r$ for $0 \le r \le \lfloor(k-1)/2\rfloor$ and forms
   $$\bm{M}(\gamma) = \sum_{r=0}^{\lfloor(k-1)/2\rfloor} \gamma^{-r} \left(\frac{1}{r!}\bm{D}_2^r\right),$$
   storing its transpose `MT(in_deg, out_deg)` so 3D action $\bm{M}(\gamma)^{\otimes 3} s$ maps directly to a single allocation-free `fast_transform(s, MT, tmp, work)` call.

---

### 3.3 Universal 1D Stencil Cache & Exact 3D Sub-Block Norms (`windowed_ns_stencil.h`)

`WindowedStencilCache1D` delegates to a process-wide thread-safe table (`get_or_build_universal_stencil_entries`) keyed solely by $(k, W, M, \epsilon, \gamma_{\text{ref}}(0,0), \gamma_{\text{ref}}(W-1,M-1))$. Consequently:
- Across all $40\text{--}80$ `BSHOperatorPtr3D(mu_i)` constructions per molecular SCF run (where orbital energies $\varepsilon_i$ and $\mu_i = \sqrt{-2\varepsilon_i}$ change every iteration), **zero** 1D quadrature builds and **zero** LAPACK `dgesvd` calls occur after the initial $W \times M$ table construction (`< 0.5 ms`).
- Each `StencilEntry1D` queries `GaussianConvolution1DCache<double>::get(k, gamma_ref, 0, LatticeRange())->nonstandard(0, disp)` for $\text{disp} = l_{\text{target}} - l_{\text{source}} \in [-2, +2]$ and stores:
  - Contiguous `(dim_in, dim_out)` matrices oriented directly for `mTxmq`: `T` $(k \times k)$, `R` $(2k \times 2k)$, `R_2k_to_k` $(2k \times k)$, and `R_k_to_2k` $(k \times 2k)$.
  - Exact 1D SVDs (`StencilSVD1D`) for all four block shapes, including **rank-major contiguous prefix matrices `U_by_rank[r]`** of shape `(d_in, r)` (`ldb == r` with zero column stride padding!) and exact Frobenius tail norms `tail_normf[r]` $= (\sum_{m \ge r} \sigma_m^2)^{1/2}$.
  - Exact 3D Frobenius norms via orthogonal block subtraction—e.g. for Pass 4b-1 (`(0, d) -> s`):
    $$\|\bm{R}_{2k \to k}(\bm{d})\big|_{(0,d)\to s}\|_F = \sqrt{\prod_{\alpha \in \{x,y,z\}} \|R_{2k\to k}(d_\alpha)\|_F^2 - \prod_{\alpha \in \{x,y,z\}} \|R_{ss}(d_\alpha)\|_F^2},$$
    analogously for `opnorm_4b2_leaf` (`(s, 0) -> d`), `opnorm_4b2_int` (`(s, d) -> d`), and `opnorm_4b2_wav` (`(0, d) -> d`), aggregated into per-shell maxima `shell_max_*[sh]` and shell $L_2$ norms `shell_l2_*[sh]` for $\text{sh} = \|\bm{d}\|_\infty \in \{0, 1, 2\}$.
  - Exact 3D constant-scaling-field wavelet residual norm `const_wav_3d_norm[m_bw]` $= \|\bm{Q}^j (\sum_{\|\bm{d}\|_\infty \le m_{\text{bw}}} \bm{R}(\bm{d})^{\otimes 3} \bm{e}_{000})\|_F$ for $m_{\text{bw}} \in \{0, 1, 2\}$.

---

### 3.4 Five-Stage Execution Pipeline (`windowed_ns_halo.h`, `windowed_ns_sumfact.h`)

Given an input `Function<T, 3>` ($T \in \{\text{double}, \text{double\_complex}\}$):

#### Stage 1 — Exact NS Analysis & Fine-Driven Interior Virtual Halo (`WindowedNSHaloWorkspace`)
1. **Stage 1a (`prepare_windowed_ns_input_3d`)**: Transitions the input function to `TreeState::nonstandard_with_leaves` (if `compressed`, `make_nonstandard(true)` populates internal $(2k)^3$ NS cubes while retaining leaf $k^3$ scaling cubes; `was_compressed` is saved so Stage 5b restores the input in zero FLOPs via `standard()`).
2. **Parallel Local Node Ingestion & Global Norm Reduction**:
   - Real nodes on each rank are partitioned by level $j \in [0, n_{\max}]$ (`WindowedOwnedRealNode<T>`), extracting read-only scaling cubes `s_in`, average coefficient `s000 = s_in(0,0,0)`, oscillatory norm squared `s_osc_norm2` $= \|s_{\text{in}}\|_F^2 - |s_{000}|^2$, and (for internal nodes) pure wavelet cubes `wav_coeff` (with octant $[0..k)^3$ zeroed) in parallel chunks.
   - Per-level maxima (`max_s_real`, `max_s_leaf`, `max_d_int`, `max_ns_int`, `has_internal`) are reduced globally across `world` (`gop.max`) to form rank-invariant `WindowedGlobalLevelStats`, including `max_coarse_leaf_above[j]` $= \max_{b < j} \text{max\_s\_leaf}[b]$ (which rigorously bounds any prolonged virtual box $\|s_{\text{virt}}\|_F$ by orthogonality of $(h_0, h_1)$).
3. **Active Bandwidth & Per-Node Shell-2 Threshold Classification (`radii_builder`)**:
   - Using `loc_tol_j = f_impl->truncate_tol(eff_tol, Key<3>(j))`, each candidate term in `raw_T_terms[j]` and `raw_gamma_terms[j]` selects active half-bandwidths `bw_T` $\in \{-1, 1, 2\}$, `bw_4b1` $\in \{-1, 0, 1, 2\}$, `bw_leaf` $\in \{-1, 1, 2\}$ (**never `0`**, preserving `[-1,+1]^3` neighbor polynomial cancellation), and `bw_int >= bw_leaf` $\in \{-1, 0, 1, 2\}$.
   - Shell-1 uses tight tolerance $0.01 \cdot \text{loc\_tol}_j / \sqrt{N_{\text{terms}}}$, while Shell-2 ($\|\bm{d}\|_\infty = 2$) uses Legacy MADNESS's threshold $\text{loc\_tol}_j / (14.14 \sqrt{N_{\text{terms}}})$ and computes per-node source-norm thresholds `sh2_s_thresh`, `sh2_4b1_thresh`, and `sh2_refine_thresh` so diffuse valence nodes only emit Shell-1 (`[-1,+1]^3`) halos even on levels where a nuclear cusp requires `m = 2`.
4. **Fine-Driven Halo Discovery & Single-Pass Interior Prolongation (`virt_cache`)**:
   - Real nodes emit neighbor candidates tagged with bitmask flags:
     - `kHaloH1T`: primary virtual neighbors around $\mathcal{T}_j$ needed for Stage 4(a) $\bm{T}^{j_0}$ + Theorem 3.1 subtraction.
     - `kHaloH1Gamma`: primary virtual neighbors around $\mathcal{T}_j^{\text{int}}$ receiving Pass 4b-1 scaling spillover `ds_corr`.
     - `kHaloH1Refine`: virtual neighbors around $\mathcal{T}_j^{\text{int}}$ eligible for dynamic wavelet output promotion in Pass 4b-2.
     - `kHaloH2Ghost`: secondary read-only virtual scaling sources completing stencils around targets.
   - Any neighbor key already present in $\mathcal{T}_j$ is NOT virtual; if it is a real leaf adjacent to an internal node, its `near_internal_4b1` / `near_internal_4b2` flags are set. Any neighbor key finer than $\mathcal{T}$ (descendant of an internal node at $\ge j$) has no coarse leaf ancestor and is skipped. Every remaining candidate strictly lies inside a unique coarser real leaf $B \in \mathcal{L}$ ($\text{level}(B) < j$).
   - **Serial (`world.size() == 1`)**: `build_serial` runs lock-free in parallel over 16-node slices using 64-bit packed coordinates `(lx << 42) | (ly << 21) | lz` at $j \le 20$.
   - **Distributed MPI (`world.size() > 1`)**: `build_distributed` executes **3 `WorldObject` collectives** via `WindowedNSHaloExchange<T>`: (1) route candidate `(Key<3>, flags)` to `pmap->owner(key)` to classify owned real leaves vs owned `H1` virtuals; (2) walk parent chains of needed virtual/real stencil sources and send deduplicated `KeyQueryMsg` batches; (3) reply with read-only real leaf `s_in` and real internal `ns_in` tensors.
   - **Single-Pass Top-Down Prolongation (`populate_prolonged_virtual_and_sources`)**: A single bottom-up pass ($j = n_{\max} \dots 2$) inserts intermediate virtual parents up to `coarse_leaf`; then a single parallel top-down pass ($j = 1 \dots n_{\max}$) prolongs each virtual box once via `general_fast_transform` using octant filters `h_oct_[child_octant][axis]`.

#### Stage 2 — Cross-Level Flattened Leaf Gauss–Christoffel Pass
All owned real leaves ($B \in \mathcal{L}$ with $\|s_B\|_F > 0$) across $j = 0 \dots n_{\max}$ are flattened into a single vector `gc_leaves` and processed via `parallel_for_chunks`:
$$\hat{s}_B \mathrel{+}= \sum_{q=1}^{Q_{\text{eff}}(n)} \widetilde{W}_{n,q}\, \bm{M}(\tilde{\gamma}_{n,q})^{\otimes 3} s_B, \qquad n = \text{level}(B),$$
using per-chunk L1-resident $k^3$ scratch cubes (`fast_transform` + `aligned_axpy`, zero heap allocations per leaf).

#### Stage 3 — Root ($j = 0$) Constant Low-Tail & Ultra-Diffuse Pass
On the rank owning root key `(0, {0,0,0})`:
1. If $w_{\text{const}}^{\text{scaled}} \ne 0$ ($\mu = 0$), adds $w_{\text{const}}^{\text{scaled}} \, s_{0,(0,0,0)}$ directly to $\hat{s}_{0,(0,0,0)}$ (exact since $P_0 \equiv 1$ is orthogonal to all higher modes and wavelets).
2. Applies $W_p^{\text{scaled}} (\bm{T}^0_{p,00})^{\otimes 3} s_0$ via `fast_transform` + `aligned_axpy` for every ultra-diffuse term with $n_{\text{diff}}(p) < 0$.

#### Stage 4(a) — Base Scaling Injection $\bm{T}^{j_0}$ & Theorem 3.1 Halo Subtraction
For each level $j \in [1, n_{\max}]$ with non-empty `active_T_terms[j]`:
1. **3-Tier Deterministic Source-Norm Bucketing**: Active scaling sources (`ld.sources`: real nodes + `H1`/`H2` virtuals) are partitioned by $\|s\|_F / M_s$ ($M_s = \text{max\_s\_source}[j]$) into up to 3 disjoint sorted buckets:
   - Bucket 0 (cusp core): $\|s\|_F > 0.02\,M_s$ (`bref = M_s`)
   - Bucket 1 (intermediate): $10^{-4} M_s < \|s\|_F \le 0.02\,M_s$ (`bref = 0.02 M_s`)
   - Bucket 2 (diffuse tail): $\|s\|_F \le 10^{-4} M_s$ (`bref = 10^{-4} M_s`)
   Each bucket independently tightens active bandwidth `bw_T_b` and 1D SVD ranks `SlotOp1D`, groups terms by effective bandwidth $m_{\text{bw}} \in \{1, 2\}$, builds one shared `SumFactGraph3D` per $(b_i, m_{\text{bw}})$, and accumulates $W_p^{\text{scaled}} \bm{T}^j_p s_{\text{src}}$ onto $\mathcal{T}_j \cup \partial\mathcal{H}_j^{(1,T)}$ (`s_out` on real nodes, `ds_corr` on `kHaloH1T` virtuals).
2. **Theorem 3.1 Boundary-Halo Commutation Subtraction**: Strictly *once* per `ActiveBaseTerm` (outside the source-norm bucket loop!), every owned virtual node $V \in \partial\mathcal{H}_j^{(1,T)}$ subtracts the interior local moment already accounted for on $V$'s coarser real leaf ancestor $B$ in Stage 2:
   $$\Delta\hat{s}_{\text{corr}}[V] \mathrel{-}= W_p^{\text{scaled}}\, \bm{M}(\gamma_{p,j})^{\otimes 3} s_{\text{virt}}[V].$$
   Where $f$ is locally polynomial inside $B$, $\sum_{\text{src}} \bm{T}^j_p s_{\text{src}} - \bm{M}(\gamma_{p,j})^{\otimes 3} s_{\text{virt}}[V] = 0$ identically; across genuine fine features adjacent to $B$, $\Delta\hat{s}_{\text{corr}}[V]$ captures the exact smooth tail spilling from $\mathcal{T}_j$ into $B$.

#### Stage 4(b) — Non-Standard Wavelet Decoupling $\bm{\Gamma}^j$ & Dynamic Output Refinement
Executed only at levels $j \in [0, n_{\max}-1]$ where `ld.global_has_internal == true`:
1. **Pass 4b-1 (`2k -> k`, `(0, d^j) -> \hat{s}^j`, Cancellation-Free)**:
   - Targets: internal real nodes (`rn.is_internal`), boundary real leaves (`rn.near_internal_4b1`), and `kHaloH1Gamma` virtual nodes (`vn.ds_corr`).
   - Sources: pure wavelet cubes `se.wav_coeff` of internal nodes (`se.is_internal`), partitioned into the 3 deterministic norm buckets by `d_normf`.
   - Self-box terms (`eff_bw == 0`) run via `execute_self_box_2k_to_k_3d`; neighbor terms (`eff_bw in {1, 2}`) run via `execute_sumfact_3d` with `R_2k_to_k` / `svd_R_2k_to_k`.
2. **Constant-Subtracted Refinement Candidate Pre-Screening (`candidate_exceeds_tol`)**:
   - Internal real nodes always participate in Pass 4b-2. Boundary real leaves (`!rn.is_internal && rn.near_internal_4b2`) and virtual `kHaloH1Refine` boxes only need `(2k)^3` wavelet evaluation if their output wavelet norm can exceed $\tau_{\text{loc}}(U) = \text{truncate\_tol}(\text{eff\_thresh}, U) \cdot \text{dynamic\_refine\_factor}$.
   - `candidate_exceeds_tol(U, c0, tol_loc)` separates the nearest-neighbor shell $\text{sh} \le 1$ (`[-1,+1]^3`, 27 boxes) from outer Shell 2 ($\text{sh} = 2$). Whenever all 27 `sh <= 1` slots lie in-domain and exist in `ld.sources` (`full_sh1_stencil == true`—**never** invalidated by pruned `mLeaf_h2 == 1` ghosts at `sh == 2`!), constant polynomial cancellation around $c_0 = s_{000}(U)$ gives the tight bound:
     $$\text{ub}_{\text{sh1}} = |c_0|\,\eta_{\text{const,sh1}} + \sum_{\|\bm{d}\|_\infty \le 1} \Big[ K_{\text{leaf,sh1}}(\bm{d}) \sqrt{\|s_{\text{osc}}\|_2^2 + |s_{000} - c_0|^2} + \mathbb{I}_{\text{int}}\big(K_{\text{outer,sh1}}(\bm{d})\|s\|_F + K_{\text{wav}}(\bm{d})\|d\|_F\big) \Big],$$
     plus unsubtracted Shell-2 bounds $\text{ub}_{\text{sh2}}$ (using $\min(\text{max\_coarse\_leaf\_above}, \text{sh2\_s\_thresh})$ for any omitted `mLeaf_h2 == 1` ghost slot). Smooth valence jump candidates are rejected (`>90%` hit rate) before allocating `(2k)^3` output cubes.
3. **Pass 4b-2 (`2k -> 2k` Internal + `k -> 2k` Leaf/Virtual -> $\hat{d}^j$)**:
   - Bucketed internal sources `se.ns_coeff` contract via `R` / `svd_R` (`bw_int in {0, 1, 2}`); bucketed leaf & virtual sources `se.s_coeff` contract via `R_k_to_2k` / `svd_R_k_to_2k` (`bw_leaf in {1, 2}`, never `0`).
   - Afterward, `zero_s0_octant` clears `[0..k)^3` on `ns_out` / `ns_virt_out` (leaving the 7 wavelet octants $\hat{d}^j$), and any boundary real leaf or `kHaloH1Refine` virtual node with $\|\hat{d}^j\|_F > \tau_{\text{loc}}$ sets `promoted_leaf = true` / `promoted = true`.
4. **Exact Shared-Tree Linearity Mode (`dynamic_refine == false`)**:
   - Replaces data-dependent level norms with `1.0` in bandwidth selection, disables per-source norm bucketing/screening, and skips output tree promotion so the discrete operator is an exact linear map ($\|\mathcal{G}(\alpha f + \beta g) - \alpha\mathcal{G}(f) - \beta\mathcal{G}(g)\|_2 = 0$ to machine precision on shared octrees).

#### Stage 5 — Dynamic Subtree Promotion, Bottom-Up Adjoint Restriction & Reconstruction
1. **Stage 5a (`assemble_output_nonstandard_tree`)**:
   - **Fast Unpromoted Path**: When no virtual node or real leaf is promoted (`restrict_all_virtual_to_coarse_leaves`), virtual `ds_corr` cubes are restricted bottom-up ($j = n_{\max} \dots 1$) one dyadic step at a time using sorted flat-vector parent slots (`ParentRestrictSlot`) and exact adjoint filters `hT_oct_[child_octant][axis]` (`general_fast_transform` + `aligned_axpy`), then accumulated onto owning real coarse leaves (sorting MPI `CoarseLeafDeltaMsg` packets by `(leaf_key, min_source_key)` for deterministic rank-independent summation).
   - **General Dynamic Promotion Path (`restrict_and_split_promoted_virtual`)**: Virtual items are grouped by `coarse_leaf_key`. Unpromoted coarse-leaf groups restrict in parallel across `ThreadPool` chunks; promoted coarse-leaf groups mark all ancestors from `bkey` to every promoted `virt_key` as `split_internals`, instantiate their direct unsplit children as `split_leaves`, restrict any deeper unpromoted `ds_corr` entries onto `split_leaves`, emit $(2k)^3$ NS internal nodes and $k^3$ leaf nodes into `result`'s `WorldContainer`, and set `result` to `TreeState::nonstandard_after_apply`.
2. **Stage 5b (`restore_windowed_ns_input_3d` + `result.reconstruct`)**:
   - Restores input `f` (`clear` if `destructive()`, `standard()` if originally `compressed`, or `reconstruct()`), and calls `result.reconstruct(fence)` to sum two-scale parent scaling + wavelet octants down to standard reconstructed leaf scaling coefficients.

---

## 4. High-Performance Kernel & Parallel Execution Optimizations

### 4.1 100% Tier-1 `mTxmq` (`libxsmm` JIT `s_dense_table`) Lowering & Z-Sweep Weight Folding

Every 3D tensor contraction in Stages 1–5 maps to MADNESS optimized tensor kernels with zero scalar C++ tensor loops:

| Stage / Operation | Kernel Lowering |
| :--- | :--- |
| Stage 1b Prolongation & Stage 5a Restriction | `general_fast_transform` (3 cyclic `mTxmq` passes) + `aligned_axpy` |
| Stage 2 Leaf GC (`Q <= 5`) & Stage 4(a) `-M(gamma)` | `fast_transform` (3 cyclic `mTxmq` passes) + `aligned_axpy` |
| Stage 3 & 4(a) Root (`j = 0`) | `fast_transform` + `aligned_axpy` |
| Stage 4(b) Self-Box (`m_bw == 0`, `2k -> k` & `2k -> 2k`) | 6-step Tier-1 `mTxmq` Tucker contraction `(2k)^3 -> r^3 -> d_out^3` using dense-packed `U_by_rank[r]` (`ldb == r`), or `general_fast_transform` / `fast_transform` when high-rank |
| Stage 4(a) & 4(b) Spatial 3D Sweeps (`m_bw in {1, 2}`) | `execute_sumfact_3d` -> `apply_slot_edge_1d`: Tier-1 `mTxmq` + L1 `fast_transpose` + `aligned_add<T, T>` |

Inside `execute_sumfact_3d` (`windowed_ns_sumfact.h`), two low-level properties are critical for throughput:
1. **Guaranteed Tier-1 `s_dense_table` Dispatch**: MADNESS's `mTxmq(dimi, dimj, dimk, c, a, b)` dispatches to its pre-JIT'd `libxsmm` `s_dense_table` kernel only when `dimi <= 400`, `dimj <= 24`, and `dimk <= 24`. In `apply_slot_edge_1d`, each cyclic sweep contracts `(d_in, dimi)` with `U_by_rank[r]` of contiguous shape `(d_in, r)` (`dimi = d_in^2 <= 256` for $k \le 8$, `dimj = r <= 16 <= 24`, `dimk = d_in <= 16 <= 24`):
   - Step 1: `mTxmq(dimi, rank, d_in, tmp_svd, src, u_ptr)` -> `[dimi, rank]` (100% Tier-1 JIT hit).
   - Step 2: If `rank == 1`, `[dimi, 1]` and `[1, dimi]` have identical memory layout, so the transpose is **elided completely**; for `rank >= 2`, L1-resident `fast_transpose(dimi, rank, tmp_svd, tmp_svd_T)` produces `[rank, dimi]`.
   - Step 3: `mTxmq(dimi, d_out, rank, target, tmp_svd_T, svt_ptr)` -> `[dimi, d_out]` (100% Tier-1 JIT hit), accumulated via `aligned_add<T, T>`.
2. **Z-Sweep Weight Folding**: By multilinearity, $\alpha (A_x \otimes A_y \otimes A_z) = A_x \otimes A_y \otimes (\alpha A_z)$. Scaling a stack-allocated 64B-aligned copy of the $\le 5$ active Z-sweep 1D factors (`<= 1,280` doubles total) by `weight` once per `execute_sumfact_3d` call allows Sweep 3 (Z) to accumulate directly into `tgt_out_ptrs[ti]`, eliminating per-term $(2k)^3 N_{\text{tgt}}$ temporary Z-buffer zeroing and `aligned_axpy(4096)` memory traffic.

### 4.2 Dual-Direction Pruned 3D Sweep Graphs & 64-Bit Packed Spatial Keys (`j <= 20`)

- **64-Bit Packed Spatial Coordinates (`windowed_ns_halo.h`)**: At any fixed tree level $j \le 20$ ($0 \le l_x, l_y, l_z < 2^{20} < 2^{21}$), spatial coordinates are packed into a single `uint64_t`:
  $$\text{code}(l_x, l_y, l_z) = (\text{uint64\_t}(l_x) \ll 42) \;\big|\; (\text{uint64\_t}(l_y) \ll 21) \;\big|\; \text{uint64\_t}(l_z).$$
  Unsigned `<` on `code` induces the **exact same total order** as `windowed_key_less` at fixed level $j$, shrinking keys $4\times$ (`8 B` vs `32 B` `Key<3>`) and turning 4-branch comparisons and `Key<3>` hashing/construction inside `build_serial`, `build_sumfact_graph_3d_packed`, `intersect_sorted_keys`, and `candidate_exceeds_tol` into 1-cycle integer operations (with automatic fallback to full `Key<3>` at $j > 20$).
- **Dual-Direction Graph Pruning (`build_sumfact_graph_3d`)**: Automatically selects forward candidate generation ($\mathcal{D}_{\text{src}} \to \mathcal{D}_x \to \mathcal{D}_{xy} \to \mathcal{D}_{\text{tgt}}$ with backward dead-pencil pruning) when $|\mathcal{D}_{\text{src}}| \le |\mathcal{D}_{\text{tgt}}|$ (e.g. a small cusp bucket feeding a larger target set) and backward candidate generation ($\mathcal{D}_{\text{tgt}} \to \mathcal{D}_{xy} \to \mathcal{D}_x$) otherwise.

### 4.3 LPT Batched Vector `apply`, `thread_purge()` Wakeup & Idle-Thread Donation

In `apply_windowed_batch_impl` (`windowed_convolution.cc`):
1. Input handles are deduplicated by `FunctionImpl*`, reconstructed/transformed to `nonstandard_with_leaves` in a single batched Stage 1a pass, and restored/reconstructed in a single batched Stage 5b pass.
2. On single-rank (`target_world.size() == 1`, `N >= 2`), indices `0..N-1` are sorted by **descending tree node count (Longest-Processing-Time-First, LPT)** so deep nuclear-cusp orbitals (e.g. `O 1s`) start immediately at $t = 0$.
3. High-priority `WindowedBatchWorkerTask` helpers are pushed to `ThreadPool::add` followed immediately by `thread_purge()` (flushing MADNESS's thread-local `prebufhi[]` when `MADNESS_DQ_USE_PREBUF=1` so all worker threads wake at $t = 0$).
4. As valence orbitals finish and decrement `g_active_windowed_batch_workers`, `parallel_for_chunks` (`kSumFactSweepGrain = 2`, `kDefaultParallelGrain = 4`) dynamically recruits the newly idle pool threads into the remaining cusp orbital's X/Y/Z sum-factorization sweeps without oversubscribing cores.

---

## 5. Verification Results (`test_windowed_convolution.cc` & Legacy MRA Suites)

`test_windowed_convolution` implements an adversarial 4-Gate verification matrix covering unit algebra, unbalanced 3D coarse-to-fine jumps, 22–25 level nuclear cusps, heterogeneous batched vector `apply`, and multi-rank MPI parity:

### 5.1 Gate-by-Gate Numerical Results

1. **Gate 1 — Commutation, 3D GC Spectral Norm, Tail Quadrature & Deep Keys (`15 / 15 PASS`)**:
   - Scale-locked stiffness drift $\max |\gamma_{\text{ref}}(w_{\text{eff}}, m) - \gamma_{p,j}|/\gamma_{\text{ref}} = \mathbf{6.22 \times 10^{-15}} \le 10^{-13}$.
   - `LegendreMomentOp` transpose identity $\|\bm{M}^{\mathsf{T}} - (\bm{M})^{\mathsf{T}}\|_F = \mathbf{0.00}$.
   - Theorem 3.1 exact boundary-halo commutation ($k \in \{6, 8, 10, 12\}$, level jumps $\Delta j \in \{1, 5, 15\}$): 1D Frobenius error $\mathbf{4.83 \times 10^{-16}}$, 3D tensor relative error $\mathbf{1.54 \times 10^{-14}}$, 3D partition-of-unity round-trip $\mathbf{5.44 \times 10^{-16}}$.
   - 3D Gauss–Christoffel spectral 2-norm error up to $n = 30$: $\mathbf{5.28 \times 10^{-14}}$ (Coulomb) and $\mathbf{8.40 \times 10^{-18}}$ (BSH).
   - Pointwise radial kernel error ($\epsilon = 10^{-9}$, $r \in [10^{-8}, \sqrt{3}]$): $\mathbf{3.18 \times 10^{-11}}$ (Coulomb, $L = 1, 8$, zero outer-tail bump) and $\mathbf{2.31 \times 10^{-12}}$ (BSH, $\mu = 0.5, 2.0$).
   - `WindowedKeyHash` & `WorldContainer` collision-free up to **`Level = 35`** (`665` extreme/midplane keys).

2. **Gate 2 — Unbalanced 3D Coarse-to-Fine Interface Jumps (`42 / 42 PASS`)**:
   - Frozen coarse octants (`level = 5`) directly abutting fine subtrees (`level = 8, 11, 15`, $\Delta j = 3, 6, 10$) across faces, edges, and corners on $L = 1$ and $L = 8$:
     - $k = 8$: Global $L_2$ error **`5.86e-10 – 6.12e-10`**, face/edge/corner pointwise error **`3.10e-9 – 3.74e-9`**, 41-point micro-grid 2nd-difference curvature residual **`8.08e-11 – 5.94e-10`** (zero `erf` interface ringing).
     - $k = 10$: Global $L_2$ error **`2.19e-10 – 2.29e-10`**, face/edge/corner pointwise error **`2.19e-10 – 5.83e-10`**.

3. **Gate 3 — 22–25 Level All-Electron Slater Cusps, Dynamic Midplane Promotion & Batched Linearity (`49 / 49 PASS`)**:
   - All-electron Slater density $\rho(\bm{r}) = \frac{Z^3}{\pi} e^{-2Z|\bm{r}-\bm{R}_A|}$ at $\bm{R}_A = (0.5000004, 0.5, 0.5)$ ($\epsilon = 10^{-8}$, analytic potential $V(\bm{R}_A) = Z$ and Hartree energy $E_H = \frac{5}{16} Z$):

     | Case | Octree Depth | Output $V$ Nodes | Global Rel $L_2$ Error | Cusp $|V(\bm{R}_A)-Z|/Z$ | Hartree $|E_H - \frac{5}{16}Z|/(\frac{5}{16}Z)$ |
     | :--- | :---: | :---: | :---: | :---: | :---: |
     | `Z=20, k=8, L=4` | `23` | `12,681` | `2.34e-09` | `2.87e-09` | `8.09e-11` |
     | `Z=20, k=10, L=4` | `22` | `4,873` | `2.68e-09` | `2.55e-10` | `4.22e-11` |
     | `Z=50, k=8, L=4` | `25` | `25,769` | `2.59e-09` | `2.91e-09` | `1.17e-10` |
     | `Z=50, k=10, L=4` | `24` | `7,273` | `2.66e-09` | `2.56e-10` | `3.27e-11` |
     | `Z=50, k=10, L=1` | `23` | `6,601` | `3.62e-10` | `2.57e-10` | `5.81e-11` |

   - Dynamic output promotion across $x = 0.5$ ($\alpha = 1.2 \times 10^{11}$ at $(0.50003, 0.5, 0.5)$ with $x < 0.5$ forced coarse at `depth = 3` in input): promoted output to **`depth = 16`** (`1.87e-10` near-midplane error, `4.73e-09` global $L_2$ error).
   - `special_level = 23` cusp projection regression test (`project_refine_op` fix in `mraimpl.h`): **`0.26 s`**.
   - Shared-tree `Function<double_complex, 3>` linearity: **`0.000000e+00 <= 1e-13`**.
   - Heterogeneous 6-element LPT batched vector `apply` (`alpha = 36.0` sharp Gaussian + 4 valence Gaussians + duplicate handle) across single-op, multi-op, and `enable_windowed_ns(true)` bridges vs scalar `apply`: **`0.000000e+00 <= 1e-13`**.

4. **Gate 4 — Multi-Rank MPI Parity, Zero-Leak Audit & Legacy Regressions (`ALL PASS`)**:
   - `np = 1 --full --write-parity` (`MAD_NUM_THREADS=4`): **`112 / 112 PASSED, 0 FAILED`** (`51.9 s`, `0` leaked functions / `0` leaked `WorldObject`s, `34` observables recorded).
   - `np = 2 --full --check-parity` (`MAD_NUM_THREADS=4`): **`113 / 113 PASSED, 0 FAILED`** (`128.7 s`, `34 / 34` observables match `np = 1` with max relative/absolute difference **`2.84e-16 <= 1e-13`**, `0` leaks).
   - `np = 4 --full --check-parity` (`MAD_NUM_THREADS=1`): **`113 / 113 PASSED, 0 FAILED`** (`131.7 s`, `34 / 34` observables match `np = 1` with max difference **`1.63e-16 <= 1e-13`**, `0` leaks).
   - Legacy MRA regression suites (`testbsh`, `testopdir`, `testsuite --small`, `test_memory_measurement`) all **PASS**, and both GCC 15.2 (`-DMADNESS_WERROR=ON`) and Clang 21.1 (`-Werror -Wall -Wextra`) compile with **`0` warnings**.

---

## 6. Performance & Timing Results

All benchmarks were executed on Linux under identical compiler flags (`-O2`), identical binary, and identical thread configuration (`np = 1`, `MAD_NUM_THREADS=4`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`) with `protocol = [1e-4, 1e-6]` (`k = 6` in Stage 1, `k = 8` in Stage 2) and `lo = 1e-10`.

### 6.1 All-Electron Hartree–Fock `madqc --geometry=water` (`H2O`, 5 MOs, `L = 51.43` Bohr)

| Metric / Phase | Default Legacy (`NS=0`) | Windowed Initial (`NS=1`) | **Windowed Optimized (`NS=1`)** | **Optimized vs Default / vs Initial** |
| :--- | :---: | :---: | :---: | :--- |
| **Total Wall Time** | `37.1 s` | `152.5 s` | **`27.5 s`** | **`1.35x` faster vs Default (`-25.9%`)** / **`5.55x` vs Initial** |
| **Total SCF Energy (`Ha`)** | `-76.067934304` | `-76.067934793` | **`-76.067934802`** | `\|ΔE\| = 4.98e-7 Ha` vs Default (`9e-9` vs Initial) |
| **SCF Iterations (`k=6` + `k=8`)** | `8 + 2 = 10` | `8 + 2 = 10` | **`8 + 2 = 10`** | Identical convergence |
| **Final MO Leaves (`O 1s` core)** | `5,784` (`1,100`) | `4,632` (`1,000`) | **`4,632` (`1,000`)** | **`-19.9%` fewer total MO leaves** |
| `guess Coulomb potn` | `0.31 s` | `0.94 s` | **`0.32 s`** | Parity vs Default / **`2.94x` vs Initial** |
| `Coulomb` (S1 + S2 = Total) | `1.01 + 1.18 = 2.19 s` | `3.50 + 5.18 = 8.68 s` | **`1.44 + 1.50 = 2.94 s`** | `0.74x` vs Default / **`2.95x` vs Initial** |
| **`HF exchange` (S1 + S2 = Total)** | `11.89 + 10.78 = 22.67 s` | `40.42 + 57.18 = 97.60 s` | **`6.02 + 7.21 = 13.23 s`** | **`1.71x` faster vs Default** / **`7.38x` vs Initial** |
| **`Apply BSH` (S1 + S2 = Total)** | `2.90 + 2.69 = 5.59 s` | `17.34 + 21.75 = 39.09 s` | **`2.46 + 2.40 = 4.86 s`** | **`1.15x` faster vs Default** / **`8.04x` vs Initial** |
| Peak RSS | `1.12 GiB` | `1.11 GiB` | **`1.53 GiB`** | Within `1.6 GiB` |

#### Optimization Pass Progression on `madqc --geometry=water`
- **Initial reference implementation**: `152.5 s` (`HF exchange 97.60 s`, `Apply BSH 39.09 s`, `Coulomb 8.68 s`)
- **Pass 1** (`fast_transform` GC, exact 3D sub-block norms, flat sorted vectors, batched `apply`): `99.8 s` (`1.53x`)
- **Pass 2** (universal `(k, W, M, eps)` stencil cache, `U_by_rank` + `mTxmq` + `aligned_add`, LPT scheduling): `65.4 s` (`2.33x`)
- **Pass 3** (Tier-1 `s_dense_table` `mTxmq` + L1 `fast_transpose`, 3-tier norm bucketing, Shell-2 pruning): `54.1 s` (`2.82x`)
- **Pass 4** (`sh <= 1` constant-subtraction fix in `candidate_exceeds_tol`, Z-sweep weight folding, 64-bit packed spatial keys, `kSumFactSweepGrain = 2`, and `thread_purge()` batch pool wakeup): **`27.5 s` (`5.55x` vs Initial, `1.35x` faster than Default `37.1 s`)**

### 6.2 All-Electron Hartree–Fock `madqc --geometry=water2` (`(H2O)2` Dimer, 10 MOs)

| Metric / Phase | Default Legacy (`NS=0`) | Windowed Initial (`NS=1`) | **Windowed Optimized (`NS=1`)** | **Optimized vs Default / vs Initial** |
| :--- | :---: | :---: | :---: | :--- |
| **Total Wall Time** | `96.4 s` | `496.5 s` | **`74.6 s`** | **`1.29x` faster vs Default (`-22.6%`)** / **`6.66x` vs Initial** |
| **Total SCF Energy (`Ha`)** | `-152.141501178` | `-152.141502456` | **`-152.141502495`** | `\|ΔE\| = 1.32e-6 Ha` vs Default (`3.9e-8` vs Initial) |
| **SCF Iterations (`k=6` + `k=8`)** | `9 + 2 = 11` | `8 + 2 = 10` | **`8 + 2 = 10`** | **`-1` Stage-1 SCF iter vs Default** |
| **Final MO Leaves (`O 1s` cores)** | `14,008` (`3,392`) | `11,752` (`2,128`) | **`11,752` (`2,128`)** | **`-16.1%` total leaves** (**`-37.3%` `O 1s` core leaves**) |
| `guess Coulomb potn` | `0.45 s` | `1.41 s` | **`0.59 s`** | `0.76x` vs Default / **`2.39x` vs Initial** |
| `Coulomb` (S1 + S2 = Total) | `1.45 + 2.26 = 3.71 s` | `8.32 + 6.19 = 14.51 s` | **`2.61 + 3.04 = 5.65 s`** | `0.66x` vs Default / **`2.57x` vs Initial** |
| **`HF exchange` (S1 + S2 = Total)** | `34.12 + 30.79 = 64.91 s` | `196.49 + 156.40 = 352.89 s` | **`18.21 + 25.06 = 43.27 s`** | **`1.50x` faster vs Default** / **`8.16x` vs Initial** |
| **`Apply BSH` (S1 + S2 = Total)** | `7.26 + 5.99 = 13.25 s` | `88.19 + 25.92 = 114.11 s` | **`5.89 + 6.12 = 12.01 s`** | **`1.10x` faster vs Default** / **`9.50x` vs Initial** |
| Peak RSS | `1.96 GiB` | `1.92 GiB` | **`2.37 GiB`** | Within `2.4 GiB` |

### 6.3 Gate 4c Deep 22-Level `Z = 30` Slater Cusp Benchmark (`depth = 22, eps = 1e-7, rho_nodes = 7,625`)

| Operator & MPI Configuration | Legacy MADNESS | Windowed Initial | **Windowed Optimized** | **Speedup vs Legacy / vs Initial & Tree Compression** |
| :--- | :---: | :---: | :---: | :--- |
| `Coulomb` (`np = 1, MAD_NUM_THREADS=4`) | `1.780 s` (`18,313` nodes) | `3.287 s` (`7,689` nodes) | **`0.828 s` (`7,689` nodes)** | **`2.15x` faster vs Legacy** / **`3.97x` vs Initial** (**`-58.0%` nodes**) |
| `BSH` (`mu = 1.0, np = 1, MAD_NUM_THREADS=4`) | `1.180 s` (`14,409` nodes) | `3.294 s` (`7,625` nodes) | **`0.572 s` (`7,625` nodes)** | **`2.06x` faster vs Legacy** / **`5.76x` vs Initial** (**`-47.1%` nodes**) |
| `Coulomb` (`np = 2, MAD_NUM_THREADS=4`) | `1.753 s` | — | **`1.262 s`** | **`1.39x` faster vs Legacy** |
| `BSH` (`mu = 1.0, np = 2, MAD_NUM_THREADS=4`) | `1.065 s` | — | **`0.883 s`** | **`1.21x` faster vs Legacy** |
| `Coulomb` (`np = 4, MAD_NUM_THREADS=1`) | `2.460 s` | — | **`1.558 s`** | **`1.58x` faster vs Legacy** |
| `BSH` (`mu = 1.0, np = 4, MAD_NUM_THREADS=1`) | `1.656 s` | — | **`1.116 s`** | **`1.48x` faster vs Legacy** |
