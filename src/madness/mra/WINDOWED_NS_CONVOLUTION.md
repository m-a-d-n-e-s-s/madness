# Tail-Completed, Two-Threshold Windowed Non-Standard 3D Convolution in MADNESS

This document describes the mathematical formulation, data structures, parallel/distributed execution pipeline, `mTxmq` kernel lowering, verification suite, and benchmark performance of the **Tail-Completed, Two-Threshold Windowed Non-Standard (NS) 3D Convolution** implemented in `src/madness/mra/`.

---

## 1. Motivation & Design Overview

Applying isotropic 3D integral operators represented as sums of separated Gaussians—most notably the free-space Coulomb kernel $1/r$, the bound-state Helmholtz (BSH) kernel $e^{-\mu r}/(4\pi r)$, the Yukawa kernel $e^{-\mu r}/r$, and general multi-scale Gaussian sums—on deeply refined, unbalanced octrees presents two classical bottlenecks in standard multiwavelet NS form:

1. **Gibbs (`erf`) Interface Oscillations Across Coarse-to-Fine Jumps**: When a coarse leaf box abuts a finer subtree across a face, edge, or corner without 2:1 tree grading, extrapolating high-order Legendre polynomials *outside* the coarse leaf's support creates artificial boundary discontinuities. Convolving those discontinuities with sharp Gaussians induces spurious error-function (`erf`) oscillations that force unnecessary adaptive refinement in Molecular Orbitals (MOs) and lose up to 5–6 digits of precision at intermediate-to-large exponents.
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
| `windowed_quadrature.h` | Two-threshold solver (`compute_gamma_spiky`, `compute_gamma_diff`), parity-sparse `LegendreMomentOp` ($\bm{M}(\gamma)$ / $\bm{M}(\gamma)^{\mathsf{T}}$), dyadic-locked `WindowedQuadraturePlan` (`Coulomb`, `BSH`, `YukawaRaw`, plus overflow-safe off-grid `CustomGaussian`), discrete low/high tail completion, and log-scale Stieltjes–Golub–Welsch `GCLevelRule` compression |
| `windowed_ns_stencil.h` | 64-bit SplitMix64 `WindowedKeyHash` (`Level <= 62`), universal process-wide `(k, W, M, eps)` memoized `WindowedStencilCache1D` (plus exact per-term stencils for `CustomGaussian`), `StencilEntry1D` transition blocks (`T`, `R`, `R_2k_to_k`, `R_k_to_2k`), contiguous `U_by_rank` 1D SVDs (`StencilSVD1D`), and exact 3D NS sub-block/shell/constant-wavelet norms |
| `windowed_ns_halo.h` | TBB (`tbb_arena`) & MADNESS `ThreadPool` parallel dispatch (`parallel_for_chunks`), 64-bit packed spatial key indexing (`j <= 20`, exact `Key<3>` fallback for `j > 20`), zero-mutation serial `TreeState::reconstructed` / `TreeState::compressed` ingestion (`WindowedSerialInputRange`), ungraded (`dj >= 2`) leaf-boundary gating (`has_ungraded_child_toward`), `WindowedNSHaloWorkspace` (serial & 3-fence MPI build), single-pass interior `virt_cache` prolongation, bottom-up adjoint restriction, direct single-rank `assemble_output_reconstructed_serial`, and dynamic subtree promotion |
| `windowed_ns_sumfact.h` | Smooth loose-tolerance tightening (`tighten_loose_windowed_eps`), thread-safe pooled 64B-aligned uninitialized scratch (`UninitAlignedBuf`, `SumFactScratchLease`), pruned dual-direction 3D sweep graphs (`build_sumfact_graph_3d`), Tier-1 `mTxmq` (`libxsmm` JIT) + L1 `fast_transpose` + Z-weight-folded `execute_sumfact_3d`, self-box Tucker contractions, 3-tier source-norm bucketing with $O(N)$ bandwidth-equivalent merging (`merge_matching_norm_buckets`), `candidate_exceeds_tol` constant-subtracted screening, and 5-stage `apply_windowed_ns_3d<T>` |
| `windowed_convolution.h`, `windowed_convolution.cc` | `WindowedSeparatedConvolution<double,3>` (including custom Gaussian sum constructors and `apply_pure_spiky_same_world` pure-GC fast path), factory functions (`WindowedCoulombOperator[Ptr]`, `WindowedBSHOperator[Ptr][3D]`, `WindowedYukawaOperator[Ptr]`), thread-safe plan/stencil cache, TBB/ThreadPool LPT-scheduled batched vector `apply`, and opt-in `SeparatedConvolution` bridge |
| `operator.h`, `mra.h`, `vmra.h` | Non-intrusive opt-in bridge (`SeparatedConvolution::enable_windowed_ns(bool)` and `MADNESS_WINDOWED_NS=1`) for scalar and vector `apply(...)` across `OT_G12`, `OT_BSH`, and isotropic 3D `CustomGaussian` operators, plus `OT_BSH` metadata population in `BSHOperatorPtr3D` |
| `mraimpl.h` | Standalone `project_refine_op` fix (`else if` in `newspecialpts` filtering) eliminating $O(2^{\text{special\_level}})$ special-point list duplication on deep nuclear cusps |
| `test_windowed_convolution.cc` | 4-Gate unit, unbalanced interface jump, 22–25 level nuclear cusp, sparse `vtransform` compressed-tree parity, heterogeneous batched `apply`, and multi-rank MPI parity suite (`--small`, `--full`, `--gate=1..4`, `--write-parity`, `--check-parity`) |
| `src/examples/large_exp_3d.cc` | Standalone 3D Gaussian convolution & octree inspection benchmark sweeping 41 kernel scales ($b \in [0.5, 5.50 \times 10^{11}]$) across diffuse, intermediate windowed, and ultra-spiky `pure-GC` regimes |

### 2.2 Usage & Runtime Opt-In

1. **Direct C++ API (`windowed_convolution.h`)**:
   ```cpp
   auto pop = WindowedCoulombOperatorPtr(world, eps);          // 1 / r
   auto bop = WindowedBSHOperatorPtr(world, mu, eps);          // exp(-mu*r) / (4*pi*r)
   auto yop = WindowedYukawaOperatorPtr(world, mu, eps);       // exp(-mu*r) / r
   WindowedSeparatedConvolution<double, 3> gop(world, coeff, expnt, eps); // sum_p coeff[p]*exp(-expnt[p]*r^2)

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
   Eligible 3D free-space cubic-cell `OT_G12` (Coulomb), `OT_BSH` (BSH), and arbitrary isotropic 3D Gaussian-sum `SeparatedConvolution<double, 3>` operators automatically delegate both scalar `apply(op, f)` (`mra.h`) and batched `apply(world, op, vec_f)` / `apply(world, vec_ops, vec_f)` (`vmra.h`) to `WindowedSeparatedConvolution<double, 3>`, while periodic, range-restricted, modified, or non-3D operators remain on the legacy code path.

---

## 3. Detailed Algorithm as Implemented

### 3.1 Two-Threshold Scale Separation & Loose-Tolerance Tightening (`windowed_quadrature.h`, `windowed_ns_sumfact.h`)

Let $L > 0$ be the cubic simulation cell width (`FunctionDefaults<3>::get_cell_width()`), $\mu_{\text{sim}} = \mu_{\text{user}} L$, and $k$ be the Legendre scaling order on each octree box. In unit simulation coordinates $\bm{x} \in [0, 1]^3$, a 3D Gaussian convolution term $c_p e^{-\alpha_{\text{sim},p} \|\bm{x}-\bm{y}\|^2}$ with $\alpha_{\text{sim},p} = e^{2t_p}$ has physical 3D operator weight
$$W_p^{\text{scaled}} = \begin{cases} \text{op\_scale} \cdot c_p \left(\dfrac{\pi}{\alpha_{\text{sim},p}}\right)^{3/2}, & \text{Coulomb / YukawaRaw } (\text{op\_scale} = L^2) \text{ or BSH } (\text{op\_scale} = L^2 / (4\pi)), \\[6pt] c_{\text{user},p} \left(\dfrac{\pi}{\alpha_{\text{sim},p}}\right)^{3/2}, & \text{CustomGaussian } (\alpha_{\text{sim},p} = L^2 \alpha_{\text{user},p}). \end{cases}$$
Every term also stores `log_weight_scaled` $= \ln|W_p^{\text{scaled}}|$ directly in log-space ($\ln|c_{\text{user},p}| + \frac{3}{2}\ln\pi - 3 t_p$) so ultra-spiky custom Gaussians ($\alpha_{\text{sim},p} \gtrsim 10^{15}$) never underflow IEEE-754 `double` before level-rescaled GC compression.

At octree level $j$ (box width $2^{-j}$), dimensionless stiffness is $\gamma_{p,j} = \alpha_{\text{sim},p} 4^{-j} = e^{2t_p - 2j\ln 2}$. Two exact thresholds partition all scales:

1. **Spiky Threshold $\gamma_{\text{spiky}}(k)$ (`WindowedQuadraturePlan::compute_gamma_spiky`)**:
   $$\gamma_{\text{spiky}}(k) = \frac{k^2 (k + 1)^2}{30}.$$
   Whenever $\gamma_{p,n} \ge \gamma_{\text{spiky}}(k)$ on a leaf box at level $n$, the Gaussian's support is strictly narrower than the degree-$(k-1)$ Legendre grid resolution, and its action on smooth polynomial scaling coefficients reduces (up to $\lesssim \epsilon$) to the truncated asymptotic differential moment series $\bm{M}(\gamma_{p,n})^{\otimes 3}$.

2. **Diffuse Threshold $\gamma_{\text{diff}}(k, \epsilon)$ (`WindowedQuadraturePlan::compute_gamma_diff`)**:
   Defined as the unique root of the exact 3D spherical-Bessel / Plancherel wavelet projection bound:
   $$F_k(\gamma) = 2 \sup_{z > 0} \left[ e^{-z^2 / \gamma} \sqrt{\sum_{n=k}^{\infty} (2n + 1) j_n(z)^2} \right] = \epsilon.$$
   - The tail ratio $S_k(z) = \sum_{n=k}^{\infty} (2n+1) j_n(z)^2 / \sum_{n=0}^{\infty} (2n+1) j_n(z)^2$ (where the denominator equals $1$ by Parseval's identity) is evaluated in `spherical_bessel_tail_sum(z, k)` via **Miller's backward recurrence** starting at $n_{\text{top}} = \max(k+30, \lceil z + 30 + 5 z^{1/3}\rceil)$ with $10^{-100}$ dynamic overflow rescaling, avoiding catastrophic cancellation as $z \to 0^+$.
   - `plancherel_bound_value(gamma, k)` brackets the global maximum around $z_{\text{peak}} \approx \sqrt{k\gamma/2}$ on a 40-interval log-grid $[ \ln z_{\text{peak}} - 3, \ln z_{\text{peak}} + 2.5 ]$ followed by 60 golden-section refinement steps.
   - `compute_gamma_diff(k, eps)` seeds from the closed-form double-factorial asymptotic root and runs 64 bisection steps in $\ln \gamma$, memoizing `(k, eps) -> gamma_diff` in a static thread-safe cache.

3. **Loose-Threshold Screening Tightening (`tighten_loose_windowed_eps`)**:
   In loose multi-center SCF protocols (`thresh >= 1e-4` with `truncate_mode == 1` on wide boxes $L \gtrsim 100$ Bohr), `truncate_tol` does not scale down on coarse octaves. Whenever `eps > 1e-5` (and no explicit `eps_loc` override is set), `tighten_loose_windowed_eps(eps)` tightens internal quadrature/stencil/screening tolerance smoothly to $\min(5\times 10^{-5},\, \max(10^{-5},\, 0.1\epsilon))$, leaving `eps <= 1e-5` (`1e-6`, `1e-8`, ...) completely untouched.

For each quadrature term $p$, the transition levels are:
$$n_{\text{spiky}}(p) = \left\lfloor \frac{2t_p - \ln\gamma_{\text{spiky}}(k)}{2\ln 2} \right\rfloor, \qquad n_{\text{diff}}(p) = \left\lfloor \frac{2t_p - \ln\gamma_{\text{diff}}(k,\epsilon)}{2\ln 2} \right\rfloor,$$
giving the active NS window $j \in [j_0(p), j_1(p)]$ with $j_0(p) = \max(0, n_{\text{spiky}}(p) + 1)$ and $j_1(p) = \min(n_{\max} - 1, n_{\text{diff}}(p))$. For $j > n_{\text{diff}}(p)$, wavelet outputs vanish below $\epsilon$ and scaling outputs follow automatically from two-scale reconstruction of coarser levels.

---

### 3.2 Dyadic-Locked Quadrature, Tail Completion & Overflow-Free GC Compression (`WindowedQuadraturePlan`)

1. **Dyadic Scale-Locking (`select_oversampling_M`)**:
   Using the exponential substitution $t = \ln s$ in $\frac{e^{-\mu r}}{r} = \frac{2}{\sqrt{\pi}} \int_{-\infty}^{\infty} \exp\!\big(-e^{2t} r^2 - \frac{\mu^2}{4} e^{-2t} + t\big)\,dt$ with step size
   $$\Delta t = \frac{\ln 2}{M}, \qquad \beta = e^{2\Delta t} = 4^{1/M},$$
   the smallest integer $M \in [2, 8]$ is chosen so Proposition 4.1's trapezoidal discretization error $2/\sqrt{\cosh(\pi^2 / \Delta t)} \le 0.1\epsilon$ (typically $M = 3\text{--}5$). Decomposing $p = q_{\text{oct}} M + m$ with Euclidean residue $m \in \{0, \dots, M-1\}$ (`divmod_euclidean`), we precompute exact integer base offsets:
   $$s_{\text{base,spiky}}[m] = \left\lfloor \frac{m}{M} - \frac{\ln\gamma_{\text{spiky}}}{2\ln 2} \right\rfloor, \qquad s_{\text{base,diff}}[m] = \left\lfloor \frac{m}{M} - \frac{\ln\gamma_{\text{diff}}}{2\ln 2} \right\rfloor.$$
   Then $n_{\text{spiky}}(p) = q_{\text{oct}} + s_{\text{base,spiky}}[m]$ and $n_{\text{diff}}(p) = q_{\text{oct}} + s_{\text{base,diff}}[m]$ hold with **zero floating-point floor drift** across all octaves $q_{\text{oct}} \in \mathbb{Z}$, and at every active level $j \in [j_0(p), j_1(p)]$:
   $$\gamma_{p,j} = \gamma_{\text{ref}}(w_{\text{eff}}, m) = \exp\!\Big(2\ln 2 \big(\tfrac{m}{M} - (s_{\text{base,spiky}}[m] + 1 + w_{\text{eff}})\big)\Big), \qquad w_{\text{eff}} = j - (n_{\text{spiky}}(p) + 1) \in [0, W - 1],$$
   where $W = \max_{0 \le m < M}(s_{\text{base,diff}}[m] - s_{\text{base,spiky}}[m])$. For arbitrary off-grid `CustomGaussian` operators, each user term stores its exact active-window stiffnesses in `ref_grid_gamma_` while sharing the identical two-threshold and GC machinery.

2. **Discrete Low-Tail Completion at $\mu = 0$**:
   On $[0, 1]^3$ ($r_{\max} = \sqrt{3}$), truncating the $t_p = p\Delta t$ sum cleanly at $p_{\min} = \lfloor \frac{1}{7\Delta t} \ln(\frac{525\sqrt{\pi}\epsilon_{\text{quad}}}{4 r_{\max}^7}) \rfloor$ ($\epsilon_{\text{quad}} = 0.1\epsilon$) leaves an infinite geometric series over $p < p_{\min}$ whose exact discrete moments in $r^2$ are:
   $$M_r = \frac{2\Delta t}{\sqrt{\pi}} \frac{e^{(2r+1) t_{p_{\min}}}}{\operatorname{expm1}((2r+1)\Delta t)}, \qquad r \in \{0, 1, 2\}.$$
   Matching $M_0, M_1, M_2$ exactly via a single off-grid Gaussian $\alpha_{\text{low}} = M_2 / M_1$, $c_{\text{low}} = M_1^2 / M_2$ plus a constant kernel shift $w_{\text{const,sim}} = M_0 - c_{\text{low}}$ (`w_const_scaled = op_scale * w_const_sim`) reproduces $1/r$ across the entire domain $[10^{-8}, \sqrt{3}]$ to $O(\epsilon_{\text{quad}})$ without any continuous-vs-discrete Euler–Maclaurin mismatch. For $\mu_{\text{sim}} > 0$, super-exponential decay $\exp(-\frac{1}{4}\mu_{\text{sim}}^2 e^{-2t_p})$ allows direct truncation at $p_{\min}$ with $w_{\text{const}} = 0$.

3. **High-$\alpha$ Cusp Completion & Log-Scale Stieltjes–Golub–Welsch Compression (`build_gc_rule_for_level`)**:
   Explicit dyadic terms are generated up to $p_{\text{tail}} = p_{\max} + 10M$, where $p_{\max} + 1$ is the first grid index spiky at $n_{\max}$. At each tree level $n \ge 0$, all spiky terms ($n_{\text{spiky}}(p) \ge n$) have inverse stiffness $u_{p,n} = \gamma_{p,n}^{-1} = 4^n e^{-2t_p} \le \gamma_{\text{spiky}}(k)^{-1}$ and level-normalized positive weights $\widehat{W}_{p,n} = 4^n W_p^{\text{scaled}}$ (or log-shifted `exp(log_weight_scaled + 2*n*ln2 - max_log_wh)` for arbitrary `CustomGaussian` sums, remaining completely free of overflow/underflow at arbitrary $n$ and $\alpha$). When $\mu_{\text{sim}} = 0$ and $n_{\text{spiky}}(p_{\text{tail}}) \ge n$, the exact infinite geometric tail $p > p_{\text{tail}}$ is appended as a single remainder spike:
   $$\widehat{W}_{\text{rem},n} = \frac{\widehat{W}_{p_{\text{tail}},n}}{\beta - 1}, \qquad u_{\text{rem},n} = \frac{u_{p_{\text{tail}},n}}{\beta + 1}.$$
   Rescaling $s_i = u_i / u_{\max} \in (0, 1]$ and $w_i = \widehat{W}_i / \sum_l \widehat{W}_l$, discrete Stieltjes orthogonalization builds the $Q_{\text{eff}} \times Q_{\text{eff}}$ symmetric tridiagonal Jacobi matrix $\bm{J}$ ($Q = 4$ for $k \le 6$, $Q = 5$ for $7 \le k \le 14$), diagonalized via `madness::syev`. Eigenvalues $s_q$ and top eigenvector components $V_{0,q}^2$ yield $Q_{\text{eff}} \le 5$ compressed nodes $\tilde{\gamma}_{n,q} = (u_{\max} s_q)^{-1}$ and weights $\widetilde{W}_{n,q}$, each storing its precomputed $k \times k$ lower-triangular matrix $\bm{M}(\tilde{\gamma}_{n,q})^{\mathsf{T}}$.

4. **Parity-Sparse 1D Local Moment Operator (`LegendreMomentOp`)**:
   On orthonormalized Legendre polynomials $\tilde{\phi}_i(x) = \sqrt{2i+1}\,P_i(2x-1)$ on $[0, 1]$, the exact second-derivative operator has entries
   $$D_2(i, q) = \frac{1}{2}\sqrt{(2i+1)(2q+1)}\,(q - i)(i + q + 1) \quad \text{for } q > i,\; (q - i) \text{ even},$$
   and zero otherwise. `LegendreMomentOp` precomputes $\frac{1}{r!}\bm{D}_2^r$ for $0 \le r \le \lfloor(k-1)/2\rfloor$ and forms
   $$\bm{M}(\gamma) = \sum_{r=0}^{\lfloor(k-1)/2\rfloor} \gamma^{-r} \left(\frac{1}{r!}\bm{D}_2^r\right),$$
   storing its transpose `MT(in_deg, out_deg)` so 3D action $\bm{M}(\gamma)^{\otimes 3} s$ maps directly to a single allocation-free `fast_transform(s, MT, tmp, work)` call.

---

### 3.3 Universal 1D Stencil Cache & Exact 3D Sub-Block Norms (`windowed_ns_stencil.h`)

For dyadic-locked `Coulomb`, `BSH`, and `YukawaRaw` plans, `WindowedStencilCache1D` delegates to a process-wide thread-safe table (`get_or_build_universal_stencil_entries`) keyed solely by $(k, W, M, \epsilon, \gamma_{\text{ref}}(0,0), \gamma_{\text{ref}}(W-1,M-1))$ (while `CustomGaussian` plans build exact 1D stencils only for each term's active window levels). Consequently:
- Across all $40\text{--}80$ `BSHOperatorPtr3D(mu_i)` constructions per molecular SCF run (where orbital energies $\varepsilon_i$ and $\mu_i = \sqrt{-2\varepsilon_i}$ change every iteration), **zero** 1D quadrature builds and **zero** LAPACK `dgesvd` calls occur after the initial $W \times M$ table construction (`< 0.5 ms`).
- Each `StencilEntry1D` queries `GaussianConvolution1DCache<double>::get(k, gamma_ref, 0, LatticeRange())->nonstandard(0, disp)` for $\text{disp} = l_{\text{target}} - l_{\text{source}} \in [-2, +2]$ and stores:
  - Contiguous `(dim_in, dim_out)` matrices oriented directly for `mTxmq`: `T` $(k \times k)$, `R` $(2k \times 2k)$, `R_2k_to_k` $(2k \times k)$, and `R_k_to_2k` $(k \times 2k)$.
  - Exact 1D SVDs (`StencilSVD1D`) for all four block shapes, including **rank-major contiguous prefix matrices `U_by_rank[r]`** of shape `(d_in, r)` (`ldb == r` with zero column stride padding!) and exact Frobenius tail norms `tail_normf[r]` $= (\sum_{m \ge r} \sigma_m^2)^{1/2}$.
  - Exact 3D Frobenius norms via orthogonal block subtraction—e.g. for Pass 4b-1 (`(0, d) -> s`):
    $$\|\bm{R}_{2k \to k}(\bm{d})\big|_{(0,d)\to s}\|_F = \sqrt{\prod_{\alpha \in \{x,y,z\}} \|R_{2k\to k}(d_\alpha)\|_F^2 - \prod_{\alpha \in \{x,y,z\}} \|R_{ss}(d_\alpha)\|_F^2},$$
    analogously for `opnorm_4b2_leaf` (`(s, 0) -> d`), `opnorm_4b2_int` (`(s, d) -> d`), and `opnorm_4b2_wav` (`(0, d) -> d`), aggregated into per-shell maxima `shell_max_*[sh]` and shell $L_2$ norms `shell_l2_*[sh]` for $\text{sh} = \|\bm{d}\|_\infty \in \{0, 1, 2\}$.
  - Exact 3D constant-scaling-field wavelet residual norm `const_wav_3d_norm[m_bw]` $= \|\bm{Q}^j (\sum_{\|\bm{d}\|_\infty \le m_{\text{bw}}} \bm{R}(\bm{d})^{\otimes 3} \bm{e}_{000})\|_F$ for $m_{\text{bw}} \in \{0, 1, 2\}$.

---

### 3.4 Five-Stage Execution Pipeline & Fast Paths (`windowed_ns_halo.h`, `windowed_ns_sumfact.h`, `windowed_convolution.cc`)

Given an input `Function<T, 3>` ($T \in \{\text{double}, \text{double\_complex}\}$):

#### Pure-Spiky (`pure-GC`) Fast Path (`apply_pure_spiky_same_world`)
Whenever every term in `plan` is spiky across the entire input tree ($n_{\text{spiky}}(p) \ge d_{\max}$, so there are zero active $\bm{T}^j$ or $\bm{\Gamma}^j$ terms and $w_{\text{const}} = 0$), `apply_pure_spiky_same_world` bypasses `make_nonstandard`, `WindowedNSHaloWorkspace`, `virt_cache`, and `result.reconstruct()` altogether: it evaluates the compressed leaf `GCLevelRule` (`fast_transform` + `aligned_axpy`) across leaf scaling cubes in parallel chunks and writes the result directly in `TreeState::reconstructed`.

#### Stage 1 — Active-Window Input Ingestion & Fine-Driven Interior Virtual Halo (`WindowedNSHaloWorkspace`)
1. **Stage 1a (`prepare_windowed_ns_input_3d`) & Zero-Mutation Single-Rank Ingestion (`WindowedSerialInputRange`)**:
   - **Single-Rank (`world.size() == 1`)**: Inputs in `TreeState::reconstructed`, `TreeState::compressed`, or `TreeState::nonstandard_with_leaves` are ingested **read-only with zero `WorldContainer` mutations** (no `make_nonstandard` or `standard` round-trip). Prior to workspace construction, `apply_windowed_ns_3d_prepared` computes `WindowedSerialInputRange` (`[min_s_int_level, max_s_int_level]` and `[min_gamma_level, max_gamma_level]`):
     - `ingest_local_real_nodes_serial_reconstructed`: extracts leaf scaling cubes directly from `coeffs`, restricts internal scaling cubes bottom-up strictly over `[min_s_int_level, n_max_tree - 1]`, and synthesizes `(2k)^3` NS/wavelet tensors (`twoscale_get_ns_and_wav`) strictly inside the active window `[min_gamma_level, max_gamma_level]` (skipping coarser untouched levels completely on narrow high-level windows).
     - `ingest_local_real_nodes_serial_compressed`: expands all 8 octree children of every internal node top-down from root `(0, {0,0,0})` in parallel (matching `FunctionImpl::reconstruct_op` 1:1 even on sparse `vtransform` compressed trees where small leaf nodes are omitted from `coeffs`, with exact reconstructed depth pre-scanned via `compute_serial_max_tree_level`), then populates active-window internal scaling/NS tensors.
   - **Distributed MPI (`world.size() > 1`)**: Transitions the input function to `TreeState::nonstandard_with_leaves` (`was_compressed` saved so Stage 5b restores compressed inputs in zero FLOPs via `standard()`).
2. **Global Level-Norm Reduction**:
   - Per-node scaling cubes `s_in`, average coefficient `s000 = s_in(0,0,0)`, oscillatory norm squared `s_osc_norm2` $= \|s_{\text{in}}\|_F^2 - |s_{000}|^2$, and internal pure wavelet cubes `wav_coeff` (octant $[0..k)^3$ zeroed) are reduced globally (`gop.max`) to form rank-invariant `WindowedGlobalLevelStats` (`max_s_real`, `max_s_leaf`, `max_d_int`, `max_ns_int`, `has_internal`, and `max_coarse_leaf_above[j]` $= \max_{b < j} \text{max\_s\_leaf}[b]$).
3. **Active Bandwidth & Per-Node Shell-2 Threshold Classification (`radii_builder`)**:
   - Using `loc_tol_j = f_impl->truncate_tol(eff_tol, Key<3>(j))`, each candidate term in `raw_T_terms[j]` and `raw_gamma_terms[j]` selects active half-bandwidths `bw_T` $\in \{-1, 1, 2\}$, `bw_4b1` $\in \{-1, 0, 1, 2\}$, `bw_leaf` $\in \{-1, 1, 2\}$ (**never `0`**, preserving `[-1,+1]^3` neighbor polynomial cancellation), and `bw_int >= bw_leaf` $\in \{-1, 0, 1, 2\}$.
   - Shell-1 uses tight tolerance $0.01 \cdot \text{loc\_tol}_j / \sqrt{N_{\text{terms}}}$, while Shell-2 ($\|\bm{d}\|_\infty = 2$) uses Legacy MADNESS's threshold $\text{loc\_tol}_j / (14.14 \sqrt{N_{\text{terms}}})$ and computes per-node source-norm thresholds `sh2_s_thresh`, `sh2_4b1_thresh`, and `sh2_refine_thresh` so diffuse valence nodes only emit Shell-1 (`[-1,+1]^3`) halos even on levels where a nuclear cusp requires `m = 2`.
4. **Fine-Driven Halo Discovery, Ungraded Leaf-Boundary Gating & Single-Pass Interior Prolongation (`virt_cache`)**:
   - Real nodes emit neighbor candidates tagged with bitmask flags `kHaloH1T` (Stage 4(a) $\bm{T}^{j_0}$ + Theorem 3.1 subtraction), `kHaloH1Gamma` (Pass 4b-1 `ds_corr`), `kHaloH1Refine` (Pass 4b-2 virtual wavelet promotion), and `kHaloH2Ghost` (read-only stencil completion).
   - **Ungraded (`dj >= 2`) Leaf-Boundary Gating (`ungraded_leaf_refine_only`)**: Each internal node records an 8-bit `internal_child_mask` indicating which of its 8 children at $j+1$ are themselves internal ($\Delta j \ge 2$ jump). When `eff_thresh >= 1e-6` (`ungraded_leaf_refine_only == true`), `has_ungraded_child_toward(mask, dx, dy, dz)` gates `near_internal_4b2` on adjacent real leaves so Pass 4b-2 boundary promotion probes only run across genuine ungraded ($\Delta j \ge 2$) coarse-to-fine jumps, eliminating thousands of redundant graded-boundary leaf probes on wide multi-center molecules while preserving full 2:1 outward leaf spreading at tight tolerances (`eff_thresh < 1e-6`).
   - **Serial (`world.size() == 1`)**: `build_serial` runs lock-free in parallel over 16-node slices using 64-bit packed coordinates `(lx << 42) | (ly << 21) | lz` at $j \le 20$ (merging sorted `owned_real` and `virt_at_level` into `ld.sources` via an $O(N)$ linear two-pointer merge whenever `rem_j` is empty).
   - **Distributed MPI (`world.size() > 1`)**: `build_distributed` executes **3 `WorldObject` collectives** via `WindowedNSHaloExchange<T>`: (1) route candidate `(Key<3>, flags)` to `pmap->owner(key)`; (2) walk parent chains of needed virtual/real stencil sources and send deduplicated `KeyQueryMsg` batches; (3) reply with read-only real leaf `s_in` and real internal `ns_in` tensors.
   - **Single-Pass Top-Down Prolongation (`populate_prolonged_virtual_and_sources`)**: A single bottom-up pass ($j = n_{\max} \dots 2$) inserts intermediate virtual parents up to `coarse_leaf`; then a single parallel top-down pass ($j = 1 \dots n_{\max}$) prolongs each virtual box once via `general_fast_transform` using octant filters `h_oct_[child_octant][axis]`.

#### Stage 2 — Cross-Level Flattened Leaf Gauss–Christoffel Pass
All owned real leaves ($B \in \mathcal{L}$ with $\|s_B\|_F > 0$) across $j = 0 \dots n_{\max}$ are flattened into a single vector `gc_leaves` and processed via `parallel_for_chunks`:
$$\hat{s}_B \mathrel{+}= \sum_{q=1}^{Q_{\text{eff}}(n)} \widetilde{W}_{n,q}\, \bm{M}(\tilde{\gamma}_{n,q})^{\otimes 3} s_B, \qquad n = \text{level}(B),$$
using per-chunk L1-resident $k^3$ scratch cubes (`fast_transform` + `aligned_axpy`, zero heap allocations per leaf).

#### Stage 3 — Root ($j = 0$) Constant Low-Tail & Ultra-Diffuse Pass
On the rank owning root key `(0, {0,0,0})`:
1. If $w_{\text{const}}^{\text{scaled}} \ne 0$ ($\mu = 0$), adds $w_{\text{const}}^{\text{scaled}} \, s_{0,(0,0,0)}$ directly to $\hat{s}_{0,(0,0,0)}$ (exact since $P_0 \equiv 1$ is orthogonal to all higher modes and wavelets).
2. Applies $W_p^{\text{scaled}} (\bm{T}^0_{p,00})^{\otimes 3} s_0$ via `fast_transform` + `aligned_axpy` for every ultra-diffuse term with $n_{\text{diff}}(p) < 0$.

#### Stage 4(a) — Base Scaling Injection $\bm{T}^{j_0}$, Norm-Bucket Merging & Theorem 3.1 Halo Subtraction
For each level $j \in [1, n_{\max}]$ with non-empty `active_T_terms[j]`:
1. **3-Tier Deterministic Source-Norm Bucketing & Bandwidth-Equivalent Merging (`merge_matching_norm_buckets`)**: Active scaling sources (`ld.sources`) are partitioned by $\|s\|_F / M_s$ ($M_s = \text{max\_s\_source}[j]$) into up to 3 disjoint sorted buckets (`> 0.02 M_s`, `1e-4 M_s .. 0.02 M_s`, `<= 1e-4 M_s`). Whenever active terms select identical 1D stencil bandwidths across two or more buckets, `merge_matching_norm_buckets` merges their sorted source arrays in $O(N)$ (`std::merge`) and uses the looser bucket's `bref` for SVD rank selection, eliminating redundant `SumFactGraph3D` builds and duplicate Y/Z `mTxmq` sweeps.
2. **Theorem 3.1 Boundary-Halo Commutation Subtraction**: Strictly *once* per `ActiveBaseTerm` (outside the source-norm bucket loop!), every owned virtual node $V \in \partial\mathcal{H}_j^{(1,T)}$ subtracts the interior local moment already accounted for on $V$'s coarser real leaf ancestor $B$ in Stage 2:
   $$\Delta\hat{s}_{\text{corr}}[V] \mathrel{-}= W_p^{\text{scaled}}\, \bm{M}(\gamma_{p,j})^{\otimes 3} s_{\text{virt}}[V].$$
   Where $f$ is locally polynomial inside $B$, $\sum_{\text{src}} \bm{T}^j_p s_{\text{src}} - \bm{M}(\gamma_{p,j})^{\otimes 3} s_{\text{virt}}[V] = 0$ identically; across genuine fine features adjacent to $B$, $\Delta\hat{s}_{\text{corr}}[V]$ captures the exact smooth tail spilling from $\mathcal{T}_j$ into $B$.

#### Stage 4(b) — Non-Standard Wavelet Decoupling $\bm{\Gamma}^j$ & Dynamic Output Refinement
Executed only at levels $j \in [0, n_{\max}-1]$ where `ld.global_has_internal == true`:
1. **Pass 4b-1 (`2k -> k`, `(0, d^j) -> \hat{s}^j`, Cancellation-Free)**:
   - Targets: internal real nodes (`rn.is_internal`), boundary real leaves (`rn.near_internal_4b1`), and `kHaloH1Gamma` virtual nodes (`vn.ds_corr`).
   - Sources: pure wavelet cubes `se.wav_coeff` of internal nodes (`se.is_internal`), partitioned into the 3 norm buckets and merged via `merge_matching_norm_buckets`.
   - Self-box terms (`eff_bw == 0`) run via `execute_self_box_2k_to_k_3d`; neighbor terms (`eff_bw in {1, 2}`) run via `execute_sumfact_3d` with `R_2k_to_k` / `svd_R_2k_to_k`.
2. **Constant-Subtracted Refinement Candidate Pre-Screening (`candidate_exceeds_tol`)**:
   - Internal real nodes always participate in Pass 4b-2. Boundary real leaves (`!rn.is_internal && rn.near_internal_4b2`) and virtual `kHaloH1Refine` boxes only need `(2k)^3` wavelet evaluation if their output wavelet norm can exceed $\tau_{\text{loc}}(U) = \text{truncate\_tol}(\text{eff\_thresh}, U) \cdot \text{dynamic\_refine\_factor}$.
   - `candidate_exceeds_tol(U, c0, tol_loc)` separates the nearest-neighbor shell $\text{sh} \le 1$ (`[-1,+1]^3`, 27 boxes) from outer Shell 2 ($\text{sh} = 2$). Whenever all 27 `sh <= 1` slots lie in-domain and exist in `ld.sources` (`full_sh1_stencil == true`—**never** invalidated by pruned `mLeaf_h2 == 1` ghosts at `sh == 2`!), constant polynomial cancellation around $c_0 = s_{000}(U)$ gives the tight bound:
     $$\text{ub}_{\text{sh1}} = |c_0|\,\eta_{\text{const,sh1}} + \sum_{\|\bm{d}\|_\infty \le 1} \Big[ K_{\text{leaf,sh1}}(\bm{d}) \sqrt{\|s_{\text{osc}}\|_2^2 + |s_{000} - c_0|^2} + \mathbb{I}_{\text{int}}\big(K_{\text{outer,sh1}}(\bm{d})\|s\|_F + K_{\text{wav}}(\bm{d})\|d\|_F\big) \Big],$$
     plus unsubtracted Shell-2 bounds $\text{ub}_{\text{sh2}}$. Smooth valence jump candidates are rejected (`>90%` hit rate) before allocating `(2k)^3` output cubes.
3. **Pass 4b-2 (`2k -> 2k` Internal + `k -> 2k` Leaf/Virtual -> $\hat{d}^j$)**:
   - Merged-bucket internal sources `se.ns_coeff` contract via `R` / `svd_R` (`bw_int in {0, 1, 2}`); merged-bucket leaf & virtual sources `se.s_coeff` contract via `R_k_to_2k` / `svd_R_k_to_2k` (`bw_leaf in {1, 2}`, never `0`).
   - Afterward, `zero_s0_octant` clears `[0..k)^3` on `ns_out` / `ns_virt_out` (leaving the 7 wavelet octants $\hat{d}^j$), and any boundary real leaf or `kHaloH1Refine` virtual node with $\|\hat{d}^j\|_F > \tau_{\text{loc}}$ sets `promoted_leaf = true` / `promoted = true`.
4. **Exact Shared-Tree Linearity Mode (`dynamic_refine == false`)**:
   - Replaces data-dependent level norms with `1.0` in bandwidth selection, disables per-source norm bucketing/screening, and skips output tree promotion so the discrete operator is an exact linear map ($\|\mathcal{G}(\alpha f + \beta g) - \alpha\mathcal{G}(f) - \beta\mathcal{G}(g)\|_2 = 0$ to machine precision on shared octrees).

#### Stage 5 — Direct Single-Rank Reconstructed Synthesis or Dynamic Subtree Promotion
1. **Stage 5a (`assemble_output_reconstructed_serial` / `assemble_output_nonstandard_tree`)**:
   - **Direct Single-Rank Reconstructed Synthesis (`assemble_output_reconstructed_serial`)**: When `world.size() == 1` and no virtual node or real leaf is promoted (`!any_promoted`), `restrict_all_virtual_unpromoted` skips zero/untouched `vn.ds_corr` entries and restricts non-zero halo deltas bottom-up onto coarse leaves; `assemble_output_reconstructed_serial` then runs top-down two-scale synthesis ($j = 0 \dots n_{\max}-1$) in parallel solely across active levels (`has_s_out || j >= min_gamma_level`) and inserts empty internal `FunctionNode`s plus leaf $k^3$ scaling tensors directly into `result` with `TreeState::reconstructed`, skipping `TreeState::nonstandard_after_apply` and `result.reconstruct()` completely.
   - **General Nonstandard / Distributed / Promoted Path (`assemble_output_nonstandard_tree`)**: Restricts unpromoted virtual groups in parallel, splits promoted coarse-leaf subtrees (`restrict_and_split_promoted_virtual`), deterministic-sorts MPI `CoarseLeafDeltaMsg` packets by `(leaf_key, min_source_key)`, emits $(2k)^3$ NS internal nodes and $k^3$ leaf nodes into `result` (`TreeState::nonstandard_after_apply`), and calls `result.reconstruct(fence)`.
2. **Stage 5b (`restore_windowed_ns_input_3d`)**:
   - On `world.size() == 1`, non-destructive `reconstructed` and `compressed` inputs were never modified and require zero work; on `world.size() > 1`, restores `f` (`clear` if `destructive()`, `standard()` if originally `compressed`, or `reconstruct()`).

---

## 4. High-Performance Kernel, Memory & Parallel Execution Optimizations

### 4.1 100% Tier-1 `mTxmq` (`libxsmm` JIT `s_dense_table`) Lowering & Z-Sweep Weight Folding

Every 3D tensor contraction in Stages 1–5 maps to MADNESS optimized tensor kernels with zero scalar C++ tensor loops:

| Stage / Operation | Kernel Lowering |
| :--- | :--- |
| Stage 1b Prolongation & Stage 5a Restriction / Synthesis | `general_fast_transform` (3 cyclic `mTxmq` passes) + `aligned_axpy` |
| Pure-GC Fast Path, Stage 2 Leaf GC (`Q <= 5`) & Stage 4(a) `-M(gamma)` | `fast_transform` (3 cyclic `mTxmq` passes) + `aligned_axpy` |
| Stage 3 & 4(a) Root (`j = 0`) | `fast_transform` + `aligned_axpy` |
| Stage 4(b) Self-Box (`m_bw == 0`, `2k -> k` & `2k -> 2k`) | 6-step Tier-1 `mTxmq` Tucker contraction `(2k)^3 -> r^3 -> d_out^3` using dense-packed `U_by_rank[r]` (`ldb == r`), or `general_fast_transform` / `fast_transform` when high-rank |
| Stage 4(a) & 4(b) Spatial 3D Sweeps (`m_bw in {1, 2}`) | `execute_sumfact_3d` -> `apply_slot_edge_1d`: Tier-1 `mTxmq` + L1 `fast_transpose` + `aligned_add<T, T>` |

Inside `execute_sumfact_3d` (`windowed_ns_sumfact.h`), two low-level properties are critical for throughput:
1. **Guaranteed Tier-1 `s_dense_table` Dispatch**: MADNESS's `mTxmq(dimi, dimj, dimk, c, a, b)` dispatches to its pre-JIT'd `libxsmm` `s_dense_table` kernel only when `dimi <= 400`, `dimj <= 24`, and `dimk <= 24`. In `apply_slot_edge_1d`, each cyclic sweep contracts `(d_in, dimi)` with `U_by_rank[r]` of contiguous shape `(d_in, r)` (`dimi = d_in^2 <= 400` for $k \le 10$, `dimj = r <= 20 <= 24`, `dimk = d_in <= 20 <= 24`):
   - Step 1: `mTxmq(dimi, rank, d_in, tmp_svd, src, u_ptr)` -> `[dimi, rank]` (100% Tier-1 JIT hit).
   - Step 2: If `rank == 1`, `[dimi, 1]` and `[1, dimi]` have identical memory layout, so the transpose is **elided completely**; for `rank >= 2`, L1-resident `fast_transpose(dimi, rank, tmp_svd, tmp_svd_T)` produces `[rank, dimi]`.
   - Step 3: `mTxmq(dimi, d_out, rank, target, tmp_svd_T, svt_ptr)` -> `[dimi, d_out]` (100% Tier-1 JIT hit), accumulated via `aligned_add<T, T>`.
2. **Z-Sweep Weight Folding**: By multilinearity, $\alpha (A_x \otimes A_y \otimes A_z) = A_x \otimes A_y \otimes (\alpha A_z)$. Scaling a stack-allocated 64B-aligned copy of the $\le 5$ active Z-sweep 1D factors (`<= 1,280` doubles total) by `weight` once per `execute_sumfact_3d` call allows Sweep 3 (Z) to accumulate directly into `tgt_out_ptrs[ti]`, eliminating per-term $(2k)^3 N_{\text{tgt}}$ temporary Z-buffer zeroing and `aligned_axpy(4096)` memory traffic.

### 4.2 Dual-Direction Pruned 3D Sweep Graphs, 64-Bit Packed Keys & Pooled Uninitialized Scratch

- **64-Bit Packed Spatial Coordinates (`windowed_ns_halo.h`)**: At any fixed tree level $j \le 20$ ($0 \le l_x, l_y, l_z < 2^{20} < 2^{21}$), spatial coordinates are packed into a single `uint64_t`:
  $$\text{code}(l_x, l_y, l_z) = (\text{uint64\_t}(l_x) \ll 42) \;\big|\; (\text{uint64\_t}(l_y) \ll 21) \;\big|\; \text{uint64\_t}(l_z).$$
  Unsigned `<` on `code` induces the **exact same total order** as `windowed_key_less` at fixed level $j$, shrinking keys $4\times$ (`8 B` vs `32 B` `Key<3>`) and turning 4-branch comparisons and `Key<3>` hashing/construction inside `build_serial`, `build_sumfact_graph_3d_packed`, `intersect_sorted_keys`, and `candidate_exceeds_tol` into 1-cycle integer operations (with automatic fallback to full `Key<3>` at $j > 20$).
- **Dual-Direction Graph Pruning (`build_sumfact_graph_3d`)**: Automatically selects forward candidate generation ($\mathcal{D}_{\text{src}} \to \mathcal{D}_x \to \mathcal{D}_{xy} \to \mathcal{D}_{\text{tgt}}$ with backward dead-pencil pruning) when $|\mathcal{D}_{\text{src}}| \le |\mathcal{D}_{\text{tgt}}|$ and backward candidate generation ($\mathcal{D}_{\text{tgt}} \to \mathcal{D}_{xy} \to \mathcal{D}_x$) otherwise.
- **Pooled 64-Byte-Aligned Uninitialized `SumFactScratch3D` Workspaces (`UninitAlignedBuf<T>` + `SumFactScratchLease<T>`)**: 3D sum-factorization intermediate pencil buffers (`buf_x` and `buf_xy`, up to `30–130 MB` on large multi-center grids) are backed by 64-byte-aligned `::operator new(..., std::align_val_t(64))` storage without `std::vector` constructor zero-fill (only the active slice is zeroed inside parallel chunks right before `mTxmq` accumulation) and leased from a thread-safe 16-slot pool (`trim_if_over` capping idle retention at `64 MiB` per slot). This eliminates repeated `mmap` / kernel page-fault zeroing / `munmap` TLB shootdowns across consecutive applies and orbitals.

### 4.3 Native Intel TBB Arena Dispatch & LPT Batched Vector `apply`

1. **Native Intel TBB (`HAVE_INTEL_TBB`) & MADNESS `ThreadPool` Support**:
   - When configured with `-DMADNESS_TASK_BACKEND=TBB`, both `parallel_for_chunks` (`windowed_ns_halo.h`) and multi-orbital batch execution (`apply_windowed_batch_impl` in `windowed_convolution.cc`) execute directly inside `ThreadPool::tbb_arena->execute(...)` via `tbb::parallel_for`, enabling lock-free work stealing across all TBB worker threads without nested task-queue starvation.
   - When configured with MADNESS's Pthreads task backend, high-priority `PoolTaskInterface` tasks are pushed to `ThreadPool::add` followed immediately by `thread_purge()`, with `g_active_windowed_batch_workers` dynamically donating newly idle threads to remaining deep-cusp X/Y/Z sweeps.
2. **LPT-Scheduled Batched Vector `apply` (`apply_windowed_batch_impl`)**:
   - Unique input handles (`FunctionImpl*`) are deduplicated so shared handles never race, and indices `0..N-1` on single-rank batches (`N >= 2`) are sorted by **descending `tree_size()` (Longest-Processing-Time-First, LPT)** so deep nuclear-cusp core orbitals launch immediately at $t = 0$.

---

## 5. Verification Results (`test_windowed_convolution.cc` & Legacy MRA Suites)

`test_windowed_convolution` implements an adversarial 4-Gate verification matrix covering unit algebra, unbalanced 3D coarse-to-fine jumps, 22–25 level nuclear cusps, sparse `TreeState::compressed` parity, heterogeneous batched vector `apply`, and multi-rank MPI parity:

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

3. **Gate 3 — 22–25 Level All-Electron Slater Cusps, Sparse Compressed Trees & Batched Linearity (`51 / 51 PASS`)**:
   - All-electron Slater density $\rho(\bm{r}) = \frac{Z^3}{\pi} e^{-2Z|\bm{r}-\bm{R}_A|}$ at $\bm{R}_A = (0.5000004, 0.5, 0.5)$ ($\epsilon = 10^{-8}$, analytic potential $V(\bm{R}_A) = Z$ and Hartree energy $E_H = \frac{5}{16} Z$):

     | Case | Octree Depth | Output $V$ Nodes | Global Rel $L_2$ Error | Cusp $|V(\bm{R}_A)-Z|/Z$ | Hartree $|E_H - \frac{5}{16}Z|/(\frac{5}{16}Z)$ |
     | :--- | :---: | :---: | :---: | :---: | :---: |
     | `Z=20, k=8, L=4` | `23` | `12,681` | `2.34e-09` | `2.87e-09` | `8.09e-11` |
     | `Z=20, k=10, L=4` | `22` | `4,873` | `2.68e-09` | `2.55e-10` | `4.22e-11` |
     | `Z=50, k=8, L=4` | `25` | `25,769` | `2.59e-09` | `2.91e-09` | `1.17e-10` |
     | `Z=50, k=10, L=4` | `24` | `7,273` | `2.66e-09` | `2.56e-10` | `3.27e-11` |
     | `Z=50, k=10, L=1` | `23` | `6,601` | `3.62e-10` | `2.57e-10` | `5.81e-11` |

   - Dynamic output promotion across $x = 0.5$ ($\alpha = 1.2 \times 10^{11}$ at $(0.50003, 0.5, 0.5)$ with $x < 0.5$ forced coarse at `depth = 3` in input): promoted output to **`depth = 16`** (`1.87e-10` near-midplane error, `4.73e-09` global $L_2$ error).
   - Sparse `vtransform(..., 1e-4)` `TreeState::compressed` BSH parity vs `TreeState::reconstructed` (where internal nodes omit screened leaf entries from `coeffs`): **`2.09e-15 <= 1e-13`** with input `TreeState::compressed` preserved unmutated.
   - `special_level = 23` cusp projection regression test (`project_refine_op` fix in `mraimpl.h`): **`0.26 s`**.
   - Shared-tree `Function<double_complex, 3>` linearity and heterogeneous 6-element LPT batched vector `apply` vs scalar `apply`: **`0.000000e+00 <= 1e-13`**.

4. **Gate 4 — Fast-Math + TBB Regression & Multi-Rank MPI Parity (`ALL PASS`)**:
   - `-march=native -ffast-math` + Intel TBB (`--small`, `MAD_NUM_THREADS=11`): **`63 / 63 PASSED, 0 FAILED`** in **`20.14 s`** (`0` leaked functions / `0` leaked `WorldObject`s).
   - Multi-rank MPI parity (`--full --write-parity` / `--check-parity`): `np = 1`, `np = 2`, and `np = 4` all **PASS** with `34 / 34` observables matching `np = 1` within **`1.63e-16 – 2.84e-16 <= 1e-13`** and `0` leaks.
   - Legacy MRA regression suites (`testbsh`, `testopdir`, `testsuite --small`, `test_memory_measurement`) all **PASS**, and both GCC 15.2 and Clang 21.1 compile with **`0` warnings**.

---

## 6. Performance & Timing Results

All production benchmarks below were built from scratch (including `libxsmm`'s 411 AVX2/FMA JIT microkernels) with:
```bash
cmake -GNinja \
  -DMADNESS_TARGET_ARCH=native \
  -DMADNESS_RELAXED_MATH=fast \
  -DMADNESS_TASK_BACKEND=TBB \
  -DENABLE_TCMALLOC_MINIMAL=ON
```
and executed sequentially on `1` MPI process with `MAD_NUM_THREADS=11` (`10 + main = 11` Intel TBB threads, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, Intel MKL `HAVE_MTXMQ`).

### 6.1 `madqc` All-Electron Hartree–Fock 10-Structure Benchmark Suite

- **Before**: `MADNESS_WINDOWED_NS=0` (legacy `SeparatedConvolution`)
- **After**: `MADNESS_WINDOWED_NS=1` (`WindowedSeparatedConvolution` for Coulomb, HF exchange, and BSH)
- Default `madqc --wf=scf` two-stage protocol: **`1e-4` (`k = 6`)** followed by **`1e-6` (`k = 8`)** (`--dft="restart=none; save=false"`, with `aobasis=3-21g` on `sr` since `6-31g` only tabulates $Z \le 30$).

#### End-to-End Summary (`10 / 10` Structures Faster After: `1.15x` – `4.11x`)

| Structure | $N_{\text{MO}}$ | $E_{\text{tot}}$ Before (Ha) | $E_{\text{tot}}$ After (Ha) | $|\Delta E|$ (Ha) | Iters (`1e-4`+`1e-6`) Bef $\rightarrow$ Aft | Total Wall Bef (s) | Total Wall Aft (s) | Total Speedup | Final MO Boxes (`1e-6`) Bef $\rightarrow$ Aft | MO Box Change | Peak RSS Bef $\rightarrow$ Aft (MB) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`water`** | 5 | -76.06793430 | -76.06793479 | `4.88e-7` | `8+2` $\rightarrow$ `8+2` | 28.2 | **16.3** | **1.73x** | 5,784 $\rightarrow$ **4,632** | **-19.9%** | 1016 $\rightarrow$ 1350 |
| **`water2`** | 10 | -152.14150118 | -152.14150248 | `1.30e-6` | `8+2` $\rightarrow$ `8+2` | 72.3 | **50.0** | **1.45x** | 14,008 $\rightarrow$ **11,752** | **-16.1%** | 1654 $\rightarrow$ 3211 |
| **`water5`** | 25 | -380.38040849 | -380.38041448 | `5.98e-6` | `8+2` $\rightarrow$ `8+2` | 281.4 | **244.2** | **1.15x** | 44,504 $\rightarrow$ **37,576** | **-15.6%** | 4065 $\rightarrow$ 5988 |
| **`water9`** | 45 | -684.69877329 | -684.69876323 | `1.01e-5` | `8+3` $\rightarrow$ `9+2` | 970.0 | **743.9** | **1.30x** | 86,384 $\rightarrow$ **79,760** | **-7.7%** | 8621 $\rightarrow$ 9385 |
| **`mg`** | 6 | -199.61453508 | -199.61453631 | `1.23e-6` | `5+2` $\rightarrow$ `5+2` | 32.4 | **14.7** | **2.20x** | 6,384 $\rightarrow$ **4,656** | **-27.1%** | 1096 $\rightarrow$ 1360 |
| **`ca`** | 10 | -676.75807709 | -676.75807987 | `2.79e-6` | `5+3` $\rightarrow$ `5+2` | 118.3 | **36.0** | **3.29x** | 14,816 $\rightarrow$ **9,032** | **-39.0%** | 1984 $\rightarrow$ 2270 |
| **`sr`** | 19 | -3131.54555063 | -3131.54557463 | `2.40e-5` | `12+5` $\rightarrow$ `10+3` | 900.6 | **218.9** | **4.11x** | 42,680 $\rightarrow$ **21,552** | **-49.5%** | 4818 $\rightarrow$ **3158** |
| **`c2h4`** | 8 | -78.06952528 | -78.06952574 | `4.63e-7` | `8+2` $\rightarrow$ `8+2` | 55.8 | **38.0** | **1.47x** | 11,328 $\rightarrow$ **10,368** | **-8.5%** | 1421 $\rightarrow$ 2669 |
| **`ne`** | 5 | -128.54700005 | -128.54700036 | `3.08e-7` | `6+2` $\rightarrow$ `6+2` | 24.3 | **11.7** | **2.08x** | 4,232 $\rightarrow$ **3,688** | **-12.9%** | 968 $\rightarrow$ 1112 |
| **`ar`** | 9 | -526.81740356 | -526.81740804 | `4.49e-6` | `6+2` $\rightarrow$ `6+2` | 75.0 | **31.6** | **2.37x** | 13,376 $\rightarrow$ **7,784** | **-41.8%** | 1511 $\rightarrow$ 2049 |

#### Protocol 1 Breakdown (`thresh = 1e-4, k = 6`)

| Structure | Iters Bef $\rightarrow$ Aft | Proto Span Bef $\rightarrow$ Aft (s) | Proto Speedup | Coulomb Bef $\rightarrow$ Aft (s) | HF Exchange Bef $\rightarrow$ Aft (s) | Apply BSH Bef $\rightarrow$ Aft (s) | MO Leaf Boxes Bef $\rightarrow$ Aft |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`water`** | `8` $\rightarrow$ `8` | 14.70 $\rightarrow$ **9.60** | **1.53x** | 0.84 $\rightarrow$ **0.47** | 9.26 $\rightarrow$ **4.81** (**1.93x**) | 2.16 $\rightarrow$ **1.94** | 3,240 $\rightarrow$ 3,240 |
| **`water2`** | `8` $\rightarrow$ `8` | 35.00 $\rightarrow$ **27.30** | **1.28x** | 1.18 $\rightarrow$ **0.83** | 23.82 $\rightarrow$ **16.10** (**1.48x**) | 4.65 $\rightarrow$ 5.23 | 6,696 $\rightarrow$ **6,664** |
| **`water5`** | `8` $\rightarrow$ `8` | 129.20 $\rightarrow$ **127.90** | **1.01x** | 2.32 $\rightarrow$ **2.07** | 96.26 $\rightarrow$ **88.15** (**1.09x**) | 13.80 $\rightarrow$ 21.29 | 18,136 $\rightarrow$ **18,112** |
| **`water9`** | `8` $\rightarrow$ `9` | 329.50 $\rightarrow$ 399.10 | 0.83x | 3.65 $\rightarrow$ 4.38 | 259.86 $\rightarrow$ 292.18 | 27.45 $\rightarrow$ 61.73 | 34,472 $\rightarrow$ **34,384** |
| **`mg`** | `5` $\rightarrow$ `5` | 13.20 $\rightarrow$ **7.50** | **1.76x** | 0.57 $\rightarrow$ **0.25** | 9.10 $\rightarrow$ **4.17** (**2.18x**) | 1.60 $\rightarrow$ **1.20** | 3,888 $\rightarrow$ 3,888 |
| **`ca`** | `5` $\rightarrow$ `5` | 31.70 $\rightarrow$ **17.30** | **1.83x** | 0.62 $\rightarrow$ **0.32** | 24.84 $\rightarrow$ **11.53** (**2.15x**) | 3.07 $\rightarrow$ **2.38** | 6,736 $\rightarrow$ 6,736 |
| **`sr`** | `12` $\rightarrow$ `10` | 264.70 $\rightarrow$ **113.00** | **2.34x** | 2.10 $\rightarrow$ **0.87** | 227.98 $\rightarrow$ **89.34** (**2.55x**) | 21.78 $\rightarrow$ **12.18** | 13,720 $\rightarrow$ 13,720 |
| **`c2h4`** | `8` $\rightarrow$ `8` | 27.40 $\rightarrow$ **20.40** | **1.34x** | 1.07 $\rightarrow$ **0.70** | 18.41 $\rightarrow$ **11.48** (**1.60x**) | 3.45 $\rightarrow$ 3.78 | 5,760 $\rightarrow$ **5,696** |
| **`ne`** | `6` $\rightarrow$ `6` | 11.20 $\rightarrow$ **6.60** | **1.70x** | 0.61 $\rightarrow$ **0.29** | 7.28 $\rightarrow$ **3.34** (**2.18x**) | 1.52 $\rightarrow$ **1.18** | 3,048 $\rightarrow$ 3,048 |
| **`ar`** | `6` $\rightarrow$ `6` | 32.20 $\rightarrow$ **17.00** | **1.89x** | 0.76 $\rightarrow$ **0.36** | 25.02 $\rightarrow$ **11.15** (**2.24x**) | 3.42 $\rightarrow$ **2.43** | 6,024 $\rightarrow$ 6,024 |

#### Protocol 2 Breakdown (`thresh = 1e-6, k = 8`: `1.32x` – `6.18x` Faster Across All 10 Structures)

| Structure | Iters Bef $\rightarrow$ Aft | Proto Span Bef $\rightarrow$ Aft (s) | Proto Speedup | Coulomb Bef $\rightarrow$ Aft (s) | HF Exchange Bef $\rightarrow$ Aft (s) | Apply BSH Bef $\rightarrow$ Aft (s) | MO Leaf Boxes Bef $\rightarrow$ Aft |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`water`** | `2` $\rightarrow$ `2` | 12.10 $\rightarrow$ **5.40** | **2.24x** | 0.84 $\rightarrow$ **0.34** (**2.47x**) | 8.18 $\rightarrow$ **3.08** (**2.66x**) | 1.90 $\rightarrow$ **0.94** (**2.02x**) | 5,784 $\rightarrow$ **4,632** (-19.9%) |
| **`water2`** | `2` $\rightarrow$ `2` | 34.60 $\rightarrow$ **20.60** | **1.68x** | 1.45 $\rightarrow$ **0.68** (**2.13x**) | 25.89 $\rightarrow$ **13.98** (**1.85x**) | 4.47 $\rightarrow$ **3.49** (**1.28x**) | 14,008 $\rightarrow$ **11,752** (-16.1%) |
| **`water5`** | `2` $\rightarrow$ `2` | 145.96 $\rightarrow$ **110.18** | **1.32x** | 3.41 $\rightarrow$ **2.03** (**1.68x**) | 118.47 $\rightarrow$ **84.59** (**1.40x**) | 14.21 $\rightarrow$ 14.55 | 44,504 $\rightarrow$ **37,576** (-15.6%) |
| **`water9`** | `3` $\rightarrow$ `2` | 625.85 $\rightarrow$ **332.32** | **1.88x** | 8.69 $\rightarrow$ **4.15** (**2.09x**) | 516.00 $\rightarrow$ **266.10** (**1.94x**) | 69.20 $\rightarrow$ **39.42** (**1.76x**) | 86,384 $\rightarrow$ **79,760** (-7.7%) |
| **`mg`** | `2` $\rightarrow$ `2` | 17.40 $\rightarrow$ **5.60** | **3.11x** | 0.93 $\rightarrow$ **0.22** (**4.23x**) | 12.69 $\rightarrow$ **3.52** (**3.61x**) | 2.50 $\rightarrow$ **0.89** (**2.81x**) | 6,384 $\rightarrow$ **4,656** (-27.1%) |
| **`ca`** | `3` $\rightarrow$ `2` | 84.30 $\rightarrow$ **16.40** | **5.14x** | 2.63 $\rightarrow$ **0.37** (**7.11x**) | 66.37 $\rightarrow$ **12.12** (**5.48x**) | 11.34 $\rightarrow$ **2.07** (**5.48x**) | 14,816 $\rightarrow$ **9,032** (-39.0%) |
| **`sr`** | `5` $\rightarrow$ `3` | 631.90 $\rightarrow$ **102.20** | **6.18x** | 9.23 $\rightarrow$ **0.69** (**13.38x**) | 522.59 $\rightarrow$ **86.46** (**6.04x**) | 82.05 $\rightarrow$ **9.13** (**8.99x**) | 42,680 $\rightarrow$ **21,552** (-49.5%) |
| **`c2h4`** | `2` $\rightarrow$ `2` | 25.99 $\rightarrow$ **15.32** | **1.70x** | 1.13 $\rightarrow$ **0.54** (**2.09x**) | 19.14 $\rightarrow$ **10.15** (**1.89x**) | 3.50 $\rightarrow$ **2.49** (**1.41x**) | 11,328 $\rightarrow$ **10,368** (-8.5%) |
| **`ne`** | `2` $\rightarrow$ `2` | 11.90 $\rightarrow$ **4.06** | **2.93x** | 0.81 $\rightarrow$ **0.22** (**3.68x**) | 8.26 $\rightarrow$ **2.40** (**3.44x**) | 1.90 $\rightarrow$ **0.67** (**2.84x**) | 4,232 $\rightarrow$ **3,688** (-12.9%) |
| **`ar`** | `2` $\rightarrow$ `2` | 40.70 $\rightarrow$ **12.80** | **3.18x** | 1.47 $\rightarrow$ **0.33** (**4.45x**) | 32.93 $\rightarrow$ **9.25** (**3.56x**) | 4.02 $\rightarrow$ **1.73** (**2.32x**) | 13,376 $\rightarrow$ **7,784** (-41.8%) |

---

### 6.2 3D Multi-Scale Gaussian Convolution & Octree Benchmark (`src/examples/large_exp_3d.cc`)

Extracted from `large_exp.cc::test6` on $[-11, +11]^3$ (`width = 22.0`, $k = 10$, `thresh = 1e-8`, `truncate_mode = 0`), convolving normalized 3D Gaussian $f(\bm{r}) = (a/\pi)^{3/2} e^{-a|\bm{r}-\bm{r}_0|^2}$ ($a = 16.0$, $\bm{r}_0 = (0.1, 0, 0)$, `652` input leaves on levels `[2..7]`, $h_{\text{user}} = 0.171875$) with $g_b(\bm{r}) = (b/\pi)^{3/2} e^{-b r^2}$ across `iscale = 0 .. 40` ($b = 2^{\text{iscale}-1} \in [0.5, 5.50 \times 10^{11}]$):

| `isc` | `userh/sig` | `expnt_b` | Window Regime | Exact Leaves | Bef Raw Leaves | Bef Trunc Leaves | Bef $L_2$ Err | Bef Time (ms) | Aft Raw Leaves | Aft Trunc Leaves | Aft $L_2$ Err | Aft Time (ms) | Raw Leaf Ratio | Speedup |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0 | `1.72e-1` | `5.00e-1` | `win[0..3]` | 176 | 512 | 176 | `1.44e-9` | 9.92 | 652 | **176** | `3.79e-8` | **9.29** | 0.79x | **1.07x** |
| 2 | `3.44e-1` | `2.00e+0` | `win[1..4]` | 344 | 904 | 344 | `2.32e-9` | 11.43 | **764** | **344** | `2.97e-8` | 12.82 | **1.18x** | 0.89x |
| 4 | `6.88e-1` | `8.00e+0` | `win[2..5]` | 568 | 1,520 | 568 | `2.43e-9` | 13.55 | **792** | **568** | `1.71e-8` | 15.76 | **1.92x** | 0.86x |
| 6 | `1.38e+0` | `3.20e+1` | `win[3..6]` | 512 | 2,220 | 624 | `4.74e-8` | 24.04 | **652** | **512** | `8.16e-9` | **20.35** | **3.40x** | **1.18x** |
| 8 | `2.75e+0` | `1.28e+2` | `win[4..7]` | 596 | 2,640 | 1,100 | `8.86e-6` | 20.87 | **652** | **596** | `2.02e-8` | 20.91 | **4.05x** | **1.00x** |
| 10 | `5.50e+0` | `5.12e+2` | `win[5..7]` | 652 | 2,864 | 1,268 | `1.56e-4` | 21.43 | **652** | **652** | `2.42e-8` | **17.35** | **4.39x** | **1.24x** |
| 12 | `1.10e+1` | `2.05e+3` | `win[6..7]` | 652 | 2,836 | 1,380 | `5.79e-4` | 20.49 | **652** | **652** | `2.55e-8` | **16.00** | **4.35x** | **1.28x** |
| 13 | `1.56e+1` | `4.10e+3` | `win[7..7]` | 652 | 2,892 | 1,380 | `7.74e-4` | 20.94 | **652** | **652** | `3.06e-8` | **9.49** | **4.44x** | **2.21x** |
| 14 | `2.20e+1` | `8.19e+3` | `win[7..7]` | 652 | 2,864 | 1,380 | `1.86e-3` | 20.73 | **652** | **652** | `2.70e-8` | **9.36** | **4.39x** | **2.21x** |
| 15 | `3.11e+1` | `1.64e+4` | `pure-GC` | 652 | 2,808 | 1,380 | `2.92e-3` | 20.10 | **652** | **652** | `2.96e-8` | **3.30** | **4.31x** | **6.09x** |
| 17 | `6.22e+1` | `6.55e+4` | `pure-GC` | 652 | 2,696 | 1,380 | `3.48e-3` | 19.16 | **652** | **652** | `2.56e-8` | **3.39** | **4.13x** | **5.65x** |
| 24 | `7.04e+2` | `8.39e+6` | `pure-GC` | 652 | 2,276 | 1,100 | `6.36e-4` | 18.01 | **652** | **652** | `2.54e-8` | **3.58** | **3.49x** | **5.03x** |
| 32 | `1.13e+4` | `2.15e+9` | `pure-GC` | 652 | 1,520 | 988 | `4.23e-5` | 13.76 | **652** | **652** | `2.54e-8` | **3.54** | **2.33x** | **3.89x** |
| 40 | `1.80e+5` | `5.50e+11` | `pure-GC` | 652 | 1,492 | 876 | `2.65e-6` | 11.58 | **652** | **652** | `2.54e-8` | **3.63** | **2.29x** | **3.19x** |
| **SUM / WORST (`0..40`)** | — | — | **41 scales** | **24,912** | **83,432** | **41,348** | **`3.48e-3`** | **684.15** | **27,320** | **24,912** | **`3.79e-8`** | **328.21** | **3.05x** | **2.08x** |

- **Zero Large-Exponent Loss of Precision**: Legacy `SeparatedConvolution` degrades to `3.48e-3` $L_2$ error once `userh/sigma > 1` (`iscale >= 6`), whereas `WindowedSeparatedConvolution` maintains **`2.14e-9 – 3.79e-8` $L_2$ error across all 41 scales** (`91,800x` lower worst-case error).
- **Exact Analytical Output Trees**: Legacy raw trees bloat up to `4.44x` (`2,892` vs `652` leaves) and retain `1,380` unprunable leaves (`+768` spurious `L7` interface boxes) even after `.truncate()`. Windowed NS matches the exact analytical truncated leaf count on **all 41/41 scales (`24,912` vs `24,912` total leaves)**.
- **Smooth Intermediate-to-Spiky Speedup**: Active-window read-only ingestion, direct `TreeState::reconstructed` output synthesis, pooled uninitialized scratch, and norm-bucket merging accelerate `iscale = 13..14` (`win[7..7]`) to **`2.21x`**, `iscale = 6..12` to **`1.00x–1.28x`**, and `iscale >= 15` (`pure-GC`) to **`3.19x–6.09x`** (`2.08x` across the entire 41-scale sweep).

---

### 6.3 Deep 22-Level `Z = 30` Slater Cusp Micro-Benchmark (`depth = 22, eps = 1e-7, rho_nodes = 7,625`)

| Operator & Runtime Configuration | Legacy MADNESS | **Windowed NS** | **Speedup vs Legacy & Tree Compression** |
| :--- | :---: | :---: | :--- |
| `Coulomb` (`np = 1, TBB MAD_NUM_THREADS=11, -ffast-math`) | `1.227 s` (`18,313` nodes) | **`0.569 s` (`7,689` nodes)** | **`2.16x` faster** (**`-58.0%` output nodes**) |
| `BSH` (`mu = 1.0, np = 1, TBB MAD_NUM_THREADS=11, -ffast-math`) | `0.816 s` (`14,409` nodes) | **`0.505 s` (`7,625` nodes)** | **`1.62x` faster** (**`-47.1%` output nodes**) |
| `Coulomb` (`np = 1, Pthreads MAD_NUM_THREADS=4, -O2`) | `1.780 s` (`18,313` nodes) | **`0.828 s` (`7,689` nodes)** | **`2.15x` faster** (**`-58.0%` output nodes**) |
| `BSH` (`mu = 1.0, np = 1, Pthreads MAD_NUM_THREADS=4, -O2`) | `1.180 s` (`14,409` nodes) | **`0.572 s` (`7,625` nodes)** | **`2.06x` faster** (**`-47.1%` output nodes**) |
| `Coulomb` (`np = 2, MAD_NUM_THREADS=4, -O2`) | `1.753 s` | **`1.262 s`** | **`1.39x` faster** |
| `BSH` (`mu = 1.0, np = 2, MAD_NUM_THREADS=4, -O2`) | `1.065 s` | **`0.883 s`** | **`1.21x` faster** |
| `Coulomb` (`np = 4, MAD_NUM_THREADS=1, -O2`) | `2.460 s` | **`1.558 s`** | **`1.58x` faster** |
| `BSH` (`mu = 1.0, np = 4, MAD_NUM_THREADS=1, -O2`) | `1.656 s` | **`1.116 s`** | **`1.48x` faster** |
