# Tail-Completed, Two-Threshold Windowed Non-Standard 3D Convolution in MADNESS

This document describes the mathematical formulation, data structures, parallel/distributed execution pipeline, `mTxmq` kernel lowering, verification suite, and benchmark performance of the **Tail-Completed, Two-Threshold Windowed Non-Standard (NS) 3D Convolution** implemented in `src/madness/mra/`.

---

## 1. Motivation & Design Overview

Applying isotropic 3D integral operators represented as sums of separated Gaussians—most notably the free-space Coulomb kernel $1/r$, the bound-state Helmholtz (BSH) kernel $e^{-\mu r}/(4\pi r)$, the Yukawa kernel $e^{-\mu r}/r$, and general multi-scale Gaussian sums—on deeply refined, unbalanced octrees presents two classical bottlenecks in standard multiwavelet NS form:

1. **Gibbs (`erf`) Interface Oscillations Across Coarse-to-Fine Jumps**: When a coarse leaf box abuts a finer subtree across a face, edge, or corner without 2:1 tree grading, extrapolating high-order Legendre polynomials *outside* the coarse leaf's support creates artificial boundary discontinuities. Convolving those discontinuities with sharp Gaussians induces spurious error-function (`erf`) oscillations that force unnecessary adaptive refinement in Molecular Orbitals (MOs) and lose up to 5–6 digits of precision at intermediate-to-large exponents.
2. **Deep-Tree Cusp Quadrature Cost**: Resolving an all-electron nuclear cusp of charge $Z$ to accuracy $\epsilon$ requires octree depths $n_{\max} \approx 20\text{--}26$ and $O(10^2)$ Gaussian exponents spanning $15+$ orders of magnitude. Sweeping every Gaussian across every tree level scales poorly with cusp depth and wastes work on levels where a Gaussian is either ultra-sharp ("spiky", acting as a local differential moment operator) or ultra-smooth ("diffuse", already resolved on coarser scales).

The Windowed NS algorithm resolves both bottlenecks simultaneously by enforcing **five invariants**:

- **Invariant 1 — Zero Exterior Polynomial Extrapolation**: Every virtual halo box $V$ at level $j$ is strictly a sub-box of an existing coarser real leaf $B \in \mathcal{L}$ ($\text{level}(B) < j$), prolonged *solely inside* $B$'s support one dyadic step at a time using exact two-scale Legendre filters $(h_0, h_1)$.
- **Invariant 2 — Two-Threshold Scale Separation & Constant Window Width $W$**: At tree level $j$, a Gaussian of simulation-cell exponent $\alpha_{\text{sim},p} = e^{2t_p}$ has dimensionless stiffness $\gamma_{p,j} = \alpha_{\text{sim},p} 4^{-j}$. It is active in 3D NS sum-factorization **only** inside a narrow scale window $j \in [j_0(p), j_1(p)]$ of constant width $W \approx 4\text{--}6$ levels bounded by an exact Plancherel diffuse threshold $\gamma_{\text{diff}}(k,\epsilon)$ and a closed-form spiky threshold $\gamma_{\text{spiky}}(k) = k^2(k+1)^2/30$ (floored at $270.0$ for $k < 8$ operators via `select_operator_gamma_spiky_override` so Gaussian boundary flux across leaf faces stays $< 0.5\%$ and preserves two-scale KAIN subspace linearity).
- **Invariant 3 — Overflow-Free Gauss–Christoffel Spiky Compression & Discrete Tail Completion**: All Gaussians spiky at leaf level $n$ ($n_{\text{spiky}}(p) \ge n$, including the analytic $p > p_{\text{tail}}$ cusp tail) are compressed via level-rescaled Stieltjes–Golub–Welsch quadrature into $Q_{\text{eff}} \le Q$ local Legendre moment operators $\bm{M}(\tilde{\gamma}_{n,q})^{\otimes 3}$ per leaf ($Q = 4$ for $k \le 6$, $Q = 5$ for $7 \le k \le 14$, and $\max(5, \lfloor(m_k + 4)/2\rfloor)$ with $m_k = \lfloor(k-1)/2\rfloor + 1$ for $k > 14$). At the opposite end ($\mu = 0$), the omitted $p < p_{\min}$ infrared tail is completed by matching exact discrete geometric moments $(M_0, M_1, M_2)$ with a single synthetic Gaussian plus a root $(0,0,0)$ constant shift $w_{\text{const}}$, eliminating Euler–Maclaurin tail bumps.
- **Invariant 4 — Exact Boundary-Halo Commutation (Theorem 3.1)**: When a Gaussian enters its NS window at base level $j_0(p)$, any coarser leaf $B$ ($\text{level}(B) < j_0(p)$) has *already* accumulated $p$'s local moment contribution inside Stage 2's Gauss–Christoffel leaf rule. Subtracting $W_p \bm{M}(\gamma_{p,j_0})^{\otimes 3} s_{\text{virt}}[V]$ on primary virtual halo boxes $V \in \partial\mathcal{H}_{j_0}^{(1,T)}$ before adjoint restriction cancels the interior polynomial contribution to machine precision while preserving exact smooth cross-interface flux.
- **Invariant 5 — Dyadic Scale-Locking ($\beta = 4^{1/M}$)**: Choosing quadrature spacing $\Delta t = \frac{\ln 2}{M}$ locks $\gamma_{p,j}$ to a universal table of $W \times M$ (`~10–35`) reference stiffnesses $\gamma_{\text{ref}}(w_{\text{eff}}, m)$ independent of tree depth $n_{\max}$, screening parameter $\mu$, and cubic cell width $L$.

---

## 2. Source Layout & Public Entry Points

### 2.1 File Map

| File | Role in Implementation |
| :--- | :--- |
| `windowed_quadrature.h` | Two-threshold solver (`compute_gamma_spiky`, `compute_gamma_diff`, `gamma_spiky_override`), parity-sparse `LegendreMomentOp` ($\bm{M}(\gamma)$ / $\bm{M}(\gamma)^{\mathsf{T}}$), dyadic-locked `WindowedQuadraturePlan` (`Coulomb`, `BSH`, `YukawaRaw`, plus overflow-safe off-grid `CustomGaussian`), precomputed remaining operator weight bounds `rem_weight_bounds_` ($W_{\text{rem}}(j)$), discrete low/high tail completion, and level-rescaled Stieltjes–Golub–Welsch `GCLevelRule` compression |
| `windowed_ns_stencil.h` | 64-bit SplitMix64 `WindowedKeyHash` (`Level <= 62`), universal process-wide `(k, W, M, eps)` memoized `WindowedStencilCache1D` (plus `offgrid_entries_` for off-grid terms and `root_T0_by_term_` for ultra-diffuse $n_{\text{diff}} < 0$ terms), `StencilEntry1D` transition blocks (`T`, `R`, `R_2k_to_k`, `R_k_to_2k`), contiguous `U_by_rank` 1D SVDs (`StencilSVD1D`), and exact 3D NS sub-block/shell, constant-wavelet (`const_wav_3d_norm`), and Theorem-3.1 constant-remainder (`const_T_rem_3d_norm`) norms |
| `windowed_ns_halo.h` | Worker-thread slot registry (`WindowedThreadSlotRegistry`) & native MADNESS task-model dispatch (`parallel_for_chunks`, `windowed_for_each_outer_tasks` over `madness::Range` + `WorldTaskQueue` / `tbb_arena`), 64-bit packed spatial key indexing (`j <= 20`, exact `Key<3>` fallback for `j > 20`), zero-mutation serial `TreeState::reconstructed` / `TreeState::compressed` ingestion (`WindowedSerialInputRange`) with bottom-up $H^1$-safe `subtree_d_tol` pruning, per-node operator-weighted Shell-1 (`node_mT_and_l1`, `node_m4b1`, `node_mRefine`) & Shell-2 (`h2_mT`, `h2_mLeaf`) halo gating with coarse-leaf Cauchy–Schwarz bounds (`max_s_leaf_linf`), ungraded (`dj >= 2`) leaf-boundary gating (`has_ungraded_child_toward`), slot-tagged multi-orbital `WindowedNSHaloExchange<T>` (`num_slots = B`), `WindowedNSHaloWorkspace` (serial & lockstep distributed `build_distributed_lockstep` / `assemble_output_nonstandard_lockstep`), single-pass interior `virt_cache` prolongation, bottom-up adjoint restriction, direct single-rank `assemble_output_reconstructed_serial`, and dynamic subtree promotion |
| `windowed_ns_sumfact.h` | Smooth loose-tolerance tightening (`tighten_loose_windowed_eps`), thread-safe pooled 64B-aligned uninitialized scratch (`UninitAlignedBuf`, `SumFactScratchLease`), FLOP-calibrated sweep grain (`compute_sumfact_sweep_grain`), pruned dual-direction 3D sweep graphs (`build_sumfact_graph_3d`), Tier-1 `mTxmq` (`libxsmm` JIT) + L1 `fast_transpose` + fused multi-term bandwidth-grouped `execute_sumfact_batch_3d` (`kMaxFusedTerms = 8`, `4 MiB` pencil cap) & `execute_self_box_*_batch_3d`, fused Theorem 3.1 `h1t_virts`, 3-tier source-norm bucketing with $O(N)$ bandwidth-equivalent merging (`merge_matching_norm_buckets`), `candidate_exceeds_tol` constant-subtracted screening, shared `WindowedOrbitalExecContext<T>`, and single-orbital / lockstep distributed `apply_windowed_ns_3d_batch_distributed<T>` |
| `windowed_convolution.h`, `windowed_convolution.cc` | `WindowedSeparatedConvolution<double,3>` (including `select_operator_gamma_spiky_override`, custom Gaussian sum constructors, and `apply_pure_spiky_same_world` pure-GC fast path), factory functions (`WindowedCoulombOperator[Ptr]`, `WindowedBSHOperator[Ptr][3D]`, `WindowedYukawaOperator[Ptr]`, `WindowedGaussianOperator[Ptr]`), thread-safe plan/stencil cache, distributed `prune_ex` compressed union-skeleton pre-pruning, single-rank adaptive LPT task-queued & multi-rank lockstep batched vector `apply`, and opt-in `SeparatedConvolution` bridge |
| `operator.h`, `mra.h`, `vmra.h` | Non-intrusive opt-in bridge (`SeparatedConvolution::enable_windowed_ns(bool)` and `MADNESS_WINDOWED_NS=1`) for scalar and vector `apply(...)` across `OT_G12`, `OT_BSH`, and isotropic 3D `CustomGaussian` operators, plus `OT_BSH` metadata population in `BSHOperatorPtr3D` |
| `mraimpl.h` | Standalone `project_refine_op` fix (`else if` in `newspecialpts` filtering) eliminating $O(2^{\text{special\_level}})$ special-point list duplication on deep nuclear cusps |
| `test_windowed_convolution.cc` | 4-Gate unit, unbalanced interface jump, 22–25 level nuclear cusp, sparse `vtransform` compressed-tree parity, heterogeneous batched `apply`, and multi-rank MPI parity suite (`--small`, `--full`, `--gate=1..4`, `--write-parity`, `--check-parity`) |
| `src/examples/large_exp_3d.cc` | Standalone 3D Gaussian convolution & octree inspection benchmark sweeping 41 kernel scales ($b \in [0.5, 5.50 \times 10^{11}]$) across diffuse, intermediate windowed, and ultra-spiky `pure-GC` regimes |
| `src/examples/hatom_sf_dirac.cc`, `src/examples/dirac-hatom.cc` | Relativistic hydrogenic Dirac example solvers (2-component spin-free eliminated Dirac with `BSHOperator3D` + `pVp` + `CoulombOperator`, and 4-component complex spinor Dirac with `BSHOperator<3>` + $c\,\boldsymbol{\alpha}\cdot\mathbf{p}$) instrumented with per-stage wall-clock timers, raw-vs-truncated octree leaf/node/max-depth diagnostics, Peak RSS tracking, and clean `mask.clear()` teardown before `finalize()` |

### 2.2 Usage & Runtime Configuration

1. **Direct C++ API (`windowed_convolution.h`)**:
   ```cpp
   auto pop = WindowedCoulombOperatorPtr(world, eps);               // 1 / r
   auto bop = WindowedBSHOperatorPtr(world, mu, eps);               // exp(-mu*r) / (4*pi*r)
   auto yop = WindowedYukawaOperatorPtr(world, mu, eps);            // exp(-mu*r) / r
   auto gop = WindowedGaussianOperatorPtr(world, coeff, expnt, eps); // sum_p coeff[p]*exp(-expnt[p]*r^2)

   real_function_3d V = apply(*pop, rho);                           // scalar real or complex
   vecfuncT Bpsi = apply(world, bop_ptrs, psi);                     // LPT-batched / MPI-lockstep vector apply
   ```
2. **Per-Operator Bridge on Legacy `SeparatedConvolution<double,3>`**:
   ```cpp
   SeparatedConvolution<double, 3> op = CoulombOperator(world, lo, eps);
   op.enable_windowed_ns(true);                                // delegates scalar & vector apply
   ```
3. **Global Environment Opt-In & Optional Tuning Knobs (Zero Application Code Changes)**:
   ```bash
   MADNESS_WINDOWED_NS=1 ./madqc --geometry=water
   ```
   When `MADNESS_WINDOWED_NS` is set to `1` (`true`/`on`/`yes`) or `op.enable_windowed_ns(true)` is called, `SeparatedConvolution::windowed_ns_enabled()` (`operator.h`) automatically delegates scalar `apply(op, f)` (`mra.h`) and batched `apply(world, op, vec_f)` / `apply(world, vec_ops, vec_f)` (`vmra.h`) to `WindowedSeparatedConvolution<double, 3>` whenever `Q == double`, `NDIM == 3`, `f.k() == op.get_k()`, the simulation cell is cubic (`cw(0) == cw(1) == cw(2) > 0`) with free-space boundary conditions, and the operator is an unmodified, non-periodic, non-range-restricted `particle == 1` `OT_G12` ($\mu \ge 0$), `OT_BSH` ($\mu \ge 0$, without `truncate_lowexp_gaussians`), or 1D `stored_coeffs`/`stored_exponents` Gaussian sum with finite positive exponents; all other operators remain on the legacy code path.
   - `MADNESS_WINDOWED_MPI_BATCH_TILE` (default `16`, clamped to `[1, 256]`): controls the number of distinct orbitals executed simultaneously in lockstep per MPI collective batch on `world.size() > 1` (can be lowered to e.g. `4` or `8` on memory-constrained nodes with very deep trees).
   - `MADNESS_WINDOWED_BATCH_WORKERS` (optional diagnostic override): overrides the adaptive single-rank outer task tile size `outer_batch_tile` (`N` when `max_nodes < 1,024`, `min(N, 4)` when `< 16,384`, `min(N, 2)` when `< 32,768`, and `1` when `>= 32,768`).

---

## 3. Detailed Algorithm as Implemented

### 3.1 Two-Threshold Scale Separation, Boundary-Flux Floor & Loose-Tolerance Tightening (`windowed_quadrature.h`, `windowed_ns_sumfact.h`, `windowed_convolution.cc`)

Let $L > 0$ be the cubic simulation cell width (`FunctionDefaults<3>::get_cell_width()`), $\mu_{\text{sim}} = \mu_{\text{user}} L$, and $k$ be the Legendre scaling order on each octree box. In unit simulation coordinates $\bm{x} \in [0, 1]^3$, a 3D Gaussian convolution term $c_p e^{-\alpha_{\text{sim},p} \|\bm{x}-\bm{y}\|^2}$ with $\alpha_{\text{sim},p} = e^{2t_p}$ has physical 3D operator weight
$$W_p^{\text{scaled}} = \begin{cases} \text{op\_scale} \cdot c_p \left(\dfrac{\pi}{\alpha_{\text{sim},p}}\right)^{3/2}, & \text{Coulomb / YukawaRaw } (\text{op\_scale} = L^2) \text{ or BSH } (\text{op\_scale} = L^2 / (4\pi)), \\[6pt] c_{\text{user},p} \left(\dfrac{\pi}{\alpha_{\text{user},p}}\right)^{3/2}, & \text{CustomGaussian } (t_p = \frac{1}{2}\ln\alpha_{\text{user},p} + \ln L,\; c_{\text{sim},p} = L\,c_{\text{user},p}). \end{cases}$$
For dyadic terms, `c_sim` and `weight_sim` are evaluated with log-argument range guards `(-740, 705)` while Stage 2's level-$n$ GC rule computes level-rescaled weights $\widehat{W}_{p,n} = 4^n W_p^{\text{scaled}}$ directly from $u_{p,n} = \exp(2n\ln 2 - 2t_p)$ without forming $W_p^{\text{scaled}}$ first, preventing IEEE-754 `double` underflow/overflow at arbitrarily large $\alpha_{\text{sim},p}$ or deep cusp levels $n$.

At octree level $j$ (box width $2^{-j}$), dimensionless stiffness is $\gamma_{p,j} = \alpha_{\text{sim},p} 4^{-j} = e^{2t_p - 2j\ln 2}$. Two thresholds partition all scales:

1. **Spiky Threshold $\gamma_{\text{spiky}}(k)$ & Low-Order Boundary-Flux Floor (`compute_gamma_spiky`, `select_operator_gamma_spiky_override`)**:
   $$\gamma_{\text{spiky}}^{\text{interior}}(k) = \frac{k^2 (k + 1)^2}{30}, \qquad \gamma_{\text{spiky}}^{\text{op}}(k) = \begin{cases} 270.0, & k < 8, \\ \dfrac{k^2(k+1)^2}{30}, & k \ge 8. \end{cases}$$
   Whenever $\gamma_{p,n} \ge \gamma_{\text{spiky}}$ on a leaf box at level $n$, the Gaussian's action on smooth polynomial scaling coefficients reduces to the truncated asymptotic differential moment series $\bm{M}(\gamma_{p,n})^{\otimes 3}$.
   - **Why `select_operator_gamma_spiky_override(k)` floors $k < 8$ operators at $270.0$ ($k = 9$ stiffness, keeping window width $W = 5$)**: At $k = 6$ (loose SCF protocol `1e-4`), the interior Taylor remainder formula gives $\gamma_{\text{spiky}}^{\text{interior}}(6) = 58.8$, at which the Gaussian boundary flux $\frac{1}{2}\operatorname{erfc}(\sqrt{\gamma}\,d)$ escaping across a leaf face from the outermost 6-point Gauss–Legendre node ($d = 0.0338$) is still **$35.7\%$**. Because Stage 2 ($\gamma \ge \gamma_{\text{spiky}}$) is strictly box-local while Stage 4 ($\gamma < \gamma_{\text{spiky}}$) is exact neighbor-coupled two-scale Galerkin, whenever a leaf box of an HF exchange pair `ij` or `Vpsi` splits or merges across `truncate_tol(1e-4)` during SCF iterations, Gaussians in $\gamma \in [58.8, 235]$ would otherwise switch between exact neighbor-coupled $\bm{T}$ and leaky box-local $\bm{M}(\gamma)$—creating a $1.215 \times 10^{-4}$ two-scale discontinuity ($\|V_c - V_f\|_2 = 6.12 \times 10^{-4}$ on bitwise-identical functions $\|f_c - f_f\|_2 = 2.3 \times 10^{-16}$) that KAIN subspace extrapolation ($\sim 14\times$ step gain) amplifies into a $\sim 1.2\text{--}1.5 \times 10^{-3}$ BSH residual plateau on large 3D molecules (`water20`, Valinomycin). Flooring `cfg.gamma_spiky_override = 270.0` for $k < 8$ suppresses boundary leakage below $0.5\%$ and reduces two-scale BSH non-invariance $\|2 G_\mu(\Delta K)\|_2$ by **$12.4\times$** ($1.215 \times 10^{-4} \to 9.76 \times 10^{-6}$), while keeping `cfg.gamma_spiky_override == 0.0` default `WindowedQuadraturePlan` behavior untouched.

2. **Diffuse Threshold $\gamma_{\text{diff}}(k, \epsilon)$ (`WindowedQuadraturePlan::compute_gamma_diff`)**:
   Defined as the unique root of the exact 3D spherical-Bessel / Plancherel wavelet projection bound:
   $$F_k(\gamma) = 2 \sup_{z > 0} \left[ e^{-z^2 / \gamma} \sqrt{\sum_{n=k}^{\infty} (2n + 1) j_n(z)^2} \right] = \epsilon.$$
   - The tail ratio $S_k(z) = \sum_{n=k}^{\infty} (2n+1) j_n(z)^2 / \sum_{n=0}^{\infty} (2n+1) j_n(z)^2$ (where the denominator equals $1$ by Parseval's identity) is evaluated in `spherical_bessel_tail_sum(z, k)` via **Miller's backward recurrence** starting at $n_{\text{top}} = \max(k+30, \lceil z + 30 + 5 z^{1/3}\rceil)$ with $10^{-100}$ dynamic overflow rescaling, avoiding catastrophic cancellation as $z \to 0^+$.
   - `plancherel_bound_value(gamma, k)` brackets the global maximum around $z_{\text{peak}} \approx \sqrt{k\gamma/2}$ on a 40-interval log-grid $[ \ln z_{\text{peak}} - 3, \ln z_{\text{peak}} + 2.5 ]$ followed by 60 golden-section refinement steps.
   - `compute_gamma_diff(k, eps)` seeds from the closed-form double-factorial asymptotic root and runs 64 bisection steps in $\ln \gamma$, memoizing `(k, eps) -> gamma_diff` (up to `256` entries) in a static thread-safe cache.

3. **Loose-Threshold Screening Tightening (`tighten_loose_windowed_eps`) & Apply-Time Tolerance Flooring**:
   In loose multi-center SCF protocols (`thresh >= 1e-4` with `truncate_mode == 1` on wide boxes $L \gtrsim 100$ Bohr), `truncate_tol` does not scale down on coarse octaves. Whenever `eps > 1e-5` (and no explicit `eps_loc > 0` override is set), `tighten_loose_windowed_eps(eps)` tightens internal quadrature/stencil/screening tolerance smoothly to $\min(5\times 10^{-5},\, \max(10^{-5},\, 0.1\epsilon))$, leaving `eps <= 1e-5` (`1e-6`, `1e-8`, ...) untouched. At apply time (`WindowedOrbitalExecContext`), when an operator was built with much tighter accuracy than the input function (`plan.eps() << eff_thresh`, e.g., SCF `coulop` built at `0.001 * thresh`), internal stencil screening floors `min_plan_tol = max(plan.eps(), 0.1 * eff_thresh)` before applying `tighten_loose_windowed_eps` (unless `eps_loc > 0` is explicitly specified).

4. **Remaining Fine-Scale Operator Weight Bounds (`rem_weight_bounds_`)**:
   `WindowedQuadraturePlan::compute_rem_weight_bound(j)` precomputes the exact $L_1(\mathbb{R}^3)$ (Young $L_2 \to L_2$) operator norm bound for all windowed + compressed spiky terms still active on wavelet details at levels $\ge j$:
   $$W_{\text{rem}}(j) = \sum_{p :\, n_{\text{spiky}}(p) < j \le n_{\text{diff}}(p)} \lvert W_p^{\text{scaled}}\rvert + \sum_{q=1}^{Q_{\text{eff}}(j)} \lvert \widetilde{W}_{j,q}^{\text{scaled}}\rvert$$
   for $j \in [0, n_{\max} + 12]$ (`rem_weight_bound(level)`). Because terms with $n_{\text{diff}}(p) < j$ have already vanished on wavelets at levels $\ge j$, $W_{\text{rem}}(j)$ decays as $O(4^{-j} L^2)$ on deep cusp levels and enables rigorous $H^1$-safe compressed-subtree screening in Stage 1a.

For each quadrature term $p$, the transition levels are:
$$n_{\text{spiky}}(p) = \left\lfloor \frac{2t_p - \ln\gamma_{\text{spiky}}}{2\ln 2} \right\rfloor, \qquad n_{\text{diff}}(p) = \left\lfloor \frac{2t_p - \ln\gamma_{\text{diff}}(k,\epsilon)}{2\ln 2} \right\rfloor,$$
giving the active NS window $j \in [j_0(p), j_1(p)]$ with $j_0(p) = \max(0, n_{\text{spiky}}(p) + 1)$ and $j_1(p) = \min(n_{\max} - 1, n_{\text{diff}}(p))$ (note Stage 4(a) base injection at $j_0(p) = n_{\max}$ is also included when $n_{\text{spiky}}(p) + 1 = n_{\max} \le n_{\text{diff}}(p)$). For $j > n_{\text{diff}}(p)$, wavelet outputs vanish below $\epsilon$ and scaling outputs follow automatically from two-scale reconstruction of coarser levels.

---

### 3.2 Dyadic-Locked Quadrature, Tail Completion & Overflow-Free GC Compression (`WindowedQuadraturePlan`)

1. **Dyadic Scale-Locking (`select_oversampling_M`)**:
   Using the exponential substitution $t = \ln s$ in $\frac{e^{-\mu r}}{r} = \frac{2}{\sqrt{\pi}} \int_{-\infty}^{\infty} \exp\!\big(-e^{2t} r^2 - \frac{\mu^2}{4} e^{-2t} + t\big)\,dt$ with step size
   $$\Delta t = \frac{\ln 2}{M}, \qquad \beta = e^{2\Delta t} = 4^{1/M},$$
   the smallest integer $M \in [2, 8]$ is chosen so Proposition 4.1's trapezoidal discretization error $2/\sqrt{\cosh(\pi^2 / \Delta t)} \le 0.1\epsilon$ ($M = 2$ for $\epsilon \ge 2\times 10^{-5}$, $M = 3$ for $\epsilon \in [10^{-8}, 10^{-5}]$, $M = 4\text{--}5$ for tighter $\epsilon$). Decomposing $p = q_{\text{oct}} M + m$ with Euclidean residue $m \in \{0, \dots, M-1\}$ (`divmod_euclidean`), we precompute exact integer base offsets:
   $$s_{\text{base,spiky}}[m] = \left\lfloor \frac{m}{M} - \frac{\ln\gamma_{\text{spiky}}}{2\ln 2} \right\rfloor, \qquad s_{\text{base,diff}}[m] = \left\lfloor \frac{m}{M} - \frac{\ln\gamma_{\text{diff}}}{2\ln 2} \right\rfloor.$$
   Then $n_{\text{spiky}}(p) = q_{\text{oct}} + s_{\text{base,spiky}}[m]$ and $n_{\text{diff}}(p) = q_{\text{oct}} + s_{\text{base,diff}}[m]$ hold with **zero floating-point floor drift** across all octaves $q_{\text{oct}} \in \mathbb{Z}$, and at every active level $j \in [j_0(p), j_1(p)]$:
   $$\gamma_{p,j} = \gamma_{\text{ref}}(w_{\text{eff}}, m) = \exp\!\Big(2\ln 2 \big(\tfrac{m}{M} - (s_{\text{base,spiky}}[m] + 1 + w_{\text{eff}})\big)\Big), \qquad w_{\text{eff}} = j - (n_{\text{spiky}}(p) + 1) \in [0, W - 1],$$
   where $W = \max_{0 \le m < M}(s_{\text{base,diff}}[m] - s_{\text{base,spiky}}[m])$. Off-dyadic-grid terms (`on_dyadic_grid == false`, i.e., the $\mu = 0$ synthetic low-tail Gaussian and `CustomGaussian` user terms) precompute exact per-`(term_idx, level)` `StencilEntry1D` records in `WindowedStencilCache1D::offgrid_entries_` while sharing the identical two-threshold and GC machinery.

2. **Discrete Low-Tail Completion at $\mu = 0$**:
   On $[0, 1]^3$ ($r_{\max} = \sqrt{3}$), truncating the $t_p = p\Delta t$ sum cleanly at $p_{\min} = \lfloor \frac{1}{7\Delta t} \ln(\frac{525\sqrt{\pi}\epsilon_{\text{quad}}}{4 r_{\max}^7}) \rfloor$ ($\epsilon_{\text{quad}} = 0.1\epsilon$) leaves an infinite geometric series over $p < p_{\min}$ whose exact discrete moments in $r^2$ are:
   $$M_r = \frac{2\Delta t}{\sqrt{\pi}} \frac{e^{(2r+1) t_{p_{\min}}}}{\operatorname{expm1}((2r+1)\Delta t)}, \qquad r \in \{0, 1, 2\}.$$
   Matching $M_0, M_1, M_2$ exactly via a single off-grid Gaussian $\alpha_{\text{low}} = M_2 / M_1$, $c_{\text{low}} = M_1^2 / M_2$ plus a constant kernel shift $w_{\text{const,sim}} = M_0 - c_{\text{low}}$ (`w_const_scaled = op_scale * w_const_sim`) reproduces $1/r$ across the entire domain $[10^{-8}, \sqrt{3}]$ to $O(\epsilon_{\text{quad}})$ without any continuous-vs-discrete Euler–Maclaurin mismatch. For $\mu_{\text{sim}} > 0$, super-exponential decay $\exp(-\frac{1}{4}\mu_{\text{sim}}^2 e^{-2t_p})$ allows direct truncation at $p_{\min}$ with $w_{\text{const}} = 0$.

3. **High-$\alpha$ Cusp Completion & Level-Rescaled Stieltjes–Golub–Welsch Compression (`build_gc_rule_for_level`)**:
   Explicit dyadic terms are generated up to $p_{\text{tail}} = p_{\max} + 10M$, where $p_{\max} + 1$ is the first grid index spiky at $n_{\max}$. At each tree level $n \ge 0$ (precomputed lock-free for $0 \le n \le n_{\max} + 12$ and memoized on demand beyond):
   - For `Coulomb` / `BSH` / `YukawaRaw`, spiky terms ($n_{\text{spiky}}(p) \ge n$) have inverse stiffness $u_{p,n} = \gamma_{p,n}^{-1} = \exp(2n\ln 2 - 2t_p) \le \gamma_{\text{spiky}}^{-1}$ and level-normalized positive weights $\widehat{W}_{p,n} = \text{op\_scale}\cdot (2\pi\Delta t)\, u_{p,n}\, e^{-\mu_{\text{arg}}}$, completely free of floating-point underflow at arbitrary cusp depth $n$. When $\mu_{\text{sim}} = 0$ and $n_{\text{spiky}}(p_{\text{tail}}) \ge n$, the exact infinite geometric tail $p > p_{\text{tail}}$ is appended as a single remainder spike:
     $$\widehat{W}_{\text{rem},n} = \frac{\widehat{W}_{p_{\text{tail}},n}}{\beta - 1}, \qquad u_{\text{rem},n} = \frac{u_{p_{\text{tail}},n}}{\beta + 1}.$$
   - For `CustomGaussian`, if any spiky term has negative weight or at most $Q$ spiky terms exist at level $n$, exact uncompressed `GCNode`s are emitted directly (with $2t_p - 2n\ln 2$ clamped at $690.0$ so ultra-spiky $\gamma_{p,n}$ never overflows `double`); when $> Q$ positive spiky terms exist, they are compressed via the same Stieltjes–Golub–Welsch solver on $u_{p,n} = \exp(\max(-690.0, 2n\ln 2 - 2t_p))$.
   - Rescaling $s_i = u_i / u_{\max} \in (0, 1]$ and $w_i = \widehat{W}_i / \sum_l \widehat{W}_l$, discrete Stieltjes orthogonalization builds the $Q_{\text{eff}} \times Q_{\text{eff}}$ symmetric tridiagonal Jacobi matrix $\bm{J}$ (`select_gc_rank_Q`: $Q = 4$ for $k \le 6$, $Q = 5$ for $7 \le k \le 14$, $\max(5, \lfloor(m_k + 4)/2\rfloor)$ for $k > 14$), diagonalized via `madness::syev`. Eigenvalues $s_q$ and top eigenvector components $V_{0,q}^2$ yield $Q_{\text{eff}} \le Q$ compressed nodes $\tilde{\gamma}_{n,q} = (u_{\max} s_q)^{-1}$ and level-$n$ weights $\widetilde{W}_{n,q}^{\text{scaled}} = e^{-2n\ln 2}\widehat{W}_{n,q}$, each storing its precomputed $k \times k$ lower-triangular matrix $\bm{M}(\tilde{\gamma}_{n,q})^{\mathsf{T}}$.

4. **Parity-Sparse 1D Local Moment Operator (`LegendreMomentOp`)**:
   On orthonormalized Legendre polynomials $\tilde{\phi}_i(x) = \sqrt{2i+1}\,P_i(2x-1)$ on $[0, 1]$, the exact second-derivative operator has entries
   $$D_2(i, q) = \frac{1}{2}\sqrt{(2i+1)(2q+1)}\,(q - i)(i + q + 1) \quad \text{for } q > i,\; (q - i) \text{ even},$$
   and zero otherwise. `LegendreMomentOp` precomputes $\frac{1}{r!}\bm{D}_2^r$ for $0 \le r \le \lfloor(k-1)/2\rfloor$ and forms
   $$\bm{M}(\gamma) = \sum_{r=0}^{\lfloor(k-1)/2\rfloor} \gamma^{-r} \left(\frac{1}{r!}\bm{D}_2^r\right),$$
   storing its transpose `MT(in_deg, out_deg)` (`fill_MT_buffer` / `make_MT`) so 3D action $\bm{M}(\gamma)^{\otimes 3} s$ maps directly to a single allocation-free `fast_transform(s, MT, tmp, work)` call.

---

### 3.3 Universal 1D Stencil Cache & Exact 3D Sub-Block Norms (`windowed_ns_stencil.h`)

For dyadic-locked `Coulomb`, `BSH`, and `YukawaRaw` plans, `WindowedStencilCache1D` delegates to a process-wide thread-safe table (`get_or_build_universal_stencil_entries`, retaining up to `32` unique tables alongside `windowed_convolution.cc`'s `64`-entry plan+stencil cache) keyed solely by $(k, W, M, \epsilon, \gamma_{\text{ref}}(0,0), \gamma_{\text{ref}}(W-1,M-1))$, while precomputing `root_T0_by_term_` for ultra-diffuse terms ($n_{\text{diff}}(p) < 0$) and `offgrid_entries_` for off-grid terms. Consequently:
- Across all $40\text{--}300+$ `BSHOperatorPtr3D(mu_i)` constructions per molecular SCF run (where orbital energies $\varepsilon_i$ and $\mu_i = \sqrt{-2\varepsilon_i}$ change every iteration), **zero** 1D quadrature builds and **zero** LAPACK `dgesvd` calls occur after the initial $W \times M$ table construction (`< 0.5 ms`).
- Each `StencilEntry1D` queries `GaussianConvolution1DCache<double>::get(k, gamma_ref, 0, LatticeRange())->nonstandard(0, disp)` for $\text{disp} = l_{\text{target}} - l_{\text{source}} \in [-2, +2]$ and stores:
  - Contiguous `(dim_in, dim_out)` matrices oriented directly for `mTxmq`: `T` $(k \times k)$, `R` $(2k \times 2k)$, `R_2k_to_k` $(2k \times k)$, and `R_k_to_2k` $(k \times 2k)$, plus static ceiling half-bandwidths `m_max_T` and `m_max_ns` (`base_tol = 1e-4 * eps`).
  - Exact 1D SVDs (`StencilSVD1D`) for all four block shapes, including **rank-major contiguous prefix matrices `U_by_rank[r]`** of shape `(d_in, r)` (`ldb == r` with zero column stride padding!) and exact Frobenius tail norms `tail_normf[r]` $= (\sum_{m \ge r} \sigma_m^2)^{1/2}$.
  - Exact 3D Frobenius norms via orthogonal block subtraction—e.g. for Pass 4b-1 (`(0, d) -> s`):
    $$\|\bm{R}_{2k \to k}(\bm{d})\big|_{(0,d)\to s}\|_F = \sqrt{\prod_{\alpha \in \{x,y,z\}} \|R_{2k\to k}(d_\alpha)\|_F^2 - \prod_{\alpha \in \{x,y,z\}} \|R_{ss}(d_\alpha)\|_F^2},$$
    analogously for `opnorm_T`, `opnorm_4b2_leaf` (`(s, 0) -> d`), `opnorm_4b2_int` (`(s, d) -> d`), and `opnorm_4b2_wav` (`(0, d) -> d`), aggregated into per-shell maxima `shell_max_*[sh]` and shell $L_2$ norms `shell_l2_*[sh]` for $\text{sh} = \|\bm{d}\|_\infty \in \{0, 1, 2\}$.
  - Exact 3D constant-scaling-field wavelet residual norm `const_wav_3d_norm[m_bw]` $= \|\bm{Q}^j (\sum_{\|\bm{d}\|_\infty \le m_{\text{bw}}} \bm{R}(\bm{d})^{\otimes 3} \bm{e}_{000})\|_F$ and Theorem-3.1 constant-field commutation remainder norm `const_T_rem_3d_norm[m_bw]` $= \|(\sum_{\|\bm{d}\|_\infty \le m_{\text{bw}}} \bm{T}(\bm{d})^{\otimes 3} - (\bm{M}(\gamma)^{\mathsf{T}})^{\otimes 3})\bm{e}_{000}\|_F$ for $m_{\text{bw}} \in \{0, 1, 2\}$.

---

### 3.4 Five-Stage Execution Pipeline, Subtree Pruning & Multi-Rank Lockstep Execution (`windowed_ns_halo.h`, `windowed_ns_sumfact.h`, `windowed_convolution.cc`)

Given one or a batch of input `Function<T, 3>` objects ($T \in \{\text{double}, \text{double\_complex}\}$):

#### Pure-Spiky (`pure-GC`) Fast Path (`apply_pure_spiky_same_world`)
Whenever every term in `plan` is spiky across the entire input tree ($n_{\text{spiky}}(p) \ge d_{\max}$, so there are zero active $\bm{T}^j$ or $\bm{\Gamma}^j$ terms and $w_{\text{const}} = 0$), `apply_pure_spiky_same_world` bypasses `make_nonstandard`, `WindowedNSHaloWorkspace`, `virt_cache`, and `result.reconstruct()` altogether: it evaluates the compressed leaf `GCLevelRule` (`fast_transform` + `aligned_axpy`) across leaf scaling cubes in parallel chunks and writes the result directly in `TreeState::reconstructed`.

#### Stage 1 — $H^1$-Pruned Active-Window Ingestion & Fine-Driven Virtual Halo (`WindowedNSHaloWorkspace`)
1. **Stage 1a (`WindowedOrbitalExecContext`, `WindowedSerialInputRange`) & $H^1$-Safe Compressed Union-Skeleton Pruning**:
   - **$H^1$-Safe Compressed Subtree Screening (`subtree_d_tol`)**: In multi-center SCF (`SCF::compute_residual`), `Vpsi[i] -= transform(psi, fock, trantol)` is formed in `TreeState::compressed` without a subsequent `.truncate()`, causing every localized `Vpsi[i]` to inherit the deep nuclear-cusp octree union of *all* $N_{\text{MO}}$ orbitals (`~10^4` union nodes/MO vs `~10^3` active leaves). When `dynamic_refine == true`, Stage 1a defines the $H^1$-safe geometric-mean operator weight $W_{H^1}(j) = \sqrt{W_{\text{rem}}(0) \cdot W_{\text{rem}}(j)}$ and detail-norm threshold
     $$\tau_{\text{sub}}(j) = \frac{0.002 \cdot \text{truncate\_tol}(\text{eff\_tol}, \text{Key}(j))}{W_{H^1}(j)} \qquad (j \ge 2),$$
     preserving own-atom nuclear cusps in the kinetic-energy ($H^1$) norm while pruning negligible distant-atom union cusps ($j \ge 2$, never pruning $j \le 1$) whose cumulative wavelet detail norm squared $\sum_{U \subset T_v} \|d_U\|_F^2 \le \tau_{\text{sub}}(j)^2$.
   - **Single-Rank (`world.size() == 1`)**: Inputs in `TreeState::reconstructed`, `TreeState::compressed`, or `TreeState::nonstandard_with_leaves` are ingested **read-only with zero `WorldContainer` mutations**:
     - `ingest_local_real_nodes_serial_reconstructed`: extracts leaf scaling cubes directly from `coeffs`, restricts internal scaling cubes bottom-up strictly over `[min_s_int_level, n_max_tree - 1]` (via $k^3$ `restrict_step_accumulate_3d` outside the $\Gamma$ window), and synthesizes `(2k)^3` `ns_in` tensors via `fast_transform(box_2k, cdata_.hgT, rn.ns_in, work_2k)` + `wavelet_octants_normf` strictly inside `[min_gamma_level, max_gamma_level]`.
     - `ingest_local_real_nodes_serial_compressed`: runs a bottom-up compressed subtree detail-norm pass (converting pruned `j >= 2` subtrees into leaves at their coarse scaling sum while safely preserving sparse `vtransform` internal nodes whose zero children or empty ancestors are absent from `coeffs`), then expands unpruned octree children top-down in parallel (`fast_transform(..., cdata_.hg)` / `prolong_step_3d`) and populates active-window internal scaling/NS tensors.
   - **Distributed MPI (`world.size() > 1`)**: For compressed inputs (`TreeState::compressed`) with `enable_dynamic_refinement`, `apply_windowed_batch_impl` counts internal vs. seed ($j \le 1$ or $\|d\|_F > \tau_{\text{sub}}(j)$) nodes in a single `gop.sum`, activates pre-pruning whenever `orig_int >= 32` and `seed_int <= 0.8 * orig_int` ($\ge 20\%$ union-skeleton reduction), and builds slim `TreeState::compressed` working copies across all `A_c` compressed orbitals via a coalesced ancestor notification (`prune_ex`, `1` send fence + `1` completion fence) *before* batched `reconstruct(false)` + `make_nonstandard(true, false)`, leaving caller inputs unmutated.
2. **Global Level-Norm Reduction & Coarse-Leaf Cauchy–Schwarz Bounds (`WindowedGlobalLevelStats`)**:
   - Per-node scaling cubes `s_in`, average coefficient `s000 = s_in(0,0,0)`, cancellation-free oscillatory norm squared `s_osc_norm2` $= \sum_{(i,l,m)\ne(0,0,0)} \lvert s_{\text{in}}(i,l,m)\rvert^2$ (`cblas::dot`), Legendre Cauchy–Schwarz leaf $L_\infty$ bound `max_s_leaf_linf` $= \lvert s_{000}\rvert + k^3 \|s_{\text{osc}}\|_F$, and internal pure wavelet norms `wavelet_octants_normf(ns_in, k)` (three disjoint slices `(s1,all,all)`, `(s0,s1,all)`, `(s0,s0,s1)` with zero `(2k)^3` copies) are reduced globally as **6 doubles per level** (`has_internal`, `max_s_real`, `max_s_leaf`, `max_s_leaf_linf`, `max_d_int`, `max_ns_int`) in a single coalesced `gop.max` across all $B$ orbitals.
   - `compute_coarse_leaf_above_bound(j)` then tightens the prolonged virtual scaling norm bound at level $j$:
     $$\text{max\_coarse\_leaf\_above}[j] = \max_{0 \le b < j} \min\!\Big(\text{max\_s\_leaf}[b],\; 2^{-1.5(j-b)} \,\text{max\_s\_leaf\_linf}[b]\Big), \quad \text{max\_s\_source}[j] = \max\big(\text{max\_s\_real}[j],\, \text{max\_coarse\_leaf\_above}[j]\big).$$
3. **Active Bandwidth & Per-Node Operator-Weighted Shell-1 / Shell-2 Gating (`make_radii_builder`, `WindowedLevelRadii`)**:
   - Using `loc_tol_j = f_impl->truncate_tol(eff_tol, Key<3>(j))`, `tol_nn = 0.01 * loc_tol_j / sqrt(N_raw)`, and `tol_sh2 = loc_tol_j / (14.14 * sqrt(N_raw))` (`0.1 * tol_nn` when `dynamic_refine == false`), each candidate term in `raw_T_terms[j]` and `raw_gamma_terms[j]` selects active half-bandwidths `bw_T` $\in \{-1, 1, 2\}$, `bw_4b1` $\in \{-1, 0, 1, 2\}$, `bw_leaf` $\in \{-1, 1, 2\}$ (**never `0`**, preserving `[-1,+1]^3` neighbor polynomial cancellation), and `bw_int >= bw_leaf` $\in \{-1, 0, 1, 2\}$.
   - **Per-Node Operator-Weighted Shell-1 (`sh == 1`) & Shell-2 (`sh == 2`) Halo Gating**: To prevent neighbor-nucleus $-Z/r \cdot \psi_{\text{tail}}$ cusps on multi-center molecules (`Vpsi` in `Apply BSH`) from inheriting global cusp-driven halo radii, `WindowedLevelRadii` gates halo emission per node and per coarse-leaf bound:
     - `node_mT_and_l1`: checks `s_normf > sh2_s_thresh` for `m_max_T == 2`, else requires both the Theorem-3.1 constant-subtracted bound $\lvert s_{000}\rvert\,\text{sh1\_T\_const\_op} + \|s_{\text{osc}}\|_F\,\text{sh1\_T\_osc\_op} > \tau_{T,1}$ ($\tau_{T,1} = \text{loc\_tol\_j}/28.28$) and the Manhattan distance `l1 = |dx|+|dy|+|dz| in {1,2,3}` face/edge/corner bound `s_normf * sh1_T_l1_op[l1] > sh1_T_tol`.
     - `node_m4b1` & `node_mRefine`: gate Pass 4b-1 `kHaloH1Gamma` (`d_normf > sh2_4b1_thresh` / `sh1_4b1_thresh`) and Pass 4b-2 `kHaloH1Refine` (`min(s_normf * sh1_refine_s_op, |s000| * sh1_refine_const_op + ||s_osc||_F * sh1_refine_osc_op) + d_normf * sh1_refine_d_op > 0.02 * refine_tol_j`).
     - `h2_mT` & `h2_mLeaf`: gate secondary `kHaloH2Ghost` completion radii (`{-1, 1, 2}`) using `max_coarse_leaf_above` against `sh1_s_T_thresh`, `sh1_s_leaf_thresh`, and `sh2_s_thresh`, cutting `water9` `Apply BSH` virtual halo boxes per MO by **`3.22x` (`17,560` $\rightarrow$ `5,457`)** with 100% MPI rank invariance.
4. **Fine-Driven Halo Discovery, Graded `coulop` Refinement & Lockstep Prolongation (`virt_cache`)**:
   - Real nodes emit neighbor candidates tagged with bitmask flags `kHaloH1T` (Stage 4(a) $\bm{T}^{j_0}$ + Theorem 3.1 subtraction), `kHaloH1Gamma` (Pass 4b-1 `ds_corr`), `kHaloH1Refine` / `kHaloH1RefineLeafUngraded` (Pass 4b-2 wavelet promotion), and `kHaloH2Ghost` (read-only stencil completion).
   - **Ungraded (`dj >= 2`) Leaf-Boundary Gating & High-Precision Operator Override (`ungraded_leaf_refine_only`)**: Each internal node records an 8-bit `internal_child_mask` indicating which of its 8 children at $j+1$ are themselves internal ($\Delta j \ge 2$ jump). `ungraded_leaf_refine_only` is set to `(eff_thresh >= 1e-6) && !((k < 8) && (plan.eps() <= 0.01 * eff_thresh))`: on `Exchange` and `BSH` (`plan.eps() ~ eff_thresh`), Pass 4b-2 leaf promotion probes run only across genuine ungraded ($\Delta j \ge 2$) coarse-to-fine jumps, whereas when a high-precision Coulomb operator (`coulop` with `plan.eps() = 0.001 * thresh = 1e-7`) acts on a loose-threshold density (`k = 6, eff_thresh = 1e-4`), full graded `sh == 1` leaf refinement remains enabled so long-range $1/r$ potentials maintain full accuracy.
   - **Serial (`world.size() == 1`)**: `build_serial` runs two flattened lock-free parallel passes (`H1` discovery + `near_internal_4b1/4b2` marking, followed by `H2_Ghost` completion) over `kHaloChunkNodes = 16`-node slices using 64-bit packed coordinates `(lx << 42) | (ly << 21) | lz` at $j \le 20$.
   - **Lockstep Distributed MPI (`world.size() > 1`, `build_distributed_lockstep`)**: Across an entire tile of $B \le 16$ orbitals simultaneously, a single `WindowedNSHaloExchange<T>(world, B)` instance tags every message with `std::uint16_t slot` $\in [0, B-1]$ and executes **1 coalesced `gop.max` + 3–4 `world.gop.fence()` collectives for the entire $B$-orbital batch**: (1) coalesced 6-stat `gop.max`; (2) optional parent-mask exchange (`1` fence when `ungraded_leaf_refine_only` is active; skipped otherwise); (3) route `H1` `(slot, Key<3>, flags)` candidates (`1` fence); (4) walk parent chains per slot and send deduplicated `KeyQueryMsg` batches (`1` fence); (5) reply with read-only real leaf `s_in` and real internal `ns_in` tensors (`1` fence).
   - **Single-Pass Top-Down Prolongation (`populate_prolonged_virtual_and_sources`)**: Inserts intermediate virtual parents bottom-up ($j = n_{\max} \dots 2$) and prolongs each virtual box once top-down ($j = 1 \dots n_{\max}$) via `general_fast_transform` using octant filters `h_oct_[child_octant][axis]`, merging sorted `owned_real` and virtual entries in $O(N)$ and eagerly freeing `prolong_levels[sj]` per level.

#### Stage 2 — Cross-Level Flattened Leaf Gauss–Christoffel Pass
All owned real leaves ($B \in \mathcal{L}$ with $\|s_B\|_F > 0$) across $j = 0 \dots n_{\max}$ are flattened into a single vector `gc_leaves` and processed via `parallel_for_chunks`:
$$\hat{s}_B \mathrel{+}= \sum_{q=1}^{Q_{\text{eff}}(n)} \widetilde{W}_{n,q}^{\text{scaled}}\, \bm{M}(\tilde{\gamma}_{n,q})^{\otimes 3} s_B, \qquad n = \text{level}(B),$$
using per-slot L1-resident $k^3$ scratch cubes (`fast_transform` + `aligned_axpy`, zero heap allocations per leaf).

#### Stage 3 — Root ($j = 0$) Constant Low-Tail & Ultra-Diffuse Pass
On the rank owning root key `(0, {0,0,0})`:
1. If $w_{\text{const}}^{\text{scaled}} \ne 0$ ($\mu = 0$), adds $w_{\text{const}}^{\text{scaled}} \, s_{\text{in}}(0,0,0)$ directly to the $(0,0,0)$ Legendre constant mode `root_rn.s_out(0,0,0)`.
2. Applies $W_p^{\text{scaled}} (\bm{T}^0_{p,00})^{\otimes 3} s_0$ via `apply_single_box_1d_kronecker_3d` (`fast_transform` + `aligned_axpy` using `root_T0_transposed`) for every ultra-diffuse term with $n_{\text{diff}}(p) < 0$.

#### Fused Per-Level Execution of Stages 4(a) & 4(b) (`LevelInputReleaseGuard`)
Stages 4(a) and 4(b) execute in a **single fused level loop** `for (int j = 0; j <= n_max_tree; ++j)` inside `WindowedOrbitalExecContext::execute_local_stages_2_to_4` (100% rank-local, zero MPI communication!) so that as soon as level $j$ finishes its scaling and wavelet sweeps, `LevelInputReleaseGuard` invokes `ws.release_stage4_inputs_for_level(j)`—immediately freeing `ld.sources`, `ld.source_packed_keys`, `rn.s_in`, `rn.ns_in`, and `vn.s_virt` at level $j$ *before* level $j+1$ allocates any `(2k)^3` output cubes.

#### Stage 4(a) — Base Scaling Injection $\bm{T}^{j_0}$, Fused Multi-Term Batching & Fused Theorem 3.1 Halo Subtraction
At $j = 0$, any non-empty `active_T_terms[0]` applies `mat_T(0)` on the root box via `apply_single_box_1d_kronecker_3d`. For each level $j \in [1, n_{\max}]$ with non-empty `active_T_terms[j]`:
1. **3-Tier Source-Norm Bucketing, Bandwidth Merging & Fused Term Batching**: Active scaling sources passing `s_normf * max_T_single_op > 0.005 * loc_tol_j / 27` are partitioned into up to 3 disjoint sorted norm buckets (`> 0.02 M`, `> 1e-4 M`, `<= 1e-4 M`) and merged via `merge_matching_norm_buckets` whenever active bandwidths agree. Within each spatial target brick (`1,024` boxes, `run_sumfact_target_tiled_3d`), active terms sharing the same `select_active_1d_slots` half-bandwidth `eff_bw` $\in \{1, 2\}$ (`svd_fac_T = 0.025`, `min_active_bw = 1`) are grouped and executed together in **`execute_sumfact_batch_3d`** (`kMaxFusedTerms = 8`, `4 MiB` `buf_x + buf_xy` cap), sharing a single `SumFactGraph3D` traversal and accumulating all `n_sub` Z-contractions directly into `tgt_out_ptrs[ti]` in L1 cache.
2. **Fused Theorem 3.1 Boundary-Halo Commutation Subtraction**: All active base terms at level $j$ are subtracted in a **single parallel dispatch** over `h1t_virts` ($V \in \partial\mathcal{H}_j^{(1,T)}$, strictly outside source-norm buckets), reusing each virtual box's L1-resident `s_virt` cube across every active $(p, j)$ term:
   $$\Delta\hat{s}_{\text{corr}}[V] \mathrel{-}= \sum_{p \in \mathcal{A}_T(j)} W_p^{\text{scaled}}\, \bm{M}(\gamma_{p,j})^{\otimes 3} s_{\text{virt}}[V].$$

#### Stage 4(b) — Non-Standard Wavelet Decoupling $\bm{\Gamma}^j$ & Dynamic Output Refinement
Executed inside the fused level loop only at levels $j \in [0, n_{\max}-1]$ where `ld.global_has_internal == true` and `active_gamma_terms[j]` is non-empty:
1. **Pass 4b-1 (`2k -> k`, `(0, d^j) -> \hat{s}^j`, Cancellation-Free)**:
   - Level-scoped pure wavelet cubes `p1_wav_cubes` (`make_pure_wavelet_cube`) are materialized for level $j$'s screened internal sources (`d_normf * max_4b1_single_op > 0.005 * loc_tol_j / 27`) and freed immediately when Pass 4b-1's scope exits—before Pass 4b-2 allocates `(2k)^3` `ns_out` tensors.
   - Self-box terms (`eff_bw == 0`) are fused across all active terms in a single parallel pass via **`execute_self_box_2k_to_k_batch_3d`**; neighbor terms (`eff_bw in {1, 2}`) run grouped by bandwidth (`512`-target bricks) via **`execute_sumfact_batch_3d`**.
2. **Constant-Subtracted Refinement Candidate Pre-Screening (`candidate_exceeds_tol`)**:
   - Separates Shell 1 (`[-1,+1]^3`, 27 boxes) from Shell 2 ($\text{sh} = 2$) and exploits constant polynomial cancellation around $c_0 = s_{000}(U)$ (`|c0| * eta_const_sh1 + ub_sh1_sub`) whenever all 27 `sh <= 1` slots exist in `ld.sources`, while adding rigorous worst-case bounds (`missing_sh1_s_bound`, `missing_sh2_s_bound`) whenever `h2_mLeaf` clamped Shell-1 or Shell-2 `H2_Ghost` sources, rejecting $>90\%$ of smooth valence jump candidates with 100% MPI rank invariance before allocating `(2k)^3` output cubes.
3. **Pass 4b-2 (`2k -> 2k` Internal + `k -> 2k` Leaf/Virtual -> $\hat{d}^j$)**:
   - Self-box internal terms (`eff_bw == 0`) execute fused via **`execute_self_box_2k_to_2k_batch_3d`**; neighbor internal terms (`eff_bw in {1, 2}`) and leaf/virtual terms (`eff_bw in {1, 2}`, `min_active_bw = 1`, never `0`) execute grouped by bandwidth inside `256`-target bricks via **`execute_sumfact_batch_3d`**.
   - `zero_s0_octant` clears `[0..k)^3` on `ns_out` / `ns_virt_out`, and any boundary real leaf (`near_internal_4b2`) or `kHaloH1Refine` virtual node with `wavelet_octants_normf > truncate_tol(eff_thresh, key) * eff_refine_fac` sets `promoted_leaf = true` / `promoted = true`.
4. **Exact Shared-Tree Linearity Mode (`dynamic_refine == false`)**:
   - Replaces data-dependent level norms with `1.0`, tightens `svd_fac` from `0.025` to `0.003`, disables norm bucketing/source screening/subtree pruning/per-node halo gating, and skips output tree promotion so $\mathcal{G}$ is an exact linear map to machine precision on shared octrees.

#### Stage 5 — Direct Single-Rank Synthesis or Lockstep Multi-Rank Assembly
1. **Stage 5a (`assemble_output_reconstructed_serial` / `assemble_output_nonstandard_lockstep`)**:
   - **Direct Single-Rank Reconstructed Synthesis (`assemble_output_reconstructed_serial`)**: When `world.size() == 1` and `!any_promoted`, restricts non-zero halo deltas bottom-up (`restrict_all_virtual_unpromoted`), frees `owned_virt`, and runs top-down two-scale synthesis (`fast_transform(..., cdata_.hg)` / `prolong_step_3d`) directly into `result` (`TreeState::reconstructed`), falling back to `assemble_output_nonstandard_tree` + `result.reconstruct()` only when dynamic promotions occurred.
   - **Lockstep Multi-Rank Assembly (`assemble_output_nonstandard_lockstep`)**: Across all $B \le 16$ orbitals in a distributed tile simultaneously, executes **1 coalesced `gop.sum` + 1–2 `world.gop.fence()` collectives**: (1) `gop.sum(promo_counts)` across all $B$ orbitals; (2) bottom-up local virtual restriction + coalesced `CoarseLeafDeltaMsg` / `VirtualSubtreeItemMsg` transport (`1` fence); (3) deterministic `(leaf_key, sort_tag_key)` accumulation per slot + optional remote `OutputNodeInsertMsg` fence (**strictly skipped when no orbital in the tile had promotions**), leaving `result[b]` in `TreeState::nonstandard_after_apply` for asynchronous `result[b].reconstruct(false)`.
2. **Stage 5b (`restore_windowed_ns_input_3d`)**:
   - On `world.size() == 1`, non-destructive `reconstructed`/`compressed`/`nonstandard_with_leaves` inputs were never modified and require zero work; on `world.size() > 1`, slim pre-pruned compressed copies are cleared while non-compressed inputs restore their original `TreeState` asynchronously before the final tile fence.

---

## 4. High-Performance Kernel, Memory & Parallel Execution Optimizations

### 4.1 100% Tier-1 `mTxmq` (`libxsmm` JIT `s_dense_table`) Lowering & Fused Multi-Term Batching

Every 3D tensor contraction in Stages 1–5 maps to MADNESS optimized tensor kernels with zero scalar C++ tensor loops:

| Stage / Operation | Kernel Lowering |
| :--- | :--- |
| Stage 1b Prolongation & Stage 5a Restriction / Synthesis | `general_fast_transform` (3 cyclic `mTxmq` passes) + `aligned_axpy` |
| Pure-GC Fast Path, Stage 2 Leaf GC (`Q_eff <= Q`) & Fused Stage 4(a) `-M(gamma)` | `fast_transform` (3 cyclic `mTxmq` passes) + `aligned_axpy` (single pass over `h1t_virts`) |
| Stage 3 & 4(a) Root (`j = 0`) | `fast_transform` + `aligned_axpy` (`apply_single_box_1d_kronecker_3d`) |
| Stage 4(b) Fused Self-Box (`m_bw == 0`, `2k -> k` & `2k -> 2k`) | `execute_self_box_*_batch_3d`: 6-step Tier-1 `mTxmq` Tucker contraction `(2k)^3 -> r^3 -> d_out^3` using `U_by_rank[r]` (`ldb == r`), or `general_fast_transform`/`fast_transform` when dense is faster |
| Stage 4(a) & 4(b) Fused Spatial 3D Sweeps (`m_bw in {1, 2}`) | `execute_sumfact_batch_3d` (`<= 8` terms/sub-batch): Tier-1 `mTxmq` + L1 `fast_transpose` + direct L1 Z-accumulation |

Inside `execute_sumfact_batch_3d` and `execute_self_box_*_batch_3d` (`windowed_ns_sumfact.h`), three low-level properties drive middle-of-tree and deep-cusp throughput:
1. **Guaranteed Tier-1 `s_dense_table` Dispatch & Exact FLOP Break-Even**: MADNESS's `mTxmq(dimi, dimj, dimk, c, a, b)` dispatches to its pre-JIT'd `libxsmm` `s_dense_table` kernel when `dimi <= 400`, `dimj <= 24`, and `dimk <= 24`. In `apply_slot_edge_1d` (`make_slot_op_1d`), 1D SVD (`U_by_rank[r]` of contiguous `(d_in, r)` with `ldb == r`, plus `SVT`) is selected iff `r * (d_in + d_out) < d_in * d_out` (`r == 1`, eliding the transpose completely!) or `5 * r * (d_in + d_out) < 4 * d_in * d_out` (`r >= 2`, using L1-resident `fast_transpose`); similarly, self-box Tucker SVD is used iff `20 * flops_svd < 17 * flops_dense` (`2k -> k`) or `r < max(1, floor(0.65 * 2k))` (`2k -> 2k`).
2. **Multi-Term Bandwidth Fusion (`kMaxFusedTerms = 8`, `4 MiB` Pencil Cap)**: Active quadrature terms at a given level sharing the same active stencil half-bandwidth `eff_bw` (`2*eff_bw + 1` displacements per axis) traverse the **exact same** `SumFactGraph3D` (`x`, `xy`, `z` CSR tables). Sub-batching up to `8` terms (subject to `(need_x_1 + need_xy_1) * sizeof(T) * n_sub <= 4 MiB` so all intermediate pencils stay inside L3 cache) parallelizes Sweep 1 over `n_sub * num_x` and Sweep 2 over `n_sub * num_xy`, while Sweep 3 parallelizes over `num_active_tgt` (`z_grain = max(1, ceil(base_z_grain / n_sub))`) and loops `s = 0 .. n_sub - 1` inside each target box `ti`, accumulating directly into `tgt_out_ptrs[ti]` while keeping it pinned in L1 cache in exact deterministic term order.
3. **Z-Sweep Weight Folding**: Scaling a 64B-aligned copy of each term's $\le 5$ active Z-sweep 1D factors (`z_stack_buf` / `z_op_buf`) by `weight[t]` once per sub-batch eliminates per-term $(2k)^3 N_{\text{tgt}}$ temporary Z-buffer zeroing and `aligned_axpy(4096)` memory traffic.

### 4.2 Dual-Direction Pruned 3D Sweep Graphs, Target-Brick Tiling, 64-Bit Packed Keys & Bounded Scratch

- **64-Bit Packed Spatial Coordinates (`windowed_ns_halo.h`)**: At any fixed tree level $j \le 20$ ($0 \le l_x, l_y, l_z < 2^{20} < 2^{21}$), spatial coordinates are packed into a single `uint64_t` `(lx << 42) | (ly << 21) | lz` whose unsigned `<` order matches `windowed_key_less` 1:1 (with automatic fallback to full `Key<3>` at $j > 20$).
- **Dual-Direction Graph Pruning (`build_sumfact_graph_3d`)**: Selects forward candidate generation when $\lvert\mathcal{D}_{\text{src}}\rvert \le \lvert\mathcal{D}_{\text{tgt}}\rvert$ (with secondary Z-edge direction check `nxy0 <= num_tgt`) and backward generation otherwise.
- **L3-Resident Spatial Target-Brick Tiling (`run_sumfact_target_tiled_3d`)**: Partitions target sets into contiguous spatial bricks (`1,024` target boxes for `k -> k`, `512` for `2k -> k`, and `256` for `2k -> 2k` and `k -> 2k`), building one pruned `SumFactGraph3D` per `(brick, eff_bw)` pair and evaluating fused term sub-batches on that brick before advancing.
- **Pooled 64-Byte-Aligned Uninitialized `SumFactScratch3D` Workspaces (`UninitAlignedBuf<T>` + `SumFactScratchLease<T>`)**: Uses 64B-aligned uninitialized storage (`kMaxParallelChunks = 128` slots) lazily allocated per worker-thread slot `cid` (`WindowedThreadSlotRegistry`), frees prior buffers *before* larger reallocations, and pools up to `min(128, ThreadPool::size() + 2)` instances trimmed at `2 MiB` per `buf_x`/`buf_xy` buffer.

### 4.3 Native MADNESS Task-Model Over-Decomposition (`madness::Range` + `WorldTaskQueue` / `tbb_arena`)

1. **Worker-Thread Slot Registry & Over-Decomposed `madness::Range` Dispatch (`windowed_ns_halo.h`, `windowed_ns_sumfact.h`)**:
   - `WindowedThreadSlotRegistry` + `ScopedWorkerSlot` decouple per-thread scratch storage (`O(N_threads)` up to `128` slots, with automatic overflow-slot acquisition across reentrant `Future::get` / `ThreadPool::await` task execution) from task count.
   - `parallel_for_chunks` decomposes `[0, n_items)` into `kOverdecompPerThread * (ThreadPool::size() + 1)` (`4x` over-decomposition) fine-grained `madness::Range<std::size_t>` blocks (bounded below by `min_grain`, with 3D sum-factorization sweeps using FLOP-calibrated `compute_sumfact_sweep_grain` $\in [1, 16]$ targeting $\ge 48\text{k}$ FLOPs/block). Blocks execute via `world.taskq.for_each` + coarse split-leader `(blk & 3) == 0` `thread_purge()` on Pthreads/PaRSEC and `ThreadPool::tbb_arena->execute` + `tbb::parallel_for(madness::Range)` on Intel TBB, dynamically load-balancing across heterogeneous big.LITTLE (P/E) cores and dynamically clocked CPUs.
2. **Two-Priority Outer Orbital/Pair Queueing & Multi-Rank Lockstep Tiling (`windowed_convolution.cc`, `windowed_ns_sumfact.h`)**:
   - **Single-Rank / Subworld (`world.size() == 1`, including `MacroTaskExchangeRow` `ntile = 10` & `ApplyTask` `nfunc = 4`)**: Vector batches sort inputs in LPT (`node_counts` descending) order and tile outer task submission (`windowed_for_each_outer_tasks`: normal-priority `world.taskq.add` at the back of `DQueue`, or `tbb_arena` `parallel_for`) according to `outer_batch_tile`: all `N` for `max_nodes < 1,024`, `min(N, 4)` for `< 16,384`, `min(N, 2)` for `< 32,768`, and `1` for `>= 32,768` (overridable via `MADNESS_WINDOWED_BATCH_WORKERS`). While `g_unfinished_outer_tasks >= total_threads`, inner `parallel_for_chunks` sweeps run inline in L1/L2 cache with zero `ForEachTask` allocations or `DQueue` contention; as soon as the batch enters its tail (`g_unfinished_outer_tasks < total_threads`), inner `hipri` `ForEachTask`s activate automatically at the front of `DQueue` so fast P-cores help finish straggling orbitals without barrier stalls.
   - **Distributed MPI (`world.size() > 1`)**: Distinct-input batches execute in lockstep tiles of `get_multi_rank_batch_tile()` orbitals (`16` by default, configurable in `[1, 256]` via `MADNESS_WINDOWED_MPI_BATCH_TILE`) through `apply_windowed_ns_3d_batch_distributed`, coalescing Stage 1b (`4–5` collectives) and Stage 5a (`2–3` collectives) across all $B$ orbitals in the tile while dispatching rank-local Stages 2–4 across `windowed_for_each_outer_tasks` in LPT order whenever `max_local_real < 512` (e.g., distributed `HF Exchange` pair batches) and running sequentially with full `4x` inner `madness::Range` over-decomposition on larger local trees.

---

## 5. Verification Results (`test_windowed_convolution.cc` & Legacy MRA Suites)

`test_windowed_convolution` implements an adversarial 4-Gate verification matrix (`113` checks in `--full` mode, `63` checks in `--small` CTest mode):

1. **Gate 1 — Commutation, 3D GC Spectral Norm, Tail Quadrature & Deep Keys (`15 / 15 PASS`)**:
   - Scale-locked stiffness drift $\max \lvert\gamma_{\text{ref}}(w_{\text{eff}}, m) - \gamma_{p,j}\rvert / \gamma_{\text{ref}} = \mathbf{6.28 \times 10^{-15}} \le 10^{-13}$.
   - Theorem 3.1 exact boundary-halo commutation ($k \in \{6, 8, 10, 12\}$, level jumps $\Delta j \in \{1, 5, 15\}$): 1D Frobenius error $\mathbf{7.33 \times 10^{-16}}$, 3D tensor relative error $\mathbf{1.97 \times 10^{-14}}$, 3D partition-of-unity round-trip $\mathbf{3.29 \times 10^{-16}}$.
   - 3D Gauss–Christoffel spectral 2-norm error up to $n = 30$: $\mathbf{1.87 \times 10^{-14}}$ (Coulomb) and $\mathbf{1.68 \times 10^{-18}}$ (BSH).
   - Pointwise radial kernel error ($\epsilon = 10^{-9}$, $r \in [10^{-8}, \sqrt{3}]$): $\mathbf{3.18 \times 10^{-11}}$ (Coulomb, $L = 1, 8$, zero outer-tail bump) and $\mathbf{2.31 \times 10^{-12}}$ (BSH, $\mu = 0.5, 2.0$).
   - `WindowedKeyHash` & `WorldContainer` collision-free up to **`Level = 35`** (`665` extreme/midplane keys).

2. **Gate 2 — Unbalanced 3D Coarse-to-Fine Interface Jumps (`42 / 42 PASS`)**:
   - Frozen coarse octants (`level = 5`) directly abutting fine subtrees (`level = 8, 11, 15`, $\Delta j = 3, 6, 10$) across faces, edges, and corners on $L = 1$ and $L = 8$: global $L_2$ error **`2.19e-10 – 6.28e-10`**, face/edge/corner pointwise error **`2.19e-10 – 1.20e-8`**, zero `erf` interface ringing.

3. **Gate 3 — 22–25 Level All-Electron Slater Cusps, Sparse Compressed Trees & Batched Linearity (`51 / 51 PASS`)**:

   | Case | Octree Depth | Output $V$ Nodes | Global Rel $L_2$ Error | Cusp $\lvert V(\bm{R}_A)-Z\rvert / Z$ | Hartree $\lvert E_H - \frac{5}{16}Z\rvert / (\frac{5}{16}Z)$ |
   | :--- | :---: | :---: | :---: | :---: | :---: |
   | `Z=20, k=8, L=4` | `23` | `12,681` | `2.34e-09` | `2.87e-09` | `8.09e-11` |
   | `Z=20, k=10, L=4` | `22` | `4,873` | `2.68e-09` | `2.55e-10` | `4.22e-11` |
   | `Z=50, k=8, L=4` | `25` | `25,769` | `2.59e-09` | `2.91e-09` | `1.17e-10` |
   | `Z=50, k=10, L=4` | `24` | `7,273` | `2.66e-09` | `2.56e-10` | `3.27e-11` |
   | `Z=50, k=10, L=1` | `23` | `6,601` | `3.62e-10` | `2.57e-10` | `5.81e-11` |

   - Sparse `vtransform(..., 1e-4)` `TreeState::compressed` BSH parity vs `TreeState::reconstructed`: **`0.00e+00 – 2.09e-15 <= 1e-13`** with input `TreeState::compressed` preserved unmutated.
   - Shared-tree `Function<double_complex, 3>` linearity and heterogeneous 6-element batched vector `apply` vs scalar `apply` (including multi-tile `MADNESS_WINDOWED_MPI_BATCH_TILE=2` on 4 MPI ranks): **`0.000000e+00 <= 1e-13`**.

4. **Gate 4 — Fast-Math + TBB + `mimalloc` Full & Multi-Rank MPI Parity (`113 / 113 PASS`)**:
   - Full 4-Gate suite (`--full`) under `-march=native -ffast-math` + Intel TBB + `mimalloc`: **`113 / 113 PASSED, 0 FAILED`** in **`42.4 s`** on `1` MPI rank (`4.53x` Coulomb & `3.40x` BSH speedups on `Z = 30` depth-22 cusp) and **`113 / 113 PASSED, 0 FAILED`** on `4` MPI ranks, with identical 1-rank vs. 4-rank output node counts (`7,689` and `7,625`), identical cross-operator relative $L_2$ diffs (`2.426253e-8` and `3.575509e-7`), and `0` leaked `FunctionImpl` or `WorldObject` instances.

---

## 6. Performance & Timing Results

All production benchmarks below were built with:
```bash
cmake -GNinja \
  -DMADNESS_TARGET_ARCH=native \
  -DMADNESS_RELAXED_MATH=fast \
  -DMADNESS_TASK_BACKEND=TBB \
  -DENABLE_TCMALLOC_MINIMAL=OFF \
  -DENABLE_GPERFTOOLS=OFF
```
and executed sequentially with **`mimalloc v3.5.3`** (`LD_PRELOAD=libmimalloc.so.3.5`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, Intel MKL `HAVE_MTXMQ`, `--dft="restart=none; save=false"`, `aobasis=3-21g` on `sr`).

### 6.1 `madqc` All-Electron Hartree–Fock 11-Structure Benchmark Suite (`1` MPI $\times$ `11` Threads, `mimalloc v3.5.3`)

- **Before**: `MADNESS_WINDOWED_NS=0` (legacy `SeparatedConvolution`)
- **After**: `MADNESS_WINDOWED_NS=1` (`WindowedSeparatedConvolution` for Coulomb, HF exchange, and BSH)
- Default `madqc --wf=scf` two-stage protocol: **`1e-4` (`k = 6`)** followed by **`1e-6` (`k = 8`)**, spanning single atoms ($Z = 10\text{--}38$), organic molecules (`c2h4`), and 3D water clusters $(\text{H}_2\text{O})_n$ up to **`water20`** ($60$ atoms, $100$ MOs, $L = 120.56$ Bohr).

#### End-to-End Summary (`11` Structures: `1.19x` Cumulative Speedup, `11 / 11` Fewer MO Boxes, `11 / 11` Lower Peak RSS)

| Structure | $N_{\text{MO}}$ | $E_{\text{tot}}$ Before (Ha) | $E_{\text{tot}}$ After (Ha) | $\lvert\Delta E\rvert$ (Ha) | Iters (`1e-4`+`1e-6`) Bef $\rightarrow$ Aft | Total Wall Bef (s) | Total Wall Aft (s) | Speedup | Final MO Boxes (`1e-6`) Bef $\rightarrow$ Aft | MO Box Change | Peak RSS Bef $\rightarrow$ Aft (MB) | RSS Change |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`water`** | 5 | -76.06793430 | -76.06793479 | `4.88e-7` | `8+2` $\rightarrow$ `8+2` | 25.8 | **15.0** | **1.72x** | 5,784 $\rightarrow$ **4,632** | **-19.9%** | 1,430 $\rightarrow$ **1,183** | **-17.3%** |
| **`water2`** | 10 | -152.14150118 | -152.14150246 | `1.28e-6` | `8+2` $\rightarrow$ `8+2` | 62.7 | **45.2** | **1.39x** | 14,008 $\rightarrow$ **11,752** | **-16.1%** | 2,336 $\rightarrow$ **1,869** | **-20.0%** |
| **`water5`** | 25 | -380.38040855 | -380.38041462 | `6.07e-6` | `8+2` $\rightarrow$ `8+2` | 233.5 | **206.2** | **1.13x** | 44,496 $\rightarrow$ **37,488** | **-15.7%** | 5,484 $\rightarrow$ **5,045** | **-8.0%** |
| **`water9`** | 45 | -684.69877329 | -684.69877629 | `3.00e-6` | `8+3` $\rightarrow$ `8+3` | 785.4 | **717.5** | **1.09x** | 86,384 $\rightarrow$ **82,384** | **-4.6%** | 11,779 $\rightarrow$ **11,037** | **-6.3%** |
| **`water20`** | 100 | -1521.28176475 | -1521.28177241 | `7.66e-6` | `10+3` $\rightarrow$ `10+3` | 2,981.7 | 3,076.4 | 0.97x\* | 242,432 $\rightarrow$ **234,136** | **-3.4%** | 28,148 $\rightarrow$ **26,336** | **-6.4%** |
| **`c2h4`** | 8 | -78.06952528 | -78.06952574 | `4.62e-7` | `8+2` $\rightarrow$ `8+2` | 50.6 | **34.4** | **1.47x** | 11,328 $\rightarrow$ **10,368** | **-8.5%** | 2,028 $\rightarrow$ **1,719** | **-15.2%** |
| **`ne`** | 5 | -128.54700003 | -128.54700036 | `3.33e-7` | `6+2` $\rightarrow$ `6+2` | 23.2 | **10.4** | **2.23x** | 4,328 $\rightarrow$ **3,688** | **-14.8%** | 1,337 $\rightarrow$ **1,007** | **-24.7%** |
| **`mg`** | 6 | -199.61453509 | -199.61453631 | `1.21e-6` | `5+2` $\rightarrow$ `5+2` | 29.3 | **13.4** | **2.19x** | 6,328 $\rightarrow$ **4,656** | **-26.4%** | 1,514 $\rightarrow$ **1,096** | **-27.6%** |
| **`ar`** | 9 | -526.81740404 | -526.81740803 | `3.99e-6` | `6+2` $\rightarrow$ `6+2` | 67.4 | **27.8** | **2.42x** | 12,664 $\rightarrow$ **7,880** | **-37.8%** | 2,093 $\rightarrow$ **1,372** | **-34.5%** |
| **`ca`** | 10 | -676.75807488 | -676.75807986 | `4.99e-6` | `5+3` $\rightarrow$ `5+2` | 109.0 | **31.0** | **3.52x** | 17,992 $\rightarrow$ **8,992** | **-50.0%** | 2,971 $\rightarrow$ **1,487** | **-50.0%** |
| **`sr`** | 19 | -3131.54555056 | -3131.54557463 | `2.41e-5` | `12+5` $\rightarrow$ `10+3` | 810.8 | **188.2** | **4.31x** | 42,600 $\rightarrow$ **21,400** | **-49.8%** | 6,763 $\rightarrow$ **3,048** | **-54.9%** |
| **TOTAL** | **242** | — | — | — | — | **5,179.4** | **4,365.5** | **1.19x** | 488,344 $\rightarrow$ **427,376** | **-12.5%** | — | **11 / 11 Lower** |

*\*On `water20` ($100$ MOs), flooring operator `gamma_spiky = 270.0` at $k = 6$ + $H^1$-safe union-skeleton pruning + graded `coulop` refinement eliminated the `25`-iteration KAIN residual stall (`4,800.8 s` $\rightarrow$ **`3,076.4 s`**, `1.56x` faster), converging cleanly in `10+3` iterations while outperforming Legacy in the production `1e-6, k = 8` protocol (`1,821.4 s` $\rightarrow$ **`1,661.5 s`**, `1.10x` faster) and saving **`-1.81 GB` (`-6.4%`) Peak RSS**.*

#### Protocol 1 Breakdown (`thresh = 1e-4, k = 6`)

| Structure | Iters Bef $\rightarrow$ Aft | Proto Span Bef $\rightarrow$ Aft (s) | Coulomb Bef $\rightarrow$ Aft (s) | HF Exchange Bef $\rightarrow$ Aft (s) | Apply BSH Bef $\rightarrow$ Aft (s) | MO Leaf Boxes Bef $\rightarrow$ Aft |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`water`** | `8` $\rightarrow$ `8` | 13.10 $\rightarrow$ **9.40** (**1.39x**) | 0.76 $\rightarrow$ **0.47** | 8.24 $\rightarrow$ **4.57** (**1.80x**) | 1.92 $\rightarrow$ 2.29 | 3,240 $\rightarrow$ 3,240 |
| **`water2`** | `8` $\rightarrow$ `8` | 29.80 $\rightarrow$ **27.40** (**1.09x**) | 1.15 $\rightarrow$ **0.85** | 20.09 $\rightarrow$ **16.48** (**1.22x**) | 4.19 $\rightarrow$ 5.82 | 6,696 $\rightarrow$ **6,664** |
| **`water5`** | `8` $\rightarrow$ `8` | 104.70 $\rightarrow$ 116.60 (0.90x) | 2.24 $\rightarrow$ **2.05** | 75.76 $\rightarrow$ 78.97 (0.96x) | 12.60 $\rightarrow$ 21.77 | 18,144 $\rightarrow$ **18,104** |
| **`water9`** | `8` $\rightarrow$ `8` | 258.50 $\rightarrow$ 309.10 (0.84x) | 3.67 $\rightarrow$ **3.57** | 196.86 $\rightarrow$ 221.31 (0.89x) | 25.07 $\rightarrow$ 52.08 | 34,472 $\rightarrow$ **34,424** |
| **`water20`** | `10` $\rightarrow$ `10` | 1,118.10 $\rightarrow$ 1,371.80 (0.82x) | 9.38 $\rightarrow$ 9.94 | 890.43 $\rightarrow$ 1,022.05 (0.87x) | 80.70 $\rightarrow$ 202.55 | 80,320 $\rightarrow$ **80,272** |
| **`c2h4`** | `8` $\rightarrow$ `8` | 24.40 $\rightarrow$ **20.40** (**1.20x**) | 1.08 $\rightarrow$ **0.78** | 16.34 $\rightarrow$ **11.49** (**1.42x**) | 3.22 $\rightarrow$ 4.41 | 5,760 $\rightarrow$ **5,696** |
| **`ne`** | `6` $\rightarrow$ `6` | 10.80 $\rightarrow$ **6.00** (**1.80x**) | 0.68 $\rightarrow$ **0.29** | 7.13 $\rightarrow$ **3.03** (**2.35x**) | 1.53 $\rightarrow$ **1.28** | 3,048 $\rightarrow$ 3,048 |
| **`mg`** | `5` $\rightarrow$ `5` | 11.60 $\rightarrow$ **7.30** (**1.59x**) | 0.52 $\rightarrow$ **0.25** | 8.03 $\rightarrow$ **4.03** (**1.99x**) | 1.42 $\rightarrow$ **1.39** | 3,888 $\rightarrow$ 3,888 |
| **`ar`** | `6` $\rightarrow$ `6` | 28.70 $\rightarrow$ **16.20** (**1.77x**) | 0.74 $\rightarrow$ **0.30** | 22.45 $\rightarrow$ **10.70** (**2.10x**) | 3.05 $\rightarrow$ **2.68** | 6,024 $\rightarrow$ 6,024 |
| **`ca`** | `5` $\rightarrow$ `5` | 28.40 $\rightarrow$ **16.60** (**1.71x**) | 0.60 $\rightarrow$ **0.30** | 22.36 $\rightarrow$ **11.12** (**2.01x**) | 2.88 $\rightarrow$ **2.52** | 6,736 $\rightarrow$ 6,736 |
| **`sr`** | `12` $\rightarrow$ `10` | 238.00 $\rightarrow$ **110.30** (**2.16x**) | 2.12 $\rightarrow$ **0.90** | 204.82 $\rightarrow$ **87.86** (**2.33x**) | 20.57 $\rightarrow$ **12.94** | 13,720 $\rightarrow$ 13,720 |

#### Protocol 2 Breakdown (`thresh = 1e-6, k = 8`: `11 / 11` Faster, `1.10x` – `7.67x`)

| Structure | Iters Bef $\rightarrow$ Aft | Proto Span Bef $\rightarrow$ Aft (s) | Proto Speedup | Coulomb Bef $\rightarrow$ Aft (s) | HF Exchange Bef $\rightarrow$ Aft (s) | Apply BSH Bef $\rightarrow$ Aft (s) | MO Leaf Boxes Bef $\rightarrow$ Aft |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`water`** | `2` $\rightarrow$ `2` | 11.39 $\rightarrow$ **4.51** | **2.53x** | 0.86 $\rightarrow$ **0.23** (**3.74x**) | 7.66 $\rightarrow$ **2.54** (**3.02x**) | 1.93 $\rightarrow$ **0.89** (**2.17x**) | 5,784 $\rightarrow$ **4,632** (-19.9%) |
| **`water2`** | `2` $\rightarrow$ `2` | 30.74 $\rightarrow$ **15.92** | **1.93x** | 1.41 $\rightarrow$ **0.51** (**2.76x**) | 22.60 $\rightarrow$ **10.77** (**2.10x**) | 4.32 $\rightarrow$ **2.48** (**1.74x**) | 14,008 $\rightarrow$ **11,752** (-16.1%) |
| **`water5`** | `2` $\rightarrow$ `2` | 123.34 $\rightarrow$ **84.55** | **1.46x** | 3.36 $\rightarrow$ **1.48** (**2.27x**) | 97.49 $\rightarrow$ **64.31** (**1.52x**) | 13.46 $\rightarrow$ **10.34** (**1.30x**) | 44,496 $\rightarrow$ **37,488** (-15.7%) |
| **`water9`** | `3` $\rightarrow$ `3` | 515.23 $\rightarrow$ **397.79** | **1.30x** | 9.08 $\rightarrow$ **4.52** (**2.01x**) | 407.27 $\rightarrow$ **308.01** (**1.32x**) | 69.12 $\rightarrow$ **57.38** (**1.20x**) | 86,384 $\rightarrow$ **82,384** (-4.6%) |
| **`water20`** | `3` $\rightarrow$ `3` | 1,821.38 $\rightarrow$ **1,661.45** | **1.10x** | 18.31 $\rightarrow$ **10.13** (**1.81x**) | 1,501.19 $\rightarrow$ **1,333.68** (**1.13x**) | 176.77 $\rightarrow$ 194.90 (0.91x) | 242,432 $\rightarrow$ **234,136** (-3.4%) |
| **`c2h4`** | `2` $\rightarrow$ `2` | 24.05 $\rightarrow$ **12.20** | **1.97x** | 1.16 $\rightarrow$ **0.39** (**2.97x**) | 17.51 $\rightarrow$ **7.97** (**2.20x**) | 3.49 $\rightarrow$ **1.97** (**1.77x**) | 11,328 $\rightarrow$ **10,368** (-8.5%) |
| **`ne`** | `2` $\rightarrow$ `2` | 11.14 $\rightarrow$ **3.22** | **3.46x** | 0.80 $\rightarrow$ **0.12** (**6.67x**) | 7.75 $\rightarrow$ **1.85** (**4.19x**) | 1.90 $\rightarrow$ **0.60** (**3.17x**) | 4,328 $\rightarrow$ **3,688** (-14.8%) |
| **`mg`** | `2` $\rightarrow$ `2` | 16.13 $\rightarrow$ **4.63** | **3.48x** | 0.91 $\rightarrow$ **0.14** (**6.50x**) | 11.82 $\rightarrow$ **2.88** (**4.10x**) | 2.48 $\rightarrow$ **0.81** (**3.06x**) | 6,328 $\rightarrow$ **4,656** (-26.4%) |
| **`ar`** | `2` $\rightarrow$ `2` | 36.70 $\rightarrow$ **9.93** | **3.70x** | 1.41 $\rightarrow$ **0.21** (**6.71x**) | 29.75 $\rightarrow$ **7.26** (**4.10x**) | 3.94 $\rightarrow$ **1.26** (**3.13x**) | 12,664 $\rightarrow$ **7,880** (-37.8%) |
| **`ca`** | `3` $\rightarrow$ `2` | 78.30 $\rightarrow$ **12.33** | **6.35x** | 2.74 $\rightarrow$ **0.25** (**10.96x**) | 61.74 $\rightarrow$ **9.16** (**6.74x**) | 10.75 $\rightarrow$ **1.51** (**7.12x**) | 17,992 $\rightarrow$ **8,992** (-50.0%) |
| **`sr`** | `5` $\rightarrow$ `3` | 568.60 $\rightarrow$ **74.10** | **7.67x** | 9.02 $\rightarrow$ **0.51** (**17.69x**) | 466.92 $\rightarrow$ **62.71** (**7.45x**) | 78.97 $\rightarrow$ **6.43** (**12.28x**) | 42,600 $\rightarrow$ **21,400** (-49.8%) |

#### Hybrid Multi-Rank Lockstep MPI + Threads (`4` MPI Processes $\times$ `MAD_NUM_THREADS=2`, `mimalloc v3.5.3`)

| Structure ($N_{\text{MO}}$) | Total Wall Bef $\rightarrow$ Aft (s) | Speedup | `1e-4` Span (Coul / K / BSH) Aft (s) | `1e-6` Span Bef $\rightarrow$ Aft (s) | `1e-6` Speedup (`BSH` Speedup) | Final MO Boxes 4-Rank (1-Rank) | Max Rank RSS Bef $\rightarrow$ Aft (MB) | Converged $E_{\text{tot}}$ Aft (Ha) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`water2`** (10) | 165.9 $\rightarrow$ **103.2** | **1.61x** | 52.40 (`1.75` / `21.90` / `13.13`) | 74.58 $\rightarrow$ **31.32** | **2.38x** (**2.06x**) | **11,752** (11,752) | 1,404.0 $\rightarrow$ 1,570.5\* | `-152.14150246` |
| **`c2h4`** (8) | 106.6 $\rightarrow$ **68.2** | **1.56x** | 33.10 (`1.33` / `13.67` / `8.15`) | 46.47 $\rightarrow$ **22.70** | **2.05x** (**2.06x**) | **10,368** (10,368) | 1,270.9 $\rightarrow$ **1,260.3** | `-78.06952574` |
| **`water5`** (25) | 512.6 $\rightarrow$ **387.4** | **1.32x** | 189.40 (`3.97` / `92.88` / `42.86`) | 262.63 $\rightarrow$ **151.00** | **1.74x** (**1.56x**) | **37,488** (37,488) | 2,454.6 $\rightarrow$ 2,943.3\* | `-380.38041462` |

*\*With default `MADNESS_WINDOWED_MPI_BATCH_TILE=16`, up to `16` orbitals share lockstep MPI halo/restriction barriers simultaneously (`16 + 9` tiles on `water5`), trading `~165–488 MB/rank` of transient batch workspace (`41.7%` below 1-rank `water5`'s `5,045 MB`) for `~20x` fewer global MPI fences and exact 1-rank vs. 4-rank leaf-box and `< 0.01` $\mu E_h$ energy invariance (tunable via `MADNESS_WINDOWED_MPI_BATCH_TILE`).*

#### 2-Node DGX ARMv9.2 big.LITTLE Cluster (`4` MPI $\times$ `8` Pthreads Across `2` $\times$ 20-Core Cortex-X925/A725 Nodes, `256 GB` RAM, `12` Structures Up to `water47`)

Executed across `dgx2` + `dgx` (`2` MPI ranks/node $\times$ `MAD_NUM_THREADS=8`, `--bind-to none`, GCC 13.3 `-mcpu=native -ffast-math`, `ARMPL`, `libxsmm` aarch64 JIT, system `glibc malloc`):

| Structure ($N_{\text{MO}}$) | Wall Bef $\rightarrow$ Aft (s) | End-to-End Speedup | `1e-6` Span Bef $\rightarrow$ Aft (s) | **`1e-6` Speedup** | `1e-6` `K` / `BSH` Aft (Bef) (s) | Final MO Boxes Bef $\rightarrow$ Aft | Max Host RSS Bef $\rightarrow$ Aft (MB) | $\Delta E$ ($\mu E_h$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`water`** (5) | 17.3 $\rightarrow$ **14.6** | **1.18x** | 5.50 $\rightarrow$ **3.50** | **1.57x** | **1.29** (2.21) / **0.68** (1.53) | 5,784 $\rightarrow$ **4,632** (-19.9%) | 1,937 $\rightarrow$ **1,773** | `-0.49` |
| **`water2`** (10) | 29.0 $\rightarrow$ **26.4** | **1.10x** | 10.10 $\rightarrow$ **7.50** | **1.35x** | **3.33** (5.18) / **1.33** (1.37) | 14,008 $\rightarrow$ **11,752** (-16.1%) | 2,632 $\rightarrow$ 2,766 | `-1.29` |
| **`water5`** (25) | 88.2 $\rightarrow$ 89.5 | 0.99x | 34.80 $\rightarrow$ **30.40** | **1.14x** | **20.19** (23.37) / 3.45 (3.23) | 44,496 $\rightarrow$ **37,496** (-15.7%) | 4,907 $\rightarrow$ **4,824** | `-6.13` |
| **`water9`** (45) | 258.4 $\rightarrow$ 260.0 | 0.99x | 135.20 $\rightarrow$ **123.90** | **1.09x** | **93.80** (99.70) / **14.22** (18.20) | 86,384 $\rightarrow$ **82,360** (-4.7%) | 9,262 $\rightarrow$ 9,597 | `-3.01` |
| **`water20`** (100) | 1,054.1 $\rightarrow$ 1,118.8 | 0.94x | 531.60 $\rightarrow$ **510.60** | **1.04x** | **418.22** (427.11) / **43.75** (51.67) | 242,184 $\rightarrow$ **234,104** (-3.3%) | 23,800 $\rightarrow$ **21,788** | `-7.79` |
| **`water47`** (235) | 5,456.6 $\rightarrow$ 5,869.9 | 0.93x | 2,209.70 $\rightarrow$ 2,261.80 | 0.98x | 1,939.64 (1,876.53) / 146.25 (144.78) | 654,912 $\rightarrow$ **634,960** (-3.0%) | 65,899 $\rightarrow$ **60,646** | `-19.98` |
| **`c2h4`** (8) | 22.0 $\rightarrow$ **20.1** | **1.09x** | 7.90 $\rightarrow$ **5.70** | **1.39x** | **2.80** (4.14) / **0.89** (1.44) | 11,328 $\rightarrow$ **10,368** (-8.5%) | 2,446 $\rightarrow$ 2,486 | `-0.46` |
| **`ne`** (5) | 13.3 $\rightarrow$ **10.2** | **1.30x** | 5.10 $\rightarrow$ **2.80** | **1.82x** | **1.14** (2.20) / **0.53** (1.22) | 4,152 $\rightarrow$ **3,688** (-11.2%) | 1,800 $\rightarrow$ **1,439** | `-0.29` |
| **`mg`** (6) | 14.4 $\rightarrow$ **11.0** | **1.31x** | 6.20 $\rightarrow$ **3.40** | **1.82x** | **1.45** (3.29) / **0.54** (1.10) | 6,672 $\rightarrow$ **4,656** (-30.2%) | 2,076 $\rightarrow$ **1,625** | `-1.34` |
| **`ar`** (9) | 26.5 $\rightarrow$ **16.9** | **1.57x** | 11.80 $\rightarrow$ **5.60** | **2.11x** | **2.89** (7.34) / **0.91** (1.63) | 12,992 $\rightarrow$ **7,752** (-40.3%) | 3,033 $\rightarrow$ **2,109** | `-4.39` |
| **`ca`** (10) | 34.6 $\rightarrow$ **17.4** | **1.99x** | 20.60 $\rightarrow$ **6.20** | **3.32x** | **3.42** (13.05) / **0.75** (3.54) | 15,768 $\rightarrow$ **9,048** (-42.6%) | 3,720 $\rightarrow$ **2,099** | `-3.58` |
| **`sr`** (19) | 228.6 $\rightarrow$ **80.4** | **2.84x** | 128.50 $\rightarrow$ **29.10** | **4.42x** | **21.54** (93.25) / **2.73** (21.68) | 42,296 $\rightarrow$ **21,424** (-49.3%) | 7,076 $\rightarrow$ **3,354** | `-22.13` |
| **CUMULATIVE** | 7,243.0 $\rightarrow$ 7,535.2 | 0.96x | **3,107.0 $\rightarrow$ 2,990.5** | **1.04x** | — | 1,140,976 $\rightarrow$ **1,062,240** (**-6.9%**) | 245.3 $\rightarrow$ **220.1 GB** (**-10.3%**) | — |

#### Hybrid GGA `B3LYP` Benchmark Suite (`--dft="xc b3lyp"`, `4` MPI $\times$ `8` Pthreads on 2-Node DGX ARM Cluster, `11` Structures Up to `water20`)

Because hybrid GGA `B3LYP` (`0.20` exact HF exchange + `Libxc` GGA density-gradient differentiation $\nabla\rho$, $\nabla\cdot(v_\sigma\nabla\rho)$) differentiates orbital densities spatially inside `XCOperator`, Legacy's coarse-to-fine Gibbs (`erf`) interface oscillations get amplified in `DFT potential` and `V*psi`. Eliminating Gibbs oscillations via Windowed NS accelerates **both** the integral operators (`Coulomb`, `HF Exchange`, `Apply BSH`) **and** the untouched GGA `DFT potential` pass (`11 / 11` faster at `1e-6`, `1.29x` faster even at `1e-4` on `water20`) while saving `-5.52 GB` / host (`-14.0%`) Peak RSS on `water20`:

| Structure ($N_{\text{MO}}$) | Iters Bef $\rightarrow$ Aft | Wall Bef $\rightarrow$ Aft (s) | **Speedup** | `1e-6` Span Bef $\rightarrow$ Aft (s) | **`1e-6` Speedup** | `1e-6` `K` / **`DFT pot`** / `BSH` Aft (Bef) (s) | Final MO Boxes Bef $\rightarrow$ Aft | Max Host RSS Bef $\rightarrow$ Aft (MB) | $\Delta E$ ($\mu E_h$) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **`water`** (5) | `8+2` $\rightarrow$ `8+2` | 26.5 $\rightarrow$ **17.9** | **1.48x** | 7.50 $\rightarrow$ **4.70** | **1.60x** | **1.40** (2.39) / **0.58** (0.62) / **1.08** (2.33) | 5,576 $\rightarrow$ **4,760** (-14.6%) | 2,143 $\rightarrow$ **2,111** | `-0.16` |
| **`water2`** (10) | `9+2` $\rightarrow$ `9+2` | 37.0 $\rightarrow$ **33.6** | **1.10x** | 12.70 $\rightarrow$ **9.20** | **1.38x** | **3.46** (5.86) / **0.92** (1.04) / **1.79** (2.24) | 13,832 $\rightarrow$ **12,240** (-11.5%) | 3,026 $\rightarrow$ 3,627 | `-0.62` |
| **`water5`** (25) | `8+3` $\rightarrow$ `8+3` | 129.7 $\rightarrow$ **126.0** | **1.03x** | 63.60 $\rightarrow$ **58.40** | **1.09x** | **35.88** (37.97) / **2.69** (3.09) / **8.52** (10.26) | 42,696 $\rightarrow$ **41,016** (-3.9%) | 7,868 $\rightarrow$ **6,944** (-11.7%) | `-0.79` |
| **`water9`** (45) | `9+3` $\rightarrow$ `9+3` | 311.2 $\rightarrow$ 325.9 | 0.95x | 158.80 $\rightarrow$ **146.20** | **1.09x** | **105.00** (108.91) / **4.01** (4.74) / **16.95** (23.20) | 93,352 $\rightarrow$ **90,872** (-2.7%) | 13,579 $\rightarrow$ **12,138** (-10.6%) | `-1.43` |
| **`water20`** (100) | `10+3` $\rightarrow$ `10+3` | 1,249.6 $\rightarrow$ 1,323.7 | 0.94x | 647.70 $\rightarrow$ **636.70** | **1.02x** | 488.87 (482.53) / **9.75** (11.19)\* / **62.79** (71.65) | 278,464 $\rightarrow$ **273,680** (-1.7%) | 39,407 $\rightarrow$ **33,887** (**-14.0%**) | `-5.08` |
| **`c2h4`** (8) | `8+2` $\rightarrow$ `8+2` | 26.8 $\rightarrow$ **24.2** | **1.11x** | 9.70 $\rightarrow$ **7.30** | **1.33x** | **2.82** (4.19) / **0.67** (0.80) / **1.62** (2.06) | 11,648 $\rightarrow$ **10,944** (-6.0%) | 2,909 $\rightarrow$ 3,563 | `-0.44` |
| **`ne`** (5) | `7+2` $\rightarrow$ `7+2` | 17.0 $\rightarrow$ **12.7** | **1.34x** | 6.00 $\rightarrow$ **3.50** | **1.71x** | **1.00** (2.27) / **0.42** (0.52) / **0.90** (1.67) | 4,264 $\rightarrow$ **3,688** (-13.5%) | 1,986 $\rightarrow$ **1,799** (-9.4%) | `-0.32` |
| **`mg`** (6) | `6+2` $\rightarrow$ `6+2` | 18.3 $\rightarrow$ **14.2** | **1.29x** | 7.40 $\rightarrow$ **4.40** | **1.68x** | **1.54** (3.30) / **0.47** (0.62) / **0.97** (1.50) | 6,928 $\rightarrow$ **5,040** (-27.3%) | 2,197 $\rightarrow$ **1,971** (-10.3%) | `-1.06` |
| **`ar`** (9) | `7+2` $\rightarrow$ `7+2` | 32.3 $\rightarrow$ **22.2** | **1.45x** | 14.00 $\rightarrow$ **7.00** | **2.00x** | **3.00** (7.55) / **0.63** (0.91) / **1.26** (2.29) | 15,512 $\rightarrow$ **8,120** (-47.7%) | 3,290 $\rightarrow$ **2,682** (-18.5%) | `-5.80` |
| **`ca`** (10) | `6+3` $\rightarrow$ `6+2` | 44.6 $\rightarrow$ **21.8** | **2.05x** | 26.80 $\rightarrow$ **7.40** | **3.62x** | **3.45** (13.05) / **0.67** (1.40) / **1.06** (7.87) | 10,336 $\rightarrow$ **9,600** (-7.1%) | 5,327 $\rightarrow$ **2,408** (-54.8%) | `+0.57` |
| **`sr`** (19) | `13+5` $\rightarrow$ `9+3` | 232.9 $\rightarrow$ **87.4** | **2.66x** | 141.20 $\rightarrow$ **32.00** | **4.41x** | **21.88** (89.62) / **1.37** (4.07) / **3.81** (33.37) | 33,288 $\rightarrow$ **22,032** (-33.8%) | 9,408 $\rightarrow$ **3,851** (-59.1%) | `-9.86` |
| **CUMULATIVE** | — | **2,125.9 $\rightarrow$ 2,009.6** | **1.06x** | **1,095.4 $\rightarrow$ 916.8** | **1.19x** | **668.3** (757.6) / **22.18** (29.00) / **100.75** (158.44) | 515,896 $\rightarrow$ **481,992** (**-6.6%**) | 168.0 $\rightarrow$ **139.8 GB** (**-16.8%**) | — |

*\*On `water20` at `1e-4 (k = 6)`, `DFT potential` also drops `24.36 s` $\rightarrow$ **`18.94 s` (`1.29x` faster)** (`55.59 s` $\rightarrow$ **`48.29 s`, `1.15x` faster** cumulative across all 11 structures).*

#### Allocator Reference: `tcmalloc_minimal` vs. `mimalloc v3.5.3` (`10` Structures, `1` MPI $\times$ `MAD_NUM_THREADS=11`)

Comparing `gperftools` `tcmalloc_minimal` against `LD_PRELOAD=libmimalloc.so.3.5` across the original 10-structure set (`water`, `water2`, `water5`, `water9`, `mg`, `ca`, `sr`, `c2h4`, `ne`, `ar`) demonstrated a uniform **`1.19x` – `1.29x` (`1.26x` cumulative: `1,393.4 s` $\rightarrow$ `1,108.1 s`)** wall-time speedup on `WindowedSeparatedConvolution` (`1.16x` cumulative on Legacy: `2,558.3 s` $\rightarrow$ `2,197.7 s`) by eliminating `tcmalloc_minimal`'s central-freelist spinlock contention during 11-thread TBB cross-thread `Tensor<double>` allocations/frees, while trading ~`23%–41%` higher arena page retention (`purge_delay = 1000 ms`, tunable via `MIMALLOC_PURGE_DELAY`).

---

### 6.2 3D Multi-Scale Gaussian Convolution & Octree Benchmark (`src/examples/large_exp_3d.cc`, `mimalloc v3.5.3`)

Extracted from `large_exp.cc::test6` on $[-11, +11]^3$ (`width = 22.0`, $k = 10$, `thresh = 1e-8`, `truncate_mode = 0`), convolving normalized 3D Gaussian $f(\bm{r}) = (a/\pi)^{3/2} e^{-a\lvert\bm{r}-\bm{r}_0\rvert^2}$ ($a = 16.0$, $\bm{r}_0 = (0.1, 0, 0)$, `652` input leaves on levels `[2..7]`, $h_{\text{user}} = 0.171875$) with $g_b(\bm{r}) = (b/\pi)^{3/2} e^{-b r^2}$ across `iscale = 0 .. 40` ($b = 2^{\text{iscale}-1} \in [0.5, 5.50 \times 10^{11}]$):

| `isc` | `userh/sig` | `expnt_b` | Window Regime | Exact Leaves | Bef Raw Leaves | Bef Trunc Leaves | Bef $L_2$ Err | Bef Time (ms) | Aft Raw Leaves | Aft Trunc Leaves | Aft $L_2$ Err | Aft Time (ms) | Raw Leaf Ratio | Speedup |
| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| 0 | `1.72e-1` | `5.00e-1` | `win[0..3]` | 176 | 512 | 176 | `1.44e-9` | 8.63 | 652 | **176** | `3.79e-8` | **7.22** | 0.79x | **1.20x** |
| 2 | `3.44e-1` | `2.00e+0` | `win[1..4]` | 344 | 904 | 344 | `2.32e-9` | 9.65 | **764** | **344** | `2.97e-8` | **8.19** | **1.18x** | **1.18x** |
| 4 | `6.88e-1` | `8.00e+0` | `win[2..5]` | 568 | 1,520 | 568 | `2.43e-9` | 11.77 | **792** | **568** | `1.71e-8` | 12.06 | **1.92x** | 0.98x |
| 6 | `1.38e+0` | `3.20e+1` | `win[3..6]` | 512 | 2,220 | 624 | `4.74e-8` | 19.80 | **652** | **512** | `8.16e-9` | **16.33** | **3.40x** | **1.21x** |
| 8 | `2.75e+0` | `1.28e+2` | `win[4..7]` | 596 | 2,640 | 1,100 | `8.86e-6` | 17.15 | **652** | **596** | `2.02e-8` | 17.90 | **4.05x** | 0.96x |
| 10 | `5.50e+0` | `5.12e+2` | `win[5..7]` | 652 | 2,864 | 1,268 | `1.56e-4` | 16.49 | **652** | **652** | `2.42e-8` | **14.08** | **4.39x** | **1.17x** |
| 12 | `1.10e+1` | `2.05e+3` | `win[6..7]` | 652 | 2,836 | 1,380 | `5.79e-4` | 15.46 | **652** | **652** | `2.55e-8` | **12.64** | **4.35x** | **1.22x** |
| 13 | `1.56e+1` | `4.10e+3` | `win[7..7]` | 652 | 2,892 | 1,380 | `7.74e-4` | 16.25 | **652** | **652** | `3.06e-8` | **7.42** | **4.44x** | **2.19x** |
| 14 | `2.20e+1` | `8.19e+3` | `win[7..7]` | 652 | 2,864 | 1,380 | `1.86e-3` | 16.75 | **652** | **652** | `2.70e-8` | **7.36** | **4.39x** | **2.27x** |
| 15 | `3.11e+1` | `1.64e+4` | `pure-GC` | 652 | 2,808 | 1,380 | `2.92e-3` | 16.04 | **652** | **652** | `2.96e-8` | **2.26** | **4.31x** | **7.11x** |
| 17 | `6.22e+1` | `6.55e+4` | `pure-GC` | 652 | 2,696 | 1,380 | `3.48e-3` | 15.97 | **652** | **652** | `2.56e-8` | **2.24** | **4.13x** | **7.14x** |
| 24 | `7.04e+2` | `8.39e+6` | `pure-GC` | 652 | 2,276 | 1,100 | `6.36e-4` | 14.05 | **652** | **652** | `2.54e-8` | **2.31** | **3.49x** | **6.09x** |
| 32 | `1.13e+4` | `2.15e+9` | `pure-GC` | 652 | 1,520 | 988 | `4.23e-5` | 11.31 | **652** | **652** | `2.54e-8` | **2.26** | **2.33x** | **5.00x** |
| 40 | `1.80e+5` | `5.50e+11` | `pure-GC` | 652 | 1,492 | 876 | `2.65e-6` | 9.26 | **652** | **652** | `2.54e-8` | **2.23** | **2.29x** | **4.15x** |
| **SUM / WORST (`0..40`)** | — | — | **41 scales** | **24,912** | **83,432** | **41,348** | **`3.48e-3`** | **550.47** | **27,320** | **24,912** | **`3.79e-8`** | **245.25** | **3.05x** | **2.24x** |

- **Zero Large-Exponent Loss of Precision**: Legacy `SeparatedConvolution` degrades to `3.48e-3` $L_2$ error once `userh/sigma > 1` (`iscale >= 6`), whereas `WindowedSeparatedConvolution` maintains **`2.14e-9 – 3.79e-8` $L_2$ error across all 41 scales** (`91,800x` lower worst-case error).
- **Exact Analytical Output Trees**: Legacy raw trees bloat up to `4.44x` (`2,892` vs `652` leaves) and retain `1,380` unprunable leaves (`+768` spurious `L7` interface boxes) even after `.truncate()`. Windowed NS matches the exact analytical truncated leaf count on **all 41/41 scales (`24,912` vs `24,912` total leaves)**.
- **Smooth Intermediate-to-Spiky Speedup**: Fused multi-term bandwidth batching, active-window read-only ingestion, direct `TreeState::reconstructed` output synthesis, and `pure-GC` execution accelerate `iscale = 13..14` (`win[7..7]`) to **`2.19x–2.27x`**, `iscale = 6..12` to **`1.17x–1.22x`**, and `iscale >= 15` (`pure-GC`) to **`3.97x–7.14x`** (`2.24x`–`2.55x` across the entire 41-scale sweep).

---

### 6.3 Deep 22-Level `Z = 30` Slater Cusp Micro-Benchmark (`depth = 22, eps = 1e-7, rho_nodes = 7,625`)

| Operator & Runtime Configuration | Legacy MADNESS | **Windowed NS** | **Speedup vs Legacy & Tree Compression** |
| :--- | :---: | :---: | :--- |
| `Coulomb` (`np = 1, TBB MAD_NUM_THREADS=4, mimalloc, -ffast-math`) | `2.234 s` (`18,313` nodes) | **`0.493 s` (`7,689` nodes)** | **`4.53x` faster** (**`-58.0%` output nodes**) |
| `BSH` (`mu = 1.0, np = 1, TBB MAD_NUM_THREADS=4, mimalloc, -ffast-math`) | `1.466 s` (`14,409` nodes) | **`0.431 s` (`7,625` nodes)** | **`3.40x` faster** (**`-47.1%` output nodes**) |
| `Coulomb` (`np = 4, TBB MAD_NUM_THREADS=2, mimalloc, -ffast-math`) | `3.427 s` (`18,313` nodes) | **`1.748 s` (`7,689` nodes)** | **`1.96x` faster** (**`-58.0%` output nodes**) |
| `BSH` (`mu = 1.0, np = 4, TBB MAD_NUM_THREADS=2, mimalloc, -ffast-math`) | `2.814 s` (`14,409` nodes) | **`1.526 s` (`7,625` nodes)** | **`1.84x` faster** (**`-47.1%` output nodes**) |

---

### 6.4 Relativistic $Z = 80$ Hydrogenic Dirac Solvers (`src/examples/hatom_sf_dirac.cc` & `src/examples/dirac-hatom.cc`)

Relativistic all-electron solvers at heavy nuclear charge ($Z = 80$, mercury-like hydrogenic ion, exact point-nucleus $1s_{1/2}$ Dirac–Coulomb energy $E_{\text{DC}} = -3532.19209293\text{ Ha}$, $\gamma = \sqrt{1 - (Z/c)^2} \approx 0.811906$) provide an extreme stress test for Gibbs (`erf`) interface noise because **unbounded spatial derivative operators are applied directly to integral-operator outputs on every iteration**:

1. **`hatom_sf_dirac`** (`src/examples/hatom_sf_dirac.cc`, 2-component spin-free / eliminated Dirac): Evaluates `phi_new = apply(BSH, rhs)`, differentiates `phi` twice via `pVp = D(D(phi) * V)` (`bspline1` derivatives), and applies `CoulombOperator` (`psimphi_new = apply(coulop, -pVp/(2\pi))`). Differentiating coarse-to-fine polynomial extrapolation (`erf`) artifacts around the steep $Z = 80$ potential amplifies high-frequency interface oscillations by $O(4^d)$, driving severe tree bloat in both `pVp` and subsequent `BSH`/`Coulomb` applies.
2. **`dirac-hatom`** (`src/examples/dirac-hatom.cc`, 4-component complex Dirac spinor): Scales the small spinor components by $c \approx 137.036$ prior to `BSHOperator<3>` (effectively truncating small-component BSH outputs at `thresh / 137`) and then applies $H_D = c\,\boldsymbol{\alpha}\cdot\mathbf{p} + \beta m c^2$ (`Hd(gvpsi)`), multiplying first derivatives of BSH interface noise by $c = 137.036$ and coupling them directly back into the large components (`spinor[0]` and the analytically-zero $m_j = +1/2$ component `spinor[1]`). Under Legacy NS (`MADNESS_WINDOWED_NS=0`), this triggers a runaway adaptive-refinement feedback loop across iterations.

All runs below were executed sequentially on `1` MPI process $\times$ `MAD_NUM_THREADS=11` (`10` TBB worker threads + main thread, `LD_PRELOAD=libmimalloc.so.3.5`, `OMP_NUM_THREADS=1`, `MKL_NUM_THREADS=1`, `600 s` per-run cap):

#### End-to-End Summary Across All $Z = 80$ Dirac Benchmarks

| Program & Configuration | Metric | Before (`WINDOWED_NS=0`) | **After (`WINDOWED_NS=1`)** | **Improvement** |
| :--- | :--- | ---: | ---: | :--- |
| **`hatom_sf_dirac`**<br>`Z=80, k=8, thresh=1e-6`<br>`L=0.5, rcut=1e-4` | Iterations (`rnorm < 1e-5`) | `16` (`7.08e-6`) | `16` (`6.61e-6`) | Identical convergence |
| | Converged Energy (Ha) | `-3531.844587` | `-3531.840623` | Finite-nuc (`rcut=1e-4`) |
| | **Total Wall Time** | `185.14 s` | **`28.13 s`** | **`6.58x` faster** (`8.43x` steady-state) |
| | `BSH_apply` Total / Steady | `35.51 s` (`2.64 s/it`) | **`1.96 s` (`0.107 s/it`)** | **`18.1x` / `24.7x` faster** |
| | `Coul_apply` Total / Steady | `73.67 s` (`5.71 s/it`) | **`5.71 s` (`0.362 s/it`)** | **`12.9x` / `15.8x` faster** |
| | `pVp` (`D(D(phi)*V)`) Total | `16.99 s` (`1.22 s/it`) | **`3.76 s` (`0.126 s/it`)** | **`4.52x` / `9.7x` faster** |
| | Raw `BSH` / `Coul` Leaves (Steady) | `50,576` / `64,856` | **`1,912` / `6,616`** | **`26.5x` / `9.8x` fewer raw leaves** |
| | Final `psi` / `phi` Leaves | `15,368` / `17,160` | **`1,992` / `1,544`** | **`7.7x` / `11.1x` smaller** |
| | **Peak RSS** | `9,796.5 MB` | **`3,634.6 MB`** | **`2.70x` lower (`-62.9%`)** |
| **`dirac-hatom` (Matched)**<br>`Z=80, k=8, thresh=1e-6`<br>`L=0.5, nucrad=1e-4, use_bsp=true` | Iterations (`maxiter=10`) | **Timed out in Iter 4 (`>600 s`)** | **`10 / 10` (`108.28 s`)** | **`>5.5x` total (`22.0x` at Iter 3)** |
| | Final Elec. Energy (Ha) | `-3531.7228` (at Iter 3) | **`-3531.747234`** (`res=5.7e-5`) | Smooth convergence |
| | Iter 3 Wall / `BSH_apply` | `233.40 s` / `113.31 s` | **`10.61 s` / `5.14 s`** | **`22.0x` / `22.1x` faster** |
| | 4-Spinor Leaf Sum (Peak obs.) | `344,992` (`d=15`, growing) | **`12,704` (`10,912` steady, `d=14`)** | **`31.6x` smaller vs steady** |
| | Raw `GVpsi` Leaf Sum (Iter 3) | `355,520` | **`34,360`** | **`10.3x` smaller** |
| | **Peak RSS** | `28,311.3 MB` (at Iter 3) | **`4,891.7 MB`** (all 10 iters) | **`5.79x` lower (`-82.7%`)** |
| **`dirac-hatom` (Deep Tree)**<br>`Z=80, k=7, thresh=1e-5`<br>`L=20, nucrad=1e-8` (`max_depth=30`) | Iterations (`maxiter=10`) | **Timed out in Iter 7 (`>600 s`)** | **`10 / 10` (`112.85 s`)** | **`>5.3x` total (`10.8x` at Iter 6)** |
| | Final Elec. Energy (Ha) | `-3532.191888` (at Iter 6) | **`-3532.191649`** (`res=7.7e-7`) | **`+0.44 mHa` vs exact point-DC** |
| | Iter 6 Wall / `BSH_apply` | `124.04 s` / `58.55 s` | **`11.50 s` / `4.44 s`** | **`10.8x` / `13.2x` faster** |
| | 4-Spinor Leaf Sum (Peak obs.) | `274,528` (`spinor[1]=87,432`) | **`17,376` (`14,304` steady, `d=23`)** | **`19.2x` smaller vs steady** |
| | Truncated `GVpsi` Leaf Sum | `155,488` (Iter 5) | **`11,744`** (steady) | **`13.2x` smaller** |
| | **Peak RSS** | `25,471.5 MB` (at Iter 6) | **`3,722.8 MB`** (all 10 iters) | **`6.84x` lower (`-85.4%`)** |

#### Per-Iteration Trajectory Highlights

1. **`hatom_sf_dirac` (`Z = 80, k = 8, thresh = 1e-6, L = 0.5, rcut = 1e-4`)**:
   Starting from identical over-refined guesses (`psi0 = 14,456` leaves at `d = 19`, `phi0 = 13,336` leaves at `d = 14`, `pVp_rhs = 7,752` leaves in Iter 0), Windowed NS collapses `psi` and `phi` by Iteration 2 down to their true physical support (**`1,800` and `1,480` leaves**), holding `phi_raw = 1,912` leaves (`BSH_apply = 0.104–0.108 s/iter`) and `pVp_rhs = 7,368–7,944` leaves (`Coul_apply = 0.341–0.389 s/iter`) flat through Iteration 15. Under Legacy NS, doubly-differentiating Iteration 0's `erf` boundary artifacts inflates `pVp_rhs` across iterations (`7,752` $\to$ `23,432` $\to$ `34,568` $\to$ `60,808` leaves), bloating `phi_raw` to `51,136` leaves (`26.7x` larger) and slowing steady-state iterations from `1.69 s` to `14.26–15.58 s`.

2. **`dirac-hatom` Matched (`k = 8, L = 0.5, nucrad = 1e-4, use_bsp = true`) & Deep-Tree Point-Like Nucleus (`k = 7, L = 20, nucrad = 1e-8`, `max_depth = 30`)**:
   - **Deep 30-Level `Key<3>` Execution**: At `L = 20` (`cell = 40 Bohr`) and `nucrad = 1e-8` (Gaussian nuclear exponent $\zeta = 9 \times 10^{16}$), `Hv` and `Vpsi` refine down to **`max_depth = 30`** (box width $40 \times 2^{-30} \approx 3.73 \times 10^{-8}\text{ Bohr}$), directly exercising `WindowedSeparatedConvolution`'s 64-bit `Key<3>` fallback path (`j > 20`) and overflow-free level-rescaled GC compression. Windowed NS completes all `10` iterations at `10.7–11.5 s/iter` (`3.72 GB` RSS), reaching `res = 7.71e-7` and `-3532.19164933 Ha` (`+4.44e-4 Ha` from exact point-nucleus Dirac–Coulomb).
   - **Analytically-Zero `spinor[1]` Diagnostic**: For the $1s_{1/2}$ ($m_j = +1/2$) ground state, $\psi_1$ (`spinor[1]`) is identically zero analytically; any non-zero tree nodes arise solely from numerical discretization/interface error coupled through $c\,\boldsymbol{\alpha}\cdot\mathbf{p}$. Under Legacy NS, `spinor[1]` explodes to **`87,432` leaves (`26.2M` complex coefficients)** in Iteration 5 (`k = 7`) and **`74,248` leaves (`38.0M` complex coefficients)** in Iteration 3 (`k = 8`) before hitting the `600 s` timeout. Under Windowed NS, `GVpsi_trunc[1]` stabilizes at **`72` leaves** (`norm = 2.4e-8` at `k = 8`, `4.9e-7` at `k = 7`) and `spinor[1]` remains flat at **`840` leaves (`k = 8`)** and **`1,992` leaves (`k = 7`, `d = 23`)** across all iterations.
