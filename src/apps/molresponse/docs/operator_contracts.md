# Operator Contracts — molresponse_v3

The mathematical definition of every response operation, the **metric /
sign convention** it uses, and the **legacy source line** that is the
ground truth for it.

## Why this file exists

`molresponse_v3` is a refactor of a working legacy implementation
(`src/apps/molresponse/`). The legacy code is the base truth. A refactor
is only safe if every operation still matches that truth — and the
dangerous bugs are the ones in *conventions* (a sign, a metric, a factor)
that compile fine and converge to a plausible-but-wrong number.

We shipped exactly such a bug: the RPA subspace overlap and the
response-vector normalization were written with the **Euclidean** metric
`⟨X|X⟩ + ⟨Y|Y⟩` instead of the **indefinite** RPA metric
`⟨X|X⟩ − ⟨Y|Y⟩`. It compiled, it "converged," and it passed a loose
ω-tolerance check — because nothing tied the code to the legacy formula
tightly and automatically.

**How to use this file:**
- When you touch an operation, check the code against the formula here,
  **and** check the formula here against the cited legacy line.
- When you add an operation, add a row with its ground-truth line.
- If the code and the cited legacy line disagree, the code is wrong until
  proven otherwise — legacy is the reference.

Notation: `X` = excitation block, `Y` = de-excitation block. `φ` =
occupied ground-state orbitals. `⟨a|b⟩` over vecfuncs means
`Σ_p ⟨a_p|b_p⟩` (sum over orbitals). Closed-shell unless a β block is
named.

---

## Metrics and inner products

| Operation | Formula | Convention | Ground truth (legacy) | v3 location |
|---|---|---|---|---|
| Euclidean bundle inner | `⟨a\|b⟩ = Σ_p ⟨a_p\|b_p⟩` | positive-definite | `response_functions.h` `response_space_inner` | `rs::inner` (response_space_ops.hpp) |
| **RPA subspace overlap `S`** | `S_ij = ⟨Xi\|Xj⟩ − ⟨Yi\|Yj⟩` | **INDEFINITE** (Y block negative) | `ExcitedResponse.cpp:871` | `rs::metric` |
| **Response-vector norm** | `‖v‖ = √(⟨X\|X⟩ − ⟨Y\|Y⟩)` | **INDEFINITE** | `ResponseBase.cpp:2097-2099` | `rs::metric_inner` |
| TDA norm (x-only) | `‖v‖ = √(⟨X\|X⟩)` | standard | `ResponseBase.cpp:2077-2088` | `rs::metric_inner` (no Y block) |

The single most error-prone convention in the whole module: **the Y
(de-excitation) block enters the overlap and the norm with a minus
sign.** For x-only (TDA) states there is no Y block, so the indefinite
metric collapses to the standard `⟨X|X⟩` — `rs::metric` and
`rs::metric_inner` are written so this falls out automatically (`if
constexpr` on block presence).

---

## Subspace eigenproblem

| Operation | Formula | Convention | Ground truth | v3 location |
|---|---|---|---|---|
| subspace matrix `A` | `A_ij = ⟨Xi\|Λxj⟩ + ⟨Yi\|Λyj⟩` | Euclidean (+) | `ExcitedResponse.cpp:885` | `ESSolver::step_*` (`rs::inner(roots, lambda)`) |
| generalized eig | `A U = S U Ω` | `S` indefinite, solved with `sygvp` | `ExcitedResponse.cpp:1030` (`excited_eig`) | `rs::diagonalize` |
| eigenvalue ordering | dominance-swap → phase fix → cluster-unmix → **ascending sort** | stable ascending | `ExcitedResponse.cpp:1046-1122` | `rs::diagonalize` steps 3-5.5 |

Note the asymmetry: `A` uses the **Euclidean** (+) inner, only `S` uses
the **indefinite** metric. `sygvp` (symmetric-definite) is valid because
`S = ⟨X|X⟩ − ⟨Y|Y⟩` stays positive-definite for excitation-dominated
states (‖X‖ > ‖Y‖); it only loses definiteness at a triplet instability.
The **final ascending sort** is required for slot-stable KAIN history —
without it the greedy dominance swap flips slot identity between
iterations for near-degenerate roots (`ExcitedResponse.cpp:1122`
`sort_eigenvalues`).

---

## Per-root building blocks (kernels)

| Operation | Formula | Notes | v3 location |
|---|---|---|---|
| response density | `ρ = 2 · Σ_p φ_p (x_p + y_p)` | factor 2 = spin; `(x+y)` folds X,Y | `Kernels<Full,ClosedShell>::compute_density` (full.hpp:44) |
| Coulomb-exchange `γ` | `γ_X = J[ρ]φ − c_xc(K[φ,x]φ + K[y,φ]φ)` | **cross-channel** K[y,φ] couples X↔Y | `compute_gamma` (full.hpp); legacy `ResponseKernel.hpp:113` |
| `V0·x` | `V_local·x − c_xc K[φ,φ]·x` | acts on X and Y independently | `compute_V0x` |
| `E0·x` (off-diag) | `F_offdiag · x` | diag ε absorbed into BSH | `compute_E0x` |
| `E0·x` (full) | `F · x` | full Fock; for Λ assembly | `compute_E0x_full` |
| `T0·x` | `−½∇²·x` | kinetic | `compute_T0x` |
| **θ** (BSH driver) | `θ = V0x − E0x + γ` | FD adds `+ V_p` source | `assemble_theta` (assembly.hpp); legacy `ResponseKernel.hpp:271` |
| **Λ** (subspace) | `Λ = T0x + V0x − E0x_full + γ` | full Fock (not off-diag) | `assemble_lambda` (assembly.hpp); legacy `ResponseKernel.hpp:397` |
| BSH apply | `x_new = Q(BSH(ω)·(−2(θ + shift·x)))` | paired ±ω for Full (X:+ω, Y:−ω) | `bsh_apply` |

Density factor convention: `spin_factor × y_factor`. Restricted Static =
4, restricted Full = 2, unrestricted Static = 2, unrestricted Full = 1.
`alpha_factor = −(spin_factor × y_factor)`.

---

## Tensor-layer exchange — `exchange_ctx` (FD `--fd-tensor` path)

The contract for `kernels/exchange_ctx.hpp` (doc 28; framework docs 26/27).
This is an **alternate assembly of θ** that is bit-for-bit the *same math* as the
per-op reference (`compute_V0x`/`compute_gamma`/`compute_E0x`) but builds the
shared convolution tensors **once** and assembles θ as cheap contractions. It is
the interface `parallel-runtime` distributes (the batched waves are taskq/subworld
work) and `perf-model` meters (the named phases below). Gated behind
`policy_.exchange_tensor` / `--fd-tensor`; **gate 0 (the per-op path) is the
reference oracle and stays untouched.** FD ClosedShell Static + Full only;
ES (`assemble_lambda`, +T0x) and VBC/β reuse the same layer but are not wired yet.

**Unified form** (verified `exchangeoperator.cc:182-187`, doc 27 §1):
`K[a,b](c)_k = Σ_i a_i · g(b,c)_{ik}`, with `g(b,c)_{ik} = Poisson(b_i·c_k)`
(`a`=ket/weight, `b`=bra, `c`=apply_to). The convolution `g(b,c)` (n² Poisson) is
the cost; the contraction `Σ_i a_i·g` is a cheap mul+reduce.

**Convolution tensors** (row-major `T[i*n+k]`, `n = |φ|`):

| tensor | definition | class | lifetime | v3 location |
|---|---|---|---|---|
| `g0` | `Poisson(φ_i·φ_k)` | 1 — φ-only ⇒ cacheable | **per protocol**, cached on `ResponseGroundState::g0_alpha` (built in `fd_solver::step` when `--fd-tensor`+empty; fresh-empty per protocol via `build_gs`) | `build_g0` |
| `Tx` | `Poisson(φ_i·x_k)` | 2 — has response ⇒ transient | per response-iter | `build_pair_tensor` (via `build_ctx_*`) |
| `Ty` | `Poisson(φ_i·y_k)` | 2 — transient (Full only) | per response-iter | `build_pair_tensor` (via `build_ctx_full_cs`) |
| `J` | `Poisson(ρ)` (Coulomb/Hartree) | — | per response-iter | `build_ctx_*` |

**Sharing insight** (doc 27 §4): `g(x,φ) = T_xᵀ`, so one `Tx` build serves BOTH
the V0x ground-exchange (term 1, `contract_col`) AND the γ cross term (term 5,
`contract_row`). Likewise `Ty` serves terms 2 and 6. Static Class-2 convs: 2n²→n²;
Full: 4n²→2n².

**Contractions** (cheap; NO convolutions):

| helper | reduction | assembles |
|---|---|---|
| `contract_col(a, T, n)` | `out_k = Σ_i a_i·T[i*n+k]` | groundK `K[φ,φ](x)`=`col(φ,Tx)`; γ-direct `K[φ,x](φ)`=`col(x,g0)` |
| `contract_row(a, T, n)` | `out_k = Σ_i a_i·T[k*n+i]` (transpose) | γ-cross `K[x,φ](φ)`=`row(φ,Tx)`; Full γ_X `K[y,φ](φ)`=`row(φ,Ty)` |

**θ map** (must mirror the reference *exactly*):
- Static: `θ = [V_local·x − c_xc·col(φ,Tx)] − E0x + Q( J·φ − c_xc·(col(x,g0) + row(φ,Tx)) )`
- Full γ_X = `col(x,g0) + row(φ,Ty)`; γ_Y = `col(y,g0) + row(φ,Tx)` (X/Y V0x use Tx/Ty resp.)

**Q discipline (the easiest way to break bit-identity, doc 28 §2):** `Q` (`gs.Qa`)
applies to the **γ block ONLY**. `V0x` and `−E0x` stay un-projected; `E0x` is reused
from the reference kernel (`Kernels<…>::compute_E0x` — call, do not modify).
Truncate intermediates at `vtol = thresh·0.1`, final θ at `thresh`.

**Interface (the batching/perf contract surface):**
```
build_g0(world, gs, vtol) -> vecfuncT                              // Class 1, cache on gs.g0_alpha
build_ctx_static_cs / build_ctx_full_cs(world, gs, state, rho, vtol) -> ResponseExchangeCtx{Tx,[Ty],J}
assemble_theta_static_cs / assemble_theta_full_cs(world, gs, state, g0, ctx) -> State
assemble_theta_tensor(world, gs, state, rho)                       // one-call gate-1 entry; overloaded on State
```

**A/B floor:** explicit-Poisson + dot-contraction vs MADNESS's `Exchange` operator
is the same math on a different truncation/accumulation path → agreement is
**~1e-3 RELATIVE** on a converged channel (`verify_fd_tensor.sh`, `TOL=2e-3`,
convergence-gated), NOT machine-eps. A larger gap on a converged channel is a real
bug — check Q discipline / truncation order.

**Perf meters** (Inc-3d; perf-model contract — named `PROFILE_BLOCK`s in
`exchange_ctx.hpp`, no-op unless `WORLD_PROFILE_ENABLE`, flowing into perf-model's
WorldProfile JSON via the already-wired `dump_json` call sites; see "Performance
profile schema"). Coarse phase blocks — the core ops (`apply`/`multiply`/`dot`/
`truncate`, already PROFILE-instrumented) roll up *within* each via inclusive time:
`rs_ext_g0_build` (φ·φ tensor, once/protocol), `rs_ext_ctx_build` (per-iter J +
Tx[/Ty] convolutions — dominant; `apply` count under it = # Poisson waves, rises
with `--fd-tensor-tile`), `rs_ext_assemble` (per-iter contractions + E0x + Q).

---

## RPA symmetric reduction — REMOVED (2026-06)

The `ESSolverFullRPA` symmetric-reduction solver (`(A−B)(A+B)u = ω²u`,
`u = X+Y`) and its `apply_AplusB` / `apply_AminusB` operators were
removed — not a direction we are pursuing. Full ES is the direct
paired-(X,Y) `ESSolver<Full, ClosedShell>` only.

---

## Top-of-iteration discipline

| Operation | Formula | Ground truth | v3 location |
|---|---|---|---|
| project | `Q·v` per spin (alpha→Qa, beta→Qb) on every block | `ExcitedResponse.cpp:2483-2517` | `rs::project` |
| orthonormalize | modified Gram-Schmidt in the **RPA metric** | `ResponseBase.cpp:2091` (normalize) | `rs::orthonormalize` |

`orthonormalize` MUST use `metric_inner` (indefinite), not the Euclidean
inner — so the bundle entering `diagonalize` is orthonormal in the same
metric the subspace overlap `S` uses (→ `S ≈ I`). This is the line that
had the bug.

---

## How this is enforced automatically

A contract table is a human aid; the gate is the test. See
`tests/test_v3_es_skeleton.cpp`:

1. **Convergence gate** — a run only PASSES if it actually converged
   (`!diverged` and `max_residual ≤ solver target`), not merely if ω
   landed near the reference. The metric bug stalled at a residual floor;
   this gate fails that.
2. **Legacy ω comparison** — converged ω must match the tabulated legacy
   value within `tol`. Legacy is an independent implementation, so a
   match is real signal.
3. **Two-solver cross-check** — retired with `ESSolverFullRPA` (removed
   2026-06). The remaining independent check is item 2 (legacy ω).

When adding an operation, add its row above **and** a check that ties it
to legacy — the table alone did not prevent the metric bug; the missing
piece was the automated tight comparison.

---

## Performance profile schema (v1)

Pinned by the **perf-model** thread (design: `docs/29_perf_model_design.md`).
The cross-thread contract: the machine-readable profile that
`WorldProfile::dump_json` emits, that `exchange` reports its Tx/tile counts +
phase timings into, and that `parallel-runtime` reads to settle the doc-24-vs-25
fork with measured numbers. **Stable key names — coordinate changes here.**

Emitted once per run by rank 0, **only** when built with `ENABLE_WORLD_PROFILE=ON`
**and** env `MADQC_PROFILE_JSON` is set (else absent — zero-effect contract).
Every per-phase statistic is a faithful copy of a reduced `WorldProfileEntry`
field: `{sum, min, max, pmin, pmax}` (`pmin`/`pmax` = the rank holding the min/max).

```jsonc
{
  "schema_version": 1,
  "world_size": 6,                 // P (MPI ranks)
  "total_cpu_s": 812.4,            // rank-0 cpu_time - WorldProfile::cpu_start
  "total_wall_s": 141.2,           // rank-0 wall_time - WorldProfile::wall_start
  "overhead_s_per_call": 3.0e-7,   // estimated profiling overhead/call
  "context": {                     // OPTIONAL; filled by the v3 caller, null in
                                   // core. The join key for the cost-model fit.
    "molecule": "h2o", "n_occ": 5, "k": 6, "thresh": 1e-4,
    "protocol": 0, "box_L": 30.0
  },
  "phases": [
    {
      "name": "FunctionImpl::apply",   // raw __FUNCTION__ / class::function key
      "phase": "apply",                // canonical taxonomy (PM-2); else "other"
      "count":      {"sum": 1.2e5, "min": 1.9e4, "max": 2.1e4, "pmin": 3, "pmax": 0},
      "cpu_excl_s": {"sum": 402.1, "min": 61.0, "max": 71.3, "pmin": 4, "pmax": 0},
      "cpu_incl_s": {"sum": 588.0, "min": 90.2, "max": 99.8, "pmin": 4, "pmax": 0},
      "nmsg_sent_excl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nmsg_sent_incl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nbyte_sent_excl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nbyte_sent_incl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nmsg_recv_excl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nmsg_recv_incl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nbyte_recv_excl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0},
      "nbyte_recv_incl": {"sum": 0, "min": 0, "max": 0, "pmin": 0, "pmax": 0}
    }
    // … one object per registered profile entry
  ]
}
```

Canonical `phase` values (PM-2): `apply`, `compress`, `reconstruct`, `multiply`,
`inner`, `gaxpy`, `truncate`, `exchange` (γ build, block `rs_exchange_gamma`),
`projection` (`Q·v`, block `rs_projection`), `other`. The exchange thread's
Tx/tile counters aggregate under the `rs_exchange_gamma` block so they land in the
`exchange` phase. Cost-model consumers (§9 of the companion doc): use `cpu_excl_s`
+ `count` for `T_compute`, `nmsg_*`/`nbyte_*` for `T_comm`, and `max/sum` ratios
for the imbalance factor φ; recover wall by joining `context` with the coarse
`StateMetrics.wall_s` layer.
