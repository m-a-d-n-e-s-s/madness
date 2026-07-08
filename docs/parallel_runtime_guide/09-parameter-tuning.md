# Chapter 9 — Parameter tuning, per operation

[← Operations](08-operations.md) · [Index](README.md) · [Next: Parallel function I/O →](10-io-and-hdf5.md)

This chapter is the practical payoff: what every knob does, the symptom that tells
you to touch it, and which knobs matter for each operation. Parameters fall into
three groups: **runtime** (environment variables — communication & threading),
**build** (CMake), and **numerical** (`FunctionDefaults` — accuracy vs cost &
distribution).

> Golden rule: change one knob at a time, measure with the counters in
> [Chapter 3 §3.5](03-rmi-thread.md#35-statistics-you-can-read) and
> [`tree_size()`](07-function-and-tree.md#75-how-big-is-the-tree), and keep the
> known-good small cases (H2O, CH3OH) as regression anchors.

---

## 9.1 The knobs and what they do

### Runtime (environment variables)

| Variable | Raise it when… | Lower it when… | Default |
|----------|----------------|----------------|---------|
| `MAD_NUM_THREADS` | cores idle, compute-bound op | oversubscription, memory pressure (each thread has working set) | #HW cores |
| `MAD_SEND_BUFFERS` (`nsend`) | a rank bursts many small remote tasks/AMs and stalls in 100 µs send waits | memory tight (each buffer ~`AmArg`-sized, but they accumulate) | 128 / 512 CrayXT |
| `MAD_RECV_BUFFERS` (`nrecv_`) | many senders target one rank (incast); receiver drops to huge-msg path | memory tight | 128 |
| `MAD_BUFFER_SIZE` (`max_msg_len_`) | coefficient tensors exceed it → rendezvous on every big message | shrink only to save memory on tiny-message workloads | 1.5 MB |
| `MAD_NSSEND` | a fast sender overruns a slow receiver (raise = less back-pressure); set small to throttle harder | — | = `nrecv_` |
| `MAD_BACKOFF_US` | reduce to cut receive latency for latency-bound, message-heavy phases | raise to cut RMI-thread CPU spin when communication is sparse | 2 µs |
| `MAD_MAX_REDUCEBCAST_MSG_SIZE` | rarely; lower to bound peak collective message size on fragile networks | — | `INT_MAX` |
| `MAD_BIND` | NUMA effects suspected; pin worker threads | — | platform |

### Build (CMake)

| Option | Use |
|--------|-----|
| `ENABLE_NEVER_SPIN` | debugging hangs / oversubscribed nodes — waits sleep instead of spin |
| `MADNESS_TASK_BACKEND` | `Pthreads` (default) vs `OneTBB` vs `PaRSEC`; TBB/PaRSEC change thread accounting (Ch.1 §1.1) |
| `MADNESS_DQ_USE_PREBUF` | high task-submission rates — fewer queue locks |
| `MADNESS_DQ_STATS` / `MADNESS_DQ_PRINT` | diagnose queue growth / contention |

### Numerical (`FunctionDefaults<NDIM>`, `funcdefaults.h`)

| Knob | Effect | Cost coupling |
|------|--------|---------------|
| `k` | wavelet order; tensor size `k^d` | per-node cost `∝ k^{d+1}`; memory `∝ k^d` |
| `thresh` | accuracy; sets where wavelets are kept | controls `N_leaf` (tree size) |
| `initial_level` | starting projection depth | startup work |
| `max_refine_level` | hard depth cap | bounds worst-case `N` |
| `special_level` | forced refinement at special points | local `N` near features |
| `truncate_mode` | level scaling of `truncate_tol` | how aggressively `truncate` prunes |
| `autorefine` | refine during `multiply` | accuracy vs more nodes |
| tensor type (`TT_FULL`/`TT_2D`/…) | full vs low-rank storage | `apply`/`gaxpy` cost via rank `r` |
| pmap | data distribution | imbalance `φ` and cross-rank edge count |

---

## 9.2 Per-operation tuning playbook

For each op: the dominant cost, the knobs that move it, and the symptom→action
mapping. (Costs from [Chapter 8](08-operations.md).)

### compress / reconstruct — `Θ(N·d·k^{d+1})`, edge-traffic + 1 fence

- **Dominant cost:** the filter/unfilter, and cross-rank parent↔child messages.
- **Primary knob:** the **pmap**. `SimplePmap` makes almost every edge remote;
  `LevelPmap` co-locates parent/child and cuts edge traffic sharply for these
  sweeps. If you compress/reconstruct often, prefer `LevelPmap` or an
  `LBDeuxPmap` built to keep subtrees together.
- **Secondary:** `MAD_BUFFER_SIZE` so child blocks (`(2k)^d` doubles) ride inline,
  not via rendezvous. For `k=10, d=3`: `8·20³ = 64 KB` — well under default, fine.
- **Symptom → action:** many small AMs + high `nmsg_sent` during state changes →
  switch pmap toward locality; long fence tail → reduce how often you change state.

### gaxpy / +,− — `Θ(N·r·k^d)`, `O(N)` AMs + 1 fence

- **Dominant cost:** the `send` AMs when operands have **different pmaps**.
- **Primary knob:** ensure `f` and `g` **share a pmap** (same `FunctionDefaults`),
  making the combine local. Use a low-rank tensor type to shrink `r·k^d`.
- **Symptom → action:** gaxpy unexpectedly chatty → check the two operands' pmaps
  match; if not, project/redistribute one to align.

### multiply — `Θ(N·d·k^{d+1})`, no comm + 1 fence

- **Dominant cost:** the two separable transforms per node; node count via
  `autorefine`.
- **Primary knobs:** `autorefine` (off → fewer nodes, lower accuracy; on → accurate,
  more nodes), `k`, and `thresh` of the result. Compute-bound, so
  **`MAD_NUM_THREADS`** matters most for wall time.
- **Symptom → action:** memory spike during multiply → `autorefine` is inflating
  the tree; truncate immediately after, or accept lower `k` for the product.

### inner — `Θ(N·k^d)` + 1 all-reduce

- **Dominant cost:** local trace; communication is a single reduce.
- **Primary knob:** **load balance** (`φ`) — the reduce waits for the slowest rank,
  so imbalance shows up directly. `MAD_NUM_THREADS` for the local traces.
- **Symptom → action:** reduce-bound at scale → improve the pmap balance; batch
  many inner products (e.g. a Gram matrix) to amortize the reduce.

### apply — `Θ(N·D_eff·M·d·k^{d+1})`, heavy target traffic + 1 fence

- **Dominant cost:** the most expensive op — both flops and network.
- **Primary knobs:**
  - **operator `thresh`** → sets `M` (number of separated terms) and the screening
    tolerance → directly scales both flops and `D_eff`. Loosen within accuracy
    budget.
  - **tensor type / low-rank** → replaces a `k^{d+1}` factor with `M·k·r`; large win
    at high `k`.
  - **pmap locality** → target nodes `K+δ` are nearby in index space; a
    locality-preserving pmap (LevelPmap/LBDeux) keeps many `accumulate2` tasks
    local, cutting the heaviest traffic.
  - **`MAD_SEND_BUFFERS`** → apply fans out many target tasks per source; raise if
    senders stall.
- **Symptom → action:** OOM right after apply → low-rank tensor type + truncate;
  network-bound at scale → locality pmap + raise send/recv buffers; too slow but
  memory OK → loosen operator `thresh`, lower `k` at coarse protocols.

### truncate — `Θ(N·d·k^{d+1})`, light comm + **2 fences**

- **Dominant cost:** two sweeps and two fences.
- **Primary knob:** `truncate_mode` / `truncate_tol` — more aggressive pruning →
  smaller `N_leaf` downstream (helps every later op and memory), at some accuracy
  cost. Call truncate **after** tree-growing ops (multiply, apply) to claw back
  memory.
- **Symptom → action:** memory grows across iterations → truncate more often /
  more aggressively; fence overhead noticeable at large `P` → don't truncate
  needlessly between every op.

---

## 9.3 Cross-cutting strategies

1. **Fewer fences.** Each operation ends in a fence (`O(log P)` per pass, more with
   stragglers). Chain operations that don't need an intervening global view, and
   prefer vector APIs (`vmra.h`) that fence once for many functions instead of once
   per function.
2. **Locality vs balance is the pmap trade-off.** `SimplePmap` = perfect balance,
   worst locality; `LBDeuxPmap` with a good cost functor = both, if the cost
   estimate is accurate. Tune the cost functor to the dominant op (coefficient
   count for memory; flop weight for apply-bound runs).
3. **Memory cliff at `k=10`.** Per-node bytes `8·k^d` and `N_leaf` both jump ~4.6×
   from `k=6→10`. Combined with ground-state **replication** across subworld tasks
   (`n_occ·N_leaf·8·k³`, independent of `P`), this is the documented OOM. Levers:
   node-local shared replica (`WorldDCNodeReplicatedPmap`), on-demand orbital
   loading (`on_demand` state), lower `--ntasks-per-node`. See companion doc §11.
4. **Overlap, don't serialize.** Keep worker threads busy with independent tasks
   while the single RMI thread moves data; avoid collectives inside hot tasks
   (they take the MPI mutex and stall a worker).

---

## 9.4 A tuning loop you can run

1. Pick one operation and a representative input (e.g. one BSH `apply` on a known
   molecule).
2. Record baseline: `RMI::stats` (msgs/bytes), `tree_size()` and per-rank
   `coeffs.size()` (⇒ `φ`), wall time, and MaxRSS.
3. Change **one** knob from §9.2 for that op.
4. Re-measure; keep it only if the target metric improved **and** H2O/CH3OH still
   pass.
5. Move to the next dominant cost (use the table in
   [Chapter 8 §8.8](08-operations.md#88-at-a-glance-cost--communication) to decide
   what's dominant).

This is exactly the calibration loop the companion performance-models doc
formalizes (fit `R, α, β, φ` on small cases, predict large ones).

[← Operations](08-operations.md) · [Index](README.md) · [Next: Parallel function I/O →](10-io-and-hdf5.md)
