# MADNESS Parallel Runtime and Performance Models

A ground-up, source-referenced guide to the MADWorld distributed runtime, how
multiresolution `Function` objects are distributed and operated on, and how to
build quantitative performance (time / communication / memory) models for the
core operations.

All file:line references point into `src/madness/world/` (the runtime) and
`src/madness/mra/` (the multiresolution layer) of this repository. Line numbers
are accurate as of the time of writing; treat them as anchors, not contracts.

> **Companion:** for the narrative, diagram-driven walkthrough of each subsystem
> and a per-operation parameter-tuning playbook, see the
> [Parallel Runtime Developer Guide](parallel_runtime_guide/README.md). This
> document is the quantitative companion to that guide.

> Scaling note: This document is written with the project scaling goal in mind
> (extend the response solver beyond 20 occupied orbitals). The memory model in
> §11 derives the `n_occupied × n_leaves × k³ × 8 bytes` per-task figure that
> drives the OOM behavior recorded in `CLAUDE.md`.

---

## 0. Mental model in one page

MADNESS is built in layers. From the bottom up:

```
        ┌─────────────────────────────────────────────────────────┐
  MRA   │ Function<T,NDIM>  (handle)  ──►  FunctionImpl<T,NDIM>     │
        │   adaptive 2^NDIM-tree of Key ──► FunctionNode (k^NDIM    │
        │   coefficient tensors); operations = tree sweeps          │
        ├─────────────────────────────────────────────────────────┤
        │ WorldContainer<Key,FunctionNode>  = distributed hash map  │
        │   pmap: Key ─► owner ProcessID   (this IS the data        │
        │   distribution)                                           │
  WORLD ├─────────────────────────────────────────────────────────┤
        │ WorldTaskQueue + Future<T> + DependencyInterface          │
        │   dataflow tasks; futures wire producers to consumers     │
        ├─────────────────────────────────────────────────────────┤
        │ WorldGOP (collectives, fence) | WorldAM (active messages) │
        ├─────────────────────────────────────────────────────────┤
        │ RMI server thread  +  SafeMPI (thread-safe MPI wrapper)   │
        └─────────────────────────────────────────────────────────┘
```

Five ideas carry almost all of the behavior:

1. **A `World` wraps an MPI communicator** and owns four services: `mpi`, `am`
   (active messages), `taskq` (task queue), `gop` (collectives). Subworlds are
   sub-communicators with their own `World`.
2. **Everything distributed is a `WorldContainer`** — a hash map whose keys are
   assigned to processes by a *process map* (`pmap`). For a `Function`, the
   container holds the tree nodes, so **the pmap literally defines how a function
   is spread across ranks.**
3. **Work is expressed as dataflow tasks.** A `TaskFn` holds a callable plus
   argument `Future`s; it becomes runnable when all its futures are set. Tasks
   may run locally or be shipped to the owner of a key.
4. **A single RMI thread per process owns the wire.** All point-to-point traffic
   is active messages sent through one communication thread; worker threads never
   call MPI directly for AM.
5. **`fence()` is global quiescence detection** — it terminates only when every
   task is done and the global count of messages-sent equals messages-received.
   Operations are separated by fences.

---

## 1. The `World` object and MPI substrate

### 1.1 `World`

`World` (`world.h:119-608`) is the central handle. It is constructed from a
`SafeMPI::Intracomm` (`world.cc:81-105`) and owns, in initialization order:

- `WorldMpiInterface& mpi` (`world.h:204`) — the MPI wrapper.
- `WorldAmInterface& am` (`world.h:205`) — active-message engine.
- `WorldTaskQueue& taskq` (`world.h:206`) — the dataflow task queue.
- `WorldGopInterface& gop` (`world.h:207`) — collective operations.

Each `World` has a universe-unique `unsigned long _id` (`world.h:187`) assigned
on rank 0 and broadcast (`world.cc:91-97`). A registry maps object ids to
pointers (`map_id_to_ptr`, `world.h:182`) so that an active message naming a
`(world_id, object_id)` pair can be dispatched to the right C++ object on the
receiver — this is how `WorldObject`-derived classes (containers, function impls)
receive remote method calls.

**Subworlds** are created by splitting the communicator
(`SafeMPI::Intracomm::Split` / `Create`, `safempi.h:616-654`); each gives a new
`World`. This is the mechanism behind MacroTaskQ subworld pools and the
`state_parallel` scheduling in molresponse_v2.

### 1.2 SafeMPI — the thread-safe MPI wrapper

`safempi.h` wraps MPI so it can be called under `MPI_THREAD_SERIALIZED`. Every
raw MPI call is guarded by a global mutex via the `SAFE_MPI_GLOBAL_MUTEX` macro
(`safempi.h:96-99`), i.e. a scoped lock on `SafeMPI::charon`. Practical
consequences for modeling:

- **MPI is effectively serialized across threads in a process.** Only one thread
  holds the MPI lock at a time. Concurrency comes from overlapping computation
  with the *single* RMI thread doing communication, not from many threads issuing
  MPI concurrently.
- Core primitives used: non-blocking `Isend`/`Irecv`/`Issend`
  (`safempi.h:731-753`), synchronous `Send` (`MPI_Ssend`, `safempi.h:755-759`),
  `Bcast`, `Reduce`/`Allreduce` (`safempi.h:782-798`), and request polling via
  `Test`/`Testsome`/`Testany` (`safempi.h:296-431`).

**Tag plan** (`safempi.h:107-113, 533-564`):

| Tag range | Use |
|-----------|-----|
| 1–999 | one-shot reserved tags (`unique_reserved_tag`) |
| 1000 | `DEFAULT_SEND_RECV_TAG` |
| 1001 | `MPIAR_TAG` (archive I/O) |
| 1023 | `RMI_TAG` — **all** RMI active messages |
| 1024–4095 | round-robin `unique_tag()` (period 3072) |
| 4096–8191 | huge-message exchange |

All routine point-to-point AM traffic rides a *single tag* (1023) with
`MPI_ANY_SOURCE` receives; ordering and dispatch are handled above MPI (see §3).

---

## 2. Active Messages (`WorldAM`)

An **active message** is a one-way message that names a handler function to run
on the receiver. This is the universal transport for "do X on rank r."

**`AmArg`** (`worldam.h:80-158`) is the message envelope: a fixed header
(`RMI::HEADER_LEN = 64` bytes) followed by `nbyte` of user payload, plus
bookkeeping (`worldid`, relative function pointer `func`, source rank, flags).
Total fixed overhead is ~96 bytes on a 64-bit machine. Payload is allocated
contiguously after the header by `alloc_am_arg` (`worldam.h:162-167`).

Handler signature: `typedef void (*am_handlerT)(const AmArg&)` (`worldam.h:77`).

**Send path** (`worldam.h:278-351`):

- A pool of `nsend` managed send buffers (default 128; 512 on Cray; `worldam.h:213-216`)
  provides **flow control**. A sender claims the next buffer round-robin
  (`cur_msg`), and if its previous MPI request hasn't completed it spins on
  `TestAndFree`, sleeping 100 µs between checks (`worldam.h:338-346`) to throttle
  injection. So the *maximum number of in-flight AMs from one rank is bounded by
  `nsend`* — a real parameter in any latency/throughput model.
- The actual enqueue calls `RMI::isend(arg, size, dest, handler, attr)`
  (`worldam.h:349`). Message size on the wire = `sizeof(AmArg) + payload`.

**Counters for quiescence:** `nsent` is incremented on send; `nrecv` is
incremented *after* the handler finishes running (`worldam.h:266`). `fence()`
compares global sums of these (see §4).

---

## 3. The RMI communication thread

Each process runs **one dedicated RMI thread** (`worldrmi.h`, `worldrmi.cc`)
that owns all AM traffic. Worker threads hand messages to it; it never blocks the
compute threads on MPI.

- **Posted receives:** `nrecv_` buffers (default 128, env `MAD_RECV_BUFFERS`),
  each `max_msg_len_` bytes (default 1.5 MB, env `MAD_BUFFER_SIZE`), all posted as
  `Irecv` on `RMI_TAG` with `MPI_ANY_SOURCE` (`worldrmi.cc:227-322`).
- **Server loop `process_some()`** (`worldrmi.cc:61-178`): polls up to `maxq_`
  requests with `Testsome`, up to 1000 spins, then backs off `MAD_BACKOFF_US`
  µs (default 2). For each arrived message it reads the handler pointer and
  attributes from the header and invokes the handler.
- **Ordering:** the top 16 bits of the attribute word carry a per-destination
  sequence counter. `ATTR_ORDERED` messages from a given source are delivered to
  their handlers strictly in order; out-of-order arrivals are queued and sorted
  (`worldrmi.cc:105-166`). The sequence counter occupies the high 16 bits of the
  attribute word, so up to ~2^16 ≈ 65k ordered messages per source before wrap.
- **Huge messages** (> `max_msg_len_`) use a rendezvous protocol on tags
  4096–8191: a control AM, an ack, then the bulk transfer
  (`worldrmi.cc:180-199, 325-343`).

**Modeling takeaway:** point-to-point latency has three additive parts —
(1) sender-side buffer wait (only under congestion), (2) MPI transit, and
(3) receiver-side polling latency (up to one backoff interval + handler time).
Bandwidth per link is bounded by the buffer size and the single-threaded,
mutex-serialized MPI.

---

## 4. Collectives and `fence()` (`WorldGOP`)

`WorldGopInterface` (`worldgop.h`) provides the collectives. All use a **binary
tree** over ranks (`binary_tree_info`, `safempi.h:863`):

- **Broadcast** (`worldgop.cc:173-208`): root → children down the tree; each node
  `Irecv`s from parent and `Isend`s to ≤2 children. Depth `O(log P)`, total `O(P)`
  messages, chunked by `MAD_MAX_REDUCEBCAST_MSG_SIZE`.
- **Reduce / Allreduce** (`worldgop.h:784-868`): up-tree combine then broadcast
  down. Depth `O(log P)`, `~2P` messages, element-wise op applied at each node.
- **Barrier** (`worldgop.h:702-706`): implemented as a global sum of ranks with a
  consistency check.

**`fence()`** (`worldgop.cc:50-159`) is the cornerstone of correctness and a
frequent serialization point. It performs **Dijkstra-style termination
detection** on the binary tree:

1. `taskq.fence()` waits for all local tasks to drain.
2. It reads (twice, for consistency) the local `taskq.size()`, `am.nsent`,
   `am.nrecv`.
3. Subtree sums of `nsent`/`nrecv` flow up the tree; the global pair is broadcast
   down.
4. **Termination** requires `global_nsent == global_nrecv` *and* the value is
   unchanged from the previous pass (`worldgop.cc:128-132`). Otherwise it loops.

Cost: `O(log P)` messages of 16 bytes per pass; **a clean program converges in
~2 passes**, but any straggler task that emits new messages forces another round.
Each fence is a global synchronization barrier — minimizing fences is a primary
optimization lever.

---

## 5. Tasks, Futures, and the thread pool

### 5.1 `Future<T>`

A `Future<T>` (`future.h:368-712`) is a shallow handle to a `FutureImpl<T>`
(`future.h:74-352`):

- `std::atomic<bool> assigned` (`future.h:95`) gives a **lock-free `probe()`**.
- A small inline stack of up to 4 callbacks (`future.h:83`) avoids heap allocation
  for the common case.
- For an already-assigned value the wrapper skips the `shared_ptr` allocation
  entirely (`future.h:397-409`) — cheap when futures are set on creation.
- **Remote futures:** `remote_ref(World&)` (`future.h:671-677`) produces a
  `RemoteReference<FutureImpl<T>>`. When such a future is `set`, the value is
  shipped via an AM (`set_handler`, `future.h:106-141, 253-264`) to the process
  holding the original future. This is how a result computed on rank B reaches the
  caller on rank A.

### 5.2 Dependencies and `TaskFn`

`DependencyInterface` (`dependency_interface.h:100-380`) holds an atomic
`ndepend` counter (`:105`). A task is runnable when `ndepend == 0`. `inc()`/`dec()`
adjust it; reaching zero fires the registered callbacks outside the lock
(`:262-284`).

`TaskFn` (`taskfn.h:473-861`) bundles a callable with up to 9 arguments and a
result `Future`. On construction, `check_dependencies()` (`:548-602`) inspects
each argument: for every argument that is an unset `Future`, it calls `inc()` and
registers itself as a callback. When all argument futures are set, `ndepend`
hits 0 and the task is submitted to the pool. `run()` (`:857-860`) unwraps the
arguments, calls the function, and sets the result future.

**This is the dataflow engine:** you can submit a task whose inputs are the
not-yet-computed outputs of other tasks; the runtime wires them together and runs
each task exactly when its inputs are ready.

### 5.3 The thread pool

`ThreadPool` (`thread.h:1166-1529`) is a process-wide singleton:

- `nthreads` worker threads (env `MAD_NUM_THREADS`) plus the main thread plus the
  separate RMI thread.
- The ready queue is a `DQueue<PoolTaskInterface*>` (`dqueue.h:79-349`), a
  cache-line-padded, dynamically growing deque with **thread-local prebuffers**
  (`MADNESS_DQ_USE_PREBUF`) to cut lock contention.
- Workers pop in batches (up to `nmax = 128`, `thread.h:1192, 1255`) and cap each
  pop at `~size/64` to avoid starving peers (`dqueue.h:290`).
- `ThreadPool::await(probe, dowork)` (`thread.h:1449-1505`) is the universal wait:
  while waiting on a condition (e.g. a future), a thread **steals and runs other
  tasks** rather than idling. Spin-vs-sleep is governed by the wait policy
  (`ENABLE_NEVER_SPIN` forces sleeping — important when oversubscribed/debugging).

### 5.4 End-to-end remote task

`taskq.add(dest, fn, args)` with `dest != me` (`world_task_queue.h:1005-1014`)
packs a `TaskHandlerInfo` — the result future's `RemoteReference`, the function
id, and serialized arguments — into an AM (`:427-439`) sent to `dest`. The handler
`remote_task_handler` (`:353-369`) deserializes and enqueues a local `TaskFn`. On
completion, setting the result future ships the value back to the caller's future
via the `set_handler` AM. **Remote task arguments must already be ready** (you
can't pass an unset remote future); local tasks may depend on unset futures.

---

## 6. `WorldContainer` — the distributed hash map

`WorldContainer<keyT,valueT>` (`worlddc.h:1125-2046`) is a shallow handle over
`WorldContainerImpl` (`:528-1098`), which **is** a `WorldObject` (so it can be a
target of remote method AMs). Local data lives in a per-process
`ConcurrentHashMap<keyT,valueT>` (`:542, 571`), thread-safe with per-entry locks
and accessor objects.

**Ownership = the pmap.** `owner(key)` (`:810-813`) delegates to a
`WorldDCPmapInterface::owner(key)` (`:135`). Implementations:

| pmap | rule | distribution |
|------|------|--------------|
| `WorldDCDefaultPmap` (`:248-267`) | `hash(key) % nproc` | balanced, scattered |
| `WorldDCLocalPmap` (`:272-289`) | always `me` | fully replicated |
| `WorldDCNodeReplicatedPmap` (`:295-315`) | lowest rank on node | one copy per node |

**Access costs:**

- `find(key)` (`:919-943`): **local** → immediate `Future` (no comm); **remote**
  → send request AM to owner; owner replies with success/failure AM that sets the
  caller's future. **One round trip (2 messages).**
- `insert`/`replace` (`:829-845`): local → hash insert under write lock; remote →
  one fire-and-forget AM to the owner.
- `send(key, memfun, args...)` (`:1459-1699`): invoke a method *on the value that
  owns `key`*, wherever it lives; the object is never serialized — only the
  method id + args travel. Auto-creates the item if absent.
- `task(key, memfun, ...)` (`:1701-1967`): like `send` but enqueues on the owner's
  task queue; supports future-dependent arguments for local targets.

Only **local** data is iterable (`begin()/end()`, `:899-917`); there is no global
iteration without communication. `redistribute()` (`:169-197`) is a 3-phase,
3-fence collective that moves misplaced keys to their new owners after a pmap
change.

---

## 7. The multiresolution function: representation

### 7.1 Handle vs implementation

`Function<T,NDIM>` (`mra.h:138-200`) stores exactly one member:
`std::shared_ptr<FunctionImpl<T,NDIM>> impl` (`mra.h:146`). Copies are shallow;
destruction defers cleanup to the next fence. `FunctionImpl<T,NDIM>`
(`funcimpl.h:945-1063`) holds the real state: `k`, `thresh`, the shared
`FunctionCommonData` (filters, quadrature) for that `k`, and the coefficient
container.

### 7.2 The adaptive tree

The domain is recursively bisected in every dimension, forming a `2^NDIM`-ary
tree. A node is named by `Key<NDIM>` (`key.h:70-175`): a level `n` plus an
integer translation vector `l` of length `NDIM`. Parent/child navigation:
`parent() = Key(n-1, l>>1)` (`key.h:289-297`); each parent has `2^NDIM` children.

`FunctionNode<T,NDIM>` (`funcimpl.h:126-495`) stores a coefficient tensor
`coeffT _coeffs` (a `GenTensor`, possibly low-rank), a `has_children` flag, and
norms. Tensor sizes:

- **reconstructed leaf:** `k^NDIM` scaling coefficients.
- **compressed internal node:** `(2k)^NDIM` (scaling + wavelet block).

The coefficient container is
`WorldContainer<Key<NDIM>, FunctionNode<T,NDIM>>` (`funcimpl.h:957, 989`).
**Therefore the function's data distribution is exactly the container's pmap over
keys.**

### 7.3 Tree states

`TreeState` (`funcdefaults.h:59-69`) tracks where coefficients live:
`reconstructed` (scaling at leaves), `compressed` (wavelets in internal nodes),
`nonstandard` / `redundant` variants used by products and operator application,
and `on_demand` (no stored coefficients — computed lazily from a functor).
Operations require specific input states, which is why `compress`/`reconstruct`
appear so often between them.

### 7.4 How the tree is spread across ranks

The pmap used for function keys (`FunctionDefaults<NDIM>::get_pmap()`,
`funcdefaults.h:401-447`) chooses the layout:

- **`SimplePmap`** (`funcimpl.h:84-101`): root on rank 0, otherwise
  `key.hash() % nproc`. Uniform but ignores parent/child locality — a parent and
  its children typically land on *different* ranks, so tree sweeps generate
  cross-rank traffic.
- **`LevelPmap`** (`funcimpl.h:103-122`): co-locates odd-level children with their
  even-level parents to improve locality of two-scale transforms.
- **`LBDeuxPmap`** (`lbdeux.h:57-91`): an explicit `map<Key,ProcessID>` produced by
  load balancing; lookups walk up to the nearest mapped ancestor, so whole
  subtrees sit on one rank.

`tree_size()` (`mraimpl.h:1882-1887`) returns the global node count via
`coeffs.size()` + `gop.sum`.

### 7.5 Load balancing (`LoadBalanceDeux`, `lbdeux.h:227-397`)

Builds a parallel cost tree (`LBNodeDeux`) using a user cost functor
`double(const Key&, const FunctionNode&)` — commonly `node.coeff().size()`
(`vmra.h`). Algorithm: (1) sum costs bottom-up; (2) descend, cutting subtrees
whose cost exceeds the per-rank target; (3) greedily bin-pack subtree roots onto
the least-loaded rank; (4) emit an `LBDeuxPmap`, then `redistribute`. This is the
knob that controls load imbalance in the performance model of §11.

---

## 8. Core operations: algorithm, parallel structure, and per-node cost

Notation used throughout:

- `k` — wavelet order (coefficients per dimension per box).
- `d ≡ NDIM` — dimensionality.
- `L` — number of occupied tree levels (depth).
- `N` — total number of tree nodes (global); `N_leaf` — leaf count.
- `P` — number of MPI ranks; `M` — separated-operator term count;
  `r` — low-rank (SVD) rank, `r ≤ k^{d/2}`.
- A **separable transform** (apply a `k×k` or `2k×2k` matrix along each of the
  `d` axes of a size-`k^d` tensor) costs `Θ(d · k^{d+1})` flops. This is the
  recurring per-node primitive.

### 8.1 `compress` / `reconstruct` — two-scale transforms

**What:** convert between scaling-at-leaves (`reconstructed`) and wavelet
(`compressed`) form by applying the two-scale filter up or down the tree.

**Parallel structure:**
- `compress` (`mraimpl.h:1500-1514`, `compress_spawn` `:3268-3328`,
  `compress_op` `:1668-1718`) is **bottom-up**: a node spawns
  `compress_spawn` tasks on the *owners of its `2^NDIM` children*
  (`:3283-3284`), collects their sum-coefficient futures, assembles a `(2k)^d`
  tensor, applies `filter()` (`:1684`), keeps the wavelet block, and returns the
  scaling block upward.
- `reconstruct` (`mraimpl.h:1468-1490`, `reconstruct_op` `:2079-2130`) is
  **top-down**: a node `unfilter()`s its `(2k)^d` block (`:2108`) and ships the
  child scaling blocks to the *owners of its children* (`:2111-2116`).

**Communication:** one task message per parent→child (or child→parent) edge whose
endpoints differ in rank — i.e. proportional to the number of tree edges crossing
rank boundaries, which the pmap controls (`SimplePmap`: almost all edges cross;
`LevelPmap`/`LBDeux`: far fewer). Each operation ends with a fence.

**Per-node cost:** one separable filter/unfilter ⇒ `Θ(d · k^{d+1})` flops; `O(k^d)`
temporary memory.

**Operation cost:** `T ≈ N · d · k^{d+1} / (P · throughput)` + (cross-rank edges) ·
(message latency) + fence.

### 8.2 `gaxpy` / `+` / `-` — linear combination

**What:** `c = α·f + β·g`, both in the same (compressed or reconstructed) state.

**Parallel structure** (`funcimpl.h:1212-1300`): iterate `g`'s local nodes; for
each, `coeffs.send(key, &FunctionNode::gaxpy_inplace, α, node, β)` to the owner of
that key in `f` (`:1275`), creating the node if absent. One fence at the end.

**Cost:** per node a tensor axpy, `Θ(r · k^d)` for low-rank (`O(k^d)` full-rank);
`O(N)` AMs (only for nodes present in `g` but mapped to a different rank than where
they're produced). Tree alignment: nodes present in only one operand are copied.

### 8.3 `multiply` — pointwise product

**What:** `h(x) = f(x)·g(x)`. Requires reconstructed inputs because products are
local in value space, not coefficient space.

**Parallel structure** (`multiply_op`, `funcimpl.h:3492-3597`, descent at
`:3654-3680`): each node, on its owner, pulls the relevant parent coefficients of
`f` and `g` (`parent_to_child`), screens whether it can be a leaf, and if so:
`coeffs2values` (transform `k^d` coeffs → values on the `npt^d` grid),
multiply pointwise, `values2coeffs` back, then rank-reduce. **No inter-node
communication** for the multiply itself; refinement may add nodes.

**Per-node cost:** two separable transforms `Θ(d · k^{d+1})` (with `npt ≈ k`) plus
`Θ(npt^d)` pointwise ⇒ dominated by `Θ(d · k^{d+1})`.

### 8.4 `inner` — inner product

**What:** scalar `⟨f|g⟩ = ∫ f g`.

**Parallel structure** (`inner_local`, `funcimpl.h:5624-5680`): iterate local
nodes, and for keys present in both, accumulate `trace_conj` of the two coefficient
tensors (`O(k^d)` each). Local accumulation via the task queue, then a single
global `gop.sum`.

**Cost:** `Θ(N_local · k^d)` compute + **one all-reduce** `O(log P)`. This is the
cheapest-communication core op — only one scalar reduction.

### 8.5 `apply` — separated convolution (the expensive one)

**What:** `g = ∫ K(x,y) f(y) dy` with `K = Σ_{μ=1}^{M}` products of 1-D kernels
(BSH, Coulomb/Poisson, etc.). Output is in non-standard form.

**Parallel structure** (`operator.h:1407-1503`; `do_apply_directed_screening`,
`funcimpl.h:5116-5220`):

- For each source node, loop over a **precomputed displacement list** ordered by
  distance (`Displacements::get_disp`, `:5156`). In 3-D the full stencil is up to
  `(2·bmax+1)^3` displacements.
- Each displacement maps the source to a target node `neighbor_in_volume(key,disp)`
  that may live on **another rank**; the per-term application `muopxv_fast2`
  applies the `M` separated 1-D operators successively across dimensions
  (`:1480-1481`), with SVD rank truncation between terms.
- **Screening** (`:5175-5184`) skips displacements whose contribution falls below
  tolerance; because displacements are distance-ordered, the loop breaks early.
- Results are accumulated at the destination via
  `coeffs.task(dest, &FunctionNode::accumulate2, ...)` (`:5783`) — a task to the
  target owner. A finalizing fence follows.

**Per-source cost:** `Θ(D_eff · M · d · k^{d+1})` where `D_eff` is the number of
displacements surviving screening (a small constant shell in practice, e.g.
~27–125 in 3-D before screening, far fewer after). Low-rank application replaces
one `k^{d+1}` factor with `M · k · r`.

**Communication:** one task per (source, surviving-displacement) whose target is
remote — the heaviest traffic of any op, carrying coefficient tensors. This is why
operator application dominates both time and network in response/SCF iterations.

### 8.6 `truncate` — norm-based pruning

**What:** drop nodes whose wavelet norm is below tolerance, promoting parent
scaling coefficients to leaves.

**Parallel structure** (`mraimpl.h:1545-1658`): pass 1 `norm_tree` computes
`sqrt(Σ child_norm²)` bottom-up; pass 2 `truncate_reconstructed_spawn/op`
descends, filters child sums, computes the wavelet-block norm, and if below
`truncate_tol(tol,key)=tol·2^{level}` deletes children and keeps the scaling block.

**Cost:** per node a filter `Θ(d·k^{d+1})` + norm `O(k^d)`; two tree sweeps, two
fences; light messages (sum coefficients only).

### 8.7 Cost summary

| Op | Per-node flops | Sweep / parallel shape | Communication |
|----|----------------|------------------------|---------------|
| compress | `d·k^{d+1}` | bottom-up | tasks on cross-rank child edges + fence |
| reconstruct | `d·k^{d+1}` | top-down | tasks on cross-rank child edges + fence |
| gaxpy | `r·k^d` | per-node local | `O(N)` AMs + 1 fence |
| multiply | `d·k^{d+1}` | per-node local | none + fence |
| inner | `k^d` | local + reduce | **1 all-reduce** `O(log P)` |
| apply | `D_eff·M·d·k^{d+1}` | source→target stencil | `O(N·D_eff)` tasks + fence |
| truncate | `d·k^{d+1}` | 2 sweeps | child-sum tasks + 2 fences |

---

## 9. Putting it together: a performance-model framework

Model any operation as the sum of three terms, taken over the **critical path**
(the rank/level that finishes last, not the average):

```
T_op  =  T_compute  +  T_comm  +  T_sync
```

**Compute term.** Per-node flop count `c_node` (from §8.7) times nodes processed
on the busiest rank, divided by per-core flop rate `R` times useful cores
`(MAD_NUM_THREADS − overhead)`:

```
T_compute ≈ (N_busiest · c_node) / (R · n_threads_eff)
```

with `c_node` scaling as `k^{d+1}` for transform-type ops. Because `k` jumps from
6 to 10 at the final protocol and node count `N` grows with it, the final protocol
dominates: for `k=10` vs `k=6`, the per-node transform cost rises by
`(10/6)^{d+1}` and (per `CLAUDE.md`) `N_leaf` rises ~4.6×.

**Communication term.** `T_comm = n_msg · α + V_bytes · β`, where `α` is per-message
latency (MPI transit + receiver poll/backoff up to `MAD_BACKOFF_US`, + sender
buffer wait under congestion when in-flight > `nsend`), and `β` is inverse
bandwidth (bounded by the single mutex-serialized RMI thread and `MAD_BUFFER_SIZE`
chunking). Message counts come from §8.7 and the pmap: cross-rank tree edges for
sweeps, `N·D_eff` for apply, one all-reduce for inner.

**Sync term.** `T_sync ≈ (#fences) · (fence_cost)`, fence_cost `≈ 2 · log₂(P) · α`
in the clean case, more if stragglers re-trigger passes. Count fences from the
operation sequence — every `compress`/`reconstruct`/`truncate`/`apply` boundary is
one.

**Load-balance factor.** Multiply `T_compute` by `N_max_rank · P / N` (the
imbalance ratio). `SimplePmap` randomizes nodes (good balance, bad locality);
`LBDeux` balances by cost functor (good locality, balance only as good as the cost
estimate). This factor is often the difference between modeled and observed time.

---

## 10. Worked example: one response solver iteration

A frequency-response iteration applies a BSH operator to each of `n_occ`
perturbed orbitals per spatial direction, with surrounding arithmetic. Sketch:

```
per orbital component:
   multiply (potential × orbital)   →  Θ(N · d·k^{d+1})  compute, fence
   compress / reconstruct pair      →  2 sweeps,          2 fences
   apply (BSH)                      →  Θ(N · D_eff·M·d·k^{d+1}) compute,
                                       O(N·D_eff) remote tasks, fence
   gaxpy (accumulate, orthogonalize)→  Θ(N · r·k^d),       fences
   inner (norms, overlaps)          →  Θ(N · k^d) + all-reduce
```

The `apply` term dominates compute; the fences (≈4–6 per component) dominate
small-`P` latency; the cross-rank coefficient traffic dominates network at large
`P`. Per-iteration time then scales roughly as:

```
T_iter ≈ n_occ · n_dir · [ apply_compute / (P·R) · imbalance
                           + (apply_msgs · α + apply_bytes · β)
                           + n_fences · 2 log₂P · α ]
```

To validate a model against a run, instrument the four quantities directly
(next section) rather than trusting the constants.

---

## 11. Memory model and the ground-state replication bottleneck

This connects the runtime to the OOM behavior in `CLAUDE.md`.

**Per-node bytes.** A reconstructed leaf stores `k^d` doubles = `8·k^d` bytes
(plus `GenTensor`/node overhead). For `d=3`: `8·k³` bytes per node.

**Per-function bytes (one rank).** With `N_leaf` leaves distributed over `P`
ranks and imbalance factor `φ ≥ 1`:

```
mem_per_function_per_rank ≈ φ · (N_leaf / P) · 8 · k³
```

**Ground-state replication.** In `state_parallel "on"` (subworld) mode, every MPI
task holds a **complete replica of all `n_occ` ground-state orbitals** at the
current protocol — they are not distributed across the subworld, they are
replicated. So the dominant resident term is:

```
mem_ground_per_task ≈ n_occ · N_leaf(k) · 8 · k³   (+ response overhead)
```

This is exactly the `n_occupied × n_leaves × k³ × 8 bytes` figure in `CLAUDE.md`.
It is independent of `P` (replication, not distribution), which is why adding
ranks does **not** relieve the per-task OOM, and why the documented workaround is
to halve `--ntasks-per-node` (fewer replicas per node) rather than add nodes.

**Why `k=10` is the cliff.** Going `k=6 → k=10` raises `8·k³` by `(10/6)³≈4.6×`
*and* `N_leaf` by ~4.6× (empirical), so the ground-state replica grows by roughly
an order of magnitude right at the final protocol — consistent with C6H6 (21 occ)
and naphthalene (34 occ) dying with exit 137 immediately after `PROTOCOL_POLICY`.

**Levers the architecture already offers** (mapping to `CLAUDE.md` priorities):
- A **pre-flight estimate** = `n_occ · N_leaf · 8 · k³` per task vs
  `/proc/meminfo`, printed before subworld allocation — cheap, prevents silent
  SIGKILL.
- **Node-local shared replica** via MPI shared-memory windows or
  `WorldDCNodeReplicatedPmap` (`worlddc.h:295-315`) would cut the multiplier from
  `ntasks_per_node` to 1 for the ground-state term.
- **Lazy/on-demand orbital loading** maps onto the `on_demand` tree state plus
  `ResponseIO` archives — keep orbitals evictable rather than resident.

---

## 12. How to measure (instrument the model)

- **Per-task RSS:** `worldmem.h` exposes RSS queries; emit
  `MEMORY_HWM rank=N protocol=N rss_GB=X` at protocol boundaries in `FrequencyLoop`
  (a `CLAUDE.md` priority) instead of post-mortem `sacct`.
- **Message/byte counts:** `RMI::stats` tracks `nmsg_recv`, `nbyte_recv`
  (`worldrmi.cc:96-98`), and send-side counts; `World::print_stats`
  (`world.cc:236-421`) dumps them. These give `n_msg`, `V_bytes` for `T_comm`.
- **Task counts / queue depth:** `DQStats` (`dqueue.h:59-68`, enable
  `MADNESS_DQ_STATS`) and `taskq.size()`.
- **Tree size / node distribution:** `FunctionImpl::tree_size()`
  (`mraimpl.h:1882`) and per-rank `coeffs.size()` give `N`, `N_busiest`, hence the
  imbalance factor `φ`.
- **Fence count:** instrument `gop.fence()` call sites along an operation sequence;
  combine with `P` for `T_sync`.
- **Tunables that enter the model:** `MAD_NUM_THREADS`, `MAD_RECV_BUFFERS`,
  `MAD_BUFFER_SIZE`, `MAD_SEND_BUFFERS`, `MAD_BACKOFF_US`,
  `MAD_MAX_REDUCEBCAST_MSG_SIZE`, `ENABLE_NEVER_SPIN`, and the chosen pmap /
  load-balance cost functor.

A practical calibration loop: run H2O and CH3OH (known-good), fit `R`, `α`, `β`,
and `φ` from the measured counts above, then *predict* C6H6 and naphthalene
memory and time before launching — and check the prediction against the
`mul_sparse` table in `CLAUDE.md`.

---

## 13. File index

Runtime (`src/madness/world/`): `world.h/.cc`, `safempi.h`, `worldmpi.h`,
`worldam.h/.cc`, `worldrmi.h/.cc`, `worldgop.h/.cc`, `world_task_queue.h`,
`taskfn.h`, `future.h/.cc`, `dependency_interface.h`, `worldref.h`, `thread.h`,
`dqueue.h`, `worlddc.h`, `worldhashmap.h`, `worldmem.h`.

MRA (`src/madness/mra/`): `mra.h`, `funcimpl.h`, `mraimpl.h`, `key.h`,
`funcdefaults.h`, `function_common_data.h`, `operator.h`, `vmra.h`, `lbdeux.h`,
`twoscale.h`, `macrotaskq.h`.
