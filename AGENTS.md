# Agent Instructions for MADNESS

Canonical guidance for coding agents working on this repo. For installing
dependencies and walking through a first build, see `INSTALL.md`. This file
covers what an agent needs beyond that: the source layout, how tests are
organized, and the project-specific gotchas that are easy to get wrong.

`CLAUDE.md` is a symlink to this file; there is no vendor-specific variant.

**How to use this file.** Treat it as a starting map, not authority. The code is the source of truth — when this file disagrees with what `grep`/`Read` show you, trust the code and flag the discrepancy so this file can be updated.

_Last reviewed against master:_ 2026-04-23.

## Platforms

MADNESS targets POSIX only — Linux and macOS on x86-64, ARM64, and IBM
Power. Windows is not supported. Do not add `#ifdef _WIN32` branches, path
handling that assumes backslashes, or dependencies that are Windows-only.
POSIX facilities (`pthread`, `mmap`, `fork`, `dlopen`, `sys/stat.h`, …) may
be used directly without feature detection beyond the existing
`HAVE_*` macros.

## Source tree

- `src/madness/world` — parallel runtime, a.k.a. **MADWorld** (tasks, futures, MPI wrappers).
- `src/madness/tensor` — dense tensors and BLAS/LAPACK wrappers.
- `src/madness/mra` — multiresolution analysis (`Function<T,NDIM>`, operators).
- `src/madness/chem` — electronic-structure building blocks (SCF, XC, …).
- `src/madness/misc`, `src/madness/external` — utilities, bundled deps.
- `src/apps/*` — end-user applications (`moldft`, `nemo`, `cc2`, `mp2`, …).
- `src/examples`, `src/pymadness` — samples and Python bindings.

## Build directories

Prefer reusing an existing out-of-tree build directory over creating a new
one from scratch, so you don't re-pay CMake configure time. Common local
conventions are `cmake-build-{debug,release,relwithdebinfo}` and `build*`.
If several exist side-by-side, check their
cached CMake options (`grep … CMakeCache.txt`) — configure state often
diverges between them.

For focused iteration, these scope flags drop rebuild time substantially:

- `-DMADNESS_BUILD_MADWORLD_ONLY=ON` — build only the runtime, skipping
  `tensor`/`mra`/`chem`/`apps`. Right for work confined to `src/madness/world`.
- `-DMADNESS_BUILD_LIBRARIES_ONLY=ON` — build the libraries, skip `src/apps`.
  Right for changes to numerical layers when you don't need `moldft`/`nemo`/
  `cc2`.

## Tests

Tests live next to the sources they exercise and are registered through
`add_unittests(component sources libs labels)` (see
`cmake/modules/AddUnittests.cmake`). The macro emits:

- one CTest entry per source file at `madness/test/<component>/<name>/run`,
  with `<name>` derived from the source filename, and
- a per-component `madness/test/<component>/build` entry that builds the
  unit-test targets on demand.

Both entries carry the `labels` passed to the macro. There is no filename
convention — most test sources follow `test*.cc` / `test_*.cc` by habit,
but the macro accepts any source name (e.g. `src/madness/misc/interp3.cc`,
`src/examples/periodic/erfcr.cc`).

The `check-short-madness` target runs everything labeled `short` or `medium`
via `ctest -L "short|medium"`. To run a subset without the full suite:

```
ctest -L short -R "madness/test/mra"         # all mra tests in the short set
ctest -R "madness/test/mra/test_cloud/run"   # a single test by name
```

Each test is also a normal binary under
`<build-dir>/<path-of-source-dir>/<name>` (mirroring the source tree, so
`src/madness/mra` tests land at `<build-dir>/src/madness/mra/<name>`, and
`src/examples/periodic` tests at `<build-dir>/src/examples/periodic/<name>`)
and can be invoked directly for debugging.

### Smoke test

`testsuite` in `src/madness/mra` is the broad numerical-library regression
harness. From a configured build directory:

```
ninja testsuite
./src/madness/mra/testsuite
```

Success: the program prints `testsuite passed:  true` and exits 0. For the
MPI variant:

```
MAD_NUM_THREADS=2 mpiexec -np 2 ./src/madness/mra/testsuite
```

### moldft smoke

```
ninja moldft
./src/apps/moldft/moldft --geometry=water --dft="xc=lda"
```

Exit code 0 on success. Reference inputs for regression checks are in
`src/apps/moldft/tests/*.in`.

## Project-specific gotchas

- **BLAS must be sequential.** MADNESS owns parallelism via its own task
  pool; a threaded BLAS inside a task oversubscribes cores. Use
  single-threaded MKL, Accelerate, or OpenBLAS-serial.
- **Static libs are the default because ASLR is assumed.** Building shared
  libs on a system without ASLR requires
  `-DMADNESS_ASSUMES_ASLR_DISABLED=ON`. Do not flip `BUILD_SHARED_LIBS`
  casually.
- **MPI is required by default.** `-DENABLE_MPI=OFF` is rarely used but
  still supported and maintained — in that mode `stubmpi.h` stands in for
  the subset of MPI that `SafeMPI` uses. Any new code that reaches for MPI
  must go through `SafeMPI` / `World::mpi` (see runtime invariants below)
  so this path continues to compile.
- **`-DENABLE_NEVER_SPIN=ON`** is a workaround for VMs/containers where
  spinlocks burn CPU. Leave it off on bare metal.
- **`-DENABLE_GENTENSOR=ON`** swaps in a compressed-tensor code path used
  for 6D work (MP2, CC2). Off by default and tested in only one CI cell,
  so changes in the tensor layer should be checked with it on when
  relevant.

## Runtime & deployment

These are concerns for launch scripts and performance tuning, not the source
itself — and they are largely undiscoverable by reading the code.

### Task backend

`-DMADNESS_TASK_BACKEND={Pthreads,TBB,PaRSEC}` (CMake, default `Pthreads`):

- **Pthreads** — MADNESS's own pool. Pick this when you need execution on
  just the main thread; TBB and PaRSEC always spawn at least one worker.
- **PaRSEC** — preferred for production. Disable PaRSEC's own thread
  binding or it will fight the launcher:
  `export PARSEC_MCA_bind_threads=0`.
- **TBB** — Intel TBB task pool; autodetected when available.

### Thread count and the comm thread

Whenever `nranks > 1`, each rank spawns a dedicated MADNESS communication
thread **in addition to** the `MAD_NUM_THREADS` workers. MPI implementations
typically want a core for their own progress thread too, so the budget is:

- Comm-heavy workloads (MRA integral operators, TA tensor contractions):
  `MAD_NUM_THREADS = hwthreads_per_rank − 2` (one hwthread each for the
  MADNESS and MPI comm threads).
- Small ranks (≈16 hwthreads) or light comm: `... − 1` is usually fine —
  the two comm threads co-exist on one core without meaningful contention.
- Never oversubscribe.

### Rank pinning

For MRA workloads, pin each rank to an L3 slice (or NUMA domain on systems
without L3-sized domains). Launch flags are deployment-specific:

- OpenMPI: `mpirun --map-by l3cache` or `--map-by numa`.
- Other MPI vendors: flags vary, some lack pinning entirely.
- SLURM: `--cpu-bind` / `--distribution` rather than the MPI flags.

Document the chosen command with the application, not here.

### Communication buffers

One-sided sends/receives use fixed-size buffers; 2-sided is supported but
much slower, so the 1-sided path should fit the common message.

- `MAD_BUFFER_SIZE` (env var, bytes) — default `3*512*1024` ≈ 1.5 MiB,
  defined at `src/madness/world/worldrmi.h`.
- MRA rule of thumb: a tree-node payload is
  `sizeof(T) · (2k)^NDIM + metadata`. The default comfortably holds
  several nodes at typical `k ≤ 10, NDIM ≤ 3`; raise it for larger `k`,
  6D work, or TiledArray-style bulk tensor transfers.
- MPI vendors have their own eager/rendezvous cutoffs (MVAPICH in
  particular) — consult vendor docs when tuning comm perf.

### Apple Silicon + Homebrew OpenMPI

Homebrew's OpenMPI on Apple M-series links HWLOC, which probes hardware
via OpenCL and crashes at startup if launched outside `mpiexec`. Either:

- Always launch via `mpiexec`, even for `-np 1`, or
- `export HWLOC_COMPONENTS=darwin,no_os,-opencl` to skip the OpenCL probe.

## Coding conventions

### Error handling

MADNESS ships several macros with distinct intent (`madness_exception.h`).
Pick by *why* you're checking, not by what's short to type:

- **`MADNESS_CHECK(cond)` / `MADNESS_CHECK_THROW(cond, msg)`** — always-on
  checks: state-machine preconditions, collective-ordering requirements,
  user input validation. Use these for invariants that must hold in release
  builds.
- **`MADNESS_ASSERT(cond)`** — debug-only, expensive invariants. Its runtime
  behavior (throw / abort / std::assert / disabled) depends on the
  `ASSERTION_TYPE` CMake variable; release builds often compile these out.
  Do **not** use for checks that must run in production.
- **`MADNESS_EXCEPTION(msg, value)`** — explicit throw sites (unreachable
  branches, unsupported configurations). Not a conditional.
- **`MADNESS_ASSERT_NOEXCEPT`** — use inside `noexcept` functions where a
  throwing `MADNESS_ASSERT` would call `std::terminate`.

### CMake harness

MADNESS is designed to be consumed by parent projects both as an installed
package *and* as a subproject (via `FetchContent` / `add_subdirectory`, e.g.
from TiledArray and MPQC). Maintain both invariants when editing CMake.

- **No generic target names.** `check`, `doc`, `unittests`, `applications`
  would collide with a parent project's targets. Define subproject-scoped
  targets with `add_custom_target_subproject(madness <name> …)` from
  `cmake/modules/AddCustomTargetSubproject.cmake`; it creates `<name>-madness`
  and hooks it into `<name>` if the parent already defined it.
- **Use the `add_mad_*` wrappers**, not raw `add_library` / `add_executable`.
  `add_mad_library` / `add_mad_hdr_library` / `add_mad_executable` apply
  project-wide include paths, install rules, and the
  `$<BUILD_INTERFACE:…>` / `$<INSTALL_INTERFACE:…>` generator expressions
  that make targets consumable from both build and install trees.
- **Consumable from both install tree and build tree.** Preserved via
  `install(EXPORT madness …)` plus
  `export(EXPORT madness FILE ${PROJECT_BINARY_DIR}/madness-targets.cmake)`
  in the top-level `CMakeLists.txt`, and the genex pair in `add_mad_library`.
  Breaks in either path (missing install rules on new targets, hard-coded
  paths outside the genex) are silent regressions for downstream consumers.
- **Optional dependencies follow the `FindOrFetch*` pattern** (e.g.
  `FindOrFetchCereal.cmake`, `FindOrFetchPARSEC.cmake`): try `find_package`
  first, fall back to `FetchContent`. Harness mantra: use preexisting deps
  when available, build from source otherwise — except for the most basic
  (MPI, BLAS), which must be found on the system.
- **Component libraries** (`MADworld`, `MADmra`, `MADtensor`, `MADchem`, …)
  are the public link targets. The aggregate `madness` library is composed
  from per-component object libraries (`MADworld-obj`, etc.), which (a) lets
  downstream projects pull a lean subset via `MADNESS_BUILD_MADWORLD_ONLY` /
  `MADNESS_BUILD_LIBRARIES_ONLY`, and (b) lets developers build a single
  component quickly without its prerequisites — useful for shorter build
  tests during iteration.
- **`CMAKE_FIND_NO_INSTALL_PREFIX=ON`** is set as a safety net so
  `find_package` doesn't pick up a stale installed copy during development.

### Commit messages

Write plain commit messages describing the change. **Do not append
`Co-Authored-By: …` trailers crediting AI tooling** (Claude, Copilot,
Codex, etc.). This project does not credit tooling in commit metadata —
the same way it doesn't credit the editor, compiler, or the vast body
of prior work the tooling was trained on.

## Consumer-facing invariants

MADNESS is consumed by external projects — TiledArray leans on the parallel
runtime (`src/madness/world`), MPQC leans on MRA (`src/madness/mra`). When
editing in either area, preserve these invariants; they are not enforced at
compile time and breakage surfaces far downstream.

### Runtime (`src/madness/world`)

- `madness::initialize` must bracket all `World` / MPI use; `madness::finalize`
  mirrors it. Creating a `World` or submitting tasks outside that window is
  undefined.
- `World` lifetimes must nest — destroy in reverse creation order.
- `WorldObject` and `WorldContainer` require collective construction: every
  rank in the enclosing `World` must create the same object in the same order.
- Active-message handlers (`World::am`) run on runtime threads and must not
  block waiting on further tasks or collective ops; post continuations instead.
- **MPI is initialized at `MPI_THREAD_MULTIPLE`** — the only level tested
  today, and required because tasks can spawn further (possibly remote)
  tasks. This also constrains MPI-implementation choice: some MPI stacks
  have performance cliffs under `THREAD_MULTIPLE`.
- **Use `SafeMPI` / `World::mpi` for every MPI call; never raw `MPI_*`.**
  Three reasons: automated error checking; optional serialization fallback
  when the MPI implementation is not reliably thread-safe; and the
  `ENABLE_MPI=OFF` build replaces MPI with `stubmpi.h`, which only declares
  the subset of MPI that `SafeMPI` uses — raw MPI calls won't compile there.
- **Never call MPI blocking collectives.** They park the calling thread
  outside the task scheduler and can deadlock progress. Use MADNESS global
  ops (`World::gop`) — non-blocking, composable with the runtime, and they
  serialize user types via `madness::archive`.
- `madness::archive` is the public serialization extension point — downstream
  code teaches `Future`, `WorldContainer`, and AM payloads to transport custom
  types by specializing archive traits, not by reaching into internals.
- Outside tasks and the main thread, assume nothing in the runtime is
  thread-safe unless the header documents it.
- **Prefer explicit `World&` parameters** over `World::get_default()` in any
  new API. TiledArray builds with `-DMADNESS_DISABLE_WORLD_GET_DEFAULT=ON`
  to surface code that implicitly assumes the default World; MADNESS itself
  is laxer, but APIs written to accept an explicit `World&` work in both
  settings. App-level code may use the default where it is clearly the top
  level.

### MRA (`src/madness/mra`)

- `madness::startup(world, argc, argv, quiet)` is the **MRA** initializer —
  loads twoscale coefficients, autocorrelation data, and populates
  `FunctionDefaults`. It must be called after `madness::initialize` and
  before constructing any `Function<T,NDIM>` or MRA operator. Runtime-only
  (MADWorld) code does not need it.
- `Function<T,NDIM>` is a **shallow-copy handle** to a distributed
  implementation — Python-like reference semantics, *not* copy-on-write.
  `Function b = a;` shares the impl; mutations through `b` are observed
  through `a`.
- `Function` operations enqueue **async** tasks. Results are not observable
  until an implicit fence (an observing op such as `norm2`, `trace`, `inner`)
  or an explicit `world.gop.fence()`.
- `Function` carries a `TreeState` (`reconstructed`, `compressed`,
  `nonstandard`, `redundant`, …). Many operations assert on state — use
  `change_tree_state` / `reconstruct` / `compress` rather than assuming.
- Prefer the vectorized forms in `src/madness/mra/vmra.h`
  (`std::vector<Function>` arithmetic, `apply`, inner products, …) over
  hand-written loops over individual `Function`s — this is the canonical
  idiom and the overhead is negligible.
- `FunctionDefaults<NDIM>` is process-global mutable state (`k`, `thresh`,
  cell, boundary conditions). Set it before constructing any `Function<NDIM>`;
  changing it mid-run will not retroactively update existing functions.
- `SeparatedConvolution` operators are tied to a specific `k`, box size, and
  threshold — reuse one across `Function`s built with different defaults and
  results are undefined.
- Drive mutation from tasks or the main thread; don't mutate `Function`s
  from arbitrary user threads.
- **Checkpointing.** `Function::save` / `Function::load` go through
  `ParallelOutputArchive` / `ParallelInputArchive` (binary fstream under the
  hood). `FunctionDefaults<NDIM>` must match between save and load — same
  cell, `k`, `thresh`, and boundary conditions — or `load` throws. There
  is no migration path.

For API depth, see the Doxygen groups rather than restating signatures here.

## Lint

A `.clang-tidy` configuration is checked in. Run `clang-tidy` against
changed translation units before finishing non-trivial C++ edits.

## Debugging aids

- Memory profiling of `FunctionImpl` objects is always available via
  `MemoryMeasurer::measure_and_print(world)`.
- `-DENABLE_MEM_STATS=ON` and `-DENABLE_TENSOR_INSTANCE_COUNT=ON` provide
  coarser, always-on counters for leak hunting (see `INSTALL.md`).

## Reviewing a PR

Guidance for agents doing PR review — focused on regressions CI won't
catch because its matrix is narrow. Flag these; don't nitpick style
(`.clang-tidy` is the source of truth).

- **Raw MPI calls.** Anything added outside `SafeMPI` / `World::mpi`
  breaks `-DENABLE_MPI=OFF` and bypasses error checking / thread-safety
  serialization.
- **Blocking MPI collectives** (`MPI_Barrier`, `MPI_Allreduce`, …).
  Push toward `World::gop`.
- **Implicit default World.** `World::get_default()` in library code
  breaks downstream builds (notably TiledArray's
  `-DMADNESS_DISABLE_WORLD_GET_DEFAULT=ON`). Push for an explicit
  `World&` parameter.
- **Non-collective construction** of `WorldObject` / `WorldContainer`,
  or any code path where ranks can diverge on construction order.
- **Threaded-BLAS assumptions** (parallel loops around `gemm`, relying
  on BLAS threading for speed).
- **`Function` mutation from user threads** (not tasks or the main
  thread) or mis-handled `TreeState` preconditions.
- **Checkpoint-format breaks.** Any change to what `Function::save`
  writes, or to `FunctionDefaults` fields consumed by `load`, is a
  silent compat break — call it out.
- **Tensor-layer changes without GENTENSOR coverage.** CI only tests
  one GENTENSOR cell; ask the author to verify with
  `-DENABLE_GENTENSOR=ON` locally if the change could affect the
  compressed-tensor path.

## Documentation

Developer API docs are generated by Doxygen and published at
<https://m-a-d-n-e-s-s.github.io/madness/api-doc/>. User-facing docs live at
<https://madness.readthedocs.io/en/latest/>.
