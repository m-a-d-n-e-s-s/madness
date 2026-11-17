/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680
*/

#ifndef MADNESS_MRA_WINDOWED_NS_HALO_H__INCLUDED
#define MADNESS_MRA_WINDOWED_NS_HALO_H__INCLUDED

/// \file windowed_ns_halo.h
/// \brief Fine-driven single-pass virtual boundary halo cache (`virt_cache`, Stage 1b)
///        and single-pass bottom-up virtual restriction / dynamic refinement assembly
///        (Stage 5a) for Windowed NS 3D convolutions (serial & multi-rank MPI).
///
/// Enforces Invariant 1 (Zero Exterior Polynomial Extrapolation): every virtual node
/// `V` is strictly a sub-box of an existing coarser real leaf `B in L` (`level(B) < level(V)`),
/// prolonged top-down one dyadic step at a time via exact interior Legendre filters
/// `(h0, h1)` (`general_fast_transform`) and restricted bottom-up via exact adjoint
/// scaling filters `(h0T, h1T)` (`general_fast_transform` + `aligned_axpy`).

#include <madness/mra/windowed_ns_stencil.h>
#include <madness/mra/mra.h>
#include <madness/world/range.h>
#include <madness/world/thread.h>
#include <madness/world/world_object.h>
#include <madness/world/world_task_queue.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <exception>
#include <functional>
#include <mutex>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#ifdef HAVE_INTEL_TBB
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/partitioner.h>
#endif

namespace madness {

/// Thread-local flag set to `true` when executing inside a function-level batched
/// task (e.g., vector `apply`), preventing nested outer batch launches.
inline thread_local bool g_in_windowed_batch_task = false;

/// RAII scope setting `g_in_windowed_batch_task = true` on a batch worker thread.
struct WindowedBatchTaskScope {
    bool prev;

    WindowedBatchTaskScope() noexcept
        : prev(g_in_windowed_batch_task) {
        g_in_windowed_batch_task = true;
    }

    ~WindowedBatchTaskScope() {
        g_in_windowed_batch_task = prev;
    }

    WindowedBatchTaskScope(const WindowedBatchTaskScope&) = delete;
    WindowedBatchTaskScope& operator=(const WindowedBatchTaskScope&) = delete;
};

namespace windowed_parallel_detail {

inline constexpr std::size_t kMaxParallelChunks = 128;
inline constexpr std::size_t kDefaultParallelGrain = 16;
inline constexpr std::size_t kHeavyParallelGrain = 4;
inline constexpr std::size_t kOverdecompPerThread = 4;
inline thread_local std::uint32_t g_windowed_task_depth = 0;
inline std::atomic<std::size_t> g_unfinished_outer_tasks{0};

struct WindowedTaskDepthScope {
    WindowedTaskDepthScope() noexcept { ++g_windowed_task_depth; }
    ~WindowedTaskDepthScope() { --g_windowed_task_depth; }
    WindowedTaskDepthScope(const WindowedTaskDepthScope&) = delete;
    WindowedTaskDepthScope& operator=(const WindowedTaskDepthScope&) = delete;
};

/// Compact worker-thread slot registry decoupling per-thread scratch storage (`O(N_threads)`)
/// from over-decomposed `madness::Range` task counts (`O(8 * N_threads)`), with full reentrancy
/// safety whenever a thread nested-executes tasks during `Future::get` / `ThreadPool::await`.
class WindowedThreadSlotRegistry {
public:
    static std::size_t acquire_slot() {
        auto& st = state();
        std::lock_guard<std::mutex> lock(st.mtx);
        for (std::size_t w = 0; w < 2u; ++w) {
            const std::uint64_t free_bits = ~st.used_mask[w];
            if (free_bits != 0ULL) {
                const std::size_t bit = static_cast<std::size_t>(__builtin_ctzll(free_bits));
                const std::size_t slot = w * 64u + bit;
                if (slot < kMaxParallelChunks) {
                    st.used_mask[w] |= (1ULL << bit);
                    std::size_t cur_hi = st.high_water.load(std::memory_order_relaxed);
                    while (slot + 1u > cur_hi &&
                           !st.high_water.compare_exchange_weak(
                               cur_hi, slot + 1u, std::memory_order_release, std::memory_order_relaxed)) {
                    }
                    return slot;
                }
            }
        }
        MADNESS_EXCEPTION("WindowedThreadSlotRegistry: exceeded kMaxParallelChunks (128) active slots", 1);
        return 0u;
    }

    static void release_slot(std::size_t slot) noexcept {
        if (slot >= kMaxParallelChunks) return;
        auto& st = state();
        std::lock_guard<std::mutex> lock(st.mtx);
        st.used_mask[slot / 64u] &= ~(1ULL << (slot % 64u));
    }

    static std::size_t slot_bound() noexcept {
        return std::min<std::size_t>(
            kMaxParallelChunks, state().high_water.load(std::memory_order_acquire));
    }

private:
    struct State {
        std::mutex mtx;
        std::uint64_t used_mask[2] = {0ULL, 0ULL};
        std::atomic<std::size_t> high_water{1u};
    };

    static State& state() {
        static State st;
        return st;
    }
};

struct TlsPrimarySlot {
    std::size_t primary_slot = SIZE_MAX;
    std::uint32_t active_depth = 0;

    std::size_t get_or_init() {
        if (primary_slot == SIZE_MAX) {
            primary_slot = WindowedThreadSlotRegistry::acquire_slot();
        }
        return primary_slot;
    }

    ~TlsPrimarySlot() {
        if (primary_slot != SIZE_MAX) {
            WindowedThreadSlotRegistry::release_slot(primary_slot);
            primary_slot = SIZE_MAX;
        }
    }
};

inline TlsPrimarySlot& this_thread_slot_state() {
    static thread_local TlsPrimarySlot tls;
    return tls;
}

class ScopedWorkerSlot {
public:
    ScopedWorkerSlot() {
        auto& tls = this_thread_slot_state();
        if (tls.active_depth == 0u) {
            slot_ = tls.get_or_init();
            is_overflow_ = false;
        } else {
            slot_ = WindowedThreadSlotRegistry::acquire_slot();
            is_overflow_ = true;
        }
        ++tls.active_depth;
    }

    ~ScopedWorkerSlot() {
        auto& tls = this_thread_slot_state();
        --tls.active_depth;
        if (is_overflow_) {
            WindowedThreadSlotRegistry::release_slot(slot_);
        }
    }

    ScopedWorkerSlot(const ScopedWorkerSlot&) = delete;
    ScopedWorkerSlot& operator=(const ScopedWorkerSlot&) = delete;

    std::size_t slot() const noexcept { return slot_; }

private:
    std::size_t slot_ = 0;
    bool is_overflow_ = false;
};

inline std::size_t slot_bound() noexcept {
    return WindowedThreadSlotRegistry::slot_bound();
}

template <typename Fn>
struct RangeBlockOp {
    struct SharedCtx {
        const Fn* fn = nullptr;
        std::size_t n_items = 0;
        std::size_t grain = 1;
        std::atomic<bool> has_exc{false};
        std::mutex exc_mtx;
        std::exception_ptr captured_exc;
    };

    SharedCtx* ctx = nullptr;

    bool operator()(std::size_t blk) const {
#if !defined(HAVE_INTEL_TBB) && !defined(HAVE_PARSEC)
        // Flush binary-split `ForEachTask` children on coarse split leaders (`blk % 4 == 0`)
        // so idle pool threads wake up immediately at `blk == 0` while cutting `DQueue` mutex
        // lock frequency 4x across leaf blocks.
        if ((blk & 3u) == 0u) {
            thread_purge();
        }
#endif
        if (ctx->has_exc.load(std::memory_order_relaxed)) return true;
        const std::size_t b = blk * ctx->grain;
        const std::size_t e = std::min<std::size_t>(ctx->n_items, b + ctx->grain);
        WindowedTaskDepthScope depth_scope;
        ScopedWorkerSlot slot_guard;
        try {
            (*ctx->fn)(b, e, slot_guard.slot());
        } catch (...) {
            std::lock_guard<std::mutex> lock(ctx->exc_mtx);
            if (!ctx->captured_exc) {
                ctx->captured_exc = std::current_exception();
            }
            ctx->has_exc.store(true, std::memory_order_release);
        }
        return true;
    }
};

/// Over-decomposed MADNESS native task executor over `[0, n_items)`:
///   - Splits `[0, n_items)` into `~4 * N_threads` fine-grained `madness::Range` blocks (`>= min_grain`)
///     so heterogeneous big.LITTLE (P/E) cores and dynamically clocked cores naturally self-balance.
///   - Routes execution through `world.taskq.for_each(madness::Range)` + `thread_purge()` + `Future::get()`
///     on Pthreads/PaRSEC and `ThreadPool::tbb_arena->execute` + `tbb::parallel_for(madness::Range)` on TBB.
///   - Assigns each executing worker thread a reentrancy-safe `tid in [0, slot_bound())` (`O(N_threads)` scratch).
template <typename Fn>
inline void parallel_for_chunks(World& world,
                                std::size_t n_items,
                                std::size_t min_grain,
                                const Fn& fn) {
    if (n_items == 0) return;
    const std::size_t eff_grain = std::max<std::size_t>(min_grain, 1u);
    const std::size_t pool_sz = ThreadPool::size();
    const std::size_t total_threads = pool_sz + 1u;
    if (pool_sz == 0u || n_items < 2u * eff_grain) {
        ScopedWorkerSlot slot_guard;
        fn(0u, n_items, slot_guard.slot());
        return;
    }

    // When executing inside an outer `windowed_for_each_outer_tasks` batch (e.g., an `HF Exchange`
    // pair tile or `ApplyTask` orbital batch) while `>= total_threads` outer tasks are still
    // unfinished, every worker thread already has a full outer task in L1/L2 cache; run inner
    // sweeps inline with zero `ForEachTask` allocations or `DQueue` mutex contention, and
    // automatically enable inner `madness::Range` over-decomposition as soon as the batch enters
    // its tail (`g_unfinished_outer_tasks < total_threads`) so fast P-cores help finish stragglers.
    if (g_windowed_task_depth > 0u &&
        g_unfinished_outer_tasks.load(std::memory_order_relaxed) >= total_threads) {
        ScopedWorkerSlot slot_guard;
        fn(0u, n_items, slot_guard.slot());
        return;
    }

    const std::size_t target_blocks = kOverdecompPerThread * total_threads;
    const std::size_t grain = std::max<std::size_t>(
        eff_grain, (n_items + target_blocks - 1u) / target_blocks);
    if (n_items <= grain) {
        ScopedWorkerSlot slot_guard;
        fn(0u, n_items, slot_guard.slot());
        return;
    }

#ifdef HAVE_INTEL_TBB
    (void)world;
    if (ThreadPool::tbb_arena) {
        std::atomic<bool> has_exc{false};
        std::mutex exc_mtx;
        std::exception_ptr captured_exc;
        const int tbb_grain =
            static_cast<int>(std::min<std::size_t>(grain, static_cast<std::size_t>(0x3fffffff)));
        const madness::Range<std::size_t> item_range(0u, n_items, tbb_grain);
        ThreadPool::tbb_arena->execute([&]() {
            tbb::parallel_for(
                item_range,
                [&](const madness::Range<std::size_t>& r) {
                    if (has_exc.load(std::memory_order_relaxed)) return;
                    WindowedTaskDepthScope depth_scope;
                    ScopedWorkerSlot slot_guard;
                    try {
                        fn(r.begin(), r.end(), slot_guard.slot());
                    } catch (...) {
                        std::lock_guard<std::mutex> lock(exc_mtx);
                        if (!captured_exc) {
                            captured_exc = std::current_exception();
                        }
                        has_exc.store(true, std::memory_order_release);
                    }
                },
                tbb::auto_partitioner());
        });
        if (captured_exc) {
            std::rethrow_exception(captured_exc);
        }
    } else {
        ScopedWorkerSlot slot_guard;
        fn(0u, n_items, slot_guard.slot());
    }
#else
    const std::size_t n_blocks = (n_items + grain - 1u) / grain;
    typename RangeBlockOp<Fn>::SharedCtx shared_ctx;
    shared_ctx.fn = &fn;
    shared_ctx.n_items = n_items;
    shared_ctx.grain = grain;

    const madness::Range<std::size_t> block_range(0u, n_blocks, 1);
    Future<bool> done = world.taskq.for_each(block_range, RangeBlockOp<Fn>{&shared_ctx});
    thread_purge();
#ifdef HAVE_PARSEC
    const bool allow_dowork = (g_windowed_task_depth == 0u);
#else
    constexpr bool allow_dowork = true;
#endif
    (void)done.get(allow_dowork);
    if (shared_ctx.captured_exc) {
        std::rethrow_exception(shared_ctx.captured_exc);
    }
#endif
}

template <typename Fn>
inline void windowed_for_each_range(World& world,
                                    std::size_t n_items,
                                    std::size_t min_grain,
                                    const Fn& fn) {
    parallel_for_chunks(world, n_items, min_grain, fn);
}

/// Submit `n_tasks` outer orbital/pair tasks (`fn(idx)`) to MADNESS's `WorldTaskQueue` / `tbb_arena`
/// without holding a leaf `ScopedWorkerSlot` across `Future::get` / nested `Range` sweeps:
///   - On Pthreads/PaRSEC, submits normal-priority `WorldTaskQueue::add` tasks (pushed to the back of
///     `DQueue`) so inner `TaskAttributes::hipri()` `ForEachTask` children (`push_front`) always drain
///     depth-first before starting additional outer tasks, bounding concurrent workspaces to `O(N_threads)`.
///   - On Intel TBB, dispatches `madness::Range<std::size_t>(0, n_tasks, 1)` inside `ThreadPool::tbb_arena`.
template <typename TaskFn>
inline void windowed_for_each_outer_tasks(World& world,
                                          std::size_t n_tasks,
                                          const TaskFn& fn) {
    if (n_tasks == 0u) return;
    if (n_tasks == 1u || ThreadPool::size() == 0u) {
        WindowedTaskDepthScope depth_scope;
        for (std::size_t idx = 0; idx < n_tasks; ++idx) {
            fn(idx);
        }
        return;
    }

    std::atomic<bool> has_exc{false};
    std::mutex exc_mtx;
    std::exception_ptr captured_exc;

    struct UnfinishedOuterTaskDec {
        ~UnfinishedOuterTaskDec() {
            g_unfinished_outer_tasks.fetch_sub(1u, std::memory_order_release);
        }
    };

#ifdef HAVE_INTEL_TBB
    (void)world;
    if (ThreadPool::tbb_arena) {
        g_unfinished_outer_tasks.fetch_add(n_tasks, std::memory_order_relaxed);
        const madness::Range<std::size_t> task_range(0u, n_tasks, 1);
        ThreadPool::tbb_arena->execute([&]() {
            tbb::parallel_for(
                task_range,
                [&](const madness::Range<std::size_t>& r) {
                    for (std::size_t idx = r.begin(); idx != r.end(); ++idx) {
                        UnfinishedOuterTaskDec dec_guard;
                        if (has_exc.load(std::memory_order_relaxed)) continue;
                        WindowedTaskDepthScope depth_scope;
                        try {
                            fn(idx);
                        } catch (...) {
                            std::lock_guard<std::mutex> lock(exc_mtx);
                            if (!captured_exc) {
                                captured_exc = std::current_exception();
                            }
                            has_exc.store(true, std::memory_order_release);
                        }
                    }
                },
                tbb::auto_partitioner());
        });
    } else {
        WindowedTaskDepthScope depth_scope;
        for (std::size_t idx = 0; idx < n_tasks; ++idx) {
            fn(idx);
        }
    }
#else
    struct OuterSharedCtx {
        const TaskFn* fn = nullptr;
        std::atomic<bool>* has_exc = nullptr;
        std::mutex* exc_mtx = nullptr;
        std::exception_ptr* captured_exc = nullptr;
    } shared_ctx{&fn, &has_exc, &exc_mtx, &captured_exc};

    auto outer_op = [ctx = &shared_ctx](std::size_t idx) -> bool {
        UnfinishedOuterTaskDec dec_guard;
        if (ctx->has_exc->load(std::memory_order_relaxed)) return true;
        WindowedTaskDepthScope depth_scope;
        try {
            (*ctx->fn)(idx);
        } catch (...) {
            std::lock_guard<std::mutex> lock(*ctx->exc_mtx);
            if (!(*ctx->captured_exc)) {
                *ctx->captured_exc = std::current_exception();
            }
            ctx->has_exc->store(true, std::memory_order_release);
        }
        return true;
    };

    g_unfinished_outer_tasks.fetch_add(n_tasks, std::memory_order_relaxed);
    std::vector<Future<bool>> done;
    done.reserve(n_tasks);
    for (std::size_t idx = 0; idx < n_tasks; ++idx) {
        done.push_back(world.taskq.add(outer_op, idx));
    }
    thread_purge();
#ifdef HAVE_PARSEC
    const bool allow_dowork = (g_windowed_task_depth == 0u);
#else
    constexpr bool allow_dowork = true;
#endif
    for (auto& fut : done) {
        (void)fut.get(allow_dowork);
    }
#endif

    if (captured_exc) {
        std::rethrow_exception(captured_exc);
    }
}

} // namespace windowed_parallel_detail

/// Bitmask flags classifying virtual halo nodes at level `j`.
enum WindowedHaloFlags : std::uint8_t {
    kHaloNone = 0u,
    kHaloH1T = 1u << 0,                  ///< Primary virtual halo dH_j^{(1,T)} around real nodes T_j (Stage 4a).
    kHaloH1Gamma = 1u << 1,              ///< Primary virtual halo dH_j^{(1,Gamma)} around T_j^{int} (Pass 4b-1 ds_corr).
    kHaloH1Refine = 1u << 2,             ///< Virtual dynamic-refinement candidate around T_j^{int} (Pass 4b-2).
    kHaloH2Ghost = 1u << 3,              ///< Secondary read-only virtual ghost dH_j^{(2)} completing scaling stencils.
    kHaloH1RefineLeafUngraded = 1u << 4  ///< Real-leaf dynamic-refinement gate: ungraded jump (dj >= 2) across neighbor.
};

namespace windowed_halo_detail {

inline double scalar_abs2(double x) noexcept { return x * x; }

inline double scalar_abs2(const double_complex& z) noexcept {
    const double re = z.real();
    const double im = z.imag();
    return re * re + im * im;
}

} // namespace windowed_halo_detail

/// Per-level active operator half-bandwidths driving minimal fine-driven halo construction.
struct WindowedLevelRadii {
    bool has_T_terms = false;
    int m_max_T = -1;          ///< Active Stage 4(a) half-bandwidth in `{-1, 1, 2}`.
    bool has_4b1_terms = false;
    int m_max_4b1 = -1;        ///< Active Pass 4b-1 `(0, d) -> s` half-bandwidth in `{-1, 0, 1, 2}`.
    bool has_4b2_int_terms = false;
    int m_max_4b2_int = -1;    ///< Active Pass 4b-2 `(s, d) -> d` half-bandwidth in `{-1, 0, 1, 2}`.
    bool has_4b2_leaf_terms = false;
    int m_max_4b2_leaf = -1;   ///< Active Pass 4b-2 `(s, 0) -> d` half-bandwidth in `{-1, 1, 2}`.

    double sh1_s_T_thresh = 0.0;              ///< Coarse-leaf `||s_virt||_F` threshold needed for shell-1 (`sh == 1`) `H2_Ghost` around `H1_T`.
    std::array<double, 4> sh1_T_l1_op{};      ///< Exact 3D `T` operator norm by Manhattan distance `l1 = |dx|+|dy|+|dz| in {1,2,3}` (face, edge, corner).
    double sh1_T_const_op = 0.0;              ///< Level-`j` Theorem-3.1 constant-mode remainder norm for `kHaloH1T`.
    double sh1_T_osc_op = 0.0;                ///< Level-`j` shell-1 oscillatory + gradient-jump operator norm bound for `kHaloH1T`.
    double sh1_T_tol = 0.0;                   ///< Per-node tolerance gating `kHaloH1T` emission.
    double sh1_4b1_thresh = 0.0;              ///< Per-node `||d_int||_F` threshold needed for shell-1 (`sh == 1`) Pass 4b-1 `kHaloH1Gamma` halo.
    double sh1_s_leaf_thresh = 0.0;           ///< Coarse-leaf `||s_virt||_F` threshold needed for shell-1 (`sh == 1`) Pass 4b-2 leaf `H2_Ghost`.
    double sh1_refine_s_op = 0.0;             ///< Level-`j` shell-1 unsubtracted scaling-to-wavelet operator norm bound for `kHaloH1Refine`.
    double sh1_refine_const_op = 0.0;         ///< Level-`j` shell-1 constant-mode wavelet leakage norm for `kHaloH1Refine`.
    double sh1_refine_osc_op = 0.0;           ///< Level-`j` shell-1 oscillatory + gradient-jump scaling-to-wavelet norm for `kHaloH1Refine`.
    double sh1_refine_d_op = 0.0;             ///< Level-`j` shell-1 wavelet-to-wavelet operator norm bound for `kHaloH1Refine`.
    double sh1_refine_tol = 0.0;              ///< Per-node wavelet output tolerance gating `kHaloH1Refine` / `near_internal_4b2` emission.

    double sh2_s_thresh = 0.0;       ///< Per-node `||s||_F` threshold needed for shell-2 (`sh == 2`) scaling halo.
    double sh2_4b1_thresh = 0.0;     ///< Per-node `||d_int||_F` threshold needed for shell-2 (`sh == 2`) Pass 4b-1 halo.
    double sh2_refine_thresh = 0.0;  ///< Per-node `||(s, d)||_F` threshold needed for shell-2 (`sh == 2`) Pass 4b-2 refine halo.

    int max_refine_cand_bw() const noexcept {
        return std::max(m_max_4b2_int, m_max_4b2_leaf);
    }

    /// Returns `{mT_inf, max_l1}` for Stage 4(a) `kHaloH1T` emission around real node `rn`,
    /// where `mT_inf in {-1, 1, 2}` is the Chebyshev box radius and `max_l1 in {0, 1, 2, 3, 6}`
    /// bounds `|dx| + |dy| + |dz|` (pruning weak edge/corner neighbors when only faces exceed `sh1_T_tol`).
    template <typename ScalarT>
    std::pair<int, int> node_mT_and_l1(bool dyn_refine,
                                       double s_normf,
                                       ScalarT s000,
                                       double s_osc_norm2) const noexcept {
        if (!has_T_terms || m_max_T < 1) return {-1, 0};
        if (!dyn_refine) return {m_max_T, 3 * m_max_T};
        if (m_max_T >= 2 && s_normf > sh2_s_thresh) return {m_max_T, 3 * m_max_T};

        const double s_dc = std::sqrt(windowed_halo_detail::scalar_abs2(s000));
        const double s_osc = std::sqrt(std::max(0.0, s_osc_norm2));
        const double csub_ub = s_dc * sh1_T_const_op + s_osc * sh1_T_osc_op;
        if (!(csub_ub > sh1_T_tol)) return {-1, 0};

        int max_l1 = 0;
        for (int l1 = 3; l1 >= 1; --l1) {
            if (s_normf * sh1_T_l1_op[static_cast<std::size_t>(l1)] > sh1_T_tol) {
                max_l1 = l1;
                break;
            }
        }
        if (max_l1 == 0) return {-1, 0};
        return {1, max_l1};
    }

    int node_m4b1(bool dyn_refine, int m1_lvl, double d_normf) const noexcept {
        if (m1_lvl < 0) return -1;
        if (!dyn_refine) return m1_lvl;
        if (m1_lvl >= 2 && d_normf > sh2_4b1_thresh) return m1_lvl;
        if (m1_lvl >= 1 && d_normf > sh1_4b1_thresh) return 1;
        return 0;
    }

    template <typename ScalarT>
    int node_mRefine(bool dyn_refine,
                     int m2_lvl,
                     double s_normf,
                     ScalarT s000,
                     double s_osc_norm2,
                     double d_normf) const noexcept {
        if (m2_lvl < 1) return -1;
        if (!dyn_refine) return m2_lvl;
        const double s_dc = std::sqrt(windowed_halo_detail::scalar_abs2(s000));
        const double s_osc = std::sqrt(std::max(0.0, s_osc_norm2));
        const double s_contrib =
            std::min(s_normf * sh1_refine_s_op,
                     s_dc * sh1_refine_const_op + s_osc * sh1_refine_osc_op);
        if (s_contrib + d_normf * sh1_refine_d_op > sh1_refine_tol) {
            return m2_lvl;
        }
        return -1;
    }

    int h2_mT(bool dyn_refine, double max_cl_above) const noexcept {
        if (!has_T_terms || m_max_T < 1) return -1;
        if (!dyn_refine) return m_max_T;
        if (m_max_T >= 2 && max_cl_above > sh2_s_thresh) return m_max_T;
        if (max_cl_above > sh1_s_T_thresh) return 1;
        return -1;
    }

    int h2_mLeaf(bool dyn_refine, bool has_int, double max_cl_above) const noexcept {
        if (!has_int || !has_4b2_leaf_terms || m_max_4b2_leaf < 1) return -1;
        if (!dyn_refine) return m_max_4b2_leaf;
        if (m_max_4b2_leaf >= 2 && max_cl_above > sh2_s_thresh) return m_max_4b2_leaf;
        if (max_cl_above > sh1_s_leaf_thresh) return 1;
        return -1;
    }
};

/// Globally reduced, rank-independent per-level input coefficient norm bounds.
struct WindowedGlobalLevelStats {
    bool has_internal = false;
    double max_s_real = 0.0;             ///< Max `||s_in||_F` over all real nodes at level `j`.
    double max_s_leaf = 0.0;             ///< Max `||s_in||_F` over real leaves at level `j`.
    double max_s_leaf_linf = 0.0;        ///< Exact Legendre Cauchy-Schwarz `L_inf` bound `|s000| + k^3 ||s_osc||_F` over real leaves at level `j`.
    double max_coarse_leaf_above = 0.0;  ///< Exact volume-scaled bound on prolonged virtual `||s_virt||_F` at level `j`.
    double max_s_source = 0.0;           ///< Exact bound `max(max_s_real[j], max_coarse_leaf_above[j])` covering all sources.
    double max_d_int = 0.0;              ///< Max `||(0, d_in)||_F` over internal nodes at level `j`.
    double max_ns_int = 0.0;             ///< Max `||(s_in, d_in)||_F` over internal nodes at level `j`.
};

namespace windowed_halo_detail {

/// Deterministic strict weak ordering over `Key<3>` independent of hash bucket layout.
inline bool windowed_key_less(const Key<3>& a, const Key<3>& b) noexcept {
    if (a.level() != b.level()) return a.level() < b.level();
    const auto& la = a.translation();
    const auto& lb = b.translation();
    if (la[0] != lb[0]) return la[0] < lb[0];
    if (la[1] != lb[1]) return la[1] < lb[1];
    return la[2] < lb[2];
}

/// Maximum tree level `j` for 64-bit packed spatial key indexing (`21` bits/axis -> `0 <= l[d] < 2^20 < 2^21`).
/// Levels `j > 20` fall back cleanly to full `Key<3>` comparisons and `shift_key_in_domain`.
inline constexpr int kMaxPackedKeyLevel = 20;
inline constexpr std::uint64_t kPackedCoordMask = (std::uint64_t(1) << 21) - 1u;

/// Pack `(lx, ly, lz)` (`0 <= lx, ly, lz < 2^21`) into an unsigned 64-bit integer whose `<` order
/// is 100% identical to `windowed_key_less` at any fixed tree level `j <= kMaxPackedKeyLevel`.
inline std::uint64_t pack_spatial_coords(Translation lx,
                                         Translation ly,
                                         Translation lz) noexcept {
    return (static_cast<std::uint64_t>(lx) << 42) |
           (static_cast<std::uint64_t>(ly) << 21) |
           static_cast<std::uint64_t>(lz);
}

inline std::uint64_t pack_spatial_key_same_level(const Key<3>& key) noexcept {
    const auto& l = key.translation();
    return pack_spatial_coords(l[0], l[1], l[2]);
}

inline Key<3> unpack_spatial_key_same_level(Level n, std::uint64_t code) noexcept {
    Vector<Translation, 3> nl;
    nl[0] = static_cast<Translation>((code >> 42) & kPackedCoordMask);
    nl[1] = static_cast<Translation>((code >> 21) & kPackedCoordMask);
    nl[2] = static_cast<Translation>(code & kPackedCoordMask);
    return Key<3>(n, nl);
}

/// Fast zero-`Key<3>`-construction neighbor coordinate shift at level `j <= kMaxPackedKeyLevel` (`dom = 2^j`).
inline bool shift_packed_coords_in_domain(Translation dom,
                                          Translation lx,
                                          Translation ly,
                                          Translation lz,
                                          int dx,
                                          int dy,
                                          int dz,
                                          const array_of_bools<3>& is_periodic,
                                          std::uint64_t& out_code) noexcept {
    Translation cx = lx + static_cast<Translation>(dx);
    if (cx < 0 || cx >= dom) {
        if (!is_periodic[0]) return false;
        cx %= dom;
        if (cx < 0) cx += dom;
    }
    Translation cy = ly + static_cast<Translation>(dy);
    if (cy < 0 || cy >= dom) {
        if (!is_periodic[1]) return false;
        cy %= dom;
        if (cy < 0) cy += dom;
    }
    Translation cz = lz + static_cast<Translation>(dz);
    if (cz < 0 || cz >= dom) {
        if (!is_periodic[2]) return false;
        cz %= dom;
        if (cz < 0) cz += dom;
    }
    out_code = pack_spatial_coords(cx, cy, cz);
    return true;
}

inline bool shift_packed_code_in_domain(Translation dom,
                                        std::uint64_t code,
                                        int dx,
                                        int dy,
                                        int dz,
                                        const array_of_bools<3>& is_periodic,
                                        std::uint64_t& out_code) noexcept {
    const Translation lx = static_cast<Translation>((code >> 42) & kPackedCoordMask);
    const Translation ly = static_cast<Translation>((code >> 21) & kPackedCoordMask);
    const Translation lz = static_cast<Translation>(code & kPackedCoordMask);
    return shift_packed_coords_in_domain(dom, lx, ly, lz, dx, dy, dz, is_periodic, out_code);
}

struct PackedKeyFlag {
    std::uint64_t code = 0;
    std::uint8_t flags = 0;
};

inline void sort_and_merge_packed_flags(std::vector<PackedKeyFlag>& vec) {
    if (vec.size() <= 1) return;
    std::sort(vec.begin(), vec.end(), [](const PackedKeyFlag& a, const PackedKeyFlag& b) noexcept {
        return a.code < b.code;
    });
    std::size_t write_idx = 0;
    for (std::size_t read_idx = 1; read_idx < vec.size(); ++read_idx) {
        if (vec[read_idx].code == vec[write_idx].code) {
            vec[write_idx].flags =
                static_cast<std::uint8_t>(vec[write_idx].flags | vec[read_idx].flags);
        } else {
            ++write_idx;
            if (write_idx != read_idx) {
                vec[write_idx] = vec[read_idx];
            }
        }
    }
    vec.resize(write_idx + 1);
}

inline bool contains_sorted_packed_flag(const std::vector<PackedKeyFlag>& vec,
                                        std::uint64_t code) noexcept {
    const auto it = std::lower_bound(
        vec.begin(), vec.end(), code, [](const PackedKeyFlag& pf, std::uint64_t c) noexcept {
            return pf.code < c;
        });
    return (it != vec.end() && it->code == code);
}

inline void sort_unique_packed_codes(std::vector<std::uint64_t>& codes) {
    if (codes.size() <= 1) return;
    std::sort(codes.begin(), codes.end());
    codes.erase(std::unique(codes.begin(), codes.end()), codes.end());
}

/// Overflow-safe neighbor computation at level `j` (`0 <= j < MAXLEVEL`).
/// Returns `true` and populates `out_key` iff `key + (dx, dy, dz)` lies inside `[0, 1]^3`
/// (wrapping modulo `2^j` along periodic axes).
inline bool shift_key_in_domain(const Key<3>& key,
                                int dx,
                                int dy,
                                int dz,
                                const array_of_bools<3>& is_periodic,
                                Key<3>& out_key) noexcept {
    const Level n = key.level();
    const Translation dom = (Translation(1) << n);
    const auto& l = key.translation();
    const int deltas[3] = {dx, dy, dz};
    Vector<Translation, 3> nl;

    for (std::size_t d = 0; d < 3; ++d) {
        Translation c = l[d] + static_cast<Translation>(deltas[d]);
        if (c < 0 || c >= dom) {
            if (!is_periodic[d]) {
                return false;
            }
            c %= dom;
            if (c < 0) c += dom;
        }
        nl[d] = c;
    }
    out_key = Key<3>(n, nl);
    return true;
}

inline int child_octant_index(const Key<3>& child_key) noexcept {
    const auto& l = child_key.translation();
    return (static_cast<int>(l[0] & Translation(1)) << 2) |
           (static_cast<int>(l[1] & Translation(1)) << 1) |
           static_cast<int>(l[2] & Translation(1));
}

/// Return true iff an internal node at level `j` with 8-bit `internal_child_mask` (bit `oct` set iff
/// child `oct` at `j + 1` is internal) has at least one internal child touching the `sh == 1` neighbor
/// at displacement `(dx, dy, dz) in {-1, 0, +1}^3 \ {(0, 0, 0)}` (i.e. an ungraded jump `dj >= 2`).
inline bool has_ungraded_child_toward(std::uint8_t internal_child_mask,
                                      int dx,
                                      int dy,
                                      int dz) noexcept {
    std::uint8_t mask = internal_child_mask;
    if (mask == 0u) return false;
    if (dx < 0) mask &= 0x0Fu;
    else if (dx > 0) mask &= 0xF0u;
    if (dy < 0) mask &= 0x33u;
    else if (dy > 0) mask &= 0xCCu;
    if (dz < 0) mask &= 0x55u;
    else if (dz > 0) mask &= 0xAAu;
    return (mask != 0u);
}

/// Compute `s000 = s(0,0,0)`, cancellation-free `s_osc_norm2 = \sum_{(i,j,l) != (0,0,0)} |s(i,j,l)|^2`,
/// and `s_normf = ||s||_F` 100% read-only on a contiguous `k^3` tensor via BLAS `cblas::dot`.
template <typename T>
inline void compute_scaling_norms_readonly(const Tensor<T>& s_cube,
                                           T& out_s000,
                                           double& out_s_osc_norm2,
                                           double& out_s_normf) {
    const long sz = s_cube.size();
    if (sz <= 0) {
        out_s000 = T(0);
        out_s_osc_norm2 = 0.0;
        out_s_normf = 0.0;
        return;
    }
    MADNESS_CHECK(s_cube.iscontiguous());
    const T* p = s_cube.ptr();
    out_s000 = p[0];
    double osc_sq = 0.0;
    if (sz > 1) {
        if constexpr (std::is_same_v<T, double>) {
            osc_sq = madness::cblas::dot(static_cast<integer>(sz - 1), p + 1, 1, p + 1, 1);
        } else {
            const double* p_re = reinterpret_cast<const double*>(p + 1);
            osc_sq = madness::cblas::dot(static_cast<integer>(2L * (sz - 1)), p_re, 1, p_re, 1);
        }
    }
    out_s_osc_norm2 = std::max(0.0, osc_sq);
    out_s_normf = std::sqrt(scalar_abs2(out_s000) + out_s_osc_norm2);
}

/// Fast BLAS `cblas::dot` check returning true iff contiguous `s_cube` has non-zero Frobenius norm.
template <typename T>
inline bool is_cube_nonzero_readonly(const Tensor<T>& s_cube) noexcept {
    const long sz = s_cube.size();
    if (sz <= 0) return false;
    const T* p = s_cube.ptr();
    double sq = 0.0;
    if constexpr (std::is_same_v<T, double>) {
        sq = madness::cblas::dot(static_cast<integer>(sz), p, 1, p, 1);
    } else {
        const double* p_re = reinterpret_cast<const double*>(p);
        sq = madness::cblas::dot(static_cast<integer>(2L * sz), p_re, 1, p_re, 1);
    }
    return sq > 0.0;
}

/// Extract a contiguous `k x k x k` scaling tensor from a `(2k)^3` or `k^3` tensor using `Slice`.
template <typename T>
inline Tensor<T> extract_scaling_cube(const Tensor<T>& full_coeff, int k) {
    const long kk = static_cast<long>(k);
    MADNESS_CHECK(full_coeff.ndim() == 3);
    if (full_coeff.dim(0) == kk && full_coeff.dim(1) == kk && full_coeff.dim(2) == kk &&
        full_coeff.iscontiguous()) {
        return copy(full_coeff);
    }
    MADNESS_CHECK(full_coeff.dim(0) >= kk && full_coeff.dim(1) >= kk && full_coeff.dim(2) >= kk);
    const Slice s0(0, kk - 1);
    return copy(full_coeff(s0, s0, s0));
}

/// Build `(0, d_int)` (`(2k)^3` tensor with `s0 = [0..k)^3` zeroed) and return `||(0, d_int)||_F`.
template <typename T>
inline Tensor<T> make_pure_wavelet_cube(const Tensor<T>& ns_coeff, int k, double* out_d_normf = nullptr) {
    const long kk = static_cast<long>(k);
    const long twok = 2L * kk;
    MADNESS_CHECK(ns_coeff.ndim() == 3 &&
                  ns_coeff.dim(0) == twok &&
                  ns_coeff.dim(1) == twok &&
                  ns_coeff.dim(2) == twok);
    const Slice s0(0, kk - 1);
    Tensor<T> wav = copy(ns_coeff);
    wav(s0, s0, s0) = T(0);
    if (out_d_normf != nullptr) {
        *out_d_normf = wav.normf();
    }
    return wav;
}

/// Single-step 3D interior prolongation (`parent -> child`) via `general_fast_transform`.
/// Zero exterior polynomial extrapolation; zero heap allocations; zero C++ element loops.
template <typename T>
inline void prolong_step_3d(const Tensor<double> (&h_oct)[8][3],
                            const Tensor<T>& parent_s,
                            const Key<3>& child_key,
                            Tensor<T>& child_s_out,
                            Tensor<T>& work1) {
    const int oct = child_octant_index(child_key);
    general_fast_transform(parent_s, h_oct[oct], child_s_out, work1);
}

/// Single-step 3D adjoint scaling restriction (`child -> parent`) via `general_fast_transform`
/// + MKL `aligned_axpy` (`parent_ds_accum += (H^T)^{\otimes 3} child_ds`).
template <typename T>
inline void restrict_step_accumulate_3d(long kk3,
                                        const Tensor<double> (&hT_oct)[8][3],
                                        const Tensor<T>& child_ds,
                                        const Key<3>& child_key,
                                        Tensor<T>& parent_ds_accum,
                                        Tensor<T>& work1,
                                        Tensor<T>& work2) {
    const int oct = child_octant_index(child_key);
    general_fast_transform(child_ds, hT_oct[oct], work2, work1);
    aligned_axpy(kk3, parent_ds_accum.ptr(), work2.ptr(), 1.0);
}

/// Compute Frobenius norm of the 7 wavelet octants (`[0..2k)^3 \ [0..k)^3`) of a `(2k)^3` NS tensor
/// via 3 disjoint slice norms (zero heap allocation, zero `s0` subtraction cancellation).
template <typename T>
inline double wavelet_octants_normf(const Tensor<T>& ns_tensor, int k) {
    const long kk = static_cast<long>(k);
    const long twok = 2L * kk;
    const Slice all(0, twok - 1);
    const Slice s0(0, kk - 1);
    const Slice s1(kk, twok - 1);
    const double n1 = ns_tensor(s1, all, all).normf();
    const double n2 = ns_tensor(s0, s1, all).normf();
    const double n3 = ns_tensor(s0, s0, s1).normf();
    return std::sqrt(n1 * n1 + n2 * n2 + n3 * n3);
}

inline void sort_and_merge_key_flags(std::vector<std::pair<Key<3>, std::uint8_t>>& vec) {
    if (vec.size() <= 1) return;
    std::sort(vec.begin(), vec.end(), [](const auto& a, const auto& b) noexcept {
        return windowed_key_less(a.first, b.first);
    });
    std::size_t write_idx = 0;
    for (std::size_t read_idx = 1; read_idx < vec.size(); ++read_idx) {
        if (vec[read_idx].first == vec[write_idx].first) {
            vec[write_idx].second =
                static_cast<std::uint8_t>(vec[write_idx].second | vec[read_idx].second);
        } else {
            ++write_idx;
            if (write_idx != read_idx) {
                vec[write_idx] = vec[read_idx];
            }
        }
    }
    vec.resize(write_idx + 1);
}

inline bool contains_sorted_key_flag(const std::vector<std::pair<Key<3>, std::uint8_t>>& vec,
                                     const Key<3>& key) noexcept {
    const auto it = std::lower_bound(
        vec.begin(), vec.end(), key, [](const auto& kv, const Key<3>& k) noexcept {
            return windowed_key_less(kv.first, k);
        });
    return (it != vec.end() && it->first == key);
}

inline void sort_unique_keys(std::vector<Key<3>>& keys) {
    if (keys.size() <= 1) return;
    std::sort(keys.begin(), keys.end(), windowed_key_less);
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
}

// --- MPI Active-Message Payload Types ---

struct HaloCandidateMsg {
    Key<3> key;
    std::uint8_t flags = 0;
    std::uint16_t slot = 0;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & key & flags & slot;
    }
};

struct KeyQueryMsg {
    Key<3> key;
    std::uint8_t need_internal_ns = 0;
    std::uint16_t slot = 0;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & key & need_internal_ns & slot;
    }
};

template <typename T>
struct KeyReplyMsg {
    Key<3> key;
    std::uint8_t is_internal = 0;
    std::uint16_t slot = 0;
    Tensor<T> coeff;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & key & is_internal & slot & coeff;
    }
};

template <typename T>
struct CoarseLeafDeltaMsg {
    Key<3> leaf_key;
    Key<3> sort_tag_key;
    Tensor<T> delta_s;
    std::uint16_t slot = 0;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & leaf_key & sort_tag_key & delta_s & slot;
    }
};

template <typename T>
struct VirtualSubtreeItemMsg {
    Key<3> coarse_leaf_key;
    Key<3> virt_key;
    std::uint8_t is_promoted = 0;
    Tensor<T> ds_corr;
    Tensor<T> d_virt_ns;
    std::uint16_t slot = 0;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & coarse_leaf_key & virt_key & is_promoted & ds_corr & d_virt_ns & slot;
    }
};

template <typename T>
struct OutputNodeInsertMsg {
    Key<3> key;
    std::uint8_t has_children = 0;
    Tensor<T> coeff;
    std::uint16_t slot = 0;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & key & has_children & coeff & slot;
    }
};

/// Dedicated `WorldObject` handling single-function and multi-orbital lockstep batched
/// multi-rank halo discovery, ghost/ancestor fetches, and Stage 5a restriction/promotion assembly.
template <typename T>
class WindowedNSHaloExchange : public WorldObject<WindowedNSHaloExchange<T>> {
public:
    using baseT = WorldObject<WindowedNSHaloExchange<T>>;

    explicit WindowedNSHaloExchange(World& world, std::size_t num_slots = 1)
        : baseT(world)
        , num_slots_(std::max<std::size_t>(1, num_slots))
        , halo_candidates_(num_slots_)
        , key_replies_(num_slots_)
        , coarse_deltas_(num_slots_)
        , subtree_items_(num_slots_)
        , output_inserts_(num_slots_) {
        this->process_pending();
    }

    ~WindowedNSHaloExchange() override = default;

    std::size_t num_slots() const noexcept { return num_slots_; }

    void recv_halo_candidates(const std::vector<HaloCandidateMsg>& msgs) {
        std::lock_guard<std::mutex> lock(mtx_);
        for (const auto& m : msgs) {
            if (m.slot < num_slots_) {
                halo_candidates_[m.slot].push_back(m);
            }
        }
    }

    void recv_key_queries(ProcessID requester, const std::vector<KeyQueryMsg>& queries) {
        std::lock_guard<std::mutex> lock(mtx_);
        key_queries_.emplace_back(requester, queries);
    }

    void recv_key_replies(const std::vector<KeyReplyMsg<T>>& replies) {
        std::lock_guard<std::mutex> lock(mtx_);
        for (const auto& r : replies) {
            if (r.slot < num_slots_) {
                key_replies_[r.slot].push_back(r);
            }
        }
    }

    void recv_coarse_deltas(const std::vector<CoarseLeafDeltaMsg<T>>& deltas) {
        std::lock_guard<std::mutex> lock(mtx_);
        for (const auto& d : deltas) {
            if (d.slot < num_slots_) {
                coarse_deltas_[d.slot].push_back(d);
            }
        }
    }

    void recv_virtual_subtree_items(const std::vector<VirtualSubtreeItemMsg<T>>& items) {
        std::lock_guard<std::mutex> lock(mtx_);
        for (const auto& it : items) {
            if (it.slot < num_slots_) {
                subtree_items_[it.slot].push_back(it);
            }
        }
    }

    void recv_output_inserts(const std::vector<OutputNodeInsertMsg<T>>& inserts) {
        std::lock_guard<std::mutex> lock(mtx_);
        for (const auto& ins : inserts) {
            if (ins.slot < num_slots_) {
                output_inserts_[ins.slot].push_back(ins);
            }
        }
    }

    std::vector<HaloCandidateMsg> take_halo_candidates(std::size_t slot = 0) {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<HaloCandidateMsg> out;
        if (slot < num_slots_) out.swap(halo_candidates_[slot]);
        return out;
    }

    std::vector<std::pair<ProcessID, std::vector<KeyQueryMsg>>> take_key_queries() {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<std::pair<ProcessID, std::vector<KeyQueryMsg>>> out;
        out.swap(key_queries_);
        return out;
    }

    std::vector<KeyReplyMsg<T>> take_key_replies(std::size_t slot = 0) {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<KeyReplyMsg<T>> out;
        if (slot < num_slots_) out.swap(key_replies_[slot]);
        return out;
    }

    std::vector<CoarseLeafDeltaMsg<T>> take_coarse_deltas(std::size_t slot = 0) {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<CoarseLeafDeltaMsg<T>> out;
        if (slot < num_slots_) out.swap(coarse_deltas_[slot]);
        return out;
    }

    std::vector<VirtualSubtreeItemMsg<T>> take_subtree_items(std::size_t slot = 0) {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<VirtualSubtreeItemMsg<T>> out;
        if (slot < num_slots_) out.swap(subtree_items_[slot]);
        return out;
    }

    std::vector<OutputNodeInsertMsg<T>> take_output_inserts(std::size_t slot = 0) {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<OutputNodeInsertMsg<T>> out;
        if (slot < num_slots_) out.swap(output_inserts_[slot]);
        return out;
    }

private:
    std::size_t num_slots_ = 1;
    std::mutex mtx_;
    std::vector<std::vector<HaloCandidateMsg>> halo_candidates_;
    std::vector<std::pair<ProcessID, std::vector<KeyQueryMsg>>> key_queries_;
    std::vector<std::vector<KeyReplyMsg<T>>> key_replies_;
    std::vector<std::vector<CoarseLeafDeltaMsg<T>>> coarse_deltas_;
    std::vector<std::vector<VirtualSubtreeItemMsg<T>>> subtree_items_;
    std::vector<std::vector<OutputNodeInsertMsg<T>>> output_inserts_;
};

} // namespace windowed_halo_detail

/// Owned real octree node record on this MPI rank at level `j`.
template <typename T>
struct WindowedOwnedRealNode {
    Key<3> key;
    bool is_internal = false;               ///< True iff `node.has_children()` in input tree.
    bool near_internal_4b1 = false;         ///< For leaves: within `m_max_4b1 >= 1` of `T_j^{int}`.
    bool near_internal_4b2 = false;         ///< For leaves: abuts ungraded (`dj >= 2`) boundary of `T_j^{int}`.
    bool promoted_leaf = false;             ///< For leaves: true iff dynamically refined in Stage 4(b).
    std::uint8_t internal_child_mask = 0u;  ///< For internals: bit `oct` set iff child `oct` at `j+1` is internal.
    T s000 = T(0);                   ///< `s_in(0, 0, 0)` constant moment.
    double s_osc_norm2 = 0.0;        ///< `||s_in - s000 e_{000}||_F^2` (cancellation-free).
    double s_normf = 0.0;            ///< `||s_in||_F`.
    double d_normf = 0.0;            ///< `||(0, d_in)||_F` (nonzero iff `is_internal`).
    double ns_normf = 0.0;           ///< `||(s_in, d_in)||_F`.
    Tensor<T> s_in;                  ///< `k^3` input scaling coefficients.
    Tensor<T> ns_in;                 ///< `(2k)^3` input NS `(s, d)` (populated iff `is_internal`).
    Tensor<T> s_out;                 ///< `k^3` accumulated scaling output `\hat{s}^j`.
    Tensor<T> ns_out;                ///< `(2k)^3` accumulated NS wavelet output `(0, \hat{d}^j)`.
};

/// Owned active `H1` virtual boundary halo node record (`dH_j^{(1,T)} \cup dH_j^{(1,Gamma)}`).
template <typename T>
struct WindowedOwnedVirtNode {
    Key<3> key;
    Key<3> coarse_leaf_key;          ///< Unique real coarse leaf `B in L` strictly containing `key`.
    std::uint8_t flags = 0;          ///< Combination of `kHaloH1T`, `kHaloH1Gamma`, `kHaloH1Refine`.
    bool promoted = false;           ///< Set true if `||d_virt||_F > eps_loc` under `dynamic_refine`.
    T s000 = T(0);                   ///< `s_virt(0, 0, 0)` constant moment.
    double s_osc_norm2 = 0.0;        ///< `||s_virt - s000 e_{000}||_F^2`.
    double s_normf = 0.0;            ///< `||s_virt||_F`.
    Tensor<T> s_virt;                ///< `k^3` prolonged virtual input scaling coefficients.
    Tensor<T> ds_corr;               ///< `k^3` spillover accumulator `\Delta\hat{s}^j_{V,corr}` (if H1T|H1Gamma).
    Tensor<T> ns_virt_out;           ///< `(2k)^3` wavelet indicator `(0, \hat{d}^j_{U,virt})` (allocated lazily).
};

/// Read-only stencil source lookup entry at level `j` (stored in a flat sorted `std::vector`).
template <typename T>
struct WindowedSourceEntry {
    Key<3> key;
    bool is_internal = false;
    bool is_virtual = false;         ///< True iff prolonged virtual halo/ghost node (`false` for real nodes).
    T s000 = T(0);
    double s_osc_norm2 = 0.0;
    double s_normf = 0.0;
    double d_normf = 0.0;
    double ns_normf = 0.0;
    Tensor<T> s_coeff;               ///< `k^3` scaling coefficients (real or prolonged virtual).
    Tensor<T> ns_coeff;              ///< `(2k)^3` NS `(s, d)` (iff `is_internal`; else empty).
};

/// Binary search `key` in sorted `sources` (`std::vector<WindowedSourceEntry<T>>`), returning `nullptr` if absent.
template <typename T>
inline const WindowedSourceEntry<T>* find_windowed_source(
    const std::vector<WindowedSourceEntry<T>>& sources,
    const Key<3>& key) noexcept {
    const auto it = std::lower_bound(
        sources.begin(), sources.end(), key,
        [](const WindowedSourceEntry<T>& se, const Key<3>& k) noexcept {
            return windowed_halo_detail::windowed_key_less(se.key, k);
        });
    if (it != sources.end() && it->key == key) {
        return &(*it);
    }
    return nullptr;
}

/// Per-level execution workspace constructed by `WindowedNSHaloWorkspace<T>`.
template <typename T>
struct WindowedLevelData {
    int level = 0;
    bool global_has_internal = false;  ///< True iff `T_j^{int}` is globally non-empty across `world`.
    std::vector<WindowedOwnedRealNode<T>> owned_real;
    std::vector<WindowedOwnedVirtNode<T>> owned_virt;
    std::vector<WindowedSourceEntry<T>> sources;        ///< Sorted by `windowed_key_less`.
    std::vector<std::uint64_t> real_packed_keys;        ///< 64-bit packed `owned_real` keys (`level <= 20`).
    std::vector<std::uint64_t> source_packed_keys;      ///< 64-bit packed `sources` keys (`level <= 20`).
};

/// Fast 64-bit packed binary search in `ld.source_packed_keys` (`ld.level <= kMaxPackedKeyLevel`).
template <typename T>
inline const WindowedSourceEntry<T>* find_windowed_source_packed(
    const WindowedLevelData<T>& ld,
    std::uint64_t code) noexcept {
    const auto& pvec = ld.source_packed_keys;
    const auto it = std::lower_bound(pvec.begin(), pvec.end(), code);
    if (it != pvec.end() && *it == code) {
        return &(ld.sources[static_cast<std::size_t>(it - pvec.begin())]);
    }
    return nullptr;
}

/// Active internal-node level bounds for zero-mutation serial (`world.size() == 1`) input ingestion:
///   - `min_s_int_level` / `max_s_int_level`: smallest / largest level `j in [0, n_max_tree - 1]`
///     that reads `s_in` (or `(s_in, d_in)`) on internal nodes (`> max_s_int_level` / `< 0` if none).
///   - `min_gamma_level` / `max_gamma_level`: smallest / largest level `j in [0, n_max_tree - 1]`
///     that reads `(2k)^3` `ns_in` / `wav_in` on internal nodes (`> max_gamma_level` / `< 0` if none).
struct WindowedSerialInputRange {
    int min_s_int_level = 0;
    int max_s_int_level = 1000;
    int min_gamma_level = 0;
    int max_gamma_level = 1000;
    const std::vector<double>* subtree_d_tol = nullptr;
};

/// Complete Stage 1b + Stage 5a distributed halo workspace for one `apply_windowed_ns_3d` call
/// (or one slot inside a multi-orbital lockstep distributed batch).
template <typename T>
class WindowedNSHaloWorkspace {
public:
    using RadiiBuilderFn = std::function<void(const std::vector<WindowedGlobalLevelStats>&,
                                              std::vector<WindowedLevelRadii>&)>;

    /// Deferred-distributed tag constructor for lockstep multi-orbital batch execution:
    /// initializes local structures without registering a per-orbital `WorldObject` or running collectives.
    struct DeferredDistributedTag {};

    WindowedNSHaloWorkspace(DeferredDistributedTag,
                            World& world,
                            const Function<T, 3>& f_in,
                            int n_max_tree,
                            bool dynamic_refine,
                            int max_refine_level,
                            bool ungraded_leaf_refine_only,
                            windowed_halo_detail::WindowedNSHaloExchange<T>* shared_exchange,
                            std::uint16_t slot)
        : world_(world)
        , k_(f_in.k())
        , n_max_tree_(n_max_tree)
        , dynamic_refine_(dynamic_refine)
        , max_refine_level_(max_refine_level)
        , ungraded_leaf_refine_only_(ungraded_leaf_refine_only)
        , cdata_(FunctionCommonData<T, 3>::get(f_in.k()))
        , is_periodic_(FunctionDefaults<3>::get_bc().is_periodic())
        , pmap_(f_in.get_pmap())
        , exchange_ptr_(shared_exchange)
        , slot_(slot) {
        MADNESS_CHECK(world_.size() > 1);
        MADNESS_CHECK(f_in.get_impl()->get_tree_state() == nonstandard_with_leaves);
        MADNESS_CHECK(exchange_ptr_ != nullptr);
        init_octant_filters();
        levels_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        global_stats_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        radii_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        coarse_leaves_by_level_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        for (int j = 0; j <= n_max_tree_; ++j) {
            levels_[static_cast<std::size_t>(j)].level = j;
        }
    }

    WindowedNSHaloWorkspace(World& world,
                            const Function<T, 3>& f_in,
                            int n_max_tree,
                            bool dynamic_refine,
                            int max_refine_level,
                            const RadiiBuilderFn& radii_builder,
                            const WindowedSerialInputRange& serial_range = WindowedSerialInputRange{},
                            bool ungraded_leaf_refine_only = false)
        : world_(world)
        , k_(f_in.k())
        , n_max_tree_(n_max_tree)
        , dynamic_refine_(dynamic_refine)
        , max_refine_level_(max_refine_level)
        , ungraded_leaf_refine_only_(ungraded_leaf_refine_only)
        , serial_range_(serial_range)
        , cdata_(FunctionCommonData<T, 3>::get(f_in.k()))
        , is_periodic_(FunctionDefaults<3>::get_bc().is_periodic())
        , pmap_(f_in.get_pmap()) {
        const TreeState in_state = f_in.get_impl()->get_tree_state();
        if (world_.size() > 1) {
            MADNESS_CHECK(in_state == nonstandard_with_leaves);
        } else {
            MADNESS_CHECK(in_state == nonstandard_with_leaves ||
                          in_state == reconstructed ||
                          in_state == compressed);
        }
        init_octant_filters();
        levels_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        global_stats_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        radii_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        coarse_leaves_by_level_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        for (int j = 0; j <= n_max_tree_; ++j) {
            levels_[static_cast<std::size_t>(j)].level = j;
        }

        if (world_.size() > 1) {
            owned_exchange_ =
                std::make_unique<windowed_halo_detail::WindowedNSHaloExchange<T>>(world_, 1u);
            exchange_ptr_ = owned_exchange_.get();
            slot_ = 0u;
            build_distributed(f_in, radii_builder);
        } else {
            build_serial(f_in, radii_builder);
        }
    }

    WindowedNSHaloWorkspace(World& world,
                            const Function<T, 3>& f_ns,
                            const std::vector<WindowedLevelRadii>& fixed_radii,
                            int n_max_tree,
                            bool dynamic_refine,
                            int max_refine_level)
        : WindowedNSHaloWorkspace(
              world,
              f_ns,
              n_max_tree,
              dynamic_refine,
              max_refine_level,
              [&fixed_radii](const std::vector<WindowedGlobalLevelStats>&,
                             std::vector<WindowedLevelRadii>& out_radii) {
                  out_radii = fixed_radii;
              }) {}

    ~WindowedNSHaloWorkspace() {
        if (owned_exchange_) {
            world_.gop.fence();
            owned_exchange_.reset();
            exchange_ptr_ = nullptr;
        }
    }

    WindowedNSHaloWorkspace(const WindowedNSHaloWorkspace&) = delete;
    WindowedNSHaloWorkspace& operator=(const WindowedNSHaloWorkspace&) = delete;

    World& world() const noexcept { return world_; }
    int k() const noexcept { return k_; }
    int n_max_tree() const noexcept { return n_max_tree_; }
    const array_of_bools<3>& is_periodic() const noexcept { return is_periodic_; }
    const FunctionCommonData<T, 3>& cdata() const noexcept { return cdata_; }
    const std::vector<WindowedGlobalLevelStats>& global_stats() const noexcept { return global_stats_; }
    const std::vector<WindowedLevelRadii>& radii() const noexcept { return radii_; }

    WindowedLevelData<T>& level_data(int j) {
        return levels_[static_cast<std::size_t>(j)];
    }

    const WindowedLevelData<T>& level_data(int j) const {
        return levels_[static_cast<std::size_t>(j)];
    }

    /// Eagerly release input scaling/NS tensors and sorted source tables at level `j`
    /// once both Stage 4(a) and Stage 4(b) have finished consuming level `j`.
    void release_stage4_inputs_for_level(int j) {
        if (j < 0 || j > n_max_tree_) return;
        auto& ld = levels_[static_cast<std::size_t>(j)];
        {
            std::vector<WindowedSourceEntry<T>>().swap(ld.sources);
            std::vector<std::uint64_t>().swap(ld.source_packed_keys);
        }
        for (auto& rn : ld.owned_real) {
            rn.s_in = Tensor<T>();
            rn.ns_in = Tensor<T>();
        }
        for (auto& vn : ld.owned_virt) {
            vn.s_virt = Tensor<T>();
        }
    }

    /// Fast single-rank (`world.size() == 1`) Stage 5a + Stage 5b synthesis directly into
    /// `TreeState::reconstructed` when no dynamic leaf/virtual promotion occurred:
    ///   - Restricts all virtual boundary corrections in flat arrays (`restrict_all_virtual_unpromoted`).
    ///   - Propagates non-zero internal `(s_out, ns_out)` top-down (`j = 0 .. n_max_tree - 1`) in parallel
    ///     flat arrays via `fast_transform(..., cdata_.hg)` (or `prolong_step_3d` when `ns_out` is empty),
    ///     completely skipping untouched/zero internal nodes and inactive coarse levels.
    ///   - Populates `result`'s `WorldContainer` ONCE directly in `TreeState::reconstructed` (zero
    ///     `result.reconstruct()` task-queue walk).
    /// Falls back cleanly to `assemble_output_nonstandard_tree(result)` (`return false`) when
    /// `world.size() > 1` or when any leaf/virtual node was dynamically promoted.
    bool assemble_output_reconstructed_serial(Function<T, 3>& result) {
        if (world_.size() > 1) {
            assemble_output_nonstandard_tree(result);
            return false;
        }

        long local_any_promoted = 0;
        for (int j = 0; j <= n_max_tree_; ++j) {
            for (const auto& rn : levels_[static_cast<std::size_t>(j)].owned_real) {
                if (rn.promoted_leaf && j < max_refine_level_) ++local_any_promoted;
            }
            if (j >= 1) {
                for (const auto& vn : levels_[static_cast<std::size_t>(j)].owned_virt) {
                    if (vn.promoted && j < max_refine_level_) ++local_any_promoted;
                }
            }
        }
        if (local_any_promoted > 0) {
            assemble_output_nonstandard_tree(result);
            return false;
        }

        restrict_all_virtual_unpromoted();
        for (int j = 1; j <= n_max_tree_; ++j) {
            std::vector<WindowedOwnedVirtNode<T>>().swap(levels_[static_cast<std::size_t>(j)].owned_virt);
        }

        const long kk = static_cast<long>(k_);
        const long kk3 = kk * kk * kk;
        const long twok = 2L * kk;
        const long twok3 = twok * twok * twok;
        const long k_dims[3] = {kk, kk, kk};
        const long twok_dims[3] = {twok, twok, twok};
        const Slice s0(0, kk - 1);

        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_full_2k;
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_unf_2k;
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_work_2k;
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_child_k3;
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_work_k3;

        std::vector<std::size_t> active_ints;
        for (int j = 0; j < n_max_tree_; ++j) {
            auto& ld = levels_[static_cast<std::size_t>(j)];
            auto& child_ld = levels_[static_cast<std::size_t>(j + 1)];
            if (ld.owned_real.empty()) continue;

            active_ints.clear();
            for (std::size_t ri = 0; ri < ld.owned_real.size(); ++ri) {
                const auto& rn = ld.owned_real[ri];
                if (!rn.is_internal) continue;
                if (rn.ns_out.size() == twok3 ||
                    windowed_halo_detail::is_cube_nonzero_readonly(rn.s_out)) {
                    active_ints.push_back(ri);
                }
            }
            if (active_ints.empty()) continue;

            windowed_parallel_detail::parallel_for_chunks(
                world_,
                active_ints.size(),
                2,
                [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
                    for (std::size_t ai = b; ai < e_hi; ++ai) {
                        const auto& rn = ld.owned_real[active_ints[ai]];
                        if (rn.ns_out.size() == twok3) {
                            if (tls_full_2k[cid].size() != twok3) {
                                tls_full_2k[cid] = Tensor<T>(3L, twok_dims, false);
                                tls_unf_2k[cid] = Tensor<T>(3L, twok_dims, false);
                                tls_work_2k[cid] = Tensor<T>(3L, twok_dims, false);
                            }
                            Tensor<T>& full_ns = tls_full_2k[cid];
                            Tensor<T>& unf_2k = tls_unf_2k[cid];
                            Tensor<T>& work_2k = tls_work_2k[cid];

                            std::memcpy(full_ns.ptr(),
                                        rn.ns_out.ptr(),
                                        static_cast<std::size_t>(twok3) * sizeof(T));
                            full_ns(s0, s0, s0) = rn.s_out;
                            fast_transform(full_ns, cdata_.hg, unf_2k, work_2k);

                            for (KeyChildIterator<3> kit(rn.key); kit; ++kit) {
                                const Key<3>& ch = kit.key();
                                const std::size_t cidx = find_owned_real_index(j + 1, ch);
                                MADNESS_CHECK(cidx != SIZE_MAX);
                                const auto& cl = ch.translation();
                                child_ld.owned_real[cidx].s_out +=
                                    unf_2k(cdata_.s[cl[0] & 1],
                                           cdata_.s[cl[1] & 1],
                                           cdata_.s[cl[2] & 1]);
                            }
                        } else {
                            if (tls_child_k3[cid].size() != kk3) {
                                tls_child_k3[cid] = Tensor<T>(3L, k_dims, false);
                                tls_work_k3[cid] = Tensor<T>(3L, k_dims, false);
                            }
                            Tensor<T>& child_tmp = tls_child_k3[cid];
                            Tensor<T>& work_k3 = tls_work_k3[cid];

                            for (KeyChildIterator<3> kit(rn.key); kit; ++kit) {
                                const Key<3>& ch = kit.key();
                                const std::size_t cidx = find_owned_real_index(j + 1, ch);
                                MADNESS_CHECK(cidx != SIZE_MAX);
                                windowed_halo_detail::prolong_step_3d(
                                    h_oct_, rn.s_out, ch, child_tmp, work_k3);
                                aligned_axpy(kk3,
                                             child_ld.owned_real[cidx].s_out.ptr(),
                                             child_tmp.ptr(),
                                             1.0);
                            }
                        }
                    }
                });
            for (std::size_t ai : active_ints) {
                auto& rn = ld.owned_real[ai];
                rn.ns_out = Tensor<T>();
                rn.s_out = Tensor<T>();
            }
        }

        auto* res_impl = result.get_impl().get();
        auto& res_coeffs = res_impl->get_coeffs();
        const TensorArgs res_targs = res_impl->get_tensor_args();
        const GenTensor<T> empty_coeff;
        for (int j = 0; j <= n_max_tree_; ++j) {
            auto& ld = levels_[static_cast<std::size_t>(j)];
            for (auto& rn : ld.owned_real) {
                if (rn.is_internal) {
                    res_coeffs.replace(rn.key, FunctionNode<T, 3>(empty_coeff, true));
                } else {
                    GenTensor<T> coeff(rn.s_out, res_targs);
                    rn.s_out = Tensor<T>();
                    res_coeffs.replace(rn.key, FunctionNode<T, 3>(coeff, false));
                }
            }
        }
        res_impl->set_tree_state(reconstructed);
        return true;
    }

    /// Lockstep Stage 1b distributed halo construction across `B` orbitals simultaneously,
    /// coalescing all `gop.max` reductions and active-message halo exchanges across `0 .. B-1`
    /// into 4–5 global collectives total for the entire batch.
    static void build_distributed_lockstep(
        World& world,
        const std::vector<WindowedNSHaloWorkspace<T>*>& ws_list,
        const std::vector<const Function<T, 3>*>& f_list,
        const std::vector<RadiiBuilderFn>& radii_builders,
        windowed_halo_detail::WindowedNSHaloExchange<T>& ex) {
        const std::size_t B = ws_list.size();
        if (B == 0) return;
        MADNESS_CHECK(world.size() > 1);
        MADNESS_CHECK(f_list.size() == B && radii_builders.size() == B);
        MADNESS_CHECK(ex.num_slots() >= B);

        const ProcessID me = world.rank();
        const int nproc = world.size();

        // 1. Local real-node ingestion + ONE coalesced `gop.max` across all `B` orbitals.
        std::size_t total_stat_doubles = 0;
        std::vector<std::size_t> stat_offsets(B, 0u);
        for (std::size_t b = 0; b < B; ++b) {
            auto* ws = ws_list[b];
            MADNESS_CHECK(ws != nullptr && ws->slot_ == static_cast<std::uint16_t>(b));
            ws->ingest_local_real_nodes_local_only(*f_list[b]);
            stat_offsets[b] = total_stat_doubles;
            total_stat_doubles += 6u * static_cast<std::size_t>(ws->n_max_tree_ + 1);
        }

        std::vector<double> gop_buf(total_stat_doubles, 0.0);
        for (std::size_t b = 0; b < B; ++b) {
            const auto* ws = ws_list[b];
            const std::size_t n_lvls = static_cast<std::size_t>(ws->n_max_tree_ + 1);
            const std::size_t base = stat_offsets[b];
            for (std::size_t j = 0; j < n_lvls; ++j) {
                const auto& st = ws->global_stats_[j];
                gop_buf[base + 6 * j + 0] = st.has_internal ? 1.0 : 0.0;
                gop_buf[base + 6 * j + 1] = st.max_s_real;
                gop_buf[base + 6 * j + 2] = st.max_s_leaf;
                gop_buf[base + 6 * j + 3] = st.max_s_leaf_linf;
                gop_buf[base + 6 * j + 4] = st.max_d_int;
                gop_buf[base + 6 * j + 5] = st.max_ns_int;
            }
        }
        world.gop.max(gop_buf.data(), gop_buf.size());

        bool any_need_parent_mask = false;
        for (std::size_t b = 0; b < B; ++b) {
            auto* ws = ws_list[b];
            const std::size_t n_lvls = static_cast<std::size_t>(ws->n_max_tree_ + 1);
            const std::size_t base = stat_offsets[b];
            for (std::size_t j = 0; j < n_lvls; ++j) {
                auto& st = ws->global_stats_[j];
                st.has_internal = (gop_buf[base + 6 * j + 0] > 0.5);
                st.max_s_real = gop_buf[base + 6 * j + 1];
                st.max_s_leaf = gop_buf[base + 6 * j + 2];
                st.max_s_leaf_linf = gop_buf[base + 6 * j + 3];
                st.max_d_int = gop_buf[base + 6 * j + 4];
                st.max_ns_int = gop_buf[base + 6 * j + 5];
            }
            ws->finalize_ingested_global_stats(true);
            radii_builders[b](ws->global_stats_, ws->radii_);
            if (ws->dynamic_refine_ && ws->ungraded_leaf_refine_only_ && ws->n_max_tree_ >= 1) {
                any_need_parent_mask = true;
            }
        }

        // 2. Coalesced parent-mask exchange (1 fence across all `B` orbitals; skipped completely if none need it).
        if (any_need_parent_mask) {
            std::vector<std::vector<windowed_halo_detail::HaloCandidateMsg>> parent_outbox(
                static_cast<std::size_t>(nproc));
            std::vector<KeyFlagVec> pvecs(static_cast<std::size_t>(nproc));

            for (std::size_t b = 0; b < B; ++b) {
                auto* ws = ws_list[b];
                if (!(ws->dynamic_refine_ && ws->ungraded_leaf_refine_only_ && ws->n_max_tree_ >= 1)) {
                    continue;
                }
                for (int p = 0; p < nproc; ++p) pvecs[static_cast<std::size_t>(p)].clear();
                for (int j = 1; j <= ws->n_max_tree_; ++j) {
                    for (const auto& ch : ws->levels_[static_cast<std::size_t>(j)].owned_real) {
                        if (!ch.is_internal) continue;
                        const Key<3> par = ch.key.parent();
                        const int oct = windowed_halo_detail::child_octant_index(ch.key);
                        const ProcessID owner = ws->pmap_->owner(par);
                        pvecs[static_cast<std::size_t>(owner)].emplace_back(
                            par, static_cast<std::uint8_t>(1u << oct));
                    }
                }
                const std::uint16_t slot = ws->slot_;
                for (int p = 0; p < nproc; ++p) {
                    auto& vec = pvecs[static_cast<std::size_t>(p)];
                    if (vec.empty()) continue;
                    windowed_halo_detail::sort_and_merge_key_flags(vec);
                    auto& dst = parent_outbox[static_cast<std::size_t>(p)];
                    for (const auto& kv : vec) {
                        dst.push_back(
                            windowed_halo_detail::HaloCandidateMsg{kv.first, kv.second, slot});
                    }
                }
            }

            for (int p = 0; p < nproc; ++p) {
                auto& batch = parent_outbox[static_cast<std::size_t>(p)];
                if (batch.empty()) continue;
                if (p == me) {
                    ex.recv_halo_candidates(batch);
                } else {
                    ex.send(static_cast<ProcessID>(p),
                            &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_halo_candidates,
                            batch);
                }
            }
            world.gop.fence();

            for (std::size_t b = 0; b < B; ++b) {
                auto* ws = ws_list[b];
                for (const auto& msg : ex.take_halo_candidates(ws->slot_)) {
                    const int pj = static_cast<int>(msg.key.level());
                    const std::size_t pidx = ws->find_owned_real_index(pj, msg.key);
                    if (pidx != SIZE_MAX) {
                        auto& prn = ws->levels_[static_cast<std::size_t>(pj)].owned_real[pidx];
                        prn.internal_child_mask =
                            static_cast<std::uint8_t>(prn.internal_child_mask | msg.flags);
                    }
                }
            }
        }

        // 3. Coalesced H1 candidate exchange (1 fence across all `B` orbitals).
        {
            std::vector<std::vector<windowed_halo_detail::HaloCandidateMsg>> h1_outbox(
                static_cast<std::size_t>(nproc));
            std::vector<KeyFlagVec> out_cand_vecs(static_cast<std::size_t>(nproc));

            for (std::size_t b = 0; b < B; ++b) {
                auto* ws = ws_list[b];
                for (int p = 0; p < nproc; ++p) out_cand_vecs[static_cast<std::size_t>(p)].clear();

                auto emit_candidate = [&](const Key<3>& nb, std::uint8_t flag_bits) {
                    const ProcessID owner = ws->pmap_->owner(nb);
                    out_cand_vecs[static_cast<std::size_t>(owner)].emplace_back(nb, flag_bits);
                };

                for (int j = 1; j <= ws->n_max_tree_; ++j) {
                    const auto& ld = ws->levels_[static_cast<std::size_t>(j)];
                    const WindowedLevelRadii& rad = ws->radii_[static_cast<std::size_t>(j)];
                    const bool can_refine_j = ws->dynamic_refine_ && (j < ws->max_refine_level_);

                    if (rad.has_T_terms && rad.m_max_T >= 1) {
                        for (const auto& rn : ld.owned_real) {
                            const auto [mT_rn, max_l1_T] = rad.node_mT_and_l1(
                                ws->dynamic_refine_, rn.s_normf, rn.s000, rn.s_osc_norm2);
                            if (mT_rn < 1) continue;
                            for (int dx = -mT_rn; dx <= mT_rn; ++dx) {
                                for (int dy = -mT_rn; dy <= mT_rn; ++dy) {
                                    for (int dz = -mT_rn; dz <= mT_rn; ++dz) {
                                        const int l1 = std::abs(dx) + std::abs(dy) + std::abs(dz);
                                        if (l1 == 0 || l1 > max_l1_T) continue;
                                        Key<3> nb;
                                        if (windowed_halo_detail::shift_key_in_domain(
                                                rn.key, dx, dy, dz, ws->is_periodic_, nb)) {
                                            emit_candidate(nb, kHaloH1T);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    if (ld.global_has_internal) {
                        const int m1 = rad.has_4b1_terms ? rad.m_max_4b1 : -1;
                        const int m2 = (can_refine_j && rad.max_refine_cand_bw() >= 1) ? 1 : -1;
                        if (std::max(m1, m2) >= 1) {
                            for (const auto& rn : ld.owned_real) {
                                if (!rn.is_internal) continue;
                                const int m1_rn = rad.node_m4b1(ws->dynamic_refine_, m1, rn.d_normf);
                                const int m2_rn = rad.node_mRefine(ws->dynamic_refine_,
                                                                   m2,
                                                                   rn.s_normf,
                                                                   rn.s000,
                                                                   rn.s_osc_norm2,
                                                                   rn.d_normf);
                                const int m_emit_rn = std::max(m1_rn, m2_rn);
                                if (m_emit_rn < 1) continue;
                                for (int dx = -m_emit_rn; dx <= m_emit_rn; ++dx) {
                                    for (int dy = -m_emit_rn; dy <= m_emit_rn; ++dy) {
                                        for (int dz = -m_emit_rn; dz <= m_emit_rn; ++dz) {
                                            if (dx == 0 && dy == 0 && dz == 0) continue;
                                            const int sh = std::max(
                                                {std::abs(dx), std::abs(dy), std::abs(dz)});
                                            Key<3> nb;
                                            if (windowed_halo_detail::shift_key_in_domain(
                                                    rn.key, dx, dy, dz, ws->is_periodic_, nb)) {
                                                std::uint8_t bits = 0u;
                                                if (sh <= m1_rn) {
                                                    bits = static_cast<std::uint8_t>(
                                                        bits | kHaloH1Gamma);
                                                }
                                                if (sh <= m2_rn) {
                                                    bits = static_cast<std::uint8_t>(
                                                        bits | kHaloH1Refine);
                                                    if (!ws->ungraded_leaf_refine_only_ ||
                                                        windowed_halo_detail::has_ungraded_child_toward(
                                                            rn.internal_child_mask, dx, dy, dz)) {
                                                        bits = static_cast<std::uint8_t>(
                                                            bits | kHaloH1RefineLeafUngraded);
                                                    }
                                                }
                                                if (bits != 0u) emit_candidate(nb, bits);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                const std::uint16_t slot = ws->slot_;
                for (int p = 0; p < nproc; ++p) {
                    auto& vec = out_cand_vecs[static_cast<std::size_t>(p)];
                    if (vec.empty()) continue;
                    windowed_halo_detail::sort_and_merge_key_flags(vec);
                    auto& dst = h1_outbox[static_cast<std::size_t>(p)];
                    for (const auto& kv : vec) {
                        dst.push_back(
                            windowed_halo_detail::HaloCandidateMsg{kv.first, kv.second, slot});
                    }
                }
            }

            for (int p = 0; p < nproc; ++p) {
                auto& batch = h1_outbox[static_cast<std::size_t>(p)];
                if (batch.empty()) continue;
                if (p == me) {
                    ex.recv_halo_candidates(batch);
                } else {
                    ex.send(static_cast<ProcessID>(p),
                            &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_halo_candidates,
                            batch);
                }
            }
            world.gop.fence();
        }

        std::vector<std::vector<KeyFlagVec>> virt_flags_all(B);
        std::vector<std::vector<std::vector<Key<3>>>> needed_vc_all(B);
        std::vector<std::vector<windowed_halo_detail::KeyQueryMsg>> query_outbox(
            static_cast<std::size_t>(nproc));

        // 4. Process H1 candidates & coalesce ancestor/source key queries across all `B` orbitals.
        {
            std::vector<KeyFlagVec> query_vecs(static_cast<std::size_t>(nproc));
            for (std::size_t b = 0; b < B; ++b) {
                auto* ws = ws_list[b];
                auto& virt_flags = virt_flags_all[b];
                virt_flags.resize(static_cast<std::size_t>(ws->n_max_tree_ + 1));

                for (const auto& cand : ex.take_halo_candidates(ws->slot_)) {
                    const int j = static_cast<int>(cand.key.level());
                    const std::size_t ridx = ws->find_owned_real_index(j, cand.key);
                    if (ridx != SIZE_MAX) {
                        auto& rn = ws->levels_[static_cast<std::size_t>(j)].owned_real[ridx];
                        if (!rn.is_internal) {
                            if ((cand.flags & kHaloH1Gamma) != 0) rn.near_internal_4b1 = true;
                            if ((cand.flags & kHaloH1RefineLeafUngraded) != 0) {
                                rn.near_internal_4b2 = true;
                            }
                        }
                    } else {
                        const std::uint8_t vflags = static_cast<std::uint8_t>(
                            cand.flags & static_cast<std::uint8_t>(~kHaloH1RefineLeafUngraded));
                        if (vflags != 0u) {
                            virt_flags[static_cast<std::size_t>(j)].emplace_back(cand.key, vflags);
                        }
                    }
                }
                for (int j = 1; j <= ws->n_max_tree_; ++j) {
                    windowed_halo_detail::sort_and_merge_key_flags(
                        virt_flags[static_cast<std::size_t>(j)]);
                }

                auto& needed_virt_capable = needed_vc_all[b];
                needed_virt_capable.resize(static_cast<std::size_t>(ws->n_max_tree_ + 1));
                std::vector<std::vector<Key<3>>> needed_real_only(
                    static_cast<std::size_t>(ws->n_max_tree_ + 1));

                auto add_ring_to = [&](std::vector<Key<3>>& dst, const Key<3>& center, int radius) {
                    if (radius < 0) return;
                    for (int dx = -radius; dx <= radius; ++dx) {
                        for (int dy = -radius; dy <= radius; ++dy) {
                            for (int dz = -radius; dz <= radius; ++dz) {
                                Key<3> nb;
                                if (windowed_halo_detail::shift_key_in_domain(
                                        center, dx, dy, dz, ws->is_periodic_, nb)) {
                                    dst.push_back(nb);
                                }
                            }
                        }
                    }
                };

                for (int j = 1; j <= ws->n_max_tree_; ++j) {
                    const std::size_t sj = static_cast<std::size_t>(j);
                    const auto& ld = ws->levels_[sj];
                    const WindowedLevelRadii& rad = ws->radii_[sj];
                    const auto& vf = virt_flags[sj];
                    const bool can_refine_j = ws->dynamic_refine_ && (j < ws->max_refine_level_);
                    const double max_cl_above = ws->global_stats_[sj].max_coarse_leaf_above;
                    auto& vc_vec = needed_virt_capable[sj];
                    auto& ro_vec = needed_real_only[sj];

                    if (rad.has_T_terms && rad.m_max_T >= 1) {
                        const int mT_h2 = rad.h2_mT(ws->dynamic_refine_, max_cl_above);
                        for (const auto& rn : ld.owned_real) {
                            const auto [mT_rn, max_l1_T] = rad.node_mT_and_l1(
                                ws->dynamic_refine_, rn.s_normf, rn.s000, rn.s_osc_norm2);
                            (void)max_l1_T;
                            if (mT_rn >= 1) {
                                add_ring_to(vc_vec, rn.key, mT_rn);
                            }
                            if (rad.m_max_T > std::max(0, mT_rn)) {
                                add_ring_to(ro_vec, rn.key, rad.m_max_T);
                            }
                        }
                        for (const auto& kv : vf) {
                            if ((kv.second & kHaloH1T) != 0) {
                                if (mT_h2 >= 1) {
                                    add_ring_to(vc_vec, kv.first, mT_h2);
                                }
                                if (rad.m_max_T > std::max(0, mT_h2)) {
                                    add_ring_to(ro_vec, kv.first, rad.m_max_T);
                                }
                            }
                        }
                    }

                    if (ld.global_has_internal) {
                        if (rad.has_4b1_terms && rad.m_max_4b1 >= 1) {
                            const int m1 = rad.m_max_4b1;
                            for (const auto& rn : ld.owned_real) {
                                if (rn.is_internal || rn.near_internal_4b1) {
                                    add_ring_to(ro_vec, rn.key, m1);
                                }
                            }
                            for (const auto& kv : vf) {
                                if ((kv.second & kHaloH1Gamma) != 0) {
                                    add_ring_to(ro_vec, kv.first, m1);
                                }
                            }
                        }

                        if (rad.has_4b2_int_terms && rad.m_max_4b2_int >= 1) {
                            const int mInt = rad.m_max_4b2_int;
                            for (const auto& rn : ld.owned_real) {
                                if (rn.is_internal || (can_refine_j && rn.near_internal_4b2)) {
                                    add_ring_to(ro_vec, rn.key, mInt);
                                }
                            }
                            if (can_refine_j) {
                                for (const auto& kv : vf) {
                                    if ((kv.second & kHaloH1Refine) != 0) {
                                        add_ring_to(ro_vec, kv.first, mInt);
                                    }
                                }
                            }
                        }

                        if (rad.has_4b2_leaf_terms && rad.m_max_4b2_leaf >= 1) {
                            const int mLeaf = rad.m_max_4b2_leaf;
                            const int mLeaf_h2 =
                                rad.h2_mLeaf(ws->dynamic_refine_, true, max_cl_above);
                            for (const auto& rn : ld.owned_real) {
                                if (rn.is_internal || (can_refine_j && rn.near_internal_4b2)) {
                                    if (mLeaf_h2 >= 1) {
                                        add_ring_to(vc_vec, rn.key, mLeaf_h2);
                                    }
                                    if (mLeaf > std::max(0, mLeaf_h2)) {
                                        add_ring_to(ro_vec, rn.key, mLeaf);
                                    }
                                }
                            }
                            if (can_refine_j) {
                                for (const auto& kv : vf) {
                                    if ((kv.second & kHaloH1Refine) != 0) {
                                        if (mLeaf_h2 >= 1) {
                                            add_ring_to(vc_vec, kv.first, mLeaf_h2);
                                        }
                                        if (mLeaf > std::max(0, mLeaf_h2)) {
                                            add_ring_to(ro_vec, kv.first, mLeaf);
                                        }
                                    }
                                }
                            }
                        }
                    }

                    windowed_halo_detail::sort_unique_keys(vc_vec);
                    windowed_halo_detail::sort_unique_keys(ro_vec);
                }

                for (int p = 0; p < nproc; ++p) query_vecs[static_cast<std::size_t>(p)].clear();

                auto register_query_chain = [&](const Key<3>& start_key) {
                    const int j_start = static_cast<int>(start_key.level());
                    const ProcessID s_owner = ws->pmap_->owner(start_key);
                    if (s_owner == me && ws->find_owned_real_index(j_start, start_key) != SIZE_MAX) {
                        return;
                    }
                    if (s_owner != me) {
                        query_vecs[static_cast<std::size_t>(s_owner)].emplace_back(start_key, 1u);
                    }

                    Key<3> curr = start_key;
                    while (curr.level() > 0) {
                        curr = curr.parent();
                        const int clvl = static_cast<int>(curr.level());
                        const ProcessID owner = ws->pmap_->owner(curr);
                        if (owner == me) {
                            if (ws->find_owned_real_index(clvl, curr) != SIZE_MAX) {
                                break;
                            }
                        } else {
                            query_vecs[static_cast<std::size_t>(owner)].emplace_back(curr, 0u);
                        }
                    }
                };

                for (int j = 1; j <= ws->n_max_tree_; ++j) {
                    for (const Key<3>& skey : needed_virt_capable[static_cast<std::size_t>(j)]) {
                        register_query_chain(skey);
                    }
                    for (const auto& kv : virt_flags[static_cast<std::size_t>(j)]) {
                        register_query_chain(kv.first);
                    }
                    for (const Key<3>& rkey : needed_real_only[static_cast<std::size_t>(j)]) {
                        const ProcessID owner = ws->pmap_->owner(rkey);
                        if (owner != me) {
                            query_vecs[static_cast<std::size_t>(owner)].emplace_back(rkey, 1u);
                        }
                    }
                }

                const std::uint16_t slot = ws->slot_;
                for (int p = 0; p < nproc; ++p) {
                    auto& qv = query_vecs[static_cast<std::size_t>(p)];
                    if (qv.empty()) continue;
                    windowed_halo_detail::sort_and_merge_key_flags(qv);
                    auto& dst = query_outbox[static_cast<std::size_t>(p)];
                    for (const auto& kv : qv) {
                        dst.push_back(
                            windowed_halo_detail::KeyQueryMsg{kv.first, kv.second, slot});
                    }
                }
            }

            for (int p = 0; p < nproc; ++p) {
                auto& batch = query_outbox[static_cast<std::size_t>(p)];
                if (batch.empty()) continue;
                ex.send(static_cast<ProcessID>(p),
                        &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_key_queries,
                        me,
                        batch);
            }
            world.gop.fence();
        }

        // 5. Coalesced key replies across all `B` orbitals (1 fence total).
        {
            std::vector<std::vector<windowed_halo_detail::KeyReplyMsg<T>>> reply_outbox(
                static_cast<std::size_t>(nproc));
            for (const auto& [requester, queries] : ex.take_key_queries()) {
                auto& replies = reply_outbox[static_cast<std::size_t>(requester)];
                for (const auto& q : queries) {
                    if (q.slot >= B) continue;
                    const auto* ws = ws_list[q.slot];
                    const int qlvl = static_cast<int>(q.key.level());
                    if (qlvl < 0 || qlvl > ws->n_max_tree_) continue;
                    const std::size_t ridx = ws->find_owned_real_index(qlvl, q.key);
                    if (ridx == SIZE_MAX) continue;

                    const auto& rn = ws->levels_[static_cast<std::size_t>(qlvl)].owned_real[ridx];
                    if (!rn.is_internal) {
                        replies.push_back(
                            windowed_halo_detail::KeyReplyMsg<T>{rn.key, 0u, q.slot, rn.s_in});
                    } else if (q.need_internal_ns != 0) {
                        replies.push_back(
                            windowed_halo_detail::KeyReplyMsg<T>{rn.key, 1u, q.slot, rn.ns_in});
                    }
                }
            }
            for (int p = 0; p < nproc; ++p) {
                auto& batch = reply_outbox[static_cast<std::size_t>(p)];
                if (batch.empty()) continue;
                ex.send(static_cast<ProcessID>(p),
                        &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_key_replies,
                        batch);
            }
            world.gop.fence();
        }

        // 6. Ingest remote sources, complete `kHaloH2Ghost`, and prolong virtual halos per orbital.
        for (std::size_t b = 0; b < B; ++b) {
            auto* ws = ws_list[b];
            auto& virt_flags = virt_flags_all[b];
            const auto& needed_virt_capable = needed_vc_all[b];

            std::vector<WindowedSourceEntry<T>> remote_real_same_level;
            for (auto& rep : ex.take_key_replies(ws->slot_)) {
                const int rlvl = static_cast<int>(rep.key.level());
                if (rep.is_internal == 0) {
                    if (rlvl >= 0 && rlvl <= ws->n_max_tree_) {
                        ws->coarse_leaves_by_level_[static_cast<std::size_t>(rlvl)].emplace_back(
                            rep.key, rep.coeff);
                    }
                    WindowedSourceEntry<T> se;
                    se.key = rep.key;
                    se.is_internal = false;
                    se.s_coeff = rep.coeff.iscontiguous() ? rep.coeff : copy(rep.coeff);
                    windowed_halo_detail::compute_scaling_norms_readonly(
                        se.s_coeff, se.s000, se.s_osc_norm2, se.s_normf);
                    se.ns_normf = se.s_normf;
                    remote_real_same_level.push_back(std::move(se));
                } else {
                    WindowedSourceEntry<T> se;
                    se.key = rep.key;
                    se.is_internal = true;
                    se.s_coeff = windowed_halo_detail::extract_scaling_cube(rep.coeff, ws->k_);
                    windowed_halo_detail::compute_scaling_norms_readonly(
                        se.s_coeff, se.s000, se.s_osc_norm2, se.s_normf);
                    se.ns_coeff = rep.coeff;
                    se.d_normf = windowed_halo_detail::wavelet_octants_normf(rep.coeff, ws->k_);
                    se.ns_normf = std::sqrt(se.s_normf * se.s_normf + se.d_normf * se.d_normf);
                    remote_real_same_level.push_back(std::move(se));
                }
            }

            for (int j = 0; j <= ws->n_max_tree_; ++j) {
                auto& cl = ws->coarse_leaves_by_level_[static_cast<std::size_t>(j)];
                std::sort(cl.begin(), cl.end(), [](const auto& a, const auto& b) noexcept {
                    return windowed_halo_detail::windowed_key_less(a.first, b.first);
                });
                cl.erase(
                    std::unique(cl.begin(), cl.end(), [](const auto& a, const auto& b) noexcept {
                        return a.first == b.first;
                    }),
                    cl.end());
            }

            std::sort(remote_real_same_level.begin(),
                      remote_real_same_level.end(),
                      [](const auto& a, const auto& b) noexcept {
                          return windowed_halo_detail::windowed_key_less(a.key, b.key);
                      });
            remote_real_same_level.erase(
                std::unique(remote_real_same_level.begin(),
                            remote_real_same_level.end(),
                            [](const auto& a, const auto& b) noexcept { return a.key == b.key; }),
                remote_real_same_level.end());

            for (int j = 1; j <= ws->n_max_tree_; ++j) {
                auto& vf = virt_flags[static_cast<std::size_t>(j)];
                std::vector<Key<3>> h2_add;
                for (const Key<3>& skey : needed_virt_capable[static_cast<std::size_t>(j)]) {
                    if (ws->find_owned_real_index(j, skey) != SIZE_MAX) continue;
                    if (find_windowed_source(remote_real_same_level, skey) != nullptr) continue;
                    if (!windowed_halo_detail::contains_sorted_key_flag(vf, skey)) {
                        h2_add.push_back(skey);
                    }
                }
                if (!h2_add.empty()) {
                    for (const Key<3>& k_add : h2_add) {
                        vf.emplace_back(k_add, kHaloH2Ghost);
                    }
                    windowed_halo_detail::sort_and_merge_key_flags(vf);
                }
            }

            ws->populate_prolonged_virtual_and_sources(virt_flags, remote_real_same_level);
        }
    }

    /// Lockstep Stage 5a distributed restriction + nonstandard output tree assembly across `B`
    /// orbitals simultaneously, coalescing promotion count reduction, coarse-leaf / virtual-subtree
    /// transport, and remote output-node inserts into 2–3 global collectives total for the batch.
    static void assemble_output_nonstandard_lockstep(
        World& world,
        const std::vector<WindowedNSHaloWorkspace<T>*>& ws_list,
        const std::vector<Function<T, 3>*>& res_list,
        windowed_halo_detail::WindowedNSHaloExchange<T>& ex) {
        const std::size_t B = ws_list.size();
        if (B == 0) return;
        MADNESS_CHECK(world.size() > 1);
        MADNESS_CHECK(res_list.size() == B && ex.num_slots() >= B);

        const ProcessID me = world.rank();
        const int nproc = world.size();
        const TensorArgs full_targs(-1.0, TT_FULL);

        // 1. Coalesced global promotion check (1 `gop.sum` across all `B` orbitals).
        std::vector<long> promo_counts(B, 0L);
        for (std::size_t b = 0; b < B; ++b) {
            promo_counts[b] = ws_list[b]->count_local_promotions();
        }
        world.gop.sum(promo_counts.data(), promo_counts.size());

        // 2. Local unpromoted bottom-up restriction + coalesced coarse-delta / virtual-subtree exchange (1 fence).
        std::vector<std::vector<windowed_halo_detail::CoarseLeafDeltaMsg<T>>> delta_outbox(
            static_cast<std::size_t>(nproc));
        std::vector<std::vector<windowed_halo_detail::VirtualSubtreeItemMsg<T>>> subtree_outbox(
            static_cast<std::size_t>(nproc));

        bool any_slot_promoted = false;
        for (std::size_t b = 0; b < B; ++b) {
            auto* ws = ws_list[b];
            const std::uint16_t slot = ws->slot_;
            const long kk3 =
                static_cast<long>(ws->k_) * static_cast<long>(ws->k_) * static_cast<long>(ws->k_);
            if (promo_counts[b] == 0) {
                auto leaf_contribs = ws->compute_local_leaf_contribs_unpromoted();
                for (auto& item : leaf_contribs) {
                    const ProcessID owner = ws->pmap_->owner(item.key);
                    delta_outbox[static_cast<std::size_t>(owner)].push_back(
                        windowed_halo_detail::CoarseLeafDeltaMsg<T>{
                            item.key, item.min_source_key, std::move(item.ds), slot});
                }
            } else {
                any_slot_promoted = true;
                for (int j = 1; j <= ws->n_max_tree_; ++j) {
                    for (auto& vn : ws->levels_[static_cast<std::size_t>(j)].owned_virt) {
                        if (!vn.promoted && vn.ds_corr.size() != kk3) continue;
                        const ProcessID owner = ws->pmap_->owner(vn.coarse_leaf_key);
                        subtree_outbox[static_cast<std::size_t>(owner)].push_back(
                            windowed_halo_detail::VirtualSubtreeItemMsg<T>{
                                vn.coarse_leaf_key,
                                vn.key,
                                static_cast<std::uint8_t>(vn.promoted ? 1u : 0u),
                                vn.ds_corr,
                                vn.ns_virt_out,
                                slot});
                    }
                }
            }
        }

        for (int p = 0; p < nproc; ++p) {
            auto& db = delta_outbox[static_cast<std::size_t>(p)];
            if (!db.empty()) {
                if (p == me) {
                    ex.recv_coarse_deltas(db);
                } else {
                    ex.send(static_cast<ProcessID>(p),
                            &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_coarse_deltas,
                            db);
                }
            }
            auto& sb = subtree_outbox[static_cast<std::size_t>(p)];
            if (!sb.empty()) {
                if (p == me) {
                    ex.recv_virtual_subtree_items(sb);
                } else {
                    ex.send(
                        static_cast<ProcessID>(p),
                        &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_virtual_subtree_items,
                        sb);
                }
            }
        }
        world.gop.fence();

        // 3. Apply received deltas / promoted subtrees + emit owned real nodes; coalesce all remote output inserts!
        std::vector<std::vector<windowed_halo_detail::OutputNodeInsertMsg<T>>> insert_outbox(
            static_cast<std::size_t>(nproc));

        for (std::size_t b = 0; b < B; ++b) {
            auto* ws = ws_list[b];
            auto* res_impl = res_list[b]->get_impl().get();
            auto& res_coeffs = res_impl->get_coeffs();
            const long kk3 =
                static_cast<long>(ws->k_) * static_cast<long>(ws->k_) * static_cast<long>(ws->k_);

            if (promo_counts[b] == 0) {
                auto received = ex.take_coarse_deltas(ws->slot_);
                std::sort(received.begin(),
                          received.end(),
                          [](const auto& a, const auto& b_msg) noexcept {
                              if (a.leaf_key != b_msg.leaf_key) {
                                  return windowed_halo_detail::windowed_key_less(
                                      a.leaf_key, b_msg.leaf_key);
                              }
                              return windowed_halo_detail::windowed_key_less(
                                  a.sort_tag_key, b_msg.sort_tag_key);
                          });
                for (const auto& msg : received) {
                    const int blvl = static_cast<int>(msg.leaf_key.level());
                    const std::size_t ridx = ws->find_owned_real_index(blvl, msg.leaf_key);
                    MADNESS_CHECK(ridx != SIZE_MAX);
                    T* dst =
                        ws->levels_[static_cast<std::size_t>(blvl)].owned_real[ridx].s_out.ptr();
                    aligned_axpy(kk3, dst, msg.delta_s.ptr(), 1.0);
                }
            } else {
                auto local_items = ex.take_subtree_items(ws->slot_);
                ws->apply_promoted_subtree_items(
                    local_items, res_coeffs, full_targs, &insert_outbox);
            }

            for (int j = 1; j <= ws->n_max_tree_; ++j) {
                std::vector<WindowedOwnedVirtNode<T>>().swap(
                    ws->levels_[static_cast<std::size_t>(j)].owned_virt);
            }

            ws->emit_owned_real_nodes_nonstandard(res_coeffs, full_targs, &insert_outbox);
        }

        // 4. Single remote output-node insert fence (strictly skipped when no orbital had promotions!).
        if (any_slot_promoted) {
            for (int p = 0; p < nproc; ++p) {
                auto& batch = insert_outbox[static_cast<std::size_t>(p)];
                if (batch.empty()) continue;
                if (p == me) {
                    ex.recv_output_inserts(batch);
                } else {
                    ex.send(static_cast<ProcessID>(p),
                            &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_output_inserts,
                            batch);
                }
            }
            world.gop.fence();

            for (std::size_t b = 0; b < B; ++b) {
                if (promo_counts[b] == 0) continue;
                auto* ws = ws_list[b];
                auto& res_coeffs = res_list[b]->get_impl()->get_coeffs();
                for (auto& msg : ex.take_output_inserts(ws->slot_)) {
                    GenTensor<T> coeff(msg.coeff, full_targs);
                    res_coeffs.replace(msg.key, FunctionNode<T, 3>(coeff, msg.has_children != 0));
                }
            }
        }

        for (std::size_t b = 0; b < B; ++b) {
            res_list[b]->get_impl()->set_tree_state(nonstandard_after_apply);
        }
    }

    /// Stage 5a + Stage 5b assembly:
    ///   1. Restrict all non-promoted virtual boundary corrections `\Delta\hat{s}^j_{V,corr}`
    ///      bottom-up (`j = n_max_tree -> 1`) onto their real coarse leaves `B in L`
    ///      (or split subtrees when dynamic output refinement promotes virtual/leaf nodes).
    ///   2. Populate `result` in `TreeState::nonstandard_after_apply` so a subsequent
    ///      `result.reconstruct(fence)` performs top-down NS synthesis cleanly.
    void assemble_output_nonstandard_tree(Function<T, 3>& result) {
        if (world_.size() > 1) {
            MADNESS_CHECK(exchange_ptr_ != nullptr);
            assemble_output_nonstandard_lockstep(world_, {this}, {&result}, *exchange_ptr_);
            return;
        }

        auto* res_impl = result.get_impl().get();
        auto& res_coeffs = res_impl->get_coeffs();
        const TensorArgs full_targs(-1.0, TT_FULL);

        const long local_any_promoted = count_local_promotions();
        if (local_any_promoted == 0) {
            restrict_all_virtual_unpromoted();
        } else {
            const long kk3 = static_cast<long>(k_) * static_cast<long>(k_) * static_cast<long>(k_);
            std::vector<windowed_halo_detail::VirtualSubtreeItemMsg<T>> local_items;
            for (int j = 1; j <= n_max_tree_; ++j) {
                for (auto& vn : levels_[static_cast<std::size_t>(j)].owned_virt) {
                    if (!vn.promoted && vn.ds_corr.size() != kk3) continue;
                    local_items.push_back(windowed_halo_detail::VirtualSubtreeItemMsg<T>{
                        vn.coarse_leaf_key,
                        vn.key,
                        static_cast<std::uint8_t>(vn.promoted ? 1u : 0u),
                        vn.ds_corr,
                        vn.ns_virt_out,
                        slot_});
                }
            }
            apply_promoted_subtree_items(local_items, res_coeffs, full_targs, nullptr);
        }
        for (int j = 1; j <= n_max_tree_; ++j) {
            std::vector<WindowedOwnedVirtNode<T>>().swap(
                levels_[static_cast<std::size_t>(j)].owned_virt);
        }

        emit_owned_real_nodes_nonstandard(res_coeffs, full_targs, nullptr);
        res_impl->set_tree_state(nonstandard_after_apply);
    }

private:
    using KeyFlagVec = std::vector<std::pair<Key<3>, std::uint8_t>>;

    struct RestrictAccum {
        Key<3> key;
        Key<3> coarse_leaf;
        Key<3> min_source_key;
        Tensor<T> ds;
    };

    long count_local_promotions() const noexcept {
        long count = 0;
        for (int j = 0; j <= n_max_tree_; ++j) {
            for (const auto& rn : levels_[static_cast<std::size_t>(j)].owned_real) {
                if (rn.promoted_leaf && j < max_refine_level_) ++count;
            }
            if (j >= 1) {
                for (const auto& vn : levels_[static_cast<std::size_t>(j)].owned_virt) {
                    if (vn.promoted && j < max_refine_level_) ++count;
                }
            }
        }
        return count;
    }

    void emit_owned_real_nodes_nonstandard(
        typename FunctionImpl<T, 3>::dcT& res_coeffs,
        const TensorArgs& full_targs,
        std::vector<std::vector<windowed_halo_detail::OutputNodeInsertMsg<T>>>* remote_outbox) {
        const long kk = static_cast<long>(k_);
        const long twok = 2L * kk;
        const Slice s0(0, kk - 1);
        const ProcessID me = world_.rank();

        for (int j = 0; j <= n_max_tree_; ++j) {
            auto& ld = levels_[static_cast<std::size_t>(j)];
            for (auto& rn : ld.owned_real) {
                if (!split_coarse_leaves_emitted_.empty() &&
                    std::binary_search(split_coarse_leaves_emitted_.begin(),
                                       split_coarse_leaves_emitted_.end(),
                                       rn.key,
                                       windowed_halo_detail::windowed_key_less)) {
                    rn.ns_out = Tensor<T>();
                    rn.s_out = Tensor<T>();
                    continue;
                }

                if (rn.is_internal) {
                    Tensor<T> full_ns = rn.ns_out;
                    if (full_ns.size() != twok * twok * twok) {
                        full_ns = Tensor<T>(twok, twok, twok);
                    }
                    full_ns(s0, s0, s0) = rn.s_out;
                    rn.ns_out = Tensor<T>();
                    rn.s_out = Tensor<T>();
                    GenTensor<T> coeff(full_ns, full_targs);
                    res_coeffs.replace(rn.key, FunctionNode<T, 3>(coeff, true));
                } else if (rn.promoted_leaf && j < max_refine_level_) {
                    Tensor<T> full_ns = rn.ns_out;
                    if (full_ns.size() != twok * twok * twok) {
                        full_ns = Tensor<T>(twok, twok, twok);
                    }
                    full_ns(s0, s0, s0) = rn.s_out;
                    rn.ns_out = Tensor<T>();
                    rn.s_out = Tensor<T>();
                    GenTensor<T> coeff(full_ns, full_targs);
                    res_coeffs.replace(rn.key, FunctionNode<T, 3>(coeff, true));

                    for (KeyChildIterator<3> kit(rn.key); kit; ++kit) {
                        const Key<3>& child = kit.key();
                        Tensor<T> zero_s(kk, kk, kk);
                        const ProcessID owner = (world_.size() == 1) ? me : pmap_->owner(child);
                        if (owner == me || remote_outbox == nullptr) {
                            GenTensor<T> c_coeff(zero_s, full_targs);
                            res_coeffs.replace(child, FunctionNode<T, 3>(c_coeff, false));
                        } else {
                            (*remote_outbox)[static_cast<std::size_t>(owner)].push_back(
                                windowed_halo_detail::OutputNodeInsertMsg<T>{
                                    child, 0u, std::move(zero_s), slot_});
                        }
                    }
                } else {
                    GenTensor<T> coeff(rn.s_out, full_targs);
                    rn.s_out = Tensor<T>();
                    res_coeffs.replace(rn.key, FunctionNode<T, 3>(coeff, false));
                }
            }
        }
    }

    void init_octant_filters() {
        for (int oct = 0; oct < 8; ++oct) {
            h_oct_[oct][0] = ((oct & 4) == 0) ? cdata_.h0 : cdata_.h1;
            h_oct_[oct][1] = ((oct & 2) == 0) ? cdata_.h0 : cdata_.h1;
            h_oct_[oct][2] = ((oct & 1) == 0) ? cdata_.h0 : cdata_.h1;

            hT_oct_[oct][0] = ((oct & 4) == 0) ? cdata_.h0T : cdata_.h1T;
            hT_oct_[oct][1] = ((oct & 2) == 0) ? cdata_.h0T : cdata_.h1T;
            hT_oct_[oct][2] = ((oct & 1) == 0) ? cdata_.h0T : cdata_.h1T;
        }
    }

    std::size_t find_owned_real_packed(int j, std::uint64_t code) const noexcept {
        const auto& pvec = levels_[static_cast<std::size_t>(j)].real_packed_keys;
        const auto it = std::lower_bound(pvec.begin(), pvec.end(), code);
        if (it != pvec.end() && *it == code) {
            return static_cast<std::size_t>(it - pvec.begin());
        }
        return SIZE_MAX;
    }

    std::size_t find_owned_real_index(int j, const Key<3>& key) const noexcept {
        if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
            return find_owned_real_packed(
                j, windowed_halo_detail::pack_spatial_key_same_level(key));
        }
        const auto& vec = levels_[static_cast<std::size_t>(j)].owned_real;
        const auto it = std::lower_bound(
            vec.begin(), vec.end(), key,
            [](const WindowedOwnedRealNode<T>& rn, const Key<3>& k) noexcept {
                return windowed_halo_detail::windowed_key_less(rn.key, k);
            });
        if (it != vec.end() && it->key == key) {
            return static_cast<std::size_t>(it - vec.begin());
        }
        return SIZE_MAX;
    }

    const Tensor<T>* find_coarse_leaf(int j, const Key<3>& key) const noexcept {
        if (world_.size() == 1) {
            const std::size_t idx = find_owned_real_index(j, key);
            if (idx != SIZE_MAX) {
                const auto& rn = levels_[static_cast<std::size_t>(j)].owned_real[idx];
                if (!rn.is_internal) {
                    return &rn.s_in;
                }
            }
            return nullptr;
        }
        const auto& vec = coarse_leaves_by_level_[static_cast<std::size_t>(j)];
        const auto it = std::lower_bound(
            vec.begin(), vec.end(), key,
            [](const std::pair<Key<3>, Tensor<T>>& kv, const Key<3>& k) noexcept {
                return windowed_halo_detail::windowed_key_less(kv.first, k);
            });
        if (it != vec.end() && it->first == key) {
            return &(it->second);
        }
        return nullptr;
    }

    void ingest_local_real_nodes_local_only(const Function<T, 3>& f_ns) {
        const long kk = static_cast<long>(k_);
        const auto& coeffs = f_ns.get_impl()->get_coeffs();

        std::vector<std::pair<Key<3>, const FunctionNode<T, 3>*>> raw_nodes;
        raw_nodes.reserve(coeffs.size());
        for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
            const int j = static_cast<int>(it->first.level());
            MADNESS_CHECK(j >= 0 && j <= n_max_tree_);
            raw_nodes.emplace_back(it->first, &(it->second));
        }

        std::vector<WindowedOwnedRealNode<T>> built_nodes(raw_nodes.size());
        windowed_parallel_detail::parallel_for_chunks(
            world_,
            raw_nodes.size(),
            windowed_parallel_detail::kDefaultParallelGrain,
            [&](std::size_t b, std::size_t e_hi, std::size_t) {
                for (std::size_t idx = b; idx < e_hi; ++idx) {
                    const Key<3>& key = raw_nodes[idx].first;
                    const FunctionNode<T, 3>& node = *raw_nodes[idx].second;

                    WindowedOwnedRealNode<T>& rn = built_nodes[idx];
                    rn.key = key;
                    rn.is_internal = node.has_children();
                    rn.s_out = Tensor<T>(kk, kk, kk);

                    const auto& gen_c = node.coeff();
                    const bool is_contig_full =
                        (gen_c.tensor_type() == TT_FULL && gen_c.get_tensor().iscontiguous());
                    const Tensor<T> full_c =
                        is_contig_full ? gen_c.get_tensor() : gen_c.full_tensor_copy();
                    if (rn.is_internal) {
                        rn.ns_in = full_c;
                        rn.s_in = windowed_halo_detail::extract_scaling_cube(full_c, k_);
                        windowed_halo_detail::compute_scaling_norms_readonly(
                            rn.s_in, rn.s000, rn.s_osc_norm2, rn.s_normf);
                        rn.d_normf = windowed_halo_detail::wavelet_octants_normf(full_c, k_);
                        rn.ns_normf = std::sqrt(rn.s_normf * rn.s_normf + rn.d_normf * rn.d_normf);
                    } else {
                        if (full_c.dim(0) == kk && full_c.dim(1) == kk && full_c.dim(2) == kk &&
                            full_c.iscontiguous()) {
                            rn.s_in = full_c;
                        } else {
                            rn.s_in = windowed_halo_detail::extract_scaling_cube(full_c, k_);
                        }
                        windowed_halo_detail::compute_scaling_norms_readonly(
                            rn.s_in, rn.s000, rn.s_osc_norm2, rn.s_normf);
                        rn.ns_normf = rn.s_normf;
                    }
                }
            });

        std::vector<std::size_t> lvl_counts(static_cast<std::size_t>(n_max_tree_ + 1), 0u);
        for (const auto& kv : raw_nodes) {
            ++lvl_counts[static_cast<std::size_t>(kv.first.level())];
        }
        for (int j = 0; j <= n_max_tree_; ++j) {
            levels_[static_cast<std::size_t>(j)].owned_real.reserve(
                lvl_counts[static_cast<std::size_t>(j)]);
        }

        const double k3_fac = static_cast<double>(k_ * k_ * k_);
        for (auto& rn : built_nodes) {
            const std::size_t sj = static_cast<std::size_t>(rn.key.level());
            auto& st = global_stats_[sj];
            st.max_s_real = std::max(st.max_s_real, rn.s_normf);
            if (rn.is_internal) {
                st.has_internal = true;
                st.max_d_int = std::max(st.max_d_int, rn.d_normf);
                st.max_ns_int = std::max(st.max_ns_int, rn.ns_normf);
            } else {
                st.max_s_leaf = std::max(st.max_s_leaf, rn.s_normf);
                const double leaf_linf =
                    std::sqrt(windowed_halo_detail::scalar_abs2(rn.s000)) +
                    k3_fac * std::sqrt(std::max(0.0, rn.s_osc_norm2));
                st.max_s_leaf_linf = std::max(st.max_s_leaf_linf, leaf_linf);
            }
            levels_[sj].owned_real.push_back(std::move(rn));
        }
    }

    double compute_coarse_leaf_above_bound(int j) const noexcept {
        double bound = 0.0;
        for (int b = 0; b < j; ++b) {
            const auto& sb = global_stats_[static_cast<std::size_t>(b)];
            if (!(sb.max_s_leaf > 0.0)) continue;
            const int shift3 = 3 * (j - b);
            const double vol_sqrt = std::ldexp(1.0, -(shift3 / 2)) *
                                    ((shift3 & 1) ? 0.7071067811865476 : 1.0);
            const double v_bound = std::min(sb.max_s_leaf, vol_sqrt * sb.max_s_leaf_linf);
            bound = std::max(bound, v_bound);
        }
        return bound;
    }

    void finalize_ingested_global_stats(bool populate_coarse_leaf_cache) {
        for (int j = 0; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            auto& st = global_stats_[sj];
            const double coarse_leaf_max = compute_coarse_leaf_above_bound(j);
            st.max_coarse_leaf_above = coarse_leaf_max;
            st.max_s_source = std::max(st.max_s_real, coarse_leaf_max);

            auto& ld = levels_[sj];
            ld.global_has_internal = st.has_internal;
            std::sort(ld.owned_real.begin(), ld.owned_real.end(), [](const auto& a, const auto& b) noexcept {
                return windowed_halo_detail::windowed_key_less(a.key, b.key);
            });
            if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                ld.real_packed_keys.resize(ld.owned_real.size());
                for (std::size_t ri = 0; ri < ld.owned_real.size(); ++ri) {
                    ld.real_packed_keys[ri] =
                        windowed_halo_detail::pack_spatial_key_same_level(ld.owned_real[ri].key);
                }
            }

            if (populate_coarse_leaf_cache) {
                auto& cl_vec = coarse_leaves_by_level_[sj];
                for (const auto& rn : ld.owned_real) {
                    if (!rn.is_internal) {
                        cl_vec.emplace_back(rn.key, rn.s_in);
                    }
                }
            }
        }
    }

    void ingest_local_real_nodes(const Function<T, 3>& f_ns) {
        ingest_local_real_nodes_local_only(f_ns);

        if (world_.size() > 1) {
            const std::size_t n_lvls = static_cast<std::size_t>(n_max_tree_ + 1);
            std::vector<double> gop_buf(6 * n_lvls, 0.0);
            for (std::size_t j = 0; j < n_lvls; ++j) {
                const auto& st = global_stats_[j];
                gop_buf[6 * j + 0] = st.has_internal ? 1.0 : 0.0;
                gop_buf[6 * j + 1] = st.max_s_real;
                gop_buf[6 * j + 2] = st.max_s_leaf;
                gop_buf[6 * j + 3] = st.max_s_leaf_linf;
                gop_buf[6 * j + 4] = st.max_d_int;
                gop_buf[6 * j + 5] = st.max_ns_int;
            }
            world_.gop.max(gop_buf.data(), gop_buf.size());
            for (std::size_t j = 0; j < n_lvls; ++j) {
                auto& st = global_stats_[j];
                st.has_internal = (gop_buf[6 * j + 0] > 0.5);
                st.max_s_real = gop_buf[6 * j + 1];
                st.max_s_leaf = gop_buf[6 * j + 2];
                st.max_s_leaf_linf = gop_buf[6 * j + 3];
                st.max_d_int = gop_buf[6 * j + 4];
                st.max_ns_int = gop_buf[6 * j + 5];
            }
        }

        finalize_ingested_global_stats(world_.size() > 1);

        if (world_.size() == 1 && dynamic_refine_ && ungraded_leaf_refine_only_) {
            populate_serial_internal_child_masks();
        }
    }

    void populate_serial_internal_child_masks() {
        for (int j = 1; j <= n_max_tree_; ++j) {
            for (const auto& ch : levels_[static_cast<std::size_t>(j)].owned_real) {
                if (!ch.is_internal) continue;
                const Key<3> par = ch.key.parent();
                const std::size_t pidx = find_owned_real_index(j - 1, par);
                if (pidx != SIZE_MAX) {
                    const int oct = windowed_halo_detail::child_octant_index(ch.key);
                    auto& prn = levels_[static_cast<std::size_t>(j - 1)].owned_real[pidx];
                    prn.internal_child_mask = static_cast<std::uint8_t>(
                        prn.internal_child_mask | static_cast<std::uint8_t>(1u << oct));
                }
            }
        }
    }

    void finalize_serial_stats_and_child_masks() {
        const double k3_fac = static_cast<double>(k_ * k_ * k_);
        for (int j = 0; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            auto& ld = levels_[sj];
            auto& st = global_stats_[sj];
            for (const auto& rn : ld.owned_real) {
                st.max_s_real = std::max(st.max_s_real, rn.s_normf);
                if (rn.is_internal) {
                    st.has_internal = true;
                    st.max_d_int = std::max(st.max_d_int, rn.d_normf);
                    st.max_ns_int = std::max(st.max_ns_int, rn.ns_normf);
                } else {
                    st.max_s_leaf = std::max(st.max_s_leaf, rn.s_normf);
                    const double leaf_linf =
                        std::sqrt(windowed_halo_detail::scalar_abs2(rn.s000)) +
                        k3_fac * std::sqrt(std::max(0.0, rn.s_osc_norm2));
                    st.max_s_leaf_linf = std::max(st.max_s_leaf_linf, leaf_linf);
                }
            }
            const double coarse_leaf_max = compute_coarse_leaf_above_bound(j);
            st.max_coarse_leaf_above = coarse_leaf_max;
            st.max_s_source = std::max(st.max_s_real, coarse_leaf_max);
            ld.global_has_internal = st.has_internal;
        }
        if (dynamic_refine_ && ungraded_leaf_refine_only_) {
            populate_serial_internal_child_masks();
        }
    }

    /// Zero-mutation read-only serial ingestion directly from `TreeState::reconstructed`:
    ///   - Reads leaf `k^3` tensors in O(1) shallow handles (zero `f.make_nonstandard()`, zero `f.reconstruct()`).
    ///   - Restricts internal nodes bottom-up strictly down to `serial_range_.min_s_int_level`, skipping all
    ///     coarser levels `j < min_s_int_level` completely (and skipping ALL internal nodes when
    ///     `serial_range_.max_s_int_level < 0`, e.g., when `j_0 == n_max_tree` and `raw_gamma_terms` is empty).
    ///   - Computes `(2k)^3` `ns_in` + `wav_in` strictly on `[min_gamma_level, max_gamma_level]`, using
    ///     lightweight `k^3` `restrict_step_accumulate_3d` (`24 k^4` flops, zero `(2k)^3` allocs) elsewhere.
    void ingest_local_real_nodes_serial_reconstructed(const Function<T, 3>& f_rec) {
        const long kk = static_cast<long>(k_);
        const long kk3 = kk * kk * kk;
        const long twok = 2L * kk;
        const long twok3 = twok * twok * twok;
        const long k_dims[3] = {kk, kk, kk};
        const long twok_dims[3] = {twok, twok, twok};
        const Slice s0(0, kk - 1);
        const auto& coeffs = f_rec.get_impl()->get_coeffs();

        std::vector<std::pair<Key<3>, const FunctionNode<T, 3>*>> raw_nodes;
        raw_nodes.reserve(coeffs.size());
        std::vector<std::size_t> lvl_counts(static_cast<std::size_t>(n_max_tree_ + 1), 0u);
        for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
            const int j = static_cast<int>(it->first.level());
            MADNESS_CHECK(j >= 0 && j <= n_max_tree_);
            raw_nodes.emplace_back(it->first, &(it->second));
            ++lvl_counts[static_cast<std::size_t>(j)];
        }

        if (raw_nodes.empty()) {
            WindowedOwnedRealNode<T> root_rn;
            root_rn.key = cdata_.key0;
            root_rn.is_internal = false;
            root_rn.s_in = Tensor<T>(kk, kk, kk);
            root_rn.s_out = Tensor<T>(kk, kk, kk);
            levels_[0].owned_real.push_back(std::move(root_rn));
            levels_[0].real_packed_keys = {
                windowed_halo_detail::pack_spatial_key_same_level(cdata_.key0)};
            finalize_serial_stats_and_child_masks();
            return;
        }

        std::vector<WindowedOwnedRealNode<T>> built_nodes(raw_nodes.size());
        windowed_parallel_detail::parallel_for_chunks(
            world_,
            raw_nodes.size(),
            windowed_parallel_detail::kDefaultParallelGrain,
            [&](std::size_t b, std::size_t e_hi, std::size_t) {
                for (std::size_t idx = b; idx < e_hi; ++idx) {
                    const Key<3>& key = raw_nodes[idx].first;
                    const FunctionNode<T, 3>& node = *raw_nodes[idx].second;

                    WindowedOwnedRealNode<T>& rn = built_nodes[idx];
                    rn.key = key;
                    rn.is_internal = node.has_children();
                    rn.s_out = Tensor<T>(kk, kk, kk);

                    if (!rn.is_internal) {
                        if (node.has_coeff()) {
                            const auto& gen_c = node.coeff();
                            const bool is_contig_full =
                                (gen_c.tensor_type() == TT_FULL && gen_c.get_tensor().iscontiguous());
                            const Tensor<T> full_c =
                                is_contig_full ? gen_c.get_tensor() : gen_c.full_tensor_copy();
                            if (full_c.dim(0) == kk && full_c.dim(1) == kk && full_c.dim(2) == kk &&
                                full_c.iscontiguous()) {
                                rn.s_in = full_c;
                            } else {
                                rn.s_in = windowed_halo_detail::extract_scaling_cube(full_c, k_);
                            }
                        } else {
                            rn.s_in = Tensor<T>(kk, kk, kk);
                        }
                        windowed_halo_detail::compute_scaling_norms_readonly(
                            rn.s_in, rn.s000, rn.s_osc_norm2, rn.s_normf);
                        rn.ns_normf = rn.s_normf;
                    }
                }
            });

        for (int j = 0; j <= n_max_tree_; ++j) {
            levels_[static_cast<std::size_t>(j)].owned_real.reserve(
                lvl_counts[static_cast<std::size_t>(j)]);
        }
        for (auto& rn : built_nodes) {
            levels_[static_cast<std::size_t>(rn.key.level())].owned_real.push_back(std::move(rn));
        }

        for (int j = 0; j <= n_max_tree_; ++j) {
            auto& ld = levels_[static_cast<std::size_t>(j)];
            std::sort(ld.owned_real.begin(), ld.owned_real.end(), [](const auto& a, const auto& b) noexcept {
                return windowed_halo_detail::windowed_key_less(a.key, b.key);
            });
            if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                ld.real_packed_keys.resize(ld.owned_real.size());
                for (std::size_t ri = 0; ri < ld.owned_real.size(); ++ri) {
                    ld.real_packed_keys[ri] =
                        windowed_halo_detail::pack_spatial_key_same_level(ld.owned_real[ri].key);
                }
            }
        }

        if (n_max_tree_ >= 1 &&
            serial_range_.max_s_int_level >= 0 &&
            serial_range_.min_s_int_level <= serial_range_.max_s_int_level) {
            const int j_stop = std::max(0, serial_range_.min_s_int_level);
            std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_box_2k;
            std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_work_2k;
            std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_w1_k3;
            std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_w2_k3;
            std::vector<std::size_t> int_indices;

            for (int j = n_max_tree_ - 1; j >= j_stop; --j) {
                auto& ld = levels_[static_cast<std::size_t>(j)];
                const auto& child_ld = levels_[static_cast<std::size_t>(j + 1)];
                int_indices.clear();
                for (std::size_t ri = 0; ri < ld.owned_real.size(); ++ri) {
                    if (ld.owned_real[ri].is_internal) {
                        int_indices.push_back(ri);
                    }
                }
                if (int_indices.empty()) continue;

                const bool need_gamma_j =
                    (j >= serial_range_.min_gamma_level && j <= serial_range_.max_gamma_level);

                windowed_parallel_detail::parallel_for_chunks(
                    world_,
                    int_indices.size(),
                    2,
                    [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
                        if (need_gamma_j) {
                            if (tls_box_2k[cid].size() != twok3) {
                                tls_box_2k[cid] = Tensor<T>(3L, twok_dims, false);
                                tls_work_2k[cid] = Tensor<T>(3L, twok_dims, false);
                            }
                            Tensor<T>& box_2k = tls_box_2k[cid];
                            Tensor<T>& work_2k = tls_work_2k[cid];

                            for (std::size_t ii = b; ii < e_hi; ++ii) {
                                auto& rn = ld.owned_real[int_indices[ii]];
                                for (KeyChildIterator<3> kit(rn.key); kit; ++kit) {
                                    const Key<3>& ch = kit.key();
                                    const std::size_t cidx = find_owned_real_index(j + 1, ch);
                                    MADNESS_CHECK(cidx != SIZE_MAX);
                                    const auto& cl = ch.translation();
                                    box_2k(cdata_.s[cl[0] & 1],
                                           cdata_.s[cl[1] & 1],
                                           cdata_.s[cl[2] & 1]) = child_ld.owned_real[cidx].s_in;
                                }
                                rn.ns_in = Tensor<T>(3L, twok_dims, false);
                                fast_transform(box_2k, cdata_.hgT, rn.ns_in, work_2k);
                                rn.s_in = copy(rn.ns_in(s0, s0, s0));
                                windowed_halo_detail::compute_scaling_norms_readonly(
                                    rn.s_in, rn.s000, rn.s_osc_norm2, rn.s_normf);
                                rn.d_normf =
                                    windowed_halo_detail::wavelet_octants_normf(rn.ns_in, k_);
                                rn.ns_normf =
                                    std::sqrt(rn.s_normf * rn.s_normf + rn.d_normf * rn.d_normf);
                            }
                        } else {
                            if (tls_w1_k3[cid].size() != kk3) {
                                tls_w1_k3[cid] = Tensor<T>(3L, k_dims, false);
                                tls_w2_k3[cid] = Tensor<T>(3L, k_dims, false);
                            }
                            Tensor<T>& w1_k3 = tls_w1_k3[cid];
                            Tensor<T>& w2_k3 = tls_w2_k3[cid];

                            for (std::size_t ii = b; ii < e_hi; ++ii) {
                                auto& rn = ld.owned_real[int_indices[ii]];
                                rn.s_in = Tensor<T>(kk, kk, kk);
                                for (KeyChildIterator<3> kit(rn.key); kit; ++kit) {
                                    const Key<3>& ch = kit.key();
                                    const std::size_t cidx = find_owned_real_index(j + 1, ch);
                                    MADNESS_CHECK(cidx != SIZE_MAX);
                                    windowed_halo_detail::restrict_step_accumulate_3d(
                                        kk3,
                                        hT_oct_,
                                        child_ld.owned_real[cidx].s_in,
                                        ch,
                                        rn.s_in,
                                        w1_k3,
                                        w2_k3);
                                }
                                windowed_halo_detail::compute_scaling_norms_readonly(
                                    rn.s_in, rn.s000, rn.s_osc_norm2, rn.s_normf);
                                rn.ns_normf = rn.s_normf;
                            }
                        }
                    });
            }
        }

        finalize_serial_stats_and_child_masks();
    }

    /// Zero-mutation read-only serial ingestion directly from `TreeState::compressed`:
    ///   - Reads root `(s^0, d^0)` and internal `(0, d^j)` tensors read-only and synthesizes `s_in`
    ///     top-down (`0 .. n_max_tree`) in flat parallel arrays via `fast_transform(..., cdata_.hg)`
    ///     with zero `WorldContainer` mutations (`f` stays in `TreeState::compressed` untouched!).
    ///   - Matches `FunctionImpl::reconstruct_op` 1:1 on sparse compressed trees produced by
    ///     `vtransform` / `gaxpy_inplace` (where leaf siblings or empty intermediate ancestors may
    ///     be absent from `coeffs`), materializing all 8 octree children of every internal node.
    void ingest_local_real_nodes_serial_compressed(const Function<T, 3>& f_comp) {
        const long kk = static_cast<long>(k_);
        const long kk3 = kk * kk * kk;
        const long twok = 2L * kk;
        const long twok3 = twok * twok * twok;
        const long k_dims[3] = {kk, kk, kk};
        const long twok_dims[3] = {twok, twok, twok};
        const Slice s0(0, kk - 1);
        const auto& coeffs = f_comp.get_impl()->get_coeffs();

        struct RawCompEntry {
            Key<3> key;
            const FunctionNode<T, 3>* node = nullptr;
            double sub_d_norm2 = 0.0;
            bool keep_internal = false;
        };
        std::vector<std::vector<RawCompEntry>> raw_by_level(
            static_cast<std::size_t>(n_max_tree_ + 1));
        for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
            const int j = static_cast<int>(it->first.level());
            if (j >= 0 && j <= n_max_tree_) {
                raw_by_level[static_cast<std::size_t>(j)].push_back(
                    RawCompEntry{it->first, &(it->second), 0.0, false});
            }
        }
        for (int j = 0; j <= n_max_tree_; ++j) {
            auto& vec = raw_by_level[static_cast<std::size_t>(j)];
            std::sort(vec.begin(), vec.end(), [](const RawCompEntry& a, const RawCompEntry& b) noexcept {
                return windowed_halo_detail::windowed_key_less(a.key, b.key);
            });
        }

        auto extract_node_tensor = [](const FunctionNode<T, 3>* fn) -> Tensor<T> {
            if (fn == nullptr || !fn->has_coeff()) return Tensor<T>();
            const auto& gen_c = fn->coeff();
            if (!gen_c.has_data()) return Tensor<T>();
            const bool is_contig_full =
                (gen_c.tensor_type() == TT_FULL && gen_c.get_tensor().iscontiguous());
            return is_contig_full ? gen_c.get_tensor() : gen_c.full_tensor_copy();
        };

        auto raw_node_is_internal = [&](int j, const FunctionNode<T, 3>* fn) -> bool {
            if (fn == nullptr || j >= n_max_tree_) return false;
            if (fn->has_children()) return true;
            if (fn->has_coeff()) {
                const auto& gen_c = fn->coeff();
                if (gen_c.has_data() && gen_c.ndim() == 3 && gen_c.dim(0) == twok) {
                    if (j > 0) return true;
                    const Tensor<T> full_c = extract_node_tensor(fn);
                    return windowed_halo_detail::wavelet_octants_normf(full_c, k_) > 0.0;
                }
            }
            return false;
        };

        const bool use_subtree_prune =
            (serial_range_.subtree_d_tol != nullptr &&
             serial_range_.subtree_d_tol->size() >= static_cast<std::size_t>(n_max_tree_ + 1));

        if (use_subtree_prune && n_max_tree_ >= 1) {
            // Ensure every present node at level j >= 1 has its parent key in `raw_by_level[j - 1]`
            // (handles sparse compressed trees from `vtransform` with omitted intermediate nodes).
            std::vector<Key<3>> missing_pars;
            for (int j = n_max_tree_; j >= 1; --j) {
                const auto& cur_vec = raw_by_level[static_cast<std::size_t>(j)];
                if (cur_vec.empty()) continue;
                auto& par_vec = raw_by_level[static_cast<std::size_t>(j - 1)];
                missing_pars.clear();
                for (const auto& ce : cur_vec) {
                    if (ce.node != nullptr && !raw_node_is_internal(j, ce.node)) continue;
                    const Key<3> par = ce.key.parent();
                    const auto pit = std::lower_bound(
                        par_vec.begin(), par_vec.end(), par,
                        [](const RawCompEntry& e, const Key<3>& k) noexcept {
                            return windowed_halo_detail::windowed_key_less(e.key, k);
                        });
                    if (pit == par_vec.end() || pit->key != par) {
                        missing_pars.push_back(par);
                    }
                }
                if (!missing_pars.empty()) {
                    windowed_halo_detail::sort_unique_keys(missing_pars);
                    std::vector<RawCompEntry> merged;
                    merged.reserve(par_vec.size() + missing_pars.size());
                    std::size_t ia = 0, ib = 0;
                    while (ia < par_vec.size() && ib < missing_pars.size()) {
                        if (windowed_halo_detail::windowed_key_less(par_vec[ia].key, missing_pars[ib])) {
                            merged.push_back(par_vec[ia++]);
                        } else {
                            merged.push_back(RawCompEntry{missing_pars[ib++], nullptr, 0.0, false});
                        }
                    }
                    while (ia < par_vec.size()) merged.push_back(par_vec[ia++]);
                    while (ib < missing_pars.size()) {
                        merged.push_back(RawCompEntry{missing_pars[ib++], nullptr, 0.0, false});
                    }
                    par_vec.swap(merged);
                }
            }

            // Bottom-up exact wavelet-detail L2 norm accumulation & operator-weighted pruning.
            for (int j = n_max_tree_ - 1; j >= 0; --j) {
                auto& vec = raw_by_level[static_cast<std::size_t>(j)];
                if (vec.empty()) continue;
                windowed_parallel_detail::parallel_for_chunks(
                    world_,
                    vec.size(),
                    16,
                    [&](std::size_t b, std::size_t e_hi, std::size_t) {
                        for (std::size_t idx = b; idx < e_hi; ++idx) {
                            auto& e = vec[idx];
                            if (e.node != nullptr && e.node->has_coeff()) {
                                const auto& gen_c = e.node->coeff();
                                if (gen_c.has_data() && gen_c.ndim() == 3 && gen_c.dim(0) == twok) {
                                    const Tensor<T> full_c = extract_node_tensor(e.node);
                                    const double dn =
                                        windowed_halo_detail::wavelet_octants_normf(full_c, k_);
                                    e.sub_d_norm2 += dn * dn;
                                }
                            }
                            const bool can_be_int =
                                (j < n_max_tree_) &&
                                ((e.node == nullptr) || raw_node_is_internal(j, e.node));
                            if (j <= 1) {
                                // Match `FunctionImpl::truncate_op` (`key.level() > 1`): never prune `j <= 1`.
                                e.keep_internal = e.keep_internal || can_be_int;
                            } else if (!e.keep_internal && can_be_int) {
                                const double dtol =
                                    (*serial_range_.subtree_d_tol)[static_cast<std::size_t>(j)];
                                if (e.sub_d_norm2 > dtol * dtol) {
                                    e.keep_internal = true;
                                }
                            }
                        }
                    });

                if (j >= 1) {
                    auto& par_vec = raw_by_level[static_cast<std::size_t>(j - 1)];
                    for (const auto& ce : vec) {
                        if (!ce.keep_internal && !(ce.sub_d_norm2 > 0.0)) continue;
                        const Key<3> par = ce.key.parent();
                        auto pit = std::lower_bound(
                            par_vec.begin(), par_vec.end(), par,
                            [](const RawCompEntry& e, const Key<3>& k) noexcept {
                                return windowed_halo_detail::windowed_key_less(e.key, k);
                            });
                        if (pit != par_vec.end() && pit->key == par) {
                            pit->sub_d_norm2 += ce.sub_d_norm2;
                            if (ce.keep_internal) {
                                pit->keep_internal = true;
                            }
                        }
                    }
                }
            }
        }

        auto find_raw_entry = [&](int j, const Key<3>& key) noexcept -> const RawCompEntry* {
            if (j < 0 || j > n_max_tree_) return nullptr;
            const auto& vec = raw_by_level[static_cast<std::size_t>(j)];
            const auto it = std::lower_bound(
                vec.begin(),
                vec.end(),
                key,
                [](const RawCompEntry& e, const Key<3>& k) noexcept {
                    return windowed_halo_detail::windowed_key_less(e.key, k);
                });
            if (it != vec.end() && it->key == key) return &(*it);
            return nullptr;
        };

        auto find_raw = [&](int j, const Key<3>& key) noexcept -> const FunctionNode<T, 3>* {
            const RawCompEntry* re = find_raw_entry(j, key);
            return (re != nullptr) ? re->node : nullptr;
        };

        auto is_comp_internal = [&](int j,
                                    const Key<3>& key,
                                    const FunctionNode<T, 3>* fn,
                                    const Tensor<T>& full_c) -> bool {
            if (j >= n_max_tree_) return false;
            if (use_subtree_prune && n_max_tree_ >= 1) {
                const RawCompEntry* re = find_raw_entry(j, key);
                return (re != nullptr) && re->keep_internal;
            }
            if (fn == nullptr) return false;
            if (fn->has_children()) return true;
            if (full_c.ndim() == 3 && full_c.dim(0) == twok) {
                if (j > 0) return true;
                return windowed_halo_detail::wavelet_octants_normf(full_c, k_) > 0.0;
            }
            return false;
        };

        auto sort_and_pack_level = [&](int j) {
            auto& ld = levels_[static_cast<std::size_t>(j)];
            std::sort(ld.owned_real.begin(),
                      ld.owned_real.end(),
                      [](const auto& a, const auto& b) noexcept {
                          return windowed_halo_detail::windowed_key_less(a.key, b.key);
                      });
            if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                ld.real_packed_keys.resize(ld.owned_real.size());
                for (std::size_t ri = 0; ri < ld.owned_real.size(); ++ri) {
                    ld.real_packed_keys[ri] =
                        windowed_halo_detail::pack_spatial_key_same_level(ld.owned_real[ri].key);
                }
            }
        };

        // Level 0 (root): holds `(s^0, d^0)` when internal, or `s^0` (in `k^3` or `(2k)^3` s0) when leaf.
        {
            const Key<3>& k0 = cdata_.key0;
            const FunctionNode<T, 3>* fn0 = find_raw(0, k0);
            Tensor<T> c0 = extract_node_tensor(fn0);

            WindowedOwnedRealNode<T> root_rn;
            root_rn.key = k0;
            root_rn.is_internal = is_comp_internal(0, k0, fn0, c0);
            root_rn.s_out = Tensor<T>(kk, kk, kk);
            root_rn.s_in = (c0.size() > 0) ? windowed_halo_detail::extract_scaling_cube(c0, k_)
                                           : Tensor<T>(kk, kk, kk);
            windowed_halo_detail::compute_scaling_norms_readonly(
                root_rn.s_in, root_rn.s000, root_rn.s_osc_norm2, root_rn.s_normf);

            if (root_rn.is_internal) {
                Tensor<T> unf_2k(3L, twok_dims, false);
                Tensor<T> work_2k(3L, twok_dims, false);
                if (c0.size() == twok3) {
                    fast_transform(c0, cdata_.hg, unf_2k, work_2k);
                } else if (c0.size() == kk3) {
                    Tensor<T> full_ns(twok, twok, twok);
                    full_ns(s0, s0, s0) = root_rn.s_in;
                    fast_transform(full_ns, cdata_.hg, unf_2k, work_2k);
                } else {
                    unf_2k = T(0);
                }

                levels_[1].owned_real.reserve(8u);
                for (KeyChildIterator<3> kit(k0); kit; ++kit) {
                    const Key<3>& ch = kit.key();
                    const auto& cl = ch.translation();
                    const FunctionNode<T, 3>* fn_ch = find_raw(1, ch);
                    Tensor<T> c_ch = extract_node_tensor(fn_ch);

                    WindowedOwnedRealNode<T> ch_rn;
                    ch_rn.key = ch;
                    ch_rn.is_internal = is_comp_internal(1, ch, fn_ch, c_ch);
                    ch_rn.s_in = copy(unf_2k(cdata_.s[cl[0] & 1],
                                             cdata_.s[cl[1] & 1],
                                             cdata_.s[cl[2] & 1]));
                    ch_rn.s_out = Tensor<T>(kk, kk, kk);
                    if (ch_rn.is_internal && c_ch.size() == twok3) {
                        ch_rn.ns_in = c_ch;
                    }
                    levels_[1].owned_real.push_back(std::move(ch_rn));
                }
                sort_and_pack_level(1);

                const bool need_gamma_0 =
                    (0 >= serial_range_.min_gamma_level && 0 <= serial_range_.max_gamma_level);
                if (need_gamma_0) {
                    if (c0.size() == twok3) {
                        root_rn.ns_in = c0;
                        root_rn.d_normf = windowed_halo_detail::wavelet_octants_normf(c0, k_);
                        root_rn.ns_normf = std::sqrt(
                            root_rn.s_normf * root_rn.s_normf + root_rn.d_normf * root_rn.d_normf);
                    } else {
                        root_rn.ns_in = Tensor<T>(twok, twok, twok);
                        root_rn.ns_in(s0, s0, s0) = root_rn.s_in;
                        root_rn.d_normf = 0.0;
                        root_rn.ns_normf = root_rn.s_normf;
                    }
                } else {
                    root_rn.ns_normf = root_rn.s_normf;
                }
            } else {
                root_rn.ns_normf = root_rn.s_normf;
            }

            levels_[0].owned_real.push_back(std::move(root_rn));
            levels_[0].real_packed_keys = {windowed_halo_detail::pack_spatial_key_same_level(k0)};
        }

        // Top-down pass across levels `1 .. n_max_tree_`, expanding all 8 children of every
        // internal node in parallel without requiring explicit leaf entries in `f_comp.coeffs`.
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_full_2k;
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_unf_2k;
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_work_2k;
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> tls_w1_k3;

        for (int j = 1; j <= n_max_tree_; ++j) {
            auto& ld = levels_[static_cast<std::size_t>(j)];
            if (ld.owned_real.empty()) break;

            std::vector<std::size_t> int_slot(ld.owned_real.size(), SIZE_MAX);
            std::size_t n_int = 0;
            if (j < n_max_tree_) {
                for (std::size_t ri = 0; ri < ld.owned_real.size(); ++ri) {
                    if (ld.owned_real[ri].is_internal) {
                        int_slot[ri] = n_int++;
                    }
                }
            }
            std::vector<WindowedOwnedRealNode<T>> next_nodes(n_int * 8u);

            const bool need_gamma_j =
                (j < n_max_tree_ &&
                 j >= serial_range_.min_gamma_level &&
                 j <= serial_range_.max_gamma_level);

            windowed_parallel_detail::parallel_for_chunks(
                world_,
                ld.owned_real.size(),
                windowed_parallel_detail::kDefaultParallelGrain,
                [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
                    for (std::size_t ri = b; ri < e_hi; ++ri) {
                        auto& rn = ld.owned_real[ri];
                        windowed_halo_detail::compute_scaling_norms_readonly(
                            rn.s_in, rn.s000, rn.s_osc_norm2, rn.s_normf);

                        if (!rn.is_internal || j >= n_max_tree_) {
                            rn.is_internal = false;
                            rn.ns_in = Tensor<T>();
                            rn.ns_normf = rn.s_normf;
                            continue;
                        }

                        const std::size_t base_out = int_slot[ri] * 8u;
                        if (rn.ns_in.size() == twok3) {
                            if (tls_full_2k[cid].size() != twok3) {
                                tls_full_2k[cid] = Tensor<T>(3L, twok_dims, false);
                                tls_unf_2k[cid] = Tensor<T>(3L, twok_dims, false);
                                tls_work_2k[cid] = Tensor<T>(3L, twok_dims, false);
                            }
                            Tensor<T>& full_ns = tls_full_2k[cid];
                            Tensor<T>& unf_2k = tls_unf_2k[cid];
                            Tensor<T>& work_2k = tls_work_2k[cid];

                            std::memcpy(full_ns.ptr(),
                                        rn.ns_in.ptr(),
                                        static_cast<std::size_t>(twok3) * sizeof(T));
                            full_ns(s0, s0, s0) = rn.s_in;
                            fast_transform(full_ns, cdata_.hg, unf_2k, work_2k);

                            std::size_t oct = 0;
                            for (KeyChildIterator<3> kit(rn.key); kit; ++kit, ++oct) {
                                const Key<3>& ch = kit.key();
                                const auto& cl = ch.translation();
                                const FunctionNode<T, 3>* fn_ch = find_raw(j + 1, ch);
                                Tensor<T> c_ch = extract_node_tensor(fn_ch);

                                WindowedOwnedRealNode<T>& ch_rn = next_nodes[base_out + oct];
                                ch_rn.key = ch;
                                ch_rn.is_internal = is_comp_internal(j + 1, ch, fn_ch, c_ch);
                                ch_rn.s_in = copy(unf_2k(cdata_.s[cl[0] & 1],
                                                         cdata_.s[cl[1] & 1],
                                                         cdata_.s[cl[2] & 1]));
                                ch_rn.s_out = Tensor<T>(kk, kk, kk);
                                if (ch_rn.is_internal && c_ch.size() == twok3) {
                                    ch_rn.ns_in = c_ch;
                                }
                            }

                            if (need_gamma_j) {
                                rn.d_normf =
                                    windowed_halo_detail::wavelet_octants_normf(full_ns, k_);
                                rn.ns_in = copy(full_ns);
                                rn.ns_normf = std::sqrt(
                                    rn.s_normf * rn.s_normf + rn.d_normf * rn.d_normf);
                            } else {
                                rn.ns_in = Tensor<T>();
                                rn.ns_normf = rn.s_normf;
                            }
                        } else {
                            if (tls_w1_k3[cid].size() != kk3) {
                                tls_w1_k3[cid] = Tensor<T>(3L, k_dims, false);
                            }
                            Tensor<T>& w1_k3 = tls_w1_k3[cid];

                            std::size_t oct = 0;
                            for (KeyChildIterator<3> kit(rn.key); kit; ++kit, ++oct) {
                                const Key<3>& ch = kit.key();
                                const FunctionNode<T, 3>* fn_ch = find_raw(j + 1, ch);
                                Tensor<T> c_ch = extract_node_tensor(fn_ch);

                                WindowedOwnedRealNode<T>& ch_rn = next_nodes[base_out + oct];
                                ch_rn.key = ch;
                                ch_rn.is_internal = is_comp_internal(j + 1, ch, fn_ch, c_ch);
                                ch_rn.s_in = Tensor<T>(3L, k_dims, false);
                                windowed_halo_detail::prolong_step_3d(
                                    h_oct_, rn.s_in, ch, ch_rn.s_in, w1_k3);
                                ch_rn.s_out = Tensor<T>(kk, kk, kk);
                                if (ch_rn.is_internal && c_ch.size() == twok3) {
                                    ch_rn.ns_in = c_ch;
                                }
                            }

                            if (need_gamma_j) {
                                rn.ns_in = Tensor<T>(twok, twok, twok);
                                rn.ns_in(s0, s0, s0) = rn.s_in;
                                rn.d_normf = 0.0;
                                rn.ns_normf = rn.s_normf;
                            } else {
                                rn.ns_in = Tensor<T>();
                                rn.ns_normf = rn.s_normf;
                            }
                        }
                    }
                });

            if (n_int > 0) {
                levels_[static_cast<std::size_t>(j + 1)].owned_real = std::move(next_nodes);
                sort_and_pack_level(j + 1);
            }
        }

        finalize_serial_stats_and_child_masks();
    }

    // =========================================================================
    // Serial Single-Rank Build (`world.size() == 1`) — Parallel + 64-Bit Packed Keys
    // =========================================================================
    void build_serial(const Function<T, 3>& f_in, const RadiiBuilderFn& radii_builder) {
        const TreeState in_state = f_in.get_impl()->get_tree_state();
        if (in_state == reconstructed) {
            ingest_local_real_nodes_serial_reconstructed(f_in);
        } else if (in_state == compressed) {
            ingest_local_real_nodes_serial_compressed(f_in);
        } else {
            ingest_local_real_nodes(f_in);
        }
        radii_builder(global_stats_, radii_);

        const std::size_t n_lvls = static_cast<std::size_t>(n_max_tree_ + 1);
        std::vector<KeyFlagVec> virt_flags(n_lvls);
        if (n_max_tree_ < 1) {
            populate_prolonged_virtual_and_sources(virt_flags, {});
            return;
        }

        constexpr std::size_t kHaloChunkNodes = 16;
        struct LevelChunkSlice {
            int j = 0;
            std::size_t begin = 0;
            std::size_t end = 0;
        };

        // Pass 1: Parallel cross-level flattened H1 discovery + `near_internal_4b1/4b2` marking.
        std::vector<LevelChunkSlice> h1_slices;
        for (int j = 1; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const auto& ld = levels_[sj];
            if (ld.owned_real.empty()) continue;
            const WindowedLevelRadii& rad = radii_[sj];
            const auto& st = global_stats_[sj];
            const bool can_refine_j = dynamic_refine_ && (j < max_refine_level_);
            const bool need_T = (rad.has_T_terms && rad.m_max_T >= 1);
            const int m1 = (ld.global_has_internal && rad.has_4b1_terms)
                               ? rad.node_m4b1(dynamic_refine_, rad.m_max_4b1, st.max_d_int)
                               : -1;
            const int m2 =
                (ld.global_has_internal && can_refine_j && rad.max_refine_cand_bw() >= 1) ? 1 : -1;
            if (!need_T && std::max(m1, m2) < 1) continue;

            for (std::size_t b = 0; b < ld.owned_real.size(); b += kHaloChunkNodes) {
                h1_slices.push_back(LevelChunkSlice{
                    j, b, std::min(b + kHaloChunkNodes, ld.owned_real.size())});
            }
        }

        std::array<std::vector<std::vector<windowed_halo_detail::PackedKeyFlag>>,
                   windowed_parallel_detail::kMaxParallelChunks>
            h1_packed_tls;
        std::array<std::vector<KeyFlagVec>, windowed_parallel_detail::kMaxParallelChunks> h1_key_tls;
        std::array<std::vector<std::vector<std::pair<std::uint32_t, std::uint8_t>>>,
                   windowed_parallel_detail::kMaxParallelChunks>
            near_real_tls;

        if (!h1_slices.empty()) {
            windowed_parallel_detail::parallel_for_chunks(
                world_, h1_slices.size(), 1, [&](std::size_t sb, std::size_t se, std::size_t cid) {
                    auto& my_h1_packed = h1_packed_tls[cid];
                    auto& my_h1_keys = h1_key_tls[cid];
                    auto& my_near = near_real_tls[cid];
                    if (my_h1_packed.empty()) {
                        my_h1_packed.resize(n_lvls);
                        my_h1_keys.resize(n_lvls);
                        my_near.resize(n_lvls);
                    }

                    for (std::size_t si = sb; si < se; ++si) {
                        const LevelChunkSlice sl = h1_slices[si];
                        const int j = sl.j;
                        const std::size_t sj = static_cast<std::size_t>(j);
                        const auto& ld = levels_[sj];
                        const WindowedLevelRadii& rad = radii_[sj];
                        const bool can_refine_j = dynamic_refine_ && (j < max_refine_level_);
                        const bool use_packed = (j <= windowed_halo_detail::kMaxPackedKeyLevel);
                        const Translation dom =
                            use_packed ? (Translation(1) << j) : Translation(0);

                        if (rad.has_T_terms && rad.m_max_T >= 1) {
                            for (std::size_t ri = sl.begin; ri < sl.end; ++ri) {
                                const auto& rn = ld.owned_real[ri];
                                const auto [mT_rn, max_l1_T] = rad.node_mT_and_l1(
                                    dynamic_refine_, rn.s_normf, rn.s000, rn.s_osc_norm2);
                                if (mT_rn < 1) continue;
                                if (use_packed) {
                                    const auto& l = rn.key.translation();
                                    auto& out_p = my_h1_packed[sj];
                                    for (int dx = -mT_rn; dx <= mT_rn; ++dx) {
                                        for (int dy = -mT_rn; dy <= mT_rn; ++dy) {
                                            for (int dz = -mT_rn; dz <= mT_rn; ++dz) {
                                                const int l1 =
                                                    std::abs(dx) + std::abs(dy) + std::abs(dz);
                                                if (l1 == 0 || l1 > max_l1_T) continue;
                                                std::uint64_t nb_code = 0;
                                                if (!windowed_halo_detail::shift_packed_coords_in_domain(
                                                        dom,
                                                        l[0],
                                                        l[1],
                                                        l[2],
                                                        dx,
                                                        dy,
                                                        dz,
                                                        is_periodic_,
                                                        nb_code)) {
                                                    continue;
                                                }
                                                if (find_owned_real_packed(j, nb_code) == SIZE_MAX) {
                                                    out_p.push_back(
                                                        windowed_halo_detail::PackedKeyFlag{
                                                            nb_code, kHaloH1T});
                                                }
                                            }
                                        }
                                    }
                                } else {
                                    auto& out_k = my_h1_keys[sj];
                                    for (int dx = -mT_rn; dx <= mT_rn; ++dx) {
                                        for (int dy = -mT_rn; dy <= mT_rn; ++dy) {
                                            for (int dz = -mT_rn; dz <= mT_rn; ++dz) {
                                                const int l1 =
                                                    std::abs(dx) + std::abs(dy) + std::abs(dz);
                                                if (l1 == 0 || l1 > max_l1_T) continue;
                                                Key<3> nb;
                                                if (!windowed_halo_detail::shift_key_in_domain(
                                                        rn.key, dx, dy, dz, is_periodic_, nb)) {
                                                    continue;
                                                }
                                                if (find_owned_real_index(j, nb) == SIZE_MAX) {
                                                    out_k.emplace_back(nb, kHaloH1T);
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }

                        if (ld.global_has_internal) {
                            const int m1 = rad.has_4b1_terms ? rad.m_max_4b1 : -1;
                            const int m2 =
                                (can_refine_j && rad.max_refine_cand_bw() >= 1) ? 1 : -1;
                            if (std::max(m1, m2) >= 1) {
                                auto& out_near = my_near[sj];
                                for (std::size_t ri = sl.begin; ri < sl.end; ++ri) {
                                    const auto& rn = ld.owned_real[ri];
                                    if (!rn.is_internal) continue;
                                    const int m1_rn =
                                        rad.node_m4b1(dynamic_refine_, m1, rn.d_normf);
                                    const int m2_rn = rad.node_mRefine(dynamic_refine_,
                                                                       m2,
                                                                       rn.s_normf,
                                                                       rn.s000,
                                                                       rn.s_osc_norm2,
                                                                       rn.d_normf);
                                    const int m_emit_rn = std::max(m1_rn, m2_rn);
                                    if (m_emit_rn < 1) continue;

                                    if (use_packed) {
                                        const auto& l = rn.key.translation();
                                        auto& out_p = my_h1_packed[sj];
                                        for (int dx = -m_emit_rn; dx <= m_emit_rn; ++dx) {
                                            for (int dy = -m_emit_rn; dy <= m_emit_rn; ++dy) {
                                                for (int dz = -m_emit_rn; dz <= m_emit_rn; ++dz) {
                                                    if (dx == 0 && dy == 0 && dz == 0) continue;
                                                    const int sh = std::max(
                                                        {std::abs(dx), std::abs(dy), std::abs(dz)});
                                                    std::uint64_t nb_code = 0;
                                                    if (!windowed_halo_detail::shift_packed_coords_in_domain(
                                                            dom,
                                                            l[0],
                                                            l[1],
                                                            l[2],
                                                            dx,
                                                            dy,
                                                            dz,
                                                            is_periodic_,
                                                            nb_code)) {
                                                        continue;
                                                    }
                                                    const std::size_t ridx =
                                                        find_owned_real_packed(j, nb_code);
                                                    if (ridx != SIZE_MAX) {
                                                        if (!ld.owned_real[ridx].is_internal) {
                                                            std::uint8_t nbits = 0u;
                                                            if (sh <= m1_rn) nbits |= 1u;
                                                            if (sh <= m2_rn &&
                                                                (!ungraded_leaf_refine_only_ ||
                                                                 windowed_halo_detail::has_ungraded_child_toward(
                                                                     rn.internal_child_mask, dx, dy, dz))) {
                                                                nbits |= 2u;
                                                            }
                                                            if (nbits != 0u) {
                                                                out_near.emplace_back(
                                                                    static_cast<std::uint32_t>(ridx),
                                                                    nbits);
                                                            }
                                                        }
                                                    } else {
                                                        std::uint8_t bits = 0u;
                                                        if (sh <= m1_rn) {
                                                            bits = static_cast<std::uint8_t>(
                                                                bits | kHaloH1Gamma);
                                                        }
                                                        if (sh <= m2_rn) {
                                                            bits = static_cast<std::uint8_t>(
                                                                bits | kHaloH1Refine);
                                                        }
                                                        if (bits != 0u) {
                                                            out_p.push_back(
                                                                windowed_halo_detail::PackedKeyFlag{
                                                                    nb_code, bits});
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    } else {
                                        auto& out_k = my_h1_keys[sj];
                                        for (int dx = -m_emit_rn; dx <= m_emit_rn; ++dx) {
                                            for (int dy = -m_emit_rn; dy <= m_emit_rn; ++dy) {
                                                for (int dz = -m_emit_rn; dz <= m_emit_rn; ++dz) {
                                                    if (dx == 0 && dy == 0 && dz == 0) continue;
                                                    const int sh = std::max(
                                                        {std::abs(dx), std::abs(dy), std::abs(dz)});
                                                    Key<3> nb;
                                                    if (!windowed_halo_detail::shift_key_in_domain(
                                                            rn.key, dx, dy, dz, is_periodic_, nb)) {
                                                        continue;
                                                    }
                                                    const std::size_t ridx =
                                                        find_owned_real_index(j, nb);
                                                    if (ridx != SIZE_MAX) {
                                                        if (!ld.owned_real[ridx].is_internal) {
                                                            std::uint8_t nbits = 0u;
                                                            if (sh <= m1_rn) nbits |= 1u;
                                                            if (sh <= m2_rn &&
                                                                (!ungraded_leaf_refine_only_ ||
                                                                 windowed_halo_detail::has_ungraded_child_toward(
                                                                     rn.internal_child_mask, dx, dy, dz))) {
                                                                nbits |= 2u;
                                                            }
                                                            if (nbits != 0u) {
                                                                out_near.emplace_back(
                                                                    static_cast<std::uint32_t>(ridx),
                                                                    nbits);
                                                            }
                                                        }
                                                    } else {
                                                        std::uint8_t bits = 0u;
                                                        if (sh <= m1_rn) {
                                                            bits = static_cast<std::uint8_t>(
                                                                bits | kHaloH1Gamma);
                                                        }
                                                        if (sh <= m2_rn) {
                                                            bits = static_cast<std::uint8_t>(
                                                                bits | kHaloH1Refine);
                                                        }
                                                        if (bits != 0u) {
                                                            out_k.emplace_back(nb, bits);
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
        }

        std::vector<std::vector<windowed_halo_detail::PackedKeyFlag>> h1_packed_lvls(n_lvls);
        std::vector<KeyFlagVec> h1_key_lvls(n_lvls);

        const std::size_t h1_slot_lim = windowed_parallel_detail::slot_bound();
        for (int j = 1; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            auto& ld = levels_[sj];
            for (std::size_t cid = 0; cid < h1_slot_lim; ++cid) {
                if (!near_real_tls[cid].empty()) {
                    for (const auto& [ridx, nbits] : near_real_tls[cid][sj]) {
                        auto& trn = ld.owned_real[ridx];
                        if ((nbits & 1u) != 0u) trn.near_internal_4b1 = true;
                        if ((nbits & 2u) != 0u) trn.near_internal_4b2 = true;
                    }
                }
            }
            if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                std::size_t total_p = 0;
                for (std::size_t cid = 0; cid < h1_slot_lim; ++cid) {
                    if (!h1_packed_tls[cid].empty()) total_p += h1_packed_tls[cid][sj].size();
                }
                auto& h1p = h1_packed_lvls[sj];
                h1p.reserve(total_p);
                for (std::size_t cid = 0; cid < h1_slot_lim; ++cid) {
                    if (!h1_packed_tls[cid].empty()) {
                        const auto& src = h1_packed_tls[cid][sj];
                        h1p.insert(h1p.end(), src.begin(), src.end());
                    }
                }
                windowed_halo_detail::sort_and_merge_packed_flags(h1p);
            } else {
                std::size_t total_k = 0;
                for (std::size_t cid = 0; cid < h1_slot_lim; ++cid) {
                    if (!h1_key_tls[cid].empty()) total_k += h1_key_tls[cid][sj].size();
                }
                auto& h1k = h1_key_lvls[sj];
                h1k.reserve(total_k);
                for (std::size_t cid = 0; cid < h1_slot_lim; ++cid) {
                    if (!h1_key_tls[cid].empty()) {
                        const auto& src = h1_key_tls[cid][sj];
                        h1k.insert(h1k.end(), src.begin(), src.end());
                    }
                }
                windowed_halo_detail::sort_and_merge_key_flags(h1k);
            }
        }

        // Pass 2: Parallel cross-level flattened H2_Ghost completion.
        std::vector<std::vector<std::pair<std::uint64_t, int>>> h2_centers_packed(n_lvls);
        std::vector<std::vector<std::pair<Key<3>, int>>> h2_centers_keys(n_lvls);
        std::vector<LevelChunkSlice> h2_slices;

        for (int j = 1; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const auto& ld = levels_[sj];
            const WindowedLevelRadii& rad = radii_[sj];
            const bool can_refine_j = dynamic_refine_ && (j < max_refine_level_);
            const double max_cl_above = global_stats_[sj].max_coarse_leaf_above;

            const int mT_h2 = rad.h2_mT(dynamic_refine_, max_cl_above);
            const int mLeaf_h2 =
                rad.h2_mLeaf(dynamic_refine_, ld.global_has_internal, max_cl_above);
            if (std::max(mT_h2, mLeaf_h2) < 1) continue;

            if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                auto& cvec = h2_centers_packed[sj];
                const auto& h1p = h1_packed_lvls[sj];
                cvec.reserve(ld.owned_real.size() + h1p.size());
                if (mLeaf_h2 >= 1) {
                    for (std::size_t ri = 0; ri < ld.owned_real.size(); ++ri) {
                        const auto& rn = ld.owned_real[ri];
                        if (rn.is_internal || (can_refine_j && rn.near_internal_4b2)) {
                            cvec.emplace_back(ld.real_packed_keys[ri], mLeaf_h2);
                        }
                    }
                }
                if (mT_h2 >= 1 || (can_refine_j && mLeaf_h2 >= 1)) {
                    for (const auto& pf : h1p) {
                        int r_h2 = -1;
                        if (mT_h2 >= 1 && (pf.flags & kHaloH1T) != 0) {
                            r_h2 = std::max(r_h2, mT_h2);
                        }
                        if (can_refine_j && mLeaf_h2 >= 1 && (pf.flags & kHaloH1Refine) != 0) {
                            r_h2 = std::max(r_h2, mLeaf_h2);
                        }
                        if (r_h2 >= 1) {
                            cvec.emplace_back(pf.code, r_h2);
                        }
                    }
                }
                for (std::size_t b = 0; b < cvec.size(); b += kHaloChunkNodes) {
                    h2_slices.push_back(
                        LevelChunkSlice{j, b, std::min(b + kHaloChunkNodes, cvec.size())});
                }
            } else {
                auto& cvec = h2_centers_keys[sj];
                const auto& h1k = h1_key_lvls[sj];
                cvec.reserve(ld.owned_real.size() + h1k.size());
                if (mLeaf_h2 >= 1) {
                    for (const auto& rn : ld.owned_real) {
                        if (rn.is_internal || (can_refine_j && rn.near_internal_4b2)) {
                            cvec.emplace_back(rn.key, mLeaf_h2);
                        }
                    }
                }
                if (mT_h2 >= 1 || (can_refine_j && mLeaf_h2 >= 1)) {
                    for (const auto& kv : h1k) {
                        int r_h2 = -1;
                        if (mT_h2 >= 1 && (kv.second & kHaloH1T) != 0) {
                            r_h2 = std::max(r_h2, mT_h2);
                        }
                        if (can_refine_j && mLeaf_h2 >= 1 && (kv.second & kHaloH1Refine) != 0) {
                            r_h2 = std::max(r_h2, mLeaf_h2);
                        }
                        if (r_h2 >= 1) {
                            cvec.emplace_back(kv.first, r_h2);
                        }
                    }
                }
                for (std::size_t b = 0; b < cvec.size(); b += kHaloChunkNodes) {
                    h2_slices.push_back(
                        LevelChunkSlice{j, b, std::min(b + kHaloChunkNodes, cvec.size())});
                }
            }
        }

        std::array<std::vector<std::vector<std::uint64_t>>,
                   windowed_parallel_detail::kMaxParallelChunks>
            h2_packed_tls;
        std::array<std::vector<std::vector<Key<3>>>, windowed_parallel_detail::kMaxParallelChunks>
            h2_key_tls;

        if (!h2_slices.empty()) {
            windowed_parallel_detail::parallel_for_chunks(
                world_, h2_slices.size(), 1, [&](std::size_t sb, std::size_t se, std::size_t cid) {
                    auto& my_h2_packed = h2_packed_tls[cid];
                    auto& my_h2_keys = h2_key_tls[cid];
                    if (my_h2_packed.empty()) {
                        my_h2_packed.resize(n_lvls);
                        my_h2_keys.resize(n_lvls);
                    }

                    for (std::size_t si = sb; si < se; ++si) {
                        const LevelChunkSlice sl = h2_slices[si];
                        const int j = sl.j;
                        const std::size_t sj = static_cast<std::size_t>(j);
                        if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                            const Translation dom = (Translation(1) << j);
                            const auto& cvec = h2_centers_packed[sj];
                            const auto& h1p = h1_packed_lvls[sj];
                            auto& out_h2 = my_h2_packed[sj];
                            for (std::size_t ci = sl.begin; ci < sl.end; ++ci) {
                                const auto [center_code, radius] = cvec[ci];
                                for (int dx = -radius; dx <= radius; ++dx) {
                                    for (int dy = -radius; dy <= radius; ++dy) {
                                        for (int dz = -radius; dz <= radius; ++dz) {
                                            if (dx == 0 && dy == 0 && dz == 0) continue;
                                            std::uint64_t nb_code = 0;
                                            if (!windowed_halo_detail::shift_packed_code_in_domain(
                                                    dom,
                                                    center_code,
                                                    dx,
                                                    dy,
                                                    dz,
                                                    is_periodic_,
                                                    nb_code)) {
                                                continue;
                                            }
                                            if (find_owned_real_packed(j, nb_code) == SIZE_MAX &&
                                                !windowed_halo_detail::contains_sorted_packed_flag(
                                                    h1p, nb_code)) {
                                                out_h2.push_back(nb_code);
                                            }
                                        }
                                    }
                                }
                            }
                        } else {
                            const auto& cvec = h2_centers_keys[sj];
                            const auto& h1k = h1_key_lvls[sj];
                            auto& out_h2 = my_h2_keys[sj];
                            for (std::size_t ci = sl.begin; ci < sl.end; ++ci) {
                                const auto& [center_key, radius] = cvec[ci];
                                for (int dx = -radius; dx <= radius; ++dx) {
                                    for (int dy = -radius; dy <= radius; ++dy) {
                                        for (int dz = -radius; dz <= radius; ++dz) {
                                            if (dx == 0 && dy == 0 && dz == 0) continue;
                                            Key<3> nb;
                                            if (!windowed_halo_detail::shift_key_in_domain(
                                                    center_key, dx, dy, dz, is_periodic_, nb)) {
                                                continue;
                                            }
                                            if (find_owned_real_index(j, nb) == SIZE_MAX &&
                                                !windowed_halo_detail::contains_sorted_key_flag(
                                                    h1k, nb)) {
                                                out_h2.push_back(nb);
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                });
        }

        // Merge sorted unique `(H1, H2_Ghost)` per level, unpacking `Key<3>` ONCE per unique virtual node.
        const std::size_t h2_slot_lim = windowed_parallel_detail::slot_bound();
        for (int j = 1; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            auto& vf = virt_flags[sj];

            if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                const Level lvl = static_cast<Level>(j);
                const auto& h1p = h1_packed_lvls[sj];
                std::size_t total_h2 = 0;
                for (std::size_t cid = 0; cid < h2_slot_lim; ++cid) {
                    if (!h2_packed_tls[cid].empty()) total_h2 += h2_packed_tls[cid][sj].size();
                }
                std::vector<std::uint64_t> h2p;
                h2p.reserve(total_h2);
                for (std::size_t cid = 0; cid < h2_slot_lim; ++cid) {
                    if (!h2_packed_tls[cid].empty()) {
                        const auto& src = h2_packed_tls[cid][sj];
                        h2p.insert(h2p.end(), src.begin(), src.end());
                    }
                }
                windowed_halo_detail::sort_unique_packed_codes(h2p);

                vf.reserve(h1p.size() + h2p.size());
                std::size_t i1 = 0;
                std::size_t i2 = 0;
                while (i1 < h1p.size() && i2 < h2p.size()) {
                    if (h1p[i1].code < h2p[i2]) {
                        vf.emplace_back(
                            windowed_halo_detail::unpack_spatial_key_same_level(lvl, h1p[i1].code),
                            h1p[i1].flags);
                        ++i1;
                    } else if (h2p[i2] < h1p[i1].code) {
                        vf.emplace_back(
                            windowed_halo_detail::unpack_spatial_key_same_level(lvl, h2p[i2]),
                            kHaloH2Ghost);
                        ++i2;
                    } else {
                        vf.emplace_back(
                            windowed_halo_detail::unpack_spatial_key_same_level(lvl, h1p[i1].code),
                            h1p[i1].flags);
                        ++i1;
                        ++i2;
                    }
                }
                while (i1 < h1p.size()) {
                    vf.emplace_back(
                        windowed_halo_detail::unpack_spatial_key_same_level(lvl, h1p[i1].code),
                        h1p[i1].flags);
                    ++i1;
                }
                while (i2 < h2p.size()) {
                    vf.emplace_back(
                        windowed_halo_detail::unpack_spatial_key_same_level(lvl, h2p[i2]),
                        kHaloH2Ghost);
                    ++i2;
                }
            } else {
                const auto& h1k = h1_key_lvls[sj];
                std::size_t total_h2 = 0;
                for (std::size_t cid = 0; cid < h2_slot_lim; ++cid) {
                    if (!h2_key_tls[cid].empty()) total_h2 += h2_key_tls[cid][sj].size();
                }
                std::vector<Key<3>> h2k;
                h2k.reserve(total_h2);
                for (std::size_t cid = 0; cid < h2_slot_lim; ++cid) {
                    if (!h2_key_tls[cid].empty()) {
                        const auto& src = h2_key_tls[cid][sj];
                        h2k.insert(h2k.end(), src.begin(), src.end());
                    }
                }
                windowed_halo_detail::sort_unique_keys(h2k);

                vf.reserve(h1k.size() + h2k.size());
                std::size_t i1 = 0;
                std::size_t i2 = 0;
                while (i1 < h1k.size() && i2 < h2k.size()) {
                    if (windowed_halo_detail::windowed_key_less(h1k[i1].first, h2k[i2])) {
                        vf.push_back(h1k[i1++]);
                    } else if (windowed_halo_detail::windowed_key_less(h2k[i2], h1k[i1].first)) {
                        vf.emplace_back(h2k[i2++], kHaloH2Ghost);
                    } else {
                        vf.push_back(h1k[i1++]);
                        ++i2;
                    }
                }
                while (i1 < h1k.size()) vf.push_back(h1k[i1++]);
                while (i2 < h2k.size()) vf.emplace_back(h2k[i2++], kHaloH2Ghost);
            }
        }

        populate_prolonged_virtual_and_sources(virt_flags, {});
    }

    // =========================================================================
    // Multi-Rank Distributed Build (`world.size() > 1`)
    // =========================================================================
    void build_distributed(const Function<T, 3>& f_ns, const RadiiBuilderFn& radii_builder) {
        MADNESS_CHECK(exchange_ptr_ != nullptr);
        build_distributed_lockstep(world_, {this}, {&f_ns}, {radii_builder}, *exchange_ptr_);
    }

    // =========================================================================
    // Shared Top-Down Single-Step Interior Prolongation (`virt_cache`)
    // =========================================================================
    void populate_prolonged_virtual_and_sources(
        const std::vector<KeyFlagVec>& virt_flags,
        const std::vector<WindowedSourceEntry<T>>& remote_real_nodes) {
        const long kk = static_cast<long>(k_);

        struct VirtProlongItem {
            Key<3> key;
            Key<3> coarse_leaf;
            T s000 = T(0);
            double s_osc_norm2 = 0.0;
            double s_normf = 0.0;
            Tensor<T> s_virt;
        };

        std::vector<std::vector<VirtProlongItem>> prolong_levels(
            static_cast<std::size_t>(n_max_tree_ + 1));

        for (int j = 1; j <= n_max_tree_; ++j) {
            const auto& vf = virt_flags[static_cast<std::size_t>(j)];
            auto& pl = prolong_levels[static_cast<std::size_t>(j)];
            pl.reserve(vf.size());
            for (const auto& kv : vf) {
                pl.push_back(VirtProlongItem{kv.first, Key<3>(), T(0), 0.0, 0.0, Tensor<T>()});
            }
        }

        // Single bottom-up pass (`j = n_max_tree .. 2`) inserting intermediate virtual ancestors.
        std::vector<Key<3>> needed_parents;
        for (int j = n_max_tree_; j >= 2; --j) {
            const auto& pl_j = prolong_levels[static_cast<std::size_t>(j)];
            if (pl_j.empty()) continue;

            needed_parents.clear();
            for (const auto& item : pl_j) {
                const Key<3> par = item.key.parent();
                if (find_coarse_leaf(j - 1, par) == nullptr) {
                    needed_parents.push_back(par);
                }
            }
            if (needed_parents.empty()) continue;
            windowed_halo_detail::sort_unique_keys(needed_parents);

            auto& pl_par = prolong_levels[static_cast<std::size_t>(j - 1)];
            std::vector<VirtProlongItem> merged;
            merged.reserve(pl_par.size() + needed_parents.size());
            std::size_t i_exist = 0;
            std::size_t i_new = 0;
            while (i_exist < pl_par.size() && i_new < needed_parents.size()) {
                if (windowed_halo_detail::windowed_key_less(pl_par[i_exist].key, needed_parents[i_new])) {
                    merged.push_back(std::move(pl_par[i_exist++]));
                } else if (windowed_halo_detail::windowed_key_less(
                               needed_parents[i_new], pl_par[i_exist].key)) {
                    merged.push_back(VirtProlongItem{
                        needed_parents[i_new++], Key<3>(), T(0), 0.0, 0.0, Tensor<T>()});
                } else {
                    merged.push_back(std::move(pl_par[i_exist++]));
                    ++i_new;
                }
            }
            while (i_exist < pl_par.size()) merged.push_back(std::move(pl_par[i_exist++]));
            while (i_new < needed_parents.size()) {
                merged.push_back(VirtProlongItem{
                    needed_parents[i_new++], Key<3>(), T(0), 0.0, 0.0, Tensor<T>()});
            }
            pl_par.swap(merged);
        }

        // Single top-down pass (`j = 1 .. n_max_tree`) prolonging each virtual box once via `general_fast_transform`.
        const long k_dims[3] = {kk, kk, kk};
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> prolong_work1;
        for (int j = 1; j <= n_max_tree_; ++j) {
            auto& pl_j = prolong_levels[static_cast<std::size_t>(j)];
            if (pl_j.empty()) continue;
            const auto& pl_par = prolong_levels[static_cast<std::size_t>(j - 1)];
            windowed_parallel_detail::parallel_for_chunks(
                world_,
                pl_j.size(),
                8u,
                [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
                    if (prolong_work1[cid].size() != kk * kk * kk) {
                        prolong_work1[cid] = Tensor<T>(3L, k_dims, false);
                    }
                    Tensor<T>& w1 = prolong_work1[cid];
                    for (std::size_t idx = b; idx < e_hi; ++idx) {
                        auto& item = pl_j[idx];
                        const Key<3> par = item.key.parent();
                        const Tensor<T>* parent_s_ptr = find_coarse_leaf(j - 1, par);
                        if (parent_s_ptr != nullptr) {
                            item.coarse_leaf = par;
                        } else {
                            const auto pit = std::lower_bound(
                                pl_par.begin(), pl_par.end(), par,
                                [](const VirtProlongItem& vp, const Key<3>& k) noexcept {
                                    return windowed_halo_detail::windowed_key_less(vp.key, k);
                                });
                            MADNESS_CHECK(pit != pl_par.end() && pit->key == par);
                            item.coarse_leaf = pit->coarse_leaf;
                            parent_s_ptr = &(pit->s_virt);
                        }

                        item.s_virt = Tensor<T>(3L, k_dims, false);
                        windowed_halo_detail::prolong_step_3d(
                            h_oct_, *parent_s_ptr, item.key, item.s_virt, w1);
                        windowed_halo_detail::compute_scaling_norms_readonly(
                            item.s_virt, item.s000, item.s_osc_norm2, item.s_normf);
                    }
                });
        }

        // Partition remote real nodes by level for O(N) sorted merging.
        std::vector<std::vector<WindowedSourceEntry<T>>> remote_by_level(
            static_cast<std::size_t>(n_max_tree_ + 1));
        for (const auto& se : remote_real_nodes) {
            const int j = static_cast<int>(se.key.level());
            if (j >= 0 && j <= n_max_tree_) {
                remote_by_level[static_cast<std::size_t>(j)].push_back(se);
            }
        }

        // Populate per-level `owned_virt` (already sorted) and sorted `sources` vectors.
        for (int j = 0; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            auto& ld = levels_[sj];
            const auto& vf = virt_flags[sj];
            const auto& rem_j = remote_by_level[sj];

            if (world_.size() == 1 && vf.empty()) {
                const auto& rad = radii_[sj];
                const bool need_sources_j =
                    (j >= 1 && rad.has_T_terms) ||
                    rad.has_4b1_terms ||
                    rad.has_4b2_int_terms ||
                    rad.has_4b2_leaf_terms;
                if (!need_sources_j) {
                    continue;
                }
            }

            ld.sources.reserve(ld.owned_real.size() + vf.size() + rem_j.size());

            auto make_real_source = [](const WindowedOwnedRealNode<T>& rn) {
                WindowedSourceEntry<T> se;
                se.key = rn.key;
                se.is_internal = rn.is_internal;
                se.s000 = rn.s000;
                se.s_osc_norm2 = rn.s_osc_norm2;
                se.s_normf = rn.s_normf;
                se.d_normf = rn.d_normf;
                se.ns_normf = rn.ns_normf;
                se.s_coeff = rn.s_in;
                if (rn.is_internal) {
                    se.ns_coeff = rn.ns_in;
                }
                return se;
            };

            if (j >= 1 && !vf.empty()) {
                const auto& pl_j = prolong_levels[sj];
                ld.owned_virt.reserve(vf.size());
                std::size_t p_cursor = 0;
                std::size_t r_idx = 0;

                auto emit_virt_entry = [&](const Key<3>& vkey, std::uint8_t flags) {
                    while (p_cursor < pl_j.size() &&
                           windowed_halo_detail::windowed_key_less(pl_j[p_cursor].key, vkey)) {
                        ++p_cursor;
                    }
                    MADNESS_CHECK(p_cursor < pl_j.size() && pl_j[p_cursor].key == vkey);
                    const VirtProlongItem& pitem = pl_j[p_cursor];

                    WindowedSourceEntry<T> se;
                    se.key = vkey;
                    se.is_internal = false;
                    se.is_virtual = true;
                    se.s000 = pitem.s000;
                    se.s_osc_norm2 = pitem.s_osc_norm2;
                    se.s_normf = pitem.s_normf;
                    se.ns_normf = pitem.s_normf;
                    se.s_coeff = pitem.s_virt;
                    ld.sources.push_back(std::move(se));

                    const std::uint8_t h1_bits =
                        static_cast<std::uint8_t>(flags & (kHaloH1T | kHaloH1Gamma | kHaloH1Refine));
                    if (h1_bits != 0) {
                        WindowedOwnedVirtNode<T> vn;
                        vn.key = vkey;
                        vn.coarse_leaf_key = pitem.coarse_leaf;
                        vn.flags = h1_bits;
                        vn.s000 = pitem.s000;
                        vn.s_osc_norm2 = pitem.s_osc_norm2;
                        vn.s_normf = pitem.s_normf;
                        vn.s_virt = pitem.s_virt;
                        if ((h1_bits & (kHaloH1T | kHaloH1Gamma)) != 0) {
                            vn.ds_corr = Tensor<T>(kk, kk, kk);
                        }
                        ld.owned_virt.push_back(std::move(vn));
                    }
                };

                // O(N) linear merge of already-sorted `ld.owned_real` and `vf` (zero `std::sort` when `rem_j.empty()`).
                std::size_t v_idx = 0;
                while (r_idx < ld.owned_real.size() && v_idx < vf.size()) {
                    if (windowed_halo_detail::windowed_key_less(
                            ld.owned_real[r_idx].key, vf[v_idx].first)) {
                        ld.sources.push_back(make_real_source(ld.owned_real[r_idx++]));
                    } else {
                        emit_virt_entry(vf[v_idx].first, vf[v_idx].second);
                        ++v_idx;
                    }
                }
                while (r_idx < ld.owned_real.size()) {
                    ld.sources.push_back(make_real_source(ld.owned_real[r_idx++]));
                }
                while (v_idx < vf.size()) {
                    emit_virt_entry(vf[v_idx].first, vf[v_idx].second);
                    ++v_idx;
                }
            } else {
                for (const auto& rn : ld.owned_real) {
                    ld.sources.push_back(make_real_source(rn));
                }
            }

            if (!rem_j.empty()) {
                for (const auto& se : rem_j) {
                    ld.sources.push_back(se);
                }
                std::sort(ld.sources.begin(), ld.sources.end(), [](const auto& a, const auto& b) noexcept {
                    return windowed_halo_detail::windowed_key_less(a.key, b.key);
                });
            }
            if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                ld.source_packed_keys.resize(ld.sources.size());
                for (std::size_t si = 0; si < ld.sources.size(); ++si) {
                    ld.source_packed_keys[si] =
                        windowed_halo_detail::pack_spatial_key_same_level(ld.sources[si].key);
                }
            }
            std::vector<VirtProlongItem>().swap(prolong_levels[sj]);
        }
        {
            std::vector<std::vector<std::pair<Key<3>, Tensor<T>>>>().swap(coarse_leaves_by_level_);
        }
    }

    // =========================================================================
    // Stage 5a Fast Path: Bottom-Up Virtual Restriction (Flat Sorted Vectors)
    // =========================================================================
    std::vector<RestrictAccum> compute_local_leaf_contribs_unpromoted() {
        const long kk = static_cast<long>(k_);
        const long kk3 = kk * kk * kk;

        auto sort_and_merge_accum = [kk3](std::vector<RestrictAccum>& vec) {
            if (vec.size() <= 1) return;
            std::sort(vec.begin(), vec.end(), [](const RestrictAccum& a, const RestrictAccum& b) noexcept {
                return windowed_halo_detail::windowed_key_less(a.key, b.key);
            });
            std::size_t w = 0;
            for (std::size_t r = 1; r < vec.size(); ++r) {
                if (vec[r].key == vec[w].key) {
                    if (windowed_halo_detail::windowed_key_less(
                            vec[r].min_source_key, vec[w].min_source_key)) {
                        vec[w].min_source_key = vec[r].min_source_key;
                    }
                    aligned_axpy(kk3, vec[w].ds.ptr(), vec[r].ds.ptr(), 1.0);
                } else {
                    ++w;
                    if (w != r) {
                        vec[w] = std::move(vec[r]);
                    }
                }
            }
            vec.resize(w + 1);
        };

        std::vector<std::vector<RestrictAccum>> accum_levels(
            static_cast<std::size_t>(n_max_tree_ + 1));

        for (int j = 1; j <= n_max_tree_; ++j) {
            auto& acc_j = accum_levels[static_cast<std::size_t>(j)];
            const auto& ov = levels_[static_cast<std::size_t>(j)].owned_virt;
            acc_j.reserve(ov.size());
            for (const auto& vn : ov) {
                if (vn.ds_corr.size() == kk3 &&
                    windowed_halo_detail::is_cube_nonzero_readonly(vn.ds_corr)) {
                    acc_j.push_back(RestrictAccum{
                        vn.key,
                        vn.coarse_leaf_key,
                        vn.key,
                        (j == n_max_tree_) ? vn.ds_corr : copy(vn.ds_corr)});
                }
            }
        }

        std::vector<RestrictAccum> leaf_contribs;
        const long k_dims[3] = {kk, kk, kk};
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> r_work1;
        std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> r_work2;

        struct ParentRestrictSlot {
            Tensor<T>* dst = nullptr;
            std::uint8_t count = 0;
            std::array<std::uint32_t, 8> children{};
        };

        for (int j = n_max_tree_; j >= 1; --j) {
            auto& acc_j = accum_levels[static_cast<std::size_t>(j)];
            if (acc_j.empty()) continue;
            sort_and_merge_accum(acc_j);

            std::vector<RestrictAccum> new_leaf_parents;
            std::vector<RestrictAccum> new_virt_parents;
            for (const auto& item : acc_j) {
                const Key<3> par = item.key.parent();
                if (par == item.coarse_leaf) {
                    new_leaf_parents.push_back(
                        RestrictAccum{par, par, item.min_source_key, Tensor<T>()});
                } else {
                    new_virt_parents.push_back(
                        RestrictAccum{par, item.coarse_leaf, item.min_source_key, Tensor<T>()});
                }
            }

            auto ensure_parent_slots = [kk](std::vector<RestrictAccum>& target,
                                            std::vector<RestrictAccum>& candidates) {
                if (candidates.empty()) return;
                std::sort(candidates.begin(),
                          candidates.end(),
                          [](const RestrictAccum& a, const RestrictAccum& b) noexcept {
                              return windowed_halo_detail::windowed_key_less(a.key, b.key);
                          });
                std::size_t w = 0;
                for (std::size_t r = 1; r < candidates.size(); ++r) {
                    if (candidates[r].key == candidates[w].key) {
                        if (windowed_halo_detail::windowed_key_less(
                                candidates[r].min_source_key, candidates[w].min_source_key)) {
                            candidates[w].min_source_key = candidates[r].min_source_key;
                        }
                    } else {
                        ++w;
                        if (w != r) candidates[w] = std::move(candidates[r]);
                    }
                }
                candidates.resize(w + 1);

                std::sort(target.begin(),
                          target.end(),
                          [](const RestrictAccum& a, const RestrictAccum& b) noexcept {
                              return windowed_halo_detail::windowed_key_less(a.key, b.key);
                          });

                std::vector<RestrictAccum> merged;
                merged.reserve(target.size() + candidates.size());
                std::size_t it_t = 0;
                std::size_t it_c = 0;
                while (it_t < target.size() && it_c < candidates.size()) {
                    if (windowed_halo_detail::windowed_key_less(target[it_t].key, candidates[it_c].key)) {
                        merged.push_back(std::move(target[it_t++]));
                    } else if (windowed_halo_detail::windowed_key_less(
                                   candidates[it_c].key, target[it_t].key)) {
                        candidates[it_c].ds = Tensor<T>(kk, kk, kk);
                        merged.push_back(std::move(candidates[it_c++]));
                    } else {
                        if (windowed_halo_detail::windowed_key_less(
                                candidates[it_c].min_source_key, target[it_t].min_source_key)) {
                            target[it_t].min_source_key = candidates[it_c].min_source_key;
                        }
                        merged.push_back(std::move(target[it_t++]));
                        ++it_c;
                    }
                }
                while (it_t < target.size()) merged.push_back(std::move(target[it_t++]));
                while (it_c < candidates.size()) {
                    candidates[it_c].ds = Tensor<T>(kk, kk, kk);
                    merged.push_back(std::move(candidates[it_c++]));
                }
                target.swap(merged);
            };

            ensure_parent_slots(leaf_contribs, new_leaf_parents);
            if (j > 1) {
                ensure_parent_slots(accum_levels[static_cast<std::size_t>(j - 1)], new_virt_parents);
            }

            const std::size_t n_leaf_p = leaf_contribs.size();
            const std::size_t n_virt_p =
                (j > 1) ? accum_levels[static_cast<std::size_t>(j - 1)].size() : 0u;
            std::vector<ParentRestrictSlot> pslots(n_leaf_p + n_virt_p);

            for (std::size_t ci = 0; ci < acc_j.size(); ++ci) {
                const RestrictAccum& item = acc_j[ci];
                const Key<3> par = item.key.parent();
                if (par == item.coarse_leaf) {
                    auto it = std::lower_bound(
                        leaf_contribs.begin(), leaf_contribs.end(), par,
                        [](const RestrictAccum& a, const Key<3>& k) noexcept {
                            return windowed_halo_detail::windowed_key_less(a.key, k);
                        });
                    MADNESS_CHECK(it != leaf_contribs.end() && it->key == par);
                    const std::size_t pid = static_cast<std::size_t>(it - leaf_contribs.begin());
                    ParentRestrictSlot& ps = pslots[pid];
                    ps.dst = &(it->ds);
                    MADNESS_CHECK(ps.count < 8u);
                    ps.children[ps.count++] = static_cast<std::uint32_t>(ci);
                } else {
                    auto& acc_par = accum_levels[static_cast<std::size_t>(j - 1)];
                    auto it = std::lower_bound(
                        acc_par.begin(), acc_par.end(), par,
                        [](const RestrictAccum& a, const Key<3>& k) noexcept {
                            return windowed_halo_detail::windowed_key_less(a.key, k);
                        });
                    MADNESS_CHECK(it != acc_par.end() && it->key == par);
                    const std::size_t pid =
                        n_leaf_p + static_cast<std::size_t>(it - acc_par.begin());
                    ParentRestrictSlot& ps = pslots[pid];
                    ps.dst = &(it->ds);
                    MADNESS_CHECK(ps.count < 8u);
                    ps.children[ps.count++] = static_cast<std::uint32_t>(ci);
                }
            }

            std::vector<ParentRestrictSlot> active_pslots;
            active_pslots.reserve(acc_j.size());
            for (const ParentRestrictSlot& ps : pslots) {
                if (ps.count > 0u) active_pslots.push_back(ps);
            }

            windowed_parallel_detail::parallel_for_chunks(
                world_,
                active_pslots.size(),
                windowed_parallel_detail::kHeavyParallelGrain,
                [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
                    if (r_work1[cid].size() != kk3) {
                        r_work1[cid] = Tensor<T>(3L, k_dims, false);
                        r_work2[cid] = Tensor<T>(3L, k_dims, false);
                    }
                    Tensor<T>& w1 = r_work1[cid];
                    Tensor<T>& w2 = r_work2[cid];
                    for (std::size_t pi = b; pi < e_hi; ++pi) {
                        const ParentRestrictSlot& ps = active_pslots[pi];
                        for (std::uint8_t c = 0; c < ps.count; ++c) {
                            const RestrictAccum& item = acc_j[ps.children[c]];
                            windowed_halo_detail::restrict_step_accumulate_3d(
                                kk3, hT_oct_, item.ds, item.key, *ps.dst, w1, w2);
                        }
                    }
                });
        }

        return leaf_contribs;
    }

    void restrict_all_virtual_unpromoted() {
        const long kk = static_cast<long>(k_);
        const long kk3 = kk * kk * kk;
        auto leaf_contribs = compute_local_leaf_contribs_unpromoted();
        MADNESS_CHECK(world_.size() == 1);
        for (const auto& item : leaf_contribs) {
            const int blvl = static_cast<int>(item.key.level());
            const std::size_t ridx = find_owned_real_index(blvl, item.key);
            MADNESS_CHECK(ridx != SIZE_MAX);
            T* dst = levels_[static_cast<std::size_t>(blvl)].owned_real[ridx].s_out.ptr();
            aligned_axpy(kk3, dst, item.ds.ptr(), 1.0);
        }
    }

    // =========================================================================
    // Stage 5a General Path: Restriction + Dynamic Virtual Subtree Promotion
    // =========================================================================
    void apply_promoted_subtree_items(
        std::vector<windowed_halo_detail::VirtualSubtreeItemMsg<T>>& local_items,
        typename FunctionImpl<T, 3>::dcT& res_coeffs,
        const TensorArgs& full_targs,
        std::vector<std::vector<windowed_halo_detail::OutputNodeInsertMsg<T>>>* remote_outbox) {
        const long kk = static_cast<long>(k_);
        const long kk3 = kk * kk * kk;
        const long twok = 2L * kk;
        const Slice s0(0, kk - 1);
        const ProcessID me = world_.rank();

        std::sort(local_items.begin(),
                  local_items.end(),
                  [](const auto& a, const auto& b) noexcept {
                      if (a.coarse_leaf_key != b.coarse_leaf_key) {
                          return windowed_halo_detail::windowed_key_less(
                              a.coarse_leaf_key, b.coarse_leaf_key);
                      }
                      return windowed_halo_detail::windowed_key_less(a.virt_key, b.virt_key);
                  });

        struct LeafGroup {
            Key<3> bkey;
            std::size_t begin = 0;
            std::size_t end = 0;
            int blvl = 0;
            std::size_t ridx = SIZE_MAX;
        };
        std::vector<LeafGroup> unpromoted_groups;
        std::vector<LeafGroup> promoted_groups;
        {
            std::size_t g_begin = 0;
            while (g_begin < local_items.size()) {
                const Key<3> bkey = local_items[g_begin].coarse_leaf_key;
                std::size_t g_end = g_begin + 1;
                while (g_end < local_items.size() &&
                       local_items[g_end].coarse_leaf_key == bkey) {
                    ++g_end;
                }
                bool leaf_has_promoted = false;
                for (std::size_t idx = g_begin; idx < g_end; ++idx) {
                    if (local_items[idx].is_promoted != 0 &&
                        static_cast<int>(local_items[idx].virt_key.level()) < max_refine_level_) {
                        leaf_has_promoted = true;
                        break;
                    }
                }
                const int blvl = static_cast<int>(bkey.level());
                const std::size_t ridx = find_owned_real_index(blvl, bkey);
                MADNESS_CHECK(ridx != SIZE_MAX);
                const auto& coarse_rn = levels_[static_cast<std::size_t>(blvl)].owned_real[ridx];
                if (coarse_rn.promoted_leaf && blvl < max_refine_level_) {
                    leaf_has_promoted = true;
                }
                const LeafGroup grp{bkey, g_begin, g_end, blvl, ridx};
                if (leaf_has_promoted) {
                    promoted_groups.push_back(grp);
                } else {
                    unpromoted_groups.push_back(grp);
                }
                g_begin = g_end;
            }
        }

        if (!unpromoted_groups.empty()) {
            const long k_dims[3] = {kk, kk, kk};
            std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> u_work1;
            std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> u_work2;
            windowed_parallel_detail::parallel_for_chunks(
                world_,
                unpromoted_groups.size(),
                8,
                [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
                    if (u_work1[cid].size() != kk3) {
                        u_work1[cid] = Tensor<T>(3L, k_dims, false);
                        u_work2[cid] = Tensor<T>(3L, k_dims, false);
                    }
                    Tensor<T>& w1 = u_work1[cid];
                    Tensor<T>& w2 = u_work2[cid];
                    std::vector<std::pair<Key<3>, Tensor<T>>> cur_lvl;
                    std::vector<std::pair<Key<3>, Tensor<T>>> par_lvl;
                    std::vector<Key<3>> par_cands;

                    for (std::size_t gi = b; gi < e_hi; ++gi) {
                        const LeafGroup& grp = unpromoted_groups[gi];
                        auto& coarse_rn =
                            levels_[static_cast<std::size_t>(grp.blvl)].owned_real[grp.ridx];

                        int max_lvl = grp.blvl;
                        for (std::size_t idx = grp.begin; idx < grp.end; ++idx) {
                            if (local_items[idx].ds_corr.size() == kk3) {
                                max_lvl = std::max(
                                    max_lvl, static_cast<int>(local_items[idx].virt_key.level()));
                            }
                        }
                        if (max_lvl <= grp.blvl) continue;

                        cur_lvl.clear();
                        for (std::size_t idx = grp.begin; idx < grp.end; ++idx) {
                            const auto& item = local_items[idx];
                            if (item.ds_corr.size() == kk3 &&
                                static_cast<int>(item.virt_key.level()) == max_lvl) {
                                cur_lvl.emplace_back(item.virt_key, copy(item.ds_corr));
                            }
                        }

                        for (int lvl = max_lvl; lvl > grp.blvl; --lvl) {
                            par_cands.clear();
                            for (const auto& kv : cur_lvl) {
                                const Key<3> par = kv.first.parent();
                                if (par != grp.bkey) par_cands.push_back(par);
                            }
                            windowed_halo_detail::sort_unique_keys(par_cands);

                            par_lvl.clear();
                            if (lvl - 1 > grp.blvl) {
                                std::vector<std::pair<Key<3>, Tensor<T>>> direct_par;
                                for (std::size_t idx = grp.begin; idx < grp.end; ++idx) {
                                    const auto& item = local_items[idx];
                                    if (item.ds_corr.size() == kk3 &&
                                        static_cast<int>(item.virt_key.level()) == lvl - 1) {
                                        direct_par.emplace_back(item.virt_key, copy(item.ds_corr));
                                    }
                                }
                                par_lvl.reserve(direct_par.size() + par_cands.size());
                                std::size_t id = 0;
                                std::size_t ic = 0;
                                while (id < direct_par.size() && ic < par_cands.size()) {
                                    if (windowed_halo_detail::windowed_key_less(
                                            direct_par[id].first, par_cands[ic])) {
                                        par_lvl.push_back(std::move(direct_par[id++]));
                                    } else if (windowed_halo_detail::windowed_key_less(
                                                   par_cands[ic], direct_par[id].first)) {
                                        par_lvl.emplace_back(par_cands[ic++], Tensor<T>(kk, kk, kk));
                                    } else {
                                        par_lvl.push_back(std::move(direct_par[id++]));
                                        ++ic;
                                    }
                                }
                                while (id < direct_par.size()) {
                                    par_lvl.push_back(std::move(direct_par[id++]));
                                }
                                while (ic < par_cands.size()) {
                                    par_lvl.emplace_back(par_cands[ic++], Tensor<T>(kk, kk, kk));
                                }
                            }

                            for (const auto& kv : cur_lvl) {
                                const Key<3>& vk = kv.first;
                                const Key<3> par = vk.parent();
                                if (par == grp.bkey) {
                                    windowed_halo_detail::restrict_step_accumulate_3d(
                                        kk3, hT_oct_, kv.second, vk, coarse_rn.s_out, w1, w2);
                                } else {
                                    auto pit = std::lower_bound(
                                        par_lvl.begin(),
                                        par_lvl.end(),
                                        par,
                                        [](const auto& a, const Key<3>& k) noexcept {
                                            return windowed_halo_detail::windowed_key_less(
                                                a.first, k);
                                        });
                                    MADNESS_CHECK(pit != par_lvl.end() && pit->first == par);
                                    windowed_halo_detail::restrict_step_accumulate_3d(
                                        kk3, hT_oct_, kv.second, vk, pit->second, w1, w2);
                                }
                            }
                            cur_lvl.swap(par_lvl);
                        }
                    }
                });
        }

        Tensor<T> work1(kk, kk, kk);
        Tensor<T> work2(kk, kk, kk);

        auto emit_node = [&](const Key<3>& key, bool has_children, const Tensor<T>& coeff_tensor) {
            const ProcessID owner = (world_.size() == 1) ? me : pmap_->owner(key);
            if (owner == me || remote_outbox == nullptr) {
                GenTensor<T> gt(coeff_tensor, full_targs);
                res_coeffs.replace(key, FunctionNode<T, 3>(gt, has_children));
            } else {
                (*remote_outbox)[static_cast<std::size_t>(owner)].push_back(
                    windowed_halo_detail::OutputNodeInsertMsg<T>{
                        key, static_cast<std::uint8_t>(has_children ? 1u : 0u), coeff_tensor, slot_});
            }
        };

        for (const LeafGroup& grp : promoted_groups) {
            const Key<3>& bkey = grp.bkey;
            const int blvl = grp.blvl;
            auto& coarse_rn = levels_[static_cast<std::size_t>(blvl)].owned_real[grp.ridx];

            std::unordered_set<Key<3>, WindowedKeyHash> split_internals;
            if (coarse_rn.promoted_leaf && blvl < max_refine_level_) {
                split_internals.insert(bkey);
            }
            std::unordered_map<Key<3>, Tensor<T>, WindowedKeyHash> promoted_ns_map;
            std::unordered_map<Key<3>, Tensor<T>, WindowedKeyHash> ds_map;
            int max_v_lvl = blvl;

            for (std::size_t idx = grp.begin; idx < grp.end; ++idx) {
                const auto& item = local_items[idx];
                const int vlvl = static_cast<int>(item.virt_key.level());
                max_v_lvl = std::max(max_v_lvl, vlvl);
                if (item.ds_corr.size() == kk3) {
                    ds_map.emplace(item.virt_key, copy(item.ds_corr));
                }
                if (item.d_virt_ns.size() == twok * twok * twok) {
                    promoted_ns_map.emplace(item.virt_key, item.d_virt_ns);
                }

                if (item.is_promoted != 0 && vlvl < max_refine_level_) {
                    Key<3> anc = item.virt_key;
                    while (anc.level() >= bkey.level()) {
                        split_internals.insert(anc);
                        if (anc == bkey) break;
                        anc = anc.parent();
                    }
                }
            }

            std::unordered_set<Key<3>, WindowedKeyHash> split_leaves;
            for (const Key<3>& ikey : split_internals) {
                for (KeyChildIterator<3> kit(ikey); kit; ++kit) {
                    const Key<3>& ch = kit.key();
                    if (split_internals.count(ch) == 0) {
                        split_leaves.insert(ch);
                    }
                }
            }

            for (int lvl = max_v_lvl; lvl > blvl; --lvl) {
                std::vector<Key<3>> deeper_keys;
                for (const auto& kv : ds_map) {
                    const Key<3>& vk = kv.first;
                    if (static_cast<int>(vk.level()) == lvl &&
                        split_internals.count(vk) == 0 &&
                        split_leaves.count(vk) == 0) {
                        deeper_keys.push_back(vk);
                    }
                }
                std::sort(deeper_keys.begin(), deeper_keys.end(), windowed_halo_detail::windowed_key_less);
                for (const Key<3>& vk : deeper_keys) {
                    const Key<3> par = vk.parent();
                    auto [pit, ins] = ds_map.emplace(par, Tensor<T>(kk, kk, kk));
                    (void)ins;
                    windowed_halo_detail::restrict_step_accumulate_3d(
                        kk3, hT_oct_, ds_map.at(vk), vk, pit->second, work1, work2);
                }
            }

            for (const Key<3>& ikey : split_internals) {
                Tensor<T> ns_box(twok, twok, twok);

                auto p_it = promoted_ns_map.find(ikey);
                if (p_it != promoted_ns_map.end() && p_it->second.size() == twok * twok * twok) {
                    ns_box = copy(p_it->second);
                }

                if (ikey == bkey) {
                    if (coarse_rn.ns_out.size() == twok * twok * twok) {
                        ns_box = copy(coarse_rn.ns_out);
                    }
                    ns_box(s0, s0, s0) = coarse_rn.s_out;
                } else {
                    auto d_it = ds_map.find(ikey);
                    if (d_it != ds_map.end()) {
                        ns_box(s0, s0, s0) = d_it->second;
                    } else {
                        ns_box(s0, s0, s0) = T(0);
                    }
                }
                emit_node(ikey, true, ns_box);
            }

            for (const Key<3>& lkey : split_leaves) {
                Tensor<T> s_box;
                auto d_it = ds_map.find(lkey);
                if (d_it != ds_map.end()) {
                    s_box = copy(d_it->second);
                } else {
                    s_box = Tensor<T>(kk, kk, kk);
                }
                emit_node(lkey, false, s_box);
            }

            split_coarse_leaves_emitted_.push_back(bkey);
        }

        windowed_halo_detail::sort_unique_keys(split_coarse_leaves_emitted_);
    }

    World& world_;
    int k_ = 0;
    int n_max_tree_ = 0;
    bool dynamic_refine_ = true;
    int max_refine_level_ = 30;
    bool ungraded_leaf_refine_only_ = false;
    WindowedSerialInputRange serial_range_{};
    const FunctionCommonData<T, 3>& cdata_;
    array_of_bools<3> is_periodic_;
    std::shared_ptr<WorldDCPmapInterface<Key<3>>> pmap_;
    Tensor<double> h_oct_[8][3];
    Tensor<double> hT_oct_[8][3];

    std::vector<WindowedLevelData<T>> levels_;
    std::vector<WindowedGlobalLevelStats> global_stats_;
    std::vector<WindowedLevelRadii> radii_;
    std::vector<std::vector<std::pair<Key<3>, Tensor<T>>>> coarse_leaves_by_level_;
    std::vector<Key<3>> split_coarse_leaves_emitted_;
    std::unique_ptr<windowed_halo_detail::WindowedNSHaloExchange<T>> owned_exchange_;
    windowed_halo_detail::WindowedNSHaloExchange<T>* exchange_ptr_ = nullptr;
    std::uint16_t slot_ = 0;
};

} // namespace madness

#endif // MADNESS_MRA_WINDOWED_NS_HALO_H__INCLUDED
