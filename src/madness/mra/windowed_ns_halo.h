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
#include <madness/world/thread.h>
#include <madness/world/world_object.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <functional>
#include <mutex>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace madness {

/// Thread-local flag set to `true` when executing inside a function-level batched
/// `ThreadPool` worker task (e.g., vector `apply`), preventing nested outer batch launches.
inline thread_local bool g_in_windowed_batch_task = false;

/// Count of outer function-batch worker threads actively executing `apply_windowed_ns_3d_prepared`
/// (`0` outside batched vector `apply`). Used by `parallel_for_chunks` to recruit idle pool threads
/// automatically when `N < total_threads` or when tail orbitals remain after peers finish.
inline std::atomic<int> g_active_windowed_batch_workers{0};

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

/// RAII guard incrementing/decrementing `g_active_windowed_batch_workers` around active batch execution.
struct WindowedActiveBatchWorkerGuard {
    WindowedActiveBatchWorkerGuard() noexcept {
        g_active_windowed_batch_workers.fetch_add(1, std::memory_order_acq_rel);
    }

    ~WindowedActiveBatchWorkerGuard() {
        g_active_windowed_batch_workers.fetch_sub(1, std::memory_order_acq_rel);
    }

    WindowedActiveBatchWorkerGuard(const WindowedActiveBatchWorkerGuard&) = delete;
    WindowedActiveBatchWorkerGuard& operator=(const WindowedActiveBatchWorkerGuard&) = delete;
};

namespace windowed_parallel_detail {

inline constexpr std::size_t kMaxParallelChunks = 16;
inline constexpr std::size_t kDefaultParallelGrain = 4;
inline std::atomic<int> g_active_inner_helper_tasks{0};
inline thread_local bool g_in_inner_sumfact_chunk = false;

struct InnerChunkFlagGuard {
    bool prev;

    InnerChunkFlagGuard() noexcept
        : prev(g_in_inner_sumfact_chunk) {
        g_in_inner_sumfact_chunk = true;
    }

    ~InnerChunkFlagGuard() {
        g_in_inner_sumfact_chunk = prev;
    }

    InnerChunkFlagGuard(const InnerChunkFlagGuard&) = delete;
    InnerChunkFlagGuard& operator=(const InnerChunkFlagGuard&) = delete;
};

struct ReservedInnerHelpersGuard {
    int n = 0;

    explicit ReservedInnerHelpersGuard(int count) noexcept
        : n(count) {}

    ~ReservedInnerHelpersGuard() {
        if (n > 0) {
            g_active_inner_helper_tasks.fetch_sub(n, std::memory_order_acq_rel);
        }
    }

    ReservedInnerHelpersGuard(const ReservedInnerHelpersGuard&) = delete;
    ReservedInnerHelpersGuard& operator=(const ReservedInnerHelpersGuard&) = delete;
};

template <typename Fn>
class SumFactChunkTask final : public PoolTaskInterface {
public:
    SumFactChunkTask(const Fn* fn,
                     std::size_t begin,
                     std::size_t end,
                     std::size_t chunk_id,
                     std::atomic<int>* remaining,
                     std::atomic<bool>* failed)
        : PoolTaskInterface(TaskAttributes::hipri())
        , fn_(fn)
        , begin_(begin)
        , end_(end)
        , chunk_id_(chunk_id)
        , remaining_(remaining)
        , failed_(failed) {}

    ~SumFactChunkTask() override {
        if (remaining_ != nullptr) {
            remaining_->fetch_sub(1, std::memory_order_release);
        }
    }

    void run(const TaskThreadEnv&) override {
        InnerChunkFlagGuard inner_guard;
        try {
            (*fn_)(begin_, end_, chunk_id_);
        } catch (...) {
            if (failed_ != nullptr) {
                failed_->store(true, std::memory_order_release);
            }
        }
    }

private:
    void get_id(std::pair<void*, unsigned short>& id) const override {
        PoolTaskInterface::make_id(id, &SumFactChunkTask::run);
    }

    const Fn* fn_;
    std::size_t begin_;
    std::size_t end_;
    std::size_t chunk_id_;
    std::atomic<int>* remaining_;
    std::atomic<bool>* failed_;
};

/// Idle-thread-aware deterministic chunk executor over `[0, n_items)`:
/// recruits up to `avail = (pool_sz + 1) - max(1, g_active_windowed_batch_workers) - g_active_inner_helper_tasks`
/// high-priority `ThreadPool` helper tasks when `n_items >= 2 * min_grain` and `!g_in_inner_sumfact_chunk`,
/// saturating all cores across single-function calls, small batches (`N < nthreads`), and tail heavy orbitals
/// with zero oversubscription and exact serial floating-point rounding order per item.
template <typename Fn>
inline void parallel_for_chunks(std::size_t n_items, std::size_t min_grain, const Fn& fn) {
    if (n_items == 0) return;
    const std::size_t eff_grain = std::max<std::size_t>(min_grain, 1);
    const std::size_t pool_sz = ThreadPool::size();
    if (g_in_inner_sumfact_chunk || pool_sz == 0 || n_items < 2 * eff_grain) {
        fn(0, n_items, 0);
        return;
    }

    const std::size_t max_by_grain = n_items / eff_grain;
    if (max_by_grain <= 1) {
        fn(0, n_items, 0);
        return;
    }

    const int max_wanted_helpers = static_cast<int>(
        std::min<std::size_t>(max_by_grain - 1, kMaxParallelChunks - 1));
    if (max_wanted_helpers <= 0) {
        fn(0, n_items, 0);
        return;
    }

    const int total_threads = static_cast<int>(pool_sz) + 1;
    int reserved_helpers = 0;
    int cur_inner = g_active_inner_helper_tasks.load(std::memory_order_relaxed);
    for (;;) {
        const int busy_outer =
            std::max(1, g_active_windowed_batch_workers.load(std::memory_order_relaxed));
        const int avail = total_threads - busy_outer - cur_inner;
        if (avail <= 0) {
            reserved_helpers = 0;
            break;
        }
        const int take = std::min(max_wanted_helpers, avail);
        if (g_active_inner_helper_tasks.compare_exchange_weak(
                cur_inner, cur_inner + take, std::memory_order_acq_rel, std::memory_order_relaxed)) {
            reserved_helpers = take;
            break;
        }
    }

    if (reserved_helpers <= 0) {
        fn(0, n_items, 0);
        return;
    }

    ReservedInnerHelpersGuard reserve_guard(reserved_helpers);
    const std::size_t n_chunks = static_cast<std::size_t>(reserved_helpers) + 1u;
    const std::size_t base = n_items / n_chunks;
    const std::size_t rem_items = n_items % n_chunks;
    auto chunk_bounds = [base, rem_items](std::size_t cid) noexcept -> std::pair<std::size_t, std::size_t> {
        const std::size_t b = cid * base + std::min<std::size_t>(cid, rem_items);
        const std::size_t e = b + base + (cid < rem_items ? 1u : 0u);
        return {b, e};
    };

    std::atomic<int> remaining{reserved_helpers};
    std::atomic<bool> failed{false};

    for (std::size_t cid = 1; cid < n_chunks; ++cid) {
        const auto [b, e] = chunk_bounds(cid);
        ThreadPool::add(new SumFactChunkTask<Fn>(&fn, b, e, cid, &remaining, &failed));
    }
    thread_purge();

    {
        InnerChunkFlagGuard chunk0_guard;
        const auto [b0, e0] = chunk_bounds(0);
        fn(b0, e0, 0);
    }

    ThreadPool::await([&remaining]() { return remaining.load(std::memory_order_acquire) == 0; }, true);
    MADNESS_CHECK(!failed.load(std::memory_order_acquire));
}

} // namespace windowed_parallel_detail

/// Bitmask flags classifying virtual halo nodes at level `j`.
enum WindowedHaloFlags : std::uint8_t {
    kHaloNone = 0u,
    kHaloH1T = 1u << 0,       ///< Primary virtual halo dH_j^{(1,T)} around real nodes T_j (Stage 4a).
    kHaloH1Gamma = 1u << 1,   ///< Primary virtual halo dH_j^{(1,Gamma)} around T_j^{int} (Pass 4b-1 ds_corr).
    kHaloH1Refine = 1u << 2,  ///< Virtual dynamic-refinement candidate around T_j^{int} (Pass 4b-2).
    kHaloH2Ghost = 1u << 3    ///< Secondary read-only virtual ghost dH_j^{(2)} completing scaling stencils.
};

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

    double sh2_s_thresh = 0.0;       ///< Per-node `||s||_F` threshold needed for shell-2 (`sh == 2`) scaling halo.
    double sh2_4b1_thresh = 0.0;     ///< Per-node `||d_int||_F` threshold needed for shell-2 (`sh == 2`) Pass 4b-1 halo.
    double sh2_refine_thresh = 0.0;  ///< Per-node `||(s, d)||_F` threshold needed for shell-2 (`sh == 2`) Pass 4b-2 refine halo.

    int max_refine_cand_bw() const noexcept {
        return std::max(m_max_4b2_int, m_max_4b2_leaf);
    }
};

/// Globally reduced, rank-independent per-level input coefficient norm bounds.
struct WindowedGlobalLevelStats {
    bool has_internal = false;
    double max_s_real = 0.0;             ///< Max `||s_in||_F` over all real nodes at level `j`.
    double max_s_leaf = 0.0;             ///< Max `||s_in||_F` over real leaves at level `j`.
    double max_coarse_leaf_above = 0.0;  ///< Exact bound `max_{b < j} max_s_leaf[b]` on virtual halo `||s_virt||_F` at level `j`.
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

inline double scalar_abs2(double x) noexcept { return x * x; }

inline double scalar_abs2(const double_complex& z) noexcept {
    const double re = z.real();
    const double im = z.imag();
    return re * re + im * im;
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

/// Compute Frobenius norm of the 7 wavelet octants (`[0..2k)^3 \ [0..k)^3`) of a `(2k)^3` NS tensor.
template <typename T>
inline double wavelet_octants_normf(const Tensor<T>& ns_tensor, int k) {
    const long kk = static_cast<long>(k);
    const Slice s0(0, kk - 1);
    const double full_n = ns_tensor.normf();
    const double s0_n = ns_tensor(s0, s0, s0).normf();
    return std::sqrt(std::max(0.0, full_n * full_n - s0_n * s0_n));
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

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & key & flags;
    }
};

struct KeyQueryMsg {
    Key<3> key;
    std::uint8_t need_internal_ns = 0;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & key & need_internal_ns;
    }
};

template <typename T>
struct KeyReplyMsg {
    Key<3> key;
    std::uint8_t is_internal = 0;
    Tensor<T> coeff;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & key & is_internal & coeff;
    }
};

template <typename T>
struct CoarseLeafDeltaMsg {
    Key<3> leaf_key;
    Key<3> sort_tag_key;
    Tensor<T> delta_s;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & leaf_key & sort_tag_key & delta_s;
    }
};

template <typename T>
struct VirtualSubtreeItemMsg {
    Key<3> coarse_leaf_key;
    Key<3> virt_key;
    std::uint8_t is_promoted = 0;
    Tensor<T> ds_corr;
    Tensor<T> d_virt_ns;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & coarse_leaf_key & virt_key & is_promoted & ds_corr & d_virt_ns;
    }
};

template <typename T>
struct OutputNodeInsertMsg {
    Key<3> key;
    std::uint8_t has_children = 0;
    Tensor<T> coeff;

    template <typename Archive>
    void serialize(Archive& ar) {
        ar & key & has_children & coeff;
    }
};

/// Dedicated `WorldObject` handling batched multi-rank halo discovery, ghost/ancestor
/// fetches, and Stage 5a coarse-leaf restriction / promotion assembly.
template <typename T>
class WindowedNSHaloExchange : public WorldObject<WindowedNSHaloExchange<T>> {
public:
    using baseT = WorldObject<WindowedNSHaloExchange<T>>;

    explicit WindowedNSHaloExchange(World& world)
        : baseT(world) {
        this->process_pending();
    }

    ~WindowedNSHaloExchange() override = default;

    void recv_halo_candidates(const std::vector<HaloCandidateMsg>& msgs) {
        std::lock_guard<std::mutex> lock(mtx_);
        halo_candidates_.insert(halo_candidates_.end(), msgs.begin(), msgs.end());
    }

    void recv_key_queries(ProcessID requester, const std::vector<KeyQueryMsg>& queries) {
        std::lock_guard<std::mutex> lock(mtx_);
        key_queries_.emplace_back(requester, queries);
    }

    void recv_key_replies(const std::vector<KeyReplyMsg<T>>& replies) {
        std::lock_guard<std::mutex> lock(mtx_);
        key_replies_.insert(key_replies_.end(), replies.begin(), replies.end());
    }

    void recv_coarse_deltas(const std::vector<CoarseLeafDeltaMsg<T>>& deltas) {
        std::lock_guard<std::mutex> lock(mtx_);
        coarse_deltas_.insert(coarse_deltas_.end(), deltas.begin(), deltas.end());
    }

    void recv_virtual_subtree_items(const std::vector<VirtualSubtreeItemMsg<T>>& items) {
        std::lock_guard<std::mutex> lock(mtx_);
        subtree_items_.insert(subtree_items_.end(), items.begin(), items.end());
    }

    void recv_output_inserts(const std::vector<OutputNodeInsertMsg<T>>& inserts) {
        std::lock_guard<std::mutex> lock(mtx_);
        output_inserts_.insert(output_inserts_.end(), inserts.begin(), inserts.end());
    }

    std::vector<HaloCandidateMsg> take_halo_candidates() {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<HaloCandidateMsg> out;
        out.swap(halo_candidates_);
        return out;
    }

    std::vector<std::pair<ProcessID, std::vector<KeyQueryMsg>>> take_key_queries() {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<std::pair<ProcessID, std::vector<KeyQueryMsg>>> out;
        out.swap(key_queries_);
        return out;
    }

    std::vector<KeyReplyMsg<T>> take_key_replies() {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<KeyReplyMsg<T>> out;
        out.swap(key_replies_);
        return out;
    }

    std::vector<CoarseLeafDeltaMsg<T>> take_coarse_deltas() {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<CoarseLeafDeltaMsg<T>> out;
        out.swap(coarse_deltas_);
        return out;
    }

    std::vector<VirtualSubtreeItemMsg<T>> take_subtree_items() {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<VirtualSubtreeItemMsg<T>> out;
        out.swap(subtree_items_);
        return out;
    }

    std::vector<OutputNodeInsertMsg<T>> take_output_inserts() {
        std::lock_guard<std::mutex> lock(mtx_);
        std::vector<OutputNodeInsertMsg<T>> out;
        out.swap(output_inserts_);
        return out;
    }

private:
    std::mutex mtx_;
    std::vector<HaloCandidateMsg> halo_candidates_;
    std::vector<std::pair<ProcessID, std::vector<KeyQueryMsg>>> key_queries_;
    std::vector<KeyReplyMsg<T>> key_replies_;
    std::vector<CoarseLeafDeltaMsg<T>> coarse_deltas_;
    std::vector<VirtualSubtreeItemMsg<T>> subtree_items_;
    std::vector<OutputNodeInsertMsg<T>> output_inserts_;
};

} // namespace windowed_halo_detail

/// Owned real octree node record on this MPI rank at level `j`.
template <typename T>
struct WindowedOwnedRealNode {
    Key<3> key;
    bool is_internal = false;        ///< True iff `node.has_children()` in input tree.
    bool near_internal_4b1 = false;  ///< For leaves: within `m_max_4b1 >= 1` of `T_j^{int}`.
    bool near_internal_4b2 = false;  ///< For leaves: within `max_refine_cand_bw >= 1` of `T_j^{int}`.
    bool promoted_leaf = false;      ///< For leaves: true iff dynamically refined in Stage 4(b).
    T s000 = T(0);                   ///< `s_in(0, 0, 0)` constant moment.
    double s_osc_norm2 = 0.0;        ///< `||s_in - s000 e_{000}||_F^2` (cancellation-free).
    double s_normf = 0.0;            ///< `||s_in||_F`.
    double d_normf = 0.0;            ///< `||(0, d_in)||_F` (nonzero iff `is_internal`).
    double ns_normf = 0.0;           ///< `||(s_in, d_in)||_F`.
    Tensor<T> s_in;                  ///< `k^3` input scaling coefficients.
    Tensor<T> ns_in;                 ///< `(2k)^3` input NS `(s, d)` (populated iff `is_internal`).
    Tensor<T> wav_in;                ///< `(2k)^3` pure wavelet `(0, d)` (populated iff `is_internal`).
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
    Tensor<T> wav_coeff;             ///< `(2k)^3` pure wavelet `(0, d)` (iff `is_internal`; else empty).
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

/// Complete Stage 1b + Stage 5a distributed halo workspace for one `apply_windowed_ns_3d` call.
template <typename T>
class WindowedNSHaloWorkspace {
public:
    using RadiiBuilderFn = std::function<void(const std::vector<WindowedGlobalLevelStats>&,
                                              std::vector<WindowedLevelRadii>&)>;

    WindowedNSHaloWorkspace(World& world,
                            const Function<T, 3>& f_ns,
                            int n_max_tree,
                            bool dynamic_refine,
                            int max_refine_level,
                            const RadiiBuilderFn& radii_builder)
        : world_(world)
        , k_(f_ns.k())
        , n_max_tree_(n_max_tree)
        , dynamic_refine_(dynamic_refine)
        , max_refine_level_(max_refine_level)
        , cdata_(FunctionCommonData<T, 3>::get(f_ns.k()))
        , is_periodic_(FunctionDefaults<3>::get_bc().is_periodic())
        , pmap_(f_ns.get_pmap()) {
        MADNESS_CHECK(f_ns.get_impl()->get_tree_state() == nonstandard_with_leaves);
        init_octant_filters();
        levels_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        global_stats_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        radii_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        coarse_leaves_by_level_.resize(static_cast<std::size_t>(n_max_tree_ + 1));
        for (int j = 0; j <= n_max_tree_; ++j) {
            levels_[static_cast<std::size_t>(j)].level = j;
        }

        if (world_.size() > 1) {
            exchange_ = std::make_unique<windowed_halo_detail::WindowedNSHaloExchange<T>>(world_);
            build_distributed(f_ns, radii_builder);
        } else {
            build_serial(f_ns, radii_builder);
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
        if (exchange_) {
            world_.gop.fence();
            exchange_.reset();
        }
    }

    WindowedNSHaloWorkspace(const WindowedNSHaloWorkspace&) = delete;
    WindowedNSHaloWorkspace& operator=(const WindowedNSHaloWorkspace&) = delete;

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

    /// Stage 5a + Stage 5b assembly:
    ///   1. Restrict all non-promoted virtual boundary corrections `\Delta\hat{s}^j_{V,corr}`
    ///      bottom-up (`j = n_max_tree -> 1`) onto their real coarse leaves `B in L`
    ///      (or split subtrees when dynamic output refinement promotes virtual/leaf nodes).
    ///   2. Populate `result` in `TreeState::nonstandard_after_apply` so a subsequent
    ///      `result.reconstruct(fence)` performs top-down NS synthesis cleanly.
    void assemble_output_nonstandard_tree(Function<T, 3>& result) {
        const long kk = static_cast<long>(k_);
        const long twok = 2L * kk;
        const Slice s0(0, kk - 1);
        auto* res_impl = result.get_impl().get();
        auto& res_coeffs = res_impl->get_coeffs();
        const TensorArgs full_targs(-1.0, TT_FULL);

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
        long global_any_promoted = local_any_promoted;
        if (world_.size() > 1) {
            world_.gop.sum(global_any_promoted);
        }

        if (global_any_promoted == 0) {
            restrict_all_virtual_unpromoted();
        } else {
            restrict_and_split_promoted_virtual(res_coeffs, full_targs);
        }

        std::vector<std::vector<windowed_halo_detail::OutputNodeInsertMsg<T>>> remote_child_inserts(
            static_cast<std::size_t>(world_.size()));

        for (int j = 0; j <= n_max_tree_; ++j) {
            auto& ld = levels_[static_cast<std::size_t>(j)];
            for (auto& rn : ld.owned_real) {
                if (!split_coarse_leaves_emitted_.empty() &&
                    std::binary_search(split_coarse_leaves_emitted_.begin(),
                                       split_coarse_leaves_emitted_.end(),
                                       rn.key,
                                       windowed_halo_detail::windowed_key_less)) {
                    continue;
                }

                if (rn.is_internal) {
                    Tensor<T> full_ns = rn.ns_out;
                    if (full_ns.size() != twok * twok * twok) {
                        full_ns = Tensor<T>(twok, twok, twok);
                    }
                    full_ns(s0, s0, s0) = rn.s_out;
                    GenTensor<T> coeff(full_ns, full_targs);
                    res_coeffs.replace(rn.key, FunctionNode<T, 3>(coeff, true));
                } else if (rn.promoted_leaf && j < max_refine_level_) {
                    Tensor<T> full_ns = rn.ns_out;
                    if (full_ns.size() != twok * twok * twok) {
                        full_ns = Tensor<T>(twok, twok, twok);
                    }
                    full_ns(s0, s0, s0) = rn.s_out;
                    GenTensor<T> coeff(full_ns, full_targs);
                    res_coeffs.replace(rn.key, FunctionNode<T, 3>(coeff, true));

                    for (KeyChildIterator<3> kit(rn.key); kit; ++kit) {
                        const Key<3>& child = kit.key();
                        Tensor<T> zero_s(kk, kk, kk);
                        const ProcessID owner = pmap_->owner(child);
                        if (owner == world_.rank()) {
                            GenTensor<T> c_coeff(zero_s, full_targs);
                            res_coeffs.replace(child, FunctionNode<T, 3>(c_coeff, false));
                        } else {
                            remote_child_inserts[static_cast<std::size_t>(owner)].push_back(
                                windowed_halo_detail::OutputNodeInsertMsg<T>{child, 0, std::move(zero_s)});
                        }
                    }
                } else {
                    GenTensor<T> coeff(rn.s_out, full_targs);
                    res_coeffs.replace(rn.key, FunctionNode<T, 3>(coeff, false));
                }
            }
        }

        if (world_.size() > 1) {
            for (ProcessID p = 0; p < world_.size(); ++p) {
                auto& batch = remote_child_inserts[static_cast<std::size_t>(p)];
                if (!batch.empty()) {
                    exchange_->send(p,
                                    &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_output_inserts,
                                    batch);
                }
            }
            world_.gop.fence();
            for (auto& msg : exchange_->take_output_inserts()) {
                GenTensor<T> coeff(msg.coeff, full_targs);
                res_coeffs.replace(msg.key, FunctionNode<T, 3>(coeff, msg.has_children != 0));
            }
            world_.gop.fence();
        }

        res_impl->set_tree_state(nonstandard_after_apply);
    }

private:
    using KeyFlagVec = std::vector<std::pair<Key<3>, std::uint8_t>>;

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

    void ingest_local_real_nodes(const Function<T, 3>& f_ns) {
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
                        rn.wav_in =
                            windowed_halo_detail::make_pure_wavelet_cube(full_c, k_, &rn.d_normf);
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
            }
            levels_[sj].owned_real.push_back(std::move(rn));
        }

        if (world_.size() > 1) {
            const std::size_t n_lvls = static_cast<std::size_t>(n_max_tree_ + 1);
            std::vector<double> gop_buf(5 * n_lvls, 0.0);
            for (std::size_t j = 0; j < n_lvls; ++j) {
                const auto& st = global_stats_[j];
                gop_buf[5 * j + 0] = st.has_internal ? 1.0 : 0.0;
                gop_buf[5 * j + 1] = st.max_s_real;
                gop_buf[5 * j + 2] = st.max_s_leaf;
                gop_buf[5 * j + 3] = st.max_d_int;
                gop_buf[5 * j + 4] = st.max_ns_int;
            }
            world_.gop.max(gop_buf.data(), gop_buf.size());
            for (std::size_t j = 0; j < n_lvls; ++j) {
                auto& st = global_stats_[j];
                st.has_internal = (gop_buf[5 * j + 0] > 0.5);
                st.max_s_real = gop_buf[5 * j + 1];
                st.max_s_leaf = gop_buf[5 * j + 2];
                st.max_d_int = gop_buf[5 * j + 3];
                st.max_ns_int = gop_buf[5 * j + 4];
            }
        }

        double coarse_leaf_max = 0.0;
        for (int j = 0; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            auto& st = global_stats_[sj];
            st.max_coarse_leaf_above = coarse_leaf_max;
            st.max_s_source = std::max(st.max_s_real, coarse_leaf_max);
            coarse_leaf_max = std::max(coarse_leaf_max, st.max_s_leaf);

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

            if (world_.size() > 1) {
                auto& cl_vec = coarse_leaves_by_level_[sj];
                for (const auto& rn : ld.owned_real) {
                    if (!rn.is_internal) {
                        cl_vec.emplace_back(rn.key, rn.s_in);
                    }
                }
            }
        }
    }

    // =========================================================================
    // Serial Single-Rank Build (`world.size() == 1`) — Parallel + 64-Bit Packed Keys
    // =========================================================================
    void build_serial(const Function<T, 3>& f_ns, const RadiiBuilderFn& radii_builder) {
        ingest_local_real_nodes(f_ns);
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
            const bool can_refine_j = dynamic_refine_ && (j < max_refine_level_);
            const bool need_T = (rad.has_T_terms && rad.m_max_T >= 1);
            const int m1 = (ld.global_has_internal && rad.has_4b1_terms) ? rad.m_max_4b1 : -1;
            const int m2 = (ld.global_has_internal && can_refine_j) ? rad.max_refine_cand_bw() : -1;
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
                h1_slices.size(), 1, [&](std::size_t sb, std::size_t se, std::size_t cid) {
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
                                const int mT_rn =
                                    (rad.m_max_T >= 2 &&
                                     (!dynamic_refine_ || rn.s_normf > rad.sh2_s_thresh))
                                        ? rad.m_max_T
                                        : std::min(rad.m_max_T, 1);
                                if (use_packed) {
                                    const auto& l = rn.key.translation();
                                    auto& out_p = my_h1_packed[sj];
                                    for (int dx = -mT_rn; dx <= mT_rn; ++dx) {
                                        for (int dy = -mT_rn; dy <= mT_rn; ++dy) {
                                            for (int dz = -mT_rn; dz <= mT_rn; ++dz) {
                                                if (dx == 0 && dy == 0 && dz == 0) continue;
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
                                                if (dx == 0 && dy == 0 && dz == 0) continue;
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
                            const int m2 = can_refine_j ? rad.max_refine_cand_bw() : -1;
                            if (std::max(m1, m2) >= 1) {
                                auto& out_near = my_near[sj];
                                for (std::size_t ri = sl.begin; ri < sl.end; ++ri) {
                                    const auto& rn = ld.owned_real[ri];
                                    if (!rn.is_internal) continue;
                                    const int m1_rn =
                                        (m1 >= 2 &&
                                         (!dynamic_refine_ || rn.d_normf > rad.sh2_4b1_thresh))
                                            ? m1
                                            : std::min(m1, 1);
                                    const int m2_rn =
                                        (m2 >= 2 &&
                                         (!dynamic_refine_ || rn.ns_normf > rad.sh2_refine_thresh))
                                            ? m2
                                            : std::min(m2, 1);
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
                                                            if (sh <= m2_rn) nbits |= 2u;
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
                                                            if (sh <= m2_rn) nbits |= 2u;
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

        for (int j = 1; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            auto& ld = levels_[sj];
            for (std::size_t cid = 0; cid < windowed_parallel_detail::kMaxParallelChunks; ++cid) {
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
                for (std::size_t cid = 0; cid < windowed_parallel_detail::kMaxParallelChunks; ++cid) {
                    if (!h1_packed_tls[cid].empty()) total_p += h1_packed_tls[cid][sj].size();
                }
                auto& h1p = h1_packed_lvls[sj];
                h1p.reserve(total_p);
                for (std::size_t cid = 0; cid < windowed_parallel_detail::kMaxParallelChunks; ++cid) {
                    if (!h1_packed_tls[cid].empty()) {
                        const auto& src = h1_packed_tls[cid][sj];
                        h1p.insert(h1p.end(), src.begin(), src.end());
                    }
                }
                windowed_halo_detail::sort_and_merge_packed_flags(h1p);
            } else {
                std::size_t total_k = 0;
                for (std::size_t cid = 0; cid < windowed_parallel_detail::kMaxParallelChunks; ++cid) {
                    if (!h1_key_tls[cid].empty()) total_k += h1_key_tls[cid][sj].size();
                }
                auto& h1k = h1_key_lvls[sj];
                h1k.reserve(total_k);
                for (std::size_t cid = 0; cid < windowed_parallel_detail::kMaxParallelChunks; ++cid) {
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

            const int mT_h2 =
                (rad.has_T_terms && rad.m_max_T >= 1)
                    ? ((rad.m_max_T >= 2 && (!dynamic_refine_ || max_cl_above > rad.sh2_s_thresh))
                           ? rad.m_max_T
                           : std::min(rad.m_max_T, 1))
                    : -1;
            const int mLeaf_h2 =
                (ld.global_has_internal && rad.has_4b2_leaf_terms && rad.m_max_4b2_leaf >= 1)
                    ? ((rad.m_max_4b2_leaf >= 2 &&
                        (!dynamic_refine_ || max_cl_above > rad.sh2_s_thresh))
                           ? rad.m_max_4b2_leaf
                           : std::min(rad.m_max_4b2_leaf, 1))
                    : -1;
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
                h2_slices.size(), 1, [&](std::size_t sb, std::size_t se, std::size_t cid) {
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
        for (int j = 1; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            auto& vf = virt_flags[sj];

            if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                const Level lvl = static_cast<Level>(j);
                const auto& h1p = h1_packed_lvls[sj];
                std::size_t total_h2 = 0;
                for (std::size_t cid = 0; cid < windowed_parallel_detail::kMaxParallelChunks; ++cid) {
                    if (!h2_packed_tls[cid].empty()) total_h2 += h2_packed_tls[cid][sj].size();
                }
                std::vector<std::uint64_t> h2p;
                h2p.reserve(total_h2);
                for (std::size_t cid = 0; cid < windowed_parallel_detail::kMaxParallelChunks; ++cid) {
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
                for (std::size_t cid = 0; cid < windowed_parallel_detail::kMaxParallelChunks; ++cid) {
                    if (!h2_key_tls[cid].empty()) total_h2 += h2_key_tls[cid][sj].size();
                }
                std::vector<Key<3>> h2k;
                h2k.reserve(total_h2);
                for (std::size_t cid = 0; cid < windowed_parallel_detail::kMaxParallelChunks; ++cid) {
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
        const ProcessID me = world_.rank();
        const int nproc = world_.size();

        ingest_local_real_nodes(f_ns);
        radii_builder(global_stats_, radii_);

        // 2. Collective Exchange 1: Discover owned H1 virtual halo nodes & near_internal real leaves.
        std::vector<KeyFlagVec> out_cand_vecs(static_cast<std::size_t>(nproc));

        auto emit_candidate = [&](const Key<3>& nb, std::uint8_t flag_bits) {
            const ProcessID owner = pmap_->owner(nb);
            out_cand_vecs[static_cast<std::size_t>(owner)].emplace_back(nb, flag_bits);
        };

        for (int j = 1; j <= n_max_tree_; ++j) {
            const auto& ld = levels_[static_cast<std::size_t>(j)];
            const WindowedLevelRadii& rad = radii_[static_cast<std::size_t>(j)];
            const bool can_refine_j = dynamic_refine_ && (j < max_refine_level_);

            if (rad.has_T_terms && rad.m_max_T >= 1) {
                for (const auto& rn : ld.owned_real) {
                    const int mT_rn =
                        (rad.m_max_T >= 2 && (!dynamic_refine_ || rn.s_normf > rad.sh2_s_thresh))
                            ? rad.m_max_T
                            : std::min(rad.m_max_T, 1);
                    for (int dx = -mT_rn; dx <= mT_rn; ++dx) {
                        for (int dy = -mT_rn; dy <= mT_rn; ++dy) {
                            for (int dz = -mT_rn; dz <= mT_rn; ++dz) {
                                if (dx == 0 && dy == 0 && dz == 0) continue;
                                Key<3> nb;
                                if (windowed_halo_detail::shift_key_in_domain(
                                        rn.key, dx, dy, dz, is_periodic_, nb)) {
                                    emit_candidate(nb, kHaloH1T);
                                }
                            }
                        }
                    }
                }
            }

            if (ld.global_has_internal) {
                const int m1 = rad.has_4b1_terms ? rad.m_max_4b1 : -1;
                const int m2 = can_refine_j ? rad.max_refine_cand_bw() : -1;
                if (std::max(m1, m2) >= 1) {
                    for (const auto& rn : ld.owned_real) {
                        if (!rn.is_internal) continue;
                        const int m1_rn =
                            (m1 >= 2 && (!dynamic_refine_ || rn.d_normf > rad.sh2_4b1_thresh))
                                ? m1
                                : std::min(m1, 1);
                        const int m2_rn =
                            (m2 >= 2 && (!dynamic_refine_ || rn.ns_normf > rad.sh2_refine_thresh))
                                ? m2
                                : std::min(m2, 1);
                        const int m_emit_rn = std::max(m1_rn, m2_rn);
                        if (m_emit_rn < 1) continue;
                        for (int dx = -m_emit_rn; dx <= m_emit_rn; ++dx) {
                            for (int dy = -m_emit_rn; dy <= m_emit_rn; ++dy) {
                                for (int dz = -m_emit_rn; dz <= m_emit_rn; ++dz) {
                                    if (dx == 0 && dy == 0 && dz == 0) continue;
                                    const int sh = std::max({std::abs(dx), std::abs(dy), std::abs(dz)});
                                    Key<3> nb;
                                    if (windowed_halo_detail::shift_key_in_domain(
                                            rn.key, dx, dy, dz, is_periodic_, nb)) {
                                        std::uint8_t bits = 0u;
                                        if (sh <= m1_rn) bits = static_cast<std::uint8_t>(bits | kHaloH1Gamma);
                                        if (sh <= m2_rn) bits = static_cast<std::uint8_t>(bits | kHaloH1Refine);
                                        if (bits != 0u) emit_candidate(nb, bits);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for (int p = 0; p < nproc; ++p) {
            auto& vec = out_cand_vecs[static_cast<std::size_t>(p)];
            if (vec.empty()) continue;
            windowed_halo_detail::sort_and_merge_key_flags(vec);
            std::vector<windowed_halo_detail::HaloCandidateMsg> batch;
            batch.reserve(vec.size());
            for (const auto& kv : vec) {
                batch.push_back(windowed_halo_detail::HaloCandidateMsg{kv.first, kv.second});
            }
            if (p == me) {
                exchange_->recv_halo_candidates(batch);
            } else {
                exchange_->send(static_cast<ProcessID>(p),
                                &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_halo_candidates,
                                batch);
            }
        }
        world_.gop.fence();

        std::vector<KeyFlagVec> virt_flags(static_cast<std::size_t>(n_max_tree_ + 1));
        for (const auto& cand : exchange_->take_halo_candidates()) {
            const int j = static_cast<int>(cand.key.level());
            const std::size_t ridx = find_owned_real_index(j, cand.key);
            if (ridx != SIZE_MAX) {
                auto& rn = levels_[static_cast<std::size_t>(j)].owned_real[ridx];
                if (!rn.is_internal) {
                    if ((cand.flags & kHaloH1Gamma) != 0) rn.near_internal_4b1 = true;
                    if ((cand.flags & kHaloH1Refine) != 0) rn.near_internal_4b2 = true;
                }
            } else {
                virt_flags[static_cast<std::size_t>(j)].emplace_back(cand.key, cand.flags);
            }
        }
        for (int j = 1; j <= n_max_tree_; ++j) {
            windowed_halo_detail::sort_and_merge_key_flags(virt_flags[static_cast<std::size_t>(j)]);
        }

        // 3. Determine needed scaling+virtual stencil sources (`needed_virt_capable`) vs
        //    strictly real-node stencil sources (`needed_real_only`) at each level j >= 1.
        std::vector<std::vector<Key<3>>> needed_virt_capable(static_cast<std::size_t>(n_max_tree_ + 1));
        std::vector<std::vector<Key<3>>> needed_real_only(static_cast<std::size_t>(n_max_tree_ + 1));

        auto add_ring_to = [&](std::vector<Key<3>>& dst, const Key<3>& center, int radius) {
            if (radius < 0) return;
            for (int dx = -radius; dx <= radius; ++dx) {
                for (int dy = -radius; dy <= radius; ++dy) {
                    for (int dz = -radius; dz <= radius; ++dz) {
                        Key<3> nb;
                        if (windowed_halo_detail::shift_key_in_domain(
                                center, dx, dy, dz, is_periodic_, nb)) {
                            dst.push_back(nb);
                        }
                    }
                }
            }
        };

        for (int j = 1; j <= n_max_tree_; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const auto& ld = levels_[sj];
            const WindowedLevelRadii& rad = radii_[sj];
            const auto& vf = virt_flags[sj];
            const bool can_refine_j = dynamic_refine_ && (j < max_refine_level_);
            const double max_cl_above = global_stats_[sj].max_coarse_leaf_above;
            auto& vc_vec = needed_virt_capable[sj];
            auto& ro_vec = needed_real_only[sj];

            if (rad.has_T_terms && rad.m_max_T >= 1) {
                const int mT_h2 =
                    (rad.m_max_T >= 2 && (!dynamic_refine_ || max_cl_above > rad.sh2_s_thresh))
                        ? rad.m_max_T
                        : std::min(rad.m_max_T, 1);
                for (const auto& rn : ld.owned_real) {
                    const int mT_rn =
                        (rad.m_max_T >= 2 && (!dynamic_refine_ || rn.s_normf > rad.sh2_s_thresh))
                            ? rad.m_max_T
                            : std::min(rad.m_max_T, 1);
                    add_ring_to(vc_vec, rn.key, mT_rn);
                    if (rad.m_max_T > mT_rn) {
                        add_ring_to(ro_vec, rn.key, rad.m_max_T);
                    }
                }
                for (const auto& kv : vf) {
                    if ((kv.second & kHaloH1T) != 0) {
                        add_ring_to(vc_vec, kv.first, mT_h2);
                        if (rad.m_max_T > mT_h2) {
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
                        (mLeaf >= 2 && (!dynamic_refine_ || max_cl_above > rad.sh2_s_thresh))
                            ? mLeaf
                            : std::min(mLeaf, 1);
                    for (const auto& rn : ld.owned_real) {
                        if (rn.is_internal || (can_refine_j && rn.near_internal_4b2)) {
                            add_ring_to(vc_vec, rn.key, mLeaf_h2);
                            if (mLeaf > mLeaf_h2) {
                                add_ring_to(ro_vec, rn.key, mLeaf);
                            }
                        }
                    }
                    if (can_refine_j) {
                        for (const auto& kv : vf) {
                            if ((kv.second & kHaloH1Refine) != 0) {
                                add_ring_to(vc_vec, kv.first, mLeaf_h2);
                                if (mLeaf > mLeaf_h2) {
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

        std::vector<KeyFlagVec> query_vecs(static_cast<std::size_t>(nproc));

        auto register_query_chain = [&](const Key<3>& start_key) {
            const int j_start = static_cast<int>(start_key.level());
            if (pmap_->owner(start_key) == me && find_owned_real_index(j_start, start_key) != SIZE_MAX) {
                return;
            }
            if (pmap_->owner(start_key) != me) {
                query_vecs[static_cast<std::size_t>(pmap_->owner(start_key))].emplace_back(start_key, 1u);
            }

            Key<3> curr = start_key;
            while (curr.level() > 0) {
                curr = curr.parent();
                const int clvl = static_cast<int>(curr.level());
                const ProcessID owner = pmap_->owner(curr);
                if (owner == me) {
                    if (find_owned_real_index(clvl, curr) != SIZE_MAX) {
                        break;
                    }
                } else {
                    query_vecs[static_cast<std::size_t>(owner)].emplace_back(curr, 0u);
                }
            }
        };

        for (int j = 1; j <= n_max_tree_; ++j) {
            for (const Key<3>& skey : needed_virt_capable[static_cast<std::size_t>(j)]) {
                register_query_chain(skey);
            }
            for (const auto& kv : virt_flags[static_cast<std::size_t>(j)]) {
                register_query_chain(kv.first);
            }
            for (const Key<3>& rkey : needed_real_only[static_cast<std::size_t>(j)]) {
                const ProcessID owner = pmap_->owner(rkey);
                if (owner != me) {
                    query_vecs[static_cast<std::size_t>(owner)].emplace_back(rkey, 1u);
                }
            }
        }

        for (int p = 0; p < nproc; ++p) {
            auto& qv = query_vecs[static_cast<std::size_t>(p)];
            if (qv.empty()) continue;
            windowed_halo_detail::sort_and_merge_key_flags(qv);
            std::vector<windowed_halo_detail::KeyQueryMsg> batch;
            batch.reserve(qv.size());
            for (const auto& kv : qv) {
                batch.push_back(windowed_halo_detail::KeyQueryMsg{kv.first, kv.second});
            }
            exchange_->send(static_cast<ProcessID>(p),
                            &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_key_queries,
                            me,
                            batch);
        }
        world_.gop.fence();

        for (const auto& [requester, queries] : exchange_->take_key_queries()) {
            std::vector<windowed_halo_detail::KeyReplyMsg<T>> replies;
            replies.reserve(queries.size());
            for (const auto& q : queries) {
                const int qlvl = static_cast<int>(q.key.level());
                if (qlvl < 0 || qlvl > n_max_tree_) continue;
                const std::size_t ridx = find_owned_real_index(qlvl, q.key);
                if (ridx == SIZE_MAX) continue;

                const auto& rn = levels_[static_cast<std::size_t>(qlvl)].owned_real[ridx];
                if (!rn.is_internal) {
                    replies.push_back(windowed_halo_detail::KeyReplyMsg<T>{rn.key, 0u, rn.s_in});
                } else if (q.need_internal_ns != 0) {
                    replies.push_back(windowed_halo_detail::KeyReplyMsg<T>{rn.key, 1u, rn.ns_in});
                }
            }
            if (!replies.empty()) {
                exchange_->send(requester,
                                &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_key_replies,
                                replies);
            }
        }
        world_.gop.fence();

        // Ingest remote real nodes & remote coarse leaves.
        std::vector<WindowedSourceEntry<T>> remote_real_same_level;
        for (auto& rep : exchange_->take_key_replies()) {
            const int rlvl = static_cast<int>(rep.key.level());
            if (rep.is_internal == 0) {
                if (rlvl >= 0 && rlvl <= n_max_tree_) {
                    coarse_leaves_by_level_[static_cast<std::size_t>(rlvl)].emplace_back(
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
                se.s_coeff = windowed_halo_detail::extract_scaling_cube(rep.coeff, k_);
                windowed_halo_detail::compute_scaling_norms_readonly(
                    se.s_coeff, se.s000, se.s_osc_norm2, se.s_normf);
                se.ns_coeff = rep.coeff;
                se.wav_coeff =
                    windowed_halo_detail::make_pure_wavelet_cube(rep.coeff, k_, &se.d_normf);
                se.ns_normf = std::sqrt(se.s_normf * se.s_normf + se.d_normf * se.d_normf);
                remote_real_same_level.push_back(std::move(se));
            }
        }

        for (int j = 0; j <= n_max_tree_; ++j) {
            auto& cl = coarse_leaves_by_level_[static_cast<std::size_t>(j)];
            std::sort(cl.begin(), cl.end(), [](const auto& a, const auto& b) noexcept {
                return windowed_halo_detail::windowed_key_less(a.first, b.first);
            });
            cl.erase(std::unique(cl.begin(), cl.end(), [](const auto& a, const auto& b) noexcept {
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

        // Mark any needed scaling/virtual source at level j that is neither local real nor remote real as `kHaloH2Ghost`.
        for (int j = 1; j <= n_max_tree_; ++j) {
            auto& vf = virt_flags[static_cast<std::size_t>(j)];
            std::vector<Key<3>> h2_add;
            for (const Key<3>& skey : needed_virt_capable[static_cast<std::size_t>(j)]) {
                if (find_owned_real_index(j, skey) != SIZE_MAX) continue;
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

        populate_prolonged_virtual_and_sources(virt_flags, remote_real_same_level);
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
                pl_j.size(),
                windowed_parallel_detail::kDefaultParallelGrain,
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

            ld.sources.reserve(ld.owned_real.size() + vf.size() + rem_j.size());

            for (const auto& rn : ld.owned_real) {
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
                    se.wav_coeff = rn.wav_in;
                }
                ld.sources.push_back(std::move(se));
            }

            if (j >= 1 && !vf.empty()) {
                const auto& pl_j = prolong_levels[sj];
                ld.owned_virt.reserve(vf.size());
                std::size_t p_cursor = 0;

                for (const auto& [vkey, flags] : vf) {
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
                }
            }

            for (const auto& se : rem_j) {
                ld.sources.push_back(se);
            }

            std::sort(ld.sources.begin(), ld.sources.end(), [](const auto& a, const auto& b) noexcept {
                return windowed_halo_detail::windowed_key_less(a.key, b.key);
            });
            if (j <= windowed_halo_detail::kMaxPackedKeyLevel) {
                ld.source_packed_keys.resize(ld.sources.size());
                for (std::size_t si = 0; si < ld.sources.size(); ++si) {
                    ld.source_packed_keys[si] =
                        windowed_halo_detail::pack_spatial_key_same_level(ld.sources[si].key);
                }
            }
        }
    }

    // =========================================================================
    // Stage 5a Fast Path: Bottom-Up Virtual Restriction (Flat Sorted Vectors)
    // =========================================================================
    void restrict_all_virtual_unpromoted() {
        const long kk = static_cast<long>(k_);
        const long kk3 = kk * kk * kk;

        struct RestrictAccum {
            Key<3> key;
            Key<3> coarse_leaf;
            Key<3> min_source_key;
            Tensor<T> ds;
        };

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
                if (vn.ds_corr.size() == kk3) {
                    acc_j.push_back(RestrictAccum{
                        vn.key, vn.coarse_leaf_key, vn.key, copy(vn.ds_corr)});
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

            // Pre-insert unique parent slots into `leaf_contribs` and `accum_levels[j - 1]`
            // so parent-grouped restriction runs in exact sorted `vkey` order without hash tables.
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
                active_pslots.size(),
                windowed_parallel_detail::kDefaultParallelGrain,
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

        if (world_.size() == 1) {
            for (const auto& item : leaf_contribs) {
                const int blvl = static_cast<int>(item.key.level());
                const std::size_t ridx = find_owned_real_index(blvl, item.key);
                MADNESS_CHECK(ridx != SIZE_MAX);
                T* dst = levels_[static_cast<std::size_t>(blvl)].owned_real[ridx].s_out.ptr();
                aligned_axpy(kk3, dst, item.ds.ptr(), 1.0);
            }
            return;
        }

        const ProcessID me = world_.rank();
        std::vector<std::vector<windowed_halo_detail::CoarseLeafDeltaMsg<T>>> outbox(
            static_cast<std::size_t>(world_.size()));
        for (auto& item : leaf_contribs) {
            const ProcessID owner = pmap_->owner(item.key);
            outbox[static_cast<std::size_t>(owner)].push_back(
                windowed_halo_detail::CoarseLeafDeltaMsg<T>{
                    item.key, item.min_source_key, std::move(item.ds)});
        }

        for (ProcessID p = 0; p < world_.size(); ++p) {
            auto& batch = outbox[static_cast<std::size_t>(p)];
            if (batch.empty()) continue;
            if (p == me) {
                exchange_->recv_coarse_deltas(batch);
            } else {
                exchange_->send(p,
                                &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_coarse_deltas,
                                batch);
            }
        }
        world_.gop.fence();

        auto received = exchange_->take_coarse_deltas();
        std::sort(received.begin(), received.end(), [](const auto& a, const auto& b) noexcept {
            if (a.leaf_key != b.leaf_key) {
                return windowed_halo_detail::windowed_key_less(a.leaf_key, b.leaf_key);
            }
            return windowed_halo_detail::windowed_key_less(a.sort_tag_key, b.sort_tag_key);
        });

        for (const auto& msg : received) {
            const int blvl = static_cast<int>(msg.leaf_key.level());
            const std::size_t ridx = find_owned_real_index(blvl, msg.leaf_key);
            MADNESS_CHECK(ridx != SIZE_MAX);
            T* dst = levels_[static_cast<std::size_t>(blvl)].owned_real[ridx].s_out.ptr();
            aligned_axpy(kk3, dst, msg.delta_s.ptr(), 1.0);
        }
    }

    // =========================================================================
    // Stage 5a General Path: Restriction + Dynamic Virtual Subtree Promotion
    // =========================================================================
    void restrict_and_split_promoted_virtual(
        typename FunctionImpl<T, 3>::dcT& res_coeffs,
        const TensorArgs& full_targs) {
        const long kk = static_cast<long>(k_);
        const long kk3 = kk * kk * kk;
        const long twok = 2L * kk;
        const Slice s0(0, kk - 1);
        const ProcessID me = world_.rank();

        std::vector<windowed_halo_detail::VirtualSubtreeItemMsg<T>> local_items;
        if (world_.size() == 1) {
            for (int j = 1; j <= n_max_tree_; ++j) {
                for (auto& vn : levels_[static_cast<std::size_t>(j)].owned_virt) {
                    if (!vn.promoted && vn.ds_corr.size() != kk3) continue;
                    local_items.push_back(windowed_halo_detail::VirtualSubtreeItemMsg<T>{
                        vn.coarse_leaf_key,
                        vn.key,
                        static_cast<std::uint8_t>(vn.promoted ? 1u : 0u),
                        vn.ds_corr,
                        vn.ns_virt_out});
                }
            }
        } else {
            std::vector<std::vector<windowed_halo_detail::VirtualSubtreeItemMsg<T>>> outbox(
                static_cast<std::size_t>(world_.size()));
            for (int j = 1; j <= n_max_tree_; ++j) {
                for (auto& vn : levels_[static_cast<std::size_t>(j)].owned_virt) {
                    if (!vn.promoted && vn.ds_corr.size() != kk3) continue;
                    const ProcessID owner = pmap_->owner(vn.coarse_leaf_key);
                    outbox[static_cast<std::size_t>(owner)].push_back(
                        windowed_halo_detail::VirtualSubtreeItemMsg<T>{
                            vn.coarse_leaf_key,
                            vn.key,
                            static_cast<std::uint8_t>(vn.promoted ? 1u : 0u),
                            vn.ds_corr,
                            vn.ns_virt_out});
                }
            }
            for (ProcessID p = 0; p < world_.size(); ++p) {
                auto& batch = outbox[static_cast<std::size_t>(p)];
                if (batch.empty()) continue;
                if (p == me) {
                    exchange_->recv_virtual_subtree_items(batch);
                } else {
                    exchange_->send(
                        p,
                        &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_virtual_subtree_items,
                        batch);
                }
            }
            world_.gop.fence();
            local_items = exchange_->take_subtree_items();
        }

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

        // Parallel hash-free flat-vector restriction across all unpromoted coarse leaves.
        if (!unpromoted_groups.empty()) {
            const long k_dims[3] = {kk, kk, kk};
            std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> u_work1;
            std::array<Tensor<T>, windowed_parallel_detail::kMaxParallelChunks> u_work2;
            windowed_parallel_detail::parallel_for_chunks(
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

        std::vector<std::vector<windowed_halo_detail::OutputNodeInsertMsg<T>>> remote_inserts(
            static_cast<std::size_t>(world_.size()));

        auto emit_node = [&](const Key<3>& key, bool has_children, const Tensor<T>& coeff_tensor) {
            const ProcessID owner = pmap_->owner(key);
            if (owner == me) {
                GenTensor<T> gt(coeff_tensor, full_targs);
                res_coeffs.replace(key, FunctionNode<T, 3>(gt, has_children));
            } else {
                remote_inserts[static_cast<std::size_t>(owner)].push_back(
                    windowed_halo_detail::OutputNodeInsertMsg<T>{
                        key, static_cast<std::uint8_t>(has_children ? 1u : 0u), coeff_tensor});
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

        if (world_.size() > 1) {
            for (ProcessID p = 0; p < world_.size(); ++p) {
                auto& batch = remote_inserts[static_cast<std::size_t>(p)];
                if (!batch.empty()) {
                    exchange_->send(p,
                                    &windowed_halo_detail::WindowedNSHaloExchange<T>::recv_output_inserts,
                                    batch);
                }
            }
            world_.gop.fence();
            for (auto& msg : exchange_->take_output_inserts()) {
                GenTensor<T> gt(msg.coeff, full_targs);
                res_coeffs.replace(msg.key, FunctionNode<T, 3>(gt, msg.has_children != 0));
            }
            world_.gop.fence();
        }
    }

    World& world_;
    int k_ = 0;
    int n_max_tree_ = 0;
    bool dynamic_refine_ = true;
    int max_refine_level_ = 30;
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
    std::unique_ptr<windowed_halo_detail::WindowedNSHaloExchange<T>> exchange_;
};

} // namespace madness

#endif // MADNESS_MRA_WINDOWED_NS_HALO_H__INCLUDED
