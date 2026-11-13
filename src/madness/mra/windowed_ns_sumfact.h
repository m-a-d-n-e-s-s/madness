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

#ifndef MADNESS_MRA_WINDOWED_NS_SUMFACT_H__INCLUDED
#define MADNESS_MRA_WINDOWED_NS_SUMFACT_H__INCLUDED

/// \file windowed_ns_sumfact.h
/// \brief 3D spatial-offset + Legendre-basis sum-factorized execution engine and
///        complete 5-stage `apply_windowed_ns_3d` driver for Version 7 Tail-Completed
///        Windowed Non-Standard (NS) 3D convolutions.
///
/// Implements Algorithm 1 of `multiwavelet_ns_v7.tex` using 100% library-backed
/// contractions (`mTxmq`, `fast_transform`, `general_fast_transform`, `aligned_add`,
/// `aligned_axpy`, `aligned_zero`, `Slice`, `Tensor::normf`):
///   - Stage 1a/1b: Forward NS analysis (`nonstandard_with_leaves`) + single-pass
///                  fine-driven interior-prolonged virtual boundary halo (`virt_cache`)
///                  with scale/tolerance-aware per-pass dynamic bandwidths.
///   - Stage 2:     Cross-level flattened leaf-local Gauss-Christoffel spiky moment pass.
///   - Stage 3:     Root (`j = 0`) low-tail constant `w_const` & ultra-diffuse `n_diff(p) < 0`.
///   - Stage 4(a):  Base scaling `T^{j_0(p)}` injection + Theorem 3.1 `-M(\gamma)^{\otimes 3}`
///                  commutation subtraction on `\partial\mathcal{H}_j^{(1,T)}`.
///   - Stage 4(b):  Cancellation-free `\Gamma^j_p` wavelet decoupling (`Pass 4b-1`: `2k -> k`
///                  on `(0, d_int)`; `Pass 4b-2`: `2k -> 2k` / `k -> 2k` wavelet octants) +
///                  exact constant-subtracted refinement candidate pre-screening + optional
///                  dynamic output tree refinement when `||\hat{d}^j_{U,virt}||_F > \epsilon_{loc}`.
///   - Stage 5a/5b: Single-pass bottom-up adjoint restriction (`h0T, h1T`) onto coarse
///                  leaves + fast input `TreeState` restoration & top-down `reconstruct`.

#include <madness/mra/windowed_ns_halo.h>
#include <madness/tensor/aligned.h>
#include <madness/world/thread.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <type_traits>
#include <utility>
#include <vector>

namespace madness {

namespace windowed_sumfact_detail {

using windowed_parallel_detail::kDefaultParallelGrain;
using windowed_parallel_detail::kMaxParallelChunks;
using windowed_parallel_detail::parallel_for_chunks;

inline constexpr std::size_t kSumFactSweepGrain = 2;

/// Binary-search index of `code` in sorted 64-bit packed `codes`, or `UINT32_MAX` if absent.
inline std::uint32_t find_sorted_packed_index(const std::vector<std::uint64_t>& codes,
                                              std::uint64_t code) noexcept {
    const auto it = std::lower_bound(codes.begin(), codes.end(), code);
    if (it != codes.end() && *it == code) {
        return static_cast<std::uint32_t>(it - codes.begin());
    }
    return UINT32_MAX;
}

/// Binary-search index of `key` in sorted `keys`, or `UINT32_MAX` if absent.
inline std::uint32_t find_sorted_key_index(const std::vector<Key<3>>& keys, const Key<3>& key) noexcept {
    const auto it = std::lower_bound(keys.begin(), keys.end(), key, windowed_halo_detail::windowed_key_less);
    if (it != keys.end() && *it == key) {
        return static_cast<std::uint32_t>(it - keys.begin());
    }
    return UINT32_MAX;
}

/// Linear two-pointer intersection `(src_idx, tgt_idx)` of two sorted same-level `Key<3>` vectors.
inline std::vector<std::pair<std::size_t, std::size_t>>
intersect_sorted_keys(const std::vector<Key<3>>& src_keys, const std::vector<Key<3>>& tgt_keys) {
    std::vector<std::pair<std::size_t, std::size_t>> matches;
    if (src_keys.empty() || tgt_keys.empty()) return matches;
    matches.reserve(std::min(src_keys.size(), tgt_keys.size()));
    std::size_t si = 0;
    std::size_t ti = 0;
    if (static_cast<int>(src_keys[0].level()) <= windowed_halo_detail::kMaxPackedKeyLevel &&
        src_keys[0].level() == tgt_keys[0].level()) {
        std::uint64_t sc = windowed_halo_detail::pack_spatial_key_same_level(src_keys[0]);
        std::uint64_t tc = windowed_halo_detail::pack_spatial_key_same_level(tgt_keys[0]);
        for (;;) {
            if (sc == tc) {
                matches.emplace_back(si, ti);
                ++si;
                ++ti;
                if (si >= src_keys.size() || ti >= tgt_keys.size()) break;
                sc = windowed_halo_detail::pack_spatial_key_same_level(src_keys[si]);
                tc = windowed_halo_detail::pack_spatial_key_same_level(tgt_keys[ti]);
            } else if (sc < tc) {
                ++si;
                if (si >= src_keys.size()) break;
                sc = windowed_halo_detail::pack_spatial_key_same_level(src_keys[si]);
            } else {
                ++ti;
                if (ti >= tgt_keys.size()) break;
                tc = windowed_halo_detail::pack_spatial_key_same_level(tgt_keys[ti]);
            }
        }
        return matches;
    }
    while (si < src_keys.size() && ti < tgt_keys.size()) {
        if (src_keys[si] == tgt_keys[ti]) {
            matches.emplace_back(si, ti);
            ++si;
            ++ti;
        } else if (windowed_halo_detail::windowed_key_less(src_keys[si], tgt_keys[ti])) {
            ++si;
        } else {
            ++ti;
        }
    }
    return matches;
}

/// Two-pointer linear merge of two disjoint sorted `(Key<3>, T*)` streams into a single sorted `(tgt_keys, tgt_ptrs)`.
template <typename T>
inline void merge_sorted_targets(const std::vector<Key<3>>& keys_a,
                                 const std::vector<T*>& ptrs_a,
                                 const std::vector<Key<3>>& keys_b,
                                 const std::vector<T*>& ptrs_b,
                                 std::vector<Key<3>>& out_keys,
                                 std::vector<T*>& out_ptrs) {
    const std::size_t na = keys_a.size();
    const std::size_t nb = keys_b.size();
    out_keys.clear();
    out_ptrs.clear();
    out_keys.reserve(na + nb);
    out_ptrs.reserve(na + nb);

    std::size_t ia = 0;
    std::size_t ib = 0;
    while (ia < na && ib < nb) {
        if (windowed_halo_detail::windowed_key_less(keys_a[ia], keys_b[ib])) {
            out_keys.push_back(keys_a[ia]);
            out_ptrs.push_back(ptrs_a[ia]);
            ++ia;
        } else {
            out_keys.push_back(keys_b[ib]);
            out_ptrs.push_back(ptrs_b[ib]);
            ++ib;
        }
    }
    while (ia < na) {
        out_keys.push_back(keys_a[ia]);
        out_ptrs.push_back(ptrs_a[ia]);
        ++ia;
    }
    while (ib < nb) {
        out_keys.push_back(keys_b[ib]);
        out_ptrs.push_back(ptrs_b[ib]);
        ++ib;
    }
}

/// Compact CSR edge `(index, disp)` for 1D spatial stencil sweeps (`disp = l_target - l_source`).
struct StencilEdge1D {
    std::uint32_t src_idx = 0;
    std::int8_t disp = 0;
};

/// Pre-indexed 3D spatial-offset sum-factorization graph for a fixed level `j` and half-bandwidth `m_max`.
///
/// Built once per `(level, m_max)` using contiguous sorted arrays (zero hash-node allocations) and reused
/// across all active quadrature terms sharing `(level, m_max)`.
struct SumFactGraph3D {
    std::size_t num_x = 0;
    std::size_t num_xy = 0;
    std::size_t num_tgt = 0;

    std::vector<std::uint32_t> x_row_ptr;
    std::vector<StencilEdge1D> x_edges;

    std::vector<std::uint32_t> xy_row_ptr;
    std::vector<StencilEdge1D> xy_edges;

    std::vector<std::uint32_t> active_tgt_idx;
    std::vector<std::uint32_t> z_row_ptr;
    std::vector<StencilEdge1D> z_edges;

    bool empty() const noexcept { return z_edges.empty(); }
};

/// Fast 64-bit packed spatial sweep graph builder for tree levels `n <= kMaxPackedKeyLevel` (`j <= 20`).
/// Produces a 100% bit-for-bit identical `SumFactGraph3D` to `Key<3>` `windowed_key_less` with 8-byte
/// cache-line-dense integers and zero `Key<3>` constructions / hash evaluations.
inline SumFactGraph3D build_sumfact_graph_3d_packed(Level n,
                                                    const std::vector<std::uint64_t>& src_codes,
                                                    const std::vector<std::uint64_t>& tgt_codes,
                                                    int m_max,
                                                    const array_of_bools<3>& is_periodic) {
    SumFactGraph3D graph;
    graph.num_tgt = tgt_codes.size();
    if (src_codes.empty() || tgt_codes.empty() || m_max < 0) {
        graph.x_row_ptr.assign(1, 0u);
        graph.xy_row_ptr.assign(1, 0u);
        graph.z_row_ptr.assign(1, 0u);
        return graph;
    }

    const Translation dom = (Translation(1) << n);
    const std::size_t span = static_cast<std::size_t>(2 * m_max + 1);

    if (src_codes.size() <= tgt_codes.size()) {
        struct PackedEdgePair {
            std::uint64_t code = 0;
            StencilEdge1D edge;
        };
        auto packed_edge_less = [](const PackedEdgePair& a, const PackedEdgePair& b) noexcept {
            if (a.code != b.code) return a.code < b.code;
            return a.edge.disp < b.edge.disp;
        };

        std::vector<PackedEdgePair> raw_x;
        raw_x.reserve(src_codes.size() * span);
        for (std::size_t si = 0; si < src_codes.size(); ++si) {
            const std::uint64_t scode = src_codes[si];
            const Translation lx =
                static_cast<Translation>((scode >> 42) & windowed_halo_detail::kPackedCoordMask);
            const Translation ly =
                static_cast<Translation>((scode >> 21) & windowed_halo_detail::kPackedCoordMask);
            const Translation lz =
                static_cast<Translation>(scode & windowed_halo_detail::kPackedCoordMask);
            for (int dx = -m_max; dx <= m_max; ++dx) {
                std::uint64_t x_code = 0;
                if (windowed_halo_detail::shift_packed_coords_in_domain(
                        dom, lx, ly, lz, dx, 0, 0, is_periodic, x_code)) {
                    raw_x.push_back(PackedEdgePair{
                        x_code,
                        StencilEdge1D{static_cast<std::uint32_t>(si),
                                      static_cast<std::int8_t>(dx)}});
                }
            }
        }
        if (raw_x.empty()) {
            graph.x_row_ptr.assign(1, 0u);
            graph.xy_row_ptr.assign(1, 0u);
            graph.z_row_ptr.assign(1, 0u);
            return graph;
        }
        std::sort(raw_x.begin(), raw_x.end(), packed_edge_less);

        std::vector<std::uint64_t> cand_x_codes;
        std::vector<std::uint32_t> cand_x_row_ptr;
        std::vector<StencilEdge1D> cand_x_edges;
        cand_x_codes.reserve(raw_x.size());
        cand_x_row_ptr.reserve(raw_x.size() + 1);
        cand_x_edges.reserve(raw_x.size());
        cand_x_row_ptr.push_back(0u);
        for (std::size_t i = 0; i < raw_x.size();) {
            const std::uint64_t cur_c = raw_x[i].code;
            cand_x_codes.push_back(cur_c);
            while (i < raw_x.size() && raw_x[i].code == cur_c) {
                cand_x_edges.push_back(raw_x[i].edge);
                ++i;
            }
            cand_x_row_ptr.push_back(static_cast<std::uint32_t>(cand_x_edges.size()));
        }

        const std::size_t nx0 = cand_x_codes.size();
        std::vector<PackedEdgePair> raw_xy;
        raw_xy.reserve(nx0 * span);
        for (std::size_t xi = 0; xi < nx0; ++xi) {
            const std::uint64_t x_code = cand_x_codes[xi];
            const Translation lx =
                static_cast<Translation>((x_code >> 42) & windowed_halo_detail::kPackedCoordMask);
            const Translation ly =
                static_cast<Translation>((x_code >> 21) & windowed_halo_detail::kPackedCoordMask);
            const Translation lz =
                static_cast<Translation>(x_code & windowed_halo_detail::kPackedCoordMask);
            for (int dy = -m_max; dy <= m_max; ++dy) {
                std::uint64_t xy_code = 0;
                if (windowed_halo_detail::shift_packed_coords_in_domain(
                        dom, lx, ly, lz, 0, dy, 0, is_periodic, xy_code)) {
                    raw_xy.push_back(PackedEdgePair{
                        xy_code,
                        StencilEdge1D{static_cast<std::uint32_t>(xi),
                                      static_cast<std::int8_t>(dy)}});
                }
            }
        }
        if (raw_xy.empty()) {
            graph.x_row_ptr.assign(1, 0u);
            graph.xy_row_ptr.assign(1, 0u);
            graph.z_row_ptr.assign(1, 0u);
            return graph;
        }
        std::sort(raw_xy.begin(), raw_xy.end(), packed_edge_less);

        std::vector<std::uint64_t> cand_xy_codes;
        std::vector<std::uint32_t> cand_xy_row_ptr;
        std::vector<StencilEdge1D> cand_xy_edges;
        cand_xy_codes.reserve(raw_xy.size());
        cand_xy_row_ptr.reserve(raw_xy.size() + 1);
        cand_xy_edges.reserve(raw_xy.size());
        cand_xy_row_ptr.push_back(0u);
        for (std::size_t i = 0; i < raw_xy.size();) {
            const std::uint64_t cur_c = raw_xy[i].code;
            cand_xy_codes.push_back(cur_c);
            while (i < raw_xy.size() && raw_xy[i].code == cur_c) {
                cand_xy_edges.push_back(raw_xy[i].edge);
                ++i;
            }
            cand_xy_row_ptr.push_back(static_cast<std::uint32_t>(cand_xy_edges.size()));
        }

        const std::size_t nxy0 = cand_xy_codes.size();
        std::vector<std::uint8_t> xy_used(nxy0, 0u);
        graph.z_row_ptr.push_back(0u);

        if (nxy0 <= tgt_codes.size()) {
            struct TgtZEdge {
                std::uint32_t t_idx = 0;
                StencilEdge1D edge;
            };
            std::vector<TgtZEdge> raw_z;
            raw_z.reserve(nxy0 * span);
            for (std::size_t xyi = 0; xyi < nxy0; ++xyi) {
                const std::uint64_t xy_code = cand_xy_codes[xyi];
                const Translation lx = static_cast<Translation>(
                    (xy_code >> 42) & windowed_halo_detail::kPackedCoordMask);
                const Translation ly = static_cast<Translation>(
                    (xy_code >> 21) & windowed_halo_detail::kPackedCoordMask);
                const Translation lz =
                    static_cast<Translation>(xy_code & windowed_halo_detail::kPackedCoordMask);
                for (int dz = -m_max; dz <= m_max; ++dz) {
                    std::uint64_t t_code = 0;
                    if (!windowed_halo_detail::shift_packed_coords_in_domain(
                            dom, lx, ly, lz, 0, 0, dz, is_periodic, t_code)) {
                        continue;
                    }
                    const std::uint32_t ti = find_sorted_packed_index(tgt_codes, t_code);
                    if (ti != UINT32_MAX) {
                        xy_used[xyi] = 1u;
                        raw_z.push_back(TgtZEdge{
                            ti,
                            StencilEdge1D{static_cast<std::uint32_t>(xyi),
                                          static_cast<std::int8_t>(dz)}});
                    }
                }
            }
            std::sort(raw_z.begin(), raw_z.end(), [](const TgtZEdge& a, const TgtZEdge& b) noexcept {
                if (a.t_idx != b.t_idx) return a.t_idx < b.t_idx;
                return a.edge.disp < b.edge.disp;
            });
            graph.active_tgt_idx.reserve(raw_z.size());
            graph.z_row_ptr.reserve(raw_z.size() + 1);
            graph.z_edges.reserve(raw_z.size());
            for (std::size_t i = 0; i < raw_z.size();) {
                const std::uint32_t cur_t = raw_z[i].t_idx;
                graph.active_tgt_idx.push_back(cur_t);
                while (i < raw_z.size() && raw_z[i].t_idx == cur_t) {
                    graph.z_edges.push_back(raw_z[i].edge);
                    ++i;
                }
                graph.z_row_ptr.push_back(static_cast<std::uint32_t>(graph.z_edges.size()));
            }
        } else {
            graph.active_tgt_idx.reserve(graph.num_tgt);
            graph.z_row_ptr.reserve(graph.num_tgt + 1);
            graph.z_edges.reserve(graph.num_tgt * span);
            for (std::size_t t = 0; t < graph.num_tgt; ++t) {
                const std::uint64_t t_code = tgt_codes[t];
                const Translation lx = static_cast<Translation>(
                    (t_code >> 42) & windowed_halo_detail::kPackedCoordMask);
                const Translation ly = static_cast<Translation>(
                    (t_code >> 21) & windowed_halo_detail::kPackedCoordMask);
                const Translation lz =
                    static_cast<Translation>(t_code & windowed_halo_detail::kPackedCoordMask);
                const std::size_t e_start = graph.z_edges.size();
                for (int dz = -m_max; dz <= m_max; ++dz) {
                    std::uint64_t xy_code = 0;
                    if (!windowed_halo_detail::shift_packed_coords_in_domain(
                            dom, lx, ly, lz, 0, 0, -dz, is_periodic, xy_code)) {
                        continue;
                    }
                    const std::uint32_t xyi = find_sorted_packed_index(cand_xy_codes, xy_code);
                    if (xyi != UINT32_MAX) {
                        xy_used[xyi] = 1u;
                        graph.z_edges.push_back(StencilEdge1D{xyi, static_cast<std::int8_t>(dz)});
                    }
                }
                if (graph.z_edges.size() > e_start) {
                    graph.active_tgt_idx.push_back(static_cast<std::uint32_t>(t));
                    graph.z_row_ptr.push_back(static_cast<std::uint32_t>(graph.z_edges.size()));
                }
            }
        }

        if (graph.z_edges.empty()) {
            graph.x_row_ptr.assign(1, 0u);
            graph.xy_row_ptr.assign(1, 0u);
            graph.z_row_ptr.assign(1, 0u);
            return graph;
        }

        std::vector<std::uint32_t> xy_remap(nxy0, UINT32_MAX);
        std::vector<std::uint8_t> x_used(nx0, 0u);
        std::uint32_t kept_xy = 0;
        for (std::size_t xyi = 0; xyi < nxy0; ++xyi) {
            if (!xy_used[xyi]) continue;
            xy_remap[xyi] = kept_xy++;
            for (std::uint32_t e = cand_xy_row_ptr[xyi]; e < cand_xy_row_ptr[xyi + 1]; ++e) {
                x_used[cand_xy_edges[e].src_idx] = 1u;
            }
        }

        std::vector<std::uint32_t> x_remap(nx0, UINT32_MAX);
        graph.x_row_ptr.reserve(nx0 + 1);
        graph.x_edges.reserve(cand_x_edges.size());
        graph.x_row_ptr.push_back(0u);
        for (std::size_t xi = 0; xi < nx0; ++xi) {
            if (!x_used[xi]) continue;
            x_remap[xi] = static_cast<std::uint32_t>(graph.num_x++);
            for (std::uint32_t e = cand_x_row_ptr[xi]; e < cand_x_row_ptr[xi + 1]; ++e) {
                graph.x_edges.push_back(cand_x_edges[e]);
            }
            graph.x_row_ptr.push_back(static_cast<std::uint32_t>(graph.x_edges.size()));
        }

        graph.num_xy = static_cast<std::size_t>(kept_xy);
        graph.xy_row_ptr.reserve(graph.num_xy + 1);
        graph.xy_edges.reserve(cand_xy_edges.size());
        graph.xy_row_ptr.push_back(0u);
        for (std::size_t xyi = 0; xyi < nxy0; ++xyi) {
            if (!xy_used[xyi]) continue;
            for (std::uint32_t e = cand_xy_row_ptr[xyi]; e < cand_xy_row_ptr[xyi + 1]; ++e) {
                StencilEdge1D ed = cand_xy_edges[e];
                ed.src_idx = x_remap[ed.src_idx];
                graph.xy_edges.push_back(ed);
            }
            graph.xy_row_ptr.push_back(static_cast<std::uint32_t>(graph.xy_edges.size()));
        }

        for (StencilEdge1D& zed : graph.z_edges) {
            zed.src_idx = xy_remap[zed.src_idx];
        }
        return graph;
    }

    // Backward candidate discovery (`|src| > |tgt|`) on 64-bit packed codes.
    std::vector<std::uint64_t> cand_xy;
    cand_xy.reserve(tgt_codes.size() * span);
    for (std::uint64_t t_code : tgt_codes) {
        const Translation lx =
            static_cast<Translation>((t_code >> 42) & windowed_halo_detail::kPackedCoordMask);
        const Translation ly =
            static_cast<Translation>((t_code >> 21) & windowed_halo_detail::kPackedCoordMask);
        const Translation lz =
            static_cast<Translation>(t_code & windowed_halo_detail::kPackedCoordMask);
        for (int dz = -m_max; dz <= m_max; ++dz) {
            std::uint64_t xy_code = 0;
            if (windowed_halo_detail::shift_packed_coords_in_domain(
                    dom, lx, ly, lz, 0, 0, -dz, is_periodic, xy_code)) {
                cand_xy.push_back(xy_code);
            }
        }
    }
    windowed_halo_detail::sort_unique_packed_codes(cand_xy);

    std::vector<std::uint64_t> cand_x;
    cand_x.reserve(cand_xy.size() * span);
    for (std::uint64_t xy_code : cand_xy) {
        const Translation lx =
            static_cast<Translation>((xy_code >> 42) & windowed_halo_detail::kPackedCoordMask);
        const Translation ly =
            static_cast<Translation>((xy_code >> 21) & windowed_halo_detail::kPackedCoordMask);
        const Translation lz =
            static_cast<Translation>(xy_code & windowed_halo_detail::kPackedCoordMask);
        for (int dy = -m_max; dy <= m_max; ++dy) {
            std::uint64_t x_code = 0;
            if (windowed_halo_detail::shift_packed_coords_in_domain(
                    dom, lx, ly, lz, 0, -dy, 0, is_periodic, x_code)) {
                cand_x.push_back(x_code);
            }
        }
    }
    windowed_halo_detail::sort_unique_packed_codes(cand_x);

    std::vector<std::uint64_t> valid_x_codes;
    valid_x_codes.reserve(cand_x.size());
    graph.x_row_ptr.reserve(cand_x.size() + 1);
    graph.x_edges.reserve(cand_x.size() * span);
    graph.x_row_ptr.push_back(0u);

    for (std::uint64_t x_code : cand_x) {
        const Translation lx =
            static_cast<Translation>((x_code >> 42) & windowed_halo_detail::kPackedCoordMask);
        const Translation ly =
            static_cast<Translation>((x_code >> 21) & windowed_halo_detail::kPackedCoordMask);
        const Translation lz =
            static_cast<Translation>(x_code & windowed_halo_detail::kPackedCoordMask);
        const std::size_t edge_start = graph.x_edges.size();
        for (int dx = -m_max; dx <= m_max; ++dx) {
            std::uint64_t s_code = 0;
            if (!windowed_halo_detail::shift_packed_coords_in_domain(
                    dom, lx, ly, lz, -dx, 0, 0, is_periodic, s_code)) {
                continue;
            }
            const std::uint32_t s_idx = find_sorted_packed_index(src_codes, s_code);
            if (s_idx != UINT32_MAX) {
                graph.x_edges.push_back(StencilEdge1D{s_idx, static_cast<std::int8_t>(dx)});
            }
        }
        if (graph.x_edges.size() > edge_start) {
            ++graph.num_x;
            valid_x_codes.push_back(x_code);
            graph.x_row_ptr.push_back(static_cast<std::uint32_t>(graph.x_edges.size()));
        }
    }

    if (graph.num_x == 0) {
        graph.xy_row_ptr.assign(1, 0u);
        graph.z_row_ptr.assign(1, 0u);
        return graph;
    }

    std::vector<std::uint64_t> valid_xy_codes;
    valid_xy_codes.reserve(cand_xy.size());
    graph.xy_row_ptr.reserve(cand_xy.size() + 1);
    graph.xy_edges.reserve(cand_xy.size() * span);
    graph.xy_row_ptr.push_back(0u);

    for (std::uint64_t xy_code : cand_xy) {
        const Translation lx =
            static_cast<Translation>((xy_code >> 42) & windowed_halo_detail::kPackedCoordMask);
        const Translation ly =
            static_cast<Translation>((xy_code >> 21) & windowed_halo_detail::kPackedCoordMask);
        const Translation lz =
            static_cast<Translation>(xy_code & windowed_halo_detail::kPackedCoordMask);
        const std::size_t edge_start = graph.xy_edges.size();
        for (int dy = -m_max; dy <= m_max; ++dy) {
            std::uint64_t x_code = 0;
            if (!windowed_halo_detail::shift_packed_coords_in_domain(
                    dom, lx, ly, lz, 0, -dy, 0, is_periodic, x_code)) {
                continue;
            }
            const std::uint32_t x_idx = find_sorted_packed_index(valid_x_codes, x_code);
            if (x_idx != UINT32_MAX) {
                graph.xy_edges.push_back(StencilEdge1D{x_idx, static_cast<std::int8_t>(dy)});
            }
        }
        if (graph.xy_edges.size() > edge_start) {
            ++graph.num_xy;
            valid_xy_codes.push_back(xy_code);
            graph.xy_row_ptr.push_back(static_cast<std::uint32_t>(graph.xy_edges.size()));
        }
    }

    if (graph.num_xy == 0) {
        graph.z_row_ptr.assign(1, 0u);
        return graph;
    }

    graph.active_tgt_idx.reserve(graph.num_tgt);
    graph.z_row_ptr.reserve(graph.num_tgt + 1);
    graph.z_edges.reserve(graph.num_tgt * span);
    graph.z_row_ptr.push_back(0u);
    for (std::size_t t = 0; t < graph.num_tgt; ++t) {
        const std::uint64_t t_code = tgt_codes[t];
        const Translation lx =
            static_cast<Translation>((t_code >> 42) & windowed_halo_detail::kPackedCoordMask);
        const Translation ly =
            static_cast<Translation>((t_code >> 21) & windowed_halo_detail::kPackedCoordMask);
        const Translation lz =
            static_cast<Translation>(t_code & windowed_halo_detail::kPackedCoordMask);
        const std::size_t edge_start = graph.z_edges.size();
        for (int dz = -m_max; dz <= m_max; ++dz) {
            std::uint64_t xy_code = 0;
            if (!windowed_halo_detail::shift_packed_coords_in_domain(
                    dom, lx, ly, lz, 0, 0, -dz, is_periodic, xy_code)) {
                continue;
            }
            const std::uint32_t xy_idx = find_sorted_packed_index(valid_xy_codes, xy_code);
            if (xy_idx != UINT32_MAX) {
                graph.z_edges.push_back(StencilEdge1D{xy_idx, static_cast<std::int8_t>(dz)});
            }
        }
        if (graph.z_edges.size() > edge_start) {
            graph.active_tgt_idx.push_back(static_cast<std::uint32_t>(t));
            graph.z_row_ptr.push_back(static_cast<std::uint32_t>(graph.z_edges.size()));
        }
    }

    return graph;
}

/// Build a pruned, deterministically ordered 3D spatial sweep graph `D_src -> D_x -> D_xy -> D_tgt`.
/// Requires BOTH `src_keys` and `tgt_keys` to be sorted according to `windowed_key_less`.
/// Uses 64-bit packed spatial coordinates for `j <= 20` and full `Key<3>` fallback for `j > 20`.
inline SumFactGraph3D build_sumfact_graph_3d(const std::vector<Key<3>>& src_keys,
                                             const std::vector<Key<3>>& tgt_keys,
                                             int m_max,
                                             const array_of_bools<3>& is_periodic) {
    if (!src_keys.empty() && !tgt_keys.empty() && m_max >= 0 &&
        static_cast<int>(src_keys[0].level()) <= windowed_halo_detail::kMaxPackedKeyLevel &&
        src_keys[0].level() == tgt_keys[0].level()) {
        std::vector<std::uint64_t> src_codes(src_keys.size());
        std::vector<std::uint64_t> tgt_codes(tgt_keys.size());
        for (std::size_t i = 0; i < src_keys.size(); ++i) {
            src_codes[i] = windowed_halo_detail::pack_spatial_key_same_level(src_keys[i]);
        }
        for (std::size_t i = 0; i < tgt_keys.size(); ++i) {
            tgt_codes[i] = windowed_halo_detail::pack_spatial_key_same_level(tgt_keys[i]);
        }
        return build_sumfact_graph_3d_packed(
            src_keys[0].level(), src_codes, tgt_codes, m_max, is_periodic);
    }

    SumFactGraph3D graph;
    graph.num_tgt = tgt_keys.size();
    if (src_keys.empty() || tgt_keys.empty() || m_max < 0) {
        graph.x_row_ptr.assign(1, 0u);
        graph.xy_row_ptr.assign(1, 0u);
        graph.z_row_ptr.assign(1, 0u);
        return graph;
    }

    const std::size_t span = static_cast<std::size_t>(2 * m_max + 1);

    if (src_keys.size() <= tgt_keys.size()) {
        struct KeyEdgePair {
            Key<3> key;
            StencilEdge1D edge;
        };
        auto key_edge_less = [](const KeyEdgePair& a, const KeyEdgePair& b) noexcept {
            if (a.key != b.key) {
                return windowed_halo_detail::windowed_key_less(a.key, b.key);
            }
            return a.edge.disp < b.edge.disp;
        };

        std::vector<KeyEdgePair> raw_x;
        raw_x.reserve(src_keys.size() * span);
        for (std::size_t si = 0; si < src_keys.size(); ++si) {
            const Key<3>& skey = src_keys[si];
            for (int dx = -m_max; dx <= m_max; ++dx) {
                Key<3> x_key;
                if (windowed_halo_detail::shift_key_in_domain(skey, dx, 0, 0, is_periodic, x_key)) {
                    raw_x.push_back(KeyEdgePair{
                        x_key, StencilEdge1D{static_cast<std::uint32_t>(si), static_cast<std::int8_t>(dx)}});
                }
            }
        }
        if (raw_x.empty()) {
            graph.x_row_ptr.assign(1, 0u);
            graph.xy_row_ptr.assign(1, 0u);
            graph.z_row_ptr.assign(1, 0u);
            return graph;
        }
        std::sort(raw_x.begin(), raw_x.end(), key_edge_less);

        std::vector<Key<3>> cand_x_keys;
        std::vector<std::uint32_t> cand_x_row_ptr;
        std::vector<StencilEdge1D> cand_x_edges;
        cand_x_keys.reserve(raw_x.size());
        cand_x_row_ptr.reserve(raw_x.size() + 1);
        cand_x_edges.reserve(raw_x.size());
        cand_x_row_ptr.push_back(0u);
        for (std::size_t i = 0; i < raw_x.size();) {
            const Key<3> cur_k = raw_x[i].key;
            cand_x_keys.push_back(cur_k);
            while (i < raw_x.size() && raw_x[i].key == cur_k) {
                cand_x_edges.push_back(raw_x[i].edge);
                ++i;
            }
            cand_x_row_ptr.push_back(static_cast<std::uint32_t>(cand_x_edges.size()));
        }

        const std::size_t nx0 = cand_x_keys.size();
        std::vector<KeyEdgePair> raw_xy;
        raw_xy.reserve(nx0 * span);
        for (std::size_t xi = 0; xi < nx0; ++xi) {
            const Key<3>& x_key = cand_x_keys[xi];
            for (int dy = -m_max; dy <= m_max; ++dy) {
                Key<3> xy_key;
                if (windowed_halo_detail::shift_key_in_domain(x_key, 0, dy, 0, is_periodic, xy_key)) {
                    raw_xy.push_back(KeyEdgePair{
                        xy_key, StencilEdge1D{static_cast<std::uint32_t>(xi), static_cast<std::int8_t>(dy)}});
                }
            }
        }
        if (raw_xy.empty()) {
            graph.x_row_ptr.assign(1, 0u);
            graph.xy_row_ptr.assign(1, 0u);
            graph.z_row_ptr.assign(1, 0u);
            return graph;
        }
        std::sort(raw_xy.begin(), raw_xy.end(), key_edge_less);

        std::vector<Key<3>> cand_xy_keys;
        std::vector<std::uint32_t> cand_xy_row_ptr;
        std::vector<StencilEdge1D> cand_xy_edges;
        cand_xy_keys.reserve(raw_xy.size());
        cand_xy_row_ptr.reserve(raw_xy.size() + 1);
        cand_xy_edges.reserve(raw_xy.size());
        cand_xy_row_ptr.push_back(0u);
        for (std::size_t i = 0; i < raw_xy.size();) {
            const Key<3> cur_k = raw_xy[i].key;
            cand_xy_keys.push_back(cur_k);
            while (i < raw_xy.size() && raw_xy[i].key == cur_k) {
                cand_xy_edges.push_back(raw_xy[i].edge);
                ++i;
            }
            cand_xy_row_ptr.push_back(static_cast<std::uint32_t>(cand_xy_edges.size()));
        }

        const std::size_t nxy0 = cand_xy_keys.size();
        std::vector<std::uint8_t> xy_used(nxy0, 0u);
        graph.z_row_ptr.push_back(0u);

        if (nxy0 <= tgt_keys.size()) {
            struct TgtZEdge {
                std::uint32_t t_idx = 0;
                StencilEdge1D edge;
            };
            std::vector<TgtZEdge> raw_z;
            raw_z.reserve(nxy0 * span);
            for (std::size_t xyi = 0; xyi < nxy0; ++xyi) {
                const Key<3>& xy_key = cand_xy_keys[xyi];
                for (int dz = -m_max; dz <= m_max; ++dz) {
                    Key<3> tkey;
                    if (!windowed_halo_detail::shift_key_in_domain(xy_key, 0, 0, dz, is_periodic, tkey)) {
                        continue;
                    }
                    const std::uint32_t ti = find_sorted_key_index(tgt_keys, tkey);
                    if (ti != UINT32_MAX) {
                        xy_used[xyi] = 1u;
                        raw_z.push_back(TgtZEdge{
                            ti, StencilEdge1D{static_cast<std::uint32_t>(xyi), static_cast<std::int8_t>(dz)}});
                    }
                }
            }
            std::sort(raw_z.begin(), raw_z.end(), [](const TgtZEdge& a, const TgtZEdge& b) noexcept {
                if (a.t_idx != b.t_idx) return a.t_idx < b.t_idx;
                return a.edge.disp < b.edge.disp;
            });
            graph.active_tgt_idx.reserve(raw_z.size());
            graph.z_row_ptr.reserve(raw_z.size() + 1);
            graph.z_edges.reserve(raw_z.size());
            for (std::size_t i = 0; i < raw_z.size();) {
                const std::uint32_t cur_t = raw_z[i].t_idx;
                graph.active_tgt_idx.push_back(cur_t);
                while (i < raw_z.size() && raw_z[i].t_idx == cur_t) {
                    graph.z_edges.push_back(raw_z[i].edge);
                    ++i;
                }
                graph.z_row_ptr.push_back(static_cast<std::uint32_t>(graph.z_edges.size()));
            }
        } else {
            graph.active_tgt_idx.reserve(graph.num_tgt);
            graph.z_row_ptr.reserve(graph.num_tgt + 1);
            graph.z_edges.reserve(graph.num_tgt * span);
            for (std::size_t t = 0; t < graph.num_tgt; ++t) {
                const Key<3>& tkey = tgt_keys[t];
                const std::size_t e_start = graph.z_edges.size();
                for (int dz = -m_max; dz <= m_max; ++dz) {
                    Key<3> xy_key;
                    if (!windowed_halo_detail::shift_key_in_domain(tkey, 0, 0, -dz, is_periodic, xy_key)) {
                        continue;
                    }
                    const std::uint32_t xyi = find_sorted_key_index(cand_xy_keys, xy_key);
                    if (xyi != UINT32_MAX) {
                        xy_used[xyi] = 1u;
                        graph.z_edges.push_back(StencilEdge1D{xyi, static_cast<std::int8_t>(dz)});
                    }
                }
                if (graph.z_edges.size() > e_start) {
                    graph.active_tgt_idx.push_back(static_cast<std::uint32_t>(t));
                    graph.z_row_ptr.push_back(static_cast<std::uint32_t>(graph.z_edges.size()));
                }
            }
        }

        if (graph.z_edges.empty()) {
            graph.x_row_ptr.assign(1, 0u);
            graph.xy_row_ptr.assign(1, 0u);
            graph.z_row_ptr.assign(1, 0u);
            return graph;
        }

        std::vector<std::uint32_t> xy_remap(nxy0, UINT32_MAX);
        std::vector<std::uint8_t> x_used(nx0, 0u);
        std::uint32_t kept_xy = 0;
        for (std::size_t xyi = 0; xyi < nxy0; ++xyi) {
            if (!xy_used[xyi]) continue;
            xy_remap[xyi] = kept_xy++;
            for (std::uint32_t e = cand_xy_row_ptr[xyi]; e < cand_xy_row_ptr[xyi + 1]; ++e) {
                x_used[cand_xy_edges[e].src_idx] = 1u;
            }
        }

        std::vector<std::uint32_t> x_remap(nx0, UINT32_MAX);
        graph.x_row_ptr.reserve(nx0 + 1);
        graph.x_edges.reserve(cand_x_edges.size());
        graph.x_row_ptr.push_back(0u);
        for (std::size_t xi = 0; xi < nx0; ++xi) {
            if (!x_used[xi]) continue;
            x_remap[xi] = static_cast<std::uint32_t>(graph.num_x++);
            for (std::uint32_t e = cand_x_row_ptr[xi]; e < cand_x_row_ptr[xi + 1]; ++e) {
                graph.x_edges.push_back(cand_x_edges[e]);
            }
            graph.x_row_ptr.push_back(static_cast<std::uint32_t>(graph.x_edges.size()));
        }

        graph.num_xy = static_cast<std::size_t>(kept_xy);
        graph.xy_row_ptr.reserve(graph.num_xy + 1);
        graph.xy_edges.reserve(cand_xy_edges.size());
        graph.xy_row_ptr.push_back(0u);
        for (std::size_t xyi = 0; xyi < nxy0; ++xyi) {
            if (!xy_used[xyi]) continue;
            for (std::uint32_t e = cand_xy_row_ptr[xyi]; e < cand_xy_row_ptr[xyi + 1]; ++e) {
                StencilEdge1D ed = cand_xy_edges[e];
                ed.src_idx = x_remap[ed.src_idx];
                graph.xy_edges.push_back(ed);
            }
            graph.xy_row_ptr.push_back(static_cast<std::uint32_t>(graph.xy_edges.size()));
        }

        for (StencilEdge1D& zed : graph.z_edges) {
            zed.src_idx = xy_remap[zed.src_idx];
        }
        return graph;
    }

    // Backward candidate discovery from D_tgt -> candidate D_xy -> candidate D_x using flat vectors.
    std::vector<Key<3>> cand_xy;
    cand_xy.reserve(tgt_keys.size() * span);
    for (const Key<3>& tkey : tgt_keys) {
        for (int dz = -m_max; dz <= m_max; ++dz) {
            Key<3> xy_key;
            if (windowed_halo_detail::shift_key_in_domain(tkey, 0, 0, -dz, is_periodic, xy_key)) {
                cand_xy.push_back(xy_key);
            }
        }
    }
    std::sort(cand_xy.begin(), cand_xy.end(), windowed_halo_detail::windowed_key_less);
    cand_xy.erase(std::unique(cand_xy.begin(), cand_xy.end()), cand_xy.end());

    std::vector<Key<3>> cand_x;
    cand_x.reserve(cand_xy.size() * span);
    for (const Key<3>& xy_key : cand_xy) {
        for (int dy = -m_max; dy <= m_max; ++dy) {
            Key<3> x_key;
            if (windowed_halo_detail::shift_key_in_domain(xy_key, 0, -dy, 0, is_periodic, x_key)) {
                cand_x.push_back(x_key);
            }
        }
    }
    std::sort(cand_x.begin(), cand_x.end(), windowed_halo_detail::windowed_key_less);
    cand_x.erase(std::unique(cand_x.begin(), cand_x.end()), cand_x.end());

    // Prune D_x to keys with >= 1 valid source in sorted `src_keys`.
    std::vector<Key<3>> valid_x_keys;
    valid_x_keys.reserve(cand_x.size());
    graph.x_row_ptr.reserve(cand_x.size() + 1);
    graph.x_edges.reserve(cand_x.size() * span);
    graph.x_row_ptr.push_back(0u);

    for (const Key<3>& x_key : cand_x) {
        const std::size_t edge_start = graph.x_edges.size();
        for (int dx = -m_max; dx <= m_max; ++dx) {
            Key<3> s_key;
            if (!windowed_halo_detail::shift_key_in_domain(x_key, -dx, 0, 0, is_periodic, s_key)) {
                continue;
            }
            const std::uint32_t s_idx = find_sorted_key_index(src_keys, s_key);
            if (s_idx != UINT32_MAX) {
                graph.x_edges.push_back(StencilEdge1D{s_idx, static_cast<std::int8_t>(dx)});
            }
        }
        if (graph.x_edges.size() > edge_start) {
            ++graph.num_x;
            valid_x_keys.push_back(x_key);
            graph.x_row_ptr.push_back(static_cast<std::uint32_t>(graph.x_edges.size()));
        }
    }

    if (graph.num_x == 0) {
        graph.xy_row_ptr.assign(1, 0u);
        graph.z_row_ptr.assign(1, 0u);
        return graph;
    }

    // Prune D_xy to keys with >= 1 valid entry in sorted `valid_x_keys`.
    std::vector<Key<3>> valid_xy_keys;
    valid_xy_keys.reserve(cand_xy.size());
    graph.xy_row_ptr.reserve(cand_xy.size() + 1);
    graph.xy_edges.reserve(cand_xy.size() * span);
    graph.xy_row_ptr.push_back(0u);

    for (const Key<3>& xy_key : cand_xy) {
        const std::size_t edge_start = graph.xy_edges.size();
        for (int dy = -m_max; dy <= m_max; ++dy) {
            Key<3> x_key;
            if (!windowed_halo_detail::shift_key_in_domain(xy_key, 0, -dy, 0, is_periodic, x_key)) {
                continue;
            }
            const std::uint32_t x_idx = find_sorted_key_index(valid_x_keys, x_key);
            if (x_idx != UINT32_MAX) {
                graph.xy_edges.push_back(StencilEdge1D{x_idx, static_cast<std::int8_t>(dy)});
            }
        }
        if (graph.xy_edges.size() > edge_start) {
            ++graph.num_xy;
            valid_xy_keys.push_back(xy_key);
            graph.xy_row_ptr.push_back(static_cast<std::uint32_t>(graph.xy_edges.size()));
        }
    }

    if (graph.num_xy == 0) {
        graph.z_row_ptr.assign(1, 0u);
        return graph;
    }

    // Populate compact Z-sweep edges into active targets (preserving caller's target index order).
    graph.active_tgt_idx.reserve(graph.num_tgt);
    graph.z_row_ptr.reserve(graph.num_tgt + 1);
    graph.z_edges.reserve(graph.num_tgt * span);
    graph.z_row_ptr.push_back(0u);
    for (std::size_t t = 0; t < graph.num_tgt; ++t) {
        const Key<3>& tkey = tgt_keys[t];
        const std::size_t edge_start = graph.z_edges.size();
        for (int dz = -m_max; dz <= m_max; ++dz) {
            Key<3> xy_key;
            if (!windowed_halo_detail::shift_key_in_domain(tkey, 0, 0, -dz, is_periodic, xy_key)) {
                continue;
            }
            const std::uint32_t xy_idx = find_sorted_key_index(valid_xy_keys, xy_key);
            if (xy_idx != UINT32_MAX) {
                graph.z_edges.push_back(StencilEdge1D{xy_idx, static_cast<std::int8_t>(dz)});
            }
        }
        if (graph.z_edges.size() > edge_start) {
            graph.active_tgt_idx.push_back(static_cast<std::uint32_t>(t));
            graph.z_row_ptr.push_back(static_cast<std::uint32_t>(graph.z_edges.size()));
        }
    }

    return graph;
}

/// Per-displacement 1D operator descriptor for `execute_sumfact_3d` (dense `(d_in, d_out)` or
/// 2-step Tier-1 `s_dense_table` SVD `U_by_rank[r]` `(d_in, r)` + `SVT` `(r, d_out)`).
struct SlotOp1D {
    bool active = false;
    bool use_svd = false;
    long rank = 0;
    const double* dense_ptr = nullptr;
    const double* u_ptr = nullptr;
    const double* svt_ptr = nullptr;
};

inline SlotOp1D make_slot_op_1d(int r,
                                long d_in,
                                long d_out,
                                const Tensor<double>& dense_mat,
                                const StencilSVD1D& svd) noexcept {
    SlotOp1D op;
    if (r <= 0 || !(svd.normf > 0.0)) return op;
    op.active = true;
    const long rl = static_cast<long>(r);
    const bool svd_faster =
        (rl == 1L) ? (rl * (d_in + d_out) < d_in * d_out)
                   : (rl * (d_in + d_out) * 5L < d_in * d_out * 4L);
    if (svd_faster &&
        static_cast<std::size_t>(r) < svd.U_by_rank.size() &&
        svd.SVT.dim(0) >= rl && svd.SVT.dim(1) == d_out) {
        op.use_svd = true;
        op.rank = rl;
        op.u_ptr = svd.U_by_rank[static_cast<std::size_t>(r)].ptr();
        op.svt_ptr = svd.SVT.ptr();
    } else {
        op.use_svd = false;
        op.dense_ptr = dense_mat.ptr();
    }
    return op;
}

/// Per-displacement 1D stencil screening & low-rank SVD selection using exact `StencilSVD1D`
/// Frobenius tail norms with symmetric `-disp/+disp` parity pairing.
///
/// Populates `out_ops[disp + kMaxDisp]` and returns `max_active_disp in [-1, m_bw]`. When
/// `min_active_bw >= 1` and `max_active_disp >= 0`, preserves `disp in [-min_active_bw, +min_active_bw]`
/// (at minimum rank 1) whenever `svd_arr[slot].normf > 0` so Theorem 3.1 boundary commutation and
/// constant-polynomial neighbor cancellation remain algebraically intact.
inline int select_active_1d_slots(
    int m_bw,
    int min_active_bw,
    long d_in,
    long d_out,
    const std::array<Tensor<double>, StencilEntry1D::kNumDisp>& dense_arr,
    const std::array<StencilSVD1D, StencilEntry1D::kNumDisp>& svd_arr,
    double abs_weight,
    double term_svd_tol,
    std::array<SlotOp1D, StencilEntry1D::kNumDisp>& out_ops) noexcept {
    out_ops = {};
    if (m_bw < 0 || !(abs_weight > 0.0)) return -1;

    double max_1d = 0.0;
    for (int d = -m_bw; d <= m_bw; ++d) {
        const std::size_t slot = static_cast<std::size_t>(d + StencilEntry1D::kMaxDisp);
        max_1d = std::max(max_1d, svd_arr[slot].normf);
    }
    if (!(max_1d > 0.0)) return -1;

    const double tol_1d = term_svd_tol / (3.0 * abs_weight * max_1d * max_1d);
    int max_active = -1;

    constexpr std::size_t slot0 = static_cast<std::size_t>(StencilEntry1D::kMaxDisp);
    const int r0 = svd_arr[slot0].required_rank(tol_1d);
    if (r0 > 0) {
        out_ops[slot0] = make_slot_op_1d(r0, d_in, d_out, dense_arr[slot0], svd_arr[slot0]);
        if (out_ops[slot0].active) max_active = 0;
    }

    for (int d = 1; d <= m_bw; ++d) {
        const std::size_t slot_neg = static_cast<std::size_t>(-d + StencilEntry1D::kMaxDisp);
        const std::size_t slot_pos = static_cast<std::size_t>(+d + StencilEntry1D::kMaxDisp);
        const int r = std::max(svd_arr[slot_neg].required_rank(tol_1d),
                               svd_arr[slot_pos].required_rank(tol_1d));
        if (r > 0) {
            out_ops[slot_neg] = make_slot_op_1d(r, d_in, d_out, dense_arr[slot_neg], svd_arr[slot_neg]);
            out_ops[slot_pos] = make_slot_op_1d(r, d_in, d_out, dense_arr[slot_pos], svd_arr[slot_pos]);
            if (out_ops[slot_neg].active || out_ops[slot_pos].active) {
                max_active = d;
            }
        }
    }

    if (max_active >= 0 && max_active < min_active_bw) {
        const int req_hi = std::min(m_bw, min_active_bw);
        for (int d = 1; d <= req_hi; ++d) {
            const std::size_t slot_neg = static_cast<std::size_t>(-d + StencilEntry1D::kMaxDisp);
            const std::size_t slot_pos = static_cast<std::size_t>(+d + StencilEntry1D::kMaxDisp);
            if (!out_ops[slot_neg].active && svd_arr[slot_neg].normf > 0.0) {
                out_ops[slot_neg] = make_slot_op_1d(1, d_in, d_out, dense_arr[slot_neg], svd_arr[slot_neg]);
            }
            if (!out_ops[slot_pos].active && svd_arr[slot_pos].normf > 0.0) {
                out_ops[slot_pos] = make_slot_op_1d(1, d_in, d_out, dense_arr[slot_pos], svd_arr[slot_pos]);
            }
            if (out_ops[slot_neg].active || out_ops[slot_pos].active) {
                max_active = std::max(max_active, d);
            }
        }
    }

    return max_active;
}

/// Reusable contiguous scratch buffers for `execute_sumfact_3d` and self-box fast paths.
template <typename T>
struct SumFactScratch3D {
    std::vector<T> buf_x;
    std::vector<T> buf_xy;
    std::vector<double> z_op_buf;
    std::array<std::vector<T>, kMaxParallelChunks> tmp_blocks;
    std::array<std::vector<T>, kMaxParallelChunks> svd_blocks;
    std::array<std::vector<T>, kMaxParallelChunks> svd_t_blocks;
    std::array<Tensor<T>, kMaxParallelChunks> self_work1;
    std::array<Tensor<T>, kMaxParallelChunks> self_work2;

    void ensure_self_cubes(long d) {
        const long dims[3] = {d, d, d};
        for (std::size_t cid = 0; cid < kMaxParallelChunks; ++cid) {
            if (self_work1[cid].ndim() != 3 || self_work1[cid].dim(0) < d ||
                self_work1[cid].dim(1) < d || self_work1[cid].dim(2) < d) {
                self_work1[cid] = Tensor<T>(3L, dims, false);
                self_work2[cid] = Tensor<T>(3L, dims, false);
            }
        }
    }
};

/// Contract one spatial stencil edge `target(i, j) = sum_q src(q, i) * A_disp(q, j)` using
/// 100% Tier-1 `s_dense_table` `mTxmq` (`dimi <= 400, dimj <= 24, dimk <= 24`):
///   - `sop.use_svd`: `mTxmq(dimi, rank, d_in)` -> L1 `fast_transpose(dimi, rank)` (elided when `rank == 1`)
///                    -> `mTxmq(dimi, d_out, rank)`.
///   - `!sop.use_svd`: `mTxmq(dimi, d_out, d_in)`.
template <typename T>
inline void apply_slot_edge_1d(long dimi,
                               long d_in,
                               long d_out,
                               long out_sz,
                               const SlotOp1D& sop,
                               const T* src,
                               T* MADNESS_RESTRICT dst,
                               T* MADNESS_RESTRICT tmp_out,
                               T* MADNESS_RESTRICT tmp_svd,
                               T* MADNESS_RESTRICT tmp_svd_T,
                               bool& wrote) {
    if (!sop.active) return;
    T* MADNESS_RESTRICT target = wrote ? tmp_out : dst;
    if (sop.use_svd) {
        // tmp_svd[i * rank + p] = sum_{q=0}^{d_in-1} src[q * dimi + i] * u_ptr[q * rank + p]
        mTxmq(dimi, sop.rank, d_in, tmp_svd, src, sop.u_ptr);
        if (sop.rank == 1L) {
            // Shape `(dimi, 1)` is identical in memory to `(1, dimi)` — zero transpose needed!
            mTxmq(dimi, d_out, 1L, target, tmp_svd, sop.svt_ptr);
        } else {
            fast_transpose(dimi, sop.rank, tmp_svd, tmp_svd_T);
            mTxmq(dimi, d_out, sop.rank, target, tmp_svd_T, sop.svt_ptr);
        }
    } else {
        mTxmq(dimi, d_out, d_in, target, src, sop.dense_ptr);
    }
    if (!wrote) {
        wrote = true;
    } else {
        aligned_add<T, T>(out_sz, dst, tmp_out);
    }
}

/// Execute one 3D spatial-offset + cyclic `mTxmq` + `aligned_add<T, T>` sum-factorized sweep:
///   `tgt_out[t] += sum_{s} (A(tx - sx) \otimes A(ty - sy) \otimes (weight * A(tz - sz))) src[s]`
/// for arbitrary `(d_in, d_out)` (`(k, k)`, `(2k, k)`, `(k, 2k)`, or `(2k, 2k)`) with 100% Tier-1
/// `s_dense_table` `mTxmq` hits, per-displacement 1D SVD acceleration, and direct weighted Z-sweep
/// accumulation into `tgt_out_ptrs[ti]` (eliminating `uz` staging and `aligned_axpy` passes).
template <typename T>
inline void execute_sumfact_3d(
    const SumFactGraph3D& graph,
    long d_in,
    long d_out,
    const std::array<SlotOp1D, StencilEntry1D::kNumDisp>& slot_ops,
    double weight,
    const std::vector<const T*>& src_ptrs,
    const std::vector<T*>& tgt_out_ptrs,
    SumFactScratch3D<T>& scratch) {
    if (graph.empty() || weight == 0.0) return;

    const long dimi_x = d_in * d_in;
    const long size_x = dimi_x * d_out;

    const long dimi_y = d_in * d_out;
    const long size_y = dimi_y * d_out;

    const long dimi_z = d_out * d_out;
    const long size_z = dimi_z * d_out;

    const std::size_t need_x = graph.num_x * static_cast<std::size_t>(size_x);
    const std::size_t need_xy = graph.num_xy * static_cast<std::size_t>(size_y);
    const std::size_t max_blk = static_cast<std::size_t>(std::max({size_x, size_y, size_z}));
    const std::size_t max_svd_blk = static_cast<std::size_t>(
        std::max({dimi_x, dimi_y, dimi_z}) * std::max(d_in, d_out));

    if (scratch.buf_x.size() < need_x) scratch.buf_x.resize(need_x);
    if (scratch.buf_xy.size() < need_xy) scratch.buf_xy.resize(need_xy);
    for (std::size_t cid = 0; cid < kMaxParallelChunks; ++cid) {
        if (scratch.tmp_blocks[cid].size() < max_blk) {
            scratch.tmp_blocks[cid].resize(max_blk);
        }
        if (scratch.svd_blocks[cid].size() < max_svd_blk) {
            scratch.svd_blocks[cid].resize(max_svd_blk);
        }
        if (scratch.svd_t_blocks[cid].size() < max_svd_blk) {
            scratch.svd_t_blocks[cid].resize(max_svd_blk);
        }
    }

    // Fold `weight` into a tiny L1-resident copy of the active 1D Z operators (`<= 5` matrices)
    // so Sweep 3 accumulates directly into `tgt_out_ptrs[ti]` with zero `uz` buffer / `aligned_axpy`.
    const std::size_t slot_stride = static_cast<std::size_t>(std::max(d_in, d_out) * d_out);
    constexpr std::size_t kStackSlotStride = 576;
    alignas(64) std::array<double, StencilEntry1D::kNumDisp * kStackSlotStride> z_stack_buf;
    double* z_buf_base = nullptr;
    if (slot_stride <= kStackSlotStride) {
        z_buf_base = z_stack_buf.data();
    } else {
        const std::size_t need_z_op = StencilEntry1D::kNumDisp * slot_stride;
        if (scratch.z_op_buf.size() < need_z_op) scratch.z_op_buf.resize(need_z_op);
        z_buf_base = scratch.z_op_buf.data();
    }

    std::array<SlotOp1D, StencilEntry1D::kNumDisp> z_slot_ops = slot_ops;
    for (std::size_t slot = 0; slot < StencilEntry1D::kNumDisp; ++slot) {
        SlotOp1D& zop = z_slot_ops[slot];
        if (!zop.active) continue;
        double* MADNESS_RESTRICT dst_op = z_buf_base + slot * slot_stride;
        if (zop.use_svd) {
            const std::size_t n_el = static_cast<std::size_t>(zop.rank * d_out);
            for (std::size_t i = 0; i < n_el; ++i) {
                dst_op[i] = weight * zop.svt_ptr[i];
            }
            zop.svt_ptr = dst_op;
        } else {
            const std::size_t n_el = static_cast<std::size_t>(d_in * d_out);
            for (std::size_t i = 0; i < n_el; ++i) {
                dst_op[i] = weight * zop.dense_ptr[i];
            }
            zop.dense_ptr = dst_op;
        }
    }

    T* const ux_base = scratch.buf_x.data();
    T* const uxy_base = scratch.buf_xy.data();

    // Sweep 1 (X): (qx, qy, qz) -> (qy, qz, ix)
    parallel_for_chunks(
        graph.num_x, kSumFactSweepGrain, [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
            T* MADNESS_RESTRICT tmp = scratch.tmp_blocks[cid].data();
            T* MADNESS_RESTRICT tmp_svd = scratch.svd_blocks[cid].data();
            T* MADNESS_RESTRICT tmp_svd_T = scratch.svd_t_blocks[cid].data();
            for (std::size_t xi = b; xi < e_hi; ++xi) {
                T* MADNESS_RESTRICT dst = ux_base + xi * static_cast<std::size_t>(size_x);
                const std::uint32_t e_begin = graph.x_row_ptr[xi];
                const std::uint32_t e_end = graph.x_row_ptr[xi + 1];
                bool wrote = false;
                for (std::uint32_t e = e_begin; e < e_end; ++e) {
                    const StencilEdge1D edge = graph.x_edges[e];
                    const SlotOp1D& sop = slot_ops[static_cast<std::size_t>(
                        static_cast<int>(edge.disp) + StencilEntry1D::kMaxDisp)];
                    if (!sop.active) continue;
                    apply_slot_edge_1d(dimi_x,
                                       d_in,
                                       d_out,
                                       size_x,
                                       sop,
                                       src_ptrs[edge.src_idx],
                                       dst,
                                       tmp,
                                       tmp_svd,
                                       tmp_svd_T,
                                       wrote);
                }
                if (!wrote) {
                    aligned_zero(size_x, dst);
                }
            }
        });

    // Sweep 2 (Y): (qy, qz, ix) -> (qz, ix, iy)
    parallel_for_chunks(
        graph.num_xy, kSumFactSweepGrain, [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
            T* MADNESS_RESTRICT tmp = scratch.tmp_blocks[cid].data();
            T* MADNESS_RESTRICT tmp_svd = scratch.svd_blocks[cid].data();
            T* MADNESS_RESTRICT tmp_svd_T = scratch.svd_t_blocks[cid].data();
            for (std::size_t xyi = b; xyi < e_hi; ++xyi) {
                T* MADNESS_RESTRICT dst = uxy_base + xyi * static_cast<std::size_t>(size_y);
                const std::uint32_t e_begin = graph.xy_row_ptr[xyi];
                const std::uint32_t e_end = graph.xy_row_ptr[xyi + 1];
                bool wrote = false;
                for (std::uint32_t e = e_begin; e < e_end; ++e) {
                    const StencilEdge1D edge = graph.xy_edges[e];
                    const SlotOp1D& sop = slot_ops[static_cast<std::size_t>(
                        static_cast<int>(edge.disp) + StencilEntry1D::kMaxDisp)];
                    if (!sop.active) continue;
                    const T* src =
                        ux_base + static_cast<std::size_t>(edge.src_idx) * static_cast<std::size_t>(size_x);
                    apply_slot_edge_1d(
                        dimi_y, d_in, d_out, size_y, sop, src, dst, tmp, tmp_svd, tmp_svd_T, wrote);
                }
                if (!wrote) {
                    aligned_zero(size_y, dst);
                }
            }
        });

    // Sweep 3 (Z): (qz, ix, iy) -> (ix, iy, iz) accumulated directly into `tgt_out_ptrs[ti]`
    const std::size_t num_active_tgt = graph.active_tgt_idx.size();
    parallel_for_chunks(
        num_active_tgt, kSumFactSweepGrain, [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
            T* MADNESS_RESTRICT tmp = scratch.tmp_blocks[cid].data();
            T* MADNESS_RESTRICT tmp_svd = scratch.svd_blocks[cid].data();
            T* MADNESS_RESTRICT tmp_svd_T = scratch.svd_t_blocks[cid].data();
            for (std::size_t ai = b; ai < e_hi; ++ai) {
                const std::uint32_t ti = graph.active_tgt_idx[ai];
                T* MADNESS_RESTRICT dst = tgt_out_ptrs[ti];
                const std::uint32_t e_begin = graph.z_row_ptr[ai];
                const std::uint32_t e_end = graph.z_row_ptr[ai + 1];
                bool wrote = true;
                for (std::uint32_t e = e_begin; e < e_end; ++e) {
                    const StencilEdge1D edge = graph.z_edges[e];
                    const SlotOp1D& sop = z_slot_ops[static_cast<std::size_t>(
                        static_cast<int>(edge.disp) + StencilEntry1D::kMaxDisp)];
                    if (!sop.active) continue;
                    const T* src = uxy_base +
                                   static_cast<std::size_t>(edge.src_idx) * static_cast<std::size_t>(size_y);
                    apply_slot_edge_1d(
                        dimi_z, d_in, d_out, size_z, sop, src, dst, tmp, tmp_svd, tmp_svd_T, wrote);
                }
            }
        });
}

/// Self-box (`m_bw == 0`, `disp = (0, 0, 0)`) fast path for Pass 4b-1 (`2k -> k`):
/// uses 6-step Tier-1 `s_dense_table` `U_by_rank[r]` + `SVT` Tucker contractions (`(2k)^3 -> r^3 -> k^3`)
/// when low-rank `r` saves FLOPs, or `general_fast_transform` + `aligned_axpy` otherwise.
template <typename T>
inline void execute_self_box_2k_to_k_3d(
    const std::vector<std::pair<std::size_t, std::size_t>>& matches,
    const std::vector<const Tensor<T>*>& src_tensors,
    const std::vector<T*>& tgt_ptrs,
    long kk,
    const StencilEntry1D& entry,
    double weight,
    double term_svd_tol,
    SumFactScratch3D<T>& scratch) {
    if (matches.empty() || weight == 0.0) return;

    constexpr std::size_t slot0 = static_cast<std::size_t>(StencilEntry1D::kMaxDisp);
    const StencilSVD1D& svd0 = entry.svd_R_2k_to_k[slot0];
    const double n0 = svd0.normf;
    if (!(n0 > 0.0)) return;

    const long twok = 2L * kk;
    const long kk3 = kk * kk * kk;
    const double tol_1d = term_svd_tol / (3.0 * std::abs(weight) * n0 * n0);
    const long r = static_cast<long>(svd0.required_rank(tol_1d));
    if (r <= 0) return;

    const long flops_dense = twok * twok * twok * kk + twok * twok * kk * kk + twok * kk * kk * kk;
    const long flops_svd =
        (twok * twok * twok * r + twok * twok * r * r + twok * r * r * r) +
        (r * r * r * kk + r * r * kk * kk + r * kk * kk * kk);
    const bool use_svd =
        (r < kk) &&
        (static_cast<std::size_t>(r) < svd0.U_by_rank.size()) &&
        (flops_svd * 20L < flops_dense * 17L);

    scratch.ensure_self_cubes(twok);
    const Tensor<double>& mat_0 = entry.R_2k_to_k[slot0];
    const Tensor<double> c3[3] = {mat_0, mat_0, mat_0};
    const double* const u_ptr = use_svd ? svd0.U_by_rank[static_cast<std::size_t>(r)].ptr() : nullptr;
    const double* const svt_ptr = use_svd ? svd0.SVT.ptr() : nullptr;

    parallel_for_chunks(
        matches.size(), kSumFactSweepGrain, [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
            Tensor<T>& w_res = scratch.self_work1[cid];
            Tensor<T>& w_tmp = scratch.self_work2[cid];
            T* MADNESS_RESTRICT w1 = w_res.ptr();
            T* MADNESS_RESTRICT w2 = w_tmp.ptr();
            for (std::size_t idx = b; idx < e_hi; ++idx) {
                const auto [s_i, t_i] = matches[idx];
                if (use_svd) {
                    mTxmq(twok * twok, r, twok, w1, src_tensors[s_i]->ptr(), u_ptr);
                    mTxmq(twok * r, r, twok, w2, w1, u_ptr);
                    mTxmq(r * r, r, twok, w1, w2, u_ptr);
                    mTxmq(r * r, kk, r, w2, w1, svt_ptr);
                    mTxmq(r * kk, kk, r, w1, w2, svt_ptr);
                    mTxmq(kk * kk, kk, r, w2, w1, svt_ptr);
                    aligned_axpy(kk3, tgt_ptrs[t_i], w2, weight);
                } else {
                    general_fast_transform(*src_tensors[s_i], c3, w_res, w_tmp);
                    aligned_axpy(kk3, tgt_ptrs[t_i], w_res.ptr(), weight);
                }
            }
        });
}

/// Self-box (`m_bw == 0`, `disp = (0, 0, 0)`) fast path for Pass 4b-2 internal (`2k -> 2k`):
/// uses 6-step Tier-1 `s_dense_table` `U_by_rank[r]` + `SVT` Tucker contractions (`(2k)^3 -> r^3 -> (2k)^3`)
/// with dense-packed `ldb == dimj` when `1 <= r < 0.65 * (2k)`, or `fast_transform` + `aligned_axpy` otherwise.
template <typename T>
inline void execute_self_box_2k_to_2k_3d(
    const std::vector<std::pair<std::size_t, std::size_t>>& matches,
    const std::vector<const Tensor<T>*>& src_tensors,
    const std::vector<T*>& tgt_ptrs,
    long twok,
    const StencilEntry1D& entry,
    double weight,
    double term_svd_tol,
    SumFactScratch3D<T>& scratch) {
    if (matches.empty() || weight == 0.0) return;

    constexpr std::size_t slot0 = static_cast<std::size_t>(StencilEntry1D::kMaxDisp);
    const StencilSVD1D& svd0 = entry.svd_R[slot0];
    const double n0 = svd0.normf;
    if (!(n0 > 0.0)) return;

    const double tol_1d = term_svd_tol / (3.0 * std::abs(weight) * n0 * n0);
    const long r = static_cast<long>(svd0.required_rank(tol_1d));
    if (r <= 0) return;

    const long break_even = std::max(1L, static_cast<long>(0.65 * static_cast<double>(twok)));
    const bool use_svd =
        (r < break_even) && (static_cast<std::size_t>(r) < svd0.U_by_rank.size());

    const long twok3 = twok * twok * twok;
    scratch.ensure_self_cubes(twok);
    const Tensor<double>& mat_R0 = entry.R[slot0];
    const double* const u_ptr = use_svd ? svd0.U_by_rank[static_cast<std::size_t>(r)].ptr() : nullptr;
    const double* const svt_ptr = use_svd ? svd0.SVT.ptr() : nullptr;

    parallel_for_chunks(
        matches.size(), kSumFactSweepGrain, [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
            Tensor<T>& work1 = scratch.self_work1[cid];
            Tensor<T>& work2 = scratch.self_work2[cid];
            T* MADNESS_RESTRICT w1 = work1.ptr();
            T* MADNESS_RESTRICT w2 = work2.ptr();
            for (std::size_t idx = b; idx < e_hi; ++idx) {
                const auto [s_i, t_i] = matches[idx];
                if (use_svd) {
                    mTxmq(twok * twok, r, twok, w1, src_tensors[s_i]->ptr(), u_ptr);
                    mTxmq(twok * r, r, twok, w2, w1, u_ptr);
                    mTxmq(r * r, r, twok, w1, w2, u_ptr);
                    mTxmq(r * r, twok, r, w2, w1, svt_ptr);
                    mTxmq(r * twok, twok, r, w1, w2, svt_ptr);
                    mTxmq(twok * twok, twok, r, w2, w1, svt_ptr);
                    aligned_axpy(twok3, tgt_ptrs[t_i], w2, weight);
                } else {
                    fast_transform(*src_tensors[s_i], mat_R0, work1, work2);
                    aligned_axpy(twok3, tgt_ptrs[t_i], work1.ptr(), weight);
                }
            }
        });
}

/// Apply `(A \otimes A \otimes A) src` on a single root box `(k, k, k)` and accumulate `weight * (...)` into `dst`.
template <typename T>
inline void apply_single_box_1d_kronecker_3d(long kk3,
                                             const Tensor<double>& mat_1d,
                                             double weight,
                                             const Tensor<T>& src,
                                             Tensor<T>& dst,
                                             Tensor<T>& work1,
                                             Tensor<T>& work2) {
    if (weight == 0.0) return;
    fast_transform(src, mat_1d, work2, work1);
    aligned_axpy(kk3, dst.ptr(), work2.ptr(), weight);
}

/// Zero the scaling octant `s0 = [0..k)^3` of a `(2k)^3` NS tensor via `Slice` so only the 7 wavelet octants remain.
template <typename T>
inline void zero_s0_octant(Tensor<T>& ns_tensor, int k) {
    const long kk = static_cast<long>(k);
    const long twok = 2L * kk;
    if (ns_tensor.size() != twok * twok * twok) return;
    const Slice s0(0, kk - 1);
    ns_tensor(s0, s0, s0) = T(0);
}

struct RawTermRef {
    std::size_t term_idx = 0;
    const StencilEntry1D* entry = nullptr;
    double weight_scaled = 0.0;
};

struct ActiveBaseTerm {
    std::size_t term_idx = 0;
    const StencilEntry1D* entry = nullptr;
    double weight_scaled = 0.0;
    int bw_T = 1;
};

struct ActiveGammaTerm {
    std::size_t term_idx = 0;
    const StencilEntry1D* entry = nullptr;
    double weight_scaled = 0.0;
    int bw_4b1 = -1;
    int bw_leaf = -1;
    int bw_int = -1;
};

} // namespace windowed_sumfact_detail

/// Stage 1a helper: prepare an input `Function<T, 3>` in `TreeState::nonstandard_with_leaves`,
/// returning `true` iff `f_work` entered in `TreeState::compressed` (so `restore_windowed_ns_input_3d`
/// can restore it via zero-FLOP `standard()`).
template <typename T>
bool prepare_windowed_ns_input_3d(Function<T, 3>& f_work, bool fence = true) {
    MADNESS_CHECK_THROW(f_work.is_initialized(),
                        "prepare_windowed_ns_input_3d: input Function must be initialized");
    MADNESS_CHECK_THROW(!f_work.is_on_demand(),
                        "prepare_windowed_ns_input_3d: on-demand input Function is not supported");
    mTxmq_init();

    const TreeState init_state = f_work.get_impl()->get_tree_state();
    const bool was_compressed = (init_state == compressed);
    if (init_state != reconstructed && init_state != nonstandard_with_leaves) {
        f_work.reconstruct(true);
    }
    if (f_work.get_impl()->get_tree_state() != nonstandard_with_leaves) {
        f_work.make_nonstandard(true, fence);
    }
    return was_compressed;
}

/// Stage 5b helper: restore `f_work` after `apply_windowed_ns_3d_prepared` using `clear(fence)` when
/// `destructive_input == true`, zero-FLOP `standard(fence)` when `was_compressed == true`, or `reconstruct(fence)`.
template <typename T>
void restore_windowed_ns_input_3d(Function<T, 3>& f_work,
                                  bool was_compressed,
                                  bool destructive_input,
                                  bool fence = true) {
    if (!f_work.is_initialized()) return;
    if (destructive_input) {
        f_work.clear(fence);
    } else if (was_compressed) {
        f_work.standard(fence);
    } else {
        f_work.reconstruct(fence);
    }
}

/// Core Stage 1b through Stage 5a Windowed NS execution on an already-prepared `f_ns`
/// (`TreeState::nonstandard_with_leaves`) and already-allocated empty `result`
/// (`result.get_impl()->get_coeffs().process_pending()` completed by caller),
/// leaving `result` in `TreeState::nonstandard_after_apply`.
///
/// When `world.size() == 1`, performs zero `WorldObject` registrations and zero `world.gop`
/// collectives, making it safe to execute concurrently across disjoint `(f_ns[i], result[i])`
/// pairs on `ThreadPool` worker threads under `WindowedBatchTaskScope`.
template <typename T>
void apply_windowed_ns_3d_prepared(World& world,
                                   const Function<T, 3>& f_ns,
                                   Function<T, 3>& result,
                                   const WindowedQuadraturePlan& plan,
                                   const WindowedStencilCache1D& stencil_cache,
                                   double thresh,
                                   bool dynamic_refine = true,
                                   double dynamic_refine_factor = 1.0) {
    using namespace windowed_sumfact_detail;

    MADNESS_CHECK(f_ns.is_initialized() && result.is_initialized());
    MADNESS_CHECK(f_ns.get_impl()->get_tree_state() == nonstandard_with_leaves);
    MADNESS_CHECK(f_ns.k() == plan.k() && f_ns.k() == stencil_cache.k());

    const int k = f_ns.k();
    const long kk = static_cast<long>(k);
    const long kk3 = kk * kk * kk;
    const long twok = 2L * kk;
    const long twok3 = twok * twok * twok;
    const double eff_thresh = (thresh > 0.0) ? thresh : f_ns.thresh();
    const double eff_tol = (eff_thresh > 0.0) ? std::min(plan.eps(), eff_thresh) : plan.eps();
    const double eff_refine_fac = (dynamic_refine_factor > 0.0) ? dynamic_refine_factor : 1.0;

    const int n_max_tree = (world.size() == 1) ? static_cast<int>(f_ns.max_local_depth())
                                               : static_cast<int>(f_ns.max_depth());
    MADNESS_CHECK_THROW(
        n_max_tree <= plan.n_max(),
        "apply_windowed_ns_3d_prepared: input function depth exceeds WindowedQuadraturePlan::n_max()");

    int max_refine_level = FunctionDefaults<3>::get_max_refine_level();
    if (max_refine_level <= 0) {
        max_refine_level = plan.n_max();
    } else {
        max_refine_level = std::min(max_refine_level, plan.n_max());
    }
    max_refine_level = std::min(max_refine_level, static_cast<int>(MAXLEVEL) - 1);

    // Pre-bucket candidate quadrature terms per tree level `j in [0, n_max_tree]`.
    std::vector<std::vector<RawTermRef>> raw_T_terms(static_cast<std::size_t>(n_max_tree + 1));
    std::vector<std::vector<RawTermRef>> raw_gamma_terms(static_cast<std::size_t>(n_max_tree + 1));
    std::vector<std::pair<std::size_t, double>> root_ultra_diffuse_terms;

    const auto& all_terms = plan.terms();
    for (std::size_t t_idx = 0; t_idx < all_terms.size(); ++t_idx) {
        const WindowedQuadratureTerm& term = all_terms[t_idx];
        if (term.weight_scaled == 0.0) continue;

        if (term.n_diff < 0) {
            root_ultra_diffuse_terms.emplace_back(t_idx, term.weight_scaled);
            continue;
        }

        if (term.j0 <= n_max_tree) {
            const StencilEntry1D& e_j0 = stencil_cache.entry_for_term(t_idx, term, term.j0);
            raw_T_terms[static_cast<std::size_t>(term.j0)].push_back(
                RawTermRef{t_idx, &e_j0, term.weight_scaled});
        }

        const int j_gamma_top = std::min(n_max_tree - 1, term.n_diff);
        for (int j = term.j0; j <= j_gamma_top; ++j) {
            const StencilEntry1D& e_j = stencil_cache.entry_for_term(t_idx, term, j);
            raw_gamma_terms[static_cast<std::size_t>(j)].push_back(
                RawTermRef{t_idx, &e_j, term.weight_scaled});
        }
    }

    const auto* f_impl = f_ns.get_impl().get();
    std::vector<std::vector<ActiveBaseTerm>> active_T_terms(static_cast<std::size_t>(n_max_tree + 1));
    std::vector<std::vector<ActiveGammaTerm>> active_gamma_terms(static_cast<std::size_t>(n_max_tree + 1));

    // Dynamic bandwidth classification:
    // - When `dynamic_refine == false`, uses unit reference norms (`1.0`) so operator bandwidths
    //   depend solely on `(plan, eff_tol, tree_topology)` and shared-tree linearity holds to ~1e-15.
    // - When `dynamic_refine == true`, uses global level norms and Legacy MADNESS's shell-2
    //   threshold `loc_tol_j / (14.14 * sqrt(N_terms))` while keeping shell-0/1 at `0.01 * loc_tol_j / sqrt(N)`.
    auto radii_builder = [&](const std::vector<WindowedGlobalLevelStats>& gstats,
                             std::vector<WindowedLevelRadii>& out_radii) {
        out_radii.assign(static_cast<std::size_t>(n_max_tree + 1), WindowedLevelRadii{});

        for (int j = 0; j <= n_max_tree; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const WindowedGlobalLevelStats& st = gstats[sj];
            const double loc_tol_j = f_impl->truncate_tol(eff_tol, Key<3>(static_cast<Level>(j)));

            double sum_sh2_s_T = 0.0;
            double max_sh2_s_T = 0.0;
            double sum_sh2_4b1 = 0.0;
            double max_sh2_4b1 = 0.0;
            double sum_sh2_s_leaf = 0.0;
            double max_sh2_s_leaf = 0.0;
            double sum_sh2_refine = 0.0;
            double max_sh2_refine = 0.0;

            // Stage 4(a) active bandwidth classification
            const auto& rT = raw_T_terms[sj];
            if (!rT.empty()) {
                if (j == 0) {
                    if (!dynamic_refine || st.max_s_real > 0.0) {
                        for (const RawTermRef& rt : rT) {
                            active_T_terms[sj].push_back(
                                ActiveBaseTerm{rt.term_idx, rt.entry, rt.weight_scaled, 0});
                        }
                    }
                } else if (!dynamic_refine || st.max_s_source > 0.0) {
                    const double sqrt_nT =
                        std::sqrt(static_cast<double>(std::max<std::size_t>(1, rT.size())));
                    const double tol_nn_T = 0.01 * loc_tol_j / sqrt_nT;
                    const double tol_sh2_T =
                        dynamic_refine ? (loc_tol_j / (14.14 * sqrt_nT)) : (0.1 * tol_nn_T);
                    const double ref_s_src = dynamic_refine ? st.max_s_source : 1.0;

                    for (const RawTermRef& rt : rT) {
                        const double aw = std::abs(rt.weight_scaled);
                        const int bw_T = rt.entry->select_bw_T(aw, ref_s_src, tol_nn_T, tol_sh2_T);
                        if (bw_T >= 1) {
                            active_T_terms[sj].push_back(
                                ActiveBaseTerm{rt.term_idx, rt.entry, rt.weight_scaled, bw_T});
                            out_radii[sj].has_T_terms = true;
                            out_radii[sj].m_max_T = std::max(out_radii[sj].m_max_T, bw_T);
                            if (rt.entry->m_max_T >= 2) {
                                const double op2 = aw * rt.entry->shell_max_T[2];
                                sum_sh2_s_T += op2;
                                max_sh2_s_T = std::max(max_sh2_s_T, sqrt_nT * op2);
                            }
                        }
                    }
                }
            }

            // Stage 4(b) active bandwidth classification (strictly skipped if level has no internal nodes)
            if (j < n_max_tree && st.has_internal && !raw_gamma_terms[sj].empty()) {
                const auto& rG = raw_gamma_terms[sj];
                const double sqrt_nG =
                    std::sqrt(static_cast<double>(std::max<std::size_t>(1, rG.size())));
                const double tol_nn_G = 0.01 * loc_tol_j / sqrt_nG;
                const double tol_sh2_G =
                    dynamic_refine ? (loc_tol_j / (14.14 * sqrt_nG)) : (0.1 * tol_nn_G);
                const double ref_s_src = dynamic_refine ? st.max_s_source : 1.0;
                const double ref_s_real = dynamic_refine ? st.max_s_real : 1.0;
                const double ref_d_int = dynamic_refine ? st.max_d_int : 1.0;

                for (const RawTermRef& rg : rG) {
                    const double aw = std::abs(rg.weight_scaled);
                    const int bw_4b1 =
                        rg.entry->select_bw_4b1(aw, ref_d_int, tol_nn_G, tol_sh2_G);
                    const int bw_leaf =
                        rg.entry->select_bw_4b2_leaf(aw, ref_s_src, tol_nn_G, tol_sh2_G);
                    const int bw_int = rg.entry->select_bw_4b2_int(
                        aw, ref_s_real, ref_d_int, tol_nn_G, tol_sh2_G, bw_leaf);

                    if (bw_4b1 >= 0) {
                        out_radii[sj].has_4b1_terms = true;
                        out_radii[sj].m_max_4b1 = std::max(out_radii[sj].m_max_4b1, bw_4b1);
                    }
                    if (bw_leaf >= 1) {
                        out_radii[sj].has_4b2_leaf_terms = true;
                        out_radii[sj].m_max_4b2_leaf = std::max(out_radii[sj].m_max_4b2_leaf, bw_leaf);
                    }
                    if (bw_int >= 0) {
                        out_radii[sj].has_4b2_int_terms = true;
                        out_radii[sj].m_max_4b2_int = std::max(out_radii[sj].m_max_4b2_int, bw_int);
                    }
                    if (bw_4b1 >= 0 || bw_leaf >= 1 || bw_int >= 0) {
                        active_gamma_terms[sj].push_back(ActiveGammaTerm{
                            rg.term_idx, rg.entry, rg.weight_scaled, bw_4b1, bw_leaf, bw_int});
                        if (rg.entry->m_max_ns >= 2) {
                            const double op_4b1 = aw * rg.entry->shell_max_4b1[2];
                            const double op_leaf = aw * rg.entry->shell_max_4b2_leaf[2];
                            const double op_ref =
                                aw * std::max(rg.entry->shell_max_4b2_leaf[2],
                                              rg.entry->shell_max_4b2_int[2]);
                            sum_sh2_4b1 += op_4b1;
                            max_sh2_4b1 = std::max(max_sh2_4b1, sqrt_nG * op_4b1);
                            sum_sh2_s_leaf += op_leaf;
                            max_sh2_s_leaf = std::max(max_sh2_s_leaf, sqrt_nG * op_leaf);
                            sum_sh2_refine += op_ref;
                            max_sh2_refine = std::max(max_sh2_refine, sqrt_nG * op_ref);
                        }
                    }
                }
            }

            if (dynamic_refine) {
                const double sh2_tol = loc_tol_j / 14.14;
                const double eff_s_op = std::max(
                    {1e-30, 2.5 * (sum_sh2_s_T + sum_sh2_s_leaf), max_sh2_s_T, max_sh2_s_leaf});
                const double eff_4b1_op = std::max({1e-30, sum_sh2_4b1, max_sh2_4b1});
                const double eff_ref_op = std::max({1e-30, sum_sh2_refine, max_sh2_refine});
                out_radii[sj].sh2_s_thresh = sh2_tol / eff_s_op;
                out_radii[sj].sh2_4b1_thresh = sh2_tol / eff_4b1_op;
                out_radii[sj].sh2_refine_thresh = sh2_tol / eff_ref_op;
            }
        }
    };

    // Inner scope guarantees `WindowedNSHaloWorkspace`'s internal `WorldObject` (when `world.size() > 1`)
    // is fenced & destroyed before `apply_windowed_ns_3d_prepared` returns.
    {
        WindowedNSHaloWorkspace<T> ws(
            world, f_ns, n_max_tree, dynamic_refine, max_refine_level, radii_builder);

        const LegendreMomentOp& mop = plan.moment_op();
        const auto& gstats = ws.global_stats();
        const auto& radii = ws.radii();

        const long k_dims[3] = {kk, kk, kk};
        Tensor<T> work1_k3(3L, k_dims, false);
        Tensor<T> work2_k3(3L, k_dims, false);
        Tensor<double> mt_mat(kk, kk);
        SumFactScratch3D<T> sf_scratch;

        // =====================================================================
        // Stage 2: Cross-Level Flattened Gauss-Christoffel Spiky Moment Pass
        // =====================================================================
        {
            struct GCLeafRef {
                int level = 0;
                WindowedOwnedRealNode<T>* rn = nullptr;
            };

            std::size_t total_real_est = 0;
            for (int j = 0; j <= n_max_tree; ++j) {
                total_real_est += ws.level_data(j).owned_real.size();
            }

            std::vector<GCLeafRef> gc_leaves;
            gc_leaves.reserve(total_real_est);
            for (int j = 0; j <= n_max_tree; ++j) {
                auto& ld = ws.level_data(j);
                for (auto& rn : ld.owned_real) {
                    if (!rn.is_internal && rn.s_normf > 0.0) {
                        gc_leaves.push_back(GCLeafRef{j, &rn});
                    }
                }
            }

            if (!gc_leaves.empty()) {
                std::array<Tensor<T>, kMaxParallelChunks> gc_tmp_k3;
                std::array<Tensor<T>, kMaxParallelChunks> gc_w1_k3;
                std::array<Tensor<T>, kMaxParallelChunks> gc_w2_k3;
                parallel_for_chunks(
                    gc_leaves.size(),
                    kDefaultParallelGrain,
                    [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
                        for (std::size_t idx = b; idx < e_hi; ++idx) {
                            const GCLeafRef& item = gc_leaves[idx];
                            plan.accumulate_gc_leaf_3d(item.level,
                                                       item.rn->s_in,
                                                       item.rn->s_out,
                                                       gc_tmp_k3[cid],
                                                       gc_w1_k3[cid],
                                                       gc_w2_k3[cid]);
                        }
                    });
            }
        }

        // =====================================================================
        // Stage 3: Root (j = 0) Constant Low-Tail `w_const` & `n_diff < 0` Pass
        // =====================================================================
        {
            auto& ld0 = ws.level_data(0);
            if (!ld0.owned_real.empty()) {
                auto& root_rn = ld0.owned_real[0];
                const double w_const = plan.w_const_scaled();
                if (w_const != 0.0) {
                    root_rn.s_out(0L, 0L, 0L) +=
                        static_cast<T>(w_const) * root_rn.s_in(0L, 0L, 0L);
                }
                for (const auto& [t_idx, w_scaled] : root_ultra_diffuse_terms) {
                    const Tensor<double>& root_T0 = stencil_cache.root_T0_transposed(t_idx);
                    apply_single_box_1d_kronecker_3d(
                        kk3, root_T0, w_scaled, root_rn.s_in, root_rn.s_out, work1_k3, work2_k3);
                }
            }
        }

        // =====================================================================
        // Stage 4(a): Base `T^{j_0}` Injection + Theorem 3.1 Halo Subtraction
        // =====================================================================
        constexpr std::size_t kMaxNormBuckets = 3;
        auto select_norm_bucket = [dynamic_refine](double norm, double max_M) noexcept -> std::size_t {
            if (!dynamic_refine) return 0u;
            if (norm > 0.02 * max_M) return 0u;
            if (norm > 1e-4 * max_M) return 1u;
            return 2u;
        };

        std::array<Tensor<T>, kMaxParallelChunks> h1t_tmp_k3;
        std::array<Tensor<T>, kMaxParallelChunks> h1t_w1_k3;

        for (int j = 0; j <= n_max_tree; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const auto& t_list = active_T_terms[sj];
            if (t_list.empty()) continue;

            auto& ld = ws.level_data(j);
            if (j == 0) {
                if (!ld.owned_real.empty()) {
                    auto& root_rn = ld.owned_real[0];
                    for (const ActiveBaseTerm& tr : t_list) {
                        apply_single_box_1d_kronecker_3d(kk3,
                                                         tr.entry->mat_T(0),
                                                         tr.weight_scaled,
                                                         root_rn.s_in,
                                                         root_rn.s_out,
                                                         work1_k3,
                                                         work2_k3);
                    }
                }
                continue;
            }

            // Collect sorted targets: all owned real nodes in T_j + all owned H1_T virtual nodes.
            std::vector<Key<3>> real_tkeys;
            std::vector<T*> real_tptrs;
            real_tkeys.reserve(ld.owned_real.size());
            real_tptrs.reserve(ld.owned_real.size());
            for (auto& rn : ld.owned_real) {
                real_tkeys.push_back(rn.key);
                real_tptrs.push_back(rn.s_out.ptr());
            }

            std::vector<Key<3>> virt_tkeys;
            std::vector<T*> virt_tptrs;
            std::vector<WindowedOwnedVirtNode<T>*> h1t_virts;
            virt_tkeys.reserve(ld.owned_virt.size());
            virt_tptrs.reserve(ld.owned_virt.size());
            h1t_virts.reserve(ld.owned_virt.size());
            for (auto& vn : ld.owned_virt) {
                if ((vn.flags & kHaloH1T) != 0) {
                    virt_tkeys.push_back(vn.key);
                    virt_tptrs.push_back(vn.ds_corr.ptr());
                    h1t_virts.push_back(&vn);
                }
            }

            std::vector<Key<3>> tgt_keys;
            std::vector<T*> tgt_ptrs;
            merge_sorted_targets(real_tkeys, real_tptrs, virt_tkeys, virt_tptrs, tgt_keys, tgt_ptrs);
            if (tgt_keys.empty()) continue;

            const double loc_tol_j = f_impl->truncate_tol(eff_tol, Key<3>(static_cast<Level>(j)));
            const double sqrt_nT_raw = std::sqrt(
                static_cast<double>(std::max<std::size_t>(1, raw_T_terms[sj].size())));
            const double sqrt_nT =
                std::sqrt(static_cast<double>(std::max<std::size_t>(1, t_list.size())));
            const double tol_nn_T = 0.01 * loc_tol_j / sqrt_nT_raw;
            const double tol_sh2_T =
                dynamic_refine ? (loc_tol_j / (14.14 * sqrt_nT_raw)) : (0.1 * tol_nn_T);
            const double svd_fac_T = dynamic_refine ? 0.025 : 0.003;
            const std::size_t n_buckets_T = dynamic_refine ? kMaxNormBuckets : 1u;

            double max_T_single_op = 0.0;
            for (const ActiveBaseTerm& tr : t_list) {
                const double aw = std::abs(tr.weight_scaled);
                for (int sh = 0; sh <= tr.bw_T; ++sh) {
                    max_T_single_op += aw * tr.entry->shell_max_T[static_cast<std::size_t>(sh)];
                }
            }
            const double t_src_screen_tol = dynamic_refine ? (0.005 * loc_tol_j / 27.0) : 0.0;
            const double max_s_T_M =
                dynamic_refine ? std::max(1e-30, gstats[sj].max_s_source) : 1.0;
            const std::array<double, kMaxNormBuckets> t_brefs = {
                max_s_T_M, 0.02 * max_s_T_M, 1e-4 * max_s_T_M};

            std::array<std::vector<Key<3>>, kMaxNormBuckets> b_T_keys;
            std::array<std::vector<const T*>, kMaxNormBuckets> b_T_ptrs;
            for (const auto& se : ld.sources) {
                if (se.s_normf > 0.0 && se.s_coeff.size() == kk3 &&
                    (!dynamic_refine || se.s_normf * max_T_single_op > t_src_screen_tol)) {
                    const std::size_t bi = select_norm_bucket(se.s_normf, max_s_T_M);
                    b_T_keys[bi].push_back(se.key);
                    b_T_ptrs[bi].push_back(se.s_coeff.ptr());
                }
            }

            struct PreparedBaseTerm {
                const ActiveBaseTerm* tr = nullptr;
                int eff_bw = -1;
                std::array<SlotOp1D, StencilEntry1D::kNumDisp> slot_ops{};
            };
            std::vector<PreparedBaseTerm> prep_T;
            prep_T.reserve(t_list.size());

            for (std::size_t bi = 0; bi < n_buckets_T; ++bi) {
                if (b_T_keys[bi].empty()) continue;
                const double bref = dynamic_refine ? std::max(1e-30, t_brefs[bi]) : 1.0;
                const double term_svd_tol_T = (svd_fac_T * loc_tol_j) / (sqrt_nT * bref);

                prep_T.clear();
                for (const ActiveBaseTerm& tr : t_list) {
                    const double aw = std::abs(tr.weight_scaled);
                    const int bw_T_b = dynamic_refine
                                           ? std::min(tr.bw_T,
                                                      tr.entry->select_bw_T(
                                                          aw, bref, tol_nn_T, tol_sh2_T))
                                           : tr.bw_T;
                    if (bw_T_b < 1) continue;
                    PreparedBaseTerm pt;
                    pt.tr = &tr;
                    pt.eff_bw = select_active_1d_slots(bw_T_b,
                                                       1,
                                                       kk,
                                                       kk,
                                                       tr.entry->T,
                                                       tr.entry->svd_T,
                                                       aw,
                                                       term_svd_tol_T,
                                                       pt.slot_ops);
                    if (pt.eff_bw >= 1) {
                        prep_T.push_back(pt);
                    }
                }

                for (int m_bw = 1; m_bw <= StencilEntry1D::kMaxDisp; ++m_bw) {
                    bool any_bw = false;
                    for (const PreparedBaseTerm& pt : prep_T) {
                        if (pt.eff_bw == m_bw) {
                            any_bw = true;
                            break;
                        }
                    }
                    if (!any_bw) continue;

                    const SumFactGraph3D graph_T =
                        build_sumfact_graph_3d(b_T_keys[bi], tgt_keys, m_bw, ws.is_periodic());
                    if (graph_T.empty()) continue;

                    for (const PreparedBaseTerm& pt : prep_T) {
                        if (pt.eff_bw != m_bw) continue;
                        execute_sumfact_3d(graph_T,
                                           kk,
                                           kk,
                                           pt.slot_ops,
                                           pt.tr->weight_scaled,
                                           b_T_ptrs[bi],
                                           tgt_ptrs,
                                           sf_scratch);
                    }
                }
            }

            // Theorem 3.1 exact commutation subtraction (executed ONCE per ActiveBaseTerm,
            // strictly outside source-norm buckets!):
            //   ds_corr[V] -= W_p * M(gamma_{p, j})^{\otimes 3} s_virt[V]  for all V in dH_j^{(1,T)}.
            if (!h1t_virts.empty()) {
                for (const ActiveBaseTerm& tr : t_list) {
                    mop.fill_MT_buffer(tr.entry->gamma, mt_mat.ptr());
                    const double neg_w = -tr.weight_scaled;
                    parallel_for_chunks(
                        h1t_virts.size(),
                        kDefaultParallelGrain,
                        [&](std::size_t b, std::size_t e_hi, std::size_t cid) {
                            if (h1t_tmp_k3[cid].size() != kk3) {
                                h1t_tmp_k3[cid] = Tensor<T>(3L, k_dims, false);
                                h1t_w1_k3[cid] = Tensor<T>(3L, k_dims, false);
                            }
                            Tensor<T>& t_cube = h1t_tmp_k3[cid];
                            Tensor<T>& w_cube = h1t_w1_k3[cid];
                            for (std::size_t vi = b; vi < e_hi; ++vi) {
                                WindowedOwnedVirtNode<T>* vn_ptr = h1t_virts[vi];
                                if (!(vn_ptr->s_normf > 0.0)) continue;
                                fast_transform(vn_ptr->s_virt, mt_mat, t_cube, w_cube);
                                aligned_axpy(kk3, vn_ptr->ds_corr.ptr(), t_cube.ptr(), neg_w);
                            }
                        });
                }
            }
        }

        // =====================================================================
        // Stage 4(b): `\Gamma^j` Wavelet Decoupling & Dynamic Output Refinement
        // =====================================================================

        for (int j = 0; j < n_max_tree; ++j) {
            const std::size_t sj = static_cast<std::size_t>(j);
            const auto& g_list = active_gamma_terms[sj];
            auto& ld = ws.level_data(j);
            if (g_list.empty() || !ld.global_has_internal) continue;

            const WindowedLevelRadii& rad = radii[sj];
            const WindowedGlobalLevelStats& st = gstats[sj];
            const double loc_tol_j = f_impl->truncate_tol(eff_tol, Key<3>(static_cast<Level>(j)));
            const double sqrt_nG_raw = std::sqrt(
                static_cast<double>(std::max<std::size_t>(1, raw_gamma_terms[sj].size())));
            const double sqrt_nG =
                std::sqrt(static_cast<double>(std::max<std::size_t>(1, g_list.size())));
            const double tol_nn_G = 0.01 * loc_tol_j / sqrt_nG_raw;
            const double tol_sh2_G =
                dynamic_refine ? (loc_tol_j / (14.14 * sqrt_nG_raw)) : (0.1 * tol_nn_G);
            const double svd_fac_G = dynamic_refine ? 0.025 : 0.003;
            const std::size_t n_buckets = dynamic_refine ? kMaxNormBuckets : 1u;
            const bool can_refine_level = dynamic_refine && (j < max_refine_level);

            // -----------------------------------------------------------------
            // Pass 4b-1 (`2k -> k`, cancellation-free):
            //   `\hat{s}^j += W_p * ((P^j)^{\otimes 3} R^j_p (0, d^j))`
            // -----------------------------------------------------------------
            if (rad.has_4b1_terms && rad.m_max_4b1 >= 0) {
                double max_4b1_single_op = 0.0;
                for (const ActiveGammaTerm& gr : g_list) {
                    if (gr.bw_4b1 < 0) continue;
                    const double aw = std::abs(gr.weight_scaled);
                    for (int sh = 0; sh <= gr.bw_4b1; ++sh) {
                        max_4b1_single_op += aw * gr.entry->shell_max_4b1[static_cast<std::size_t>(sh)];
                    }
                }
                const double p1_src_screen_tol = dynamic_refine ? (0.005 * loc_tol_j / 27.0) : 0.0;

                std::vector<Key<3>> p1_real_tkeys;
                std::vector<T*> p1_real_tptrs;
                p1_real_tkeys.reserve(ld.owned_real.size());
                p1_real_tptrs.reserve(ld.owned_real.size());
                for (auto& rn : ld.owned_real) {
                    if (rn.is_internal || (rad.m_max_4b1 >= 1 && rn.near_internal_4b1)) {
                        p1_real_tkeys.push_back(rn.key);
                        p1_real_tptrs.push_back(rn.s_out.ptr());
                    }
                }

                std::vector<Key<3>> p1_virt_tkeys;
                std::vector<T*> p1_virt_tptrs;
                if (rad.m_max_4b1 >= 1) {
                    p1_virt_tkeys.reserve(ld.owned_virt.size());
                    p1_virt_tptrs.reserve(ld.owned_virt.size());
                    for (auto& vn : ld.owned_virt) {
                        if ((vn.flags & kHaloH1Gamma) != 0) {
                            p1_virt_tkeys.push_back(vn.key);
                            p1_virt_tptrs.push_back(vn.ds_corr.ptr());
                        }
                    }
                }

                std::vector<Key<3>> p1_tgt_keys;
                std::vector<T*> p1_tgt_ptrs;
                merge_sorted_targets(
                    p1_real_tkeys, p1_real_tptrs, p1_virt_tkeys, p1_virt_tptrs, p1_tgt_keys, p1_tgt_ptrs);

                if (!p1_tgt_keys.empty()) {
                    const double max_d_M = dynamic_refine ? std::max(1e-30, st.max_d_int) : 1.0;
                    const std::array<double, kMaxNormBuckets> p1_brefs = {
                        max_d_M, 0.02 * max_d_M, 1e-4 * max_d_M};

                    std::array<std::vector<Key<3>>, kMaxNormBuckets> b_src_keys;
                    std::array<std::vector<const T*>, kMaxNormBuckets> b_src_ptrs;
                    std::array<std::vector<const Tensor<T>*>, kMaxNormBuckets> b_src_tensors;

                    for (const auto& se : ld.sources) {
                        if (se.is_internal && se.d_normf > 0.0 && se.wav_coeff.size() == twok3 &&
                            (!dynamic_refine || se.d_normf * max_4b1_single_op > p1_src_screen_tol)) {
                            const std::size_t bi = select_norm_bucket(se.d_normf, max_d_M);
                            b_src_keys[bi].push_back(se.key);
                            b_src_ptrs[bi].push_back(se.wav_coeff.ptr());
                            b_src_tensors[bi].push_back(&se.wav_coeff);
                        }
                    }

                    struct PreparedP1Term {
                        const ActiveGammaTerm* gr = nullptr;
                        int eff_bw = -1;
                        std::array<SlotOp1D, StencilEntry1D::kNumDisp> slot_ops{};
                    };
                    std::vector<PreparedP1Term> prep_p1;
                    prep_p1.reserve(g_list.size());

                    for (std::size_t bi = 0; bi < n_buckets; ++bi) {
                        if (b_src_keys[bi].empty()) continue;
                        const double bref = dynamic_refine ? std::max(1e-30, p1_brefs[bi]) : 1.0;
                        const double term_svd_tol_p1 = (svd_fac_G * loc_tol_j) / (sqrt_nG * bref);

                        prep_p1.clear();
                        bool any_self_p1 = false;
                        for (const ActiveGammaTerm& gr : g_list) {
                            if (gr.bw_4b1 < 0) continue;
                            const double aw = std::abs(gr.weight_scaled);
                            const int bw_b = dynamic_refine
                                                 ? std::min(gr.bw_4b1,
                                                            gr.entry->select_bw_4b1(
                                                                aw, bref, tol_nn_G, tol_sh2_G))
                                                 : gr.bw_4b1;
                            if (bw_b < 0) continue;
                            PreparedP1Term pp;
                            pp.gr = &gr;
                            pp.eff_bw = select_active_1d_slots(bw_b,
                                                               0,
                                                               twok,
                                                               kk,
                                                               gr.entry->R_2k_to_k,
                                                               gr.entry->svd_R_2k_to_k,
                                                               aw,
                                                               term_svd_tol_p1,
                                                               pp.slot_ops);
                            if (pp.eff_bw >= 0) {
                                if (pp.eff_bw == 0) any_self_p1 = true;
                                prep_p1.push_back(pp);
                            }
                        }

                        if (any_self_p1) {
                            const auto p1_matches =
                                intersect_sorted_keys(b_src_keys[bi], p1_tgt_keys);
                            if (!p1_matches.empty()) {
                                for (const PreparedP1Term& pp : prep_p1) {
                                    if (pp.eff_bw != 0) continue;
                                    execute_self_box_2k_to_k_3d(p1_matches,
                                                                b_src_tensors[bi],
                                                                p1_tgt_ptrs,
                                                                kk,
                                                                *pp.gr->entry,
                                                                pp.gr->weight_scaled,
                                                                term_svd_tol_p1,
                                                                sf_scratch);
                                }
                            }
                        }

                        for (int m_bw = 1; m_bw <= StencilEntry1D::kMaxDisp; ++m_bw) {
                            bool any_bw = false;
                            for (const PreparedP1Term& pp : prep_p1) {
                                if (pp.eff_bw == m_bw) {
                                    any_bw = true;
                                    break;
                                }
                            }
                            if (!any_bw) continue;

                            const SumFactGraph3D graph_p1 = build_sumfact_graph_3d(
                                b_src_keys[bi], p1_tgt_keys, m_bw, ws.is_periodic());
                            if (graph_p1.empty()) continue;

                            for (const PreparedP1Term& pp : prep_p1) {
                                if (pp.eff_bw != m_bw) continue;
                                execute_sumfact_3d(graph_p1,
                                                   twok,
                                                   kk,
                                                   pp.slot_ops,
                                                   pp.gr->weight_scaled,
                                                   b_src_ptrs[bi],
                                                   p1_tgt_ptrs,
                                                   sf_scratch);
                            }
                        }
                    }
                }
            }

            // -----------------------------------------------------------------
            // Pass 4b-2 (`2k -> 2k` internal + `k -> 2k` leaf/virtual):
            //   `\hat{d}^j += W_p * Q^j (R^j_p (s^j, d^j))`
            // Includes single-pass parallel constant-subtracted pre-screening of
            // refinement-only candidates `U` (`!rn.is_internal && rn.near_internal_4b2`, `kHaloH1Refine`).
            // -----------------------------------------------------------------
            if (!rad.has_4b2_int_terms && !rad.has_4b2_leaf_terms) continue;

            constexpr int kDMax = StencilEntry1D::kMaxDisp;
            constexpr int kSpan = StencilEntry1D::kNumDisp;
            double eta_const_sh1 = 0.0;
            std::array<double, kSpan * kSpan * kSpan> K_leaf_sh1{};
            std::array<double, kSpan * kSpan * kSpan> K_leaf_sh2{};
            std::array<double, kSpan * kSpan * kSpan> K_outer_int_sh1{};
            std::array<double, kSpan * kSpan * kSpan> K_unsub_s{};
            std::array<double, kSpan * kSpan * kSpan> K_wav{};

            auto kidx = [](int dx, int dy, int dz) constexpr noexcept -> std::size_t {
                return static_cast<std::size_t>(((dx + kDMax) * kSpan + (dy + kDMax)) * kSpan +
                                                (dz + kDMax));
            };

            const int m_req_leaf = rad.has_4b2_leaf_terms ? rad.m_max_4b2_leaf : -1;
            const int m_req_int = rad.has_4b2_int_terms ? rad.m_max_4b2_int : -1;
            const int m_cand = std::max(m_req_leaf, m_req_int);
            const bool mLeaf_h2_clamped =
                dynamic_refine && (m_req_leaf >= 2) &&
                !(st.max_coarse_leaf_above > rad.sh2_s_thresh);
            const double missing_sh2_s_bound =
                mLeaf_h2_clamped ? std::min(st.max_coarse_leaf_above, rad.sh2_s_thresh) : 0.0;
            const bool use_packed_cand = (j <= windowed_halo_detail::kMaxPackedKeyLevel);
            const Translation cand_dom =
                use_packed_cand ? (Translation(1) << j) : Translation(0);

            if (can_refine_level && m_cand >= 0) {
                for (const ActiveGammaTerm& gr : g_list) {
                    const double aw = std::abs(gr.weight_scaled);
                    if (gr.bw_leaf >= 1) {
                        eta_const_sh1 += aw * gr.entry->const_wav_3d_norm[1];
                    }
                    const int b_hi = std::max(gr.bw_leaf, gr.bw_int);
                    for (int dx = -b_hi; dx <= b_hi; ++dx) {
                        for (int dy = -b_hi; dy <= b_hi; ++dy) {
                            for (int dz = -b_hi; dz <= b_hi; ++dz) {
                                const int sh = std::max({std::abs(dx), std::abs(dy), std::abs(dz)});
                                const std::size_t idx = kidx(dx, dy, dz);
                                if (sh <= 1) {
                                    if (gr.bw_leaf >= 1) {
                                        const double op_l = aw * gr.entry->opnorm_4b2_leaf(dx, dy, dz);
                                        K_leaf_sh1[idx] += op_l;
                                        K_unsub_s[idx] += op_l;
                                    } else if (sh <= gr.bw_int) {
                                        const double op_l = aw * gr.entry->opnorm_4b2_leaf(dx, dy, dz);
                                        K_outer_int_sh1[idx] += op_l;
                                        K_unsub_s[idx] += op_l;
                                    }
                                    if (sh <= gr.bw_int) {
                                        K_wav[idx] += aw * gr.entry->opnorm_4b2_wav(dx, dy, dz);
                                    }
                                } else {
                                    if (gr.bw_leaf >= 2) {
                                        const double op_l = aw * gr.entry->opnorm_4b2_leaf(dx, dy, dz);
                                        K_leaf_sh2[idx] += op_l;
                                        K_unsub_s[idx] += op_l;
                                    } else if (gr.bw_int >= 2) {
                                        K_unsub_s[idx] += aw * gr.entry->opnorm_4b2_leaf(dx, dy, dz);
                                    }
                                    if (gr.bw_int >= 2) {
                                        K_wav[idx] += aw * gr.entry->opnorm_4b2_wav(dx, dy, dz);
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Single-pass, 100% MPI-rank-invariant upper bound on `||\hat{d}^j_U||_F`:
            //   - Uses exact 27-box (`sh <= 1`) constant polynomial subtraction (`full_sh1_stencil`),
            //     never invalidated by pruned shell-2 (`sh == 2`) `H2_Ghost` sources!
            //   - Evaluates `find_windowed_source_packed` (`uint64_t`) at `j <= 20` with zero `Key<3>` allocs.
            auto candidate_exceeds_tol = [&](const Key<3>& u_key, T c0, double tol_loc) -> bool {
                if (m_cand < 0) return false;

                bool full_sh1_stencil = (m_req_leaf >= 1);
                double ub_sh1_sub = 0.0;
                double ub_sh1_unsub = 0.0;
                double ub_sh2 = 0.0;

                const auto& u_tr = u_key.translation();
                const Translation ulx = u_tr[0];
                const Translation uly = u_tr[1];
                const Translation ulz = u_tr[2];

                for (int dx = -m_cand; dx <= m_cand; ++dx) {
                    for (int dy = -m_cand; dy <= m_cand; ++dy) {
                        for (int dz = -m_cand; dz <= m_cand; ++dz) {
                            const int sh = std::max({std::abs(dx), std::abs(dy), std::abs(dz)});
                            bool in_dom = false;
                            const WindowedSourceEntry<T>* se = nullptr;
                            if (use_packed_cand) {
                                std::uint64_t s_code = 0;
                                in_dom = windowed_halo_detail::shift_packed_coords_in_domain(
                                    cand_dom, ulx, uly, ulz, -dx, -dy, -dz, ws.is_periodic(), s_code);
                                if (in_dom) {
                                    se = find_windowed_source_packed(ld, s_code);
                                }
                            } else {
                                Key<3> s_key;
                                in_dom = windowed_halo_detail::shift_key_in_domain(
                                    u_key, -dx, -dy, -dz, ws.is_periodic(), s_key);
                                if (in_dom) {
                                    se = find_windowed_source(ld.sources, s_key);
                                }
                            }

                            const std::size_t idx = kidx(dx, dy, dz);
                            if (sh <= 1) {
                                if (!in_dom || se == nullptr) {
                                    if (m_req_leaf >= 1) full_sh1_stencil = false;
                                    continue;
                                }
                                if (full_sh1_stencil) {
                                    const double diff_s = std::sqrt(
                                        se->s_osc_norm2 +
                                        windowed_halo_detail::scalar_abs2(se->s000 - c0));
                                    ub_sh1_sub += K_leaf_sh1[idx] * diff_s;
                                    if (se->is_internal && sh <= m_req_int) {
                                        ub_sh1_sub += K_outer_int_sh1[idx] * se->s_normf +
                                                      K_wav[idx] * se->d_normf;
                                    }
                                }
                                if (se->is_internal && sh <= m_req_int) {
                                    ub_sh1_unsub +=
                                        K_unsub_s[idx] * se->s_normf + K_wav[idx] * se->d_normf;
                                } else if (!se->is_internal && sh <= m_req_leaf) {
                                    ub_sh1_unsub += K_leaf_sh1[idx] * se->s_normf;
                                }
                            } else {
                                if (!in_dom) continue;
                                if (se != nullptr && se->is_internal) {
                                    if (m_req_int >= 2 || m_req_leaf >= 2) {
                                        ub_sh2 += K_unsub_s[idx] * se->s_normf +
                                                  K_wav[idx] * se->d_normf;
                                    }
                                } else if (m_req_leaf >= 2) {
                                    if (!mLeaf_h2_clamped) {
                                        if (se != nullptr) {
                                            ub_sh2 += K_leaf_sh2[idx] * se->s_normf;
                                        }
                                    } else {
                                        if (se != nullptr && !se->is_virtual) {
                                            ub_sh2 += K_leaf_sh2[idx] * se->s_normf;
                                        } else {
                                            ub_sh2 += K_leaf_sh2[idx] * missing_sh2_s_bound;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }

                const double ub_sh1 =
                    full_sh1_stencil
                        ? (std::sqrt(windowed_halo_detail::scalar_abs2(c0)) * eta_const_sh1 +
                           ub_sh1_sub)
                        : ub_sh1_unsub;
                const double ub = ub_sh1 + ub_sh2;
                return (ub * (1.0 + 1e-12) > tol_loc);
            };

            std::vector<Key<3>> p2_real_tkeys;
            std::vector<T*> p2_real_tptrs;
            if (!ld.owned_real.empty()) {
                std::vector<std::uint8_t> keep_real(ld.owned_real.size(), 0u);
                parallel_for_chunks(
                    ld.owned_real.size(),
                    kDefaultParallelGrain,
                    [&](std::size_t b, std::size_t e_hi, std::size_t) {
                        for (std::size_t ri = b; ri < e_hi; ++ri) {
                            auto& rn = ld.owned_real[ri];
                            bool keep = rn.is_internal;
                            if (!keep && can_refine_level && rn.near_internal_4b2) {
                                const double tol_loc =
                                    f_impl->truncate_tol(eff_thresh, rn.key) * eff_refine_fac;
                                keep = candidate_exceeds_tol(rn.key, rn.s000, tol_loc);
                            }
                            if (keep) {
                                keep_real[ri] = 1u;
                                if (rn.ns_out.size() != twok3) {
                                    rn.ns_out = Tensor<T>(twok, twok, twok);
                                }
                            }
                        }
                    });
                p2_real_tkeys.reserve(ld.owned_real.size());
                p2_real_tptrs.reserve(ld.owned_real.size());
                for (std::size_t ri = 0; ri < ld.owned_real.size(); ++ri) {
                    if (keep_real[ri] != 0u) {
                        p2_real_tkeys.push_back(ld.owned_real[ri].key);
                        p2_real_tptrs.push_back(ld.owned_real[ri].ns_out.ptr());
                    }
                }
            }

            std::vector<Key<3>> p2_virt_tkeys;
            std::vector<T*> p2_virt_tptrs;
            if (can_refine_level && !ld.owned_virt.empty()) {
                std::vector<std::uint8_t> keep_virt(ld.owned_virt.size(), 0u);
                parallel_for_chunks(
                    ld.owned_virt.size(),
                    kDefaultParallelGrain,
                    [&](std::size_t b, std::size_t e_hi, std::size_t) {
                        for (std::size_t vi = b; vi < e_hi; ++vi) {
                            auto& vn = ld.owned_virt[vi];
                            if ((vn.flags & kHaloH1Refine) == 0) continue;
                            const double tol_loc =
                                f_impl->truncate_tol(eff_thresh, vn.key) * eff_refine_fac;
                            if (!candidate_exceeds_tol(vn.key, vn.s000, tol_loc)) continue;
                            keep_virt[vi] = 1u;
                            if (vn.ns_virt_out.size() != twok3) {
                                vn.ns_virt_out = Tensor<T>(twok, twok, twok);
                            }
                        }
                    });
                p2_virt_tkeys.reserve(ld.owned_virt.size());
                p2_virt_tptrs.reserve(ld.owned_virt.size());
                for (std::size_t vi = 0; vi < ld.owned_virt.size(); ++vi) {
                    if (keep_virt[vi] != 0u) {
                        p2_virt_tkeys.push_back(ld.owned_virt[vi].key);
                        p2_virt_tptrs.push_back(ld.owned_virt[vi].ns_virt_out.ptr());
                    }
                }
            }

            std::vector<Key<3>> p2_tgt_keys;
            std::vector<T*> p2_tgt_ptrs;
            merge_sorted_targets(
                p2_real_tkeys, p2_real_tptrs, p2_virt_tkeys, p2_virt_tptrs, p2_tgt_keys, p2_tgt_ptrs);

            if (!p2_tgt_keys.empty()) {
                double max_4b2_int_op = 0.0;
                double max_4b2_leaf_op = 0.0;
                for (const ActiveGammaTerm& gr : g_list) {
                    const double aw = std::abs(gr.weight_scaled);
                    for (int sh = 0; sh <= gr.bw_int; ++sh) {
                        max_4b2_int_op +=
                            aw * gr.entry->shell_max_4b2_int[static_cast<std::size_t>(sh)];
                    }
                    for (int sh = 0; sh <= gr.bw_leaf; ++sh) {
                        max_4b2_leaf_op +=
                            aw * gr.entry->shell_max_4b2_leaf[static_cast<std::size_t>(sh)];
                    }
                }
                const double p2_src_screen_tol = dynamic_refine ? (0.005 * loc_tol_j / 27.0) : 0.0;

                const double max_ns_M = dynamic_refine ? std::max(1e-30, st.max_ns_int) : 1.0;
                const double max_s_M = dynamic_refine ? std::max(1e-30, st.max_s_source) : 1.0;
                const std::array<double, kMaxNormBuckets> int_brefs = {
                    max_ns_M, 0.02 * max_ns_M, 1e-4 * max_ns_M};
                const std::array<double, kMaxNormBuckets> leaf_brefs = {
                    max_s_M, 0.02 * max_s_M, 1e-4 * max_s_M};

                std::array<std::vector<Key<3>>, kMaxNormBuckets> b_int_keys;
                std::array<std::vector<const T*>, kMaxNormBuckets> b_int_ptrs;
                std::array<std::vector<const Tensor<T>*>, kMaxNormBuckets> b_int_tensors;
                std::array<std::vector<Key<3>>, kMaxNormBuckets> b_leaf_keys;
                std::array<std::vector<const T*>, kMaxNormBuckets> b_leaf_ptrs;

                for (const auto& se : ld.sources) {
                    if (se.is_internal) {
                        if (se.ns_normf > 0.0 && se.ns_coeff.size() == twok3 &&
                            (!dynamic_refine || se.ns_normf * max_4b2_int_op > p2_src_screen_tol)) {
                            const std::size_t bi = select_norm_bucket(se.ns_normf, max_ns_M);
                            b_int_keys[bi].push_back(se.key);
                            b_int_ptrs[bi].push_back(se.ns_coeff.ptr());
                            b_int_tensors[bi].push_back(&se.ns_coeff);
                        }
                    } else if (se.s_normf > 0.0 && se.s_coeff.size() == kk3 &&
                               (!dynamic_refine ||
                                se.s_normf * max_4b2_leaf_op > p2_src_screen_tol)) {
                        const std::size_t bi = select_norm_bucket(se.s_normf, max_s_M);
                        b_leaf_keys[bi].push_back(se.key);
                        b_leaf_ptrs[bi].push_back(se.s_coeff.ptr());
                    }
                }

                // 1. Internal `(s, d) -> d_out` contractions (`bw_int in {0, 1, 2}`)
                if (rad.has_4b2_int_terms) {
                    struct PreparedP2IntTerm {
                        const ActiveGammaTerm* gr = nullptr;
                        int eff_bw = -1;
                        std::array<SlotOp1D, StencilEntry1D::kNumDisp> slot_ops{};
                    };
                    std::vector<PreparedP2IntTerm> prep_int;
                    prep_int.reserve(g_list.size());

                    for (std::size_t bi = 0; bi < n_buckets; ++bi) {
                        if (b_int_keys[bi].empty()) continue;
                        const double bref = dynamic_refine ? std::max(1e-30, int_brefs[bi]) : 1.0;
                        const double term_svd_tol_int = (svd_fac_G * loc_tol_j) / (sqrt_nG * bref);

                        prep_int.clear();
                        bool any_self_int = false;
                        for (const ActiveGammaTerm& gr : g_list) {
                            if (gr.bw_int < 0) continue;
                            const double aw = std::abs(gr.weight_scaled);
                            const int bw_leaf_b =
                                dynamic_refine
                                    ? std::min(gr.bw_leaf,
                                               gr.entry->select_bw_4b2_leaf(
                                                   aw, bref, tol_nn_G, tol_sh2_G))
                                    : gr.bw_leaf;
                            const int bw_int_b =
                                dynamic_refine
                                    ? std::min(gr.bw_int,
                                               gr.entry->select_bw_4b2_int(
                                                   aw, bref, bref, tol_nn_G, tol_sh2_G, bw_leaf_b))
                                    : gr.bw_int;
                            if (bw_int_b < 0) continue;

                            PreparedP2IntTerm pi;
                            pi.gr = &gr;
                            const int min_bw_int = (bw_leaf_b >= 1) ? 1 : 0;
                            pi.eff_bw = select_active_1d_slots(bw_int_b,
                                                               min_bw_int,
                                                               twok,
                                                               twok,
                                                               gr.entry->R,
                                                               gr.entry->svd_R,
                                                               aw,
                                                               term_svd_tol_int,
                                                               pi.slot_ops);
                            if (pi.eff_bw >= 0) {
                                if (pi.eff_bw == 0) any_self_int = true;
                                prep_int.push_back(pi);
                            }
                        }

                        if (any_self_int) {
                            const auto p2_int_matches =
                                intersect_sorted_keys(b_int_keys[bi], p2_tgt_keys);
                            if (!p2_int_matches.empty()) {
                                for (const PreparedP2IntTerm& pi : prep_int) {
                                    if (pi.eff_bw != 0) continue;
                                    execute_self_box_2k_to_2k_3d(p2_int_matches,
                                                                 b_int_tensors[bi],
                                                                 p2_tgt_ptrs,
                                                                 twok,
                                                                 *pi.gr->entry,
                                                                 pi.gr->weight_scaled,
                                                                 term_svd_tol_int,
                                                                 sf_scratch);
                                }
                            }
                        }

                        for (int m_bw = 1; m_bw <= StencilEntry1D::kMaxDisp; ++m_bw) {
                            bool any_bw = false;
                            for (const PreparedP2IntTerm& pi : prep_int) {
                                if (pi.eff_bw == m_bw) {
                                    any_bw = true;
                                    break;
                                }
                            }
                            if (!any_bw) continue;

                            const SumFactGraph3D graph_p2_int = build_sumfact_graph_3d(
                                b_int_keys[bi], p2_tgt_keys, m_bw, ws.is_periodic());
                            if (graph_p2_int.empty()) continue;

                            for (const PreparedP2IntTerm& pi : prep_int) {
                                if (pi.eff_bw != m_bw) continue;
                                execute_sumfact_3d(graph_p2_int,
                                                   twok,
                                                   twok,
                                                   pi.slot_ops,
                                                   pi.gr->weight_scaled,
                                                   b_int_ptrs[bi],
                                                   p2_tgt_ptrs,
                                                   sf_scratch);
                            }
                        }
                    }
                }

                // 2. Leaf & virtual `(s, 0) -> d_out` contractions (`bw_leaf in {1, 2}`, never `0`)
                if (rad.has_4b2_leaf_terms) {
                    struct PreparedP2LeafTerm {
                        const ActiveGammaTerm* gr = nullptr;
                        int eff_bw = -1;
                        std::array<SlotOp1D, StencilEntry1D::kNumDisp> slot_ops{};
                    };
                    std::vector<PreparedP2LeafTerm> prep_leaf;
                    prep_leaf.reserve(g_list.size());

                    for (std::size_t bi = 0; bi < n_buckets; ++bi) {
                        if (b_leaf_keys[bi].empty()) continue;
                        const double bref = dynamic_refine ? std::max(1e-30, leaf_brefs[bi]) : 1.0;
                        const double term_svd_tol_leaf = (svd_fac_G * loc_tol_j) / (sqrt_nG * bref);

                        prep_leaf.clear();
                        for (const ActiveGammaTerm& gr : g_list) {
                            if (gr.bw_leaf < 1) continue;
                            const double aw = std::abs(gr.weight_scaled);
                            const int bw_leaf_b =
                                dynamic_refine
                                    ? std::min(gr.bw_leaf,
                                               gr.entry->select_bw_4b2_leaf(
                                                   aw, bref, tol_nn_G, tol_sh2_G))
                                    : gr.bw_leaf;
                            if (bw_leaf_b < 1) continue;

                            PreparedP2LeafTerm pl;
                            pl.gr = &gr;
                            pl.eff_bw = select_active_1d_slots(bw_leaf_b,
                                                               1,
                                                               kk,
                                                               twok,
                                                               gr.entry->R_k_to_2k,
                                                               gr.entry->svd_R_k_to_2k,
                                                               aw,
                                                               term_svd_tol_leaf,
                                                               pl.slot_ops);
                            if (pl.eff_bw >= 1) {
                                prep_leaf.push_back(pl);
                            }
                        }

                        for (int m_bw = 1; m_bw <= StencilEntry1D::kMaxDisp; ++m_bw) {
                            bool any_bw = false;
                            for (const PreparedP2LeafTerm& pl : prep_leaf) {
                                if (pl.eff_bw == m_bw) {
                                    any_bw = true;
                                    break;
                                }
                            }
                            if (!any_bw) continue;

                            const SumFactGraph3D graph_p2_leaf = build_sumfact_graph_3d(
                                b_leaf_keys[bi], p2_tgt_keys, m_bw, ws.is_periodic());
                            if (graph_p2_leaf.empty()) continue;

                            for (const PreparedP2LeafTerm& pl : prep_leaf) {
                                if (pl.eff_bw != m_bw) continue;
                                execute_sumfact_3d(graph_p2_leaf,
                                                   kk,
                                                   twok,
                                                   pl.slot_ops,
                                                   pl.gr->weight_scaled,
                                                   b_leaf_ptrs[bi],
                                                   p2_tgt_ptrs,
                                                   sf_scratch);
                            }
                        }
                    }
                }

                // Clear `s0` octant on allocated `ns_out` / `ns_virt_out` and evaluate dynamic refinement.
                parallel_for_chunks(
                    ld.owned_real.size(),
                    kDefaultParallelGrain,
                    [&](std::size_t b, std::size_t e_hi, std::size_t) {
                        for (std::size_t ri = b; ri < e_hi; ++ri) {
                            auto& rn = ld.owned_real[ri];
                            if (rn.ns_out.size() != twok3) continue;
                            zero_s0_octant(rn.ns_out, k);
                            if (!rn.is_internal && can_refine_level && rn.near_internal_4b2) {
                                const double tol_loc =
                                    f_impl->truncate_tol(eff_thresh, rn.key) * eff_refine_fac;
                                if (windowed_halo_detail::wavelet_octants_normf(rn.ns_out, k) > tol_loc) {
                                    rn.promoted_leaf = true;
                                }
                            }
                        }
                    });
                if (can_refine_level && !ld.owned_virt.empty()) {
                    parallel_for_chunks(
                        ld.owned_virt.size(),
                        kDefaultParallelGrain,
                        [&](std::size_t b, std::size_t e_hi, std::size_t) {
                            for (std::size_t vi = b; vi < e_hi; ++vi) {
                                auto& vn = ld.owned_virt[vi];
                                if (vn.ns_virt_out.size() != twok3) continue;
                                zero_s0_octant(vn.ns_virt_out, k);
                                const double tol_loc =
                                    f_impl->truncate_tol(eff_thresh, vn.key) * eff_refine_fac;
                                if (windowed_halo_detail::wavelet_octants_normf(vn.ns_virt_out, k) >
                                    tol_loc) {
                                    vn.promoted = true;
                                }
                            }
                        });
                }
            }
        }

        // =====================================================================
        // Stage 5a: Bottom-Up Virtual Restriction & Nonstandard Tree Assembly
        // =====================================================================
        ws.assemble_output_nonstandard_tree(result);
    }
}

/// Apply a 3D isotropic Gaussian-quadrature integral operator (`WindowedQuadraturePlan`)
/// to `f` using Version 7 Tail-Completed Windowed Non-Standard (NS) form.
///
/// Supports both `T = double` and `T = double_complex`, arbitrary unbalanced octrees
/// up to `plan.n_max()` without 2:1 grading, single-rank and multi-rank MPI execution,
/// subworld `World&` contexts, and optional `destructive_input` input clearing.
template <typename T>
Function<T, 3> apply_windowed_ns_3d(
    World& world,
    const Function<T, 3>& f,
    const WindowedQuadraturePlan& plan,
    const WindowedStencilCache1D& stencil_cache,
    double thresh,
    bool dynamic_refine = true,
    double dynamic_refine_factor = 1.0,
    bool fence = true,
    bool destructive_input = false) {
    MADNESS_CHECK_THROW(f.is_initialized(), "apply_windowed_ns_3d: input Function must be initialized");
    MADNESS_CHECK_THROW(!f.is_on_demand(), "apply_windowed_ns_3d: on-demand input Function is not supported");
    MADNESS_CHECK_THROW(
        f.k() == plan.k() && f.k() == stencil_cache.k(),
        "apply_windowed_ns_3d: mismatch between Function::k(), plan.k(), and stencil_cache.k()");

    const int k = f.k();
    const double eff_thresh = (thresh > 0.0) ? thresh : f.thresh();
    const double eff_refine_fac = (dynamic_refine_factor > 0.0) ? dynamic_refine_factor : 1.0;

    if (&f.world() == &world) {
        Function<T, 3>& ff = const_cast<Function<T, 3>&>(f);
        const bool was_compressed = prepare_windowed_ns_input_3d(ff, true);

        Function<T, 3> result =
            FunctionFactory<T, 3>(world).k(k).thresh(eff_thresh).pmap(ff.get_pmap()).empty();
        result.get_impl()->get_coeffs().process_pending();
        if (world.size() > 1) {
            world.gop.fence();
        }

        apply_windowed_ns_3d_prepared(
            world, ff, result, plan, stencil_cache, eff_thresh, dynamic_refine, eff_refine_fac);

        restore_windowed_ns_input_3d(ff, was_compressed, destructive_input, false);
        result.reconstruct(fence);
        return result;
    }

    Function<T, 3> result;
    {
        Function<T, 3> f_sub = copy(world, f, true);
        const bool was_compressed = prepare_windowed_ns_input_3d(f_sub, true);

        result =
            FunctionFactory<T, 3>(world).k(k).thresh(eff_thresh).pmap(f_sub.get_pmap()).empty();
        result.get_impl()->get_coeffs().process_pending();
        if (world.size() > 1) {
            world.gop.fence();
        }

        apply_windowed_ns_3d_prepared(
            world, f_sub, result, plan, stencil_cache, eff_thresh, dynamic_refine, eff_refine_fac);

        restore_windowed_ns_input_3d(f_sub, was_compressed, true, true);
        result.reconstruct(true);
    }
    if (destructive_input) {
        const_cast<Function<T, 3>&>(f).clear(fence);
    } else if (fence) {
        world.gop.fence();
    }
    return result;
}

} // namespace madness

#endif // MADNESS_MRA_WINDOWED_NS_SUMFACT_H__INCLUDED
