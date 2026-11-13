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

#ifndef MADNESS_MRA_WINDOWED_NS_STENCIL_H__INCLUDED
#define MADNESS_MRA_WINDOWED_NS_STENCIL_H__INCLUDED

/// \file windowed_ns_stencil.h
/// \brief Deep-tree 64-bit SplitMix64 key hashing, L1/L2-cached universal
///        scale-locked 1D two-scale transition tables, exact 3D sub-block
///        Frobenius norms, and 1D SVD decompositions for Windowed NS 3D convolutions.

#include <madness/mra/windowed_quadrature.h>
#include <madness/mra/convolution1d.h>
#include <madness/mra/key.h>
#include <madness/tensor/mTxmq.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/madness_exception.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <mutex>
#include <vector>

namespace madness {

#ifndef MADNESS_WINDOWED_KEY_HASH_DEFINED
#define MADNESS_WINDOWED_KEY_HASH_DEFINED

/// Overflow-free 64-bit SplitMix64 hash over `Key<3>` `(level, l[0], l[1], l[2])`.
///
/// Supports arbitrary octree depth (`Level` up to `MAXLEVEL = 62`, full 64-bit
/// `Translation` per axis) without packing `(lx, ly, lz)` into a 64-bit Morton integer.
struct WindowedKeyHash {
    static constexpr std::uint64_t splitmix64(std::uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31);
    }

    std::size_t operator()(const Key<3>& key) const noexcept {
        const auto& l = key.translation();
        const std::uint64_t n = static_cast<std::uint64_t>(key.level());
        const std::uint64_t lx = static_cast<std::uint64_t>(l[0]);
        const std::uint64_t ly = static_cast<std::uint64_t>(l[1]);
        const std::uint64_t lz = static_cast<std::uint64_t>(l[2]);

        std::uint64_t h = splitmix64(n ^ 0x243f6a8885a308d3ULL);
        h = splitmix64(h ^ splitmix64(lx + 0x13198a2e03707344ULL));
        h = splitmix64(h ^ splitmix64(ly + 0xa4093822299f31d0ULL));
        h = splitmix64(h ^ splitmix64(lz + 0x082efa98ec4e6c89ULL));
        return static_cast<std::size_t>(h);
    }
};

#endif // MADNESS_WINDOWED_KEY_HASH_DEFINED

/// Precomputed 1D SVD decomposition `A = U * SVT` (`SVT` rows scaled by singular values)
/// with dense-packed prefix matrices `U_by_rank[r]` (`(d_in, r)`, `ldb == r` for Tier-1
/// `libxsmm` JIT `s_dense_table`) and exact Frobenius tail norms `tail_normf[r] = ||A - U_r SVT_r||_F`.
struct StencilSVD1D {
    double normf = 0.0;
    double norm_s1 = 0.0;
    Tensor<double> U;                        ///< `(d_in, rmax)` left singular vectors.
    Tensor<double> SVT;                      ///< `(rmax, d_out)` singular-value-scaled right singular vectors.
    std::vector<Tensor<double>> U_by_rank;   ///< Length `rmax + 1`: `U_by_rank[r]` is contiguous `(d_in, r)` (`ldb == r`).
    std::vector<double> tail_normf;          ///< Length `rmax + 1`: `tail_normf[r]` for rank `r in [0, rmax]`.

    int required_rank(double abs_frob_tol) const noexcept {
        if (normf <= abs_frob_tol || tail_normf.empty()) return 0;
        const int rmax = static_cast<int>(tail_normf.size()) - 1;
        for (int r = 1; r <= rmax; ++r) {
            if (tail_normf[static_cast<std::size_t>(r)] <= abs_frob_tol) {
                return r;
            }
        }
        return rmax;
    }
};

/// Pre-transposed 1D two-scale and scaling transition matrices, SVD factors, and exact 3D
/// NS sub-block Frobenius norms for a single dimensionless stiffness `gamma = alpha_sim * 4^{-j} > 0`
/// on the reference `[0, 1]` grid.
///
/// Displacement convention matches `Convolution1D::nonstandard(0, disp)`:
///   `disp = l_target - l_source` (`l_source = l_target - disp`).
///
/// All stored matrices are contiguous row-major `(dim_in, dim_out)` oriented directly
/// for `mTxmq(dimi, dim_out, dim_in, c, a, mat.ptr())`:
///   - `R[slot]`:          `(2k, 2k)` full two-scale transition (`[s, d] -> [s, d]`).
///   - `R_2k_to_k[slot]`:  `(2k,  k)` scaling-output columns `0..k-1` of `R` (`[s, d] -> s`).
///   - `R_k_to_2k[slot]`:  `( k, 2k)` scaling-input rows `0..k-1` of `R` (`s -> [s, d]`).
///   - `T[slot]`:          `( k,  k)` scaling-to-scaling transition (`s -> s`).
struct StencilEntry1D {
    static constexpr int kMaxDisp = 2;
    static constexpr int kNumDisp = 2 * kMaxDisp + 1;

    double gamma = 0.0;
    int k = 0;
    int m_max_T = 1;   ///< Static upper half-bandwidth for base scaling T^{j_0}.
    int m_max_ns = 2;  ///< Static upper half-bandwidth for NS wavelet decoupling Gamma^j.

    std::array<Tensor<double>, kNumDisp> R;
    std::array<Tensor<double>, kNumDisp> R_2k_to_k;
    std::array<Tensor<double>, kNumDisp> R_k_to_2k;
    std::array<Tensor<double>, kNumDisp> T;

    // Direct ConvolutionData1D SVD approximations (`operator.h` strided-U / VT layout)
    std::array<Tensor<double>, kNumDisp> RU;
    std::array<Tensor<double>, kNumDisp> RVT;
    std::array<Tensor<double>, kNumDisp> Rs;
    std::array<double, kNumDisp> Rnorm{};

    std::array<Tensor<double>, kNumDisp> TU;
    std::array<Tensor<double>, kNumDisp> TVT;
    std::array<Tensor<double>, kNumDisp> Ts;
    std::array<double, kNumDisp> Tnorm{};

    // Rectangular / square 1D SVD decompositions with exact Frobenius tail norms
    std::array<StencilSVD1D, kNumDisp> svd_R;
    std::array<StencilSVD1D, kNumDisp> svd_R_2k_to_k;
    std::array<StencilSVD1D, kNumDisp> svd_R_k_to_2k;
    std::array<StencilSVD1D, kNumDisp> svd_T;

    // Exact 1D sub-block Frobenius norms per displacement slot `disp + kMaxDisp`
    std::array<double, kNumDisp> normf_T{};           ///< `||T||_F`
    std::array<double, kNumDisp> normf_R_ss{};        ///< `||R[0..k, 0..k)||_F`
    std::array<double, kNumDisp> normf_R_to_s{};      ///< `||R[0..2k, 0..k)||_F` (`R_2k_to_k`)
    std::array<double, kNumDisp> normf_R_s_to_all{};  ///< `||R[0..k, 0..2k)||_F` (`R_k_to_2k`)
    std::array<double, kNumDisp> normf_R_full{};      ///< `||R[0..2k, 0..2k)||_F` (`R`)
    std::array<double, kNumDisp> normf_R_d_to_all{};  ///< `||R[k..2k, 0..2k)||_F`
    std::array<double, kNumDisp> normf_R_ds{};        ///< `||R[k..2k, 0..k)||_F`

    // Precomputed 3D shell norms (`sh = max({|dx|,|dy|,|dz|}) in {0, 1, 2}`):
    //   `shell_max_*[sh]` = max single-box 3D Frobenius norm in shell `sh`
    //   `shell_l2_*[sh]`  = sqrt(sum_{||d||_inf == sh} opnorm(d)^2) over entire shell `sh`
    std::array<double, kMaxDisp + 1> shell_max_T{};
    std::array<double, kMaxDisp + 1> shell_max_4b1{};
    std::array<double, kMaxDisp + 1> shell_max_4b2_leaf{};
    std::array<double, kMaxDisp + 1> shell_max_4b2_int{};
    std::array<double, kMaxDisp + 1> shell_max_4b2_wav{};

    std::array<double, kMaxDisp + 1> shell_l2_T{};
    std::array<double, kMaxDisp + 1> shell_l2_4b1{};
    std::array<double, kMaxDisp + 1> shell_l2_4b2_leaf{};
    std::array<double, kMaxDisp + 1> shell_l2_4b2_int{};
    std::array<double, kMaxDisp + 1> shell_l2_4b2_wav{};

    /// Exact 3D wavelet-octant Frobenius norm `||Q^j (\sum_{|d| <= m_bw} R(d)^{\otimes 3} e_{000})||_F`
    /// of a unit constant scaling field `e_{000}` summed over `[-m_bw, +m_bw]^3` (`m_bw in {0, 1, 2}`).
    std::array<double, kMaxDisp + 1> const_wav_3d_norm{};

    const Tensor<double>& mat_R(int disp) const noexcept {
        return R[static_cast<std::size_t>(disp + kMaxDisp)];
    }

    const Tensor<double>& mat_R_2k_to_k(int disp) const noexcept {
        return R_2k_to_k[static_cast<std::size_t>(disp + kMaxDisp)];
    }

    const Tensor<double>& mat_R_k_to_2k(int disp) const noexcept {
        return R_k_to_2k[static_cast<std::size_t>(disp + kMaxDisp)];
    }

    const Tensor<double>& mat_T(int disp) const noexcept {
        return T[static_cast<std::size_t>(disp + kMaxDisp)];
    }

    /// Exact 3D Frobenius norm of `T(dx) \otimes T(dy) \otimes T(dz)` (`s -> s`).
    double opnorm_T(int dx, int dy, int dz) const noexcept {
        const std::size_t sx = static_cast<std::size_t>(dx + kMaxDisp);
        const std::size_t sy = static_cast<std::size_t>(dy + kMaxDisp);
        const std::size_t sz = static_cast<std::size_t>(dz + kMaxDisp);
        return normf_T[sx] * normf_T[sy] * normf_T[sz];
    }

    /// Exact 3D Frobenius norm of Pass 4b-1 block `(0, d_in) -> s_out` (7 wavelet octants to scaling).
    double opnorm_4b1(int dx, int dy, int dz) const noexcept {
        const std::size_t sx = static_cast<std::size_t>(dx + kMaxDisp);
        const std::size_t sy = static_cast<std::size_t>(dy + kMaxDisp);
        const std::size_t sz = static_cast<std::size_t>(dz + kMaxDisp);
        const double p_to_s = normf_R_to_s[sx] * normf_R_to_s[sy] * normf_R_to_s[sz];
        const double p_ss = normf_R_ss[sx] * normf_R_ss[sy] * normf_R_ss[sz];
        return std::sqrt(std::max(0.0, p_to_s * p_to_s - p_ss * p_ss));
    }

    /// Exact 3D Frobenius norm of Pass 4b-2 leaf block `(s_in, 0) -> d_out` (scaling to 7 wavelet octants).
    double opnorm_4b2_leaf(int dx, int dy, int dz) const noexcept {
        const std::size_t sx = static_cast<std::size_t>(dx + kMaxDisp);
        const std::size_t sy = static_cast<std::size_t>(dy + kMaxDisp);
        const std::size_t sz = static_cast<std::size_t>(dz + kMaxDisp);
        const double p_s_all = normf_R_s_to_all[sx] * normf_R_s_to_all[sy] * normf_R_s_to_all[sz];
        const double p_ss = normf_R_ss[sx] * normf_R_ss[sy] * normf_R_ss[sz];
        return std::sqrt(std::max(0.0, p_s_all * p_s_all - p_ss * p_ss));
    }

    /// Exact 3D Frobenius norm of Pass 4b-2 internal block `(s_in, d_in) -> d_out`.
    double opnorm_4b2_int(int dx, int dy, int dz) const noexcept {
        const std::size_t sx = static_cast<std::size_t>(dx + kMaxDisp);
        const std::size_t sy = static_cast<std::size_t>(dy + kMaxDisp);
        const std::size_t sz = static_cast<std::size_t>(dz + kMaxDisp);
        const double p_full = normf_R_full[sx] * normf_R_full[sy] * normf_R_full[sz];
        const double p_to_s = normf_R_to_s[sx] * normf_R_to_s[sy] * normf_R_to_s[sz];
        return std::sqrt(std::max(0.0, p_full * p_full - p_to_s * p_to_s));
    }

    /// Exact 3D Frobenius norm of pure wavelet-to-wavelet block `(0, d_in) -> d_out`
    /// across ALL 7 input wavelet octants and ALL 7 output wavelet octants.
    double opnorm_4b2_wav(int dx, int dy, int dz) const noexcept {
        const std::size_t sx = static_cast<std::size_t>(dx + kMaxDisp);
        const std::size_t sy = static_cast<std::size_t>(dy + kMaxDisp);
        const std::size_t sz = static_cast<std::size_t>(dz + kMaxDisp);
        const double p_full = normf_R_full[sx] * normf_R_full[sy] * normf_R_full[sz];
        const double p_to_s = normf_R_to_s[sx] * normf_R_to_s[sy] * normf_R_to_s[sz];
        const double p_s_all = normf_R_s_to_all[sx] * normf_R_s_to_all[sy] * normf_R_s_to_all[sz];
        const double p_ss = normf_R_ss[sx] * normf_R_ss[sy] * normf_R_ss[sz];
        const double diff2 = (p_full * p_full - p_to_s * p_to_s) - (p_s_all * p_s_all - p_ss * p_ss);
        return std::sqrt(std::max(0.0, diff2));
    }

    /// Dynamic half-bandwidth in `{-1, 1, 2}` for Stage 4(a) base scaling `T`
    /// (`>= 1` whenever active so Theorem 3.1 face/edge/corner halo commutation stays intact).
    int select_bw_T(double abs_weight,
                    double ref_s_source,
                    double tol_nn,
                    double tol_sh2) const noexcept {
        const double scale = abs_weight * ref_s_source;
        if (!(scale > 0.0)) return -1;
        if (m_max_T >= 2 && scale * shell_max_T[2] > tol_sh2) return 2;
        if (scale * std::max(shell_l2_T[1], shell_l2_T[0]) > tol_nn) return 1;
        return -1;
    }

    /// Dynamic half-bandwidth in `{-1, 0, 1, 2}` for Pass 4b-1 `(0, d_int) -> s_out`.
    int select_bw_4b1(double abs_weight,
                      double ref_d_int,
                      double tol_nn,
                      double tol_sh2) const noexcept {
        const double scale = abs_weight * ref_d_int;
        if (!(scale > 0.0)) return -1;
        if (m_max_ns >= 2 && scale * shell_max_4b1[2] > tol_sh2) return 2;
        if (scale * shell_l2_4b1[1] > tol_nn) return 1;
        if (scale * shell_l2_4b1[0] > tol_nn) return 0;
        return -1;
    }

    /// Dynamic half-bandwidth in `{-1, 1, 2}` for Pass 4b-2 leaf `(s, 0) -> d_out`.
    /// NEVER returns `0`: preserves Theorem 3.1 `[-1, 1]^3` neighbor polynomial cancellation!
    int select_bw_4b2_leaf(double abs_weight,
                           double ref_s_source,
                           double tol_nn,
                           double tol_sh2) const noexcept {
        const double scale = abs_weight * ref_s_source;
        if (!(scale > 0.0)) return -1;
        if (m_max_ns >= 2 && scale * shell_max_4b2_leaf[2] > tol_sh2) return 2;
        if (scale * std::max(shell_l2_4b2_leaf[1], shell_l2_4b2_leaf[0]) > tol_nn) return 1;
        return -1;
    }

    /// Dynamic half-bandwidth in `{-1, 0, 1, 2}` for Pass 4b-2 internal `(s, d) -> d_out`,
    /// guaranteed `>= bw_leaf` so mixed internal+leaf neighborhoods share uniform cancellation radii.
    int select_bw_4b2_int(double abs_weight,
                          double ref_s_real,
                          double ref_d_int,
                          double tol_nn,
                          double tol_sh2,
                          int bw_leaf) const noexcept {
        int bw = bw_leaf;
        if (abs_weight > 0.0) {
            const double sh2_max =
                abs_weight * (shell_max_4b2_leaf[2] * ref_s_real + shell_max_4b2_wav[2] * ref_d_int);
            const double sh1_l2 =
                abs_weight * (shell_l2_4b2_leaf[1] * ref_s_real + shell_l2_4b2_wav[1] * ref_d_int);
            const double sh0_l2 =
                abs_weight * (shell_l2_4b2_leaf[0] * ref_s_real + shell_l2_4b2_wav[0] * ref_d_int);
            if (m_max_ns >= 2 && sh2_max > tol_sh2) {
                bw = std::max(bw, 2);
            } else if (sh1_l2 > tol_nn) {
                bw = std::max(bw, 1);
            } else if (sh0_l2 > tol_nn) {
                bw = std::max(bw, 0);
            }
        }
        return bw;
    }
};

namespace windowed_stencil_detail {

inline StencilSVD1D build_svd_1d(const Tensor<double>& mat) {
    StencilSVD1D out;
    out.normf = mat.normf();
    if (!(out.normf > 1e-20) || mat.ndim() != 2 || mat.dim(0) <= 0 || mat.dim(1) <= 0) {
        out.normf = 0.0;
        return out;
    }

    Tensor<double> s;
    Tensor<double> VT;
    svd(mat, out.U, s, VT);

    const long d_in = out.U.dim(0);
    const long rmax = s.dim(0);
    const long d_out = VT.dim(1);
    out.SVT = copy(VT);
    for (long r = 0; r < rmax; ++r) {
        const double sr = s(r);
        out.norm_s1 += sr;
        double* row = out.SVT.ptr() + r * d_out;
        madness::cblas::scal(static_cast<integer>(d_out), sr, row, 1);
    }

    out.U_by_rank.resize(static_cast<std::size_t>(rmax + 1));
    const Slice s_rows(0, d_in - 1);
    for (long r = 1; r <= rmax; ++r) {
        out.U_by_rank[static_cast<std::size_t>(r)] = copy(out.U(s_rows, Slice(0, r - 1)));
    }

    out.tail_normf.assign(static_cast<std::size_t>(rmax + 1), 0.0);
    double rem_sq = 0.0;
    for (long r = rmax - 1; r >= 0; --r) {
        const double sr = s(r);
        rem_sq += sr * sr;
        out.tail_normf[static_cast<std::size_t>(r)] = std::sqrt(rem_sq);
    }
    return out;
}

/// Build a complete `StencilEntry1D` from `GaussianConvolution1DCache<double>` at reference level 0.
inline StencilEntry1D build_stencil_entry_1d(int k, double gamma, double eps) {
    MADNESS_CHECK(k >= 1);
    MADNESS_CHECK(gamma > 0.0);

    StencilEntry1D entry;
    entry.gamma = gamma;
    entry.k = k;

    const long kk = static_cast<long>(k);
    const long twok = 2L * kk;
    const Slice s_lo(0, kk - 1);
    const Slice s_hi(kk, twok - 1);
    const Slice s_all(0, twok - 1);

    auto gc = GaussianConvolution1DCache<double>::get(k, gamma, 0, LatticeRange());

    for (int disp = -StencilEntry1D::kMaxDisp; disp <= StencilEntry1D::kMaxDisp; ++disp) {
        const std::size_t slot = static_cast<std::size_t>(disp + StencilEntry1D::kMaxDisp);
        const ConvolutionData1D<double>* ns = gc->nonstandard(0, static_cast<Translation>(disp));

        Tensor<double> R_mat(twok, twok);
        Tensor<double> R_2k_k_mat(twok, kk);
        Tensor<double> R_k_2k_mat(kk, twok);
        Tensor<double> T_mat(kk, kk);

        if (ns != nullptr && ns->R.size() == twok * twok) {
            R_mat = copy(ns->R);
            R_2k_k_mat = copy(R_mat(s_all, s_lo));
            R_k_2k_mat = copy(R_mat(s_lo, s_all));

            entry.normf_R_full[slot] = R_mat.normf();
            entry.normf_R_to_s[slot] = R_2k_k_mat.normf();
            entry.normf_R_s_to_all[slot] = R_k_2k_mat.normf();
            entry.normf_R_ss[slot] = R_mat(s_lo, s_lo).normf();
            entry.normf_R_d_to_all[slot] = R_mat(s_hi, s_all).normf();
            entry.normf_R_ds[slot] = R_mat(s_hi, s_lo).normf();

            entry.RU[slot] = ns->RU;
            entry.RVT[slot] = ns->RVT;
            entry.Rs[slot] = ns->Rs;
            entry.Rnorm[slot] = ns->Rnorm;
        }

        if (ns != nullptr && ns->T.size() == kk * kk) {
            T_mat = copy(ns->T);
            entry.normf_T[slot] = T_mat.normf();
            entry.TU[slot] = ns->TU;
            entry.TVT[slot] = ns->TVT;
            entry.Ts[slot] = ns->Ts;
            entry.Tnorm[slot] = ns->Tnorm;
        }

        entry.svd_R[slot] = build_svd_1d(R_mat);
        entry.svd_R_2k_to_k[slot] = build_svd_1d(R_2k_k_mat);
        entry.svd_R_k_to_2k[slot] = build_svd_1d(R_k_2k_mat);
        entry.svd_T[slot] = build_svd_1d(T_mat);

        entry.R[slot] = std::move(R_mat);
        entry.R_2k_to_k[slot] = std::move(R_2k_k_mat);
        entry.R_k_to_2k[slot] = std::move(R_k_2k_mat);
        entry.T[slot] = std::move(T_mat);
    }

    // Precompute exact 3D shell norms (`sh in {0, 1, 2}`) for all 5 NS operator blocks.
    std::array<double, StencilEntry1D::kMaxDisp + 1> sum2_T{};
    std::array<double, StencilEntry1D::kMaxDisp + 1> sum2_4b1{};
    std::array<double, StencilEntry1D::kMaxDisp + 1> sum2_leaf{};
    std::array<double, StencilEntry1D::kMaxDisp + 1> sum2_int{};
    std::array<double, StencilEntry1D::kMaxDisp + 1> sum2_wav{};

    for (int dx = -StencilEntry1D::kMaxDisp; dx <= StencilEntry1D::kMaxDisp; ++dx) {
        for (int dy = -StencilEntry1D::kMaxDisp; dy <= StencilEntry1D::kMaxDisp; ++dy) {
            for (int dz = -StencilEntry1D::kMaxDisp; dz <= StencilEntry1D::kMaxDisp; ++dz) {
                const std::size_t sh = static_cast<std::size_t>(
                    std::max({std::abs(dx), std::abs(dy), std::abs(dz)}));

                const double nT = entry.opnorm_T(dx, dy, dz);
                const double n4b1 = entry.opnorm_4b1(dx, dy, dz);
                const double nLeaf = entry.opnorm_4b2_leaf(dx, dy, dz);
                const double nInt = entry.opnorm_4b2_int(dx, dy, dz);
                const double nWav = entry.opnorm_4b2_wav(dx, dy, dz);

                entry.shell_max_T[sh] = std::max(entry.shell_max_T[sh], nT);
                entry.shell_max_4b1[sh] = std::max(entry.shell_max_4b1[sh], n4b1);
                entry.shell_max_4b2_leaf[sh] = std::max(entry.shell_max_4b2_leaf[sh], nLeaf);
                entry.shell_max_4b2_int[sh] = std::max(entry.shell_max_4b2_int[sh], nInt);
                entry.shell_max_4b2_wav[sh] = std::max(entry.shell_max_4b2_wav[sh], nWav);

                sum2_T[sh] += nT * nT;
                sum2_4b1[sh] += n4b1 * n4b1;
                sum2_leaf[sh] += nLeaf * nLeaf;
                sum2_int[sh] += nInt * nInt;
                sum2_wav[sh] += nWav * nWav;
            }
        }
    }

    for (std::size_t sh = 0; sh <= static_cast<std::size_t>(StencilEntry1D::kMaxDisp); ++sh) {
        entry.shell_l2_T[sh] = std::sqrt(sum2_T[sh]);
        entry.shell_l2_4b1[sh] = std::sqrt(sum2_4b1[sh]);
        entry.shell_l2_4b2_leaf[sh] = std::sqrt(sum2_leaf[sh]);
        entry.shell_l2_4b2_int[sh] = std::sqrt(sum2_int[sh]);
        entry.shell_l2_4b2_wav[sh] = std::sqrt(sum2_wav[sh]);
    }

    // Precompute exact 3D wavelet-octant Frobenius norm of a unit constant scaling field `e_{000}`
    // summed over `[-m_bw, +m_bw]^3` for `m_bw in {0, 1, 2}`:
    //   `v_1d[i_out] = \sum_{|disp| <= m_bw} R[disp](0, i_out)` (`i_out in [0, 2k)`),
    //   `||Q^j (v_1d^{\otimes 3})||_F = sqrt(max(0, (||v_1d||_2^2)^3 - (||v_1d[0..k)||_2^2)^3))`.
    for (int m_bw = 0; m_bw <= StencilEntry1D::kMaxDisp; ++m_bw) {
        double s_sq = 0.0;
        double all_sq = 0.0;
        for (long i_out = 0; i_out < twok; ++i_out) {
            double v = 0.0;
            for (int disp = -m_bw; disp <= m_bw; ++disp) {
                const std::size_t slot = static_cast<std::size_t>(disp + StencilEntry1D::kMaxDisp);
                v += entry.R[slot](0L, i_out);
            }
            const double v2 = v * v;
            all_sq += v2;
            if (i_out < kk) {
                s_sq += v2;
            }
        }
        const double all_cube = all_sq * all_sq * all_sq;
        const double s_cube = s_sq * s_sq * s_sq;
        entry.const_wav_3d_norm[static_cast<std::size_t>(m_bw)] =
            std::sqrt(std::max(0.0, all_cube - s_cube));
    }

    // Static ceiling half-bandwidths (finer per-level decision made by `select_bw_*`).
    const double base_tol = 1e-4 * eps;
    entry.m_max_T = (gamma >= 2.0 && entry.shell_l2_T[2] > base_tol) ? 2 : 1;
    const double ns_sh2 = std::max(entry.shell_l2_4b1[2], entry.shell_l2_4b2_int[2]);
    entry.m_max_ns = (ns_sh2 > base_tol) ? 2 : 1;
    return entry;
}

/// Build the root-level `(0, 0)` `k x k` scaling operator matrix oriented for `mTxmq`
/// (`T0_T(in_q, out_i) = r^0_0(out_i, in_q)`), avoiding child-level construction on ultra-diffuse terms.
inline Tensor<double> build_root_T0_transposed(int k, double alpha_sim) {
    MADNESS_CHECK(k >= 1);
    const long kk = static_cast<long>(k);
    Tensor<double> T0_T(kk, kk);
    if (!(alpha_sim > 0.0)) {
        return T0_T;
    }

    auto gc = GaussianConvolution1DCache<double>::get(k, alpha_sim, 0, LatticeRange());
    const Tensor<double> r0 = gc->rnlij(0, 0);
    if (r0.size() == kk * kk) {
        T0_T = transpose(r0);
    }
    return T0_T;
}

/// Thread-safe universal cache of `W * M` `StencilEntry1D` tables keyed by
/// `(k, W, M, eps, gamma_ref(0,0), gamma_ref(W-1,M-1))` — completely independent of
/// `mu_user` and `cell_width_L`, so constructing BSH operators across changing orbital
/// energies in SCF reuses the pre-SVD'd `W * M` stencil table in `O(1)` (`0` LAPACK calls).
inline std::shared_ptr<const std::vector<StencilEntry1D>>
get_or_build_universal_stencil_entries(const WindowedQuadraturePlan& plan) {
    struct UniversalStencilCacheRecord {
        int k = 0;
        int W = 0;
        int M = 0;
        double eps = 0.0;
        double gamma_00 = 0.0;
        double gamma_hi = 0.0;
        std::shared_ptr<const std::vector<StencilEntry1D>> entries;
    };

    const int k = plan.k();
    const int W = plan.window_width_W();
    const int M = plan.M();
    const double eps = plan.eps();
    const double g00 = plan.gamma_ref(0, 0);
    const double ghi = plan.gamma_ref(W - 1, M - 1);

    static std::mutex s_mtx;
    static std::vector<UniversalStencilCacheRecord> s_cache;

    std::lock_guard<std::mutex> lock(s_mtx);
    for (const auto& rec : s_cache) {
        if (rec.k == k && rec.W == W && rec.M == M && rec.eps == eps &&
            rec.gamma_00 == g00 && rec.gamma_hi == ghi) {
            return rec.entries;
        }
    }

    auto built = std::make_shared<std::vector<StencilEntry1D>>(
        static_cast<std::size_t>(W) * static_cast<std::size_t>(M));
    for (int w_eff = 0; w_eff < W; ++w_eff) {
        for (int m = 0; m < M; ++m) {
            const std::size_t idx =
                static_cast<std::size_t>(w_eff) * static_cast<std::size_t>(M) +
                static_cast<std::size_t>(m);
            (*built)[idx] = build_stencil_entry_1d(k, plan.gamma_ref(w_eff, m), eps);
        }
    }

    if (s_cache.size() >= 32) {
        s_cache.erase(s_cache.begin());
    }
    s_cache.push_back(UniversalStencilCacheRecord{k, W, M, eps, g00, ghi, built});
    return built;
}

} // namespace windowed_stencil_detail

/// Universal scale-locked 1D stencil cache for a `WindowedQuadraturePlan`.
///
/// Exploiting dyadic scale-locking `beta = 4^{1/M}`, every dyadic quadrature term `p`
/// inside its active window `j in [j_0(p), j_1(p)]` has dimensionless stiffness
///   `gamma_{p, j} = plan.gamma_ref(w_eff, m)`
/// determined solely by `w_eff = j - (term.n_spiky + 1) in [0, W - 1]` and
/// `m = term.sub_octave_m in [0, M - 1]`. Consequently only `W * M` (~15-35) 1D
/// stencil entries are needed across arbitrarily deep 3D nuclear trees (`n_max >= 30`),
/// and are shared universally across all `mu_user` and `cell_width_L` configurations.
class WindowedStencilCache1D {
public:
    explicit WindowedStencilCache1D(const WindowedQuadraturePlan& plan)
        : k_(plan.k())
        , W_(plan.window_width_W())
        , M_(plan.M())
        , eps_(plan.eps()) {
        MADNESS_CHECK(k_ >= 1);
        MADNESS_CHECK(W_ >= 1);
        MADNESS_CHECK(M_ >= 1);

        entries_ptr_ = windowed_stencil_detail::get_or_build_universal_stencil_entries(plan);

        // Precompute root T^0_{p, 00} matrices for ultra-diffuse terms (n_diff < 0)
        // and precompute any off-dyadic-grid active window entries (if ever non-empty)
        // so runtime execution is 100% lock-free and allocation-free.
        const auto& terms = plan.terms();
        root_T0_by_term_.resize(terms.size());
        const double ln2 = std::log(2.0);
        for (std::size_t idx = 0; idx < terms.size(); ++idx) {
            const WindowedQuadratureTerm& term = terms[idx];
            if (term.weight_scaled == 0.0) continue;

            if (term.n_diff < 0) {
                const double alpha_sim = std::exp(2.0 * term.t_p);
                root_T0_by_term_[idx] =
                    windowed_stencil_detail::build_root_T0_transposed(k_, alpha_sim);
            } else if (!term.on_dyadic_grid) {
                for (int j = term.j0; j <= term.n_diff; ++j) {
                    const double gamma_pj =
                        std::exp(2.0 * term.t_p - 2.0 * static_cast<double>(j) * ln2);
                    offgrid_entries_.push_back(OffGridEntry{
                        idx,
                        j,
                        windowed_stencil_detail::build_stencil_entry_1d(k_, gamma_pj, eps_)});
                }
            }
        }
    }

    int k() const noexcept { return k_; }
    int window_width_W() const noexcept { return W_; }
    int M() const noexcept { return M_; }
    double eps() const noexcept { return eps_; }

    /// Lookup universal 1D stencil entry for `(w_eff, m)` with `0 <= w_eff < W`, `0 <= m < M`.
    const StencilEntry1D& entry(int w_eff, int m) const {
        MADNESS_CHECK(w_eff >= 0 && w_eff < W_);
        MADNESS_CHECK(m >= 0 && m < M_);
        return (*entries_ptr_)[index(w_eff, m)];
    }

    /// Lookup 1D stencil entry for term `term_idx` (`plan.terms()[term_idx]`) at active level `j`.
    const StencilEntry1D& entry_for_term(std::size_t term_idx,
                                         const WindowedQuadratureTerm& term,
                                         int j) const {
        if (term.on_dyadic_grid) {
            const int w_eff = j - (term.n_spiky + 1);
            return entry(w_eff, term.sub_octave_m);
        }
        for (const auto& og : offgrid_entries_) {
            if (og.term_idx == term_idx && og.level == j) {
                return og.stencil;
            }
        }
        MADNESS_EXCEPTION("WindowedStencilCache1D::entry_for_term: missing off-grid stencil entry", 1);
        return (*entries_ptr_)[0];
    }

    /// Pre-transposed `k x k` root operator `T^0_{p, 00}` for ultra-diffuse term `term_idx` (`n_diff < 0`).
    const Tensor<double>& root_T0_transposed(std::size_t term_idx) const {
        MADNESS_CHECK(term_idx < root_T0_by_term_.size());
        return root_T0_by_term_[term_idx];
    }

private:
    struct OffGridEntry {
        std::size_t term_idx = 0;
        int level = 0;
        StencilEntry1D stencil;
    };

    std::size_t index(int w_eff, int m) const noexcept {
        return static_cast<std::size_t>(w_eff) * static_cast<std::size_t>(M_) +
               static_cast<std::size_t>(m);
    }

    int k_ = 0;
    int W_ = 0;
    int M_ = 0;
    double eps_ = 0.0;
    std::shared_ptr<const std::vector<StencilEntry1D>> entries_ptr_;
    std::vector<Tensor<double>> root_T0_by_term_;
    std::vector<OffGridEntry> offgrid_entries_;
};

} // namespace madness

#endif // MADNESS_MRA_WINDOWED_NS_STENCIL_H__INCLUDED
