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

#ifndef MADNESS_MRA_WINDOWED_QUADRATURE_H__INCLUDED
#define MADNESS_MRA_WINDOWED_QUADRATURE_H__INCLUDED

/// \file windowed_quadrature.h
/// \brief Two-threshold scale separation, dyadic-locked Coulomb/BSH/Yukawa
///        quadratures with discrete tail completion, overflow-safe log-scale
///        Stieltjes--Golub--Welsch compression, and parity-sparse Legendre
///        local moment operators for Windowed Non-Standard (NS) 3D convolutions.

#include <madness/madness_config.h>
#include <madness/constants.h>
#include <madness/tensor/aligned.h>
#include <madness/tensor/mTxmq.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>
#include <madness/world/madness_exception.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <mutex>
#include <utility>
#include <vector>

namespace madness {

/// Kernel family and physical normalization convention for Windowed NS operators.
enum class WindowedKernelKind {
    Coulomb,   ///< Physical kernel 1 / r (matches MADNESS CoulombOperator).
    BSH,       ///< Physical kernel exp(-mu * r) / (4 * pi * r) (matches BSHOperator).
    YukawaRaw  ///< Physical kernel exp(-mu * r) / r (unnormalized Yukawa).
};

/// Configuration parameters for constructing a WindowedQuadraturePlan.
struct WindowedQuadratureConfig {
    int k = 8;                                           ///< Legendre scaling order k >= 1.
    double eps = 1e-8;                                   ///< Target relative accuracy in (0, 1).
    double mu_user = 0.0;                                ///< Screening parameter in user coordinates (>= 0).
    double cell_width_L = 1.0;                           ///< Cubic simulation cell side width L > 0.
    int max_tree_level = 30;                             ///< Maximum tree level n_max for precomputation.
    int oversampling_M = 0;                              ///< Terms per dyadic octave (0 => auto from eps).
    int gc_rank_Q = 0;                                   ///< Gauss-Christoffel rank Q (0 => auto in {4, 5}).
    WindowedKernelKind kind = WindowedKernelKind::Coulomb;
};

/// Metadata and weights for a single Gaussian term in the dyadic-locked expansion.
struct WindowedQuadratureTerm {
    int p = 0;                     ///< Dyadic index (p_min - 1 for synthetic low-tail term).
    bool on_dyadic_grid = true;    ///< False only for the synthetic discrete low-tail term (mu = 0).
    int sub_octave_m = 0;          ///< Euclidean residue ((p % M) + M) % M in [0, M - 1].
    double t_p = 0.0;              ///< Log-exponent parameter: alpha_sim = exp(2 * t_p).
    double coeff_sim = 0.0;        ///< Unscaled unit-cube radial weight c_p in sum c_p exp(-alpha_sim r^2).
    double weight_sim = 0.0;       ///< Unscaled unit-cube 3D Gaussian operator weight W_p = c_p (pi/alpha_sim)^{3/2}.
    double weight_scaled = 0.0;    ///< Scaled 3D operator weight op_scale * weight_sim.
    int n_spiky = 0;               ///< floor((2*t_p - ln(gamma_spiky)) / (2*ln 2)).
    int n_diff = 0;                ///< floor((2*t_p - ln(gamma_diff)) / (2*ln 2)).
    int j0 = 0;                    ///< Base NS window level max(0, n_spiky + 1).
    int j1 = -1;                   ///< Top NS window level min(n_max - 1, n_diff).
};

/// Single compressed Gauss-Christoffel node at tree level n.
struct GCNode {
    double gamma_tilde = 0.0;      ///< Compressed dimensionless stiffness >= gamma_spiky(k).
    double weight_hat = 0.0;       ///< Overflow-free normalized weight exp(2*n*ln 2) * W_tilde_scaled.
    double weight_scaled = 0.0;    ///< Level-n 3D weight exp(-2*n*ln 2) * weight_hat (includes op_scale).
    Tensor<double> MT;             ///< k x k lower-triangular transpose(M(gamma_tilde)) for fast_transform/mTxmq.
};

/// Compressed Gauss-Christoffel rule for all spiky Gaussians (n_spiky(p) >= level) at a tree level.
struct GCLevelRule {
    int level = 0;
    std::vector<GCNode> nodes;     ///< Up to Q compressed nodes (empty only if spiky set is empty).
};

namespace windowed_quad_detail {

/// Unconditionally stable Miller backward recurrence + Parseval normalization for
/// S_k(z) = sum_{n=k}^{infty} (2n + 1) * j_n(z)^2 on z > 0, immune to small-z cancellation.
inline double spherical_bessel_tail_sum(double z, int k) noexcept {
    if (z <= 0.0) return 0.0;
    const double z_cbrt = std::cbrt(z);
    const int n_est = static_cast<int>(std::ceil(z + 30.0 + 5.0 * z_cbrt));
    const int n_top = std::max(k + 30, n_est);

    double f_next = 0.0;
    double f_curr = 1.0;
    double sum_all = static_cast<double>(2 * n_top + 1) * (f_curr * f_curr);
    double sum_tail = (n_top >= k) ? sum_all : 0.0;

    for (int n = n_top; n >= 1; --n) {
        double f_prev = (static_cast<double>(2 * n + 1) / z) * f_curr - f_next;
        if (std::abs(f_prev) > 1e100) {
            constexpr double rescale = 1e-100;
            constexpr double rescale_sq = rescale * rescale;
            f_prev *= rescale;
            f_curr *= rescale;
            sum_all *= rescale_sq;
            sum_tail *= rescale_sq;
        }
        const double term = static_cast<double>(2 * (n - 1) + 1) * (f_prev * f_prev);
        sum_all += term;
        if (n - 1 >= k) {
            sum_tail += term;
        }
        f_next = f_curr;
        f_curr = f_prev;
    }

    if (sum_all <= 0.0) return 0.0;
    const double ratio = sum_tail / sum_all;
    return std::clamp(ratio, 0.0, 1.0);
}

/// Evaluate the exact Plancherel bound (Eq. 5):
///   F_k(gamma) = 2 * sup_{z > 0} exp(-z^2 / gamma) * sqrt(sum_{n=k}^{infty} (2n+1) j_n(z)^2).
inline double plancherel_bound_value(double gamma, int k) noexcept {
    if (gamma <= 0.0) return 0.0;
    const double z_peak_est = std::sqrt(std::max(1e-300, 0.5 * static_cast<double>(k) * gamma));
    const double log_z_center = std::log(z_peak_est);
    const double u_min = log_z_center - 3.0;
    const double u_max = log_z_center + 2.5;

    auto eval_at_log_z = [gamma, k](double u) noexcept -> double {
        const double z = std::exp(u);
        const double expo = (z * z) / gamma;
        if (expo > 700.0) return 0.0;
        const double tail = spherical_bessel_tail_sum(z, k);
        if (tail <= 0.0) return 0.0;
        return 2.0 * std::exp(-expo) * std::sqrt(tail);
    };

    constexpr int n_grid = 40;
    const double du = (u_max - u_min) / static_cast<double>(n_grid);
    int best_idx = 0;
    double best_val = -1.0;
    for (int i = 0; i <= n_grid; ++i) {
        const double v = eval_at_log_z(u_min + static_cast<double>(i) * du);
        if (v > best_val) {
            best_val = v;
            best_idx = i;
        }
    }

    double lo = u_min + static_cast<double>(std::max(0, best_idx - 1)) * du;
    double hi = u_min + static_cast<double>(std::min(n_grid, best_idx + 1)) * du;
    const double inv_phi = (std::sqrt(5.0) - 1.0) * 0.5;
    double c = hi - inv_phi * (hi - lo);
    double d = lo + inv_phi * (hi - lo);
    double fc = eval_at_log_z(c);
    double fd = eval_at_log_z(d);

    for (int iter = 0; iter < 60; ++iter) {
        if (fc < fd) {
            lo = c;
            c = d;
            fc = fd;
            d = lo + inv_phi * (hi - lo);
            fd = eval_at_log_z(d);
        } else {
            hi = d;
            d = c;
            fd = fc;
            c = hi - inv_phi * (hi - lo);
            fc = eval_at_log_z(c);
        }
    }
    return std::max(best_val, std::max(fc, fd));
}

/// Euclidean integer modulo in [0, M - 1] and corresponding floor quotient q = floor(p / M).
inline std::pair<int, int> divmod_euclidean(int p, int M) noexcept {
    int m = p % M;
    if (m < 0) m += M;
    const int q = (p - m) / M;
    return {q, m};
}

} // namespace windowed_quad_detail

/// Parity-sparse 1D Legendre polynomial local moment operator
///   M(gamma) = sum_{r=0}^{m_k - 1} (gamma^{-r} / r!) D_2^r
/// and zero-allocation 3D tensor-product contractions.
///
/// Indexing convention:
///   - `make_M(gamma)` returns the k x k **upper-triangular** matrix M(out_deg, in_deg)
///     in standard column-vector notation (hat_s = M * s), nonzero only for
///     in_deg >= out_deg and (in_deg - out_deg) even.
///   - `make_MT(gamma)` returns the k x k **lower-triangular** transpose
///     MT(in_deg, out_deg) = M(out_deg, in_deg), which is the orientation expected
///     by MADNESS `transform`, `fast_transform`, and `mTxmq` (all of which contract
///     the first index of the transformation matrix).
class LegendreMomentOp {
public:
    explicit LegendreMomentOp(int k)
        : k_(k)
        , m_k_((k - 1) / 2 + 1)
        , d2_scaled_powers_(static_cast<std::size_t>(m_k_)) {
        MADNESS_CHECK(k_ >= 1);
        const long kk = static_cast<long>(k_);

        Tensor<double> D2(kk, kk);
        D2 = 0.0;
        for (int i = 0; i < k_; ++i) {
            const double wi = std::sqrt(static_cast<double>(2 * i + 1));
            for (int q = i + 2; q < k_; q += 2) {
                const double wq = std::sqrt(static_cast<double>(2 * q + 1));
                D2(i, q) = 0.5 * wi * wq * static_cast<double>((q - i) * (i + q + 1));
            }
        }

        Tensor<double> eye(kk, kk);
        eye = 0.0;
        for (int i = 0; i < k_; ++i) {
            eye(i, i) = 1.0;
        }
        d2_scaled_powers_[0] = eye;

        Tensor<double> prev = eye;
        for (int r = 1; r < m_k_; ++r) {
            Tensor<double> next(kk, kk);
            next = 0.0;
            const double inv_r = 1.0 / static_cast<double>(r);
            for (int i = 0; i < k_; ++i) {
                for (int q = i + 2 * r; q < k_; q += 2) {
                    double acc = 0.0;
                    for (int mid = i + 2 * (r - 1); mid <= q - 2; mid += 2) {
                        acc += prev(i, mid) * D2(mid, q);
                    }
                    next(i, q) = acc * inv_r;
                }
            }
            d2_scaled_powers_[static_cast<std::size_t>(r)] = next;
            prev = next;
        }
    }

    int order() const noexcept { return k_; }

    int num_even_derivs() const noexcept { return m_k_; }

    /// Fill a contiguous row-major k x k buffer with MT(in_deg, out_deg) = M(out_deg, in_deg)
    /// without heap allocation.
    void fill_MT_buffer(double gamma, double* mt_buf) const {
        MADNESS_CHECK(gamma > 0.0);
        const std::size_t n2 = static_cast<std::size_t>(k_) * static_cast<std::size_t>(k_);
        for (std::size_t idx = 0; idx < n2; ++idx) {
            mt_buf[idx] = 0.0;
        }
        for (int i = 0; i < k_; ++i) {
            mt_buf[static_cast<std::size_t>(i) * static_cast<std::size_t>(k_) + static_cast<std::size_t>(i)] = 1.0;
        }

        const double inv_gamma = 1.0 / gamma;
        double inv_gamma_pow = 1.0;
        for (int r = 1; r < m_k_; ++r) {
            inv_gamma_pow *= inv_gamma;
            const double* p_pow = d2_scaled_powers_[static_cast<std::size_t>(r)].ptr();
            for (int out_i = 0; out_i < k_; ++out_i) {
                for (int in_q = out_i + 2 * r; in_q < k_; in_q += 2) {
                    const double val = inv_gamma_pow *
                        p_pow[static_cast<std::size_t>(out_i) * static_cast<std::size_t>(k_) +
                              static_cast<std::size_t>(in_q)];
                    mt_buf[static_cast<std::size_t>(in_q) * static_cast<std::size_t>(k_) +
                           static_cast<std::size_t>(out_i)] += val;
                }
            }
        }
    }

    /// Return k x k upper-triangular M(out_deg, in_deg) in standard matrix-vector convention.
    Tensor<double> make_M(double gamma) const {
        MADNESS_CHECK(gamma > 0.0);
        const long kk = static_cast<long>(k_);
        Tensor<double> M(kk, kk);
        M = 0.0;
        for (int i = 0; i < k_; ++i) {
            M(i, i) = 1.0;
        }

        const double inv_gamma = 1.0 / gamma;
        double inv_gamma_pow = 1.0;
        for (int r = 1; r < m_k_; ++r) {
            inv_gamma_pow *= inv_gamma;
            const Tensor<double>& P = d2_scaled_powers_[static_cast<std::size_t>(r)];
            for (int i = 0; i < k_; ++i) {
                for (int q = i + 2 * r; q < k_; q += 2) {
                    M(i, q) += inv_gamma_pow * P(i, q);
                }
            }
        }
        return M;
    }

    /// Return k x k lower-triangular MT(in_deg, out_deg) = transpose(M(gamma))
    /// ready for MADNESS `transform`, `fast_transform`, and `mTxmq`.
    Tensor<double> make_MT(double gamma) const {
        const long kk = static_cast<long>(k_);
        Tensor<double> MT(kk, kk);
        fill_MT_buffer(gamma, MT.ptr());
        return MT;
    }

    /// Apply M^{\otimes 3} using a precomputed k x k lower-triangular MT Tensor via MADNESS `fast_transform`.
    template <typename T>
    void apply_MT_3d(const Tensor<double>& MT,
                     const Tensor<T>& s_in,
                     Tensor<T>& s_out,
                     Tensor<T>& work1,
                     Tensor<T>& work2) const {
        const long kk = static_cast<long>(k_);
        MADNESS_CHECK(s_in.ndim() == 3 && s_in.dim(0) == kk && s_in.dim(1) == kk && s_in.dim(2) == kk);
        MADNESS_CHECK(MT.ndim() == 2 && MT.dim(0) == kk && MT.dim(1) == kk && MT.iscontiguous());

        ensure_cube_buffer(s_out, kk);
        ensure_cube_buffer(work1, kk);

        if (s_in.iscontiguous() && s_in.ptr() != s_out.ptr() && s_in.ptr() != work1.ptr()) {
            fast_transform(s_in, MT, s_out, work1);
        } else {
            work2 = copy(s_in);
            fast_transform(work2, MT, s_out, work1);
        }
    }

    /// Apply M^{\otimes 3} using a precomputed k x k lower-triangular MT buffer via `mTxmq` / `fast_transform`.
    template <typename T>
    void apply_MT_buffer_3d(const double* mt_ptr,
                            const Tensor<T>& s_in,
                            Tensor<T>& s_out,
                            Tensor<T>& work1,
                            Tensor<T>& work2) const {
        const long kk = static_cast<long>(k_);
        MADNESS_CHECK(s_in.ndim() == 3 && s_in.dim(0) == kk && s_in.dim(1) == kk && s_in.dim(2) == kk);

        if constexpr (!std::is_same_v<T, double>) {
            const long mdims[2] = {kk, kk};
            Tensor<double> MT(2L, mdims, false);
            std::copy_n(mt_ptr, static_cast<std::size_t>(kk * kk), MT.ptr());
            apply_MT_3d(MT, s_in, s_out, work1, work2);
            return;
        } else {
            ensure_cube_buffer(s_out, kk);
            ensure_cube_buffer(work1, kk);

            const T* src_ptr = nullptr;
            if (s_in.iscontiguous() && s_in.ptr() != s_out.ptr() && s_in.ptr() != work1.ptr()) {
                src_ptr = s_in.ptr();
            } else {
                work2 = copy(s_in);
                src_ptr = work2.ptr();
            }

            const long dimi = kk * kk;
            mTxmq(dimi, kk, kk, s_out.ptr(), src_ptr, mt_ptr);
            mTxmq(dimi, kk, kk, work1.ptr(), s_out.ptr(), mt_ptr);
            mTxmq(dimi, kk, kk, s_out.ptr(), work1.ptr(), mt_ptr);
        }
    }

    /// Apply M(gamma)^{\otimes 3} to a 3D k x k x k tensor with caller-provided scratch buffers
    /// and zero heap allocations on steady-state calls.
    template <typename T>
    void apply_M_3d(double gamma,
                    const Tensor<T>& s_in,
                    Tensor<T>& s_out,
                    Tensor<T>& work1,
                    Tensor<T>& work2) const {
        if constexpr (!std::is_same_v<T, double>) {
            const Tensor<double> MT = make_MT(gamma);
            apply_MT_3d(MT, s_in, s_out, work1, work2);
        } else {
            constexpr int kStackMax = 64;
            if (k_ <= kStackMax) {
                alignas(64) std::array<double, kStackMax * kStackMax> mt_stack{};
                fill_MT_buffer(gamma, mt_stack.data());
                apply_MT_buffer_3d(mt_stack.data(), s_in, s_out, work1, work2);
            } else {
                std::vector<double> mt_heap(static_cast<std::size_t>(k_) * static_cast<std::size_t>(k_), 0.0);
                fill_MT_buffer(gamma, mt_heap.data());
                apply_MT_buffer_3d(mt_heap.data(), s_in, s_out, work1, work2);
            }
        }
    }

    /// Convenience overload returning a newly allocated k x k x k Tensor<T>.
    template <typename T>
    Tensor<T> apply_M_3d(double gamma, const Tensor<T>& s_in) const {
        const long kk = static_cast<long>(k_);
        const long dims[3] = {kk, kk, kk};
        Tensor<T> s_out(3L, dims, false);
        Tensor<T> work1(3L, dims, false);
        Tensor<T> work2;
        apply_M_3d(gamma, s_in, s_out, work1, work2);
        return s_out;
    }

    static void ensure_cube_buffer_3d(Tensor<double>& buf, long kk) {
        ensure_cube_buffer(buf, kk);
    }
    template <typename T>
    static void ensure_cube_buffer(Tensor<T>& buf, long kk) {
        if (buf.ndim() != 3 || buf.dim(0) != kk || buf.dim(1) != kk || buf.dim(2) != kk || !buf.iscontiguous()) {
            const long dims[3] = {kk, kk, kk};
            buf = Tensor<T>(3L, dims, false);
        }
    }

private:

    int k_;
    int m_k_;
    std::vector<Tensor<double>> d2_scaled_powers_;
};

/// Precomputed two-threshold Windowed NS quadrature plan, universal stiffness reference table,
/// and per-level overflow-safe Gauss-Christoffel compressed spiky operators.
class WindowedQuadraturePlan {
public:
    explicit WindowedQuadraturePlan(const WindowedQuadratureConfig& cfg)
        : k_(cfg.k)
        , eps_(cfg.eps)
        , mu_user_((cfg.kind == WindowedKernelKind::Coulomb) ? 0.0 : cfg.mu_user)
        , cell_width_L_(cfg.cell_width_L)
        , mu_sim_(mu_user_ * cfg.cell_width_L)
        , op_scale_(compute_op_scale(cfg.kind, cfg.cell_width_L))
        , gamma_spiky_(compute_gamma_spiky(cfg.k))
        , gamma_diff_(compute_gamma_diff(cfg.k, cfg.eps))
        , M_(select_oversampling_M(cfg.eps, cfg.oversampling_M))
        , gc_rank_Q_(select_gc_rank_Q(cfg.k, cfg.gc_rank_Q))
        , n_max_(cfg.max_tree_level)
        , kind_(cfg.kind)
        , moment_op_(cfg.k) {
        MADNESS_CHECK(k_ >= 1);
        MADNESS_CHECK(eps_ > 0.0 && eps_ < 1.0);
        MADNESS_CHECK(cell_width_L_ > 0.0);
        MADNESS_CHECK(mu_user_ >= 0.0);
        MADNESS_CHECK(n_max_ >= 0);

        const double ln2 = std::log(2.0);
        const double two_ln2 = 2.0 * ln2;
        delta_t_ = ln2 / static_cast<double>(M_);
        beta_ = std::exp(2.0 * delta_t_);

        const double log_gs = std::log(gamma_spiky_);
        const double log_gd = std::log(gamma_diff_);

        // Exact integer residue-class base levels so (w_eff, m) <-> gamma_ref(w_eff, m)
        // holds with zero 1-ULP floor drift across all octaves p in Z.
        s_base_spiky_.resize(static_cast<std::size_t>(M_));
        s_base_diff_.resize(static_cast<std::size_t>(M_));
        window_width_W_ = 1;
        for (int m = 0; m < M_; ++m) {
            const double frac_m = static_cast<double>(m) / static_cast<double>(M_);
            const int sb_spiky = static_cast<int>(std::floor(frac_m - log_gs / two_ln2));
            const int sb_diff = static_cast<int>(std::floor(frac_m - log_gd / two_ln2));
            s_base_spiky_[static_cast<std::size_t>(m)] = sb_spiky;
            s_base_diff_[static_cast<std::size_t>(m)] = sb_diff;
            window_width_W_ = std::max(window_width_W_, sb_diff - sb_spiky);
        }

        // Populate canonical reference stiffnesses gamma_ref(w_eff, m) for w_eff in [0, W-1], m in [0, M-1].
        gamma_ref_table_.resize(static_cast<std::size_t>(window_width_W_) * static_cast<std::size_t>(M_));
        for (int w_eff = 0; w_eff < window_width_W_; ++w_eff) {
            for (int m = 0; m < M_; ++m) {
                const double frac_m = static_cast<double>(m) / static_cast<double>(M_);
                const int sb_spiky = s_base_spiky_[static_cast<std::size_t>(m)];
                const double exponent = two_ln2 * (frac_m - static_cast<double>(sb_spiky + 1 + w_eff));
                gamma_ref_table_[static_cast<std::size_t>(w_eff) * static_cast<std::size_t>(M_) +
                                 static_cast<std::size_t>(m)] = std::exp(exponent);
            }
        }

        // Low-tail cutoff p_min and discrete low-tail completion (Eqs. 15-17).
        const double r_max = std::sqrt(3.0);
        const double eps_quad = 0.1 * eps_;
        if (mu_sim_ == 0.0) {
            const double r_max_7 = std::pow(r_max, 7);
            const double t_target = (1.0 / 7.0) *
                std::log((525.0 * constants::sqrt_pi * eps_quad) / (4.0 * r_max_7));
            p_min_ = static_cast<int>(std::floor(t_target / delta_t_));

            const double tp_min = static_cast<double>(p_min_) * delta_t_;
            const double pref = (2.0 * delta_t_) * constants::inv_sqrt_pi;
            const double M0 = pref * std::exp(1.0 * tp_min) / std::expm1(1.0 * delta_t_);
            const double M1 = pref * std::exp(3.0 * tp_min) / std::expm1(3.0 * delta_t_);
            const double M2 = pref * std::exp(5.0 * tp_min) / std::expm1(5.0 * delta_t_);

            const double alpha_low = M2 / M1;
            const double t_low = 0.5 * std::log(alpha_low);
            const double c_low = (M1 * M1) / M2;
            const double w_low_sim = c_low * std::pow(constants::pi / alpha_low, 1.5);
            w_const_sim_ = M0 - c_low;
            w_const_scaled_ = op_scale_ * w_const_sim_;

            WindowedQuadratureTerm low_term;
            low_term.p = p_min_ - 1;
            low_term.on_dyadic_grid = false;
            low_term.sub_octave_m = windowed_quad_detail::divmod_euclidean(p_min_ - 1, M_).second;
            low_term.t_p = t_low;
            low_term.coeff_sim = c_low;
            low_term.weight_sim = w_low_sim;
            low_term.weight_scaled = op_scale_ * w_low_sim;
            low_term.n_spiky = static_cast<int>(std::floor((2.0 * t_low - log_gs) / two_ln2));
            low_term.n_diff = static_cast<int>(std::floor((2.0 * t_low - log_gd) / two_ln2));
            low_term.j0 = std::max(0, low_term.n_spiky + 1);
            low_term.j1 = std::min(n_max_ - 1, low_term.n_diff);
            terms_.push_back(low_term);
        } else {
            w_const_sim_ = 0.0;
            w_const_scaled_ = 0.0;
            const double k_rmax = (mu_sim_ * r_max < 700.0) ? (std::exp(-mu_sim_ * r_max) / r_max) : 0.0;
            const double eps_abs = eps_quad * std::min(1.0 / r_max, std::max(1e-15, k_rmax));
            const double log_denom = 4.0 * std::log(std::max(100.0, 100.0 / std::min(eps_quad, eps_abs)));
            const double alpha_cut_mu = (mu_sim_ * mu_sim_) / log_denom;
            const double alpha_cut_unif = (constants::pi * eps_abs * eps_abs) / 4.0;
            const double alpha_min = std::max(alpha_cut_mu, alpha_cut_unif);
            p_min_ = static_cast<int>(std::floor(0.5 * std::log(alpha_min) / delta_t_));
        }

        // High-alpha cusp cutoff: smallest p_first such that n_spiky_grid(p_first) >= n_max_.
        int p_probe = static_cast<int>(std::floor((static_cast<double>(n_max_) * ln2 + 0.5 * log_gs) / delta_t_)) - 2;
        while (grid_n_spiky(p_probe) >= n_max_) {
            --p_probe;
        }
        while (grid_n_spiky(p_probe) < n_max_) {
            ++p_probe;
        }
        p_max_ = p_probe - 1;
        p_tail_ = p_max_ + 10 * M_;

        // Populate uniform dyadic-locked terms p in [p_min_, p_tail_].
        terms_.reserve(terms_.size() + static_cast<std::size_t>(std::max(0, p_tail_ - p_min_ + 1)));
        for (int p = p_min_; p <= p_tail_; ++p) {
            const auto [q_oct, m_sub] = windowed_quad_detail::divmod_euclidean(p, M_);
            const double tp = static_cast<double>(p) * delta_t_;

            double mu_arg = 0.0;
            if (mu_sim_ > 0.0) {
                if (-2.0 * tp > 700.0) {
                    mu_arg = std::numeric_limits<double>::infinity();
                } else {
                    mu_arg = 0.25 * mu_sim_ * mu_sim_ * std::exp(-2.0 * tp);
                }
            }

            double c_sim = 0.0;
            double w_sim = 0.0;
            if (std::isfinite(mu_arg) && mu_arg < 700.0) {
                if (tp - mu_arg > -740.0 && tp - mu_arg < 705.0) {
                    c_sim = (2.0 * delta_t_ * constants::inv_sqrt_pi) * std::exp(tp - mu_arg);
                }
                const double w_log_exp = -2.0 * tp - mu_arg;
                if (w_log_exp > -740.0 && w_log_exp < 705.0) {
                    w_sim = (2.0 * constants::pi * delta_t_) * std::exp(w_log_exp);
                }
            }

            WindowedQuadratureTerm term;
            term.p = p;
            term.on_dyadic_grid = true;
            term.sub_octave_m = m_sub;
            term.t_p = tp;
            term.coeff_sim = c_sim;
            term.weight_sim = w_sim;
            term.weight_scaled = op_scale_ * w_sim;
            term.n_spiky = q_oct + s_base_spiky_[static_cast<std::size_t>(m_sub)];
            term.n_diff = q_oct + s_base_diff_[static_cast<std::size_t>(m_sub)];
            term.j0 = std::max(0, term.n_spiky + 1);
            term.j1 = std::min(n_max_ - 1, term.n_diff);
            terms_.push_back(term);
        }

        // Precompute compressed Gauss-Christoffel rules lock-free for all levels 0 .. n_max_ + 12.
        const int precomp_top = n_max_ + 12;
        gc_rules_.resize(static_cast<std::size_t>(precomp_top + 1));
        for (int lvl = 0; lvl <= precomp_top; ++lvl) {
            gc_rules_[static_cast<std::size_t>(lvl)] = build_gc_rule_for_level(lvl);
        }
    }

    /// Exact spiky threshold gamma_spiky(k) = k^2 * (k + 1)^2 / 30 (Theorem 2.2 / Invariant 2).
    static double compute_gamma_spiky(int k) noexcept {
        const double kd = static_cast<double>(k);
        const double kp1 = static_cast<double>(k + 1);
        return (kd * kd * kp1 * kp1) / 30.0;
    }

    /// Exact spherical-Bessel / Plancherel diffuse threshold gamma_diff(k, eps) (Theorem 2.1, Eq. 5).
    static double compute_gamma_diff(int k, double eps) {
        MADNESS_CHECK(k >= 1);
        MADNESS_CHECK(eps > 0.0 && eps < 1.0);

        // Thread-safe tiny memoization cache for repeated operator constructions at same (k, eps).
        struct CacheEntry {
            int k;
            double eps;
            double val;
        };
        static std::mutex cache_mtx;
        static std::vector<CacheEntry> cache;
        {
            std::lock_guard<std::mutex> lock(cache_mtx);
            for (const auto& entry : cache) {
                if (entry.k == k && entry.eps == eps) {
                    return entry.val;
                }
            }
        }

        // Closed-form asymptotic initial seed from Eq. (5).
        double log_dfact = 0.0;
        for (int j = 1; j <= k; ++j) {
            log_dfact += std::log(static_cast<double>(2 * j + 1));
        }
        const double log_pref = std::log(2.0) + 0.5 * std::log(static_cast<double>(2 * k + 1)) - log_dfact;
        const double kd = static_cast<double>(k);
        const double e_const = std::exp(1.0);
        const double gamma_seed = (2.0 * e_const / kd) * std::exp((2.0 / kd) * (std::log(eps) - log_pref));

        double lo = std::log(std::max(1e-30, gamma_seed)) - 2.0;
        double hi = std::log(std::max(1e-30, gamma_seed)) + 2.0;
        while (windowed_quad_detail::plancherel_bound_value(std::exp(lo), k) > eps) {
            lo -= 2.0;
        }
        while (windowed_quad_detail::plancherel_bound_value(std::exp(hi), k) < eps) {
            hi += 2.0;
        }
        for (int iter = 0; iter < 64; ++iter) {
            const double mid = 0.5 * (lo + hi);
            if (windowed_quad_detail::plancherel_bound_value(std::exp(mid), k) < eps) {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        const double result = std::exp(0.5 * (lo + hi));

        {
            std::lock_guard<std::mutex> lock(cache_mtx);
            if (cache.size() < 256) {
                cache.push_back(CacheEntry{k, eps, result});
            }
        }
        return result;
    }

    int k() const noexcept { return k_; }
    double eps() const noexcept { return eps_; }
    double gamma_spiky() const noexcept { return gamma_spiky_; }
    double gamma_diff() const noexcept { return gamma_diff_; }
    int M() const noexcept { return M_; }
    int window_width_W() const noexcept { return window_width_W_; }
    int gc_rank_Q() const noexcept { return gc_rank_Q_; }
    int n_max() const noexcept { return n_max_; }
    int p_min() const noexcept { return p_min_; }
    int p_max() const noexcept { return p_max_; }
    int p_tail() const noexcept { return p_tail_; }
    double delta_t() const noexcept { return delta_t_; }
    double beta() const noexcept { return beta_; }
    double mu_user() const noexcept { return mu_user_; }
    double mu_sim() const noexcept { return mu_sim_; }
    double cell_width_L() const noexcept { return cell_width_L_; }
    double op_scale() const noexcept { return op_scale_; }
    double w_const_sim() const noexcept { return w_const_sim_; }
    double w_const_scaled() const noexcept { return w_const_scaled_; }
    WindowedKernelKind kind() const noexcept { return kind_; }

    /// Whether an analytic p > p_tail discrete remainder mass is included in GC compression.
    bool has_high_tail_remainder() const noexcept { return mu_sim_ == 0.0; }

    /// Exact analytic p > p_tail remainder spike {weight_scaled_rem, gamma_rem} at tree level `level`
    /// (zero weight when mu_sim > 0). Useful for Gate 1 uncompressed vs GC norm verification.
    std::pair<double, double> high_tail_remainder_node(int level) const noexcept {
        if (mu_sim_ != 0.0) return {0.0, 1.0};
        const double ln2 = std::log(2.0);
        const double tp_tail = static_cast<double>(p_tail_) * delta_t_;
        const double u_tail = std::exp(2.0 * static_cast<double>(level) * ln2 - 2.0 * tp_tail);
        const double wh_tail = op_scale_ * (2.0 * constants::pi * delta_t_) * u_tail;
        const double wh_rem = wh_tail / (beta_ - 1.0);
        const double u_rem = u_tail / (beta_ + 1.0);
        const double w_scaled_rem = std::exp(-2.0 * static_cast<double>(level) * ln2) * wh_rem;
        return {w_scaled_rem, 1.0 / u_rem};
    }

    /// Reference dimensionless stiffness for Milestone 2's universal (w_eff, m) 1D stencil cache:
    ///   w_eff = j - (term.n_spiky + 1) in [0, window_width_W() - 1],
    ///   m     = term.sub_octave_m      in [0, M() - 1].
    double gamma_ref(int w_eff, int m) const {
        MADNESS_CHECK(w_eff >= 0 && w_eff < window_width_W_);
        MADNESS_CHECK(m >= 0 && m < M_);
        return gamma_ref_table_[static_cast<std::size_t>(w_eff) * static_cast<std::size_t>(M_) +
                                static_cast<std::size_t>(m)];
    }

    /// Full ordered list of quadrature terms:
    ///   - If mu_sim == 0: synthetic off-grid low-tail term (p = p_min - 1) followed by p = p_min .. p_tail.
    ///   - If mu_sim > 0:  dyadic-grid terms p = p_min .. p_tail.
    const std::vector<WindowedQuadratureTerm>& terms() const noexcept { return terms_; }

    /// Compressed Gauss-Christoffel rule at tree level `level >= 0` (lock-free O(1) for level <= n_max() + 12).
    const GCLevelRule& gc_rule(int level) const {
        MADNESS_CHECK(level >= 0);
        const std::size_t idx = static_cast<std::size_t>(level);
        if (idx < gc_rules_.size()) {
            return gc_rules_[idx];
        }
        std::lock_guard<std::mutex> lock(extra_gc_mtx_);
        for (const auto& r : extra_gc_rules_) {
            if (r.level == level) return r;
        }
        extra_gc_rules_.push_back(build_gc_rule_for_level(level));
        return extra_gc_rules_.back();
    }

    const LegendreMomentOp& moment_op() const noexcept { return moment_op_; }

    /// Stage 2 leaf helper: accumulate sum_{q=1}^Q W_tilde_{n,q} M(gamma_tilde_{n,q})^{\otimes 3} s_leaf
    /// into `s_hat_accum` via MADNESS `fast_transform` + MKL `aligned_axpy` with zero heap allocations
    /// on steady-state calls.
    template <typename T>
    void accumulate_gc_leaf_3d(int level,
                               const Tensor<T>& s_leaf,
                               Tensor<T>& s_hat_accum,
                               Tensor<T>& tmp,
                               Tensor<T>& work1,
                               Tensor<T>& work2) const {
        const GCLevelRule& rule = gc_rule(level);
        if (rule.nodes.empty()) return;

        const long kk = static_cast<long>(k_);
        MADNESS_CHECK(s_leaf.ndim() == 3 &&
                      s_leaf.dim(0) == kk &&
                      s_leaf.dim(1) == kk &&
                      s_leaf.dim(2) == kk);
        MADNESS_CHECK(s_hat_accum.ndim() == 3 &&
                      s_hat_accum.dim(0) == kk &&
                      s_hat_accum.dim(1) == kk &&
                      s_hat_accum.dim(2) == kk);

        LegendreMomentOp::ensure_cube_buffer(tmp, kk);
        LegendreMomentOp::ensure_cube_buffer(work1, kk);

        const Tensor<T>* src_ptr = &s_leaf;
        if (!s_leaf.iscontiguous() || s_leaf.ptr() == tmp.ptr() || s_leaf.ptr() == work1.ptr()) {
            work2 = copy(s_leaf);
            src_ptr = &work2;
        }

        const long total_elems = kk * kk * kk;
        if (s_hat_accum.iscontiguous() && s_hat_accum.ptr() != tmp.ptr() && s_hat_accum.ptr() != work1.ptr()) {
            T* MADNESS_RESTRICT p_acc = s_hat_accum.ptr();
            for (const GCNode& node : rule.nodes) {
                if (node.weight_scaled == 0.0) continue;
                fast_transform(*src_ptr, node.MT, tmp, work1);
                aligned_axpy(total_elems, p_acc, tmp.ptr(), node.weight_scaled);
            }
        } else {
            for (const GCNode& node : rule.nodes) {
                if (node.weight_scaled == 0.0) continue;
                fast_transform(*src_ptr, node.MT, tmp, work1);
                s_hat_accum.gaxpy(T(1), tmp, T(node.weight_scaled));
            }
        }
    }

    /// Convenience overload allocating temporary buffers automatically.
    template <typename T>
    void accumulate_gc_leaf_3d(int level,
                               const Tensor<T>& s_leaf,
                               Tensor<T>& s_hat_accum) const {
        const long kk = static_cast<long>(k_);
        const long dims[3] = {kk, kk, kk};
        Tensor<T> tmp(3L, dims, false);
        Tensor<T> work1(3L, dims, false);
        Tensor<T> work2;
        accumulate_gc_leaf_3d(level, s_leaf, s_hat_accum, tmp, work1, work2);
    }

    /// Evaluate the unscaled unit-cube radial kernel sum at r_sim in (0, sqrt(3)]:
    ///   K_sim_unit(r_sim) = w_const_sim + sum_{term} c_p * exp(-alpha_sim * r_sim^2)
    /// (approximating 1 / r_sim for mu_sim == 0 and exp(-mu_sim * r_sim) / r_sim for mu_sim > 0).
    double eval_kernel_sim(double r_sim) const {
        MADNESS_CHECK(r_sim > 0.0);
        const double log_r2 = 2.0 * std::log(r_sim);
        double sum = w_const_sim_;
        for (const WindowedQuadratureTerm& term : terms_) {
            if (term.coeff_sim == 0.0) continue;
            const double log_arg = 2.0 * term.t_p + log_r2;
            if (log_arg < 45.0) {
                sum += term.coeff_sim * std::exp(-std::exp(log_arg));
            }
        }
        // If called at tiny r_sim on a shallow max_tree_level plan where exp(2*t_{p_tail})*r_sim^2 < 45,
        // include the analytic high-p dyadic grid continuation beyond p_tail.
        const double top_log_arg = 2.0 * static_cast<double>(p_tail_) * delta_t_ + log_r2;
        if (top_log_arg < 45.0) {
            for (int p = p_tail_ + 1;; ++p) {
                const double tp = static_cast<double>(p) * delta_t_;
                const double log_arg = 2.0 * tp + log_r2;
                if (log_arg >= 45.0) break;
                const double mu_arg = (mu_sim_ > 0.0 && -2.0 * tp < 700.0)
                    ? (0.25 * mu_sim_ * mu_sim_ * std::exp(-2.0 * tp))
                    : 0.0;
                if (mu_arg < 700.0) {
                    const double cp = (2.0 * delta_t_ * constants::inv_sqrt_pi) * std::exp(tp - mu_arg);
                    sum += cp * std::exp(-std::exp(log_arg));
                }
            }
        }
        return sum;
    }

    /// Evaluate the physical user-space radial kernel approximation at r_user > 0:
    ///   - Coulomb:   approx 1 / r_user
    ///   - BSH:       approx exp(-mu_user * r_user) / (4 * pi * r_user)
    ///   - YukawaRaw: approx exp(-mu_user * r_user) / r_user
    double eval_kernel_user(double r_user) const {
        MADNESS_CHECK(r_user > 0.0);
        const double r_sim = r_user / cell_width_L_;
        const double l_cubed = cell_width_L_ * cell_width_L_ * cell_width_L_;
        return (op_scale_ / l_cubed) * eval_kernel_sim(r_sim);
    }

private:
    static double compute_op_scale(WindowedKernelKind kind, double L) noexcept {
        const double l2 = L * L;
        if (kind == WindowedKernelKind::BSH) {
            return l2 / (4.0 * constants::pi);
        }
        return l2;
    }

    static int select_oversampling_M(double eps, int m_override) noexcept {
        if (m_override >= 1) return m_override;
        const double ln2 = std::log(2.0);
        const double pi2 = constants::pi * constants::pi;
        const double target = 0.1 * eps;
        for (int M = 2; M <= 8; ++M) {
            const double dt = ln2 / static_cast<double>(M);
            const double x = pi2 / dt;
            // Exact Proposition 4.1 leading error: 2 / sqrt(cosh(x)) = 2*sqrt(2)*exp(-x/2)/sqrt(1 + exp(-2x)).
            const double e_disc = 2.0 * std::sqrt(2.0) * std::exp(-0.5 * x) /
                std::sqrt(1.0 + ((2.0 * x < 700.0) ? std::exp(-2.0 * x) : 0.0));
            if (e_disc <= target) {
                return M;
            }
        }
        return 8;
    }

    static int select_gc_rank_Q(int k, int q_override) noexcept {
        if (q_override >= 1) return q_override;
        if (k <= 6) return 4;
        if (k <= 14) return 5;
        const int m_k = (k - 1) / 2 + 1;
        return std::max(5, (m_k + 4) / 2);
    }

    int grid_n_spiky(int p) const noexcept {
        const auto [q_oct, m_sub] = windowed_quad_detail::divmod_euclidean(p, M_);
        return q_oct + s_base_spiky_[static_cast<std::size_t>(m_sub)];
    }

    /// Overflow-safe log-scale Stieltjes orthogonalization + Golub-Welsch eigensolve via madness::syev.
    GCLevelRule build_gc_rule_for_level(int level) const {
        GCLevelRule rule;
        rule.level = level;

        const double ln2 = std::log(2.0);
        const double two_n_ln2 = 2.0 * static_cast<double>(level) * ln2;

        std::vector<double> u_pts;
        std::vector<double> w_hats;
        u_pts.reserve(terms_.size() + 1);
        w_hats.reserve(terms_.size() + 1);

        for (const WindowedQuadratureTerm& term : terms_) {
            if (term.n_spiky < level) continue;
            const double u_pn = std::exp(two_n_ln2 - 2.0 * term.t_p);
            if (u_pn <= 0.0) continue;

            double wh_pn = 0.0;
            if (term.on_dyadic_grid) {
                double mu_factor = 1.0;
                if (mu_sim_ > 0.0) {
                    if (-2.0 * term.t_p > 700.0) {
                        mu_factor = 0.0;
                    } else {
                        const double mu_arg = 0.25 * mu_sim_ * mu_sim_ * std::exp(-2.0 * term.t_p);
                        mu_factor = (mu_arg < 700.0) ? std::exp(-mu_arg) : 0.0;
                    }
                }
                wh_pn = op_scale_ * (2.0 * constants::pi * delta_t_) * u_pn * mu_factor;
            } else {
                // Off-grid term (only possible if n_spiky(low) >= level, never for typical eps).
                wh_pn = (-two_n_ln2 < 700.0) ? (std::exp(two_n_ln2) * term.weight_scaled) : 0.0;
            }

            if (wh_pn > 0.0) {
                u_pts.push_back(u_pn);
                w_hats.push_back(wh_pn);
            }
        }

        // Include exact analytic discrete high-tail remainder spike for p > p_tail when mu_sim == 0.
        if (mu_sim_ == 0.0) {
            const double tp_tail = static_cast<double>(p_tail_) * delta_t_;
            if (grid_n_spiky(p_tail_) >= level) {
                const double u_tail = std::exp(two_n_ln2 - 2.0 * tp_tail);
                const double wh_tail = op_scale_ * (2.0 * constants::pi * delta_t_) * u_tail;
                const double wh_rem = wh_tail / (beta_ - 1.0);
                const double u_rem = u_tail / (beta_ + 1.0);
                if (wh_rem > 0.0 && u_rem > 0.0) {
                    u_pts.push_back(u_rem);
                    w_hats.push_back(wh_rem);
                }
            }
        }

        if (u_pts.empty()) {
            return rule;
        }

        const std::size_t n_src = u_pts.size();
        int Q_eff = std::min(gc_rank_Q_, static_cast<int>(n_src));

        double u_max = 0.0;
        double w_sum = 0.0;
        for (std::size_t i = 0; i < n_src; ++i) {
            u_max = std::max(u_max, u_pts[i]);
            w_sum += w_hats[i];
        }
        if (u_max <= 0.0 || w_sum <= 0.0) {
            return rule;
        }

        // Rescale u -> s = u / u_max in (0, 1] and normalize weights to sum to 1
        // so Stieltjes/Lanczos polynomial norms remain O(1) at arbitrary k and Q.
        std::vector<double> s_pts(n_src);
        std::vector<double> w_norm(n_src);
        for (std::size_t i = 0; i < n_src; ++i) {
            s_pts[i] = u_pts[i] / u_max;
            w_norm[i] = w_hats[i] / w_sum;
        }

        std::vector<double> diag_a(static_cast<std::size_t>(Q_eff), 0.0);
        std::vector<double> sub_b(static_cast<std::size_t>(std::max(0, Q_eff - 1)), 0.0);
        std::vector<double> p_prev(n_src, 0.0);
        std::vector<double> p_curr(n_src, 1.0);
        std::vector<double> p_next(n_src, 0.0);

        double norm_sq_prev = 1.0;
        for (int j = 0; j < Q_eff; ++j) {
            double norm_sq = 0.0;
            double s_num = 0.0;
            for (std::size_t i = 0; i < n_src; ++i) {
                const double wp2 = w_norm[i] * p_curr[i] * p_curr[i];
                norm_sq += wp2;
                s_num += wp2 * s_pts[i];
            }
            if (norm_sq <= 1e-300) {
                Q_eff = j;
                break;
            }
            const double aj = s_num / norm_sq;
            diag_a[static_cast<std::size_t>(j)] = aj;

            if (j + 1 < Q_eff) {
                const double beta_sq = (j > 0) ? (norm_sq / norm_sq_prev) : 0.0;
                double norm_sq_next = 0.0;
                for (std::size_t i = 0; i < n_src; ++i) {
                    const double val = (s_pts[i] - aj) * p_curr[i] - ((j > 0) ? (beta_sq * p_prev[i]) : 0.0);
                    p_next[i] = val;
                    norm_sq_next += w_norm[i] * val * val;
                }
                if (norm_sq_next <= 1e-28 * norm_sq) {
                    Q_eff = j + 1;
                    break;
                }
                sub_b[static_cast<std::size_t>(j)] = std::sqrt(norm_sq_next / norm_sq);
                p_prev.swap(p_curr);
                p_curr.swap(p_next);
                norm_sq_prev = norm_sq;
            }
        }

        if (Q_eff <= 0) {
            return rule;
        }

        const long q_dim = static_cast<long>(Q_eff);
        Tensor<double> J(q_dim, q_dim);
        J = 0.0;
        for (int j = 0; j < Q_eff; ++j) {
            J(j, j) = diag_a[static_cast<std::size_t>(j)];
            if (j + 1 < Q_eff) {
                const double bj = sub_b[static_cast<std::size_t>(j)];
                J(j, j + 1) = bj;
                J(j + 1, j) = bj;
            }
        }

        Tensor<double> V;
        Tensor<double> evals;
        if (Q_eff == 1) {
            V = Tensor<double>(1L, 1L);
            V(0L, 0L) = 1.0;
            evals = Tensor<double>(1L);
            evals(0L) = J(0L, 0L);
        } else {
            syev(J, V, evals);
        }

        const double level_downscale = (-two_n_ln2 > -745.0) ? std::exp(-two_n_ln2) : 0.0;
        rule.nodes.reserve(static_cast<std::size_t>(Q_eff));
        for (int q = 0; q < Q_eff; ++q) {
            const double s_q = std::clamp(evals(q), 1e-300, 1.0);
            const double u_tilde = u_max * s_q;
            const double gamma_tilde = 1.0 / u_tilde;
            const double v0 = V(0L, q);
            const double wh_q = w_sum * (v0 * v0);

            GCNode node;
            node.gamma_tilde = gamma_tilde;
            node.weight_hat = wh_q;
            node.weight_scaled = level_downscale * wh_q;
            node.MT = moment_op_.make_MT(gamma_tilde);
            rule.nodes.push_back(std::move(node));
        }

        return rule;
    }

    int k_;
    double eps_;
    double mu_user_;
    double cell_width_L_;
    double mu_sim_;
    double op_scale_;
    double gamma_spiky_;
    double gamma_diff_;
    int M_;
    int window_width_W_ = 5;
    int gc_rank_Q_;
    int n_max_;
    int p_min_ = 0;
    int p_max_ = 0;
    int p_tail_ = 0;
    double delta_t_ = 0.0;
    double beta_ = 1.0;
    double w_const_sim_ = 0.0;
    double w_const_scaled_ = 0.0;
    WindowedKernelKind kind_;

    LegendreMomentOp moment_op_;
    std::vector<int> s_base_spiky_;
    std::vector<int> s_base_diff_;
    std::vector<double> gamma_ref_table_;
    std::vector<WindowedQuadratureTerm> terms_;
    std::vector<GCLevelRule> gc_rules_;

    mutable std::mutex extra_gc_mtx_;
    mutable std::vector<GCLevelRule> extra_gc_rules_;
};

/// Free-function convenience wrapper returning the k x k upper-triangular local moment matrix M(gamma).
inline Tensor<double> make_local_moment_matrix(int k, double gamma) {
    LegendreMomentOp op(k);
    return op.make_M(gamma);
}

/// Free-function convenience wrapper applying M(gamma)^{\otimes 3} to a k x k x k scaling tensor.
template <typename T>
inline Tensor<T> apply_M_3d(int k, double gamma, const Tensor<T>& s_in) {
    LegendreMomentOp op(k);
    return op.apply_M_3d(gamma, s_in);
}

} // namespace madness

#endif // MADNESS_MRA_WINDOWED_QUADRATURE_H__INCLUDED
