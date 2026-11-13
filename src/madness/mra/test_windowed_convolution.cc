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
*/

/// \file test_windowed_convolution.cc
/// \brief Independent 4-gate adversarial verification suite for Version 7
///        Tail-Completed Windowed Non-Standard (NS) 3D Convolutions.
///
/// Covers:
///   - Gate 1: Machine-precision Theorem 3.1 filter-moment commutation (< 1e-13),
///             3D Gauss-Christoffel spectral 2-norm error (< 1e-11 up to n = 30),
///             pointwise Coulomb/BSH kernel quadrature & discrete low-tail audit,
///             and deep-tree Key<3> / 64-bit hash collision audit up to Level = 35.
///   - Gate 2: Artificially unbalanced 3D coarse-to-fine jump octrees (3, 6, 10
///             level jumps across faces, edges, and corners; unit & non-unit cells)
///             verifying zero erf boundary oscillations.
///   - Gate 3: 22+ level all-electron point-nucleus Slater cusp oracle (Z = 20, 50)
///             at R_A = (0.5000004, 0.5, 0.5) next to x = 0.5, midplane dynamic
///             output refinement, BSH cross-checks, vector apply, and opt-in bridge.
///   - Gate 4: Deterministic multi-rank parity checksums (--write-parity /
///             --check-parity), MemoryMeasurer + WorldObject zero-leak check,
///             and runtime/node-count benchmark vs CoulombOperator / BSHOperator3D.

#include <madness/mra/mra.h>
#include <madness/mra/vmra.h>
#include <madness/mra/operator.h>
#include <madness/mra/memory_measurement.h>
#include <madness/mra/windowed_quadrature.h>
#include <madness/mra/windowed_convolution.h>
#if __has_include(<madness/mra/windowed_ns_halo.h>)
#include <madness/mra/windowed_ns_halo.h>
#endif
#include <madness/constants.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {

using madness::BSHOperator3D;
using madness::BSHOperatorPtr3D;
using madness::CoulombOperator;
using madness::Function;
using madness::FunctionCommonData;
using madness::FunctionDefaults;
using madness::FunctionFactory;
using madness::FunctionFunctorInterface;
using madness::FunctionNode;
using madness::GCLevelRule;
using madness::GCNode;
using madness::GenTensor;
using madness::Key;
using madness::KeyChildIterator;
using madness::LegendreMomentOp;
using madness::Level;
using madness::MemoryMeasurer;
using madness::SeparatedConvolution;
using madness::Tensor;
using madness::Translation;
using madness::TT_FULL;
using madness::Vector;
using madness::WindowedBSHOperator3D;
using madness::WindowedBSHOperatorPtr;
using madness::WindowedBSHOperatorPtr3D;
using madness::WindowedConvolutionOptions;
using madness::WindowedCoulombOperator;
using madness::WindowedCoulombOperatorPtr;
using madness::WindowedKernelKind;
using madness::WindowedQuadratureConfig;
using madness::WindowedQuadraturePlan;
using madness::WindowedQuadratureTerm;
using madness::WindowedSeparatedConvolution;
using madness::World;
using madness::WorldContainer;

using coord_3d = Vector<double, 3>;

struct SuiteOptions {
    bool small = false;
    bool run_gate1 = true;
    bool run_gate2 = true;
    bool run_gate3 = true;
    bool run_gate4 = true;
    std::string write_parity_path;
    std::string check_parity_path;
};

/// Deterministic platform-independent SplitMix64 generator & hash mixer.
class SplitMix64 {
public:
    explicit SplitMix64(std::uint64_t seed) noexcept : state_(seed) {}

    static std::uint64_t mix(std::uint64_t x) noexcept {
        x += 0x9e3779b97f4a7c15ULL;
        x = (x ^ (x >> 30U)) * 0xbf58476d1ce4e5b9ULL;
        x = (x ^ (x >> 27U)) * 0x94d049bb133111ebULL;
        return x ^ (x >> 31U);
    }

    std::uint64_t next_u64() noexcept {
        state_ += 0x9e3779b97f4a7c15ULL;
        std::uint64_t z = state_;
        z = (z ^ (z >> 30U)) * 0xbf58476d1ce4e5b9ULL;
        z = (z ^ (z >> 27U)) * 0x94d049bb133111ebULL;
        return z ^ (z >> 31U);
    }

    double next_signed_unit() noexcept {
        const std::uint64_t bits = next_u64() >> 11U;
        const double u = static_cast<double>(bits) * (1.0 / 9007199254740992.0);
        return 2.0 * u - 1.0;
    }

private:
    std::uint64_t state_ = 0;
};

struct SplitMixKeyHash {
    std::size_t operator()(const Key<3>& key) const noexcept {
        const auto& l = key.translation();
        std::uint64_t h = SplitMix64::mix(static_cast<std::uint64_t>(static_cast<std::uint32_t>(key.level())));
        h ^= SplitMix64::mix(static_cast<std::uint64_t>(l[0]) + 0x100000001b3ULL);
        h = SplitMix64::mix(h ^ static_cast<std::uint64_t>(l[1]));
        h = SplitMix64::mix(h ^ static_cast<std::uint64_t>(l[2]));
        return static_cast<std::size_t>(h);
    }
};

/// Collective MPI-safe checkpoint reporter printing 16-digit quantitative diagnostics on rank 0.
class GateReporter {
public:
    explicit GateReporter(World& world) noexcept : world_(world) {}

    void begin_gate(const std::string& title) {
        world_.gop.fence();
        gate_start_ = madness::wall_time();
        if (world_.rank() == 0) {
            std::cout << "\n========================================================================\n"
                      << title << "\n"
                      << "========================================================================\n";
        }
    }

    void end_gate(const std::string& name) {
        world_.gop.fence();
        const double dt = madness::wall_time() - gate_start_;
        if (world_.rank() == 0) {
            std::cout << "[TIME] " << name << " completed in " << std::fixed << std::setprecision(3)
                      << dt << " s\n" << std::defaultfloat;
        }
    }

    bool check_le(const std::string& label, double metric, double bound) {
        int local_fail = (!std::isfinite(metric) || !(metric <= bound)) ? 1 : 0;
        world_.gop.max(local_fail);
        ++total_checks_;
        if (local_fail != 0) {
            ++failed_checks_;
        }
        if (world_.rank() == 0) {
            std::cout << (local_fail == 0 ? "  [PASS] " : "  [FAIL] ")
                      << std::left << std::setw(54) << label
                      << " val=" << std::scientific << std::setprecision(6) << metric
                      << " <= " << std::setprecision(3) << bound
                      << std::defaultfloat << "\n";
        }
        return local_fail == 0;
    }

    bool check_true(const std::string& label, bool cond, const std::string& detail = "") {
        int local_fail = cond ? 0 : 1;
        world_.gop.max(local_fail);
        ++total_checks_;
        if (local_fail != 0) {
            ++failed_checks_;
        }
        if (world_.rank() == 0) {
            std::cout << (local_fail == 0 ? "  [PASS] " : "  [FAIL] ")
                      << std::left << std::setw(54) << label;
            if (!detail.empty()) {
                std::cout << " (" << detail << ")";
            }
            std::cout << "\n";
        }
        return local_fail == 0;
    }

    int total_checks() const noexcept { return total_checks_; }
    int failed_checks() const noexcept { return failed_checks_; }

private:
    World& world_;
    double gate_start_ = 0.0;
    int total_checks_ = 0;
    int failed_checks_ = 0;
};

/// Ledger of deterministic observables for multi-rank parity verification (np=1 vs np=2,4).
class ParityLedger {
public:
    void record(const std::string& key, double val) {
        entries_[key] = val;
    }

    const std::map<std::string, double>& entries() const noexcept { return entries_; }

    bool write_to_file(World& world, const std::string& path) const {
        int ok = 1;
        if (world.rank() == 0) {
            std::ofstream out(path, std::ios::out | std::ios::trunc);
            if (!out) {
                ok = 0;
            } else {
                out << std::setprecision(17);
                for (const auto& [k, v] : entries_) {
                    out << k << " " << v << "\n";
                }
            }
        }
        world.gop.broadcast(ok, 0);
        return ok == 1;
    }

    bool check_against_file(World& world, const std::string& path, GateReporter& rep, double tol) const {
        std::map<std::string, double> ref;
        int read_ok = 1;
        if (world.rank() == 0) {
            std::ifstream in(path);
            if (!in) {
                read_ok = 0;
            } else {
                std::string k;
                double v = 0.0;
                while (in >> k >> v) {
                    ref[k] = v;
                }
            }
        }
        world.gop.broadcast(read_ok, 0);
        if (read_ok == 0) {
            return rep.check_true("Parity reference file readable: " + path, false);
        }
        world.gop.broadcast_serializable(ref, 0);

        double max_rel_diff = 0.0;
        std::size_t matched = 0;
        for (const auto& [k, v] : entries_) {
            const auto it = ref.find(k);
            if (it == ref.end()) {
                return rep.check_true("Parity key present in reference: " + k, false);
            }
            const double denom = std::max(1.0, std::abs(it->second));
            const double rel = std::abs(v - it->second) / denom;
            max_rel_diff = std::max(max_rel_diff, rel);
            if (world.rank() == 0) {
                std::cout << "    parity[" << k << "] cur=" << std::setprecision(16) << v
                          << " ref=" << it->second << " rel_diff=" << std::scientific
                          << std::setprecision(4) << rel << std::defaultfloat << "\n";
            }
            ++matched;
        }
        rep.check_true("All parity observables matched", matched == entries_.size(),
                       "count=" + std::to_string(matched));
        return rep.check_le("Max multi-rank parity discrepancy vs reference", max_rel_diff, tol);
    }

private:
    std::map<std::string, double> entries_;
};

/// RAII saver/restorer for FunctionDefaults<3> state across sub-tests.
class ScopedMraDefaults {
public:
    explicit ScopedMraDefaults(World& world)
        : world_(world)
        , k_(FunctionDefaults<3>::get_k())
        , thresh_(FunctionDefaults<3>::get_thresh())
        , initial_level_(FunctionDefaults<3>::get_initial_level())
        , special_level_(FunctionDefaults<3>::get_special_level())
        , max_refine_level_(FunctionDefaults<3>::get_max_refine_level())
        , truncate_mode_(FunctionDefaults<3>::get_truncate_mode())
        , refine_(FunctionDefaults<3>::get_refine())
        , autorefine_(FunctionDefaults<3>::get_autorefine())
        , truncate_on_project_(FunctionDefaults<3>::get_truncate_on_project())
        , cell_(madness::copy(FunctionDefaults<3>::get_cell())) {}

    ~ScopedMraDefaults() {
        world_.gop.fence();
        FunctionDefaults<3>::set_cell(cell_);
        FunctionDefaults<3>::set_k(k_);
        FunctionDefaults<3>::set_thresh(thresh_);
        FunctionDefaults<3>::set_max_refine_level(30);
        FunctionDefaults<3>::set_initial_level(initial_level_);
        FunctionDefaults<3>::set_special_level(special_level_);
        FunctionDefaults<3>::set_max_refine_level(max_refine_level_);
        FunctionDefaults<3>::set_truncate_mode(truncate_mode_);
        FunctionDefaults<3>::set_refine(refine_);
        FunctionDefaults<3>::set_autorefine(autorefine_);
        FunctionDefaults<3>::set_truncate_on_project(truncate_on_project_);
        world_.gop.fence();
    }

    ScopedMraDefaults(const ScopedMraDefaults&) = delete;
    ScopedMraDefaults& operator=(const ScopedMraDefaults&) = delete;

private:
    World& world_;
    int k_;
    double thresh_;
    int initial_level_;
    int special_level_;
    int max_refine_level_;
    int truncate_mode_;
    bool refine_;
    bool autorefine_;
    bool truncate_on_project_;
    Tensor<double> cell_;
};

// ============================================================================
// Exact Analytical Oracles (Gaussian, Coulomb-Gaussian, Slater Cusp, BSH)
// ============================================================================

/// Normalized 3D Gaussian f(r) = (alpha / pi)^{3/2} * exp(-alpha * |r - center|^2).
class NormalizedGaussianFunctor final : public FunctionFunctorInterface<double, 3> {
public:
    NormalizedGaussianFunctor(const coord_3d& center, double alpha,
                              std::vector<coord_3d> sp = {}, int sp_lvl = 0)
        : center_(center)
        , alpha_(alpha)
        , norm_pref_(std::pow(alpha / madness::constants::pi, 1.5))
        , special_pts_(std::move(sp))
        , special_lvl_(sp_lvl) {}

    double operator()(const coord_3d& x) const override {
        const double dx = x[0] - center_[0];
        const double dy = x[1] - center_[1];
        const double dz = x[2] - center_[2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double arg = alpha_ * r2;
        if (arg > 700.0) return 0.0;
        return norm_pref_ * std::exp(-arg);
    }

    std::vector<coord_3d> special_points() const override { return special_pts_; }
    Level special_level() const override { return special_lvl_; }

private:
    coord_3d center_;
    double alpha_ = 1.0;
    double norm_pref_ = 1.0;
    std::vector<coord_3d> special_pts_;
    int special_lvl_ = 0;
};

/// Exact Coulomb potential of NormalizedGaussianFunctor: V(r) = erf(sqrt(alpha) * r) / r.
class ExactGaussianCoulombFunctor final : public FunctionFunctorInterface<double, 3> {
public:
    ExactGaussianCoulombFunctor(const coord_3d& center, double alpha,
                                std::vector<coord_3d> sp = {}, int sp_lvl = 0)
        : center_(center)
        , sqrt_alpha_(std::sqrt(alpha))
        , special_pts_(std::move(sp))
        , special_lvl_(sp_lvl) {}

    double operator()(const coord_3d& x) const override {
        const double dx = x[0] - center_[0];
        const double dy = x[1] - center_[1];
        const double dz = x[2] - center_[2];
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double u = sqrt_alpha_ * r;
        if (u < 1e-4) {
            const double u2 = u * u;
            const double series = 1.0 - (u2 / 3.0) + (u2 * u2 / 10.0) - (u2 * u2 * u2 / 42.0);
            return (2.0 * sqrt_alpha_ * madness::constants::inv_sqrt_pi) * series;
        }
        return std::erf(u) / r;
    }

    std::vector<coord_3d> special_points() const override { return special_pts_; }
    Level special_level() const override { return special_lvl_; }

private:
    coord_3d center_;
    double sqrt_alpha_ = 1.0;
    std::vector<coord_3d> special_pts_;
    int special_lvl_ = 0;
};

/// Overflow-free scaled complementary error function erfcx(x) = exp(x^2) * erfc(x).
inline double erfcx_safe(double x) noexcept {
    if (x < 0.0) {
        if (x < -26.0) return 0.0; // Guard caller handles negative x via direct exp(-mu*r)*erfc(x) when needed
        return std::exp(x * x) * std::erfc(x);
    }
    if (x < 25.0) {
        return std::exp(x * x) * std::erfc(x);
    }
    // Asymptotic continued-fraction / inverse-power expansion for x >= 25 (accurate to < 1e-16).
    const double inv_x = 1.0 / x;
    const double inv_x2 = inv_x * inv_x;
    double sum = 1.0;
    double term = 1.0;
    for (int m = 1; m <= 12; ++m) {
        term *= -0.5 * static_cast<double>(2 * m - 1) * inv_x2;
        sum += term;
    }
    return madness::constants::inv_sqrt_pi * inv_x * sum;
}

/// Exact BSH convolution (exp(-mu*r)/(4*pi*r)) * NormalizedGaussianFunctor(center, alpha).
class ExactGaussianBSHFunctor final : public FunctionFunctorInterface<double, 3> {
public:
    ExactGaussianBSHFunctor(const coord_3d& center, double alpha, double mu)
        : center_(center)
        , alpha_(alpha)
        , mu_(mu)
        , sqrt_alpha_(std::sqrt(alpha)) {}

    double operator()(const coord_3d& x) const override {
        const double dx = x[0] - center_[0];
        const double dy = x[1] - center_[1];
        const double dz = x[2] - center_[2];
        const double r2 = dx * dx + dy * dy + dz * dz;
        const double r = std::sqrt(r2);
        const double c = 0.5 * mu_ / sqrt_alpha_;
        const double u0 = (1.0 / (4.0 * madness::constants::pi)) *
                          (2.0 * sqrt_alpha_ * madness::constants::inv_sqrt_pi - mu_ * erfcx_safe(c));
        if (r < 1e-5) {
            const double f0 = std::pow(alpha_ / madness::constants::pi, 1.5);
            const double c2 = (mu_ * mu_ * u0 - f0) / 6.0;
            const double c4 = (mu_ * mu_ * mu_ * mu_ * u0 - mu_ * mu_ * f0 + 6.0 * alpha_ * f0) / 120.0;
            return u0 + c2 * r2 + c4 * r2 * r2;
        }

        const double sr = sqrt_alpha_ * r;
        const double arg_minus = c - sr;
        const double arg_plus = c + sr;

        double term_minus = 0.0;
        if (arg_minus >= 0.0) {
            term_minus = (alpha_ * r2 < 700.0) ? (std::exp(-alpha_ * r2) * erfcx_safe(arg_minus)) : 0.0;
        } else {
            // c - sr < 0: use exp(c^2 - mu*r) * erfc(c - sr), where c^2 - mu*r = (c - sr)^2 - alpha*r^2 < c^2.
            const double log_pref = c * c - mu_ * r;
            if (log_pref > -740.0) {
                term_minus = std::exp(log_pref) * std::erfc(arg_minus);
            }
        }

        const double term_plus = (alpha_ * r2 < 700.0) ? (std::exp(-alpha_ * r2) * erfcx_safe(arg_plus)) : 0.0;
        return (term_minus - term_plus) / (8.0 * madness::constants::pi * r);
    }

private:
    coord_3d center_;
    double alpha_ = 1.0;
    double mu_ = 1.0;
    double sqrt_alpha_ = 1.0;
};

/// Normalized point-nucleus 1s Slater density: rho(r) = (Z^3 / pi) * exp(-2 * Z * |r - R_A|).
class SlaterCuspDensityFunctor final : public FunctionFunctorInterface<double, 3> {
public:
    SlaterCuspDensityFunctor(const coord_3d& center, double Z, int min_cusp_level = 11)
        : center_(center)
        , Z_(Z)
        , pref_((Z * Z * Z) / madness::constants::pi)
        , min_cusp_level_(min_cusp_level) {}

    double operator()(const coord_3d& x) const override {
        const double dx = x[0] - center_[0];
        const double dy = x[1] - center_[1];
        const double dz = x[2] - center_[2];
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double expo = 2.0 * Z_ * r;
        if (expo > 700.0) return 0.0;
        return pref_ * std::exp(-expo);
    }

    std::vector<coord_3d> special_points() const override { return {center_}; }
    Level special_level() const override { return min_cusp_level_; }

private:
    coord_3d center_;
    double Z_ = 1.0;
    double pref_ = 1.0;
    int min_cusp_level_ = 11;
};

/// Exact Coulomb potential of SlaterCuspDensityFunctor:
///   V(r) = (1 - (1 + Z*r) * exp(-2*Z*r)) / r,   V(0) = Z.
class ExactSlaterCoulombFunctor final : public FunctionFunctorInterface<double, 3> {
public:
    ExactSlaterCoulombFunctor(const coord_3d& center, double Z, int min_cusp_level = 11)
        : center_(center)
        , Z_(Z)
        , min_cusp_level_(min_cusp_level) {}

    double operator()(const coord_3d& x) const override {
        const double dx = x[0] - center_[0];
        const double dy = x[1] - center_[1];
        const double dz = x[2] - center_[2];
        const double r = std::sqrt(dx * dx + dy * dy + dz * dz);
        const double u = 2.0 * Z_ * r;
        if (u < 1e-3) {
            const double u2 = u * u;
            const double u3 = u2 * u;
            const double u4 = u2 * u2;
            const double u5 = u4 * u;
            const double u6 = u3 * u3;
            return Z_ * (1.0 - (u2 / 6.0) + (u3 / 12.0) - (u4 / 40.0) + (u5 / 180.0) - (u6 / 1008.0));
        }
        if (u > 700.0) {
            return 1.0 / r;
        }
        return (1.0 - (1.0 + Z_ * r) * std::exp(-u)) / r;
    }

    std::vector<coord_3d> special_points() const override { return {center_}; }
    Level special_level() const override { return min_cusp_level_; }

private:
    coord_3d center_;
    double Z_ = 1.0;
    int min_cusp_level_ = 11;
};

/// Iterative, deterministic, MPI-safe custom octree function builder (no recursion).
Function<double, 3> build_custom_octree_function(
    World& world,
    const std::shared_ptr<FunctionFunctorInterface<double, 3>>& functor,
    int k,
    double eps,
    const std::function<bool(const Key<3>&)>& should_split) {

    Function<double, 3> f = FunctionFactory<double, 3>(world)
                                .functor(functor)
                                .k(k)
                                .thresh(eps)
                                .empty();
    auto impl = f.get_impl();
    auto& dc = impl->get_coeffs();

    std::vector<Key<3>> stack;
    stack.reserve(4096);
    stack.emplace_back(0);

    while (!stack.empty()) {
        const Key<3> key = stack.back();
        stack.pop_back();

        if (should_split(key)) {
            if (dc.is_local(key)) {
                dc.replace(key, FunctionNode<double, 3>(GenTensor<double>(), true));
            }
            for (KeyChildIterator<3> kit(key); kit; ++kit) {
                stack.push_back(kit.key());
            }
        } else {
            if (dc.is_local(key)) {
                const Tensor<double> s = impl->project(key);
                dc.replace(key, FunctionNode<double, 3>(GenTensor<double>(s, eps, TT_FULL), false));
            }
        }
    }

    dc.process_pending();
    world.gop.fence();
    f.verify_tree();
    return f;
}

// ============================================================================
// Gate 1: Spectral, Quadrature, Commutation & Deep-Key Identities
// ============================================================================

Tensor<double> matmul_2d(const Tensor<double>& A, const Tensor<double>& B) {
    const long m = A.dim(0);
    const long n = A.dim(1);
    const long p = B.dim(1);
    MADNESS_CHECK(B.dim(0) == n);
    Tensor<double> C(m, p);
    C = 0.0;
    for (long i = 0; i < m; ++i) {
        for (long r = 0; r < n; ++r) {
            const double air = A(i, r);
            for (long j = 0; j < p; ++j) {
                C(i, j) += air * B(r, j);
            }
        }
    }
    return C;
}

Tensor<double> transpose_2d(const Tensor<double>& A) {
    const long m = A.dim(0);
    const long n = A.dim(1);
    Tensor<double> AT(n, m);
    for (long i = 0; i < m; ++i) {
        for (long j = 0; j < n; ++j) {
            AT(j, i) = A(i, j);
        }
    }
    return AT;
}

Tensor<double> prolong_3d_steps(const Tensor<double>& s_in,
                                const Tensor<double>& h0,
                                const Tensor<double>& h1,
                                int delta_j,
                                Translation dlx,
                                Translation dly,
                                Translation dlz) {
    Tensor<double> cur = madness::copy(s_in);
    const Tensor<double> h[2] = {h0, h1};
    for (int step = 0; step < delta_j; ++step) {
        const int shift = delta_j - 1 - step;
        const int bx = static_cast<int>((dlx >> shift) & 1LL);
        const int by = static_cast<int>((dly >> shift) & 1LL);
        const int bz = static_cast<int>((dlz >> shift) & 1LL);
        Tensor<double> mats[3] = {h[bx], h[by], h[bz]};
        cur = madness::general_transform(cur, mats);
    }
    return cur;
}

Tensor<double> restrict_3d_steps(const Tensor<double>& s_in,
                                 const Tensor<double>& h0T,
                                 const Tensor<double>& h1T,
                                 int delta_j,
                                 Translation dlx,
                                 Translation dly,
                                 Translation dlz) {
    Tensor<double> cur = madness::copy(s_in);
    const Tensor<double> hT[2] = {h0T, h1T};
    for (int step = delta_j - 1; step >= 0; --step) {
        const int shift = delta_j - 1 - step;
        const int bx = static_cast<int>((dlx >> shift) & 1LL);
        const int by = static_cast<int>((dly >> shift) & 1LL);
        const int bz = static_cast<int>((dlz >> shift) & 1LL);
        Tensor<double> mats[3] = {hT[bx], hT[by], hT[bz]};
        cur = madness::general_transform(cur, mats);
    }
    return cur;
}

void run_gate1(World& world, const SuiteOptions& opts, GateReporter& rep) {
    rep.begin_gate("GATE 1: Theorem 3.1 Commutation, 3D GC 2-Norm, Tail Quadrature & Level-35 Keys");

    // ------------------------------------------------------------------------
    // Gate 1a: Two-threshold scale separation & Theorem 3.1 exact commutation
    // ------------------------------------------------------------------------
    const std::vector<int> k_list = opts.small ? std::vector<int>{8, 10} : std::vector<int>{6, 8, 10, 12};
    const std::vector<int> jumps = {1, 5, 15};

    double max_comm_1d_err = 0.0;
    double max_comm_3d_rel = 0.0;
    double max_pou_3d_rel = 0.0;
    double max_gamma_ref_rel = 0.0;
    double max_mt_transpose_err = 0.0;
    bool thresholds_valid = true;

    for (const int k : k_list) {
        WindowedQuadratureConfig cfg;
        cfg.k = k;
        cfg.eps = 1e-9;
        cfg.max_tree_level = 30;
        cfg.kind = WindowedKernelKind::Coulomb;
        const WindowedQuadraturePlan plan(cfg);

        const double expected_gs = static_cast<double>(k * k * (k + 1) * (k + 1)) / 30.0;
        if (std::abs(plan.gamma_spiky() - expected_gs) > 1e-12 * expected_gs ||
            !(plan.gamma_diff() > 0.01 && plan.gamma_diff() < 5.0) ||
            !(plan.gamma_spiky() / plan.gamma_diff() > 30.0) ||
            !(plan.window_width_W() >= 5 && plan.window_width_W() <= 8) ||
            plan.p_tail() != plan.p_max() + 10 * plan.M()) {
            thresholds_valid = false;
        }

        // Verify universal (w_eff, m) stiffness table consistency across all active dyadic terms.
        const double two_ln2 = 2.0 * std::log(2.0);
        for (const WindowedQuadratureTerm& term : plan.terms()) {
            if (!term.on_dyadic_grid || term.n_diff < 0) continue;
            for (int j = term.j0; j <= term.j1; ++j) {
                const int w_eff = j - (term.n_spiky + 1);
                if (w_eff < 0 || w_eff >= plan.window_width_W()) {
                    thresholds_valid = false;
                    continue;
                }
                const double gamma_pj = std::exp(2.0 * term.t_p - static_cast<double>(j) * two_ln2);
                const double gref = plan.gamma_ref(w_eff, term.sub_octave_m);
                const double rel = std::abs(gamma_pj - gref) / gref;
                max_gamma_ref_rel = std::max(max_gamma_ref_rel, rel);
            }
        }

        const LegendreMomentOp mop(k);
        const auto& cdata = FunctionCommonData<double, 3>::get(k);
        const Tensor<double> H_mat[2] = {transpose_2d(cdata.h0), transpose_2d(cdata.h1)};

        // Deterministic asymmetric k^3 scaling tensor to expose any M vs MT transpose error.
        SplitMix64 rng(0x12345678ULL + static_cast<std::uint64_t>(k) * 977ULL);
        Tensor<double> s_coarse(static_cast<long>(k), static_cast<long>(k), static_cast<long>(k));
        for (long ix = 0; ix < k; ++ix) {
            for (long iy = 0; iy < k; ++iy) {
                for (long iz = 0; iz < k; ++iz) {
                    s_coarse(ix, iy, iz) = rng.next_signed_unit() *
                        (1.0 + 0.3 * static_cast<double>(ix) + 0.7 * static_cast<double>(iy * iy));
                }
            }
        }

        for (const double g_mult : {1.0, 4.0}) {
            const double gamma_nB = g_mult * plan.gamma_spiky();
            const Tensor<double> M_nB = mop.make_M(gamma_nB);
            const Tensor<double> MT_nB = mop.make_MT(gamma_nB);
            max_mt_transpose_err = std::max(max_mt_transpose_err, (transpose_2d(M_nB) - MT_nB).normf());

            // Verify strict upper-triangular parity sparsity of M(gamma).
            for (int row = 0; row < k; ++row) {
                for (int col = 0; col < k; ++col) {
                    if (col < row || ((col - row) % 2 != 0)) {
                        max_mt_transpose_err = std::max(max_mt_transpose_err, std::abs(M_nB(row, col)));
                    }
                }
            }

            const Tensor<double> Ms_coarse = mop.apply_M_3d(gamma_nB, s_coarse);

            for (const int dj : jumps) {
                const double gamma_jstar = std::ldexp(gamma_nB, -2 * dj);
                const Tensor<double> M_jstar = mop.make_M(gamma_jstar);
                const Translation span = (Translation(1) << dj);
                const std::array<Translation, 4> sample_dl = {
                    Translation(0),
                    span - 1,
                    span / 3,
                    (2 * span) / 3
                };

                for (const Translation dl : sample_dl) {
                    Tensor<double> H_path(static_cast<long>(k), static_cast<long>(k));
                    H_path = 0.0;
                    for (int i = 0; i < k; ++i) H_path(i, i) = 1.0;
                    for (int step = 0; step < dj; ++step) {
                        const int shift = dj - 1 - step;
                        const int bit = static_cast<int>((dl >> shift) & 1LL);
                        H_path = matmul_2d(H_mat[bit], H_path);
                    }
                    const Tensor<double> lhs_1d = matmul_2d(M_jstar, H_path);
                    const Tensor<double> rhs_1d = matmul_2d(H_path, M_nB);
                    max_comm_1d_err = std::max(max_comm_1d_err, (lhs_1d - rhs_1d).normf());
                }

                // 3D tensor commutation along corner and interior sub-box paths.
                const std::array<std::array<Translation, 3>, 3> paths_3d = {{
                    {Translation(0), Translation(0), Translation(0)},
                    {span - 1, span - 1, span - 1},
                    {span / 3, (2 * span) / 3, span / 2}
                }};
                for (const auto& p3 : paths_3d) {
                    const Tensor<double> s_V = prolong_3d_steps(
                        s_coarse, cdata.h0, cdata.h1, dj, p3[0], p3[1], p3[2]);
                    const Tensor<double> path_a = mop.apply_M_3d(gamma_jstar, s_V);
                    const Tensor<double> path_b = prolong_3d_steps(
                        Ms_coarse, cdata.h0, cdata.h1, dj, p3[0], p3[1], p3[2]);
                    const double rel = (path_a - path_b).normf() / std::max(1e-30, path_b.normf());
                    max_comm_3d_rel = std::max(max_comm_3d_rel, rel);
                }
            }

            // Partition-of-unity sum over all 8^{dj} virtual sub-boxes for dj = 1, 2.
            const int max_pou_dj = opts.small ? 1 : 2;
            for (int dj = 1; dj <= max_pou_dj; ++dj) {
                const double gamma_jstar = std::ldexp(gamma_nB, -2 * dj);
                const Translation span = (Translation(1) << dj);
                Tensor<double> accum(static_cast<long>(k), static_cast<long>(k), static_cast<long>(k));
                accum = 0.0;
                for (Translation ix = 0; ix < span; ++ix) {
                    for (Translation iy = 0; iy < span; ++iy) {
                        for (Translation iz = 0; iz < span; ++iz) {
                            const Tensor<double> s_V = prolong_3d_steps(
                                s_coarse, cdata.h0, cdata.h1, dj, ix, iy, iz);
                            const Tensor<double> Ms_V = mop.apply_M_3d(gamma_jstar, s_V);
                            accum += restrict_3d_steps(
                                Ms_V, cdata.h0T, cdata.h1T, dj, ix, iy, iz);
                        }
                    }
                }
                const double rel_pou = (accum - Ms_coarse).normf() / std::max(1e-30, Ms_coarse.normf());
                max_pou_3d_rel = std::max(max_pou_3d_rel, rel_pou);
            }
        }
    }

    rep.check_true("Invariant 2 & 4: two-threshold separation & p_tail = p_max + 10*M", thresholds_valid);
    rep.check_le("Universal stencil gamma_ref(w_eff, m) max relative drift", max_gamma_ref_rel, 1e-13);
    rep.check_le("LegendreMomentOp MT == transpose(M) & parity sparsity", max_mt_transpose_err, 1e-15);
    rep.check_le("Theorem 3.1 1D commutation ||M(g_j*)H - H M(g_nB)||_F (dj=1,5,15)", max_comm_1d_err, 1e-13);
    rep.check_le("Theorem 3.1 3D tensor commutation rel error (dj=1,5,15)", max_comm_3d_rel, 1e-13);
    rep.check_le("Theorem 3.1 3D partition-of-unity round-trip rel error", max_pou_3d_rel, 1e-13);

    // ------------------------------------------------------------------------
    // Gate 1b: 3D Gauss-Christoffel spectral 2-norm error (< 1e-11 up to n=30)
    // ------------------------------------------------------------------------
    const std::vector<int> gc_levels = opts.small
        ? std::vector<int>{0, 10, 22, 25, 30}
        : std::vector<int>{0, 5, 10, 15, 20, 22, 25, 30};

    double max_gc_sigma_coul = 0.0;
    double max_gc_sigma_bsh = 0.0;
    double max_gc_helper_rel = 0.0;
    bool gc_finite_ok = true;

    for (const int k : k_list) {
        for (int mode = 0; mode < 2; ++mode) {
            WindowedQuadratureConfig cfg;
            cfg.k = k;
            cfg.eps = 1e-9;
            cfg.max_tree_level = 30;
            if (mode == 0) {
                cfg.kind = WindowedKernelKind::Coulomb;
                cfg.cell_width_L = 1.0;
            } else {
                cfg.kind = WindowedKernelKind::BSH;
                cfg.mu_user = 1.0;
                cfg.cell_width_L = 1.0;
            }
            const WindowedQuadraturePlan plan(cfg);
            const LegendreMomentOp& mop = plan.moment_op();
            const long kk = static_cast<long>(k);
            const double two_ln2 = 2.0 * std::log(2.0);

            for (const int lvl : gc_levels) {
                const GCLevelRule& rule = plan.gc_rule(lvl);
                for (const GCNode& node : rule.nodes) {
                    if (!std::isfinite(node.gamma_tilde) || !std::isfinite(node.weight_hat) ||
                        !std::isfinite(node.weight_scaled) || node.gamma_tilde < 0.99 * plan.gamma_spiky() ||
                        node.weight_hat < 0.0) {
                        gc_finite_ok = false;
                    }
                }

                // Collect uncompressed normalized spiky terms (n_spiky >= lvl) + independent 40-point high-tail continuation.
                struct TermPair {
                    double wh = 0.0;
                    Tensor<double> MT;
                    Tensor<double> M;
                };
                std::vector<TermPair> exact_spiky;
                exact_spiky.reserve(plan.terms().size() + 40);
                const double two_lvl_ln2 = static_cast<double>(lvl) * two_ln2;

                for (const WindowedQuadratureTerm& term : plan.terms()) {
                    if (term.n_spiky < lvl) continue;
                    const double gamma_pn = std::exp(2.0 * term.t_p - two_lvl_ln2);
                    const double u_pn = 1.0 / gamma_pn;
                    double mu_fac = 1.0;
                    if (plan.mu_sim() > 0.0) {
                        const double mu_arg = (-2.0 * term.t_p < 700.0)
                            ? (0.25 * plan.mu_sim() * plan.mu_sim() * std::exp(-2.0 * term.t_p))
                            : 800.0;
                        mu_fac = (mu_arg < 700.0) ? std::exp(-mu_arg) : 0.0;
                    }
                    const double wh = term.on_dyadic_grid
                        ? (plan.op_scale() * 2.0 * madness::constants::pi * plan.delta_t() * u_pn * mu_fac)
                        : (std::exp(two_lvl_ln2) * term.weight_scaled);
                    if (wh > 0.0) {
                        exact_spiky.push_back(TermPair{wh, mop.make_MT(gamma_pn), mop.make_M(gamma_pn)});
                    }
                }

                // Independent 40-term explicit dyadic continuation above p_tail for p > p_tail remainder verification.
                if (plan.mu_sim() == 0.0) {
                    for (int extra = 1; extra <= 40; ++extra) {
                        const int p_ext = plan.p_tail() + extra;
                        const double tp_ext = static_cast<double>(p_ext) * plan.delta_t();
                        const int n_sp_ext = static_cast<int>(
                            std::floor((2.0 * tp_ext - std::log(plan.gamma_spiky())) / two_ln2));
                        if (n_sp_ext < lvl) continue;
                        const double gamma_ext = std::exp(2.0 * tp_ext - two_lvl_ln2);
                        const double wh_ext = plan.op_scale() *
                            (2.0 * madness::constants::pi * plan.delta_t()) / gamma_ext;
                        exact_spiky.push_back(TermPair{wh_ext, mop.make_MT(gamma_ext), mop.make_M(gamma_ext)});
                    }
                }

                std::vector<Tensor<double>> gc_M;
                gc_M.reserve(rule.nodes.size());
                for (const GCNode& node : rule.nodes) {
                    gc_M.push_back(mop.make_M(node.gamma_tilde));
                }

                Tensor<double> tmp(kk, kk, kk);
                Tensor<double> w1(kk, kk, kk);
                Tensor<double> w2(kk, kk, kk);

                auto apply_err_op = [&](const Tensor<double>& x, bool adjoint) -> Tensor<double> {
                    Tensor<double> y(kk, kk, kk);
                    y = 0.0;
                    for (const auto& st : exact_spiky) {
                        const double* mat = adjoint ? st.M.ptr() : st.MT.ptr();
                        mop.apply_MT_buffer_3d(mat, x, tmp, w1, w2);
                        y.gaxpy(1.0, tmp, st.wh);
                    }
                    for (std::size_t q = 0; q < rule.nodes.size(); ++q) {
                        const double* mat = adjoint ? gc_M[q].ptr() : rule.nodes[q].MT.ptr();
                        mop.apply_MT_buffer_3d(mat, x, tmp, w1, w2);
                        y.gaxpy(1.0, tmp, -rule.nodes[q].weight_hat);
                    }
                    return y;
                };

                // 25-step singular-value power iteration on E_hat^T * E_hat for genuine spectral 2-norm.
                SplitMix64 rng(0xabcdef01ULL + static_cast<std::uint64_t>(k * 131 + lvl * 17 + mode));
                Tensor<double> v(kk, kk, kk);
                for (long i = 0; i < kk; ++i) {
                    for (long j = 0; j < kk; ++j) {
                        for (long l = 0; l < kk; ++l) {
                            v(i, j, l) = rng.next_signed_unit();
                        }
                    }
                }
                v.scale(1.0 / v.normf());

                double sigma_est = 0.0;
                const int n_power = opts.small ? 18 : 25;
                for (int it = 0; it < n_power; ++it) {
                    const Tensor<double> ev = apply_err_op(v, false);
                    sigma_est = ev.normf();
                    if (sigma_est <= 1e-30) {
                        sigma_est = 0.0;
                        break;
                    }
                    Tensor<double> et_ev = apply_err_op(ev, true);
                    const double nrm = et_ev.normf();
                    if (nrm <= 1e-30) break;
                    et_ev.scale(1.0 / nrm);
                    v = et_ev;
                }

                if (mode == 0) {
                    max_gc_sigma_coul = std::max(max_gc_sigma_coul, sigma_est);
                } else {
                    max_gc_sigma_bsh = std::max(max_gc_sigma_bsh, sigma_est);
                }

                // Cross-check Stage 2 accumulate_gc_leaf_3d helper at representative levels.
                if (lvl <= 22) {
                    Tensor<double> s_acc(kk, kk, kk);
                    s_acc = 0.0;
                    plan.accumulate_gc_leaf_3d(lvl, v, s_acc, tmp, w1, w2);

                    Tensor<double> s_ref(kk, kk, kk);
                    s_ref = 0.0;
                    const double scale_down = std::exp(-two_lvl_ln2);
                    for (const GCNode& node : rule.nodes) {
                        mop.apply_MT_buffer_3d(node.MT.ptr(), v, tmp, w1, w2);
                        s_ref.gaxpy(1.0, tmp, scale_down * node.weight_hat);
                    }
                    const double helper_rel = (s_acc - s_ref).normf() / std::max(1e-30, s_ref.normf());
                    max_gc_helper_rel = std::max(max_gc_helper_rel, helper_rel);
                }
            }
        }
    }

    rep.check_true("GC nodes finite, positive weights, gamma_tilde >= gamma_spiky up to n=30", gc_finite_ok);
    rep.check_le("3D GC normalized spectral 2-norm ||E_hat_Q||_2 (Coulomb, n<=30)", max_gc_sigma_coul, 1e-11);
    rep.check_le("3D GC normalized spectral 2-norm ||E_hat_Q||_2 (BSH, n<=30)", max_gc_sigma_bsh, 1e-11);
    rep.check_le("Stage 2 accumulate_gc_leaf_3d consistency rel error", max_gc_helper_rel, 1e-14);

    // ------------------------------------------------------------------------
    // Gate 1c: Pointwise kernel quadrature & discrete low-tail audit
    // ------------------------------------------------------------------------
    double max_coul_unit_rel = 0.0;
    double max_coul_outer_rel = 0.0;
    double max_coul_cell8_rel = 0.0;
    double max_bsh_scaled_err = 0.0;

    {
        WindowedQuadratureConfig cfg_u;
        cfg_u.k = 10;
        cfg_u.eps = 1e-9;
        cfg_u.cell_width_L = 1.0;
        cfg_u.max_tree_level = 25;
        cfg_u.kind = WindowedKernelKind::Coulomb;
        const WindowedQuadraturePlan plan_u(cfg_u);

        WindowedQuadratureConfig cfg_8 = cfg_u;
        cfg_8.cell_width_L = 8.0;
        const WindowedQuadraturePlan plan_8(cfg_8);

        WindowedQuadratureConfig cfg_bsh1;
        cfg_bsh1.k = 10;
        cfg_bsh1.eps = 1e-9;
        cfg_bsh1.mu_user = 0.5;
        cfg_bsh1.cell_width_L = 1.0;
        cfg_bsh1.max_tree_level = 25;
        cfg_bsh1.kind = WindowedKernelKind::BSH;
        const WindowedQuadraturePlan plan_bsh1(cfg_bsh1);

        WindowedQuadratureConfig cfg_bsh2 = cfg_bsh1;
        cfg_bsh2.mu_user = 2.0;
        cfg_bsh2.cell_width_L = 6.0;
        const WindowedQuadraturePlan plan_bsh2(cfg_bsh2);

        constexpr int n_pts = 500;
        const double log_r_min = std::log(1e-8);
        const double log_r_max = std::log(std::sqrt(3.0));
        for (int i = 0; i < n_pts; ++i) {
            const double frac = static_cast<double>(i) / static_cast<double>(n_pts - 1);
            const double r_sim = std::exp(log_r_min + frac * (log_r_max - log_r_min));

            const double k_exact_sim = 1.0 / r_sim;
            const double rel_u = std::abs(plan_u.eval_kernel_sim(r_sim) - k_exact_sim) / k_exact_sim;
            max_coul_unit_rel = std::max(max_coul_unit_rel, rel_u);
            if (r_sim >= 0.5) {
                max_coul_outer_rel = std::max(max_coul_outer_rel, rel_u);
            }

            const double r_user8 = 8.0 * r_sim;
            const double k_exact_8 = 1.0 / r_user8;
            const double rel_8 = std::abs(plan_8.eval_kernel_user(r_user8) - k_exact_8) / k_exact_8;
            max_coul_cell8_rel = std::max(max_coul_cell8_rel, rel_8);

            const double k_bsh1_ex = std::exp(-0.5 * r_sim) / (4.0 * madness::constants::pi * r_sim);
            const double err_bsh1 = std::abs(plan_bsh1.eval_kernel_user(r_sim) - k_bsh1_ex) *
                                    (4.0 * madness::constants::pi * r_sim);
            max_bsh_scaled_err = std::max(max_bsh_scaled_err, err_bsh1);

            const double r_user6 = 6.0 * r_sim;
            const double k_bsh2_ex = std::exp(-2.0 * r_user6) / (4.0 * madness::constants::pi * r_user6);
            const double err_bsh2 = std::abs(plan_bsh2.eval_kernel_user(r_user6) - k_bsh2_ex) *
                                    (4.0 * madness::constants::pi * r_user6);
            max_bsh_scaled_err = std::max(max_bsh_scaled_err, err_bsh2);
        }
    }

    rep.check_le("Coulomb kernel max rel error on r in [1e-8, sqrt(3)] (eps=1e-9)", max_coul_unit_rel, 5e-9);
    rep.check_le("Discrete low-tail outer window r in [0.5, sqrt(3)] rel error (no EM bump)", max_coul_outer_rel, 3e-9);
    rep.check_le("Coulomb non-unit cell L=8 kernel scaling max rel error", max_coul_cell8_rel, 5e-9);
    rep.check_le("BSH kernel (mu=0.5 L=1 & mu=2.0 L=6) max 4*pi*r scaled error", max_bsh_scaled_err, 5e-9);

    // ------------------------------------------------------------------------
    // Gate 1d: Deep-tree Key<3> & 64-bit hash collision audit to Level = 35
    // ------------------------------------------------------------------------
    bool deep_key_ok = true;
    std::size_t total_deep_keys = 0;
    std::unordered_map<Key<3>, std::size_t, madness::Hash<Key<3>>> mad_map;
    std::unordered_map<Key<3>, std::size_t, SplitMixKeyHash> sm_map;
    std::unordered_set<std::size_t> sm_hashes;
    WorldContainer<Key<3>, std::uint64_t> wc(world);

    for (const Level lvl : {21, 22, 25, 30, 35}) {
        const Translation mid = (Translation(1) << (lvl - 1));
        const Translation top = (Translation(1) << lvl) - 1;

        for (Translation dx = -2; dx <= 2; ++dx) {
            for (Translation dy = -2; dy <= 2; ++dy) {
                for (Translation dz = -2; dz <= 2; ++dz) {
                    const Vector<Translation, 3> l{mid + dx, mid + dy, mid + dz};
                    const Key<3> key(lvl, l);
                    if (key.level() != lvl || key.translation()[0] != l[0] ||
                        key.translation()[1] != l[1] || key.translation()[2] != l[2]) {
                        deep_key_ok = false;
                    }
                    const Key<3> par = key.parent();
                    if (par.level() != lvl - 1 || par.translation()[0] != (l[0] >> 1)) {
                        deep_key_ok = false;
                    }
                    int found_self_in_parent = 0;
                    for (KeyChildIterator<3> kit(par); kit; ++kit) {
                        if (kit.key() == key) ++found_self_in_parent;
                    }
                    if (found_self_in_parent != 1) deep_key_ok = false;

                    const std::size_t id = total_deep_keys++;
                    mad_map[key] = id;
                    sm_map[key] = id;
                    sm_hashes.insert(SplitMixKeyHash{}(key));
                    if (wc.is_local(key)) {
                        wc.replace(key, static_cast<std::uint64_t>(id));
                    }
                }
            }
        }

        // Far-corner keys near 2^lvl - 1 (35-bit translation coordinates).
        for (Translation cx = top - 1; cx <= top; ++cx) {
            for (Translation cy = top - 1; cy <= top; ++cy) {
                for (Translation cz = top - 1; cz <= top; ++cz) {
                    const Key<3> key(lvl, Vector<Translation, 3>{cx, cy, cz});
                    const std::size_t id = total_deep_keys++;
                    mad_map[key] = id;
                    sm_map[key] = id;
                    sm_hashes.insert(SplitMixKeyHash{}(key));
                    if (wc.is_local(key)) {
                        wc.replace(key, static_cast<std::uint64_t>(id));
                    }
                }
            }
        }
    }
    wc.process_pending();
    world.gop.fence();

    std::size_t local_wc_count = 0;
    for (const auto& kv : wc) {
        (void)kv;
        ++local_wc_count;
    }
    world.gop.sum(local_wc_count);

    if (mad_map.size() != total_deep_keys || sm_map.size() != total_deep_keys ||
        sm_hashes.size() != total_deep_keys || local_wc_count != total_deep_keys) {
        deep_key_ok = false;
    }

    rep.check_true("Key<3>, WorldContainer & 64-bit SplitMix hash collision-free to Level=35",
                   deep_key_ok, "distinct_keys=" + std::to_string(total_deep_keys));

    rep.end_gate("Gate 1");
}

// ============================================================================
// Gate 2: Adversarial Unbalanced 3D Coarse-to-Fine Jump Test (Zero erf Ringing)
// ============================================================================

struct JumpCase {
    std::string tag;
    double cell_lo = 0.0;
    double cell_hi = 1.0;
    coord_3d center{0.485, 0.490, 0.495};
    double alpha = 150.0;
    int k = 8;
    int n_coarse = 4;
    int delta_j = 6;
    double eps = 1e-8;
};

void run_gate2(World& world, const SuiteOptions& opts, GateReporter& rep, ParityLedger& ledger) {
    rep.begin_gate("GATE 2: Adversarial Unbalanced 3D Coarse-to-Fine Jump (Zero erf Oscillations)");

    std::vector<JumpCase> cases;
    if (opts.small) {
        cases.push_back(JumpCase{"unit_k8_dj6", 0.0, 1.0, coord_3d{0.485, 0.490, 0.495}, 150.0, 8, 5, 6, 1e-8});
        cases.push_back(JumpCase{"box8_k8_dj6", -4.0, 4.0, coord_3d{-0.10, -0.08, -0.06}, 3.0, 8, 5, 6, 1e-8});
    } else {
        for (const int dj : {3, 6, 10}) {
            cases.push_back(JumpCase{"unit_k8_dj" + std::to_string(dj),
                                     0.0, 1.0, coord_3d{0.485, 0.490, 0.495}, 150.0, 8, 5, dj, 1e-8});
        }
        cases.push_back(JumpCase{"unit_k10_dj10", 0.0, 1.0, coord_3d{0.485, 0.490, 0.495}, 150.0, 10, 5, 10, 1e-8});
        cases.push_back(JumpCase{"box8_k8_dj6", -4.0, 4.0, coord_3d{-0.10, -0.08, -0.06}, 3.0, 8, 5, 6, 1e-8});
        cases.push_back(JumpCase{"box8_k10_dj10", -4.0, 4.0, coord_3d{-0.10, -0.08, -0.06}, 3.0, 10, 5, 10, 1e-8});
    }

    for (const JumpCase& jc : cases) {
        ScopedMraDefaults guard(world);
        FunctionDefaults<3>::set_cubic_cell(jc.cell_lo, jc.cell_hi);
        FunctionDefaults<3>::set_k(jc.k);
        FunctionDefaults<3>::set_thresh(jc.eps);
        FunctionDefaults<3>::set_truncate_mode(1);

        const int target_fine_level = jc.n_coarse + jc.delta_j;
        auto src_functor = std::make_shared<NormalizedGaussianFunctor>(jc.center, jc.alpha);
        const ExactGaussianCoulombFunctor exact_pot(jc.center, jc.alpha);

        // Unbalanced octree predicate:
        //   - Uniformly split down to n_coarse (so all boxes on x_sim >= 0.5, y_sim >= 0.5, z_sim >= 0.5
        //     remain coarse leaves at level n_coarse).
        //   - Inside octant x_sim < 0.5, y_sim < 0.5, z_sim < 0.5, drill an ungraded jump spike down to
        //     n_coarse + delta_j along boxes touching the face x_sim = 0.5, edge (0.5, 0.5, z), and
        //     corner (0.5, 0.5, 0.5) simultaneously with ZERO 2:1 balancing!
        auto should_split = [n_c = jc.n_coarse, n_f = target_fine_level](const Key<3>& key) -> bool {
            const int lvl = key.level();
            if (lvl < n_c) return true;
            if (lvl >= n_f) return false;

            const Translation mid = (Translation(1) << (lvl - 1));
            const auto& l = key.translation();
            // Corner spike touching (0.5, 0.5, 0.5) across face, edge, and corner simultaneously:
            const bool on_corner_spike = (l[0] == mid - 1 && l[1] == mid - 1 && l[2] == mid - 1);
            // Face-interior spike touching x_sim = 0.5 away from edges:
            const Translation face_idx = mid - 2;
            const bool on_face_spike = (l[0] == mid - 1 && l[1] == face_idx && l[2] == face_idx);
            // Edge-interior spike touching (0.5, 0.5, z) away from corner:
            const bool on_edge_spike = (l[0] == mid - 1 && l[1] == mid - 1 && l[2] == face_idx);
            return on_corner_spike || on_face_spike || on_edge_spike;
        };

        Function<double, 3> f_unbal = build_custom_octree_function(world, src_functor, jc.k, jc.eps, should_split);
        const std::size_t max_d = f_unbal.max_depth();
        const double f_norm = f_unbal.norm2();
        const double f_proj_rel = f_unbal.err(*src_functor) / std::max(1e-30, f_norm);

        const coord_3d mid_user{
            0.5 * (jc.cell_lo + jc.cell_hi),
            0.5 * (jc.cell_lo + jc.cell_hi),
            0.5 * (jc.cell_lo + jc.cell_hi)
        };
        const double L = jc.cell_hi - jc.cell_lo;
        const Level coarse_probe_depth = f_unbal.depthpt(
            coord_3d{mid_user[0] + 0.01 * L, mid_user[1] + 0.01 * L, mid_user[2] + 0.01 * L});
        const Level fine_probe_depth = f_unbal.depthpt(
            coord_3d{mid_user[0] - 1e-6 * L, mid_user[1] - 1e-6 * L, mid_user[2] - 1e-6 * L});

        rep.check_true(jc.tag + ": unbalanced jump tree constructed",
                       static_cast<int>(max_d) == target_fine_level &&
                       coarse_probe_depth == jc.n_coarse &&
                       fine_probe_depth == target_fine_level,
                       "coarse=" + std::to_string(coarse_probe_depth) +
                       " fine=" + std::to_string(fine_probe_depth));
        rep.check_le(jc.tag + ": source projection rel L2 error", f_proj_rel, 5e-9);

        WindowedConvolutionOptions wopts;
        wopts.enable_dynamic_refinement = true;
        auto win_coul = WindowedCoulombOperator(world, jc.eps, jc.k, wopts);
        Function<double, 3> V = win_coul(f_unbal);
        V.verify_tree();

        Function<double, 3> V_ref_proj = FunctionFactory<double, 3>(world)
                                             .functor(ExactGaussianCoulombFunctor(jc.center, jc.alpha))
                                             .k(jc.k)
                                             .thresh(jc.eps);
        const double v_exact_norm = std::max(1e-30, V_ref_proj.norm2());
        const double l2_rel_err = V.err(exact_pot) / v_exact_norm;
        rep.check_le(jc.tag + ": global Coulomb rel L2 error", l2_rel_err, 50.0 * jc.eps);

        // Dense pointwise scans across coarse-to-fine FACE, EDGE, and CORNER jumps.
        const double v_mid_scale = std::max(1e-30, std::abs(exact_pot(mid_user)));
        const std::array<double, 14> log_offsets_sim = {
            -1e-2, -3e-3, -1e-3, -3e-4, -1e-4, -3e-5, -1e-5,
            +1e-5, +3e-5, +1e-4, +3e-4, +1e-3, +3e-3, +1e-2
        };

        double max_face_pt_rel = 0.0;
        double max_edge_pt_rel = 0.0;
        double max_corner_pt_rel = 0.0;

        // Offset inside the finest level-15 corner leaf [0.5 - 2^{-15}, 0.5] (~[0.4999695, 0.5]):
        const double tang_sim_off = -1.5e-5 * L;
        for (const double s_sim : log_offsets_sim) {
            const double d = s_sim * L;
            // 1) Face crossing (x_sim = 0.5)
            const coord_3d p_face{mid_user[0] + d, mid_user[1] + tang_sim_off, mid_user[2] + tang_sim_off};
            max_face_pt_rel = std::max(max_face_pt_rel, std::abs(V(p_face) - exact_pot(p_face)) / v_mid_scale);

            // 2) Edge crossing (x_sim = 0.5, y_sim = 0.5)
            const coord_3d p_edge{mid_user[0] + d, mid_user[1] + d, mid_user[2] + tang_sim_off};
            max_edge_pt_rel = std::max(max_edge_pt_rel, std::abs(V(p_edge) - exact_pot(p_edge)) / v_mid_scale);

            // 3) Corner body-diagonal crossing (0.5, 0.5, 0.5)
            const coord_3d p_corner{mid_user[0] + d, mid_user[1] + d, mid_user[2] + d};
            max_corner_pt_rel = std::max(max_corner_pt_rel, std::abs(V(p_corner) - exact_pot(p_corner)) / v_mid_scale);
        }

        // Uniform 41-point micro-grid straddling the coarse-to-fine face jump to detect any erf oscillation kink.
        constexpr int n_micro = 41;
        const double ds_user = 5e-5 * L;
        std::array<double, n_micro> v_num_line{};
        std::array<double, n_micro> v_ex_line{};
        for (int m = 0; m < n_micro; ++m) {
            const double shift = static_cast<double>(m - (n_micro / 2)) * ds_user;
            const coord_3d pt{mid_user[0] + shift, mid_user[1] + tang_sim_off, mid_user[2] + tang_sim_off};
            v_num_line[static_cast<std::size_t>(m)] = V(pt);
            v_ex_line[static_cast<std::size_t>(m)] = exact_pot(pt);
            max_face_pt_rel = std::max(
                max_face_pt_rel,
                std::abs(v_num_line[static_cast<std::size_t>(m)] - v_ex_line[static_cast<std::size_t>(m)]) / v_mid_scale);
        }

        double max_d2_osc = 0.0;
        for (int m = 1; m + 1 < n_micro; ++m) {
            const auto im = static_cast<std::size_t>(m);
            const double d2_num = v_num_line[im + 1] - 2.0 * v_num_line[im] + v_num_line[im - 1];
            const double d2_ex = v_ex_line[im + 1] - 2.0 * v_ex_line[im] + v_ex_line[im - 1];
            max_d2_osc = std::max(max_d2_osc, std::abs(d2_num - d2_ex) / v_mid_scale);
        }

        rep.check_le(jc.tag + ": face jump pointwise rel error (zero erf jump)", max_face_pt_rel, 50.0 * jc.eps);
        rep.check_le(jc.tag + ": edge jump pointwise rel error (zero erf jump)", max_edge_pt_rel, 50.0 * jc.eps);
        rep.check_le(jc.tag + ": corner jump pointwise rel error (zero erf jump)", max_corner_pt_rel, 50.0 * jc.eps);
        rep.check_le(jc.tag + ": micro-grid 2nd-diff oscillation across jump", max_d2_osc, 100.0 * jc.eps);

        ledger.record("gate2_" + jc.tag + "_V_norm2", V.norm2());
        ledger.record("gate2_" + jc.tag + "_V_mid", V(mid_user));
    }

    rep.end_gate("Gate 2");
}

// ============================================================================
// Gate 3: 22+ Level Slater Cusp, Midplane Dynamic Refinement & BSH Cross-Check
// ============================================================================

struct CuspCase {
    std::string tag;
    double cell_lo = -1.5;
    double cell_hi = 2.5; // Centered at (0.5, 0.5, 0.5) => x=0.5 is exact root octree midplane x_sim=0.5!
    double Z = 20.0;
    int k = 8;
    double eps = 1e-8;
};

void run_gate3(World& world, const SuiteOptions& opts, GateReporter& rep, ParityLedger& ledger) {
    rep.begin_gate("GATE 3: 22+ Level Slater Cusp (Z=20,50), Midplane Dynamic Refinement & BSH");

    // ------------------------------------------------------------------------
    // Gate 3a: All-electron point-nucleus Slater cusp oracle next to x = 0.5
    // ------------------------------------------------------------------------
    const coord_3d R_A{0.5000004, 0.5, 0.5};
    std::vector<CuspCase> cusp_cases;
    if (opts.small) {
        cusp_cases.push_back(CuspCase{"Z20_k8_L4", -1.5, 2.5, 20.0, 8, 1e-8});
    } else {
        cusp_cases.push_back(CuspCase{"Z20_k8_L4", -1.5, 2.5, 20.0, 8, 1e-8});
        cusp_cases.push_back(CuspCase{"Z20_k10_L4", -1.5, 2.5, 20.0, 10, 1e-8});
        cusp_cases.push_back(CuspCase{"Z50_k8_L4", -1.5, 2.5, 50.0, 8, 1e-8});
        cusp_cases.push_back(CuspCase{"Z50_k10_L4", -1.5, 2.5, 50.0, 10, 1e-8});
        cusp_cases.push_back(CuspCase{"Z50_k10_unit", 0.0, 1.0, 50.0, 10, 1e-8});
    }

    for (const CuspCase& cc : cusp_cases) {
        ScopedMraDefaults guard(world);
        FunctionDefaults<3>::set_cubic_cell(cc.cell_lo, cc.cell_hi);
        FunctionDefaults<3>::set_k(cc.k);
        FunctionDefaults<3>::set_thresh(cc.eps);
        FunctionDefaults<3>::set_truncate_mode(1);
        FunctionDefaults<3>::set_max_refine_level(30);

        // Cap initial special_level at 11 to avoid FunctionImpl::project_refine_op's 2^special_level
        // newspecialpts vector duplication; adaptive dnorm >= truncate_tol refines the Slater cusp past 22+.
        const int init_special_lvl = 11;
        SlaterCuspDensityFunctor rho_fn(R_A, cc.Z, init_special_lvl);
        ExactSlaterCoulombFunctor v_ex_fn(R_A, cc.Z, init_special_lvl);

        Function<double, 3> rho = FunctionFactory<double, 3>(world)
                                      .functor(rho_fn)
                                      .k(cc.k)
                                      .thresh(cc.eps)
                                      .special_points({R_A})
                                      .special_level(init_special_lvl)
                                      .max_refine_level(30);
        const std::size_t rho_depth = rho.max_depth();
        const std::size_t rho_nodes = rho.tree_size();
        const double total_charge = rho.trace();

        WindowedConvolutionOptions wopts;
        wopts.max_level = 35;
        wopts.enable_dynamic_refinement = true;
        auto win_coul = WindowedCoulombOperator(world, cc.eps, cc.k, wopts);

        const double t0 = madness::wall_time();
        Function<double, 3> V = win_coul(rho);
        const double elapsed = madness::wall_time() - t0;
        V.verify_tree();

        Function<double, 3> V_ref = FunctionFactory<double, 3>(world)
                                        .functor(v_ex_fn)
                                        .k(cc.k)
                                        .thresh(cc.eps)
                                        .special_points({R_A})
                                        .special_level(init_special_lvl)
                                        .max_refine_level(30);
        const double v_ref_norm = std::max(1e-30, V_ref.norm2());
        const double v_l2_rel = V.err(v_ex_fn) / v_ref_norm;

        const double v_at_cusp = V(R_A);
        const double cusp_rel_err = std::abs(v_at_cusp - cc.Z) / cc.Z;

        const double e_hartree_num = 0.5 * madness::inner(rho, V);
        const double e_hartree_exact = (5.0 / 16.0) * cc.Z;
        const double hartree_rel_err = std::abs(e_hartree_num - e_hartree_exact) / e_hartree_exact;
        const std::size_t v_nodes = V.tree_size();

        if (world.rank() == 0) {
            std::cout << "  [INFO] " << cc.tag << ": depth=" << rho_depth
                      << " rho_nodes=" << rho_nodes << " V_nodes=" << v_nodes
                      << " Q=" << std::setprecision(12) << total_charge
                      << " V(R_A)=" << v_at_cusp << " (exact " << cc.Z << ")"
                      << " E_H=" << e_hartree_num << " (exact " << e_hartree_exact << ")"
                      << " time=" << std::fixed << std::setprecision(3) << elapsed << "s\n"
                      << std::defaultfloat;
        }

        rep.check_true(cc.tag + ": octree depth >= 22 around point nucleus",
                       rho_depth >= 22, "max_depth=" + std::to_string(rho_depth));
        rep.check_le(cc.tag + ": total charge |<1|rho> - 1|", std::abs(total_charge - 1.0), 50.0 * cc.eps);
        rep.check_le(cc.tag + ": Coulomb potential global rel L2 error", v_l2_rel, 50.0 * cc.eps);
        rep.check_le(cc.tag + ": nuclear cusp |V(R_A) - Z| / Z", cusp_rel_err, 50.0 * cc.eps);
        rep.check_le(cc.tag + ": Hartree self-energy |0.5<rho|V> - 5Z/16| / (5Z/16)",
                     hartree_rel_err, 50.0 * cc.eps);

        ledger.record("gate3_" + cc.tag + "_V_cusp", v_at_cusp);
        ledger.record("gate3_" + cc.tag + "_E_hartree", e_hartree_num);
        ledger.record("gate3_" + cc.tag + "_V_norm2", V.norm2());
    }

    // ------------------------------------------------------------------------
    // Gate 3b: Dynamic output refinement across octree midplane x = 0.5
    // ------------------------------------------------------------------------
    {
        ScopedMraDefaults guard(world);
        const int k = 8;
        const double eps = 1e-8;
        FunctionDefaults<3>::set_cubic_cell(0.0, 1.0);
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(eps);
        FunctionDefaults<3>::set_truncate_mode(1);

        // Steep nuclear Gaussian at R_nuc = (0.50003, 0.5, 0.5):
        //   alpha = 1.2e11 => on ALL of x <= 0.5, alpha * |r - R_nuc|^2 >= 108 => f(r) < 1e-31!
        //   So forcing x < 0.5 to remain coarse at level <= 3 incurs ZERO (< 1e-31) source error,
        //   while true Coulomb potential at probe_left = (0.49997, 0.5, 0.5) is ~16666.67 with steep 1/r curvature!
        const coord_3d R_nuc{0.50003, 0.5, 0.5};
        const double alpha_nuc = 1.2e11;
        auto nuc_src = std::make_shared<NormalizedGaussianFunctor>(R_nuc, alpha_nuc);
        const ExactGaussianCoulombFunctor nuc_v_ex(R_nuc, alpha_nuc, {R_nuc}, 11);

        auto split_right_only = [&R_nuc](const Key<3>& key) -> bool {
            const int lvl = key.level();
            if (lvl < 3) return true;
            if (lvl >= 21) return false;
            const Translation mid_x = (Translation(1) << (lvl - 1));
            // Strictly forbid any refinement on x < 0.5 (l[0] < 2^{lvl-1})!
            if (key.translation()[0] < mid_x) return false;

            const double scale = std::ldexp(1.0, lvl);
            const Translation tx = static_cast<Translation>(std::floor(R_nuc[0] * scale));
            const Translation ty = static_cast<Translation>(std::floor(R_nuc[1] * scale));
            const Translation tz = static_cast<Translation>(std::floor(R_nuc[2] * scale));
            const auto& l = key.translation();
            return std::abs(l[0] - tx) <= 2 && std::abs(l[1] - ty) <= 2 && std::abs(l[2] - tz) <= 2;
        };

        Function<double, 3> rho_half = build_custom_octree_function(world, nuc_src, k, eps, split_right_only);
        const coord_3d probe_left_near{0.49997, 0.50001, 0.50001};
        const coord_3d probe_left_mid{0.49985, 0.50001, 0.50001};
        const Level in_left_depth = rho_half.depthpt(probe_left_near);
        const std::size_t in_max_depth = rho_half.max_depth();

        WindowedConvolutionOptions wopts_dyn;
        wopts_dyn.max_level = 30;
        wopts_dyn.enable_dynamic_refinement = true;
        auto win_coul_dyn = WindowedCoulombOperator(world, eps, k, wopts_dyn);

        Function<double, 3> V_dyn = win_coul_dyn(rho_half);
        V_dyn.verify_tree();

        const Level out_left_depth = V_dyn.depthpt(probe_left_near);
        const double v_near_num = V_dyn(probe_left_near);
        const double v_near_ex = nuc_v_ex(probe_left_near);
        const double near_rel_err = std::abs(v_near_num - v_near_ex) / std::abs(v_near_ex);

        const double v_mid_num = V_dyn(probe_left_mid);
        const double v_mid_ex = nuc_v_ex(probe_left_mid);
        const double mid_rel_err = std::abs(v_mid_num - v_mid_ex) / std::abs(v_mid_ex);

        Function<double, 3> V_nuc_ref = FunctionFactory<double, 3>(world)
                                            .functor(nuc_v_ex)
                                            .k(k)
                                            .thresh(eps)
                                            .special_points({R_nuc})
                                            .special_level(11)
                                            .max_refine_level(30);
        const double dyn_l2_rel = V_dyn.err(nuc_v_ex) / std::max(1e-30, V_nuc_ref.norm2());

        if (world.rank() == 0) {
            std::cout << "  [INFO] Midplane dynamic refinement: in_left_depth=" << in_left_depth
                      << " -> out_left_depth=" << out_left_depth
                      << " (in_max_depth=" << in_max_depth << ")"
                      << " V(probe_left)=" << std::setprecision(12) << v_near_num
                      << " exact=" << v_near_ex << "\n" << std::defaultfloat;
        }

        rep.check_true("Dynamic output refinement promoted coarse x<0.5 leaves across midplane",
                       in_left_depth == 3 && out_left_depth >= 10,
                       "in_depth=" + std::to_string(in_left_depth) +
                       " out_depth=" + std::to_string(out_left_depth));
        rep.check_le("Dynamic output refinement x<0.5 near-midplane rel error", near_rel_err, 5e-6);
        rep.check_le("Dynamic output refinement x<0.5 mid-field rel error", mid_rel_err, 1e-6);
        rep.check_le("Dynamic output refinement global rel L2 error", dyn_l2_rel, 5e-6);

        ledger.record("gate3_dyn_V_left_near", v_near_num);
        ledger.record("gate3_dyn_V_norm2", V_dyn.norm2());
    }

    // ------------------------------------------------------------------------
    // Gate 3c: BSH cross-checks, vectorized apply, & SeparatedConvolution bridge
    // ------------------------------------------------------------------------
    {
        struct BSHCase {
            std::string tag;
            double cell_lo = 0.0;
            double cell_hi = 1.0;
            coord_3d center{0.47, 0.52, 0.49};
            double alpha = 150.0;
            double mu = 1.5;
            int k = 8;
            double eps = 1e-8;
        };

        std::vector<BSHCase> bsh_cases = {
            BSHCase{"bsh_unit_mu1.5_k8", 0.0, 1.0, coord_3d{0.47, 0.52, 0.49}, 150.0, 1.5, 8, 1e-8},
            BSHCase{"bsh_box12_mu0.5_k8", -6.0, 6.0, coord_3d{0.25, -0.35, 0.15}, 3.5, 0.5, 8, 1e-8}
        };
        if (!opts.small) {
            bsh_cases.push_back(
                BSHCase{"bsh_box12_mu1.5_k10", -6.0, 6.0, coord_3d{-0.30, 0.20, 0.40}, 4.0, 1.5, 10, 1e-8});
        }

        for (const BSHCase& bc : bsh_cases) {
            ScopedMraDefaults guard(world);
            FunctionDefaults<3>::set_cubic_cell(bc.cell_lo, bc.cell_hi);
            FunctionDefaults<3>::set_k(bc.k);
            FunctionDefaults<3>::set_thresh(bc.eps);
            FunctionDefaults<3>::set_truncate_mode(1);

            NormalizedGaussianFunctor src_fn(bc.center, bc.alpha);
            ExactGaussianBSHFunctor bsh_ex_fn(bc.center, bc.alpha, bc.mu);

            Function<double, 3> f = FunctionFactory<double, 3>(world).functor(src_fn).k(bc.k).thresh(bc.eps);
            auto win_bsh = WindowedBSHOperator3D(world, bc.mu, bc.eps, bc.k);
            Function<double, 3> u_win = madness::apply(win_bsh, f);
            u_win.verify_tree();

            SeparatedConvolution<double, 3> leg_bsh = BSHOperator3D(world, bc.mu, 1e-7, bc.eps);
            Function<double, 3> u_leg = leg_bsh(f);

            const double u_norm = std::max(1e-30, u_leg.norm2());
            const double err_vs_exact = u_win.err(bsh_ex_fn) / u_norm;
            const double leg_err_vs_exact = u_leg.err(bsh_ex_fn) / u_norm;
            const double err_vs_legacy = (u_win - u_leg).norm2() / u_norm;

            if (world.rank() == 0) {
                std::cout << "  [INFO] " << bc.tag << ": win_err_exact=" << std::scientific
                          << std::setprecision(4) << err_vs_exact
                          << " leg_err_exact=" << leg_err_vs_exact
                          << " diff_vs_leg=" << err_vs_legacy << std::defaultfloat << "\n";
            }

            rep.check_le(bc.tag + ": WindowedBSH vs analytical BSH rel L2 error", err_vs_exact, 50.0 * bc.eps);
            rep.check_le(bc.tag + ": WindowedBSH vs legacy BSHOperator3D rel L2 diff", err_vs_legacy, 100.0 * bc.eps);

            ledger.record("gate3_" + bc.tag + "_norm2", u_win.norm2());
        }

        // Vectorized apply(world, op, vec_f), apply(world, vec_ops, vec_f), and opt-in SeparatedConvolution bridge.
        ScopedMraDefaults guard(world);
        const int k = 8;
        const double eps = 1e-8;
        FunctionDefaults<3>::set_cubic_cell(-6.0, 6.0);
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(eps);
        FunctionDefaults<3>::set_truncate_mode(1);

        Function<double, 3> f1 = FunctionFactory<double, 3>(world)
            .functor(NormalizedGaussianFunctor(coord_3d{-0.3, 0.1, 0.2}, 3.0)).k(k).thresh(eps);
        Function<double, 3> f2 = FunctionFactory<double, 3>(world)
            .functor(NormalizedGaussianFunctor(coord_3d{0.4, -0.2, -0.1}, 2.5)).k(k).thresh(eps);
        Function<double, 3> f3 = FunctionFactory<double, 3>(world)
            .functor(NormalizedGaussianFunctor(coord_3d{0.05, 0.15, -0.25}, 36.0)).k(k).thresh(eps);
        Function<double, 3> f4 = FunctionFactory<double, 3>(world)
            .functor(NormalizedGaussianFunctor(coord_3d{-0.7, -0.4, 0.5}, 1.4)).k(k).thresh(eps);
        Function<double, 3> f5 = FunctionFactory<double, 3>(world)
            .functor(NormalizedGaussianFunctor(coord_3d{0.6, 0.5, -0.6}, 6.0)).k(k).thresh(eps);
        // 6-element heterogeneous batch (including sharp cusp-like f3 and duplicate handle f1 at slot 5)
        // exercising LPT cost-descending scheduling, 4-thread WindowedBatchWorkerTask, and tail donation.
        std::vector<Function<double, 3>> vf = {f1, f2, f3, f4, f5, f1};

        auto win_coul_ptr = WindowedCoulombOperatorPtr(world, eps, k);
        std::vector<Function<double, 3>> vc_vec = madness::apply(world, *win_coul_ptr, vf);
        std::vector<Function<double, 3>> vc_scalar;
        vc_scalar.reserve(vf.size());
        double vec_single_op_diff = 0.0;
        double vc_batch_norm_sum = 0.0;
        for (std::size_t idx = 0; idx < vf.size(); ++idx) {
            vc_scalar.push_back((*win_coul_ptr)(vf[idx]));
            const double snorm = std::max(1e-30, vc_scalar[idx].norm2());
            vec_single_op_diff = std::max(vec_single_op_diff, (vc_vec[idx] - vc_scalar[idx]).norm2() / snorm);
            vc_batch_norm_sum += vc_vec[idx].norm2();
        }
        const Function<double, 3>& vc_s1 = vc_scalar[0];
        const Function<double, 3>& vc_s2 = vc_scalar[1];

        const std::array<double, 6> mu_batch = {0.6, 1.1, 4.5, 0.45, 1.8, 0.85};
        std::vector<std::shared_ptr<WindowedSeparatedConvolution<double, 3>>> bsh_ops;
        bsh_ops.reserve(vf.size());
        for (std::size_t idx = 0; idx < vf.size(); ++idx) {
            bsh_ops.push_back((idx % 2 == 0)
                ? WindowedBSHOperatorPtr3D(world, mu_batch[idx], eps, k)
                : WindowedBSHOperatorPtr(world, mu_batch[idx], eps, k));
        }
        std::vector<Function<double, 3>> vb_vec = madness::apply(world, bsh_ops, vf);
        std::vector<Function<double, 3>> vb_scalar;
        vb_scalar.reserve(vf.size());
        double vec_multi_op_diff = 0.0;
        double vb_batch_norm_sum = 0.0;
        for (std::size_t idx = 0; idx < vf.size(); ++idx) {
            vb_scalar.push_back((*bsh_ops[idx])(vf[idx]));
            const double snorm = std::max(1e-30, vb_scalar[idx].norm2());
            vec_multi_op_diff = std::max(vec_multi_op_diff, (vb_vec[idx] - vb_scalar[idx]).norm2() / snorm);
            vb_batch_norm_sum += vb_vec[idx].norm2();
        }

        rep.check_le("Vector apply(world, op, vec_f) LPT batch parity vs scalar apply", vec_single_op_diff, 1e-13);
        rep.check_le("Vector apply(world, vec_ops, vec_f) LPT batch parity vs scalar apply", vec_multi_op_diff, 1e-13);
        ledger.record("gate3_vec_coul_batch_norm_sum", vc_batch_norm_sum);
        ledger.record("gate3_vec_bsh_batch_norm_sum", vb_batch_norm_sum);

        // Opt-in SeparatedConvolution::enable_windowed_ns(true) bridge check (single + batched vector apply).
        auto try_enable_windowed_ns = [](auto& op, bool flag) -> bool {
            if constexpr (requires { op.enable_windowed_ns(flag); }) {
                op.enable_windowed_ns(flag);
                return true;
            } else {
                (void)op;
                (void)flag;
                return false;
            }
        };

        SeparatedConvolution<double, 3> coul_bridge = CoulombOperator(world, 1e-5, eps);
        const bool bridge_api_coul = try_enable_windowed_ns(coul_bridge, true);
        rep.check_true("SeparatedConvolution::enable_windowed_ns(bool) API present", bridge_api_coul);
        if (bridge_api_coul) {
            const Function<double, 3> vc_bridge = madness::apply(coul_bridge, f1);
            const double bridge_coul_diff = (vc_bridge - vc_s1).norm2() / std::max(1e-30, vc_s1.norm2());

            std::vector<Function<double, 3>> vc_bridge_vec = madness::apply(world, coul_bridge, vf);
            double bridge_coul_vec_diff = 0.0;
            for (std::size_t idx = 0; idx < vf.size(); ++idx) {
                const double snorm = std::max(1e-30, vc_scalar[idx].norm2());
                bridge_coul_vec_diff = std::max(
                    bridge_coul_vec_diff, (vc_bridge_vec[idx] - vc_scalar[idx]).norm2() / snorm);
            }

            std::vector<std::shared_ptr<SeparatedConvolution<double, 3>>> leg_vec_ops;
            leg_vec_ops.reserve(vf.size());
            for (std::size_t idx = 0; idx < vf.size(); ++idx) {
                std::shared_ptr<SeparatedConvolution<double, 3>> bptr(
                    BSHOperatorPtr3D(world, mu_batch[idx], 1e-5, eps));
                try_enable_windowed_ns(*bptr, true);
                leg_vec_ops.push_back(std::move(bptr));
            }
            std::vector<Function<double, 3>> vb_bridge = madness::apply(world, leg_vec_ops, vf);
            double bridge_bsh_diff = 0.0;
            for (std::size_t idx = 0; idx < vf.size(); ++idx) {
                const double snorm = std::max(1e-30, vb_scalar[idx].norm2());
                bridge_bsh_diff = std::max(
                    bridge_bsh_diff, (vb_bridge[idx] - vb_scalar[idx]).norm2() / snorm);
            }

            rep.check_le("SeparatedConvolution::enable_windowed_ns(true) Coulomb parity", bridge_coul_diff, 1e-12);
            rep.check_le("SeparatedConvolution::enable_windowed_ns(true) batched Coulomb parity",
                         bridge_coul_vec_diff, 1e-12);
            rep.check_le("SeparatedConvolution::enable_windowed_ns(true) batched vector BSH parity",
                         bridge_bsh_diff, 1e-12);
        }

        // Complex Function<double_complex, 3> linearity check:
        //   (1) Exact < 1e-13 linearity on a shared octree T(f_cplx) (fixed tree, no norm-dependent promotion).
        //   (2) Adaptive consistency <= 50*eps vs independent adaptive trees T(f1) != T(f2).
        [&]<typename CplxT>() {
            if constexpr (requires(const WindowedSeparatedConvolution<double, 3>& op,
                                   const Function<CplxT, 3>& fc) {
                              madness::apply(op, fc);
                          }) {
                const CplxT imag_unit(0.0, 1.0);
                Function<CplxT, 3> f_cplx =
                    madness::convert<double, CplxT, 3>(f1) +
                    imag_unit * madness::convert<double, CplxT, 3>(f2);
                Function<double, 3> f_re_shared = madness::real(f_cplx);
                Function<double, 3> f_im_shared = madness::imag(f_cplx);

                WindowedConvolutionOptions wopts_fixed;
                wopts_fixed.enable_dynamic_refinement = false;
                auto win_coul_fixed = WindowedCoulombOperator(world, eps, k, wopts_fixed);

                Function<CplxT, 3> vc_cplx_fixed = madness::apply(win_coul_fixed, f_cplx);
                Function<CplxT, 3> vc_shared_ref =
                    madness::convert<double, CplxT, 3>(win_coul_fixed(f_re_shared)) +
                    imag_unit * madness::convert<double, CplxT, 3>(win_coul_fixed(f_im_shared));
                const double cplx_shared_rel = (vc_cplx_fixed - vc_shared_ref).norm2() /
                                               std::max(1e-30, vc_shared_ref.norm2());
                rep.check_le("Complex Function<double_complex,3> shared-tree linearity", cplx_shared_rel, 1e-13);

                Function<CplxT, 3> vc_cplx = madness::apply(*win_coul_ptr, f_cplx);
                Function<CplxT, 3> vc_distinct_ref =
                    madness::convert<double, CplxT, 3>(vc_s1) +
                    imag_unit * madness::convert<double, CplxT, 3>(vc_s2);
                const double cplx_distinct_rel = (vc_cplx - vc_distinct_ref).norm2() /
                                                 std::max(1e-30, vc_distinct_ref.norm2());
                rep.check_le("Complex Function<double_complex,3> distinct-tree consistency",
                             cplx_distinct_rel, 50.0 * eps);
            }
        }.template operator()<::double_complex>();
    }

    rep.end_gate("Gate 3");
}

// ============================================================================
// Gate 4: Multi-Rank Parity, MemoryMeasurer Zero-Leak & Legacy Benchmark
// ============================================================================

std::pair<long, double> snapshot_function_memory(World& world) {
    world.gop.fence();
    const auto mmap = MemoryMeasurer::measure_and_print(world);
    long total_funcs = 0;
    double total_gb = 0.0;
    if (world.rank() == 0) {
        total_gb = MemoryMeasurer::total_memory(mmap);
        for (const auto& [key, info] : mmap) {
            (void)key;
            total_funcs += info.num_functions;
        }
    }
    world.gop.broadcast(total_funcs, 0);
    world.gop.broadcast(total_gb, 0);
    return {total_funcs, total_gb};
}

void run_gate4(World& world, const SuiteOptions& opts, GateReporter& rep, ParityLedger& ledger) {
    rep.begin_gate("GATE 4: Multi-Rank Checksums, MemoryMeasurer Zero-Leak & Legacy Benchmark");

    // ------------------------------------------------------------------------
    // Gate 4a: Multi-rank parity dump / verification
    // ------------------------------------------------------------------------
    if (world.rank() == 0) {
        std::cout << "  [INFO] Recorded " << ledger.entries().size() << " deterministic observables on np="
                  << world.size() << ":\n";
        for (const auto& [k, v] : ledger.entries()) {
            std::cout << "    checksum[" << k << "] = " << std::setprecision(16) << v
                      << std::defaultfloat << "\n";
        }
    }
    if (!opts.write_parity_path.empty()) {
        const bool wrote = ledger.write_to_file(world, opts.write_parity_path);
        rep.check_true("Wrote deterministic parity reference file: " + opts.write_parity_path, wrote);
    }
    if (!opts.check_parity_path.empty()) {
        ledger.check_against_file(world, opts.check_parity_path, rep, 1e-13);
    }

    // ------------------------------------------------------------------------
    // Gate 4b: MemoryMeasurer & WorldObject zero-leak audit across repeated calls
    // ------------------------------------------------------------------------
    {
        ScopedMraDefaults guard(world);
        const int k = 8;
        const double eps = 1e-7;
        FunctionDefaults<3>::set_cubic_cell(-4.0, 4.0);
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(eps);
        FunctionDefaults<3>::set_truncate_mode(1);

        Function<double, 3> f_src = FunctionFactory<double, 3>(world)
            .functor(NormalizedGaussianFunctor(coord_3d{0.15, -0.20, 0.10}, 6.0))
            .k(k)
            .thresh(eps);

        // One warmup call so static caches (FunctionCommonData, GaussianConvolution1DCache) settle.
        {
            auto warm_coul = WindowedCoulombOperator(world, eps, k);
            auto warm_bsh = WindowedBSHOperator3D(world, 0.8, eps, k);
            Function<double, 3> g1 = warm_coul(f_src);
            Function<double, 3> g2 = warm_bsh(f_src);
            (void)g1;
            (void)g2;
        }
        world.gop.fence();

        const auto [base_funcs, base_gb] = snapshot_function_memory(world);
        std::size_t base_world_objs = world.get_object_ids().size();
        world.gop.max(base_world_objs);

        const int n_iters = opts.small ? 3 : 5;
        double accum_trace = 0.0;
        for (int iter = 0; iter < n_iters; ++iter) {
            auto op_c = WindowedCoulombOperator(world, eps, k);
            auto op_b = WindowedBSHOperator3D(world, 0.8, eps, k);
            Function<double, 3> vc = op_c(f_src);
            Function<double, 3> vb = op_b(f_src);
            accum_trace += vc.trace() + vb.trace();
        }
        world.gop.fence();

        const auto [after_funcs, after_gb] = snapshot_function_memory(world);
        std::size_t after_world_objs = world.get_object_ids().size();
        world.gop.max(after_world_objs);

        const long leaked_funcs = after_funcs - base_funcs;
        const double leaked_gb = std::abs(after_gb - base_gb);
        const long leaked_world_objs = static_cast<long>(after_world_objs) - static_cast<long>(base_world_objs);

        rep.check_true("Zero leaked FunctionImpl objects after repeated Windowed NS applies",
                       leaked_funcs == 0 && std::isfinite(accum_trace),
                       "base=" + std::to_string(base_funcs) + " after=" + std::to_string(after_funcs));
        rep.check_le("MemoryMeasurer leaked memory_GB delta", leaked_gb, 1e-15);
        rep.check_true("Zero leaked WorldObject IDs (halo/context containers)",
                       leaked_world_objs == 0,
                       "base_objs=" + std::to_string(base_world_objs) +
                       " after_objs=" + std::to_string(after_world_objs));
    }

    // ------------------------------------------------------------------------
    // Gate 4b.2: Deep special_level(23) nuclear cusp projection regression
    //            (verifies FunctionImpl::project_refine_op O(1) specialpts
    //             propagation without 2^{special_level} duplication + zero leak)
    // ------------------------------------------------------------------------
    {
        ScopedMraDefaults guard(world);
        const int k = 8;
        const double eps = 1e-8;
        const double Z_sp = 20.0;
        const coord_3d R_sp{0.5000004, 0.5, 0.5};

        FunctionDefaults<3>::set_cubic_cell(-1.5, 2.5);
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(eps);
        FunctionDefaults<3>::set_truncate_mode(1);
        FunctionDefaults<3>::set_max_refine_level(30);

        {
            const auto& cdata_warm = FunctionCommonData<double, 3>::get(k);
            (void)cdata_warm;
        }
        world.gop.fence();

        const auto [sp_base_funcs, sp_base_gb] = snapshot_function_memory(world);
        std::size_t sp_base_objs = world.get_object_ids().size();
        world.gop.max(sp_base_objs);

        std::size_t sp_depth = 0;
        std::size_t sp_nodes = 0;
        double sp_norm2 = 0.0;
        double sp_charge = 0.0;
        const double exact_sp_norm2 =
            std::sqrt((Z_sp * Z_sp * Z_sp) / (8.0 * madness::constants::pi));
        double dt_sp_proj = 0.0;
        {
            SlaterCuspDensityFunctor cusp_sp23_fn(R_sp, Z_sp, 23);
            const double t_sp0 = madness::wall_time();
            Function<double, 3> rho_sp23 = FunctionFactory<double, 3>(world)
                                               .functor(cusp_sp23_fn)
                                               .k(k)
                                               .thresh(eps)
                                               .truncate_mode(1)
                                               .special_points({R_sp})
                                               .special_level(23)
                                               .max_refine_level(30);
            world.gop.fence();
            dt_sp_proj = madness::wall_time() - t_sp0;

            rho_sp23.verify_tree();
            sp_depth = rho_sp23.max_depth();
            sp_nodes = rho_sp23.tree_size();
            sp_norm2 = rho_sp23.norm2();
            sp_charge = rho_sp23.trace();
        }
        world.gop.fence();

        const auto [sp_after_funcs, sp_after_gb] = snapshot_function_memory(world);
        std::size_t sp_after_objs = world.get_object_ids().size();
        world.gop.max(sp_after_objs);

        const long sp_leaked_funcs = sp_after_funcs - sp_base_funcs;
        const double sp_leaked_gb = std::abs(sp_after_gb - sp_base_gb);
        const long sp_leaked_objs =
            static_cast<long>(sp_after_objs) - static_cast<long>(sp_base_objs);
        const double sp_norm_rel_err =
            std::abs(sp_norm2 - exact_sp_norm2) / exact_sp_norm2;
        const double sp_charge_err = std::abs(sp_charge - 1.0);

        if (world.rank() == 0) {
            std::cout << "  [INFO] special_level(23) cusp projection: depth=" << sp_depth
                      << " nodes=" << sp_nodes
                      << " norm2=" << std::setprecision(12) << sp_norm2
                      << " (exact " << exact_sp_norm2 << ")"
                      << " Q=" << sp_charge
                      << " time=" << std::fixed << std::setprecision(3) << dt_sp_proj << "s\n"
                      << std::defaultfloat;
        }

        rep.check_true("special_level(23) cusp projection: finite norm & charge",
                       std::isfinite(sp_norm2) && std::isfinite(sp_charge));
        rep.check_true("special_level(23) cusp projection: max_depth >= 23",
                       sp_depth >= 23, "max_depth=" + std::to_string(sp_depth));
        rep.check_le("special_level(23) cusp projection: L2 norm rel error",
                     sp_norm_rel_err, 50.0 * eps);
        rep.check_le("special_level(23) cusp projection: |<1|rho> - 1|",
                     sp_charge_err, 50.0 * eps);
        rep.check_true("special_level(23) cusp projection: zero leaked objects",
                       sp_leaked_funcs == 0 && sp_leaked_objs == 0,
                       "leaked_funcs=" + std::to_string(sp_leaked_funcs) +
                       " leaked_objs=" + std::to_string(sp_leaked_objs));
        rep.check_le("special_level(23) cusp projection: MemoryMeasurer leaked_GB",
                     sp_leaked_gb, 1e-15);
    }

    // ------------------------------------------------------------------------
    // Gate 4c: Runtime & tree-size benchmark vs CoulombOperator / BSHOperator3D
    // ------------------------------------------------------------------------
    {
        ScopedMraDefaults guard(world);
        const int k = 8;
        const double eps = 1e-7;
        const double Z_bench = 30.0;
        const coord_3d R_bench{0.5000004, 0.5, 0.5};

        FunctionDefaults<3>::set_cubic_cell(-1.5, 2.5);
        FunctionDefaults<3>::set_k(k);
        FunctionDefaults<3>::set_thresh(eps);
        FunctionDefaults<3>::set_truncate_mode(1);
        FunctionDefaults<3>::set_max_refine_level(30);

        SlaterCuspDensityFunctor rho_fn(R_bench, Z_bench, 11);
        Function<double, 3> rho = FunctionFactory<double, 3>(world)
                                      .functor(rho_fn)
                                      .k(k)
                                      .thresh(eps)
                                      .special_points({R_bench})
                                      .special_level(11)
                                      .max_refine_level(30);

        // Construct operators outside timed apply regions (separate setup vs apply timings).
        world.gop.fence();
        const double t_wc_setup0 = madness::wall_time();
        auto win_coul = WindowedCoulombOperator(world, eps, k);
        world.gop.fence();
        const double dt_wc_setup = madness::wall_time() - t_wc_setup0;

        const double t_wc_app0 = madness::wall_time();
        Function<double, 3> V_win = win_coul(rho);
        world.gop.fence();
        const double dt_wc_apply = madness::wall_time() - t_wc_app0;

        // Legacy CoulombOperator with lo = 1e-7 (~2e-8 cusp resolution on L=4, depth ~25).
        const double t_lc_setup0 = madness::wall_time();
        SeparatedConvolution<double, 3> leg_coul = CoulombOperator(world, 1e-7, eps);
        world.gop.fence();
        const double dt_lc_setup = madness::wall_time() - t_lc_setup0;

        const double t_lc_app0 = madness::wall_time();
        Function<double, 3> V_leg = leg_coul(rho);
        world.gop.fence();
        const double dt_lc_apply = madness::wall_time() - t_lc_app0;

        const double coul_cross_rel = (V_win - V_leg).norm2() / std::max(1e-30, V_leg.norm2());

        // BSH benchmark at mu = 1.0.
        const double t_wb_app0 = madness::wall_time();
        auto win_bsh = WindowedBSHOperator3D(world, 1.0, eps, k);
        Function<double, 3> B_win = win_bsh(rho);
        world.gop.fence();
        const double dt_wb_total = madness::wall_time() - t_wb_app0;

        const double t_lb_app0 = madness::wall_time();
        SeparatedConvolution<double, 3> leg_bsh = BSHOperator3D(world, 1.0, 1e-7, eps);
        Function<double, 3> B_leg = leg_bsh(rho);
        world.gop.fence();
        const double dt_lb_total = madness::wall_time() - t_lb_app0;

        const double bsh_cross_rel = (B_win - B_leg).norm2() / std::max(1e-30, B_leg.norm2());

        const std::size_t bench_depth = rho.max_depth();
        const std::size_t bench_rho_nodes = rho.tree_size();
        const std::size_t v_win_nodes = V_win.tree_size();
        const std::size_t v_leg_nodes = V_leg.tree_size();
        const std::size_t b_win_nodes = B_win.tree_size();
        const std::size_t b_leg_nodes = B_leg.tree_size();

        if (world.rank() == 0) {
            std::cout << "  [BENCH] Deep Slater cusp Z=" << Z_bench << " (depth=" << bench_depth
                      << ", rho_nodes=" << bench_rho_nodes << ", np=" << world.size() << "):\n"
                      << "    WindowedCoulomb: setup=" << std::fixed << std::setprecision(4) << dt_wc_setup
                      << "s apply=" << dt_wc_apply << "s out_nodes=" << v_win_nodes << "\n"
                      << "    LegacyCoulomb:   setup=" << dt_lc_setup
                      << "s apply=" << dt_lc_apply << "s out_nodes=" << v_leg_nodes
                      << " (apply speedup=" << std::setprecision(2)
                      << (dt_lc_apply / std::max(1e-9, dt_wc_apply)) << "x)\n"
                      << "    WindowedBSH(total)=" << std::setprecision(4) << dt_wb_total
                      << "s out_nodes=" << b_win_nodes
                      << " vs LegacyBSH(total)=" << dt_lb_total << "s out_nodes=" << b_leg_nodes
                      << " (speedup=" << std::setprecision(2)
                      << (dt_lb_total / std::max(1e-9, dt_wb_total)) << "x)\n"
                      << std::defaultfloat;
        }

        rep.check_le("Benchmark Z=30 Coulomb cross-operator rel L2 diff", coul_cross_rel, 100.0 * eps);
        rep.check_le("Benchmark Z=30 BSH cross-operator rel L2 diff", bsh_cross_rel, 100.0 * eps);
    }

    rep.end_gate("Gate 4");
}

SuiteOptions parse_cli(int argc, char** argv) {
    SuiteOptions opts;
    if (std::getenv("MAD_SMALL_TESTS") != nullptr) {
        opts.small = true;
    }

    bool explicit_gate = false;
    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--small") {
            opts.small = true;
        } else if (arg == "--full") {
            opts.small = false;
        } else if (arg.rfind("--write-parity=", 0) == 0) {
            opts.write_parity_path = arg.substr(std::strlen("--write-parity="));
        } else if (arg.rfind("--check-parity=", 0) == 0) {
            opts.check_parity_path = arg.substr(std::strlen("--check-parity="));
        } else if (arg.rfind("--gate=", 0) == 0) {
            if (!explicit_gate) {
                opts.run_gate1 = false;
                opts.run_gate2 = false;
                opts.run_gate3 = false;
                opts.run_gate4 = false;
                explicit_gate = true;
            }
            const std::string g = arg.substr(std::strlen("--gate="));
            if (g == "1") opts.run_gate1 = true;
            if (g == "2") opts.run_gate2 = true;
            if (g == "3") opts.run_gate3 = true;
            if (g == "4") opts.run_gate4 = true;
        }
    }
    return opts;
}

} // namespace

int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    int exit_code = 0;
    try {
        madness::startup(world, argc, argv, false);
        world.gop.fence();

        const SuiteOptions opts = parse_cli(argc, argv);
        if (world.rank() == 0) {
            std::cout << "========================================================================\n"
                      << "MADNESS Windowed Non-Standard 3D Convolution — 4-Gate Adversarial Suite\n"
                      << "  MPI ranks = " << world.size()
                      << ", mode = " << (opts.small ? "SMALL (CTest short)" : "FULL (adversarial)") << "\n"
                      << "========================================================================\n";
        }

        GateReporter rep(world);
        ParityLedger ledger;
        const double t_suite0 = madness::wall_time();

        if (opts.run_gate1) run_gate1(world, opts, rep);
        if (opts.run_gate2) run_gate2(world, opts, rep, ledger);
        if (opts.run_gate3) run_gate3(world, opts, rep, ledger);
        if (opts.run_gate4) run_gate4(world, opts, rep, ledger);

        world.gop.fence();
        const double dt_suite = madness::wall_time() - t_suite0;
        exit_code = (rep.failed_checks() == 0) ? 0 : 1;

        if (world.rank() == 0) {
            std::cout << "\n========================================================================\n"
                      << "FINAL SUMMARY: " << (rep.total_checks() - rep.failed_checks()) << " / "
                      << rep.total_checks() << " checks PASSED, " << rep.failed_checks()
                      << " FAILED in " << std::fixed << std::setprecision(3) << dt_suite << " s\n"
                      << "========================================================================\n"
                      << std::defaultfloat;
        }
    } catch (const madness::MadnessException& e) {
        std::cerr << "[RANK " << world.rank() << "] Caught MadnessException: " << e << "\n";
        exit_code = 2;
    } catch (const std::exception& e) {
        std::cerr << "[RANK " << world.rank() << "] Caught std::exception: " << e.what() << "\n";
        exit_code = 3;
    } catch (...) {
        std::cerr << "[RANK " << world.rank() << "] Caught unknown exception\n";
        exit_code = 4;
    }

    world.gop.max(exit_code);
    world.gop.fence();
    madness::finalize();
    return exit_code;
}
