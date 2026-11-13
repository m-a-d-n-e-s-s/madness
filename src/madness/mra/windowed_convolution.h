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

#ifndef MADNESS_MRA_WINDOWED_CONVOLUTION_H__INCLUDED
#define MADNESS_MRA_WINDOWED_CONVOLUTION_H__INCLUDED

/// \file windowed_convolution.h
/// \brief Public C++ operator hierarchy and factory functions for Tail-Completed,
///        Two-Threshold Windowed Non-Standard (NS) 3D Convolutions (Coulomb, BSH, Yukawa).

#include <madness/mra/mra.h>
#include <madness/mra/windowed_quadrature.h>

#include <cstddef>
#include <memory>
#include <type_traits>
#include <vector>

namespace madness {

#ifndef MADNESS_WINDOWED_CONV_FWD_DECLARED
#define MADNESS_WINDOWED_CONV_FWD_DECLARED
template <typename Q = double, std::size_t NDIM = 3>
class WindowedSeparatedConvolution;
#endif

class WindowedStencilCache1D;

/// Configuration options for `WindowedSeparatedConvolution` and its factories.
struct WindowedConvolutionOptions {
    int oversampling_M = 0;                ///< Terms per dyadic octave (0 => auto from eps).
    int gc_rank_Q = 0;                     ///< Gauss-Christoffel rank Q (0 => auto in {4, 5}).
    int max_level = 35;                    ///< Maximum precomputed octree level (extendable on demand).
    bool enable_dynamic_refinement = true; ///< Promote fine output nodes across coarse halos when needed.
    double dynamic_refine_factor = 1.0;    ///< Multiplier on function truncation threshold for dynamic refinement.
    double eps_loc = 0.0;                  ///< Optional explicit local threshold override (> 0 overrides function thresh).
};

/// Tail-Completed, Two-Threshold Windowed Non-Standard (NS) 3D Convolution Operator.
///
/// Implements Algorithm 1 of `multiwavelet_ns_v7`:
///   - Stage 1: Exact NS tree transform (`nonstandard_with_leaves`) + interior-only top-down
///              `virt_cache` halo prolongation (zero exterior polynomial extrapolation).
///   - Stage 2: Per-leaf Gauss-Christoffel compressed local moment pass (`Q <= 5` terms).
///   - Stage 3: Root-level (`n = 0`) ultra-diffuse and discrete constant low-tail pass.
///   - Stage 4: Scale-windowed (`W` levels, `<= W * M` active Gaussians/level) sum-factorized
///              NS application with Theorem 3.1 boundary-halo commutation.
///   - Stage 5: Bottom-up coarse-leaf spillover restriction + standard tree reconstruction.
template <typename Q, std::size_t NDIM>
class WindowedSeparatedConvolution {
    static_assert(NDIM == 3 && std::is_same_v<Q, double>,
                  "WindowedSeparatedConvolution currently supports real-valued 3D kernels (Q = double, NDIM = 3).");

public:
    using opT = Q;
    using keyT = Key<NDIM>;
    static constexpr std::size_t opdim = NDIM;

    bool print_timings = false;

    WindowedSeparatedConvolution(World& world,
                                 WindowedKernelKind kind,
                                 double mu_user,
                                 double eps,
                                 int k = FunctionDefaults<NDIM>::get_k(),
                                 const WindowedConvolutionOptions& opts = {});

    WindowedSeparatedConvolution(World& world,
                                 WindowedKernelKind kind,
                                 double mu_user,
                                 double eps,
                                 const WindowedConvolutionOptions& opts)
        : WindowedSeparatedConvolution(world, kind, mu_user, eps, FunctionDefaults<NDIM>::get_k(), opts) {}

    WindowedSeparatedConvolution(const WindowedSeparatedConvolution&) = default;
    WindowedSeparatedConvolution(WindowedSeparatedConvolution&&) noexcept = default;
    WindowedSeparatedConvolution& operator=(const WindowedSeparatedConvolution&) = default;
    WindowedSeparatedConvolution& operator=(WindowedSeparatedConvolution&&) noexcept = default;
    ~WindowedSeparatedConvolution() = default;

    World& get_world() const noexcept { return *world_; }
    World& world() const noexcept { return *world_; }

    int get_k() const noexcept { return k_; }
    int k() const noexcept { return k_; }

    double eps() const noexcept { return eps_; }
    double thresh() const noexcept { return eps_; }

    double mu() const noexcept { return mu_user_; }
    double mu_user() const noexcept { return mu_user_; }
    double gamma() const noexcept { return mu_user_; }

    double cell_width() const noexcept { return cell_width_L_; }
    WindowedKernelKind kind() const noexcept { return kind_; }

    const WindowedQuadraturePlan& plan() const noexcept { return *plan_; }
    const WindowedStencilCache1D& stencil_cache() const noexcept { return *stencil_cache_; }
    const WindowedConvolutionOptions& options() const noexcept { return opts_; }

    bool& destructive() noexcept { return destructive_; }
    const bool& destructive() const noexcept { return destructive_; }

    bool modified() const noexcept { return false; }

    bool dynamic_refinement() const noexcept { return opts_.enable_dynamic_refinement; }
    bool enable_dynamic_refinement() const noexcept { return opts_.enable_dynamic_refinement; }

    WindowedSeparatedConvolution& enable_dynamic_refinement(bool flag) noexcept {
        opts_.enable_dynamic_refinement = flag;
        return *this;
    }

    WindowedSeparatedConvolution& set_dynamic_refinement(bool flag) noexcept {
        opts_.enable_dynamic_refinement = flag;
        return *this;
    }

    WindowedSeparatedConvolution& set_dynamic_refine_factor(double factor) noexcept {
        opts_.dynamic_refine_factor = factor;
        return *this;
    }

    WindowedSeparatedConvolution& set_eps_loc(double tol) noexcept {
        opts_.eps_loc = tol;
        return *this;
    }

    void print_timer() const noexcept {}
    void reset_timer() const noexcept {}

    Function<double, 3> apply_func(const Function<double, 3>& f, bool fence = true) const;
    Function<double_complex, 3> apply_func(const Function<double_complex, 3>& f, bool fence = true) const;

    Function<double, 3> apply_on_world(World& world, const Function<double, 3>& f, bool fence = true) const;
    Function<double_complex, 3> apply_on_world(World& world, const Function<double_complex, 3>& f, bool fence = true) const;

    static std::vector<Function<double, 3>> apply_batch_on_world(
        World& world,
        const std::vector<const WindowedSeparatedConvolution<Q, NDIM>*>& ops,
        const std::vector<Function<double, 3>>& f,
        bool fence = true,
        bool ignore_op_destructive = false,
        bool force_destructive = false);

    static std::vector<Function<double_complex, 3>> apply_batch_on_world(
        World& world,
        const std::vector<const WindowedSeparatedConvolution<Q, NDIM>*>& ops,
        const std::vector<Function<double_complex, 3>>& f,
        bool fence = true,
        bool ignore_op_destructive = false,
        bool force_destructive = false);

    template <typename R>
    Function<TENSOR_RESULT_TYPE(Q, R), NDIM>
    operator()(const Function<R, NDIM>& f, bool fence = true) const {
        static_assert(std::is_same_v<R, double> || std::is_same_v<R, double_complex>,
                      "WindowedSeparatedConvolution supports Function<double, 3> and Function<double_complex, 3>.");
        return apply_func(f, fence);
    }

    template <typename R>
    std::vector<Function<TENSOR_RESULT_TYPE(Q, R), NDIM>>
    operator()(const std::vector<Function<R, NDIM>>& f, bool fence = true) const;

private:
    World* world_ = nullptr;
    WindowedKernelKind kind_ = WindowedKernelKind::Coulomb;
    double mu_user_ = 0.0;
    double eps_ = 1e-8;
    int k_ = 8;
    double cell_width_L_ = 1.0;
    WindowedConvolutionOptions opts_{};
    bool destructive_ = false;

    std::shared_ptr<const WindowedQuadraturePlan> plan_;
    std::shared_ptr<const WindowedStencilCache1D> stencil_cache_;
};

extern template class WindowedSeparatedConvolution<double, 3>;

/// Apply a `WindowedSeparatedConvolution` to a single 3D MRA Function.
template <typename Q, typename R, std::size_t NDIM>
inline Function<TENSOR_RESULT_TYPE(Q, R), NDIM>
apply(const WindowedSeparatedConvolution<Q, NDIM>& op,
      const Function<R, NDIM>& f,
      bool fence = true) {
    return op.apply_func(f, fence);
}

/// Apply a single `WindowedSeparatedConvolution` to a vector of 3D MRA Functions on `world`.
template <typename Q, typename R, std::size_t NDIM>
inline std::vector<Function<TENSOR_RESULT_TYPE(Q, R), NDIM>>
apply(World& world,
      const WindowedSeparatedConvolution<Q, NDIM>& op,
      const std::vector<Function<R, NDIM>>& f,
      bool fence = true) {
    std::vector<const WindowedSeparatedConvolution<Q, NDIM>*> ops(f.size(), &op);
    return WindowedSeparatedConvolution<Q, NDIM>::apply_batch_on_world(
        world, ops, f, fence, false, op.destructive());
}

/// Apply a single `WindowedSeparatedConvolution` to a vector of 3D MRA Functions on `op.get_world()`.
template <typename Q, typename R, std::size_t NDIM>
inline std::vector<Function<TENSOR_RESULT_TYPE(Q, R), NDIM>>
apply(const WindowedSeparatedConvolution<Q, NDIM>& op,
      const std::vector<Function<R, NDIM>>& f,
      bool fence = true) {
    return apply(op.get_world(), op, f, fence);
}

/// Apply a vector of `WindowedSeparatedConvolution` shared pointers pairwise to a vector of Functions (3-arg overload matching vmra.h).
template <typename Q, typename R, std::size_t NDIM>
inline std::vector<Function<TENSOR_RESULT_TYPE(Q, R), NDIM>>
apply(World& world,
      const std::vector<std::shared_ptr<WindowedSeparatedConvolution<Q, NDIM>>>& ops,
      const std::vector<Function<R, NDIM>>& f) {
    MADNESS_CHECK(ops.size() == f.size());
    std::vector<const WindowedSeparatedConvolution<Q, NDIM>*> op_ptrs(ops.size(), nullptr);
    for (std::size_t i = 0; i < ops.size(); ++i) {
        MADNESS_CHECK(static_cast<bool>(ops[i]));
        op_ptrs[i] = ops[i].get();
    }
    return WindowedSeparatedConvolution<Q, NDIM>::apply_batch_on_world(
        world, op_ptrs, f, true, false, false);
}

/// Apply a vector of `WindowedSeparatedConvolution` shared pointers pairwise to a vector of Functions (4-arg overload with fence).
template <typename Q, typename R, std::size_t NDIM>
inline std::vector<Function<TENSOR_RESULT_TYPE(Q, R), NDIM>>
apply(World& world,
      const std::vector<std::shared_ptr<WindowedSeparatedConvolution<Q, NDIM>>>& ops,
      const std::vector<Function<R, NDIM>>& f,
      bool fence) {
    MADNESS_CHECK(ops.size() == f.size());
    std::vector<const WindowedSeparatedConvolution<Q, NDIM>*> op_ptrs(ops.size(), nullptr);
    for (std::size_t i = 0; i < ops.size(); ++i) {
        MADNESS_CHECK(static_cast<bool>(ops[i]));
        op_ptrs[i] = ops[i].get();
    }
    return WindowedSeparatedConvolution<Q, NDIM>::apply_batch_on_world(
        world, op_ptrs, f, fence, false, false);
}

/// Apply a vector of `WindowedSeparatedConvolution` instances pairwise to a vector of Functions.
template <typename Q, typename R, std::size_t NDIM>
inline std::vector<Function<TENSOR_RESULT_TYPE(Q, R), NDIM>>
apply(World& world,
      const std::vector<WindowedSeparatedConvolution<Q, NDIM>>& ops,
      const std::vector<Function<R, NDIM>>& f,
      bool fence = true) {
    MADNESS_CHECK(ops.size() == f.size());
    std::vector<const WindowedSeparatedConvolution<Q, NDIM>*> op_ptrs(ops.size(), nullptr);
    for (std::size_t i = 0; i < ops.size(); ++i) {
        op_ptrs[i] = &ops[i];
    }
    return WindowedSeparatedConvolution<Q, NDIM>::apply_batch_on_world(
        world, op_ptrs, f, fence, false, false);
}

template <typename Q, std::size_t NDIM>
template <typename R>
inline std::vector<Function<TENSOR_RESULT_TYPE(Q, R), NDIM>>
WindowedSeparatedConvolution<Q, NDIM>::operator()(const std::vector<Function<R, NDIM>>& f, bool fence) const {
    return madness::apply(*world_, *this, f, fence);
}

// ============================================================================
// Factory Functions: Coulomb (1/r), BSH (exp(-mu*r)/(4*pi*r)), Yukawa (exp(-mu*r)/r)
// ============================================================================

inline WindowedSeparatedConvolution<double, 3>
WindowedCoulombOperator(World& world,
                        double eps,
                        int k = FunctionDefaults<3>::get_k(),
                        const WindowedConvolutionOptions& opts = {}) {
    return WindowedSeparatedConvolution<double, 3>(
        world, WindowedKernelKind::Coulomb, 0.0, eps, k, opts);
}

inline WindowedSeparatedConvolution<double, 3>
WindowedCoulombOperator(World& world,
                        double eps,
                        const WindowedConvolutionOptions& opts) {
    return WindowedSeparatedConvolution<double, 3>(
        world, WindowedKernelKind::Coulomb, 0.0, eps, FunctionDefaults<3>::get_k(), opts);
}

inline WindowedSeparatedConvolution<double, 3>
WindowedCoulombOperator(World& world,
                        double eps,
                        bool enable_dynamic_refinement) {
    WindowedConvolutionOptions opts;
    opts.enable_dynamic_refinement = enable_dynamic_refinement;
    return WindowedSeparatedConvolution<double, 3>(
        world, WindowedKernelKind::Coulomb, 0.0, eps, FunctionDefaults<3>::get_k(), opts);
}

inline std::shared_ptr<WindowedSeparatedConvolution<double, 3>>
WindowedCoulombOperatorPtr(World& world,
                           double eps,
                           int k = FunctionDefaults<3>::get_k(),
                           const WindowedConvolutionOptions& opts = {}) {
    return std::make_shared<WindowedSeparatedConvolution<double, 3>>(
        world, WindowedKernelKind::Coulomb, 0.0, eps, k, opts);
}

inline std::shared_ptr<WindowedSeparatedConvolution<double, 3>>
WindowedCoulombOperatorPtr(World& world,
                           double eps,
                           const WindowedConvolutionOptions& opts) {
    return std::make_shared<WindowedSeparatedConvolution<double, 3>>(
        world, WindowedKernelKind::Coulomb, 0.0, eps, FunctionDefaults<3>::get_k(), opts);
}

inline std::shared_ptr<WindowedSeparatedConvolution<double, 3>>
WindowedCoulombOperatorPtr(World& world,
                           double eps,
                           bool enable_dynamic_refinement) {
    WindowedConvolutionOptions opts;
    opts.enable_dynamic_refinement = enable_dynamic_refinement;
    return std::make_shared<WindowedSeparatedConvolution<double, 3>>(
        world, WindowedKernelKind::Coulomb, 0.0, eps, FunctionDefaults<3>::get_k(), opts);
}

template <std::size_t NDIM = 3>
inline WindowedSeparatedConvolution<double, NDIM>
WindowedBSHOperator(World& world,
                    double mu,
                    double eps,
                    int k = FunctionDefaults<NDIM>::get_k(),
                    const WindowedConvolutionOptions& opts = {}) {
    static_assert(NDIM == 3, "WindowedBSHOperator currently supports NDIM == 3.");
    return WindowedSeparatedConvolution<double, NDIM>(
        world, WindowedKernelKind::BSH, mu, eps, k, opts);
}

template <std::size_t NDIM = 3>
inline WindowedSeparatedConvolution<double, NDIM>
WindowedBSHOperator(World& world,
                    double mu,
                    double eps,
                    const WindowedConvolutionOptions& opts) {
    static_assert(NDIM == 3, "WindowedBSHOperator currently supports NDIM == 3.");
    return WindowedSeparatedConvolution<double, NDIM>(
        world, WindowedKernelKind::BSH, mu, eps, FunctionDefaults<NDIM>::get_k(), opts);
}

template <std::size_t NDIM = 3>
inline WindowedSeparatedConvolution<double, NDIM>
WindowedBSHOperator(World& world,
                    double mu,
                    double eps,
                    bool enable_dynamic_refinement) {
    static_assert(NDIM == 3, "WindowedBSHOperator currently supports NDIM == 3.");
    WindowedConvolutionOptions opts;
    opts.enable_dynamic_refinement = enable_dynamic_refinement;
    return WindowedSeparatedConvolution<double, NDIM>(
        world, WindowedKernelKind::BSH, mu, eps, FunctionDefaults<NDIM>::get_k(), opts);
}

template <std::size_t NDIM = 3>
inline std::shared_ptr<WindowedSeparatedConvolution<double, NDIM>>
WindowedBSHOperatorPtr(World& world,
                       double mu,
                       double eps,
                       int k = FunctionDefaults<NDIM>::get_k(),
                       const WindowedConvolutionOptions& opts = {}) {
    static_assert(NDIM == 3, "WindowedBSHOperatorPtr currently supports NDIM == 3.");
    return std::make_shared<WindowedSeparatedConvolution<double, NDIM>>(
        world, WindowedKernelKind::BSH, mu, eps, k, opts);
}

template <std::size_t NDIM = 3>
inline std::shared_ptr<WindowedSeparatedConvolution<double, NDIM>>
WindowedBSHOperatorPtr(World& world,
                       double mu,
                       double eps,
                       const WindowedConvolutionOptions& opts) {
    static_assert(NDIM == 3, "WindowedBSHOperatorPtr currently supports NDIM == 3.");
    return std::make_shared<WindowedSeparatedConvolution<double, NDIM>>(
        world, WindowedKernelKind::BSH, mu, eps, FunctionDefaults<NDIM>::get_k(), opts);
}

template <std::size_t NDIM = 3>
inline std::shared_ptr<WindowedSeparatedConvolution<double, NDIM>>
WindowedBSHOperatorPtr(World& world,
                       double mu,
                       double eps,
                       bool enable_dynamic_refinement) {
    static_assert(NDIM == 3, "WindowedBSHOperatorPtr currently supports NDIM == 3.");
    WindowedConvolutionOptions opts;
    opts.enable_dynamic_refinement = enable_dynamic_refinement;
    return std::make_shared<WindowedSeparatedConvolution<double, NDIM>>(
        world, WindowedKernelKind::BSH, mu, eps, FunctionDefaults<NDIM>::get_k(), opts);
}

inline WindowedSeparatedConvolution<double, 3>
WindowedBSHOperator3D(World& world,
                      double mu,
                      double eps,
                      int k = FunctionDefaults<3>::get_k(),
                      const WindowedConvolutionOptions& opts = {}) {
    return WindowedBSHOperator<3>(world, mu, eps, k, opts);
}

inline WindowedSeparatedConvolution<double, 3>
WindowedBSHOperator3D(World& world,
                      double mu,
                      double eps,
                      const WindowedConvolutionOptions& opts) {
    return WindowedBSHOperator<3>(world, mu, eps, opts);
}

inline WindowedSeparatedConvolution<double, 3>
WindowedBSHOperator3D(World& world,
                      double mu,
                      double eps,
                      bool enable_dynamic_refinement) {
    return WindowedBSHOperator<3>(world, mu, eps, enable_dynamic_refinement);
}

inline std::shared_ptr<WindowedSeparatedConvolution<double, 3>>
WindowedBSHOperatorPtr3D(World& world,
                         double mu,
                         double eps,
                         int k = FunctionDefaults<3>::get_k(),
                         const WindowedConvolutionOptions& opts = {}) {
    return WindowedBSHOperatorPtr<3>(world, mu, eps, k, opts);
}

inline std::shared_ptr<WindowedSeparatedConvolution<double, 3>>
WindowedBSHOperatorPtr3D(World& world,
                         double mu,
                         double eps,
                         const WindowedConvolutionOptions& opts) {
    return WindowedBSHOperatorPtr<3>(world, mu, eps, opts);
}

inline std::shared_ptr<WindowedSeparatedConvolution<double, 3>>
WindowedBSHOperatorPtr3D(World& world,
                         double mu,
                         double eps,
                         bool enable_dynamic_refinement) {
    return WindowedBSHOperatorPtr<3>(world, mu, eps, enable_dynamic_refinement);
}

inline WindowedSeparatedConvolution<double, 3>
WindowedYukawaOperator(World& world,
                       double mu,
                       double eps,
                       int k = FunctionDefaults<3>::get_k(),
                       const WindowedConvolutionOptions& opts = {}) {
    return WindowedSeparatedConvolution<double, 3>(
        world, WindowedKernelKind::YukawaRaw, mu, eps, k, opts);
}

inline std::shared_ptr<WindowedSeparatedConvolution<double, 3>>
WindowedYukawaOperatorPtr(World& world,
                          double mu,
                          double eps,
                          int k = FunctionDefaults<3>::get_k(),
                          const WindowedConvolutionOptions& opts = {}) {
    return std::make_shared<WindowedSeparatedConvolution<double, 3>>(
        world, WindowedKernelKind::YukawaRaw, mu, eps, k, opts);
}

} // namespace madness

#endif // MADNESS_MRA_WINDOWED_CONVOLUTION_H__INCLUDED
