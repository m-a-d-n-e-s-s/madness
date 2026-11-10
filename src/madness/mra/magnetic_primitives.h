/// \file magnetic_primitives.h
/// \brief Building blocks shared by the magnetic Green's-function kernels (magnetic_kernel.h)
///
/// All kernels are proper-time (Laplace) quadratures over heat-kernel factors built from
///   - normalized Gaussian convolutions on one or two axes of a 3D function (low-dimensional apply),
///   - one-term separated Gaussian operators with anisotropic exponents (2D or 3D),
///   - pointwise multiplication by a bilinear phase e^{i kappa (x0-o0)(x1-o1)}, applied on the
///     values at the quadrature points after refining the tree where the phase is not resolved,
///   - a tau quadrature (Gauss-Legendre or trapezoid in u = ln tau).
/// See magnetic_kernels.md for the theory and the provenance of each piece.

#ifndef MADNESS_MRA_MAGNETIC_PRIMITIVES_H__INCLUDED
#define MADNESS_MRA_MAGNETIC_PRIMITIVES_H__INCLUDED

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/legendre.h>
#include <madness/constants.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

namespace madness {
namespace magnetic {

// =====================================================================================
// Coefficients of the harmonic / magnetic heat kernels
// =====================================================================================

/// (a, s) with e^{-tau T_perp} = e^{-a pi_y^2} e^{-s pi_x^2} e^{-a pi_y^2} (free limit for |B tau| < 1e-8)
inline std::pair<double, double> splitting_coefficients(const double B, const double tau) {
    const double x = B * tau;
    if (std::abs(x) < 1e-8) return {0.25 * tau, 0.5 * tau};
    return {std::tanh(0.5 * x) / (2 * B), std::sinh(x) / (2 * B)};
}

/// Mehler mask exponent m = (w/2) tanh(w tau/2) of the harmonic oscillator with frequency w:
/// <x|e^{-tau H_HO}|x'> = sqrt(c/pi) e^{-m x^2} e^{-c (x-x')^2} e^{-m x'^2} per dimension
inline double mehler_m(const double w, const double tau) { return 0.5 * w * std::tanh(0.5 * w * tau); }

/// log c with c = w / (2 sinh(w tau)), finite for w -> 0 (c -> 1/(2 tau)) and w tau -> infinity
inline double mehler_log_c(const double w, const double tau) {
    const double bt = std::abs(w) * tau;
    return -std::log(2 * tau) + (bt == 0.0 ? 0.0 : std::log(2 * bt) - bt - std::log(-std::expm1(-2.0 * bt)));
}

// =====================================================================================
// tau quadrature
// =====================================================================================

/// int_0^tmax dtau F(tau) ~ tmin F(0) + sum_i exp(log_dw_i) F(tau_i)
struct TauQuadrature {
    struct Node {
        double tau, log_dw;     ///< node and log of its weight in tau (including the Jacobian tau)
    };
    std::vector<Node> nodes;
    double tmin = 0.0;          ///< weight of the [0, tmin] part, applied to F(0) = f
};

/// Gauss-Legendre in u = ln tau on [ln tmin, ln tmax], ascending nodes
inline TauQuadrature gauss_legendre_log(const int n, const double tmin, const double tmax) {
    MADNESS_CHECK_THROW(n > 0 && tmin > 0 && tmax > tmin, "gauss_legendre_log: need n > 0, 0 < tmin < tmax");
    std::vector<double> x(n), w(n);
    MADNESS_CHECK_THROW(gauss_legendre(n, -1.0, 1.0, x.data(), w.data()), "gauss_legendre_log: no quadrature");
    std::vector<int> idx(n);
    for (int i = 0; i < n; ++i) idx[i] = i;
    std::sort(idx.begin(), idx.end(), [&x](int a, int b) { return x[a] < x[b]; });
    const double s0 = std::log(tmin), s1 = std::log(tmax);
    TauQuadrature q;
    q.tmin = tmin;
    for (int i : idx) {
        const double t = std::exp(0.5 * (s1 - s0) * x[i] + 0.5 * (s1 + s0));
        const double ws = 0.5 * (s1 - s0) * w[i];
        q.nodes.push_back({t, std::log(ws * t)});
    }
    return q;
}

/// trapezoid rule in u = ln tau with step h from tmin to tmax; the first node carries the
/// endpoint weight 1/2 (the last one does not: the integrand has decayed there)
inline TauQuadrature trapezoid_log(const double h, const double tmin, const double tmax) {
    MADNESS_CHECK_THROW(h > 0 && tmin > 0 && tmax > tmin, "trapezoid_log: need h > 0, 0 < tmin < tmax");
    const double umin = std::log(tmin), umax = std::log(tmax);
    const long npt = long(std::floor((umax - umin) / h + 0.5)) + 1;
    TauQuadrature q;
    q.tmin = std::exp(umin);
    for (long i = 0; i < npt; ++i) {
        const double tau = std::exp(umin + i * h);
        double lw = std::log(h * tau);
        if (i == 0) lw += std::log(0.5);
        q.nodes.push_back({tau, lw});
    }
    return q;
}

// =====================================================================================
// Low-dimensional defaults and Gaussian convolutions
// =====================================================================================

/// FunctionDefaults<1> and <2> with the cell, k and thresh of the cubic FunctionDefaults<3>
/// (required by the low-dimensional apply and by the 1D/2D functions of the FT kernel)
inline void set_lowdim_defaults_from_3d() {
    const Tensor<double> cell = FunctionDefaults<3>::get_cell();
    const double lo = cell(0, 0), hi = cell(0, 1);
    for (int d = 1; d < 3; ++d)
        MADNESS_CHECK_THROW(cell(d, 0) == lo && cell(d, 1) == hi, "magnetic kernels: the 3D cell must be cubic");
    FunctionDefaults<1>::set_cubic_cell(lo, hi);
    FunctionDefaults<2>::set_cubic_cell(lo, hi);
    FunctionDefaults<1>::set_k(FunctionDefaults<3>::get_k());
    FunctionDefaults<2>::set_k(FunctionDefaults<3>::get_k());
    FunctionDefaults<1>::set_thresh(FunctionDefaults<3>::get_thresh());
    FunctionDefaults<2>::set_thresh(FunctionDefaults<3>::get_thresh());
}

/// normalized 1D Gaussian (2 pi v)^{-1/2} exp(-x^2/(2 v)) in the simulation coordinates of an axis
/// of the given width (v in user coordinates)
inline std::shared_ptr<Convolution1D<double>> gauss_kernel_1d(const int k, const double width, const double variance) {
    MADNESS_CHECK_THROW(variance > 0.0, "gauss kernel: variance must be positive");
    const double expnt = 0.5 * width * width / variance;
    return GaussianConvolution1DCache<double>::get(k, expnt, 0, LatticeRange(false));
}

/// FunctionDefaults<ODIM> must describe axes [first, first+ODIM) of FunctionDefaults<NDIM>
template <std::size_t ODIM, typename T, std::size_t NDIM>
void check_lowdim_defaults(const Function<T, NDIM>& f, const int first, const char* name) {
    const Tensor<double>& cello = FunctionDefaults<ODIM>::get_cell();
    const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();
    for (std::size_t d = 0; d < ODIM; ++d)
        MADNESS_CHECK_THROW(cello(d, 0) == cell(first + d, 0) && cello(d, 1) == cell(first + d, 1),
                            (std::string(name) + ": FunctionDefaults" + std::to_string(ODIM) +
                             "D cell must equal the function's cell along the active axes").c_str());
    MADNESS_CHECK_THROW(FunctionDefaults<ODIM>::get_k() == f.k(),
                        (std::string(name) + ": FunctionDefaults" + std::to_string(ODIM) +
                         "D k must equal the function's k").c_str());
}

/// refine significant leaves next to refined neighbors until the tree is stable (at most npass passes)
template <typename T, std::size_t NDIM>
void broaden_until_stable(const Function<T, NDIM>& f, const int npass) {
    for (int pass = 0; pass < npass; ++pass) {
        const std::size_t size0 = f.size();
        f.broaden();
        if (f.size() == size0) break;
    }
}

/// convolution of a 3D function along one axis (0 or 2) with the normalized Gaussian of variance v
template <typename T>
Function<T, 3> gauss_conv_1d(const Function<T, 3>& f, const int axis, const double variance, const int nbroaden = 0) {
    MADNESS_CHECK_THROW(axis == 0 || axis == 2, "gauss_conv_1d: axis must be the first or the last axis");
    check_lowdim_defaults<1>(f, axis, "gauss_conv_1d");
    const double width = FunctionDefaults<3>::get_cell_width()[axis];
    SeparatedConvolution<double, 1> op(f.world(),
        std::vector<ConvolutionND<double, 1>>{ConvolutionND<double, 1>(gauss_kernel_1d(f.k(), width, variance))}, f.k());
    op.set_particle(axis == 0 ? 1 : 2);
    broaden_until_stable(f, nbroaden);
    return madness::apply(op, f);
}

/// convolution of a 3D function along axes (first, first+1), first = 0 or 1, with the product of
/// normalized Gaussians of variances v0, v1; the remaining axis is a spectator
template <typename T>
Function<T, 3> gauss_conv_2d(const Function<T, 3>& f, const int first, const double v0, const double v1,
                             const int nbroaden = 0) {
    MADNESS_CHECK_THROW(first == 0 || first == 1, "gauss_conv_2d: first must be 0 or 1");
    check_lowdim_defaults<2>(f, first, "gauss_conv_2d");
    const auto& width = FunctionDefaults<3>::get_cell_width();
    ConvolutionND<double, 2> conv;
    conv.setop(0, gauss_kernel_1d(f.k(), width[first], v0));
    conv.setop(1, gauss_kernel_1d(f.k(), width[first + 1], v1));
    SeparatedConvolution<double, 2> op(f.world(), std::vector<ConvolutionND<double, 2>>{conv}, f.k());
    op.set_particle(first == 0 ? 1 : 2);
    broaden_until_stable(f, nbroaden);
    return madness::apply(op, f);
}

/// one-term separated operator prod_d exp(-a_d x_d^2) (unnormalized, user coordinates) on the
/// FunctionDefaults<NDIM> cell, a_d = exp(log_a[d]); an exponent below 1e-10 in simulation coordinates
/// (it underflows at large tau) is raised to it: the kernel is then 1 on the cell to ~1e-10
template <std::size_t NDIM>
std::shared_ptr<SeparatedConvolution<double, NDIM>> gauss_product_operator(World& world,
                                                                          const std::array<double, NDIM>& log_a) {
    const Tensor<double>& width = FunctionDefaults<NDIM>::get_cell_width();
    const int k = FunctionDefaults<NDIM>::get_k();
    const auto lattice = FunctionDefaults<NDIM>::get_bc().lattice_range();
    ConvolutionND<double, NDIM> term;
    double fac = 1.0;
    for (std::size_t d = 0; d < NDIM; ++d) {
        const double w2 = width[d] * width[d];
        const double a = std::exp(std::max(log_a[d] + std::log(w2), std::log(1e-10)));
        term.setop(d, GaussianConvolution1DCache<double>::get(k, a, 0, lattice[d]));
        fac *= std::sqrt(constants::pi / a) * width[d];
    }
    term.setfac(fac);
    return std::make_shared<SeparatedConvolution<double, NDIM>>(world, std::vector<ConvolutionND<double, NDIM>>{term}, k);
}

// =====================================================================================
// Bilinear phase e^{i kappa (x0-o0)(x1-o1)} on axes 0 and 1
// =====================================================================================

/// how the tree is refined before the phase is multiplied onto the values
enum class PhaseRefinement {
    poly_error,     ///< refine where the best degree-(k-1) fit of the phase misses tolfac*truncate_tol
                    ///  relative to the box norm (the chain kernels)
    box_phase       ///< refine boxes within rho < rext where the phase varies by more than dphi and the
                    ///  coefficient norm exceeds tol (the FT kernel)
};

struct BilinearPhase {
    double kappa = 0.0, o0 = 0.0, o1 = 0.0;   ///< theta = kappa (x0-o0)(x1-o1)
    PhaseRefinement refinement = PhaseRefinement::poly_error;
    double tolfac = 0.1;                      ///< poly_error: tolerance factor
    bool recursive = true;                    ///< poly_error: one recursive sweep (single process) or passes
    double dphi = 3.0, rext = 1e10, tol = 0.0;  ///< box_phase
    double rcut = 0.0, rtaper = 0.0;          ///< taper 0.5 erfc((rho-rcut)/rtaper) of the factor; rcut <= 0: none

    BilinearPhase() = default;
    BilinearPhase(double kappa, double o0, double o1) : kappa(kappa), o0(o0), o1(o1) {}
    BilinearPhase with_kappa(double k) const { BilinearPhase p(*this); p.kappa = k; return p; }
};

namespace detail {

/// refinement test of PhaseRefinement::poly_error: the best degree-(k-1) approximation of e^{i theta}
/// with theta varying by 2w over a box has error ~ w^k / (2^(k-1) k!)
template <typename T, std::size_t NDIM>
struct PolyErrorRefineOp {
    double kappa = 0.0, o0 = 0.0, o1 = 0.0, tolfac = 0.1;

    bool operator()(const FunctionImpl<T, NDIM>* impl, const Key<NDIM>& key, const FunctionNode<T, NDIM>& node) const {
        const double tol = tolfac * impl->truncate_tol(impl->get_thresh(), key);
        const double cnorm = node.coeff().normf();
        if (cnorm < tol) return false;
        const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();
        const auto& width = FunctionDefaults<NDIM>::get_cell_width();
        const double h = std::pow(0.5, key.level());
        double tmin = 1e300, tmax = -1e300;
        for (int c0 = 0; c0 < 2; ++c0) {
            for (int c1 = 0; c1 < 2; ++c1) {
                const double x0 = cell(0, 0) + width[0] * h * (key.translation()[0] + c0) - o0;
                const double x1 = cell(1, 0) + width[1] * h * (key.translation()[1] + c1) - o1;
                tmin = std::min(tmin, kappa * x0 * x1);
                tmax = std::max(tmax, kappa * x0 * x1);
            }
        }
        // the bilinear phase takes its extrema over a box at the corners
        const double w = 0.5 * (tmax - tmin);
        const int k = impl->get_k();
        double err = std::pow(0.5, k - 1);
        for (int i = 1; i <= k; ++i) err *= w / i;
        return cnorm * err > tol;
    }
    template <typename Archive> void serialize(Archive& ar) { ar & kappa & o0 & o1 & tolfac; }
};

/// refinement test of PhaseRefinement::box_phase
template <typename T, std::size_t NDIM>
struct BoxPhaseRefineOp {
    double kappa = 0, o0 = 0, o1 = 0, rext = 0, dphi = 0, tol = 0;

    template <typename implT, typename nodeT>
    bool operator()(implT*, const Key<NDIM>& key, const nodeT& t) const {
        const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();
        const double width = FunctionDefaults<NDIM>::get_cell_width()[0];
        const double h = width / double(1L << key.level());
        const double x0 = cell(0, 0) + key.translation()[0] * h - o0, y0 = cell(1, 0) + key.translation()[1] * h - o1;
        auto nearest = [h](double a) { return (a <= 0.0 && a + h >= 0.0) ? 0.0 : std::min(std::abs(a), std::abs(a + h)); };
        const double rx = nearest(x0), ry = nearest(y0);
        if (rx * rx + ry * ry > rext * rext) return false;
        const double mx = std::max(std::abs(x0), std::abs(x0 + h)), my = std::max(std::abs(y0), std::abs(y0 + h));
        if (std::abs(kappa) * (mx + my) * h <= dphi) return false;
        return t.coeff().normf() > tol;
    }
    template <typename Archive> void serialize(Archive& ar) { ar & kappa & o0 & o1 & rext & dphi & tol; }
};

/// values at the quadrature points times [taper(rho)] e^{i theta}
template <std::size_t NDIM>
struct PhaseValueOp {
    double kappa = 0.0, o0 = 0.0, o1 = 0.0, rcut = 0.0, rtaper = 0.0;

    void operator()(const Key<NDIM>& key, Tensor<double_complex>& values) const {
        const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();
        const auto& width = FunctionDefaults<NDIM>::get_cell_width();
        const int k = values.dim(0);
        const Tensor<double>& qx = FunctionCommonData<double, NDIM>::get(k).quad_x;
        const double h = std::pow(0.5, key.level());
        std::vector<double> x0(k), x1(k);
        for (int i = 0; i < k; ++i) {
            x0[i] = cell(0, 0) + width[0] * h * (key.translation()[0] + qx(i)) - o0;
            x1[i] = cell(1, 0) + width[1] * h * (key.translation()[1] + qx(i)) - o1;
        }
        const long rest = values.size() / (long(k) * k);    // remaining dimensions, contiguous
        double_complex* v = values.ptr();
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                const double th = kappa * x0[i] * x1[j];
                double_complex fac(std::cos(th), std::sin(th));
                if (rcut > 0.0) fac *= 0.5 * std::erfc((std::sqrt(x0[i] * x0[i] + x1[j] * x1[j]) - rcut) / rtaper);
                double_complex* p = v + (long(i) * k + j) * rest;
                for (long r = 0; r < rest; ++r) p[r] *= fac;
            }
        }
    }
    template <typename Archive> void serialize(Archive& ar) { ar & kappa & o0 & o1 & rcut & rtaper; }
};

/// out-of-place form of PhaseValueOp for madness::unary_op (values in, new values out)
template <std::size_t NDIM>
struct PhaseValueOutOp {
    typedef double_complex resultT;
    PhaseValueOp<NDIM> op;
    Tensor<double_complex> operator()(const Key<NDIM>& key, const Tensor<double_complex>& values) const {
        Tensor<double_complex> out = copy(values);
        op(key, out);
        return out;
    }
    template <typename Archive> void serialize(Archive& ar) { ar & op; }
};

/// one node of the recursive refinement: split key with FunctionImpl::refine_op and give each new
/// child its own task (one sweep without fences between levels); key must be a leaf
template <typename T, std::size_t NDIM, typename opT>
void refine_node_recursive(World* world, FunctionImpl<T, NDIM>* impl, const opT& op, const Key<NDIM>& key) {
    impl->template refine_op<opT>(op, key);
    bool split;
    {
        typename FunctionImpl<T, NDIM>::dcT::const_accessor acc;
        const bool found = impl->get_coeffs().find(acc, key);
        MADNESS_CHECK(found);
        split = acc->second.has_children();
    }
    if (split)
        for (KeyChildIterator<NDIM> kit(key); kit; ++kit)
            world->taskq.add(&refine_node_recursive<T, NDIM, opT>, world, impl, op, kit.key());
}

/// refine the reconstructed g until op holds nowhere: one recursive sweep (single-process worlds), or
/// refine_general passes (at most 30)
template <typename T, std::size_t NDIM, typename opT>
void refine_recursive(const Function<T, NDIM>& g, const opT& op, const bool recursive) {
    World& world = g.world();
    if (!recursive || world.size() > 1) {
        for (int pass = 0; pass < 30; ++pass) {
            const std::size_t size0 = g.size();
            g.refine_general(op);
            if (g.size() == size0) break;
        }
        return;
    }
    if (!g.is_reconstructed()) g.reconstruct();
    FunctionImpl<T, NDIM>* impl = g.get_impl().get();
    std::vector<Key<NDIM>> leaves;
    for (auto it = impl->get_coeffs().begin(); it != impl->get_coeffs().end(); ++it)
        if (it->second.has_coeff() && !it->second.has_children()) leaves.push_back(it->first);
    for (const auto& key : leaves) world.taskq.add(&refine_node_recursive<T, NDIM, opT>, &world, impl, op, key);
    world.gop.fence();
}

template <std::size_t NDIM>
PhaseValueOp<NDIM> value_op(const BilinearPhase& p) {
    PhaseValueOp<NDIM> vop;
    vop.kappa = p.kappa; vop.o0 = p.o0; vop.o1 = p.o1; vop.rcut = p.rcut; vop.rtaper = p.rtaper;
    return vop;
}

}  // namespace detail

/// refine the tree of f (same values) where the phase p is not resolved, according to p.refinement
template <std::size_t NDIM>
void refine_for_phase(const Function<double_complex, NDIM>& f, const BilinearPhase& p) {
    f.reconstruct();
    if (p.refinement == PhaseRefinement::poly_error) {
        MADNESS_CHECK_THROW(p.tolfac > 0.0, "refine_for_phase: tolfac must be positive");
        detail::PolyErrorRefineOp<double_complex, NDIM> rop;
        rop.kappa = p.kappa; rop.o0 = p.o0; rop.o1 = p.o1; rop.tolfac = p.tolfac;
        detail::refine_recursive(f, rop, p.recursive);
    } else {
        detail::BoxPhaseRefineOp<double_complex, NDIM> rop;
        rop.kappa = p.kappa; rop.o0 = p.o0; rop.o1 = p.o1; rop.rext = p.rext; rop.dphi = p.dphi; rop.tol = p.tol;
        for (int pass = 0; pass < 12; ++pass) {
            const std::size_t before = f.tree_size();
            f.refine_general(rop);
            if (f.tree_size() == before) break;
        }
    }
}

/// f e^{i theta} (times the taper, if any), out of place and untruncated; the tree of f is refined
/// (same values)
template <std::size_t NDIM>
Function<double_complex, NDIM> mul_bilinear_phase(const Function<double_complex, NDIM>& f, const BilinearPhase& p) {
    refine_for_phase(f, p);
    detail::PhaseValueOutOp<NDIM> vop;
    vop.op = detail::value_op<NDIM>(p);
    return madness::unary_op(f, vop);
}

/// f <- f e^{i theta} (times the taper, if any) in place, then truncated
template <std::size_t NDIM>
void mul_bilinear_phase_inplace(Function<double_complex, NDIM>& f, const BilinearPhase& p) {
    refine_for_phase(f, p);
    f.unaryop(detail::value_op<NDIM>(p));
    f.truncate();
}

}  // namespace magnetic
}  // namespace madness

#endif  // MADNESS_MRA_MAGNETIC_PRIMITIVES_H__INCLUDED
