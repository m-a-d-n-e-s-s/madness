/*
  pymadness - Python bindings for MADNESS
  Operators: SeparatedConvolution (BSH, Coulomb, Slater), Derivative.

  SeparatedConvolution and Derivative inherit from WorldObject which is
  not copyable/movable. We use shared_ptr holders and pointer factory
  functions throughout.
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/derivative.h>

namespace py = pybind11;
using namespace madness;

using SepConv3D = SeparatedConvolution<double, 3>;
using SepConv2D = SeparatedConvolution<double, 2>;
using SepConv1D = SeparatedConvolution<double, 1>;

/// Refine significant leaves next to refined neighbors until the tree is stable (at most npass passes).

/// A kernel narrower than the boxes couples neighbors through edge terms of size
/// ~f(edge)*sigma, resolved only at the source level + 1; they cancel only if both sides of
/// a box face are at the same level.  Broadening leaves the function itself unchanged.
template <typename T, std::size_t NDIM>
static void broaden_until_stable(const Function<T, NDIM>& f, const int npass) {
    for (int pass = 0; pass < npass; ++pass) {
        const std::size_t size0 = f.size();
        f.broaden();
        if (f.size() == size0) break;
    }
}

/// Normalized 1D Gaussian kernel (2 pi v)^{-1/2} exp(-x^2/(2 v)) in the simulation coordinates
/// of an axis of the given width (user coordinates)
static std::shared_ptr<Convolution1D<double>> gauss_kernel_1d(const int k, const double width,
                                                              const double variance) {
    MADNESS_CHECK_THROW(variance > 0.0, "gauss kernel: variance must be positive");
    // exponent in simulation coordinates; the cached kernel has unit integral
    const double expnt = 0.5 * width * width / variance;
    return GaussianConvolution1DCache<double>::get(k, expnt, 0, LatticeRange(false));
}

/// FunctionDefaults<ODIM> must describe axes [first, first+ODIM) of FunctionDefaults<NDIM>
template <std::size_t ODIM, typename T, std::size_t NDIM>
static void check_lowdim_defaults(const Function<T, NDIM>& f, const int first, const char* name) {
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

/// Convolution of a 2D or 3D function along one axis with the normalized Gaussian of the given
/// variance (user coordinates), (2 pi v)^{-1/2} exp(-x^2/(2 v)); the other axes are spectators.

/// A 1D SeparatedConvolution acting on particle 1 (axis 0) or 2 (the last axis) of the
/// function, i.e. the low-dimensional apply used for the exchange operator in 6D.  For a 2D
/// function both axes are reachable; for 3D only the outer ones (axis 1 is not a contiguous
/// block at either end).  FunctionDefaults<1> must describe the same interval and k as the
/// active axis.  T = double or double_complex (the kernel is real).
template <typename T, std::size_t NDIM>
static Function<T, NDIM> apply_gauss_1d(const Function<T, NDIM>& f,
                                             const int axis, const double variance,
                                             const int nbroaden) {
    MADNESS_CHECK_THROW(axis == 0 || axis == int(NDIM) - 1,
                        "apply_gauss_1d: axis must be the first or the last axis");
    check_lowdim_defaults<1>(f, axis, "apply_gauss_1d");
    const double width = FunctionDefaults<NDIM>::get_cell_width()[axis];
    SeparatedConvolution<double, 1> op(f.world(),
        std::vector<ConvolutionND<double, 1>>{
            ConvolutionND<double, 1>(gauss_kernel_1d(f.k(), width, variance))}, f.k());
    op.set_particle(axis == 0 ? 1 : 2);
    broaden_until_stable(f, nbroaden);
    return madness::apply(op, f);
}

/// Convolution of a 3D function along two adjacent axes (first, first+1) with the product of
/// normalized Gaussians of variances v0, v1 (user coordinates); the remaining axis is a spectator.

/// A 2D SeparatedConvolution on particle 1 (first = 0) or particle 2 (first = 1).
/// FunctionDefaults<2> must describe the same intervals and k as the active axes.
template <typename T>
static Function<T, 3> apply_gauss_2d(const Function<T, 3>& f, const int first,
                                          const double v0, const double v1, const int nbroaden) {
    MADNESS_CHECK_THROW(first == 0 || first == 1, "apply_gauss_2d: first must be 0 or 1");
    check_lowdim_defaults<2>(f, first, "apply_gauss_2d");
    const auto& width = FunctionDefaults<3>::get_cell_width();
    ConvolutionND<double, 2> conv;
    conv.setop(0, gauss_kernel_1d(f.k(), width[first], v0));
    conv.setop(1, gauss_kernel_1d(f.k(), width[first + 1], v1));
    SeparatedConvolution<double, 2> op(f.world(), std::vector<ConvolutionND<double, 2>>{conv}, f.k());
    op.set_particle(first == 0 ? 1 : 2);
    broaden_until_stable(f, nbroaden);
    return madness::apply(op, f);
}

/// Refinement test for mul_phase: refine a box if the phase theta = kappa (x0-o0)(x1-o1) is not
/// resolved by the polynomials of the box to the truncation tolerance, relative to the box norm.

/// The best degree-(k-1) approximation of e^{i theta} with theta varying by 2w over an interval
/// has error ~ w^k / (2^(k-1) k!); refine if that times the coefficient norm exceeds
/// tolfac * truncate_tol.
template <typename T, std::size_t NDIM>
struct PhaseRefineOp {
    double kappa = 0.0, o0 = 0.0, o1 = 0.0, tolfac = 0.1;

    bool operator()(const FunctionImpl<T, NDIM>* impl, const Key<NDIM>& key,
                    const FunctionNode<T, NDIM>& node) const {
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
        for (int i = 1; i <= k; ++i) err *= w / i;     // w^k / (2^(k-1) k!)
        return cnorm * err > tol;
    }
    template <typename Archive> void serialize(Archive& ar) { ar & kappa & o0 & o1 & tolfac; }
};

/// One node of refine_phase (recursive mode): split key with the core FunctionImpl::refine_op
/// (its own test, the impl's max_refine_level) and give each new child its own task, so the
/// whole refinement is one sweep without fences between levels.  key must be a leaf.
template <typename T, std::size_t NDIM, typename opT>
static void refine_node_recursive(World* world, FunctionImpl<T, NDIM>* impl, const opT& op,
                                  const Key<NDIM>& key) {
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

/// Refine the reconstructed g until op holds nowhere, in one sweep: a task per leaf, and a task
/// per child as soon as a box is split; one fence at the end.  refine_general splits flagged
/// leaves by one level per call (a walk from the root and a fence each), so the phase step needed
/// one call per new level plus a final one, each with two global size() reductions.
/// Preconditions: a conforming reconstructed tree, and an op that reads only the box it tests
/// (true for PhaseRefineOp); then the result is the multi-pass fixed point (unless more than 30
/// new levels are needed, where the passes stop early).  Single-process worlds only; otherwise,
/// and with recursive = false, the passes are used.
template <typename T, std::size_t NDIM, typename opT>
static void refine_phase(const Function<T, NDIM>& g, const opT& op, const bool recursive) {
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
    for (const auto& key : leaves)
        world.taskq.add(&refine_node_recursive<T, NDIM, opT>, &world, impl, op, key);
    world.gop.fence();
}

/// Multiply function values at the quadrature points by cos(theta) or sin(theta),
/// theta = kappa (x0-o0)(x1-o1)
template <std::size_t NDIM>
struct PhaseValueOp {
    double kappa = 0.0, o0 = 0.0, o1 = 0.0;
    bool use_sin = false;

    void operator()(const Key<NDIM>& key, Tensor<double>& values) const {
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
        double* v = values.ptr();
        for (int i = 0; i < k; ++i) {
            for (int j = 0; j < k; ++j) {
                const double th = kappa * x0[i] * x1[j];
                const double fac = use_sin ? std::sin(th) : std::cos(th);
                double* p = v + (long(i) * k + j) * rest;
                for (long r = 0; r < rest; ++r) p[r] *= fac;
            }
        }
    }
    template <typename Archive> void serialize(Archive& ar) { ar & kappa & o0 & o1 & use_sin; }
};

/// cos or sin of the phase on values, out-of-place form for madness::unary_op
template <std::size_t NDIM>
struct PhaseValueOutOp {
    typedef double resultT;
    PhaseValueOp<NDIM> op;
    Tensor<double> operator()(const Key<NDIM>& key, const Tensor<double>& values) const {
        Tensor<double> out = copy(values);
        op(key, out);
        return out;
    }
    template <typename Archive> void serialize(Archive& ar) { ar & op; }
};

/// (f cos theta, f sin theta) with theta = kappa (x0-o0)(x1-o1), without projecting the phase

/// The tree of f is refined (one recursive sweep, see refine_phase) where the phase is not resolved
/// to tolfac * truncate_tol (see PhaseRefineOp); f itself is refined (same values, finer tree).
/// The products are formed from the values at the quadrature points of the leaves by the
/// out-of-place madness::unary_op (2-3x faster than copy + in-place unaryop).  Truncate the
/// results afterwards.
template <std::size_t NDIM>
static std::pair<Function<double, NDIM>, Function<double, NDIM>>
mul_phase(const Function<double, NDIM>& f, const double kappa, const double o0, const double o1,
          const double tolfac, const bool recursive = true) {
    MADNESS_CHECK_THROW(tolfac > 0.0, "mul_phase: tolfac must be positive");
    f.reconstruct();
    PhaseRefineOp<double, NDIM> rop;
    rop.kappa = kappa; rop.o0 = o0; rop.o1 = o1; rop.tolfac = tolfac;
    refine_phase(f, rop, recursive);                 // f keeps the refined tree (same values)
    PhaseValueOutOp<NDIM> vop;
    vop.op.kappa = kappa; vop.op.o0 = o0; vop.op.o1 = o1;
    vop.op.use_sin = false;
    Function<double, NDIM> fc = madness::unary_op(f, vop);
    vop.op.use_sin = true;
    Function<double, NDIM> fs = madness::unary_op(f, vop);
    return {fc, fs};
}

/// Multiply complex function values at the quadrature points by e^{i theta},
/// theta = kappa (x0-o0)(x1-o1)
template <std::size_t NDIM>
struct PhaseExpiOp {
    double kappa = 0.0, o0 = 0.0, o1 = 0.0;

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
                const double_complex fac(std::cos(th), std::sin(th));
                double_complex* p = v + (long(i) * k + j) * rest;
                for (long r = 0; r < rest; ++r) p[r] *= fac;
            }
        }
    }
    template <typename Archive> void serialize(Archive& ar) { ar & kappa & o0 & o1; }
};

/// e^{i theta} on values in the out-of-place form expected by madness::unary_op
/// (coeff_value_adaptor: values in, new values out)
template <std::size_t NDIM>
struct PhaseExpiValueOp {
    typedef double_complex resultT;
    PhaseExpiOp<NDIM> op;
    Tensor<double_complex> operator()(const Key<NDIM>& key, const Tensor<double_complex>& values) const {
        Tensor<double_complex> out = copy(values);
        op(key, out);
        return out;
    }
    template <typename Archive> void serialize(Archive& ar) { ar & op; }
};

/// f e^{i theta} for a complex f, theta = kappa (x0-o0)(x1-o1), without projecting the phase:
/// the refinement of mul_phase on f itself (same values, finer tree), then the out-of-place
/// madness::unary_op on the values.  The result is untruncated.
template <std::size_t NDIM>
static Function<double_complex, NDIM>
mul_expi(const Function<double_complex, NDIM>& f, const double kappa, const double o0, const double o1,
         const double tolfac, const bool recursive = true) {
    MADNESS_CHECK_THROW(tolfac > 0.0, "mul_expi: tolfac must be positive");
    f.reconstruct();
    PhaseRefineOp<double_complex, NDIM> rop;
    rop.kappa = kappa; rop.o0 = o0; rop.o1 = o1; rop.tolfac = tolfac;
    refine_phase(f, rop, recursive);                 // f keeps the refined tree (same values)
    PhaseExpiValueOp<NDIM> vop;
    vop.op.kappa = kappa; vop.op.o0 = o0; vop.op.o1 = o1;
    return madness::unary_op(f, vop);
}

/// Test variants of f e^{i theta} (comparison only):
///   mode 0: mul_expi (copy f, refine the copy, in-place value op on the copy);
///   mode 1: refine f itself, then out-of-place madness::unary_op (unaryXXvalues) into a new tree;
///           f keeps the refined tree (same values);
///   mode 2: refine f itself, then the in-place value op on f: f is overwritten by the result.
template <std::size_t NDIM>
static Function<double_complex, NDIM>
mul_expi_variant(const Function<double_complex, NDIM>& f, const double kappa, const double o0,
                 const double o1, const double tolfac, const int mode) {
    if (mode == 0) {                                // the previous mul_expi: copy, refine, in place
        Function<double_complex, NDIM> g = copy(f);
        g.reconstruct();
        PhaseRefineOp<double_complex, NDIM> rop;
        rop.kappa = kappa; rop.o0 = o0; rop.o1 = o1; rop.tolfac = tolfac;
        refine_phase(g, rop, true);
        PhaseExpiOp<NDIM> vop;
        vop.kappa = kappa; vop.o0 = o0; vop.o1 = o1;
        g.unaryop(vop);
        return g;
    }
    MADNESS_CHECK_THROW(mode == 1 || mode == 2, "mul_expi_variant: mode must be 0, 1 or 2");
    f.reconstruct();
    PhaseRefineOp<double_complex, NDIM> rop;
    rop.kappa = kappa; rop.o0 = o0; rop.o1 = o1; rop.tolfac = tolfac;
    refine_phase(f, rop, true);
    PhaseExpiOp<NDIM> vop;
    vop.kappa = kappa; vop.o0 = o0; vop.o1 = o1;
    if (mode == 1) {
        PhaseExpiValueOp<NDIM> v;
        v.op = vop;
        return madness::unary_op(f, v);
    }
    Function<double_complex, NDIM> g = f;          // shallow: same impl
    g.unaryop(vop);
    return g;
}

void bind_operators(py::module_& m) {

    // =====================================================================
    // SeparatedConvolution — held via shared_ptr (WorldObject not copyable)
    // =====================================================================
    py::class_<SepConv3D, std::shared_ptr<SepConv3D>>(m, "SeparatedConvolution3D",
        "Integral operator represented as a separated convolution in 3D");

    py::class_<SepConv2D, std::shared_ptr<SepConv2D>>(m, "SeparatedConvolution2D",
        "Integral operator represented as a separated convolution in 2D");

    py::class_<SepConv1D, std::shared_ptr<SepConv1D>>(m, "SeparatedConvolution1D",
        "Integral operator represented as a separated convolution in 1D");

    // =====================================================================
    // Operator factory functions — use Ptr variants that return raw pointers,
    // wrap into shared_ptr for Python ownership.
    // =====================================================================

    // Coulomb operator: 1/(4*pi*r)
    m.def("CoulombOperator", [](World& world, double lo, double eps)
            -> std::shared_ptr<SepConv3D> {
        py::gil_scoped_release release;
        return std::shared_ptr<SepConv3D>(
            madness::CoulombOperatorPtr(world, lo, eps));
    }, py::arg("world"), py::arg("lo"), py::arg("eps"),
       "Create Coulomb operator 1/(4*pi*r).\n"
       "  lo: smallest length scale to resolve\n"
       "  eps: precision");

    // BSH operator: exp(-mu*r)/(4*pi*r)
    m.def("BSHOperator3D", [](World& world, double mu, double lo, double eps)
            -> std::shared_ptr<SepConv3D> {
        py::gil_scoped_release release;
        return std::shared_ptr<SepConv3D>(
            madness::BSHOperatorPtr3D(world, mu, lo, eps));
    }, py::arg("world"), py::arg("mu"), py::arg("lo"), py::arg("eps"),
       "Create BSH (bound-state Helmholtz) operator exp(-mu*r)/(4*pi*r).\n"
       "  mu: exponent (sqrt(-2*E) for bound states)\n"
       "  lo: smallest length scale\n"
       "  eps: precision");

    // BSH operator for 1D
    m.def("BSHOperator1D", [](World& world, double mu, double lo, double eps)
            -> std::shared_ptr<SepConv1D> {
        py::gil_scoped_release release;
        return std::shared_ptr<SepConv1D>(
            madness::BSHOperatorPtr<1>(world, mu, lo, eps));
    }, py::arg("world"), py::arg("mu"), py::arg("lo"), py::arg("eps"),
       "Create 1D BSH operator");

    // BSH operator for 2D: K_0(mu*r)/(2*pi), i.e. (-Laplacian + mu^2)^{-1}
    m.def("BSHOperator2D", [](World& world, double mu, double lo, double eps)
            -> std::shared_ptr<SepConv2D> {
        py::gil_scoped_release release;
        return std::shared_ptr<SepConv2D>(
            madness::BSHOperatorPtr<2>(world, mu, lo, eps));
    }, py::arg("world"), py::arg("mu"), py::arg("lo"), py::arg("eps"),
       "Create 2D BSH operator");

    // Gauss operator: exp(-mu*r*r)  (unnormalized)
    m.def("GaussOperator3D", [](World& world, double mu)
            -> std::shared_ptr<SepConv3D> {
        py::gil_scoped_release release;
        return std::shared_ptr<SepConv3D>(
            madness::GaussOperatorPtr<3>(world, mu));
    }, py::arg("world"), py::arg("mu"),
       "Create Gaussian convolution operator exp(-mu*r^2) in 3D (unnormalized).\n"
       "  mu: Gaussian exponent");

    m.def("GaussOperator2D", [](World& world, double mu)
            -> std::shared_ptr<SepConv2D> {
        py::gil_scoped_release release;
        return std::shared_ptr<SepConv2D>(
            madness::GaussOperatorPtr<2>(world, mu));
    }, py::arg("world"), py::arg("mu"),
       "Create Gaussian convolution operator exp(-mu*r^2) in 2D (unnormalized).\n"
       "  mu: Gaussian exponent");

    m.def("GaussOperator1D", [](World& world, double mu)
            -> std::shared_ptr<SepConv1D> {
        py::gil_scoped_release release;
        return std::shared_ptr<SepConv1D>(
            madness::GaussOperatorPtr<1>(world, mu));
    }, py::arg("world"), py::arg("mu"),
       "Create Gaussian convolution operator exp(-mu*x^2) in 1D (unnormalized).\n"
       "  mu: Gaussian exponent");

    // Slater operator: exp(-mu*r)
    m.def("SlaterOperator", [](World& world, double mu, double lo, double eps)
            -> std::shared_ptr<SepConv3D> {
        py::gil_scoped_release release;
        return std::shared_ptr<SepConv3D>(
            madness::SlaterOperatorPtr(world, mu, lo, eps));
    }, py::arg("world"), py::arg("mu"), py::arg("lo"), py::arg("eps"),
       "Create Slater operator exp(-mu*r)");

    // =====================================================================
    // apply() — apply operator to function
    // =====================================================================

    // apply(op, f) for 3D real
    m.def("apply", [](const SepConv3D& op, const Function<double, 3>& f) {
        py::gil_scoped_release release;
        return madness::apply(op, f);
    }, py::arg("op"), py::arg("f"),
       "Apply integral operator to function: result = op(f)");

    // apply(op, f) for 1D real
    m.def("apply", [](const SepConv1D& op, const Function<double, 1>& f) {
        py::gil_scoped_release release;
        return madness::apply(op, f);
    }, py::arg("op"), py::arg("f"),
       "Apply 1D integral operator to function");

    // apply(op, f) for 2D real
    m.def("apply", [](const SepConv2D& op, const Function<double, 2>& f) {
        py::gil_scoped_release release;
        return madness::apply(op, f);
    }, py::arg("op"), py::arg("f"),
       "Apply 2D integral operator to function");

    // 1D Gaussian convolution along one axis of a 2D function (other axis spectator)
    m.def("apply_gauss_1d", [](const Function<double, 2>& f, int axis, double variance, int nbroaden) {
        py::gil_scoped_release release;
        return apply_gauss_1d(f, axis, variance, nbroaden);
    }, py::arg("f"), py::arg("axis"), py::arg("variance"), py::arg("nbroaden") = 6,
       "Convolve a 2D function along axis 0 or 1 with the normalized Gaussian of the given\n"
       "variance; requires FunctionDefaults1D cell and k to match the 2D axis");
    m.def("apply_gauss_1d", [](const Function<double, 3>& f, int axis, double variance, int nbroaden) {
        py::gil_scoped_release release;
        return apply_gauss_1d(f, axis, variance, nbroaden);
    }, py::arg("f"), py::arg("axis"), py::arg("variance"), py::arg("nbroaden") = 6,
       "Convolve a 3D function along axis 0 or 2 with the normalized Gaussian of the given\n"
       "variance; requires FunctionDefaults1D cell and k to match the 3D axis");

    // multiplication by cos/sin of a bilinear phase
    m.def("mul_phase", [](const Function<double, 2>& f, double kappa, double o0, double o1,
                          double tolfac, bool recursive) {
        py::gil_scoped_release release;
        return mul_phase(f, kappa, o0, o1, tolfac, recursive);
    }, py::arg("f"), py::arg("kappa"), py::arg("o0") = 0.0, py::arg("o1") = 0.0,
       py::arg("tolfac") = 0.1, py::arg("recursive") = true,
       "(f cos theta, f sin theta), theta = kappa (x0-o0)(x1-o1); refines f where the phase is not\n"
       "resolved to tolfac * truncate_tol; truncate the results afterwards");
    m.def("mul_phase", [](const Function<double, 3>& f, double kappa, double o0, double o1,
                          double tolfac, bool recursive) {
        py::gil_scoped_release release;
        return mul_phase(f, kappa, o0, o1, tolfac, recursive);
    }, py::arg("f"), py::arg("kappa"), py::arg("o0") = 0.0, py::arg("o1") = 0.0,
       py::arg("tolfac") = 0.1, py::arg("recursive") = true,
       "(f cos theta, f sin theta), theta = kappa (x0-o0)(x1-o1); refines f where the phase is not\n"
       "resolved to tolfac * truncate_tol; truncate the results afterwards");

    // screened pointwise products (madness::mul_sparse, vmra.h); mul_sparse leaves its inputs
    // redundant, restore=True returns them to reconstructed form
    m.def("mul_sparse", [](const Function<double, 3>& f, const Function<double, 3>& g, double tol,
                           bool restore) {
        py::gil_scoped_release release;
        Function<double, 3> r = madness::mul_sparse(f, g, tol, true);
        if (restore) {
            f.change_tree_state(TreeState::reconstructed, false);
            g.change_tree_state(TreeState::reconstructed, false);
            f.world().gop.fence();
        }
        return r;
    }, py::arg("f"), py::arg("g"), py::arg("tol"), py::arg("restore") = true,
       "f * g, skipping box pairs whose estimated contribution is below tol (absolute accuracy of\n"
       "the product; tol = 0 multiplies exactly, as operator*)");
    m.def("mul_sparse", [](const Function<double, 3>& f, const std::vector<Function<double, 3>>& v,
                           double tol, bool restore) {
        py::gil_scoped_release release;
        std::vector<Function<double, 3>> r = madness::mul_sparse(f.world(), f, v, tol, true);
        if (restore) {
            f.change_tree_state(TreeState::reconstructed, false);
            for (const auto& g : v) g.change_tree_state(TreeState::reconstructed, false);
            f.world().gop.fence();
        }
        return r;
    }, py::arg("f"), py::arg("v"), py::arg("tol"), py::arg("restore") = true,
       "[f * g for g in v] in one traversal of f, screened as mul_sparse(f, g, tol)");

    // ---- complex functions (ComplexFunction3D) for the magnetic chain ----
    m.def("apply_gauss_1d", [](const Function<double_complex, 3>& f, int axis, double variance, int nbroaden) {
        py::gil_scoped_release release;
        return apply_gauss_1d(f, axis, variance, nbroaden);
    }, py::arg("f"), py::arg("axis"), py::arg("variance"), py::arg("nbroaden") = 6,
       "Complex 3D function: convolve along axis 0 or 2 with the normalized Gaussian (real kernel)");
    m.def("apply_gauss_2d", [](const Function<double_complex, 3>& f, int first, double v0, double v1,
                               int nbroaden) {
        py::gil_scoped_release release;
        return apply_gauss_2d(f, first, v0, v1, nbroaden);
    }, py::arg("f"), py::arg("first"), py::arg("v0"), py::arg("v1"), py::arg("nbroaden") = 6,
       "Complex 3D function: anisotropic 2D Gaussian convolution on axes (first, first+1)");
    m.def("mul_expi", [](const Function<double_complex, 3>& f, double kappa, double o0, double o1,
                         double tolfac, bool recursive) {
        py::gil_scoped_release release;
        return mul_expi(f, kappa, o0, o1, tolfac, recursive);
    }, py::arg("f"), py::arg("kappa"), py::arg("o0") = 0.0, py::arg("o1") = 0.0,
       py::arg("tolfac") = 0.1, py::arg("recursive") = true,
       "f exp(i theta), theta = kappa (x0-o0)(x1-o1), for a complex f; refines as mul_phase;\n"
       "truncate the result afterwards");
    m.def("mul_expi_variant", [](const Function<double_complex, 3>& f, double kappa, double o0, double o1,
                                 double tolfac, int mode) {
        py::gil_scoped_release release;
        return mul_expi_variant(f, kappa, o0, o1, tolfac, mode);
    }, py::arg("f"), py::arg("kappa"), py::arg("o0"), py::arg("o1"), py::arg("tolfac"), py::arg("mode"),
       "comparison variants of mul_expi: 0 copy+in-place, 1 refine f + out-of-place unary_op,\n"
       "2 refine f + in-place on f (overwrites f)");
    m.def("mul", [](const Function<double, 3>& f, const Function<double_complex, 3>& g) {
        py::gil_scoped_release release;
        return f * g;
    }, py::arg("f"), py::arg("g"), "real f times complex g");
    m.def("mul_sparse", [](const Function<double, 3>& f, const Function<double_complex, 3>& g, double tol,
                           bool restore) {
        py::gil_scoped_release release;
        Function<double_complex, 3> r = madness::mul_sparse(f, g, tol, true);
        if (restore) {
            f.change_tree_state(TreeState::reconstructed, false);
            g.change_tree_state(TreeState::reconstructed, false);
            f.world().gop.fence();
        }
        return r;
    }, py::arg("f"), py::arg("g"), py::arg("tol"), py::arg("restore") = true,
       "real f times complex g, screened (see mul_sparse for real functions)");
    m.def("real", [](const Function<double_complex, 3>& z) {
        py::gil_scoped_release release;
        return madness::real(z);
    }, py::arg("z"), "real part of a complex function");
    m.def("imag", [](const Function<double_complex, 3>& z) {
        py::gil_scoped_release release;
        return madness::imag(z);
    }, py::arg("z"), "imaginary part of a complex function");
    m.def("to_complex", [](const Function<double, 3>& u, const Function<double, 3>& v) {
        py::gil_scoped_release release;
        Function<double_complex, 3> z = convert<double, double_complex, 3>(u);
        Function<double_complex, 3> w = convert<double, double_complex, 3>(v);
        w.scale(double_complex(0.0, 1.0));
        Function<double_complex, 3> r = z + w;      // operator+ compresses both
        r.reconstruct();
        return r;
    }, py::arg("u"), py::arg("v"), "u + i v");

    // 2D anisotropic Gaussian convolution along axes (first, first+1) of a 3D function
    m.def("apply_gauss_2d", [](const Function<double, 3>& f, int first, double v0, double v1,
                               int nbroaden) {
        py::gil_scoped_release release;
        return apply_gauss_2d(f, first, v0, v1, nbroaden);
    }, py::arg("f"), py::arg("first"), py::arg("v0"), py::arg("v1"), py::arg("nbroaden") = 6,
       "Convolve a 3D function along axes (first, first+1), first = 0 or 1, with normalized\n"
       "Gaussians of variances v0, v1; requires FunctionDefaults2D cell and k to match those axes");

    // =====================================================================
    // Derivative operators — also held via shared_ptr (WorldObject)
    // =====================================================================

    // Free function: create derivative and apply in one shot
    m.def("diff", [](World& world, const Function<double, 1>& f, std::size_t axis) {
        py::gil_scoped_release release;
        Derivative<double, 1> D(world, axis);
        return D(f);
    }, py::arg("world"), py::arg("f"), py::arg("axis") = 0,
       "Differentiate a 1D function along the given axis");

    m.def("diff", [](World& world, const Function<double, 2>& f, std::size_t axis) {
        py::gil_scoped_release release;
        Derivative<double, 2> D(world, axis);
        return D(f);
    }, py::arg("world"), py::arg("f"), py::arg("axis"),
       "Differentiate a 2D function along the given axis");

    m.def("diff", [](World& world, const Function<double, 3>& f, std::size_t axis) {
        py::gil_scoped_release release;
        Derivative<double, 3> D(world, axis);
        return D(f);
    }, py::arg("world"), py::arg("f"), py::arg("axis"),
       "Differentiate a 3D function along the given axis");

    // Derivative as a reusable class via shared_ptr holder
    py::class_<Derivative<double, 1>, std::shared_ptr<Derivative<double, 1>>>(
        m, "Derivative1D", "First derivative operator in 1D")
        .def(py::init([](World& w, std::size_t axis) {
            return std::make_shared<Derivative<double, 1>>(w, axis);
        }), py::arg("world"), py::arg("axis") = 0)
        .def("__call__", [](Derivative<double, 1>& D,
                            const Function<double, 1>& f) {
            py::gil_scoped_release release;
            return D(f);
        }, py::arg("f"), "Apply derivative to function")
        .def("set_bspline1", &Derivative<double, 1>::set_bspline1,
             "Use the B-spline first-derivative stencil (k <= 18)")
        .def("set_ble1", &Derivative<double, 1>::set_ble1,
             "Use the BLE first-derivative stencil (k <= 15)");

    py::class_<Derivative<double, 2>, std::shared_ptr<Derivative<double, 2>>>(
        m, "Derivative2D", "First derivative operator in 2D")
        .def(py::init([](World& w, std::size_t axis) {
            return std::make_shared<Derivative<double, 2>>(w, axis);
        }), py::arg("world"), py::arg("axis"))
        .def("__call__", [](Derivative<double, 2>& D,
                            const Function<double, 2>& f) {
            py::gil_scoped_release release;
            return D(f);
        }, py::arg("f"), "Apply derivative to function")
        .def("set_bspline1", &Derivative<double, 2>::set_bspline1,
             "Use the B-spline first-derivative stencil (k <= 18)")
        .def("set_ble1", &Derivative<double, 2>::set_ble1,
             "Use the BLE first-derivative stencil (k <= 15)");

    py::class_<Derivative<double, 3>, std::shared_ptr<Derivative<double, 3>>>(
        m, "Derivative3D", "First derivative operator in 3D")
        .def(py::init([](World& w, std::size_t axis) {
            return std::make_shared<Derivative<double, 3>>(w, axis);
        }), py::arg("world"), py::arg("axis"))
        .def("__call__", [](Derivative<double, 3>& D,
                            const Function<double, 3>& f) {
            py::gil_scoped_release release;
            return D(f);
        }, py::arg("f"), "Apply derivative to function")
        .def("set_bspline1", &Derivative<double, 3>::set_bspline1,
             "Use the B-spline first-derivative stencil (k <= 18)")
        .def("set_ble1", &Derivative<double, 3>::set_ble1,
             "Use the BLE first-derivative stencil (k <= 15)");

    // =====================================================================
    // Convenience: gradient (returns list of 3 functions)
    // =====================================================================
    m.def("gradient", [](World& world, const Function<double, 3>& f) {
        py::gil_scoped_release release;
        std::vector<Function<double, 3>> grad(3);
        for (int axis = 0; axis < 3; ++axis) {
            Derivative<double, 3> D(world, axis);
            grad[axis] = D(f);
        }
        return grad;
    }, py::arg("world"), py::arg("f"),
       "Compute gradient of a 3D function, returns [df/dx, df/dy, df/dz]");
}
