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
using SepConv1D = SeparatedConvolution<double, 1>;

void bind_operators(py::module_& m) {

    // =====================================================================
    // SeparatedConvolution — held via shared_ptr (WorldObject not copyable)
    // =====================================================================
    py::class_<SepConv3D, std::shared_ptr<SepConv3D>>(m, "SeparatedConvolution3D",
        "Integral operator represented as a separated convolution in 3D");

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
        }, py::arg("f"), "Apply derivative to function");

    py::class_<Derivative<double, 2>, std::shared_ptr<Derivative<double, 2>>>(
        m, "Derivative2D", "First derivative operator in 2D")
        .def(py::init([](World& w, std::size_t axis) {
            return std::make_shared<Derivative<double, 2>>(w, axis);
        }), py::arg("world"), py::arg("axis"))
        .def("__call__", [](Derivative<double, 2>& D,
                            const Function<double, 2>& f) {
            py::gil_scoped_release release;
            return D(f);
        }, py::arg("f"), "Apply derivative to function");

    py::class_<Derivative<double, 3>, std::shared_ptr<Derivative<double, 3>>>(
        m, "Derivative3D", "First derivative operator in 3D")
        .def(py::init([](World& w, std::size_t axis) {
            return std::make_shared<Derivative<double, 3>>(w, axis);
        }), py::arg("world"), py::arg("axis"))
        .def("__call__", [](Derivative<double, 3>& D,
                            const Function<double, 3>& f) {
            py::gil_scoped_release release;
            return D(f);
        }, py::arg("f"), "Apply derivative to function");

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
