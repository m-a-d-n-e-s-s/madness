/*
  pymadness - Python bindings for MADNESS
  Function<T,NDIM> bindings for Phase 1: real functions in 1D, 2D, 3D.
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/vmra.h>
#include "py_functor.h"

namespace py = pybind11;
using namespace madness;

template<typename T, std::size_t NDIM>
static void bind_function_type(py::module_& m, const char* name) {
    using FuncT = Function<T, NDIM>;
    using coordT = Vector<double, NDIM>;

    py::class_<FuncT>(m, name)
        // --- Construction ---
        .def(py::init<>(), "Create an uninitialized function")

        .def(py::init([](World& world, py::object callable, int k, double thresh) {
            std::shared_ptr<FunctionFunctorInterface<T, NDIM>> functor =
                std::make_shared<PyFunctor<T, NDIM>>(callable);
            FunctionFactory<T, NDIM> factory(world);
            factory.functor(functor);
            if (k > 0) factory.k(k);
            if (thresh > 0) factory.thresh(thresh);
            // Release the GIL so MADNESS worker threads can call back
            // into Python via PyFunctor::operator() during projection.
            py::gil_scoped_release release;
            return FuncT(factory);
        }),
            py::arg("world"), py::arg("f"),
            py::arg("k") = -1, py::arg("thresh") = -1.0,
            "Create function by projecting a Python callable.\n"
            "The callable must accept a numpy array of shape (NDIM,) and return a scalar.")

        // --- Point evaluation ---
        .def("__call__", [](const FuncT& f, py::array_t<double> r) -> T {
            auto buf = r.unchecked<1>();
            if (buf.shape(0) != static_cast<py::ssize_t>(NDIM)) {
                throw std::invalid_argument(
                    "Coordinate array must have " + std::to_string(NDIM) + " elements");
            }
            coordT coord;
            for (std::size_t i = 0; i < NDIM; ++i) coord[i] = buf(i);
            // Release GIL so MADNESS worker threads can process the
            // eval task and resolve the Future.
            T result;
            {
                py::gil_scoped_release release;
                if (!f.is_reconstructed()) f.reconstruct();
                result = f.eval(coord).get();
            }
            return result;
        }, py::arg("r"), "Evaluate function at point r (numpy array of length NDIM)")

        // --- Arithmetic operators ---
        .def("__add__", [](const FuncT& a, const FuncT& b) {
            py::gil_scoped_release release;
            return a + b;
        })
        .def("__sub__", [](const FuncT& a, const FuncT& b) {
            py::gil_scoped_release release;
            return a - b;
        })
        .def("__mul__", [](const FuncT& a, const FuncT& b) {
            py::gil_scoped_release release;
            return a * b;
        }, "Pointwise multiplication of two functions")
        .def("__mul__", [](const FuncT& a, T s) {
            py::gil_scoped_release release;
            return a * s;
        })
        .def("__rmul__", [](const FuncT& a, T s) {
            py::gil_scoped_release release;
            return s * a;
        })
        .def("__neg__", [](const FuncT& a) {
            py::gil_scoped_release release;
            return a * T(-1.0);
        })
        .def("__iadd__", [](FuncT& a, const FuncT& b) -> FuncT& {
            py::gil_scoped_release release;
            a += b;
            return a;
        }, py::return_value_policy::reference_internal)
        .def("__isub__", [](FuncT& a, const FuncT& b) -> FuncT& {
            py::gil_scoped_release release;
            a -= b;
            return a;
        }, py::return_value_policy::reference_internal)
        .def("__imul__", [](FuncT& a, T s) -> FuncT& {
            py::gil_scoped_release release;
            a.scale(s);
            return a;
        }, py::return_value_policy::reference_internal)

        // --- Norms and integrals ---
        .def("norm2", [](const FuncT& f) {
            py::gil_scoped_release release;
            return f.norm2();
        }, "L2 norm: sqrt(integral |f|^2)")
        .def("trace", [](const FuncT& f) {
            py::gil_scoped_release release;
            return f.trace();
        }, "Integral of f over the simulation cell")
        .def("inner", [](const FuncT& f, const FuncT& g) {
            py::gil_scoped_release release;
            return f.inner(g);
        }, py::arg("g"), "Inner product: integral f(x)*g(x) dx")

        // --- Tree state management ---
        .def("compress", [](const FuncT& f) {
            py::gil_scoped_release release;
            f.compress();
        }, "Transform to compressed (wavelet) representation")
        .def("reconstruct", [](const FuncT& f) {
            py::gil_scoped_release release;
            f.reconstruct();
        }, "Transform to reconstructed (scaling function) representation")
        .def("truncate", [](FuncT& f, double tol) {
            py::gil_scoped_release release;
            f.truncate(tol);
        }, py::arg("tol") = 0.0,
           "Truncate small wavelet coefficients (0 = use default threshold)")

        // --- In-place operations ---
        .def("scale", [](FuncT& f, T s) {
            py::gil_scoped_release release;
            f.scale(s);
        }, py::arg("s"), "Scale function in-place by constant s")
        .def("gaxpy", [](FuncT& f, T alpha, const FuncT& g, T beta) {
            py::gil_scoped_release release;
            f.gaxpy(alpha, g, beta);
        }, py::arg("alpha"), py::arg("g"), py::arg("beta"),
           "In-place f = alpha*f + beta*g")
        .def("square", [](FuncT& f) {
            py::gil_scoped_release release;
            f.square();
        }, "Square the function in-place: f(x) = f(x)^2")
        .def("abs", [](FuncT& f) {
            py::gil_scoped_release release;
            f.abs();
        }, "Absolute value in-place: f(x) = |f(x)|")
        .def("abs_square", [](FuncT& f) {
            py::gil_scoped_release release;
            f.abs_square();
        }, "Absolute square in-place: f(x) = |f(x)|^2")

        // --- Properties ---
        .def("k", &FuncT::k, "Wavelet order")
        .def("thresh", &FuncT::thresh, "Truncation threshold")
        .def("max_depth", &FuncT::max_depth, "Maximum depth of the function tree")
        .def("tree_size", &FuncT::tree_size, "Total number of nodes in the tree")
        .def("is_initialized", &FuncT::is_initialized,
            "True if the function has been initialized")
        .def("is_compressed", &FuncT::is_compressed,
            "True if in compressed (wavelet) representation")
        .def("is_reconstructed", &FuncT::is_reconstructed,
            "True if in reconstructed (scaling function) representation")
        .def("size", &FuncT::size,
            "Number of coefficients in the function tree")

        // --- Copy ---
        .def("copy", [](const FuncT& f) {
            py::gil_scoped_release release;
            return madness::copy(f);
        }, "Create a deep copy of this function")

        // --- Save / Load ---
        .def("save", [](const FuncT& f, const std::string& filename) {
            archive::BinaryFstreamOutputArchive ar(filename.c_str());
            f.store(ar);
        }, py::arg("filename"), "Save function to binary file")
        .def_static("load", [](World& w, const std::string& filename) {
            FuncT f;
            archive::BinaryFstreamInputArchive ar(filename.c_str());
            f.load(w, ar);
            return f;
        }, py::arg("world"), py::arg("filename"),
           "Load function from binary file")

        // --- Print info ---
        .def("print_size", [](const FuncT& f, const std::string& msg) {
            f.print_size(msg);
        }, py::arg("msg") = "", "Print tree size information")

        .def("__repr__", [](const FuncT& f) {
            std::ostringstream os;
            os << "Function<" << typeid(T).name() << "," << NDIM << ">(";
            if (f.is_initialized()) {
                os << "k=" << f.k()
                   << ", thresh=" << f.thresh()
                   << ", nodes=" << f.tree_size()
                   << ", depth=" << f.max_depth();
            } else {
                os << "uninitialized";
            }
            os << ")";
            return os.str();
        });
}

// --- Free functions operating on Function objects ---
template<typename T, std::size_t NDIM>
static void bind_function_ops(py::module_& m) {
    using FuncT = Function<T, NDIM>;

    // copy() free function
    m.def("copy", [](const FuncT& f) {
        py::gil_scoped_release release;
        return madness::copy(f);
    }, py::arg("f"), "Deep copy of a function");

    // inner() free function
    m.def("inner", [](const FuncT& f, const FuncT& g) {
        py::gil_scoped_release release;
        return madness::inner(f, g);
    }, py::arg("f"), py::arg("g"), "Inner product of two functions");
}

void bind_function(py::module_& m) {
    // Phase 1: real functions in 1D, 2D, 3D
    bind_function_type<double, 1>(m, "Function1D");
    bind_function_type<double, 2>(m, "Function2D");
    bind_function_type<double, 3>(m, "Function3D");

    // Complex 3D (needed for periodic systems)
    bind_function_type<double_complex, 3>(m, "ComplexFunction3D");

    // Free functions - overloaded for each type via pybind11 dispatch
    bind_function_ops<double, 1>(m);
    bind_function_ops<double, 2>(m);
    bind_function_ops<double, 3>(m);
    bind_function_ops<double_complex, 3>(m);
}
