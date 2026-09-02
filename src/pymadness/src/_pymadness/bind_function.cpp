/*
  pymadness - Python bindings for MADNESS
  Function<T,NDIM> bindings for dimensions 1–6.
*/

#include <cstdint>
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

/// A pointwise operator over the scaling-function VALUES of several functions.
///
/// The Python analogue of chem's `xc_potential` functor: MADNESS hands it the
/// value tensors of every input function inside one box, and it returns the
/// value tensors of the outputs on that same box.  This is the primitive that
/// makes a pointwise composition -- evaluating libxc, say -- expressible
/// without a per-point MRA eval.
///
/// Two threading subtleties, both load-bearing:
///
///  * `multi_to_multi_op_values` COPIES the op into every task, on worker
///    threads that do not hold the GIL.  Copying a `py::object` there would
///    race on the refcount, so the callable is held behind a `shared_ptr`
///    (atomic to copy) whose deleter reacquires the GIL to destroy it.
///  * The caller must NOT hold the GIL across the fence, or worker threads
///    calling back into Python deadlock.  Hence the release in the binding
///    below and the acquire here.
class PyMultiOp {
public:
    PyMultiOp(py::object fn, std::size_t n_out, std::size_t k)
        : n_out_(n_out), k_(k) {
        // custom deleter: the last owner may be a worker thread without the GIL
        fn_ = std::shared_ptr<py::object>(
            new py::object(std::move(fn)),
            [](py::object* p) {
                py::gil_scoped_acquire acquire;
                delete p;
            });
    }

    std::size_t get_result_size() const { return n_out_; }

    std::vector<Tensor<double>> operator()(
            const Key<3>& key,
            const std::vector<Tensor<double>>& t) const {
        py::gil_scoped_acquire acquire;

        // This box's quadrature points, in user coordinates.  Same construction
        // fcube() uses, and the same one chem's nemo_u1_functors::append writes
        // out (xcfunctional.h:762-773) -- it is what lets an analytic functor
        // such as U1 be evaluated pointwise instead of being projected first.
        // Ordering matches the C-order flattening of the value tensors: z runs
        // fastest.
        const Tensor<double>& qx =
            FunctionCommonData<double, 3>::get(k_).quad_x;
        const long npt = qx.dim(0);
        const double h = std::pow(0.5, double(key.level()));
        const Tensor<double>& cell = FunctionDefaults<3>::get_cell();
        const Tensor<double>& cw = FunctionDefaults<3>::get_cell_width();

        py::array_t<double> coords({static_cast<py::ssize_t>(npt * npt * npt),
                                    static_cast<py::ssize_t>(3)});
        double* cp = coords.mutable_data();
        long idx = 0;
        for (long i = 0; i < npt; ++i) {
            const double x =
                cell(0, 0) + h * cw[0] * (key.translation()[0] + qx(i));
            for (long j = 0; j < npt; ++j) {
                const double y =
                    cell(1, 0) + h * cw[1] * (key.translation()[1] + qx(j));
                for (long l = 0; l < npt; ++l, ++idx) {
                    const double z =
                        cell(2, 0) + h * cw[2] * (key.translation()[2] + qx(l));
                    cp[3 * idx + 0] = x;
                    cp[3 * idx + 1] = y;
                    cp[3 * idx + 2] = z;
                }
            }
        }

        py::list inputs;
        for (const auto& x : t) {
            if (x.size() == 0) {
                inputs.append(py::none());
                continue;
            }
            // hand out a flat float64 copy; the callable sees k^3 values
            py::array_t<double> a(static_cast<py::ssize_t>(x.size()));
            std::copy(x.ptr(), x.ptr() + x.size(), a.mutable_data());
            inputs.append(std::move(a));
        }

        py::object res = (*fn_)(coords, inputs);
        py::sequence seq = py::cast<py::sequence>(res);
        if (static_cast<std::size_t>(py::len(seq)) != n_out_) {
            throw std::runtime_error(
                "pointwise op returned " + std::to_string(py::len(seq))
                + " arrays, expected " + std::to_string(n_out_));
        }

        const long n_values = static_cast<long>(k_ * k_ * k_);
        std::vector<Tensor<double>> out;
        out.reserve(n_out_);
        for (std::size_t i = 0; i < n_out_; ++i) {
            auto a = py::cast<py::array_t<double, py::array::c_style
                                          | py::array::forcecast>>(seq[i]);
            if (a.size() != n_values) {
                throw std::runtime_error(
                    "pointwise op output " + std::to_string(i) + " has "
                    + std::to_string(a.size()) + " values, expected "
                    + std::to_string(n_values));
            }
            Tensor<double> r(k_, k_, k_);
            std::copy(a.data(), a.data() + n_values, r.ptr());
            out.push_back(r);
        }
        return out;
    }

private:
    std::shared_ptr<py::object> fn_;
    std::size_t n_out_;
    std::size_t k_;
};

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
            "Two calling conventions are supported (probed at construction time):\n"
            "  - Vectorized: f(r) where r has shape (npts, NDIM) and returns\n"
            "    an array of shape (npts,).  Preferred: the GIL is acquired\n"
            "    only once per box, making projection much faster.\n"
            "  - Scalar: f(r) where r has shape (NDIM,) and returns a float.\n"
            "    Used as a fallback if the callable does not accept the\n"
            "    vectorized form.")

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
        .def("unaryop", [](FuncT& f, const std::string& op) {
            // Use SimpleUnaryOpWrapper via the function-pointer overload
            // of unaryop(T (*f)(T)). We need a plain function pointer,
            // so define static helpers.
            using fptr_t = T(*)(T);
            fptr_t fn = nullptr;
            if (op == "exp") {
                fn = +[](T x) -> T { return std::exp(x); };
            } else if (op == "log") {
                fn = +[](T x) -> T { return std::log(x); };
            } else if (op == "abs") {
                fn = +[](T x) -> T { return std::abs(x); };
            } else if (op == "sqrt") {
                fn = +[](T x) -> T { return std::sqrt(x); };
            } else {
                throw std::invalid_argument(
                    "Unknown unary op '" + op + "'. "
                    "Supported: exp, log, abs, sqrt");
            }
            py::gil_scoped_release release;
            f.unaryop(fn);
        }, py::arg("op"),
           "Apply a pointwise unary operation in-place.\n"
           "Supported: 'exp', 'log', 'abs', 'sqrt'")

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
        // Use ParallelArchive to match the C++ save()/load() free functions
        .def("save", [](const FuncT& f, const std::string& filename) {
            py::gil_scoped_release release;
            archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(f.world(), filename.c_str(), 1);
            ar & f;
        }, py::arg("filename"), "Save function to binary file")
        .def_static("load", [](World& w, const std::string& filename) {
            FuncT f;
            {
                py::gil_scoped_release release;
                archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(w, filename.c_str());
                ar & f;
            }
            return f;
        }, py::arg("world"), py::arg("filename"),
           "Load function from binary file")

        // --- Construct from C function pointer (numba @cfunc support) ---
        .def_static("from_cfunc", [](World& world, uintptr_t addr, int k, double thresh) {
            // Cast the integer address to a C function pointer
            auto fptr = reinterpret_cast<T (*)(const Vector<double, NDIM>&)>(addr);
            FunctionFactory<T, NDIM> factory(world);
            factory.f(fptr);
            if (k > 0) factory.k(k);
            if (thresh > 0) factory.thresh(thresh);
            py::gil_scoped_release release;
            return FuncT(factory);
        },
            py::arg("world"), py::arg("cfunc_address"),
            py::arg("k") = -1, py::arg("thresh") = -1.0,
            "Create function from a C function pointer (e.g. from numba @cfunc).\n"
            "The function must have signature T(const Vector<double,NDIM>&).")

        // --- Grid evaluation ---
        .def("eval_cube", [](const FuncT& f,
                             py::array_t<double, py::array::c_style> cell_arr,
                             std::vector<long> npt) {
            // cell_arr is shape (NDIM, 2): [[lo0, hi0], [lo1, hi1], ...]
            auto buf = cell_arr.unchecked<2>();
            if (buf.shape(0) != static_cast<py::ssize_t>(NDIM) || buf.shape(1) != 2) {
                throw std::invalid_argument(
                    "cell must have shape (" + std::to_string(NDIM) + ", 2)");
            }
            if (npt.size() != NDIM) {
                throw std::invalid_argument(
                    "npt must have " + std::to_string(NDIM) + " elements");
            }
            Tensor<double> cell(NDIM, 2L);
            for (std::size_t d = 0; d < NDIM; ++d) {
                cell(d, 0L) = buf(d, 0);
                cell(d, 1L) = buf(d, 1);
            }
            Tensor<T> result;
            {
                py::gil_scoped_release release;
                result = f.eval_cube(cell, npt);
            }
            // Convert Tensor<T> to numpy array
            std::vector<py::ssize_t> shape(result.ndim());
            for (long i = 0; i < result.ndim(); ++i) {
                shape[i] = result.dim(i);
            }
            py::array_t<T> arr(shape);
            T* dst = arr.mutable_data();
            const T* src = result.ptr();
            std::copy(src, src + result.size(), dst);
            return arr;
        },
            py::arg("cell"), py::arg("npt"),
            "Evaluate function on a regular grid.\n\n"
            "Args:\n"
            "    cell: numpy array of shape (NDIM, 2) with [lo, hi] per dimension\n"
            "    npt: list of ints giving number of points per dimension\n\n"
            "Returns:\n"
            "    numpy array of shape (npt[0], npt[1], ..., npt[NDIM-1])")

        // --- Print info ---
        .def("print_size", [](const FuncT& f, const std::string& msg) {
            py::gil_scoped_release release;
            f.print_size(msg);
        }, py::arg("msg") = "", "Print tree size information")

        .def("__repr__", [](const FuncT& f) {
            std::ostringstream os;
            os << "Function<" << typeid(T).name() << "," << NDIM << ">(";
            if (f.is_initialized()) {
                int k;
                double thresh;
                long nodes, depth;
                {
                    py::gil_scoped_release release;
                    k = f.k();
                    thresh = f.thresh();
                    nodes = f.tree_size();
                    depth = f.max_depth();
                }
                os << "k=" << k
                   << ", thresh=" << thresh
                   << ", nodes=" << nodes
                   << ", depth=" << depth;
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

// --- FunctionFactory<T,NDIM> Python binding ---
template<typename T, std::size_t NDIM>
static void bind_factory_type(py::module_& m, const char* name) {
    using FactoryT = FunctionFactory<T, NDIM>;

    py::class_<FactoryT>(m, name,
        "Chainable factory for configuring Function construction")
        .def(py::init<World&>(), py::arg("world"))

        // Functor / callable
        .def("functor", [](FactoryT& self, py::object callable) -> FactoryT& {
            std::shared_ptr<FunctionFunctorInterface<T, NDIM>> f =
                std::make_shared<PyFunctor<T, NDIM>>(callable);
            self.functor(f);
            return self;
        }, py::arg("f"), py::return_value_policy::reference_internal,
           "Set Python callable as functor")

        // Wavelet settings
        .def("k", [](FactoryT& self, int k) -> FactoryT& {
            return self.k(k);
        }, py::arg("k"), py::return_value_policy::reference_internal,
           "Set wavelet order")
        .def("thresh", [](FactoryT& self, double t) -> FactoryT& {
            return self.thresh(t);
        }, py::arg("thresh"), py::return_value_policy::reference_internal,
           "Set truncation threshold")

        // Refinement
        // Force deep refinement in the boxes containing given points -- the
        // only way to make a nucleus resolved beyond what the adaptive error
        // test asks for.  Without it a cusp sits inside a box the projector
        // considers converged, and every downstream derivative inherits that.
        .def("special_points", [](FactoryT& self,
                                  const std::vector<std::vector<double>>& pts)
                                  -> FactoryT& {
            std::vector<Vector<double, NDIM>> sp;
            sp.reserve(pts.size());
            for (const auto& p : pts) {
                if (p.size() != NDIM)
                    throw std::runtime_error(
                        "special_points: each point needs " +
                        std::to_string(NDIM) + " coordinates");
                Vector<double, NDIM> v;
                for (std::size_t i = 0; i < NDIM; ++i) v[i] = p[i];
                sp.push_back(v);
            }
            return self.special_points(sp);
        }, py::arg("points"),
           "Refine boxes containing these points to at least special_level")
        .def("special_level", [](FactoryT& self, int level) -> FactoryT& {
            return self.special_level(level);
        }, py::arg("level"),
           "Minimum refinement level for boxes holding a special point")
        .def("initial_level", [](FactoryT& self, int level) -> FactoryT& {
            return self.initial_level(level);
        }, py::arg("level"), py::return_value_policy::reference_internal,
           "Set initial projection level")
        .def("max_refine_level", [](FactoryT& self, int level) -> FactoryT& {
            return self.max_refine_level(level);
        }, py::arg("level"), py::return_value_policy::reference_internal,
           "Set maximum adaptive refinement level")
        .def("refine", [](FactoryT& self, bool flag) -> FactoryT& {
            return self.refine(flag);
        }, py::arg("flag") = true, py::return_value_policy::reference_internal,
           "Enable/disable adaptive refinement")
        .def("norefine", [](FactoryT& self) -> FactoryT& {
            return self.norefine(true);
        }, py::return_value_policy::reference_internal,
           "Disable adaptive refinement")
        .def("autorefine", [](FactoryT& self) -> FactoryT& {
            return self.autorefine();
        }, py::return_value_policy::reference_internal,
           "Enable auto-refinement in operations")
        .def("noautorefine", [](FactoryT& self) -> FactoryT& {
            return self.noautorefine();
        }, py::return_value_policy::reference_internal,
           "Disable auto-refinement")

        // Truncation
        .def("truncate_mode", [](FactoryT& self, int mode) -> FactoryT& {
            return self.truncate_mode(mode);
        }, py::arg("mode"), py::return_value_policy::reference_internal,
           "Set truncation mode (0, 1, or 2)")
        .def("truncate_on_project", [](FactoryT& self) -> FactoryT& {
            return self.truncate_on_project();
        }, py::return_value_policy::reference_internal,
           "Enable truncation during projection")
        .def("notruncate_on_project", [](FactoryT& self) -> FactoryT& {
            return self.notruncate_on_project();
        }, py::return_value_policy::reference_internal,
           "Disable truncation during projection")

        // Fencing
        .def("fence", [](FactoryT& self, bool flag) -> FactoryT& {
            return self.fence(flag);
        }, py::arg("flag") = true, py::return_value_policy::reference_internal,
           "Enable/disable fence")
        .def("nofence", [](FactoryT& self) -> FactoryT& {
            return self.nofence();
        }, py::return_value_policy::reference_internal,
           "Disable fence")

        // Empty function
        .def("empty", [](FactoryT& self) -> FactoryT& {
            return self.empty();
        }, py::return_value_policy::reference_internal,
           "Create an empty function (no projection)")

        // Build the function
        .def("create", [](FactoryT& self) {
            using FuncT = Function<T, NDIM>;
            py::gil_scoped_release release;
            return FuncT(self);
        }, "Build and return the Function from this factory's settings");
}

void bind_function(py::module_& m) {
    // Real functions in 1D–6D
    bind_function_type<double, 1>(m, "Function1D");
    bind_function_type<double, 2>(m, "Function2D");
    bind_function_type<double, 3>(m, "Function3D");
    bind_function_type<double, 4>(m, "Function4D");
    bind_function_type<double, 5>(m, "Function5D");
    bind_function_type<double, 6>(m, "Function6D");

    // Complex 3D (needed for periodic systems)
    bind_function_type<double_complex, 3>(m, "ComplexFunction3D");

    // Free functions - overloaded for each type via pybind11 dispatch
    bind_function_ops<double, 1>(m);
    bind_function_ops<double, 2>(m);
    bind_function_ops<double, 3>(m);
    bind_function_ops<double, 4>(m);
    bind_function_ops<double, 5>(m);
    bind_function_ops<double, 6>(m);
    bind_function_ops<double_complex, 3>(m);

    // FunctionFactory bindings
    bind_factory_type<double, 1>(m, "FunctionFactory1D");
    bind_factory_type<double, 2>(m, "FunctionFactory2D");
    bind_factory_type<double, 3>(m, "FunctionFactory3D");
    bind_factory_type<double, 4>(m, "FunctionFactory4D");
    bind_factory_type<double, 5>(m, "FunctionFactory5D");
    bind_factory_type<double, 6>(m, "FunctionFactory6D");
    bind_factory_type<double_complex, 3>(m, "ComplexFunctionFactory3D");

    // -----------------------------------------------------------------
    // Pointwise composition over a common tree
    // -----------------------------------------------------------------
    m.def("refine_to_common_level",
          [](std::vector<Function<double, 3>> vf) {
              py::gil_scoped_release release;
              if (vf.empty()) return vf;
              refine_to_common_level(vf[0].world(), vf, true);
              return vf;
          },
          py::arg("functions"),
          "Refine all functions to a common (finest) level, in place.\n"
          "Required before multiop_values.  Note this lifts EVERY function to\n"
          "the deepest one's level, which is how a single deep intermediate\n"
          "taxes a whole pipeline.");

    m.def("multiop_values",
          [](std::vector<Function<double, 3>> vin, py::object fn,
             std::size_t n_out) {
              if (vin.empty())
                  throw std::runtime_error("multiop_values needs at least one input");
              const std::size_t k = static_cast<std::size_t>(vin[0].k());
              for (std::size_t i = 1; i < vin.size(); ++i) {
                  if (static_cast<std::size_t>(vin[i].k()) != k)
                      throw std::runtime_error(
                          "multiop_values: inputs have different k (" +
                          std::to_string(k) + " vs " +
                          std::to_string(vin[i].k()) + "); the quadrature "
                          "points would not line up with the values");
              }
              PyMultiOp op(std::move(fn), n_out, k);
              py::gil_scoped_release release;
              // The op reads SCALING-FUNCTION coefficients, so every input must
              // be in reconstructed form.  Any intervening arithmetic (a
              // subtraction, an inner product) leaves a function COMPRESSED,
              // and reading wavelet coefficients as values fails silently --
              // plausible magnitudes, malformed output tree, no exception.  So
              // restore the precondition rather than trust the caller; this is
              // a no-op when the inputs are already reconstructed.
              reconstruct(vin[0].world(), vin, true);
              return multi_to_multi_op_values(op, vin, true);
          },
          py::arg("functions"), py::arg("op"), py::arg("n_out"),
          "Apply a pointwise operator to the VALUES of several functions.\n\n"
          "`op` is called once per box as op(coords, [a0, a1, ...]) ->\n"
          "[b0, ..., bN].  `coords` is an (k^3, 3) array of that box's\n"
          "quadrature points in user coordinates; each ai is a flat float64\n"
          "array of k^3 scaling-function values, in the same order, and each\n"
          "bj must be the same length.  Inputs must already be at a common\n"
          "level (see refine_to_common_level).\n\n"
          "`coords` is what lets an analytic functor be evaluated POINTWISE\n"
          "rather than projected first -- the difference between the two is\n"
          "measurable and, for a cuspy functor, large.\n\n"
          "The result is projected onto the EXISTING tree with no refinement\n"
          "and no error test -- the same non-adaptive projection MADNESS uses\n"
          "for the XC potential.  It is faithful, not automatically accurate.");
}
