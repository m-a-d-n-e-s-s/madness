/*
  pymadness - Python bindings for MADNESS
  Adapter to use Python callables as MADNESS FunctionFunctorInterface.
*/

#ifndef PYMADNESS_PY_FUNCTOR_H
#define PYMADNESS_PY_FUNCTOR_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <madness/mra/function_interface.h>

namespace py = pybind11;

/// Wraps a Python callable as a MADNESS FunctionFunctorInterface<T,NDIM>.
///
/// Supports two calling conventions:
///   - Scalar: f(r) where r has shape (NDIM,), returns scalar
///   - Vectorized: f(r) where r has shape (npts, NDIM), returns array of shape (npts,)
///
/// The vectorized path is always used (supports_vectorized() returns true).
/// On the first batch call, we probe the callable to detect which convention
/// it uses and cache the result.
template<typename T, std::size_t NDIM>
class PyFunctor : public madness::FunctionFunctorInterface<T, NDIM> {
    py::object py_callable_;
    mutable int vectorized_mode_ = 0;  // 0=unknown, 1=returns array, 2=returns scalar

public:
    explicit PyFunctor(py::object f) : py_callable_(std::move(f)) {}

    T operator()(const madness::Vector<double, NDIM>& r) const override {
        py::gil_scoped_acquire gil;

        // Build a numpy array from the coordinate vector
        py::array_t<double> arr(NDIM);
        auto buf = arr.mutable_unchecked<1>();
        for (std::size_t i = 0; i < NDIM; ++i) {
            buf(i) = r[i];
        }

        py::object result = py_callable_(arr);
        return result.cast<T>();
    }

    bool supports_vectorized() const override { return true; }

private:
    /// Core vectorized evaluation: builds numpy array from coordinate pointers,
    /// calls the Python callable once, and fills fvals.
    void eval_vectorized(const madness::Vector<double*, NDIM>& xvals, T* fvals, int npts) const {
        py::gil_scoped_acquire gil;

        // Build a numpy array of shape (npts, NDIM) from the coordinate arrays
        py::array_t<double> coords({static_cast<py::ssize_t>(npts),
                                     static_cast<py::ssize_t>(NDIM)});
        auto cbuf = coords.mutable_unchecked<2>();
        for (int i = 0; i < npts; ++i) {
            for (std::size_t d = 0; d < NDIM; ++d) {
                cbuf(i, d) = xvals[d][i];
            }
        }

        // Probe on first call to determine if callable returns array or scalar
        if (vectorized_mode_ == 0) {
            try {
                py::object result = py_callable_(coords);
                py::array_t<double> arr = result.cast<py::array_t<double>>();
                if (arr.ndim() == 1 && arr.shape(0) == npts) {
                    vectorized_mode_ = 1;  // returns array
                    auto rbuf = arr.unchecked<1>();
                    for (int i = 0; i < npts; ++i) {
                        fvals[i] = static_cast<T>(rbuf(i));
                    }
                    return;
                }
            } catch (...) {
                // Call failed or didn't return a proper array — fall through
                // Clear any pending Python exception
                PyErr_Clear();
            }
            // Fall back to per-point evaluation in Python.
            vectorized_mode_ = 2;
            eval_scalar_loop(coords, fvals, npts);
            return;
        }

        if (vectorized_mode_ == 1) {
            py::object result = py_callable_(coords);
            py::array_t<double> arr = result.cast<py::array_t<double>>();
            auto rbuf = arr.unchecked<1>();
            for (int i = 0; i < npts; ++i) {
                fvals[i] = static_cast<T>(rbuf(i));
            }
        } else {
            eval_scalar_loop(coords, fvals, npts);
        }
    }

    /// Evaluate the callable per-point using scalar convention (still 1 GIL acquisition).
    void eval_scalar_loop(const py::array_t<double>& coords, T* fvals, int npts) const {
        auto cbuf = coords.unchecked<2>();
        for (int i = 0; i < npts; ++i) {
            py::array_t<double> pt(NDIM);
            auto pbuf = pt.mutable_unchecked<1>();
            for (std::size_t d = 0; d < NDIM; ++d) {
                pbuf(d) = cbuf(i, d);
            }
            py::object result = py_callable_(pt);
            fvals[i] = result.cast<T>();
        }
    }

public:
    // Vectorized operator() overrides for NDIM = 1, 2, 3
    // These match the signatures in FunctionFunctorInterface (function_interface.h:99-121).
    // We provide all three and only the one matching our NDIM will be called.

    void operator()(const madness::Vector<double*, 1>& xvals, T* fvals, int npts) const override {
        if constexpr (NDIM == 1) {
            eval_vectorized(xvals, fvals, npts);
        }
    }

    void operator()(const madness::Vector<double*, 2>& xvals, T* fvals, int npts) const override {
        if constexpr (NDIM == 2) {
            eval_vectorized(xvals, fvals, npts);
        }
    }

    void operator()(const madness::Vector<double*, 3>& xvals, T* fvals, int npts) const override {
        if constexpr (NDIM == 3) {
            eval_vectorized(xvals, fvals, npts);
        }
    }
};

#endif  // PYMADNESS_PY_FUNCTOR_H
