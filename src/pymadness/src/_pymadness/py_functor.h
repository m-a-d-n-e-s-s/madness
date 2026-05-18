/*
  pymadness - Python bindings for MADNESS
  Adapter to use Python callables as MADNESS FunctionFunctorInterface.
*/

#ifndef PYMADNESS_PY_FUNCTOR_H
#define PYMADNESS_PY_FUNCTOR_H

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/complex.h>
#include <madness/mra/function_interface.h>
#include <atomic>
#include <stdexcept>
#include <string>
#include <thread>

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
    mutable std::atomic<int> vectorized_mode_{0};  // -1=probing, 0=unknown, 1=returns array, 2=returns scalar

public:
    explicit PyFunctor(py::object f) : py_callable_(std::move(f)) {}

    ~PyFunctor() {
        // py::object's destructor decrements Python's refcount, which requires
        // the GIL.  MADNESS may destroy the shared_ptr<PyFunctor> from a worker
        // thread that does not hold the GIL, so we must acquire it here.
        py::gil_scoped_acquire gil;
        py_callable_ = py::none();
    }

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

        // Probe on first call to determine if callable returns array or scalar.
        // Note: eval_vectorized always holds the GIL (acquired above), so only one
        // thread can execute this block at a time.  The atomic CAS still ensures
        // correct memory ordering across threads.
        int expected = 0;
        if (vectorized_mode_.compare_exchange_strong(expected, -1,
                std::memory_order_acq_rel, std::memory_order_acquire)) {
            // We won the probe; expected was 0, now set to -1 (in-progress sentinel)
            int result_mode = 2;  // default: fall back to per-point scalar loop
            try {
                py::object result = py_callable_(coords);
                py::array_t<T> arr = result.cast<py::array_t<T>>();
                if (arr.ndim() == 1 && arr.shape(0) == npts) {
                    result_mode = 1;
                    auto rbuf = arr.template unchecked<1>();
                    for (int i = 0; i < npts; ++i) {
                        fvals[i] = rbuf(i);
                    }
                }
            } catch (...) {
                // Call failed or didn't return a proper array — fall through
                // Clear any pending Python exception
                PyErr_Clear();
            }
            vectorized_mode_.store(result_mode, std::memory_order_release);
            if (result_mode == 2) {
                eval_scalar_loop(coords, fvals, npts);
            }
            return;
        }

        // If another thread is probing (-1), we must release the GIL while
        // waiting.  The probing thread needs the GIL to finish its Python call,
        // so spinning here with the GIL held would deadlock.
        int mode = expected;  // CAS failure: expected holds the current value
        while (mode == -1) {
            {
                py::gil_scoped_release release;
                std::this_thread::yield();
            }
            mode = vectorized_mode_.load(std::memory_order_acquire);
        }

        if (mode == 1) {
            py::object result = py_callable_(coords);
            py::array_t<T> arr = result.cast<py::array_t<T>>();
            if (arr.ndim() != 1 || arr.shape(0) != npts) {
                throw std::runtime_error(
                    "PyFunctor: vectorized callable returned array with wrong shape; "
                    "expected 1D array of length " + std::to_string(npts));
            }
            auto rbuf = arr.template unchecked<1>();
            for (int i = 0; i < npts; ++i) {
                fvals[i] = rbuf(i);
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
    // Vectorized operator() overrides for NDIM = 1..6
    // These match the signatures in FunctionFunctorInterface (function_interface.h:99-121).
    // We provide all six and only the one matching our NDIM will be called.

    void operator()(const madness::Vector<double*, 1>& xvals, T* fvals, int npts) const override {
        if constexpr (NDIM == 1) { eval_vectorized(xvals, fvals, npts); }
    }

    void operator()(const madness::Vector<double*, 2>& xvals, T* fvals, int npts) const override {
        if constexpr (NDIM == 2) { eval_vectorized(xvals, fvals, npts); }
    }

    void operator()(const madness::Vector<double*, 3>& xvals, T* fvals, int npts) const override {
        if constexpr (NDIM == 3) { eval_vectorized(xvals, fvals, npts); }
    }

    void operator()(const madness::Vector<double*, 4>& xvals, T* fvals, int npts) const override {
        if constexpr (NDIM == 4) { eval_vectorized(xvals, fvals, npts); }
    }

    void operator()(const madness::Vector<double*, 5>& xvals, T* fvals, int npts) const override {
        if constexpr (NDIM == 5) { eval_vectorized(xvals, fvals, npts); }
    }

    void operator()(const madness::Vector<double*, 6>& xvals, T* fvals, int npts) const override {
        if constexpr (NDIM == 6) { eval_vectorized(xvals, fvals, npts); }
    }
};

#endif  // PYMADNESS_PY_FUNCTOR_H
