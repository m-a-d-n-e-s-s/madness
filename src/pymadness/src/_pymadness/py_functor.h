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
/// The callable must accept a numpy array of shape (NDIM,) and return a scalar
/// convertible to T.
template<typename T, std::size_t NDIM>
class PyFunctor : public madness::FunctionFunctorInterface<T, NDIM> {
    py::object py_callable_;

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

    /// Optional: support vectorized evaluation for better performance.
    /// For now we rely on the scalar interface; MADNESS will call operator()
    /// for each quadrature point.
};

#endif  // PYMADNESS_PY_FUNCTOR_H
