/*
  pymadness - Python bindings for MADNESS
  Tensor<T> ↔ numpy array conversion.
*/

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <algorithm>
#include <madness/tensor/tensor.h>

namespace py = pybind11;
using namespace madness;

/// Convert a MADNESS Tensor<double> to a numpy array (copies data).
static py::array_t<double> tensor_to_numpy(const Tensor<double>& t) {
    std::vector<py::ssize_t> shape(t.ndim());
    for (long i = 0; i < t.ndim(); ++i) {
        shape[i] = t.dim(i);
    }
    // Allocate a C-contiguous numpy array
    py::array_t<double> arr(shape);
    double* dst = arr.mutable_data();
    if (t.iscontiguous()) {
        std::copy(t.ptr(), t.ptr() + t.size(), dst);
    } else {
        // Non-contiguous tensor (e.g. from slicing): make a contiguous
        // copy first, then flat-copy into the numpy buffer.
        Tensor<double> tc = madness::copy(t);
        std::copy(tc.ptr(), tc.ptr() + tc.size(), dst);
    }
    return arr;
}

/// Convert a numpy array to a MADNESS Tensor<double> (copies data).
static Tensor<double> numpy_to_tensor(py::array_t<double, py::array::c_style> arr) {
    py::buffer_info info = arr.request();
    std::vector<long> dims(info.ndim);
    for (int i = 0; i < info.ndim; ++i) {
        dims[i] = info.shape[i];
    }
    Tensor<double> t(dims);
    const double* src = static_cast<const double*>(info.ptr);
    double* dst = t.ptr();
    std::copy(src, src + t.size(), dst);
    return t;
}

void bind_tensor(py::module_& m) {
    py::class_<Tensor<double>>(m, "Tensor",
        "Multi-dimensional array (MADNESS Tensor)")
        .def(py::init<>())
        .def("ndim", &Tensor<double>::ndim, "Number of dimensions")
        .def("dim", &Tensor<double>::dim, py::arg("i"), "Size of dimension i")
        .def("size", &Tensor<double>::size, "Total number of elements")
        .def("normf", &Tensor<double>::normf, "Frobenius norm")
        .def("to_numpy", [](const Tensor<double>& t) {
            return tensor_to_numpy(t);
        }, "Convert to numpy array (copies data)")
        .def_static("from_numpy", [](py::array_t<double, py::array::c_style> arr) {
            return numpy_to_tensor(arr);
        }, py::arg("array"), "Create Tensor from numpy array (copies data)")
        .def("__repr__", [](const Tensor<double>& t) {
            std::ostringstream os;
            os << "Tensor(ndim=" << t.ndim() << ", size=" << t.size();
            if (t.ndim() > 0) {
                os << ", shape=[";
                for (long i = 0; i < t.ndim(); ++i) {
                    if (i > 0) os << ",";
                    os << t.dim(i);
                }
                os << "]";
            }
            os << ")";
            return os.str();
        });

    // Conversion helpers exposed at module level
    m.def("tensor_to_numpy", &tensor_to_numpy,
        py::arg("tensor"), "Convert MADNESS Tensor to numpy array");
    m.def("numpy_to_tensor", &numpy_to_tensor,
        py::arg("array"), "Convert numpy array to MADNESS Tensor");
}
