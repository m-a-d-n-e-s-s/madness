/*
  pymadness - Python bindings for MADNESS
  Copyright (C) 2024 MADNESS developers

  Main pybind11 module definition.
*/

#include <pybind11/pybind11.h>

namespace py = pybind11;

// Forward declarations for sub-module binders
void bind_world(py::module_& m);
void bind_defaults(py::module_& m);
void bind_tensor(py::module_& m);
void bind_function(py::module_& m);
void bind_operators(py::module_& m);

PYBIND11_MODULE(_pymadness, m) {
    m.doc() = "MADNESS: Multiresolution Adaptive Numerical Environment for Scientific Simulation";

    bind_world(m);
    bind_defaults(m);
    bind_tensor(m);
    bind_function(m);
    bind_operators(m);
}
