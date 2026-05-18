/*
  pymadness - Python bindings for MADNESS
  FunctionDefaults<NDIM> bindings for dimensions 1, 2, 3.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <madness/mra/mra.h>

namespace py = pybind11;
using namespace madness;

template<std::size_t NDIM>
static void bind_function_defaults(py::module_& m, const char* name) {
    using FD = FunctionDefaults<NDIM>;

    py::class_<FD>(m, name,
        "Static global defaults for numerical parameters of Function objects")
        // Wavelet order
        .def_static("set_k", &FD::set_k,
            py::arg("k"), "Set default wavelet order")
        .def_static("get_k", &FD::get_k,
            "Get default wavelet order")

        // Threshold
        .def_static("set_thresh", &FD::set_thresh,
            py::arg("thresh"), "Set default truncation threshold")
        .def_static("get_thresh", &FD::get_thresh,
            "Get default truncation threshold")

        // Cell / domain
        .def_static("set_cubic_cell", &FD::set_cubic_cell,
            py::arg("lo"), py::arg("hi"),
            "Set cubic simulation cell [-L, L]^NDIM")
        .def_static("get_cell", &FD::get_cell,
            "Get simulation cell as Tensor (NDIM x 2)")
        .def_static("get_cell_width", &FD::get_cell_width,
            "Get cell width along each dimension")
        .def_static("get_cell_volume", &FD::get_cell_volume,
            "Get volume of the simulation cell")

        // Refinement levels
        .def_static("set_initial_level", &FD::set_initial_level,
            py::arg("level"), "Set initial projection level")
        .def_static("get_initial_level", &FD::get_initial_level,
            "Get initial projection level")
        .def_static("set_max_refine_level", &FD::set_max_refine_level,
            py::arg("level"), "Set maximum adaptive refinement level")
        .def_static("get_max_refine_level", &FD::get_max_refine_level,
            "Get maximum adaptive refinement level")

        // Refinement control
        .def_static("set_refine", &FD::set_refine,
            py::arg("value"), "Enable/disable adaptive refinement")
        .def_static("get_refine", &FD::get_refine,
            "Get adaptive refinement flag")
        .def_static("set_autorefine", &FD::set_autorefine,
            py::arg("value"), "Enable/disable auto-refinement")
        .def_static("get_autorefine", &FD::get_autorefine,
            "Get auto-refinement flag")
        .def_static("set_truncate_on_project", &FD::set_truncate_on_project,
            py::arg("value"), "Enable/disable truncation on projection")
        .def_static("get_truncate_on_project", &FD::get_truncate_on_project,
            "Get truncate-on-project flag")
        .def_static("set_truncate_mode", &FD::set_truncate_mode,
            py::arg("mode"), "Set truncation mode (0, 1, or 2)")
        .def_static("get_truncate_mode", &FD::get_truncate_mode,
            "Get truncation mode")

        // Debug
        .def_static("set_debug", &FD::set_debug,
            py::arg("value"), "Enable/disable debug output")
        .def_static("get_debug", &FD::get_debug,
            "Get debug flag")

        // Defaults
        .def_static("set_defaults", &FD::set_defaults,
            py::arg("world"),
            "Reset all defaults (k=7, thresh=1e-5, unit cube)")
        .def_static("print", &FD::print,
            "Print current default settings");
}

void bind_defaults(py::module_& m) {
    bind_function_defaults<1>(m, "FunctionDefaults1D");
    bind_function_defaults<2>(m, "FunctionDefaults2D");
    bind_function_defaults<3>(m, "FunctionDefaults3D");
    bind_function_defaults<4>(m, "FunctionDefaults4D");
    bind_function_defaults<5>(m, "FunctionDefaults5D");
    bind_function_defaults<6>(m, "FunctionDefaults6D");
}
