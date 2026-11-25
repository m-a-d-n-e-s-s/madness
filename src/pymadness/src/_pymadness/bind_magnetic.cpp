/*
  pymadness - Python bindings for MADNESS

  Magnetic Green's-function kernels (madness/mra/magnetic_kernels.h, see magnetic_kernels.md):
  the same C++ kernels that znemo uses, for Python drivers.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <madness/mra/magnetic_kernels.h>

namespace py = pybind11;
using namespace madness;

void bind_magnetic(py::module_& m) {

    py::enum_<MagneticGauge>(m, "MagneticGauge")
        .value("symmetric", MagneticGauge::symmetric)
        .value("landau", MagneticGauge::landau);

    py::enum_<magnetic::PhaseRefinement>(m, "PhaseRefinement")
        .value("poly_error", magnetic::PhaseRefinement::poly_error)
        .value("box_phase", magnetic::PhaseRefinement::box_phase);

    // keyword construction: verbose; any other keyword sets that parameter (user-defined, e.g. origin=[x, y, z],
    // lo=, nbroaden=, ft_qextra=); the derived ones are set by make_magnetic_kernel
    py::class_<MagneticKernelParameters>(m, "MagneticKernelParameters",
        "Parameters of the magnetic kernels (input block magnetic_kernel); derived ones from B, thresh and the nuclei")
        .def(py::init([](bool verbose, py::kwargs kwargs) {
                 MagneticKernelParameters p;
                 for (auto item : kwargs) {
                     const py::handle v = item.second;
                     p.set_user_value(std::string(py::str(item.first)),
                                      py::isinstance<py::bool_>(v) ? (v.cast<bool>() ? "true" : "false")
                                                                   : std::string(py::str(v)));
                 }
                 p.verbose = verbose;
                 return p;
             }),
             py::arg("verbose") = false)
        .def_property_readonly("origin", [](const MagneticKernelParameters& p) {
                 const coord_3d o = p.origin();
                 return std::array<double, 3>{o[0], o[1], o[2]}; })
        .def_readwrite("verbose", &MagneticKernelParameters::verbose)
        .def("set_derived_values", [](MagneticKernelParameters& p, double B, py::object thresh,
                                      const std::vector<std::tuple<double, double, double, double>>& nuclei,
                                      bool cusp_free) {
                 std::vector<MagneticKernelNucleus> nuc;
                 for (const auto& [Z, x, c, rho0] : nuclei) nuc.push_back({Z, x, c, rho0});
                 p.set_derived_values(B, thresh.is_none() ? FunctionDefaults<3>::get_thresh() : thresh.cast<double>(),
                                      nuc, cusp_free);
             }, py::arg("B"), py::arg("thresh") = py::none(),
             py::arg("nuclei") = std::vector<std::tuple<double, double, double, double>>(), py::arg("cusp_free") = false,
             "derived parameters for the field B; nuclei: list of (Z, x, smoothing radius c, rho(0))")
        .def_property_readonly("lo", &MagneticKernelParameters::lo)
        .def_property_readonly("eps", &MagneticKernelParameters::eps)
        .def_property_readonly("tau_step", &MagneticKernelParameters::tau_step)
        .def_property_readonly("tmin", &MagneticKernelParameters::tmin)
        .def_property_readonly("t_lll", &MagneticKernelParameters::t_lll)
        .def_property_readonly("ft_qmax", &MagneticKernelParameters::ft_qmax)
        .def_property_readonly("ft_rcut", &MagneticKernelParameters::ft_rcut)
        .def_property_readonly("ft_rtaper", &MagneticKernelParameters::ft_rtaper)
        .def_property_readonly("ft_period", &MagneticKernelParameters::ft_period)
        .def("tmin_cusp_free", &MagneticKernelParameters::tmin_cusp_free, py::arg("C"))
        .def("__str__", [](const MagneticKernelParameters& p) { return p.print_to_string(); });

    py::class_<MagneticKernel, std::shared_ptr<MagneticKernel>>(m, "MagneticKernel",
        "(T_k - E)^{-1} for an electron in a field B along z, pi = p + A (magnetic_kernels.md)")
        .def("apply", [](const MagneticKernel& K, const complex_function_3d& f, double E) {
            py::gil_scoped_release release;
            return K.apply(f, E);
        }, py::arg("f"), py::arg("E"), "(T_k - E)^{-1} f, E < bottom()")
        .def("apply", [](const MagneticKernel& K, const std::vector<complex_function_3d>& f,
                         const std::vector<double>& E) {
            py::gil_scoped_release release;
            Tensor<double> t(long(E.size()));
            for (std::size_t i = 0; i < E.size(); ++i) t(long(i)) = E[i];
            return K.apply(f, t);
        }, py::arg("f"), py::arg("E"), "[(T_k - E_i)^{-1} f_i]")
        .def("bottom", &MagneticKernel::bottom, "inf spec T_k = |B|/2")
        .def("B", &MagneticKernel::B, "field along z")
        .def("absorbs_paramagnetic", &MagneticKernel::absorbs_paramagnetic)
        .def("absorbs_diamagnetic", &MagneticKernel::absorbs_diamagnetic)
        .def("gauge", &MagneticKernel::gauge)
        .def("name", &MagneticKernel::name)
        .def("from_symmetric_gauge", [](const MagneticKernel& K, const complex_function_3d& f) {
            py::gil_scoped_release release;
            return K.from_symmetric_gauge(f);
        }, py::arg("f"), "e^{+i chi} f, chi = (B/2) x~ y~ (identity in the symmetric gauge); untruncated")
        .def("to_symmetric_gauge", [](const MagneticKernel& K, const complex_function_3d& f) {
            py::gil_scoped_release release;
            return K.to_symmetric_gauge(f);
        }, py::arg("f"), "e^{-i chi} f; untruncated");

    m.def("make_magnetic_kernel", [](World& world, const std::string& kind, double B, py::object params,
                                     py::object thresh,
                                     const std::vector<std::tuple<double, double, double, double>>& nuclei,
                                     bool cusp_free) {
        MagneticKernelParameters p = params.is_none() ? MagneticKernelParameters() : params.cast<MagneticKernelParameters>();
        std::vector<MagneticKernelNucleus> nuc;
        for (const auto& [Z, x, c, rho0] : nuclei) nuc.push_back({Z, x, c, rho0});
        p.set_derived_values(B, thresh.is_none() ? FunctionDefaults<3>::get_thresh() : thresh.cast<double>(), nuc,
                             cusp_free);
        py::gil_scoped_release release;
        return make_magnetic_kernel(world, kind, B, p);
    }, py::arg("world"), py::arg("kind"), py::arg("B"), py::arg("params") = py::none(), py::arg("thresh") = py::none(),
       py::arg("nuclei") = std::vector<std::tuple<double, double, double, double>>(), py::arg("cusp_free") = false,
       "kernel of the given kind (ft, chain_symmetric, chain_landau, mehler) in the field B along z; the derived\n"
       "parameters of params are set for B, thresh (default FunctionDefaults3D) and the nuclei (Z, x, c, rho(0));\n"
       "FunctionDefaults3D must be set");
    m.def("magnetic_kernel_names", &magnetic_kernel_names);

    m.def("mul_bilinear_phase", [](const complex_function_3d& f, double kappa, double o0, double o1,
                                   magnetic::PhaseRefinement refinement, double tolfac, bool recursive,
                                   double dphi, double rext, double tol, double rcut, double rtaper) {
        py::gil_scoped_release release;
        magnetic::BilinearPhase p(kappa, o0, o1);
        p.refinement = refinement; p.tolfac = tolfac; p.recursive = recursive;
        p.dphi = dphi; p.rext = rext; p.tol = tol; p.rcut = rcut; p.rtaper = rtaper;
        return magnetic::mul_bilinear_phase(f, p);
    }, py::arg("f"), py::arg("kappa"), py::arg("o0") = 0.0, py::arg("o1") = 0.0,
       py::arg("refinement") = magnetic::PhaseRefinement::poly_error, py::arg("tolfac") = 0.1,
       py::arg("recursive") = true, py::arg("dphi") = 3.0, py::arg("rext") = 1e10, py::arg("tol") = 0.0,
       py::arg("rcut") = 0.0, py::arg("rtaper") = 0.0,
       "f e^{i kappa (x0-o0)(x1-o1)} (axes 0, 1), the tree of f refined where the phase is not resolved;\n"
       "untruncated (magnetic_primitives.h)");
}
