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

    // all fields read-write, keyword construction with the C++ defaults
    const MagneticKernelParameters d;
    py::class_<MagneticKernelParameters>(m, "MagneticKernelParameters",
        "Parameters of the magnetic kernels; the defaults are those of the original implementations")
        .def(py::init([](const std::string& kind, double B, double ox, double oy, bool verbose, int nodes,
                         double tmin, double eps, double nlll, double max_btau, double phase_tolfac,
                         bool phase_recursive, int nbroaden, double ft_h, double ft_tmin, double ft_tmax_factor,
                         double ft_period, double ft_qextra, double ft_rcut, double ft_rtaper, double ft_dphi,
                         double ft_fq_tol, long ft_qbatch) {
                 MagneticKernelParameters p;
                 p.kind = kind; p.B = B; p.ox = ox; p.oy = oy; p.verbose = verbose; p.nodes = nodes;
                 p.tmin = tmin; p.eps = eps; p.nlll = nlll; p.max_btau = max_btau; p.phase_tolfac = phase_tolfac;
                 p.phase_recursive = phase_recursive; p.nbroaden = nbroaden; p.ft_h = ft_h; p.ft_tmin = ft_tmin;
                 p.ft_tmax_factor = ft_tmax_factor; p.ft_period = ft_period; p.ft_qextra = ft_qextra;
                 p.ft_rcut = ft_rcut; p.ft_rtaper = ft_rtaper; p.ft_dphi = ft_dphi; p.ft_fq_tol = ft_fq_tol;
                 p.ft_qbatch = ft_qbatch;
                 return p;
             }),
             py::arg("kind") = d.kind, py::arg("B") = d.B, py::arg("ox") = d.ox, py::arg("oy") = d.oy,
             py::arg("verbose") = d.verbose, py::arg("nodes") = d.nodes, py::arg("tmin") = d.tmin,
             py::arg("eps") = d.eps, py::arg("nlll") = d.nlll, py::arg("max_btau") = d.max_btau,
             py::arg("phase_tolfac") = d.phase_tolfac, py::arg("phase_recursive") = d.phase_recursive,
             py::arg("nbroaden") = d.nbroaden, py::arg("ft_h") = d.ft_h, py::arg("ft_tmin") = d.ft_tmin,
             py::arg("ft_tmax_factor") = d.ft_tmax_factor, py::arg("ft_period") = d.ft_period,
             py::arg("ft_qextra") = d.ft_qextra, py::arg("ft_rcut") = d.ft_rcut, py::arg("ft_rtaper") = d.ft_rtaper,
             py::arg("ft_dphi") = d.ft_dphi, py::arg("ft_fq_tol") = d.ft_fq_tol, py::arg("ft_qbatch") = d.ft_qbatch)
        .def_readwrite("kind", &MagneticKernelParameters::kind)
        .def_readwrite("B", &MagneticKernelParameters::B)
        .def_readwrite("ox", &MagneticKernelParameters::ox)
        .def_readwrite("oy", &MagneticKernelParameters::oy)
        .def_readwrite("verbose", &MagneticKernelParameters::verbose)
        .def_readwrite("nodes", &MagneticKernelParameters::nodes)
        .def_readwrite("tmin", &MagneticKernelParameters::tmin)
        .def_readwrite("eps", &MagneticKernelParameters::eps)
        .def_readwrite("nlll", &MagneticKernelParameters::nlll)
        .def_readwrite("max_btau", &MagneticKernelParameters::max_btau)
        .def_readwrite("phase_tolfac", &MagneticKernelParameters::phase_tolfac)
        .def_readwrite("phase_recursive", &MagneticKernelParameters::phase_recursive)
        .def_readwrite("nbroaden", &MagneticKernelParameters::nbroaden)
        .def_readwrite("ft_h", &MagneticKernelParameters::ft_h)
        .def_readwrite("ft_tmin", &MagneticKernelParameters::ft_tmin)
        .def_readwrite("ft_tmax_factor", &MagneticKernelParameters::ft_tmax_factor)
        .def_readwrite("ft_period", &MagneticKernelParameters::ft_period)
        .def_readwrite("ft_qextra", &MagneticKernelParameters::ft_qextra)
        .def_readwrite("ft_rcut", &MagneticKernelParameters::ft_rcut)
        .def_readwrite("ft_rtaper", &MagneticKernelParameters::ft_rtaper)
        .def_readwrite("ft_dphi", &MagneticKernelParameters::ft_dphi)
        .def_readwrite("ft_fq_tol", &MagneticKernelParameters::ft_fq_tol)
        .def_readwrite("ft_qbatch", &MagneticKernelParameters::ft_qbatch);

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

    m.def("make_magnetic_kernel", [](World& world, const MagneticKernelParameters& p) {
        py::gil_scoped_release release;
        return make_magnetic_kernel(world, p);
    }, py::arg("world"), py::arg("params"),
       "kernel of kind params.kind (ft, chain_symmetric, chain_landau, mehler); FunctionDefaults3D must be set");
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
