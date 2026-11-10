/// \file test_magnetic_kernels.cc
/// \brief Consistency of the magnetic Green's-function kernels (magnetic_kernels.h)
///
/// Without arguments: self-consistency tests at B = 2 (see magnetic_kernels.md, "Verification"):
///   - ft vs chain_symmetric on an axial (m = 0), an m = -1 and an off-axis (London phase) Gaussian;
///   - chain_landau conjugated to the symmetric gauge vs chain_symmetric;
///   - mehler vs chain_symmetric on the m = 0 function (Lz = 0: the kernels coincide);
///   - the weak resolvent identity <f|(T_B - E) G f> = <f|f>.
/// With mode=apply: apply one kernel to a saved function and save the result (comparison with the
/// original implementations), e.g.
///   test_magnetic_kernels mode=apply kernel=ft B=10 E=4 k=10 thresh=1e-5 L=20 in=f out=gf

#include <madness/mra/magnetic_kernels.h>
#include <madness/mra/vmra.h>
#include <madness/world/test_utilities.h>

#include <cstdio>
#include <map>
#include <string>

using namespace madness;

namespace {

std::map<std::string, std::string> parse_args(int argc, char** argv) {
    std::map<std::string, std::string> a;
    for (int i = 1; i < argc; ++i) {
        const std::string s(argv[i]);
        const auto eq = s.find('=');
        if (eq != std::string::npos) a[s.substr(0, eq)] = s.substr(eq + 1);
    }
    return a;
}

double getd(const std::map<std::string, std::string>& a, const std::string& k, double d) {
    return a.count(k) ? std::stod(a.at(k)) : d;
}

/// (x~ - i y~)^|m| exp(-alpha |r - R|^2) exp(-i (B/2)(R_x y - R_y x))  (magnetically translated to R)
complex_function_3d gaussian(World& world, const double alpha, const int am, const coord_3d& R, const double B) {
    return complex_factory_3d(world).functor([=](const coord_3d& r) {
        const double dx = r[0] - R[0], dy = r[1] - R[1], dz = r[2] - R[2];
        double_complex v = std::exp(-alpha * (dx * dx + dy * dy + dz * dz))
                         * std::exp(double_complex(0.0, -0.5 * B * (R[0] * r[1] - R[1] * r[0])));
        for (int i = 0; i < am; ++i) v *= double_complex(dx, -dy);
        return v;
    });
}

double rel_diff(const complex_function_3d& a, const complex_function_3d& b) { return (a - b).norm2() / b.norm2(); }

/// <f|(T_B - E) g> in the symmetric gauge about r = 0, by first derivatives (weak form)
double_complex weak_resolvent(World& world, const complex_function_3d& f, const complex_function_3d& g,
                              const double B, const double E) {
    double_complex kin = 0.0;
    for (int d = 0; d < 3; ++d) {
        complex_derivative_3d D(world, d);
        const complex_function_3d Af = complex_factory_3d(world).functor([B, d](const coord_3d& r) {
            const double A[3] = {-0.5 * B * r[1], 0.5 * B * r[0], 0.0};
            return double_complex(A[d], 0.0);
        });
        // pi_d = -i d_d + A_d
        complex_function_3d pf = double_complex(0.0, -1.0) * D(f) + Af * f;
        complex_function_3d pg = double_complex(0.0, -1.0) * D(g) + Af * g;
        kin += 0.5 * inner(pf, pg);
    }
    return kin - E * inner(f, g);
}

int self_test(World& world, const bool with_ft, const int nbroaden) {
    const double B = 2.0, E = 0.0;
    test_output t("magnetic kernels, B = 2");
    t.set_cout_to_terminal();

    MagneticKernelParameters p;
    p.B = B;
    p.nodes = 50;
    p.tmin = 1e-5;
    p.nbroaden = nbroaden;
    auto make = [&](const std::string& kind) {
        MagneticKernelParameters q = p;
        q.kind = kind;
        return make_magnetic_kernel(world, q);
    };
    const auto ft = make("ft"), sym = make("chain_symmetric"), lan = make("chain_landau"), meh = make("mehler");
    const double tol = 10 * FunctionDefaults<3>::get_thresh();

    const complex_function_3d f0 = gaussian(world, 1.0, 0, {0, 0, 0}, B);
    const complex_function_3d f1 = gaussian(world, 1.0, 1, {0, 0, 0}, B);
    const complex_function_3d f2 = gaussian(world, 1.0, 0, {1.0, 0.5, 0.3}, B);

    for (const auto& [name, f] : std::vector<std::pair<std::string, complex_function_3d>>{{"m=0", f0}, {"m=-1", f1}, {"off-axis", f2}}) {
        const complex_function_3d gs = sym->apply(f, E);
        const double dft = with_ft ? rel_diff(ft->apply(f, E), gs) : 0.0;
        const complex_function_3d gl = lan->to_symmetric_gauge(lan->apply(lan->from_symmetric_gauge(f), E));
        const double dl = rel_diff(gl, gs);
        const double_complex lhs = weak_resolvent(world, f, gs, B, E), rhs = inner(f, f);
        const double dres = std::abs(lhs - rhs) / std::abs(rhs);
        if (world.rank() == 0)
            printf("  %-8s |G_sym f| = %.8f  ft-sym %.2e  landau-sym %.2e  resolvent %.2e\n", name.c_str(), gs.norm2(), dft, dl, dres);
        t.checkpoint(dft < tol, name + " ft vs chain_symmetric");
        t.checkpoint(dl < tol, name + " chain_landau vs chain_symmetric");
        t.checkpoint(dres < 100 * tol, name + " weak resolvent identity");
        if (name == "m=0") {
            const complex_function_3d gm = meh->apply(f, E);
            const double dm = rel_diff(gm, gs);
            if (world.rank() == 0) printf("  %-8s mehler-sym %.2e  mehler-landau %.2e  |G_meh f| = %.8f |G_lan f| = %.8f\n",
                                          name.c_str(), dm, rel_diff(gm, gl), gm.norm2(), gl.norm2());
            t.checkpoint(dm < tol, name + " mehler vs chain_symmetric");
        }
    }
    return t.end();
}

int apply_mode(World& world, const std::map<std::string, std::string>& a) {
    MagneticKernelParameters p;
    p.kind = a.at("kernel");
    p.B = getd(a, "B", 0.0);
    p.nodes = int(getd(a, "nodes", p.nodes));
    p.tmin = getd(a, "tmin", p.tmin);
    p.ft_rcut = getd(a, "rcut", p.ft_rcut);
    p.ft_rtaper = getd(a, "rtaper", p.ft_rtaper);
    p.verbose = true;
    const double E = getd(a, "E", 0.0);
    auto K = make_magnetic_kernel(world, p);
    complex_function_3d f;
    load(f, a.at("in"));
    const double w0 = wall_time();
    const complex_function_3d g = K->apply(f, E);
    if (world.rank() == 0)
        printf("apply %s B=%g E=%g: |f| = %.12e  |G f| = %.12e  <f|G f> = %.12e %+.12e i  wall %.1f s\n", p.kind.c_str(),
               p.B, E, f.norm2(), g.norm2(), std::real(inner(f, g)), std::imag(inner(f, g)), wall_time() - w0);
    if (a.count("out")) save(g, a.at("out"));
    return 0;
}

}  // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);
    const auto a = parse_args(argc, argv);
    const bool apply = a.count("mode") && a.at("mode") == "apply";

    const double L = getd(a, "L", apply ? 20.0 : 10.0);
    FunctionDefaults<3>::set_k(int(getd(a, "k", 8)));
    FunctionDefaults<3>::set_thresh(getd(a, "thresh", apply ? 1e-5 : 1e-4));
    FunctionDefaults<3>::set_cubic_cell(-L, L);

    int success = 0;
    {
        success = apply ? apply_mode(world, a) : self_test(world, getd(a, "ft", 1.0) != 0.0, int(getd(a, "nbroaden", 6)));
    }
    world.gop.fence();
    finalize();
    return success;
}
