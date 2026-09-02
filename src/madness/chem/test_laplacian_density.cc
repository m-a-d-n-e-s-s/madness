/// \file test_laplacian_density.cc
/// \brief accuracy of the Laplacian of the density, direct vs. log-density route
///
/// The density Laplacian is needed by the kinetic-energy potential and by any
/// Laplacian-level functional, and it is the worst-conditioned object in the
/// chain: rho has a nuclear cusp, so differentiating it twice resolves a
/// singularity, and the MRA noise of rho is amplified by 2^n per derivative.
///
/// The log-density route avoids differentiating rho itself. With
/// rho = R^2 n and zeta = log(n) (n the regularized, cusp-free nemo density),
///
///   grad rho    = R^2 ( grad n - 2 U1 n ) = rho ( grad zeta - 2 U1 )
///   laplacian rho = div [ rho ( grad zeta - 2 U1 ) ]
///
/// so grad(rho) is built analytically from a smooth zeta and an analytic U1,
/// and only one numerical divergence is taken. U1 = -grad(R)/R.
///
/// The accuracy measure is the round trip through the Poisson equation: since
/// -laplacian(u) = 4 pi f  for  u = f * (1/r),
///
///   rho_rec = -1/(4 pi) * CoulombOperator( laplacian rho )
///
/// must return rho. Two further invariants are checked: int laplacian(rho) = 0 by
/// the divergence theorem, and int rho_rec = number of electrons.
///
/// WHAT THE ROUND TRIP DOES AND DOES NOT MEASURE. (-laplacian)^{-1} has symbol
/// 1/k^2, so it amplifies long-wavelength error and damps short-wavelength error
/// by the same factor. It is therefore
///   * maximally sensitive to a spurious monopole -- an error with nonzero mean
///     becomes a m/(4 pi r) tail in rho_rec, i.e. a charge error;
///   * nearly blind to near-nucleus error, which is where the density Laplacian
///     is actually hard: at depth 21, h ~ 2.5e-5 bohr, so k ~ 4e4 and cusp-region
///     error is suppressed by ~1e-9 before it reaches rho_rec;
///   * bounded below by the Coulomb operator's own eps, which is set to thresh
///     here. With thresh 1e-6 the best routes land at 1-2e-6, i.e. at the floor,
///     so differences among them are not resolved by this measure.
/// It is a conservation check. Measured on HF/LiH at thresh 1e-6:
///
///   route                                  int(lap)    int(rho_rec)   rel
///   make_laplacian_density (3-term, old)    7.3e-04     2.5433        8.8e-04
///   make_laplacian_density (1 div, new)    -6.1e-10     4.0000025     2.0e-06
///   direct: div(grad(rho))                  4.6e-14     4.0000016     1.0e-06
///   log form: div(rho (grad zeta - 2 U1))   4.7e-14     3.9999977     1.2e-05
///
/// The old expansion assembled laplacian(rho) from three terms of norm 9.5, 110
/// and 53 whose exact sum has integral zero, so its absolute error was set by the
/// relative error of the largest -- itself a difference of two nuclear operators
/// carrying the -Z/r singularity. A divergence cannot create a monopole whatever
/// its truncation error, which is why the single-divergence form conserves.
///
/// To rank the surviving routes something else is needed: the exact identity
/// int |grad rho|^2 = -int rho laplacian(rho) (both sides dominated by the
/// near-nucleus region, so sensitive where the round trip is blind), a weak-form
/// scan int v laplacian(rho) + int grad(v).grad(rho) = 0 over smooth v of varying
/// width (which lets one choose the length scale probed), the Kato cusp
/// coefficient -4 Z rho(0) of the 1/r divergence, or a tighter-threshold
/// self-reference. None of those is implemented yet.

#include <madness.h>
#include <madness/chem/nemo.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/write_test_input.h>

using namespace madness;

/// pointwise log with a floor

/// only grad(zeta) is used, so the additive constant is irrelevant
struct logme {
    typedef double resultT;
    double floor_ = 1.e-14;
    logme() = default;
    explicit logme(double f) : floor_(f) {}
    Tensor<double> operator()(const Key<3>& key, const Tensor<double>& t) const {
        Tensor<double> r = copy(t);
        double* p = r.ptr();
        for (long i = 0; i < r.size(); ++i) p[i] = std::log(std::max(floor_, p[i]));
        return r;
    }
    template<typename Archive> void serialize(Archive& ar) {}
};

/// laplacian by two successive numerical derivatives of the argument
real_function_3d laplacian_direct(World& world, const real_function_3d& f) {
    real_function_3d fr = copy(f).refine();
    vecfuncT parts(3);
    for (int axis = 0; axis < 3; ++axis) {
        real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
        real_function_3d df = D(fr).refine();
        parts[axis] = D(df);
    }
    return sum(world, parts, true).truncate();
}

/// laplacian from a genuine second-derivative operator

/// grad_ble_two / grad_bpsline_two return the three second derivatives directly,
/// so the shipped grad -> smoothen -> diff chain is not needed. `smoothen`
/// projects to k-1 and back and averages with the original, i.e. it throws away
/// the highest polynomial mode of the first derivative before the second one is
/// taken -- exactly the mode the second derivative needs most.
real_function_3d laplacian_second(World& world, const real_function_3d& f,
                                  const std::string& variant) {
    real_function_3d fc = copy(f);          // grad_*_two refines in place
    vecfuncT d2;
    if (variant == "ble2") d2 = grad_ble_two(fc, true);
    else if (variant == "bspline2") d2 = grad_bpsline_two(fc, true);
    else MADNESS_EXCEPTION("unknown second-derivative variant", 1);
    return sum(world, d2, true).truncate();
}

/// the shipped chain: first derivative, smoothen, first derivative again
real_function_3d laplacian_smoothen(World& world, const real_function_3d& f) {
    real_function_3d fr = copy(f).refine();
    vecfuncT parts(3);
    for (int axis = 0; axis < 3; ++axis) {
        real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
        real_function_3d df = D(fr).refine();
        Nemo::smoothen(df);
        parts[axis] = D(df);
    }
    return sum(world, parts, true).truncate();
}

/// round trip: rho_rec = -1/(4 pi) * (1/r) * laplacian(rho)
real_function_3d reconstruct(World& world, const real_function_3d& lap,
                             double lo, double eps) {
    real_convolution_3d poisson = CoulombOperator(world, lo, eps);
    return (-1.0 / (4.0 * constants::pi)) * poisson(lap);
}

void report(World& world, const std::string& name, const real_function_3d& lap,
            const real_function_3d& rho, double nelectron, double lo, double eps) {
    const real_function_3d rec = reconstruct(world, lap, lo, eps);
    const real_function_3d diff = (rec - rho).truncate();
    const double err = diff.norm2();
    const double rel = err / rho.norm2();
    if (world.rank() == 0)
        printf("%-14s  tree %8lu  depth %3lu | int(lap) %11.3e | "
               "int(rec) %12.8f (N=%.1f) | ||rec-rho|| %10.3e  rel %10.3e\n",
               name.c_str(), (unsigned long) lap.tree_size(),
               (unsigned long) lap.max_depth(), lap.trace(), rec.trace(),
               nelectron, err, rel);
}

int test_laplacian(World& world) {

    FunctionDefaults<3>::set_thresh(1.e-6);
    FunctionDefaults<3>::set_k(8);
    const double thresh = FunctionDefaults<3>::get_thresh();
    const double lo = 1.e-6;

    // converged Hartree-Fock LiH through the nemo route
    CalculationParameters param;
    param.set_user_defined_value<std::vector<double>>("protocol", {1.e-4, 1.e-6});
    param.set_user_defined_value("k", 8);
    // dconv 1e-5 is not reachable here: the energy settles at 1e-8 while the BSH
    // residual plateaus, so value() reports failure and returns 0 with perfectly
    // good orbitals. Loosen it rather than iterate to maxiter for nothing.
    param.set_user_defined_value("econv", 1.e-6);
    param.set_user_defined_value("dconv", 1.e-4);
    param.set_user_defined_value("maxiter", 20);
    write_test_input test_input(param);          // LiH by default
    commandlineparser parser;
    parser.set_keyval("input", test_input.filename());

    Nemo nemo(world, parser);
    const double energy = nemo.value();
    // value() returns 0.0 when the convergence criteria are not all met; the
    // orbitals are still the last (converged-in-energy) set, so report and go on
    if (world.rank() == 0) {
        print("\nHF/LiH nemo energy", energy);
        print("nuclear correlation factor", nemo.get_nemo_param().ncf().first,
              nemo.get_nemo_param().ncf().second);
    }

    // the total densities. make_density returns the alpha density, hence the 2
    const real_function_3d n =
            (2.0 * nemo.make_density(nemo.get_calc()->get_aocc(),
                                     nemo.get_calc()->get_amo())).truncate();
    const real_function_3d R2 = nemo.R_square;
    const real_function_3d rho = (R2 * n).truncate();
    const double nelectron = rho.trace();
    if (world.rank() == 0) {
        print("regularized density  tree", n.tree_size(), " depth", n.max_depth());
        print("physical    density  tree", rho.tree_size(), " depth", rho.max_depth());
        print("number of electrons  ", nelectron);
    }

    // zeta = log(n) and its gradient; n is cusp-free, so this is the smooth object
    const real_function_3d zeta = unary_op(n, logme());
    vecfuncT gz = grad(zeta, true);
    truncate(world, gz);
    const vecfuncT U1 = nemo.ncf->U1vec();
    if (world.rank() == 0)
        print("zeta                 tree", zeta.tree_size(), " depth", zeta.max_depth());

    if (world.rank() == 0) print("");

    // ---- first, one level down: how well does the log form reproduce grad(rho)?
    //      grad rho = rho (grad zeta - 2 U1) must equal the numerical gradient of
    //      rho. If this already disagrees, the loss is in the log form; if it
    //      agrees, the loss is in the divergence that follows.
    {
        real_function_3d rho_refined = copy(rho).refine();
        vecfuncT grad_direct(3), grad_log(3);
        for (int axis = 0; axis < 3; ++axis) {
            real_derivative_3d D = free_space_derivative<double, 3>(world, axis);
            grad_direct[axis] = D(rho_refined).truncate();
            grad_log[axis] = (rho * (gz[axis] - 2.0 * U1[axis])).truncate();
        }
        if (world.rank() == 0) print("grad(rho): direct vs log form");
        for (int axis = 0; axis < 3; ++axis) {
            const double nd = grad_direct[axis].norm2();
            const double d = (grad_log[axis] - grad_direct[axis]).norm2();
            if (world.rank() == 0)
                printf("  axis %d  ||grad_direct|| %10.4e  ||log-direct|| %10.3e"
                       "  rel %10.3e  trees %lu / %lu\n", axis, nd, d, d / nd,
                       (unsigned long) grad_direct[axis].tree_size(),
                       (unsigned long) grad_log[axis].tree_size());
        }
        if (world.rank() == 0) print("");
    }

    // ---- route 1: two numerical derivatives of the *physical* density
    report(world, "direct(rho)", laplacian_direct(world, rho), rho, nelectron, lo, thresh);

    // ---- route 2: two numerical derivatives of the *regularized* density,
    //      with the R^2 part handled analytically (Nemo's own route)
    report(world, "nemo", nemo.make_laplacian_density(n), rho, nelectron, lo, thresh);

    // ---- route 2b: grad(rho) from make_ddensity, then one divergence.
    //      Same regularization as route 2 but assembled as a single divergence.
    {
        vecfuncT drho(3);
        for (int axis = 0; axis < 3; ++axis)
            drho[axis] = nemo.make_ddensity(n, axis);
        truncate(world, drho);
        report(world, "ddens-flux", div(drho, true), rho, nelectron, lo, thresh);
    }

    // ---- route 3: grad(rho) from the log form, then one divergence
    //      grad rho = rho (grad zeta - 2 U1)
    {
        vecfuncT flux(3);
        for (int axis = 0; axis < 3; ++axis)
            flux[axis] = (rho * (gz[axis] - 2.0 * U1[axis])).truncate();
        report(world, "log-flux", div(flux, true), rho, nelectron, lo, thresh);
    }

    // ---- route 4: Nemo's skeleton, but the density Laplacian from the log form
    //      laplacian(n) = div(n grad zeta) instead of sum_a D_a D_a n
    {
        NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nemo.ncf.get());
        const real_function_3d U1dot =
                real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();
        const Nuclear<double, 3> U_op(world, nemo.ncf);
        const Nuclear<double, 3> V_op(world, nemo.get_calc().get());

        vecfuncT flux(3);
        for (int axis = 0; axis < 3; ++axis)
            flux[axis] = (n * gz[axis]).truncate();

        vecfuncT terms(3);
        terms[0] = (2.0 * U1dot * n).truncate();
        terms[1] = (-4.0 * (U_op(n) - V_op(n))).truncate();
        terms[2] = div(flux, true);
        const real_function_3d inner_sum = sum(world, terms, true);
        report(world, "nemo+log", (R2 * inner_sum).truncate(), rho, nelectron, lo, thresh);
    }

    // ---- WHICH int(Delta rho) IS REAL?
    //      Algorithm A'' inside make_xc_potential reports int(Delta rho) ~ -39 for
    //      LiH; make_laplacian_density reports ~ -6e-10 on the same molecule by
    //      what should be the same recipe. Ten orders apart, so one of them is
    //      wrong. Compute every variant here on ONE density, in one run, so the
    //      difference has to be the recipe and not the context.
    {
        auto div_with = [&](vecfuncT v, bool freespace) {
            v = copy(world, v);                      // refine() mutates in place
            reconstruct(world, v);
            refine(world, v);
            std::vector<std::shared_ptr<real_derivative_3d> > D(3);
            for (int a = 0; a < 3; ++a)
                D[a] = freespace
                       ? std::shared_ptr<real_derivative_3d>(new real_derivative_3d(
                             free_space_derivative<double, 3>(world, a)))
                       : std::shared_ptr<real_derivative_3d>(new real_derivative_3d(world, a));
            vecfuncT dv(3);
            for (int a = 0; a < 3; ++a) dv[a] = apply(*D[a], v[a], false);
            world.gop.fence();
            return sum(world, dv, true);
        };

        // (i) grad(rho) exactly as A'' builds it: default-BC derivative of a
        //     refined n, then truncate
        vecfuncT g_a2(3);
        {
            real_function_3d n_ref = copy(n).refine();
            for (int a = 0; a < 3; ++a) {
                real_derivative_3d D(world, a);
                g_a2[a] = (R2 * (D(n_ref) - 2.0 * U1[a] * n)).truncate();
            }
        }
        // (ii) grad(rho) from make_ddensity (free-space derivative internally)
        vecfuncT g_md(3);
        for (int a = 0; a < 3; ++a) g_md[a] = nemo.make_ddensity(n, a);
        truncate(world, g_md);

        // (iii) A''-style grad but built from a DEEPLY refined n, mimicking what
        //       refine_to_common_level does to the xc_args before A'' runs
        vecfuncT g_deep(3);
        {
            real_function_3d n_deep = copy(n);
            n_deep.refine(); n_deep.refine(); n_deep.refine();
            for (int a = 0; a < 3; ++a) {
                real_derivative_3d D(world, a);
                g_deep[a] = (R2 * (D(n_deep) - 2.0 * U1[a] * n_deep)).truncate();
            }
        }

        if (world.rank() == 0) print("\nLAPCMP  int(Delta rho), one density, five recipes (exact value 0)");
        auto row = [&](const char* tag, const real_function_3d& f) {
            const double t = f.trace(), nn = f.norm2();
            if (world.rank() == 0)
                printf("  %-42s %14.6e   norm %11.4e  tree %8lu\n", tag, t, nn,
                       (unsigned long) f.tree_size());
        };
        row("make_laplacian_density", nemo.make_laplacian_density(n));
        row("A''grad + default-BC div", div_with(g_a2, false));
        row("A''grad + free-space div", div_with(g_a2, true));
        row("make_ddensity + default-BC div", div_with(g_md, false));
        row("make_ddensity + free-space div", div_with(g_md, true));
        row("A''grad from 3x-refined n + default-BC div", div_with(g_deep, false));
        // and the boundary flux, which is what a nonzero integral must come from
        if (world.rank() == 0) {
            const double L = FunctionDefaults<3>::get_cell_width().max() * 0.5;
            double f_a2 = 0.0, f_md = 0.0;
            const int m = 5;
            for (int i = 0; i < m; ++i) for (int j = 0; j < m; ++j) {
                const double y = -L + 2*L*i/(m-1), z = -L + 2*L*j/(m-1);
                f_a2 += std::abs(g_a2[0](coord_3d{ L, y, z})) + std::abs(g_a2[0](coord_3d{-L, y, z}));
                f_md += std::abs(g_md[0](coord_3d{ L, y, z})) + std::abs(g_md[0](coord_3d{-L, y, z}));
            }
            printf("  |grad(rho)_x| sampled on the x-walls: A''grad %10.3e   make_ddensity %10.3e\n",
                   f_a2/(2*m*m), f_md/(2*m*m));
        }
    }

    // ---- why the nemo route leaks charge: it assembles laplacian(rho) as a sum
    //      of three individually large terms that must cancel to give
    //      int laplacian(rho) = 0. The absolute accuracy of the result is then set
    //      by the relative accuracy of the largest term, not by thresh. The
    //      log-flux route computes one divergence instead, and a divergence
    //      integrates to zero by construction.
    {
        NuclearCorrelationFactor::U1_dot_U1_functor u1_dot_u1(nemo.ncf.get());
        const real_function_3d U1dot =
                real_factory_3d(world).functor(u1_dot_u1).truncate_on_project();
        const Nuclear<double, 3> U_op(world, nemo.ncf);
        const Nuclear<double, 3> V_op(world, nemo.get_calc().get());

        const real_function_3d t1 = (2.0 * U1dot * n).truncate();
        const real_function_3d t2 = (-4.0 * (U_op(n) - V_op(n))).truncate();

        if (world.rank() == 0) print("\nthe three terms of the nemo route, weighted by R^2:");
        auto term = [&](const char* nm, const real_function_3d& t) {
            const real_function_3d w = (R2 * t).truncate();
            if (world.rank() == 0)
                printf("  %-22s int %14.6e   norm %12.4e\n", nm, w.trace(), w.norm2());
            return w.trace();
        };
        const double i1 = term("2 U1^2 n", t1);
        const double i2 = term("-4 (U-V) n", t2);

        // the density Laplacian, five ways
        std::vector<std::pair<std::string, real_function_3d>> lapn;
        lapn.emplace_back("abgv+smoothen (shipped)", laplacian_smoothen(world, n));
        lapn.emplace_back("abgv plain", laplacian_direct(world, n));
        lapn.emplace_back("ble2", laplacian_second(world, n, "ble2"));
        lapn.emplace_back("bspline2", laplacian_second(world, n, "bspline2"));
        {
            vecfuncT f2(3);
            for (int axis = 0; axis < 3; ++axis) f2[axis] = (n * gz[axis]).truncate();
            lapn.emplace_back("log form div(n grad z)", div(f2, true));
        }

        if (world.rank() == 0)
            print("\nnemo route with laplacian(n) computed five ways "
                  "(int should be 0, N should be", nelectron, "):");
        for (auto& v : lapn) {
            const double i3 = term(("laplacian(n): " + v.first).c_str(), v.second);
            vecfuncT terms(3);
            terms[0] = t1; terms[1] = t2; terms[2] = v.second;
            const real_function_3d L = (R2 * sum(world, terms, true)).truncate();
            report(world, "  -> " + v.first, L, rho, nelectron, lo, thresh);
            if (world.rank() == 0)
                printf("     (term integrals %.4e + %.4e + %.4e = %.4e)\n",
                       i1, i2, i3, i1 + i2 + i3);
        }
    }

    // ---- and the direct route with a true second-derivative operator
    if (world.rank() == 0) print("\ndirect laplacian of the physical density, "
                                 "second-derivative operators:");
    report(world, "direct ble2", laplacian_second(world, rho, "ble2"), rho, nelectron, lo, thresh);
    report(world, "direct bspline2", laplacian_second(world, rho, "bspline2"), rho, nelectron, lo, thresh);

    // ---- where the error sits: profile along the internuclear axis (Li at 0,
    //      H at 1.4375 bohr). The cusps are the interesting points.
    if (world.rank() == 0) print("\nprofile along the Li-H axis (y=z=0)\n");
    const real_function_3d L_direct = laplacian_direct(world, rho);
    const real_function_3d L_log = [&]() {
        vecfuncT flux(3);
        for (int axis = 0; axis < 3; ++axis)
            flux[axis] = (rho * (gz[axis] - 2.0 * U1[axis])).truncate();
        return div(flux, true);
    }();
    const real_function_3d rec_direct = reconstruct(world, L_direct, lo, thresh);
    const real_function_3d rec_log = reconstruct(world, L_log, lo, thresh);

    const std::vector<double> xs = {-0.5, -0.2, -0.05, 0.0, 0.05, 0.2, 0.5,
                                    1.0, 1.35, 1.4375, 1.5, 2.0, 3.0};
    if (world.rank() == 0)
        printf("%9s %14s %14s %14s %12s %12s\n", "x", "rho", "rec(direct)",
               "rec(log)", "err(direct)", "err(log)");
    for (double x : xs) {
        const coord_3d r{x, 0.0, 0.0};
        const double v = rho(r), a = rec_direct(r), b = rec_log(r);
        if (world.rank() == 0)
            printf("%9.4f %14.6e %14.6e %14.6e %12.2e %12.2e\n",
                   x, v, a, b, a - v, b - v);
    }

    return 0;
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world, argc, argv);
    std::cout.precision(6);

    int result = 0;
    try {
        result = test_laplacian(world);
    } catch (std::exception& e) {
        print("caught an exception:", e.what());
        result = 1;
    }
    world.gop.fence();
    madness::finalize();
    return result;
}
