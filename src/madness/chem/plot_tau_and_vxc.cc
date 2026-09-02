/// \file plot_tau_and_vxc.cc
/// \brief tau and the meta-GGA xc potential on LiH/TPSS, four recipes each
///
/// The question is whether tau can be made smooth. With a nuclear correlation
/// factor psi = R F, and the product rule
///
///   |grad psi|^2 = R^2 ( |grad F|^2 - 2 F U1.grad F + |U1|^2 F^2 ),   U1 = -grad(R)/R
///
/// only the cusp-free F is differentiated numerically -- that much the shipped
/// code already does. What it also does is carry U1 and |U1|^2 as MRA Functions
/// and multiply. U1_x ~ x/r is non-smooth at each nucleus componentwise, so as a
/// Function it needs depth ~18, and the products inherit that depth whatever the
/// depth of their other factor. refine_to_common_level then imposes it on every
/// xc intermediate, and the 1/r pole of div(X) gets resolved into pointwise
/// excursions of order 1e5 in the potential. That is the mechanism of
/// XC-NEMO-MEMORY-INVESTIGATION.md 9.2/9.8, and phase 1 of its 11.0 plan is to
/// never form those products: evaluate U1 from its analytic functor at the
/// quadrature points of the (shallow) orbital tree instead.
///
/// This program measures what that buys, on ONE density, in one run, so that a
/// difference between two curves is the tau recipe and nothing else.
///
///   tau  (i)   U1 projected into MRA        TauU1::mra        -- shipped
///        (ii)  U1 evaluated pointwise       TauU1::pointwise  -- the proposal
///        (iii) no U1 at all                 TauU1::none       -- DIAGNOSTIC, wrong tau
///        (iv)  the moldft recipe            psi = R F differentiated directly,
///                                           no regularization anywhere
///
///   v_xc (i)   from tau (i)
///        (ii)  from tau (iii), i.e. the wrong tau without U1
///        (iii) from tau = 0
///        (iv)  from tau (iv), the moldft recipe
///
/// (ii) and (iv) are the two that must agree with (i) to within the discretization
/// if the implementation is right: they are three routes to the same physical
/// quantity. (iii) and v_xc(iii) are controls, not candidates.
///
/// A caveat on v_xc(iii): make_libxc_args floors tau from below at the von
/// Weizsaecker bound rho*chi/8, exactly, so tau = 0 does not reach libxc as zero.
/// It reaches it as tau_W. For a one-orbital system that IS the exact tau; for LiH
/// it is a lower bound, and the curve should be read as "tau replaced by its vW
/// bound", not "tau removed".
///
/// Output, all on identical sample points so the columns can be subtracted:
///   tau_line.dat  vxc_line.dat   along the Li-H axis, log-refined at both nuclei
///   tau_plane.dat vxc_plane.dat  the z=0 plane, gnuplot pm3d blocks
/// plus a census (tree size, depth, norm, pointwise min/max, integral) per curve.
///
/// Run on one rank: the sampling goes through Function::operator()(coord), which
/// is a broadcast per point.

#include <madness.h>
#include <madness/chem/nemo.h>
#include <madness/chem/SCFOperators.h>
#include <madness/chem/xcfunctional.h>
#include <madness/chem/correlationfactor.h>
#include <madness/chem/write_test_input.h>

#include <algorithm>
#include <fstream>
#include <sstream>
#include <iomanip>

using namespace madness;

namespace {

/// pointwise log with a floor -- the operation prep_xc_args applies to build zeta

/// Reproduced here rather than reached for, so the sigma comparison uses exactly
/// what the default zeta path uses.
struct logme {
    typedef double resultT;
    double floor_ = 1.e-14;
    Tensor<double> operator()(const Key<3>& key, const Tensor<double>& t) const {
        Tensor<double> r = copy(t);
        double* p = r.ptr();
        for (long i = 0; i < r.size(); ++i) p[i] = std::log(std::max(floor_, p[i]));
        return r;
    }
    template<typename Archive> void serialize(Archive& ar) {}
};

/// pointwise min/max of a Function, which Function itself does not carry
struct minmax_op {
    typedef double resultT;
    std::shared_ptr<std::pair<double, double> > acc;
    std::shared_ptr<Mutex> mtx;
    minmax_op() : acc(new std::pair<double, double>(1.e300, -1.e300)), mtx(new Mutex()) {}
    Tensor<double> operator()(const Key<3>& key, const Tensor<double>& t) const {
        double lo = 1.e300, hi = -1.e300;
        const double* p = t.ptr();
        for (long i = 0; i < t.size(); ++i) {
            lo = std::min(lo, p[i]);
            hi = std::max(hi, p[i]);
        }
        ScopedMutex<Mutex> lock(mtx.get());
        acc->first = std::min(acc->first, lo);
        acc->second = std::max(acc->second, hi);
        return copy(t);
    }
    template<typename Archive> void serialize(Archive& ar) {}
};

/// one census line: the four numbers that decide whether a function is cheap and smooth
void census(World& world, const std::string& tag, const real_function_3d& f) {
    if (not f.is_initialized()) {
        if (world.rank() == 0) printf("  %-26s  (not initialized)\n", tag.c_str());
        return;
    }
    real_function_3d g = copy(f);
    g.reconstruct();
    minmax_op mm;
    real_function_3d dummy = unary_op(g, mm);
    dummy.clear();
    world.gop.min(mm.acc->first);
    world.gop.max(mm.acc->second);
    if (world.rank() == 0)
        printf("  %-26s  norm %12.5e  int %13.6e  tree %8lu  depth %3lu  GB %7.4f"
               "  min %13.5e  max %13.5e\n",
               tag.c_str(), g.norm2(), g.trace(), (unsigned long) g.tree_size(),
               (unsigned long) g.max_depth(), double(g.size()) * 8.e-9,
               mm.acc->first, mm.acc->second);
}

/// abscissae along the Li-H axis: linear everywhere, logarithmic at both nuclei

/// the excursions live within ~1e-3 bohr of a nucleus, which a linear scan at any
/// affordable spacing steps straight over
std::vector<double> axis_abscissae(const std::vector<double>& nuclei) {
    std::vector<double> x;
    const int nlin = 400;
    for (int i = 0; i <= nlin; ++i) x.push_back(-2.0 + 5.0 * double(i) / nlin);
    for (double c : nuclei) {
        const int nlog = 80;
        for (int i = 0; i <= nlog; ++i) {
            const double s = std::pow(10.0, -6.0 + 6.3 * double(i) / nlog);
            x.push_back(c + s);
            x.push_back(c - s);
        }
        x.push_back(c);
    }
    std::sort(x.begin(), x.end());
    x.erase(std::unique(x.begin(), x.end()), x.end());
    return x;
}

/// sample every function on the same points and write one file, columns in order
void line_dump(World& world, const std::string& fname,
               const std::vector<std::string>& labels,
               const std::vector<real_function_3d>& fs,
               const std::vector<double>& xs) {
    MADNESS_CHECK(labels.size() == fs.size());
    std::vector<std::vector<double> > vals(fs.size(), std::vector<double>(xs.size()));
    for (size_t j = 0; j < fs.size(); ++j)
        for (size_t i = 0; i < xs.size(); ++i)
            vals[j][i] = fs[j](coord_3d{xs[i], 0.0, 0.0});
    if (world.rank() != 0) return;
    std::ofstream of(fname);
    of << "# column 1: x (bohr), along the Li-H axis, y=z=0\n";
    for (size_t j = 0; j < labels.size(); ++j)
        of << "# column " << j + 2 << ": " << labels[j] << "\n";
    of << std::scientific << std::setprecision(8);
    for (size_t i = 0; i < xs.size(); ++i) {
        of << std::setw(18) << xs[i];
        for (size_t j = 0; j < fs.size(); ++j) of << std::setw(18) << vals[j][i];
        of << "\n";
    }
    print("wrote", fname, "  ", xs.size(), "points x", fs.size(), "curves");
}

/// sample MRA functions and analytic functors on one set of points

/// The functor columns are evaluated directly at each sample coordinate, never
/// through an MRA function. That is the whole point: the pointwise route contracts
/// U1 against the smooth pieces at quadrature points without projecting anything,
/// so to see what it actually computes the plot has to do the same -- take the
/// smooth pieces from MRA, take U1 from its functor, and multiply on the grid.
/// A column of projected U1 would reintroduce exactly the representation error the
/// route avoids.
void line_dump_mixed(World& world, const std::string& fname,
                     const std::vector<std::string>& flabels,
                     const std::vector<real_function_3d>& fs,
                     const std::vector<std::string>& glabels,
                     const std::vector<std::shared_ptr<FunctionFunctorInterface<double,3> > >& gs,
                     const std::vector<double>& xs) {
    MADNESS_CHECK(flabels.size() == fs.size() and glabels.size() == gs.size());
    const size_t nf = fs.size(), ng = gs.size();
    std::vector<std::vector<double> > vals(nf + ng, std::vector<double>(xs.size()));
    for (size_t j = 0; j < nf; ++j)
        for (size_t i = 0; i < xs.size(); ++i)
            vals[j][i] = fs[j](coord_3d{xs[i], 0.0, 0.0});
    for (size_t j = 0; j < ng; ++j)
        for (size_t i = 0; i < xs.size(); ++i)
            vals[nf + j][i] = (*gs[j])(coord_3d{xs[i], 0.0, 0.0});
    if (world.rank() != 0) return;
    std::ofstream of(fname);
    of << "# column 1: x (bohr), along the Li-H axis, y=z=0\n";
    for (size_t j = 0; j < nf; ++j) of << "# column " << j + 2 << ": " << flabels[j] << "   [MRA]\n";
    for (size_t j = 0; j < ng; ++j) of << "# column " << nf + j + 2 << ": " << glabels[j] << "   [functor, pointwise]\n";
    of << std::scientific << std::setprecision(10);
    for (size_t i = 0; i < xs.size(); ++i) {
        of << std::setw(20) << xs[i];
        for (size_t j = 0; j < nf + ng; ++j) of << std::setw(20) << vals[j][i];
        of << "\n";
    }
    print("wrote", fname, "  ", xs.size(), "points x", nf, "MRA +", ng, "functor columns");
}

/// the z=0 plane on a regular grid, as gnuplot pm3d blocks (blank line per row)
void plane_dump(World& world, const std::string& fname,
                const std::vector<std::string>& labels,
                const std::vector<real_function_3d>& fs,
                const double xlo, const double xhi, const double ylo, const double yhi,
                const int npt) {
    MADNESS_CHECK(labels.size() == fs.size());
    std::vector<std::vector<double> > vals(fs.size(), std::vector<double>(size_t(npt) * npt));
    for (size_t j = 0; j < fs.size(); ++j) {
        size_t p = 0;
        for (int ix = 0; ix < npt; ++ix) {
            const double x = xlo + (xhi - xlo) * ix / (npt - 1);
            for (int iy = 0; iy < npt; ++iy, ++p) {
                const double y = ylo + (yhi - ylo) * iy / (npt - 1);
                vals[j][p] = fs[j](coord_3d{x, y, 0.0});
            }
        }
    }
    if (world.rank() != 0) return;
    std::ofstream of(fname);
    of << "# columns 1,2: x y (bohr), z=0\n";
    for (size_t j = 0; j < labels.size(); ++j)
        of << "# column " << j + 3 << ": " << labels[j] << "\n";
    of << std::scientific << std::setprecision(8);
    size_t p = 0;
    for (int ix = 0; ix < npt; ++ix) {
        const double x = xlo + (xhi - xlo) * ix / (npt - 1);
        for (int iy = 0; iy < npt; ++iy, ++p) {
            const double y = ylo + (yhi - ylo) * iy / (npt - 1);
            of << std::setw(18) << x << std::setw(18) << y;
            for (size_t j = 0; j < fs.size(); ++j) of << std::setw(18) << vals[j][p];
            of << "\n";
        }
        of << "\n";
    }
    print("wrote", fname, "  ", npt, "x", npt, "points x", fs.size(), "curves");
}

/// an environment knob, so one binary can drive the whole matrix of runs
std::string env_or(const char* name, const std::string& fallback) {
    const char* e = std::getenv(name);
    return (e == nullptr or std::string(e).empty()) ? fallback : std::string(e);
}

/// "1e-4,1e-6,1e-8" -> {1e-4, 1e-6, 1e-8}
std::vector<double> parse_protocol(const std::string& s) {
    std::vector<double> p;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) if (not tok.empty()) p.push_back(std::atof(tok.c_str()));
    MADNESS_CHECK_THROW(not p.empty(), "TAUPLOT_PROTOCOL parsed to nothing");
    return p;
}

}   // anonymous namespace


int run(World& world) {

    // LiH/TPSS on the nemo path. Two protocol rungs: the point is the shape of
    // tau near a nucleus, which is already fully developed at 1e-6, and the 1e-8
    // rung is what makes the divergence form unaffordable in the first place.
    // TAUPLOT_PROTOCOL  comma-separated thresholds, default "1e-4,1e-6"
    // TAUPLOT_DERIV     abgv | ble | bspline, default abgv
    // TAUPLOT_K         polynomial order, default 8
    // TAUPLOT_TAG       prefix on every output file, so runs do not overwrite
    // TAUPLOT_NOVXC     set to skip the potential section entirely
    const std::vector<double> protocol = parse_protocol(env_or("TAUPLOT_PROTOCOL", "1e-4,1e-6"));
    const std::string deriv = env_or("TAUPLOT_DERIV", "abgv");
    const int korder = std::atoi(env_or("TAUPLOT_K", "8").c_str());
    const std::string tag = env_or("TAUPLOT_TAG", "");
    const bool novxc = (std::getenv("TAUPLOT_NOVXC") != nullptr);

    CalculationParameters param;
    param.set_user_defined_value<std::vector<double> >("protocol", protocol);
    param.set_user_defined_value("k", korder);
    param.set_user_defined_value("xc", std::string("tpss"));
    param.set_user_defined_value("dft_deriv", deriv);
    param.set_user_defined_value("econv", 1.e-6);
    param.set_user_defined_value("dconv", 1.e-4);
    param.set_user_defined_value("maxiter", 25);
    write_test_input test_input(param);          // LiH: Li at 0, H at 1.4375 bohr
    commandlineparser parser;
    parser.set_keyval("input", test_input.filename());

    Nemo nemo(world, parser);
    const double energy = nemo.value();
    if (world.rank() == 0) {
        print("\nRUN  protocol", protocol, " dft_deriv", deriv, " k", korder,
              " tag", tag.empty() ? "(none)" : tag, novxc ? " [tau only]" : "");
        print("\nTPSS/LiH nemo energy", energy);
        print("nuclear correlation factor", nemo.get_nemo_param().ncf().first,
              nemo.get_nemo_param().ncf().second);
        print("thresh", FunctionDefaults<3>::get_thresh(), " k", FunctionDefaults<3>::get_k());
    }
    MADNESS_CHECK_THROW(bool(nemo.ncf), "no nuclear correlation factor -- nothing to compare");

    const vecfuncT& amo = nemo.get_calc()->amo;
    const Tensor<double>& aocc = nemo.get_calc()->aocc;
    const real_function_3d R = nemo.ncf->function();
    const double T_ref = nemo.compute_kinetic_energy(amo);

    // ---------------------------------------------------------------- tau, 4 ways
    //
    // (i)-(iii) go through set_tau on a freshly built operator each time, so that
    // no run inherits a tree that a previous one refined. (iv) is assembled here:
    // it is moldft's recipe, psi = R F differentiated directly with the cusp in
    // it and no U1 anywhere -- the independent route that the nemo_he_tpss deck
    // cross-checks against.
    auto tau_from = [&](const TauU1 mode) {
        XCOperator<double, 3> xcop(world, &nemo, 0);
        xcop.set_tau(amo, aocc, vecfuncT(), Tensor<double>(), mode);
        return copy(xcop.get_tau(0));
    };

    const real_function_3d tau_mra = tau_from(TauU1::mra);
    const real_function_3d tau_pw = tau_from(TauU1::pointwise);
    const real_function_3d tau_nou1 = tau_from(TauU1::none);

    const real_function_3d tau_moldft = [&]() {
        vecfuncT psi = mul(world, R, amo);
        truncate(world, psi);
        real_function_3d t = real_factory_3d(world).compressed();
        for (int axis = 0; axis < 3; ++axis) {
            real_derivative_3d D(world, axis);
            vecfuncT p = copy(world, psi);
            refine(world, p);
            vecfuncT dp = apply(world, D, p);
            // aocc is 1 per spin channel for a closed shell, matching set_tau's
            // weighting; fold it in the same way rather than assuming
            for (size_t i = 0; i < dp.size(); ++i)
                if (aocc(long(i)) != 1.0) dp[i] = std::sqrt(aocc(long(i))) * dp[i];
            t += dot(world, dp, dp);
        }
        return (0.5 * t).truncate();
    }();


    if (world.rank() == 0) {
        print("\nKINETIC ENERGY. T from nemo's analytic decomposition is", T_ref,
              "; for a closed shell 2*int(tau_alpha) must reproduce it.");
        print("TAU CENSUS  (2*int is the check; min/max is the question)");
    }
    census(world, "(i)   U1 as Functions", tau_mra);
    // NB on the pointwise route tau is not a Function: get_tau() assembles one by
    // projection for display, which is the very thing that route avoids. This row
    // is therefore the mra route by another name, and is kept only so the tau
    // plots stay comparable across implementations.
    census(world, "(ii)  pw, reconstructed", tau_pw);
    census(world, "(iii) no U1 [wrong]", tau_nou1);
    census(world, "(iv)  moldft recipe", tau_moldft);
    if (world.rank() == 0) {
        printf("  T check:  2*int tau_mra %14.8f   2*int tau_pw %14.8f"
               "   2*int tau_moldft %14.8f   T_ref %14.8f\n",
               2.0 * tau_mra.trace(), 2.0 * tau_pw.trace(),
               2.0 * tau_moldft.trace(), T_ref);
        printf("  ||tau_pw - tau_mra||     %12.5e   (rel %10.3e)\n",
               (tau_pw - tau_mra).norm2(), (tau_pw - tau_mra).norm2() / tau_mra.norm2());
        printf("  ||tau_moldft - tau_mra|| %12.5e   (rel %10.3e)\n",
               (tau_moldft - tau_mra).norm2(), (tau_moldft - tau_mra).norm2() / tau_mra.norm2());
    }

    // ------------------------------------------------- the pointwise op's inputs
    //
    // G = sum_i w_i F_i grad(F_i) = 1/2 grad(n) is what u1_tau_terms contracts
    // against U1 at fixed quadrature points, with no adaptive error control. If
    // the 1e-2 bohr anomaly is a sampling failure, it has to be visible in G, in
    // U1, or in their product -- so dump all of them on the same abscissae as tau.
    // Rebuilt here rather than plumbed out of set_tau: same orbitals, same
    // occupations, same derivative operator, so it is the same object.
    vecfuncT G(3);
    real_function_3d n_reg, gradf;
    {
        vecfuncT wmo;
        for (size_t i = 0; i < amo.size(); ++i) {
            const double w = aocc(long(i));
            if (w == 0.0) continue;
            wmo.push_back(w == 1.0 ? amo[i] : std::sqrt(w) * amo[i]);
        }
        n_reg = dot(world, wmo, wmo);
        gradf = real_factory_3d(world).compressed();
        for (int axis = 0; axis < 3; ++axis) {
            real_derivative_3d Dax(world, axis);
            if (deriv == "bspline") Dax.set_bspline1();
            else if (deriv == "ble") Dax.set_ble1();
            vecfuncT mo_copy = copy(world, wmo);
            refine(world, mo_copy);
            vecfuncT dmo = apply(world, Dax, mo_copy);
            G[axis] = dot(world, mo_copy, dmo);
            gradf += dot(world, dmo, dmo);
        }
        truncate(world, G);
        gradf.truncate();
    }
    const real_function_3d U1x = nemo.ncf->U1(0);
    NuclearCorrelationFactor::U1_dot_U1_functor u1dotu1(nemo.ncf.get());
    const real_function_3d U1sq =
            real_factory_3d(world).functor(u1dotu1).truncate_on_project();
    // the two U1 terms as the *shipped* route forms them: projected products
    const real_function_3d cross_mra = (-2.0 * U1x * G[0]).truncate();
    const real_function_3d diag_mra  = (U1sq * n_reg).truncate();

    if (world.rank() == 0) print("\nINPUTS TO THE POINTWISE OP");
    census(world, "n = sum w F^2", n_reg);
    census(world, "gradf = sum w |grad F|^2", gradf);
    census(world, "G_x = sum w F dF/dx", G[0]);
    census(world, "R^2", nemo.ncf->square());
    census(world, "U1_x   [projected]", U1x);
    census(world, "|U1|^2 [projected]", U1sq);
    census(world, "-2 U1_x G_x  [projected]", cross_mra);
    census(world, "|U1|^2 n     [projected]", diag_mra);

    // ---------------------------------------------------------------- tau plots
    //
    // Written here, before the potentials, on purpose: at a tight protocol the
    // multiplicative v_xc build is the expensive object (the divergence form is
    // what the memory investigation is about), so if this run dies it dies after
    // the tau record is already on disk.
    const std::vector<double> nuclei = {0.0, 1.4375};
    const std::vector<double> xs = axis_abscissae(nuclei);

    // Route order and labels are identical in every output file of this program,
    // so a given column means the same thing in tau_line, vxc_line, dfdsigma_line
    // and dfdtau_line.
    line_dump(world, tag + "tau_line.dat",
              {"U1 as shipped", "no U1 [wrong]", "U1 pointwise", "moldft"},
              {tau_mra, tau_nou1, tau_pw, tau_moldft}, xs);

    // The pointwise contraction, laid out so it can be rebuilt on gnuplot's grid:
    //   tau_pw = 1/2 R^2 ( gradf - 2(U1x Gx + U1y Gy + U1z Gz) + |U1|^2 n )
    // Smooth pieces from MRA, U1 straight from its functor at the same sample
    // points, product formed by the plotting layer. That is exactly what
    // assemble_nemo_tau does inside make_libxc_args, only at these points instead
    // of the quadrature points -- so this shows what the functional is fed,
    // without any Function of a U1 product existing anywhere.

    // The pointwise contraction, laid out so it can be rebuilt on gnuplot's grid:
    //   tau_pw = 1/2 R^2 ( gradf - 2(U1x Gx + U1y Gy + U1z Gz) + |U1|^2 n )
    // Smooth pieces from MRA, U1 straight from its functor at the same sample
    // points, product formed by the plotting layer. That is exactly what
    // assemble_nemo_tau does inside make_libxc_args, only at these points instead
    // of the quadrature points -- so this shows what the functional is fed,
    // without any Function of a U1 product existing anywhere.
    {
        typedef FunctionFunctorInterface<double, 3> functorT;
        std::vector<std::shared_ptr<functorT> > u1f;
        for (int axis = 0; axis < 3; ++axis)
            u1f.push_back(std::shared_ptr<functorT>(
                    new NuclearCorrelationFactor::U1_functor(nemo.ncf.get(), axis)));
        u1f.push_back(std::shared_ptr<functorT>(
                new NuclearCorrelationFactor::U1_dot_U1_functor(nemo.ncf.get())));

        // rho_alpha and zeta exactly as prep_xc_args builds them on the DEFAULT
        // path: zeta = grad(log rho), a numerical derivative of a kinked function.
        // Dumped alongside the regularized pieces so sigma can be formed both ways
        // on one grid -- rho^2 (zeta.zeta) as libxc gets it today, against
        // 4 R^4 [ |G|^2 - 2 n U1.G + n^2 |U1|^2 ] from the split form.
        const real_function_3d rho_a = (nemo.ncf->square() * n_reg).truncate();
        const real_function_3d logdens = unary_op(rho_a, logme());
        vecfuncT zeta = grad(logdens);
        truncate(world, zeta);

        line_dump_mixed(world, tag + "pointwise_terms.dat",
                        {"n = sum w F^2", "gradf = sum w |grad F|^2",
                         "G_x", "G_y", "G_z", "R^2",
                         "rho_alpha", "zeta_x", "zeta_y", "zeta_z"},
                        {n_reg, gradf, G[0], G[1], G[2], nemo.ncf->square(),
                         rho_a, zeta[0], zeta[1], zeta[2]},
                        {"U1_x", "U1_y", "U1_z", "|U1|^2"}, u1f, xs);
    }

    // ------------------------------------------------ the unprojected probe
    //
    // Everything vxc returns normally is an MRA function: multi_to_multi_op_values
    // evaluates libxc at the quadrature points and then calls values2coeffs on the
    // refine_to_common_level tree, with no refinement and no error test. So a plot of
    // df/dsigma shows the *projection* of df/dsigma, and the two are not the same.
    //
    // XCfunctional::vxc, though, is a pure function of value tensors. Feed it the
    // arguments evaluated at the sample points -- shaped (npt,1,1), since it takes
    // its result dims from t[0] -- and it returns libxc's own values at those points,
    // with no Function anywhere in the chain. Dumped next to the projected version
    // from the same run, so the difference is the representation and nothing else.
    {
        XCfunctional xcf;
        xcf.initialize(nemo.get_calc()->param.xc(),
                       not nemo.get_calc()->param.spin_restricted(), world);

        const long npt = xs.size();
        const long dims[3] = {npt, 1, 1};
        auto mk = [&](const real_function_3d& f) {
            Tensor<double> v(3L, dims);
            for (long i = 0; i < npt; ++i) v.ptr()[i] = f(coord_3d{xs[i], 0.0, 0.0});
            return v;
        };
        auto mkf = [&](const FunctionFunctorInterface<double,3>& g) {
            Tensor<double> v(3L, dims);
            for (long i = 0; i < npt; ++i) v.ptr()[i] = g(coord_3d{xs[i], 0.0, 0.0});
            return v;
        };

        std::vector<Tensor<double> > t(XCfunctional::number_xc_args);
        const real_function_3d rho_a2 = (nemo.ncf->square() * n_reg).truncate();
        t[XCfunctional::enum_rhoa] = mk(rho_a2);
        {   // zeta, the default route's ingredient, so the fallback path also works
            const real_function_3d ld = unary_op(rho_a2, logme());
            vecfuncT z = grad(ld); truncate(world, z);
            t[XCfunctional::enum_zetaa_x] = mk(z[0]);
            t[XCfunctional::enum_zetaa_y] = mk(z[1]);
            t[XCfunctional::enum_zetaa_z] = mk(z[2]);
        }
        t[XCfunctional::enum_nemo_R2] = mk(nemo.ncf->square());
        t[XCfunctional::enum_na]      = mk(n_reg);
        t[XCfunctional::enum_gradfa]  = mk(gradf);
        t[XCfunctional::enum_Ga_x]    = mk(G[0]);
        t[XCfunctional::enum_Ga_y]    = mk(G[1]);
        t[XCfunctional::enum_Ga_z]    = mk(G[2]);
        for (int axis = 0; axis < 3; ++axis)
            t[XCfunctional::enum_u1_x + axis] =
                    mkf(NuclearCorrelationFactor::U1_functor(nemo.ncf.get(), axis));
        t[XCfunctional::enum_u1sq] =
                mkf(NuclearCorrelationFactor::U1_dot_U1_functor(nemo.ncf.get()));

        const std::vector<Tensor<double> > res = xcf.vxc(t, 0);
        const int islot = xc_dfds_slot(xcf.is_spin_polarized(), xcf.needs_tau(),
                                       xcf.needs_sigma());
        const int tslot = xcf.is_spin_polarized() ? 7 : 4;

        if (world.rank() == 0) {
            std::ofstream of(tag + "unprojected.dat");
            // Only the unprojected pair here. The projected counterpart is column 4
            // of dfdsigma_line.dat / dfdtau_line.dat, on these same abscissae, so the
            // two are compared in the plotting layer -- and this block stays
            // independent of whether the potentials were built at all.
            of << "# column 1: x (bohr)\n"
               << "# column 2: df/dsigma, unprojected (libxc values at these points)\n"
               << "# column 3: df/dtau, unprojected\n";
            of << std::scientific << std::setprecision(10);
            for (long i = 0; i < npt; ++i) {
                of << std::setw(20) << xs[i]
                   << std::setw(20) << (islot >= 0 ? res[islot].ptr()[i] : 0.0)
                   << std::setw(20) << (xcf.needs_tau() ? res[tslot].ptr()[i] : 0.0)
                   << "\n";
            }
            print("wrote", tag + "unprojected.dat", " slot", islot, "npt", npt);
        }
    }

    line_dump(world, tag + "grad_line.dat",
              {"n = sum w F^2", "G_x = sum w F dF/dx", "G_y", "G_z",
               "U1_x", "|U1|^2", "-2 U1_x G_x [projected]", "|U1|^2 n [projected]",
               "R^2"},
              {n_reg, G[0], G[1], G[2], U1x, U1sq, cross_mra, diag_mra, nemo.ncf->square()},
              xs);

    plane_dump(world, tag + "tau_plane.dat",
               {"tau (i) U1 as Functions", "tau (ii) U1 pointwise",
                "tau (iii) no U1 [wrong]", "tau (iv) moldft recipe"},
               {tau_mra, tau_pw, tau_nou1, tau_moldft},
               -1.5, 3.0, -2.0, 2.0, 161);

    if (novxc) {
        if (world.rank() == 0) print("\nTAUPLOT_NOVXC set -- stopping before the potentials");
        return 0;
    }

    // ------------------------------------------------- the four routes, four ways
    //
    // One pass per route. Each builds a fresh XCOperator -- make_xc_potential()
    // calls refine_to_common_level, which mutates the intermediates in place, so a
    // second call on the same operator would compare histories rather than taus.
    //
    // Two ways in, and the difference matters. Driving set_tau by *mode* lets
    // make_xc_potential build the U1 functors, so the pointwise route is exercised
    // end to end: tau is assembled from values inside make_libxc_args and never
    // becomes a Function. Injecting a ready-made tau Function is the only way to ask
    // for moldft's tau, but it cannot test the pointwise route -- injecting
    // get_tau()'s projected reconstruction just re-measures the mra route under
    // another name.
    struct route {
        std::string label;
        bool by_mode;
        TauU1 mode;                  // used when by_mode
        real_function_3d tau_in;     // used otherwise
        real_function_3d tau, vxc, dfds, dfdt;
    };
    std::vector<route> routes;
    routes.push_back({"U1 as shipped", true,  TauU1::mra,       {}, {}, {}, {}, {}});
    routes.push_back({"no U1 [wrong]", true,  TauU1::none,      {}, {}, {}, {}, {}});
    routes.push_back({"U1 pointwise",  true,  TauU1::pointwise, {}, {}, {}, {}, {}});
    routes.push_back({"moldft",        false, TauU1::mra, tau_moldft, {}, {}, {}, {}});

    for (auto& r : routes) {
        XCOperator<double, 3> xcop(world, &nemo, 0);
        if (r.by_mode) xcop.set_tau(amo, aocc, vecfuncT(), Tensor<double>(), r.mode);
        else           xcop.set_tau(r.tau_in);
        r.vxc  = xcop.make_xc_potential();      // must precede get_vtau/get_dfdsigma
        r.dfdt = xcop.get_vtau();
        r.dfds = xcop.get_dfdsigma();
        // On the pointwise route get_tau() reconstructs by projection -- see the
        // header. Kept for the absolute-tau column, flagged in the doc.
        r.tau  = copy(xcop.get_tau(0));
    }

    auto report = [&](const char* what, real_function_3d route::*m) {
        if (world.rank() == 0) print("\n" + std::string(what) + " CENSUS");
        for (auto& r : routes) census(world, r.label, r.*m);
    };
    report("TAU", &route::tau);
    report("V_XC", &route::vxc);
    report("DF/DSIGMA", &route::dfds);
    report("DF/DTAU", &route::dfdt);

    if (world.rank() == 0) {
        print("\nT check: 2*int(tau) against T =", T_ref);
        for (auto& r : routes)
            if (r.tau.is_initialized())
                printf("  %-16s %16.8f   err %12.3e\n", r.label.c_str(),
                       2.0 * r.tau.trace(), 2.0 * r.tau.trace() - T_ref);
    }

    // ---------------------------------------------------- one file per quantity
    //
    // Column order is the route order above and is the same in every file, so the
    // four quantities can be read side by side for a given route.
    {
        std::vector<std::string> labels;
        for (auto& r : routes) labels.push_back(r.label);
        auto cols = [&](real_function_3d route::*m) {
            std::vector<real_function_3d> v;
            for (auto& r : routes) v.push_back(r.*m);
            return v;
        };
        // tau_line.dat is deliberately not rewritten here: the copy written before
        // the potentials is the one that survives a run that dies in the expensive
        // half, and it holds the same values.
        line_dump(world, tag + "vxc_line.dat",      labels, cols(&route::vxc),  xs);
        line_dump(world, tag + "dfdsigma_line.dat", labels, cols(&route::dfds), xs);
        line_dump(world, tag + "dfdtau_line.dat",   labels, cols(&route::dfdt), xs);

        plane_dump(world, tag + "vxc_plane.dat", labels, cols(&route::vxc),
                   -1.5, 3.0, -2.0, 2.0, 161);
    }

    return 0;
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world, argc, argv);
    std::cout.precision(8);

    int result = 0;
    try {
        result = run(world);
    } catch (std::exception& e) {
        print("caught an exception:", e.what());
        result = 1;
    }
    world.gop.fence();
    madness::finalize();
    return result;
}
