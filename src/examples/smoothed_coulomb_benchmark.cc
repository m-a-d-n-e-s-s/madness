/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
*/

/*!
  \file examples/smoothed_coulomb_benchmark.cc
  \brief Absolute-accuracy benchmark for the smoothed-nuclear-density Coulomb energy.

  The total Coulomb energy of a molecule can be evaluated two ways:

    conventional   E = E_nn + E_ne + E_ee
    smoothed       E = 1/2 <s+e; s+e> + <n-s; e> - sqrt(a/2pi) sum_i Z_i^2

  where s is a Gaussian-smeared nuclear charge density.  The two are algebraically
  identical, so comparing them against each other can only demonstrate agreement, never
  accuracy.  Measuring accuracy needs an independent exact reference.

  This driver supplies one.  It replaces the SCF electron density with a MODEL density
  that is a fixed sum of normalized Gaussians on the nuclei, normalized to -Z_i on each
  atom.  Every integral appearing in either formulation then has a closed form, so the
  exact value of the total AND of each individual term is known to machine precision, at
  any system size.  Subtracting the numerically computed value from it gives the absolute
  error directly.

  Because there is no SCF there is also no convergence criterion, no orbital truncation,
  and no variational cancellation between kinetic and Coulomb error -- all of which
  otherwise mask the quantity of interest.  What remains is exactly the error of
  evaluating the Coulomb energy from a given density, which is what the smoothed
  formulation is intended to improve.

  Usage:
      smoothed_coulomb_benchmark input                 # moldft-style geometry block
      smoothed_coulomb_benchmark input --thresh=1e-6
      smoothed_coulomb_benchmark input --thresh=1e-4,1e-6,1e-8 --k=8
      smoothed_coulomb_benchmark input --opeps=1e-5    # sweep the operator accuracy (E11)
      smoothed_coulomb_benchmark input --a=4           # exponent of the smoothed density

  A third, "consistent" route is printed below each row: the conventional assembly with the
  nuclear potential and the nuclear repulsion taken from the same operator as E_ee,
  Vnuc~ = -G s + (Vnuc + u(s)) and E_nn~ = 1/2 <s;G s> - sqrt(a/2pi) sum Z^2 + E_ov.
  Algebraically it equals the smoothed route; numerically it applies G to s and rho
  separately instead of to s - rho.

  Sign convention: all densities below are POSITIVE and the pairing is
      P(f,g) = int int f(r) g(r') / |r-r'|,
  so with n, s positive nuclear densities and rho the positive electron density,
      E_ne = -P(n,rho),  E_ee = 1/2 P(rho,rho),  <s+e;s+e> = P(s-rho, s-rho).
*/

#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
#include <madness/mra/commandlineparser.h>
#include <madness/chem/molecule.h>
#include <madness/chem/potentialmanager.h>
#include <madness/constants.h>

#include <cmath>
#include <cstdio>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;

static const double PI = madness::constants::pi;

//==============================================================================
// Closed forms
//==============================================================================

/// Coulomb interaction of two unit-charge normalized Gaussians of exponents a, b
/// whose centres are separated by R:   erf(sqrt(ab/(a+b)) R) / R
static double gauss_gauss(double a, double b, double R) {
    const double t = std::sqrt(a * b / (a + b));
    if (R < 1e-12) return 2.0 * t / std::sqrt(PI);   // coincident-centre limit
    return std::erf(t * R) / R;
}

/// Coulomb interaction of a unit point charge with a unit-charge normalized Gaussian
/// of exponent b at separation R.  This is the a -> infinity limit of gauss_gauss.
static double point_gauss(double b, double R) {
    if (R < 1e-12) return 2.0 * std::sqrt(b) / std::sqrt(PI);
    return std::erf(std::sqrt(b) * R) / R;
}

//==============================================================================
// The model electron density: a sum of normalized Gaussians, positive, charge Z_i per atom
//==============================================================================

struct GaussianDensity {
    std::vector<coord_3d> centre;   ///< one entry per primitive
    std::vector<double> expnt;      ///< exponent of each primitive
    std::vector<double> coeff;      ///< charge carried by each primitive (positive)

    size_t size() const { return expnt.size(); }
    double total_charge() const {
        double q = 0.0;
        for (double c : coeff) q += c;
        return q;
    }
};

/// Build a model density with `nprim` Gaussians per atom.  Exponents scale as Z^2, which
/// is the hydrogenic scaling, so the tightest primitive is comparable to (in fact tighter
/// than) the smoothed nuclear density itself and exercises the same refinement.
/// Charges are chosen so each atom carries exactly -Z_i, i.e. the system is locally neutral.
static GaussianDensity make_model_density(const Molecule& mol,
                                          const std::vector<double>& zeta,
                                          const std::vector<double>& weight) {
    MADNESS_CHECK(zeta.size() == weight.size());
    double wsum = 0.0;
    for (double w : weight) wsum += w;
    MADNESS_CHECK(wsum > 0.0);

    GaussianDensity g;
    for (size_t i = 0; i < mol.natom(); ++i) {
        const Atom& a = mol.get_atom(i);
        const double Z = static_cast<double>(a.atomic_number);
        for (size_t k = 0; k < zeta.size(); ++k) {
            g.centre.push_back(coord_3d{a.x, a.y, a.z});
            g.expnt.push_back(zeta[k] * Z * Z);
            g.coeff.push_back(weight[k] / wsum * Z);
        }
    }
    return g;
}

class GaussianDensityFunctor : public FunctionFunctorInterface<double, 3> {
    const GaussianDensity& g;
    const std::vector<coord_3d> pts;
public:
    explicit GaussianDensityFunctor(const GaussianDensity& g, std::vector<coord_3d> pts)
        : g(g), pts(std::move(pts)) {}

    double operator()(const coord_3d& x) const override {
        double sum = 0.0;
        for (size_t p = 0; p < g.size(); ++p) {
            const double dx = x[0] - g.centre[p][0];
            const double dy = x[1] - g.centre[p][1];
            const double dz = x[2] - g.centre[p][2];
            const double rsq = dx * dx + dy * dy + dz * dz;
            const double arg = g.expnt[p] * rsq;
            if (arg > 46.0) continue;   // exp(-46) ~ 1e-20, below any threshold we use
            sum += g.coeff[p] * std::pow(g.expnt[p] / PI, 1.5) * std::exp(-arg);
        }
        return sum;
    }

    std::vector<coord_3d> special_points() const override { return pts; }
};

//==============================================================================
// Exact values of every quantity in both formulations
//==============================================================================

struct ExactEnergies {
    double P_ss = 0.0;     ///< P(s,s), smoothed nuclear with itself (includes self-energy)
    double P_srho = 0.0;   ///< P(s,rho)
    double P_rhorho = 0.0; ///< P(rho,rho)
    double P_nrho = 0.0;   ///< P(n,rho), point nuclei with the model density
    double E_nn = 0.0;     ///< sum_{i<j} Z_i Z_j / R_ij
    double E_ne = 0.0;     ///< -P(n,rho)
    double E_ee = 0.0;     ///< 1/2 P(rho,rho)
    double nrep_corr = 0.0;///< sqrt(a/2pi) sum_i Z_i^2
    double E_ov = 0.0;     ///< sum_{i<j} Z_i Z_j erfc(sqrt(a/2) R_ij)/R_ij, the overlap of the Gaussians
    double e_approx = 0.0; ///< 1/2 P(s-rho, s-rho)
    double e_n_corr = 0.0; ///< -P(n,rho) + P(s,rho)
    double E_conv = 0.0;   ///< E_nn + E_ne + E_ee
    double E_smooth = 0.0; ///< e_approx + e_n_corr - nrep_corr + E_ov
};

static double distance(const coord_3d& a, const coord_3d& b) {
    const double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

static ExactEnergies compute_exact(const Molecule& mol, const GaussianDensity& g, double avalue) {
    ExactEnergies x;
    const size_t na = mol.natom();

    std::vector<coord_3d> R(na);
    std::vector<double> Z(na);
    for (size_t i = 0; i < na; ++i) {
        const Atom& a = mol.get_atom(i);
        R[i] = coord_3d{a.x, a.y, a.z};
        Z[i] = static_cast<double>(a.atomic_number);
    }

    // P(s,s): every pair of smoothed nuclei, including i==j (the self-energy)
    for (size_t i = 0; i < na; ++i)
        for (size_t j = 0; j < na; ++j)
            x.P_ss += Z[i] * Z[j] * gauss_gauss(avalue, avalue, distance(R[i], R[j]));

    // P(s,rho) and P(n,rho)
    for (size_t i = 0; i < na; ++i) {
        for (size_t p = 0; p < g.size(); ++p) {
            const double d = distance(R[i], g.centre[p]);
            x.P_srho += Z[i] * g.coeff[p] * gauss_gauss(avalue, g.expnt[p], d);
            x.P_nrho += Z[i] * g.coeff[p] * point_gauss(g.expnt[p], d);
        }
    }

    // P(rho,rho)
    for (size_t p = 0; p < g.size(); ++p)
        for (size_t q = 0; q < g.size(); ++q)
            x.P_rhorho += g.coeff[p] * g.coeff[q]
                        * gauss_gauss(g.expnt[p], g.expnt[q], distance(g.centre[p], g.centre[q]));

    // conventional decomposition
    for (size_t i = 0; i < na; ++i)
        for (size_t j = i + 1; j < na; ++j)
            x.E_nn += Z[i] * Z[j] / distance(R[i], R[j]);
    x.E_ne = -x.P_nrho;
    x.E_ee = 0.5 * x.P_rhorho;
    x.E_conv = x.E_nn + x.E_ne + x.E_ee;

    // smoothed decomposition
    double zsq = 0.0;
    for (size_t i = 0; i < na; ++i) zsq += Z[i] * Z[i];
    x.nrep_corr = std::sqrt(avalue / (2.0 * PI)) * zsq;
    x.e_approx = 0.5 * (x.P_ss - 2.0 * x.P_srho + x.P_rhorho);
    x.e_n_corr = -x.P_nrho + x.P_srho;
    for (size_t i = 0; i < na; ++i)
        for (size_t j = i + 1; j < na; ++j) {
            const double r = distance(R[i], R[j]);
            x.E_ov += Z[i] * Z[j] * std::erfc(std::sqrt(0.5 * avalue) * r) / r;
        }
    x.E_smooth = x.e_approx + x.e_n_corr - x.nrep_corr + x.E_ov;

    return x;
}

//==============================================================================

static std::vector<double> parse_list(const std::string& s) {
    std::vector<double> v;
    std::stringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ',')) {
        if (!tok.empty()) v.push_back(std::atof(tok.c_str()));
    }
    return v;
}

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    startup(world, argc, argv, true);
    std::cout.precision(10);

    commandlineparser parser(argc, argv);

    // ---- options -----------------------------------------------------------
    std::vector<double> threshes = {1e-4, 1e-6, 1e-8};
    if (parser.key_exists("thresh")) threshes = parse_list(parser.value("thresh"));

    int kwave = -1;                       // -1 => let set_defaults pick from thresh
    if (parser.key_exists("k")) kwave = std::atoi(parser.value("k").c_str());

    double lo = 1e-6;
    if (parser.key_exists("lo")) lo = std::atof(parser.value("lo").c_str());

    // Coulomb operator construction accuracy.  moldft hardcodes 0.001*thresh
    // (SCF.h). Exposing it here is what makes the E11 operator sweep possible:
    // a loose operator is what shelved this project in 2024.
    double opeps = -1.0;                  // negative => use 0.001*thresh
    if (parser.key_exists("opeps")) opeps = std::atof(parser.value("opeps").c_str());

    // the Gaussian exponent of the smoothed nuclear density, both sides
    double avalue = 42.0;
    if (parser.key_exists("a")) avalue = std::atof(parser.value("a").c_str());

    // model density: exponents zeta_k * Z^2, weights normalized to sum to 1
    std::vector<double> zeta = {1.0, 0.1, 0.01};
    std::vector<double> weight = {0.3, 0.3, 0.4};
    if (parser.key_exists("zeta")) zeta = parse_list(parser.value("zeta"));
    if (parser.key_exists("weight")) weight = parse_list(parser.value("weight"));

    // ---- molecule ----------------------------------------------------------
    Molecule molecule(world, parser);
    const size_t na = molecule.natom();

    // Simulation box.  Same rule CalculationParameters uses (bounding_cube + 50), so the
    // box matches what moldft would have chosen for the same geometry.  Getting this wrong
    // is silent and catastrophic: with the default unit cell an atom at the origin sits on
    // a corner and only 1/8 of its density is inside the box.
    double L = molecule.bounding_cube() + 50.0;
    if (parser.key_exists("l")) L = std::atof(parser.value("l").c_str());

    const GaussianDensity gdens = make_model_density(molecule, zeta, weight);
    const ExactEnergies exact = compute_exact(molecule, gdens, avalue);

    if (world.rank() == 0) {
        printf("\n");
        printf("  Smoothed-density Coulomb benchmark against an analytic reference\n");
        printf("  ================================================================\n\n");
        printf("  atoms                  %zu\n", na);
        printf("  box  [-L,L]            %.3f\n", L);
        printf("  primitives             %zu (%zu per atom)\n", gdens.size(), zeta.size());
        printf("  model charge           %.12f  (should be sum Z_i)\n", gdens.total_charge());
        printf("  smoothing exponent a   %.4f\n", avalue);
        printf("  exponents (zeta*Z^2)  ");
        for (double z : zeta) printf(" %.4g", z);
        printf("\n  weights               ");
        for (double w : weight) printf(" %.4g", w);
        printf("\n\n");
        printf("  EXACT values (closed form)\n");
        printf("    E_nn                 %20.12f\n", exact.E_nn);
        printf("    E_ne                 %20.12f\n", exact.E_ne);
        printf("    E_ee                 %20.12f\n", exact.E_ee);
        printf("    E_conv               %20.12f\n", exact.E_conv);
        printf("    1/2<s+e;s+e>         %20.12f\n", exact.e_approx);
        printf("    <n-s;e>              %20.12f\n", exact.e_n_corr);
        printf("    sqrt(a/2pi) sum Z^2  %20.12f\n", exact.nrep_corr);
        printf("    E_ov                 %20.12e\n", exact.E_ov);
        printf("    E_smooth             %20.12f\n", exact.E_smooth);
        printf("    identity residual    %20.3e   (E_conv - E_smooth, must be ~0)\n",
               exact.E_conv - exact.E_smooth);
        printf("\n");
        printf("  %8s %4s %14s %14s %14s %14s %14s %14s\n",
               "thresh", "k", "d(E_ne)", "d(E_ee)", "d(E_conv)",
               "d(quad)", "d(corr)", "d(E_smooth)");
        fflush(stdout);
    }

    // ---- numerical evaluation at each threshold ----------------------------
    for (double thresh : threshes) {
        FunctionDefaults<3>::set_thresh(thresh);
        if (kwave > 0) {
            FunctionDefaults<3>::set_k(kwave);
        } else {
            // same ladder moldft uses
            int kk = 8;
            if (thresh >= 1e-4) kk = 6;
            else if (thresh >= 1e-6) kk = 8;
            else if (thresh >= 1e-8) kk = 10;
            else kk = 12;
            FunctionDefaults<3>::set_k(kk);
        }
        FunctionDefaults<3>::set_refine(true);
        FunctionDefaults<3>::set_initial_level(2);
        FunctionDefaults<3>::set_autorefine(false);
        FunctionDefaults<3>::set_apply_randomize(false);
        FunctionDefaults<3>::set_project_randomize(false);
        FunctionDefaults<3>::set_cubic_cell(-L, L);
        // MUST be cleared when k or thresh changes, or the next threshold silently reuses
        // the previous one's convolution coefficients and every row of the table comes out
        // identical.
        GaussianConvolution1DCache<double>::map.clear();

        const double op_thresh = (opeps > 0.0) ? opeps : 0.001 * thresh;
        SeparatedConvolution<double, 3> coulop = CoulombOperator(world, lo, op_thresh);

        // nuclear potential and the smoothed nuclear quantities
        PotentialManager pm(molecule, "none");
        pm.make_nuclear_potential(world);
        pm.make_smoothed_nuclear_charge_density(world, avalue);
        pm.make_smoothed_nuclear_potential(world, avalue);
        const real_function_3d& vnuc = pm.vnuclear();
        const real_function_3d& sdens = pm.smooth_density();
        const real_function_3d& unuc = pm.unuclear();

        // the model electron density, POSITIVE
        real_function_3d rho =
            real_factory_3d(world)
                .functor(real_functor_3d(new GaussianDensityFunctor(
                    gdens, molecule.get_all_coords_vec())))
                .thresh(thresh)
                .truncate_on_project();
        rho.truncate();

        // conventional route, exactly as moldft computes it
        const double E_nn_num = molecule.nuclear_repulsion_energy();
        const double E_ne_num = inner(rho, vnuc);   // vnuc = -u(r;n), rho > 0, so this is negative
        const double E_ee_num = 0.5 * inner(rho, apply(coulop, rho));
        const double E_conv_num = E_nn_num + E_ne_num + E_ee_num;

        // smoothed route, exactly as SCF::solve computes it
        real_function_3d spe = sdens - rho;
        real_function_3d Gspe = apply(coulop, spe);
        const double e_approx_num = 0.5 * inner(spe, Gspe);
        spe.clear(false);
        Gspe.clear(false);
        const double e_n_corr_num = inner(vnuc + unuc, rho);
        const double nrep_corr_num = molecule.nuclear_repulsion_correction(avalue);
        const double e_ov_num = molecule.smoothed_density_overlap_energy(avalue);
        const double E_smooth_num = e_approx_num + e_n_corr_num - nrep_corr_num + e_ov_num;

        // consistent route: the nuclear potential and the nuclear repulsion from the same
        // operator as E_ee. G is applied to s and rho separately, the smoothed route applies
        // it to s - rho; by linearity the two agree up to screening and truncation.
        real_function_3d Gs = apply(coulop, sdens);
        const double E_nn_cons = 0.5 * inner(sdens, Gs) - nrep_corr_num + e_ov_num;
        const double E_ne_cons = inner(rho, vnuc + unuc - Gs);   // <rho, Vnuc~>, Vnuc~ = -G s + c
        Gs.clear(false);
        const double E_cons_num = E_nn_cons + E_ne_cons + E_ee_num;

        // density-independent diagnostics
        const double strace = sdens.trace();
        const double ss_proj = inner(sdens, unuc);
        const double rho_trace = rho.trace();

        if (world.rank() == 0) {
            printf("  %8.1e %4d %14.4e %14.4e %14.4e %14.4e %14.4e %14.4e\n",
                   thresh, FunctionDefaults<3>::get_k(),
                   E_ne_num - exact.E_ne,
                   E_ee_num - exact.E_ee,
                   E_conv_num - exact.E_conv,
                   e_approx_num - exact.e_approx,
                   e_n_corr_num - exact.e_n_corr,
                   E_smooth_num - exact.E_smooth);
            printf("           %s op_eps=%.1e  int(rho)=%.10f (exact %.10f)  "
                   "int(s)=%.10f  <s;s>_proj err=%.3e\n",
                   "         ", op_thresh, rho_trace, gdens.total_charge(),
                   strace, ss_proj - exact.P_ss);
            printf("           %s consistent: d(E_nn~)=%.4e  d(E_ne~)=%.4e  d(E_cons)=%.4e\n",
                   "         ", E_nn_cons - exact.E_nn, E_ne_cons - exact.E_ne,
                   E_cons_num - exact.E_conv);
            fflush(stdout);
        }
        world.gop.fence();
    }

    if (world.rank() == 0) printf("\n");
    finalize();
    return 0;
}
