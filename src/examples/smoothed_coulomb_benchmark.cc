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
      smoothed_coulomb_benchmark input --thresh=1e-6 --construction   # timing harness, below
      smoothed_coulomb_benchmark input --thresh=1e-6 --construction=cost   # the same, no model density

  A third, "consistent" route is printed below each row: the conventional assembly with the
  nuclear potential and the nuclear repulsion taken from the same operator as E_ee,
  Vnuc~ = -G s + (Vnuc + u(s)) and E_nn~ = 1/2 <s;G s> - sqrt(a/2pi) sum Z^2 + E_ov.
  Algebraically it equals the smoothed route; numerically it applies G to s and rho
  separately instead of to s - rho. A last line repeats it with c = Vnuc + u(s) projected
  directly (NuclearPotentialRemainderFunctor, O(natom)), as moldft does with
  nuclear_potential=convolution, and compares <rho;c> with the two O(natom^2) projections.

  --construction: a harness for the cost of the nuclear potential at a geometry (phase 3). It
  builds the potential as moldft does, once per route: projected (master's
  PotentialManager::make_nuclear_potential, O(natom^2)) and by convolution (Make sdens, Project
  cnuclear, Convolve sdens, as SCF::make_nuclear_potential_by_convolution). It times each piece,
  reports the trees, and checks both against the exact reference: d(E_conv) and d(E_cons), and
  <rho; Vnuc~ - Vnuc>, and d(E_nn~) = E_nn~ - E_nn. It skips the smoothed route and the
  O(natom^2) u(s) projection. The model density is O(natom^2) too, since every point sums all
  its primitives, and from a few hundred atoms on it costs more than everything else:
  --construction=cost leaves it out, and with it every check but d(E_nn~). --profile then builds
  the route once more, one piece at a time: each projection and its reconstruct, and the
  convolution's apply, inner product, subtraction, truncate and reconstruct. With --s_safety=X
  that second s is projected at X thresh instead of 0.1 thresh, and its d(E_nn~) (and in the
  full mode its d(E_cons)) is printed.

  Sign convention: all densities below are POSITIVE and the pairing is
      P(f,g) = int int f(r) g(r') / |r-r'|,
  so with n, s positive nuclear densities and rho the positive electron density,
      E_ne = -P(n,rho),  E_ee = 1/2 P(rho,rho),  <s+e;s+e> = P(s-rho, s-rho).
*/

#include <madness/mra/mra.h>
#include <madness/mra/lbdeux.h>
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

/// FunctionDefaults for one threshold, as moldft's protocol sets them
static void set_function_defaults(double thresh, int kwave, double L) {
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
}

/// the node cost SCF::loadbal balances with (lbcost in SCF.h, which this example does not include)
struct scf_lbcost {
    double leaf_value;
    double parent_value;
    double operator()(const Key<3>& key, const FunctionNode<double, 3>& node) const {
        if (key.level() < 1) return leaf_value + parent_value;
        return node.is_leaf() ? leaf_value : parent_value;
    }
};

/// wall time of f(), fenced on both sides so that every rank's work is counted
template <typename funcT>
static double timed(World& world, funcT&& f) {
    world.gop.fence();
    const double t0 = wall_time();
    f();
    world.gop.fence();
    return wall_time() - t0;
}

/// --construction: the cost of the nuclear potential by projection and by convolution.
/// Without \p with_density (--construction=cost) it skips the model density, whose projection
/// sums all primitives at every point and dominates the harness from a few hundred atoms on.
static void construction_harness(World& world, const Molecule& molecule, const GaussianDensity& gdens,
                                 const ExactEnergies& exact, const std::vector<double>& threshes,
                                 int kwave, double lo, double opeps, double avalue, double L,
                                 bool with_density, bool profile, double s_safety) {
    if (world.rank() == 0) {
        printf("  construction harness: wall times in s, sizes in coefficients. The route is timed with a\n");
        printf("  warm operator cache, on a process map balanced as SCF::loadbal balances it.\n");
        if (with_density) {
            printf("  t(J) is the first apply of the operator (cold cache, model density).\n");
        } else {
            printf("  --construction=cost, no model density: the map is balanced on the projected potential\n");
            printf("  alone; t(cold) is a first apply of the operator to s (cold cache); only d(E_nn~) checks\n");
            printf("  accuracy.\n");
        }
        printf("  %8s %4s %10s %10s %10s %10s %10s %10s %7s %11s %11s %11s %11s %12s %12s %12s %12s\n",
               "thresh", "k", "t(proj)", with_density ? "t(J)" : "t(cold)", "t(sdens)", "t(cnuc)", "t(conv)",
               "t(conv rt)", "ratio", "size(proj)", "size(s)", "size(c)", "size(conv)", "d(E_conv)",
               "d(E_cons)", "<rho;dV>", "d(E_nn~)");
        fflush(stdout);
    }
    for (double thresh : threshes) {
        set_function_defaults(thresh, kwave, L);
        // as the SCF constructor sets it: the harness should build the trees moldft builds (the
        // benchmark's other mode keeps the default, so its numbers stay comparable with phase 2's)
        FunctionDefaults<3>::set_truncate_mode(1);
        const double op_thresh = (opeps > 0.0) ? opeps : 0.001 * thresh;
        SeparatedConvolution<double, 3> coulop = CoulombOperator(world, lo, op_thresh);
        PotentialManager pm(molecule, "none");

        // master's route; keep its handle, since the convolution replaces pm's vnuc
        const double t_proj = timed(world, [&] { pm.make_nuclear_potential(world); });
        const real_function_3d vnuc_proj = pm.vnuclear();

        // the model density and E_ee first: the first apply of an operator fills its cache of
        // convolution coefficients. An SCF pays that once per protocol level whichever route it
        // takes (master's in its first J), so the route is timed with a warm cache.
        real_function_3d rho;
        if (with_density) {
            rho = real_factory_3d(world)
                      .functor(real_functor_3d(new GaussianDensityFunctor(gdens, molecule.get_all_coords_vec())))
                      .thresh(thresh)
                      .truncate_on_project();
            rho.truncate();
        }

        // load balance as SCF::loadbal does (vnucextra 2, loadbalparts 2), with the model density
        // in place of the SCF's densities and orbitals: the SCF's later protocol levels build their
        // potential on such a map, and an apply depends on it
        if (world.size() > 1) {
            LoadBalanceDeux<3> lb(world);
            lb.add_tree(vnuc_proj, scf_lbcost{2.0, 16.0}, false);
            if (with_density) lb.add_tree(rho, scf_lbcost{1.0, 8.0}, false);
            world.gop.fence();
            FunctionDefaults<3>::redistribute(world, lb.load_balance(2.0));
            world.gop.fence();
        }
        double t_first = 0.0, E_ee = 0.0;
        if (with_density) {
            real_function_3d Grho;
            t_first = timed(world, [&] { Grho = apply(coulop, rho); });
            E_ee = 0.5 * inner(rho, Grho);
        }

        // the default's route, the three timers of SCF::make_nuclear_potential_by_convolution
        const double t_sdens = timed(world, [&] { pm.make_smoothed_nuclear_charge_density(world, avalue); });
        const double t_cnuc = timed(world, [&] { pm.make_nuclear_potential_remainder(world, avalue); });
        if (!with_density) {
            // without the model density, s fills the cache, and its potential is thrown away
            t_first = timed(world, [&] { real_function_3d Gs = apply(coulop, pm.smooth_density()); });
        }
        const double t_conv = timed(world, [&] { pm.make_nuclear_potential_by_convolution(coulop); });
        const real_function_3d& vnuc_conv = pm.vnuclear();
        const double t_route = t_sdens + t_cnuc + t_conv;

        // accuracy against the closed forms; d(E_nn~) needs no density
        const double E_nn_cons = 0.5 * pm.smoothed_density_coulomb()
                                 - molecule.nuclear_repulsion_correction(avalue)
                                 + molecule.smoothed_density_overlap_energy(avalue);
        const double d_nn = E_nn_cons - exact.E_nn;
        double d_conv = 0.0, d_cons = 0.0, rho_dV = 0.0;
        if (with_density) {
            const double E_ne_proj = inner(rho, vnuc_proj);
            const double E_ne_conv = inner(rho, vnuc_conv);
            d_conv = molecule.nuclear_repulsion_energy() + E_ne_proj + E_ee - exact.E_conv;
            d_cons = E_nn_cons + E_ne_conv + E_ee - exact.E_conv;
            rho_dV = E_ne_conv - E_ne_proj;
        }

        // in one tree state: inner products leave their arguments compressed, which carries about
        // 8/7 times the coefficients of the reconstructed tree
        const auto sz = [](const real_function_3d& f) {
            f.reconstruct();
            return static_cast<double>(f.size());
        };
        const double s_proj = sz(vnuc_proj), s_s = sz(pm.smooth_density());
        const double s_c = sz(pm.nuclear_potential_remainder()), s_conv = sz(vnuc_conv);
        if (world.rank() == 0) {
            printf("  %8.1e %4d %10.2f %10.2f %10.2f %10.2f %10.2f %10.2f %7.2f %11.4e %11.4e %11.4e %11.4e",
                   thresh, FunctionDefaults<3>::get_k(), t_proj, t_first, t_sdens, t_cnuc, t_conv, t_route,
                   t_route / t_proj, s_proj, s_s, s_c, s_conv);
            if (with_density) {
                printf(" %12.4e %12.4e %12.4e", d_conv, d_cons, rho_dV);
            } else {
                printf(" %12s %12s %12s", "-", "-", "-");
            }
            printf(" %12.4e\n", d_nn);
            fflush(stdout);
        }

        // --profile: the route once more, one piece at a time, in the order and the tree states
        // in which PotentialManager runs them (projections at 0.1 thresh, then reconstruct; apply;
        // the inner product, which compresses s and G~s; c - G~s; truncate; reconstruct).
        // --s_safety projects this s at s_safety * thresh instead of PotentialManager's 0.1 * thresh
        if (profile) {
            const double tol = 0.1 * thresh;
            real_function_3d s2, c2, Gs, v;
            const double tp_s = timed(world, [&] {
                s2 = real_factory_3d(world)
                         .functor(real_functor_3d(new SmoothedNuclearChargeDensityFunctor(molecule, avalue)))
                         .thresh(s_safety * thresh)
                         .truncate_on_project();
            });
            const double tr_s = timed(world, [&] { s2.set_thresh(thresh); s2.reconstruct(); });
            const double tp_c = timed(world, [&] {
                c2 = real_factory_3d(world)
                         .functor(real_functor_3d(new NuclearPotentialRemainderFunctor(molecule, avalue)))
                         .thresh(tol)
                         .truncate_on_project();
            });
            const double tr_c = timed(world, [&] { c2.set_thresh(thresh); c2.reconstruct(); });
            double ss = 0.0;
            const double t_apply = timed(world, [&] { Gs = apply(coulop, s2); });
            const double t_inner = timed(world, [&] { ss = inner(s2, Gs); });
            const double t_sub = timed(world, [&] { v = c2 - Gs; });
            const double t_trunc = timed(world, [&] { v.truncate(tol); });
            const double t_recon = timed(world, [&] { v.reconstruct(); });
            // the accuracy of this s: d(E_nn~), and with the model density d(E_cons)
            const double E_nn_p = 0.5 * ss - molecule.nuclear_repulsion_correction(avalue)
                                  + molecule.smoothed_density_overlap_energy(avalue);
            const double d_cons_p = with_density ? E_nn_p + inner(rho, v) + E_ee - exact.E_conv : 0.0;
            const double s_s2 = sz(s2);
            if (world.rank() == 0) {
                printf("  profile: s %.2f (project %.2f + reconstruct %.2f); c %.2f (%.2f + %.2f); "
                       "conv %.2f (apply %.2f + inner %.2f + c-Gs %.2f + truncate %.2f + reconstruct %.2f); "
                       "<s;G~s> - route's %.3e\n",
                       tp_s + tr_s, tp_s, tr_s, tp_c + tr_c, tp_c, tr_c,
                       t_apply + t_inner + t_sub + t_trunc + t_recon, t_apply, t_inner, t_sub, t_trunc, t_recon,
                       ss - pm.smoothed_density_coulomb());
                printf("  profile: s at %.2g thresh: size(s) %.4e, d(E_nn~) %.4e", s_safety, s_s2, E_nn_p - exact.E_nn);
                if (with_density) printf(", d(E_cons) %.4e", d_cons_p);
                printf("\n");
                fflush(stdout);
            }
        }
        world.gop.fence();
    }
}

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
        printf("    eprec                %20.3e\n", molecule.get_eprec());
        printf("\n");
    }
    if (parser.key_exists("construction")) {
        const std::string mode = parser.value("construction");
        MADNESS_CHECK_THROW(mode.empty() || mode == "cost", "--construction takes no value, or =cost");
        const double s_safety = parser.key_exists("s_safety") ? std::atof(parser.value("s_safety").c_str()) : 0.1;
        construction_harness(world, molecule, gdens, exact, threshes, kwave, lo, opeps, avalue, L, mode.empty(),
                             parser.key_exists("profile"), s_safety);
        if (world.rank() == 0) printf("\n");
        finalize();
        return 0;
    }
    if (world.rank() == 0) {
        printf("  %8s %4s %14s %14s %14s %14s %14s %14s\n",
               "thresh", "k", "d(E_ne)", "d(E_ee)", "d(E_conv)",
               "d(quad)", "d(corr)", "d(E_smooth)");
        fflush(stdout);
    }

    // ---- numerical evaluation at each threshold ----------------------------
    for (double thresh : threshes) {
        set_function_defaults(thresh, kwave, L);

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

        // the same with c projected directly, O(natom), as moldft builds Vnuc~ with
        // nuclear_potential=convolution
        pm.make_nuclear_potential_remainder(world, avalue);
        const real_function_3d& cnuc = pm.nuclear_potential_remainder();
        const double c_err = inner(rho, cnuc) - e_n_corr_num;    // <rho;c> - <rho;vnuc+unuc>
        const double E_ne_cons_c = inner(rho, cnuc - Gs);
        Gs.clear(false);
        const double E_cons_num = E_nn_cons + E_ne_cons + E_ee_num;
        const double E_cons_c_num = E_nn_cons + E_ne_cons_c + E_ee_num;

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
            printf("           %s remainder c: <rho;c>-<rho;vnuc+unuc>=%.4e  d(E_cons,c)=%.4e\n",
                   "         ", c_err, E_cons_c_num - exact.E_conv);
            fflush(stdout);
        }
        world.gop.fence();
    }

    if (world.rank() == 0) printf("\n");
    finalize();
    return 0;
}
