/// \file test_population.cc
/// \brief the population schemes of population.h on matrices alone
///
/// The schemes need only S = <chi|chi> and X = <chi|phi>, so they are checked in
/// R^n: a non-orthogonal "minimal basis" chi and orthonormal occupied orbitals phi
/// that lie partly outside its span, as MRA orbitals do. The IAO route is compared
/// with a direct construction of Knizia's definition in R^n; all schemes are checked
/// for the sum rule, for invariance under rotations of the occupied orbitals and,
/// for orbitals inside the basis, against the textbook Mulliken and Loewdin formulas.

#include <madness/chem/population.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <iostream>
#include <random>
#include <string>
#include <vector>

using namespace madness;
using namespace madness::population;

namespace {

int failures = 0;
int checks = 0;

void check(const bool ok, const std::string& what, const double value = 0.0) {
    ++checks;
    if (not ok) {
        ++failures;
        print("FAILED:", what, value);
    }
}

std::mt19937 rng(7);

Tensor<double> random_matrix(const long n, const long m) {
    std::normal_distribution<double> g(0.0, 1.0);
    Tensor<double> a(n, m);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < m; ++j) a(i, j) = g(rng);
    return a;
}

/// columns made orthonormal (symmetric orthonormalization)
Tensor<double> orthonormal_columns(const Tensor<double>& a) {
    return inner(a, matrix_power(inner(transpose(a), a), -0.5));
}

double max_abs_diff(const Tensor<double>& a, const Tensor<double>& b) {
    return (a - b).absmax();
}

double sum(const Tensor<double>& a) {
    double s = 0.0;
    for (long i = 0; i < a.dim(0); ++i) s += a(i);
    return s;
}

const std::vector<int> atom_of = {0, 0, 0, 0, 0, 1, 1, 2, 2};   // minimal basis function -> atom
const long natom = 3, n = 40, m = 9, o = 4;                       // R^n, minimal basis, occupied orbitals

/// Knizia's IAO populations built directly in R^n
Tensor<double> iao_direct(const Tensor<double>& chi, const Tensor<double>& phi, const Tensor<double>& occ) {
    const Tensor<double> S = inner(transpose(chi), chi), X = inner(transpose(chi), phi);
    Tensor<double> pt = inner(chi, inner(matrix_power(S, -1.0), X));      // projections of phi onto span(chi)
    pt = orthonormal_columns(pt);                                         // the depolarized orbitals
    const Tensor<double> O = inner(phi, transpose(phi)), Ot = inner(pt, transpose(pt));
    Tensor<double> one(n, n);
    for (long i = 0; i < n; ++i) one(i, i) = 1.0;
    const Tensor<double> A = inner(inner(O, Ot) + inner(one - O, one - Ot), chi);
    const Tensor<double> Q = inner(transpose(orthonormal_columns(A)), phi);
    Tensor<double> e(natom);
    for (long rho = 0; rho < m; ++rho)
        for (long i = 0; i < o; ++i) e(atom_of[rho]) += occ(i) * Q(rho, i) * Q(rho, i);
    return e;
}

void test_iao_and_sum_rules() {
    const Tensor<double> chi = random_matrix(n, m);
    // occupied orbitals: mostly in span(chi), plus a polarization part outside it
    const Tensor<double> phi = orthonormal_columns(inner(chi, random_matrix(m, o)) + 0.3 * random_matrix(n, o));
    Tensor<double> occ(o);
    occ.fill(1.0);
    const Tensor<double> S = inner(transpose(chi), chi), X = inner(transpose(chi), phi);

    const SpinPopulation iao = iao_population(S, X, occ, atom_of, natom);
    check(std::abs(sum(iao.electrons) - o) < 1.e-12, "iao: sum rule", sum(iao.electrons) - o);
    check(iao.max_norm_error < 1.e-12, "iao: unit columns of Q", iao.max_norm_error);
    check(iao.min_w_eigenvalue > 0.0, "iao: W positive", iao.min_w_eigenvalue);
    const double dd = max_abs_diff(iao.electrons, iao_direct(chi, phi, occ));
    check(dd < 1.e-10, "iao: matrix route equals the direct construction", dd);

    // the spilling: the share of phi outside span(chi), from the projector in R^n
    const Tensor<double> proj = inner(chi, inner(matrix_power(S, -1.0), X));  // P_chi phi
    double inside = 0.0;
    for (long i = 0; i < o; ++i)
        for (long k = 0; k < n; ++k) inside += proj(k, i) * phi(k, i);
    for (const std::string scheme : {"mulliken", "lowdin"}) {
        const SpinPopulation p = projection_population(scheme, S, X, occ, atom_of, natom);
        check(std::abs(sum(p.electrons) - o) < 1.e-12, scheme + ": sum rule", sum(p.electrons) - o);
        check(std::abs(p.spilling - (1.0 - inside / o)) < 1.e-12, scheme + ": spilling",
              p.spilling - (1.0 - inside / o));
        check(p.spilling > 0.01, scheme + ": the polarization part spills", p.spilling);
    }

    // invariance under a rotation of the occupied orbitals
    const Tensor<double> U = orthonormal_columns(random_matrix(o, o));
    const Tensor<double> X2 = inner(X, U);
    for (const std::string& scheme : schemes()) {
        const double d = max_abs_diff(spin_population(scheme, S, X, occ, atom_of, natom).electrons,
                                      spin_population(scheme, S, X2, occ, atom_of, natom).electrons);
        check(d < 1.e-12, scheme + ": invariant under occupied rotations", d);
    }

    // doubly occupied (a restricted closed shell): twice the electrons
    Tensor<double> occ2(o);
    occ2.fill(2.0);
    for (const std::string& scheme : schemes()) {
        const double d = max_abs_diff(spin_population(scheme, S, X, occ2, atom_of, natom).electrons,
                                      2.0 * spin_population(scheme, S, X, occ, atom_of, natom).electrons);
        check(d < 1.e-12, scheme + ": occupations scale the populations", d);
    }
}

void test_textbook_formulas() {
    // orbitals inside span(chi): no spilling, and the textbook formulas apply
    const Tensor<double> chi = random_matrix(n, m);
    const Tensor<double> phi = orthonormal_columns(inner(chi, random_matrix(m, o)));
    Tensor<double> occ(o);
    occ.fill(1.0);
    const Tensor<double> S = inner(transpose(chi), chi), X = inner(transpose(chi), phi);
    const Tensor<double> C = inner(matrix_power(S, -1.0), X);
    const Tensor<double> Pd = inner(C, transpose(C));                       // the density matrix
    const Tensor<double> PS = inner(Pd, S);
    const Tensor<double> Sh = matrix_power(S, 0.5);
    const Tensor<double> SPS = inner(inner(Sh, Pd), Sh);
    Tensor<double> mull(natom), lowd(natom);
    for (long mu = 0; mu < m; ++mu) {
        mull(atom_of[mu]) += PS(mu, mu);
        lowd(atom_of[mu]) += SPS(mu, mu);
    }
    const SpinPopulation pm = projection_population("mulliken", S, X, occ, atom_of, natom);
    const SpinPopulation pl = projection_population("lowdin", S, X, occ, atom_of, natom);
    check(std::abs(pm.spilling) < 1.e-12, "in-span orbitals: no spilling", pm.spilling);
    check(max_abs_diff(pm.electrons, mull) < 1.e-12, "mulliken: (PS)_mumu", max_abs_diff(pm.electrons, mull));
    check(max_abs_diff(pl.electrons, lowd) < 1.e-12, "lowdin: (S^1/2 P S^1/2)_mumu", max_abs_diff(pl.electrons, lowd));
}

void test_singular_w() {
    // an occupied orbital orthogonal to the minimal basis: W is singular, iao must refuse
    const Tensor<double> chi = random_matrix(n, m);
    Tensor<double> phi = orthonormal_columns(inner(chi, random_matrix(m, o)));
    Tensor<double> v = random_matrix(n, 1);
    const Tensor<double> S = inner(transpose(chi), chi);
    v -= inner(chi, inner(matrix_power(S, -1.0), inner(transpose(chi), v)));   // remove its part in span(chi)
    for (long i = 0; i < o - 1; ++i) {                                          // and its overlap with phi
        double c = 0.0;
        for (long k = 0; k < n; ++k) c += phi(k, i) * v(k, 0);
        for (long k = 0; k < n; ++k) v(k, 0) -= c * phi(k, i);
    }
    v.scale(1.0 / v.normf());
    phi(_, o - 1) = v(_, 0);
    Tensor<double> occ(o);
    occ.fill(1.0);
    bool threw = false;
    try {
        iao_population(S, inner(transpose(chi), phi), occ, atom_of, natom);
    } catch (const MadnessException&) {
        threw = true;
    }
    check(threw, "iao: refuses a minimal basis that misses an orbital");
}

} // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    {
        test_iao_and_sum_rules();
        test_textbook_formulas();
        test_singular_w();
        if (world.rank() == 0) print("test_population:", checks, "checks,", failures, "failed");
    }
    finalize();
    return failures == 0 ? 0 : 1;
}
