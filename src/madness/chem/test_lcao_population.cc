/// \file test_lcao_population.cc
/// \brief population analysis of LCAO wavefunctions (analytic overlaps), and the cross-basis overlap
///
/// Water (RHF) and BeH (UHF doublet) in 6-31g**, converged with the LCAO SCF and its in-core integrals.
/// Checked: the overlap of a basis with itself against the SCF's own overlap; the sum rules of all
/// three schemes (charges and spins); invariance under rotations of the occupied orbitals; unit
/// columns of the IAO coefficients; and no spilling when the projection basis is the orbitals' own.

#include <madness/chem/lcao_integrals.h>
#include <madness/chem/lcao_scf.h>
#include <madness/chem/molecularbasis.h>
#include <madness/chem/population.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <random>
#include <string>
#include <vector>

using namespace madness;

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

Molecule water() {
    Molecule m;
    m.add_atom(0.0, 0.0, 0.2216649, 8.0, 8);
    m.add_atom(0.0, 1.4309006, -0.8866595, 1.0, 1);
    m.add_atom(0.0, -1.4309006, -0.8866595, 1.0, 1);
    return m;
}

Molecule beh() {
    Molecule m;
    m.add_atom(0.0, 0.0, 0.0, 4.0, 4);
    m.add_atom(0.0, 0.0, 2.54, 1.0, 1);
    return m;
}

Tensor<double> random_rotation(const long n) {
    std::mt19937 rng(11);
    std::normal_distribution<double> g(0.0, 1.0);
    Tensor<double> a(n, n);
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j) a(i, j) = g(rng);
    return inner(a, population::matrix_power(inner(transpose(a), a), -0.5));
}

std::vector<double> get(const nlohmann::json& j, const std::string& scheme, const std::string& what) {
    return j["schemes"][scheme][what].get<std::vector<double>>();
}

double maxdiff(const std::vector<double>& a, const std::vector<double>& b) {
    double d = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) d = std::max(d, std::abs(a[i] - b[i]));
    return d;
}

double sum(const std::vector<double>& a) {
    double s = 0.0;
    for (const double x : a) s += x;
    return s;
}

void test_overlap(World& world) {
    const Molecule mol = water();
    AtomicBasisSet big, min;
    big.read_file("6-31gss");
    min.read_file("sto-3g");
    const std::vector<lcao::Shell> sb = lcao::make_shells(mol, big), sm = lcao::make_shells(mol, min);
    const lcao::SeparatedGaussianIntegrals ints(sb, lcao::GaussianKernel::coulomb(1.e-4, 50.0, 1.e-3));
    const double d = (lcao::overlap(sb, sb) - ints.overlap(world)).absmax();
    check(d < 1.e-14, "overlap of a basis with itself equals the SCF's overlap", d);
    const double t = (lcao::overlap(sm, sb) - transpose(lcao::overlap(sb, sm))).absmax();
    check(t < 1.e-15, "cross overlap: S(a,b) = S(b,a)^T", t);
}

void test_molecule(World& world, const std::string& name, const Molecule& mol, const long na, const long nb) {
    AtomicBasisSet basis;
    basis.read_file("6-31gss");
    LCAOParameters lp;
    lp.set_user_defined_value("print_level", 0);
    lp.set_user_defined_value("econv", 1.e-10);
    lp.set_user_defined_value("dconv", 1.e-8);
    lcao::LCAOSCF scf(world, mol, basis, na, nb, lp, false);
    scf.solve();
    check(scf.converged(), name + ": LCAO SCF converged");
    const Tensor<double> Ca = copy(scf.coefficients(0)(_, Slice(0, na - 1)));
    const Tensor<double> Cb = copy(scf.coefficients(1)(_, Slice(0, nb - 1)));
    const std::vector<std::string> schemes = {"mulliken", "lowdin", "iao"};
    const nlohmann::json j = lcao::population_analysis(mol, basis, Ca, Cb, schemes, "6-31g", "sto-3g", 0.0, 0);
    const nlohmann::json jr = lcao::population_analysis(mol, basis, inner(Ca, random_rotation(na)),
                                                       inner(Cb, random_rotation(nb)), schemes, "6-31g", "sto-3g",
                                                       0.0, 0);
    for (const auto& s : schemes) {
        const std::vector<double> q = get(j, s, "charges");
        check(std::abs(sum(q)) < 1.e-10, name + " " + s + ": charges sum to 0", sum(q));
        const double dq = maxdiff(q, get(jr, s, "charges"));
        check(dq < 1.e-10, name + " " + s + ": invariant under occupied rotations", dq);
        if (na != nb) {
            const std::vector<double> sp = get(j, s, "spin");
            check(std::abs(sum(sp) - double(na - nb)) < 1.e-10, name + " " + s + ": spins sum to na-nb",
                  sum(sp) - double(na - nb));
            const double ds = maxdiff(sp, get(jr, s, "spin"));
            check(ds < 1.e-10, name + " " + s + ": spins invariant under occupied rotations", ds);
        }
    }
    const double qn = j["schemes"]["iao"]["max_norm_error"].get<double>();
    check(qn < 1.e-10, name + " iao: unit columns of Q", qn);
    // the orbitals' own basis holds them completely
    const nlohmann::json js = lcao::population_analysis(mol, basis, Ca, Cb, {"mulliken"}, "6-31gss", "sto-3g", 0.0, 0);
    const double sp = std::abs(js["schemes"]["mulliken"]["spilling_alpha"].get<double>());
    check(sp < 1.e-12, name + ": no spilling in the orbitals' own basis", sp);
}

} // namespace

int main(int argc, char** argv) {
    World& world = initialize(argc, argv);
    {
        startup(world, argc, argv, true);
        test_overlap(world);
        test_molecule(world, "water", water(), 5, 5);
        test_molecule(world, "BeH", beh(), 3, 2);
        if (world.rank() == 0) print("test_lcao_population:", checks, "checks,", failures, "failed");
    }
    finalize();
    return failures == 0 ? 0 : 1;
}
