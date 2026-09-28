//
// test the DFT-D3 empirical dispersion correction (chem/dispersion.h)
//
// Reference numbers were produced with the s-dftd3 command line program that
// ships with the library, on Turbomole `coord` files (which are in bohr, so the
// geometries below are the same numbers verbatim):
//
//     s-dftd3 --bj   <functional> [--atm] --grad g.txt --json out.json coord
//     s-dftd3 --zero <functional> [--atm] --grad g.txt --json out.json coord
//
// Regenerating them is a numerical claim: they pin this interface's argument
// order, units and array layout against an independent implementation of the
// same model, so a mismatch means MADNESS is calling the library wrong, not
// that the reference is stale.
//

#include <madness/mra/mra.h>
#include <madness/chem/dispersion.h>
#include <madness/chem/molecule.h>
#include <madness/world/test_utilities.h>

#include <algorithm>
#include <cmath>
#include <string>
#include <vector>

using namespace madness;

namespace {

/// water at the geometry of the scf_h2o_hf example deck (bohr)
Molecule h2o() {
    Molecule mol;
    mol.add_atom(0.0, 0.0, 0.0, 8.0, 8);
    mol.add_atom(0.0, 1.421178283, -1.064386617, 1.0, 1);
    mol.add_atom(0.0, -1.421178283, -1.064386617, 1.0, 1);
    return mol;
}

/// two neon atoms 6 bohr apart -- pure intermolecular dispersion
Molecule ne2() {
    Molecule mol;
    mol.add_atom(0.0, 0.0, 0.0, 10.0, 10);
    mol.add_atom(0.0, 0.0, 6.0, 10.0, 10);
    return mol;
}

/// LiH with a 3 bohr bond -- the geometry of the qc example deck
Molecule lih() {
    Molecule mol;
    mol.add_atom(0.0, 0.0, 0.0, 3.0, 3);
    mol.add_atom(0.0, 0.0, 3.0, 1.0, 1);
    return mol;
}

} // anonymous namespace


/// energies against s-dftd3's own command line program
int test_energies(World& world) {
    test_output t("testing D3 energies against the s-dftd3 reference values");

    // {spec, functional, atm}, molecule, reference energy in Hartree
    struct ref { const char* spec; const char* func; bool atm; Molecule mol; double e; };
    const std::vector<ref> refs = {
        {"d3bj",   "pbe",   false, h2o(),  -3.5943972560030610e-04},
        {"d3zero", "pbe",   false, h2o(),  -3.4068551849814270e-06},
        {"d3bj",   "pbe",   true,  h2o(),  -3.5943954605137060e-04},
        {"d3bj",   "b3lyp", false, h2o(),  -5.7382229713423441e-04},
        {"d3bj",   "hf",    false, h2o(),  -4.5133155058089784e-03},
        {"d3bj",   "pbe",   false, ne2(),  -8.4984069633701052e-05},
        {"d3zero", "pbe",   false, ne2(),  -1.3791110201610425e-04},
        {"d3bj",   "pbe",   false, lih(),  -2.2896520201444101e-04},
        {"d3zero", "pbe",   false, lih(),  -5.6181999915464882e-07},
    };

    for (const auto& r : refs) {
        DispersionCorrection d(r.spec, "", r.func, r.atm);
        const double e = d.energy(world, r.mol);
        t.checkpoint(e, r.e, 1.e-12,
                     d.description() + (r.atm ? " (atm)" : "") + ", "
                     + std::to_string(r.mol.natom()) + " atoms");
    }
    return t.end();
}


/// the analytic gradient, against s-dftd3 and against finite differences
int test_gradient(World& world) {
    test_output t("testing the D3 gradient");

    // s-dftd3 --bj pbe --grad, water, [3*atom + axis]
    const std::vector<double> gref_h2o = {
        0.0000000000000000e+00,  0.0000000000000000e+00, -8.9435960556183396e-07,
        0.0000000000000000e+00,  1.2222306768601004e-06,  4.4717980278091698e-07,
        0.0000000000000000e+00, -1.2222306768601004e-06,  4.4717980278091698e-07};

    {
        DispersionCorrection d("d3bj", "pbe", "", false);
        const Tensor<double> g = d.gradient(world, h2o());
        double err = 0.0;
        for (long i = 0; i < g.size(); ++i)
            err = std::max(err, std::fabs(g[i] - gref_h2o[i]));
        t.checkpoint(err, 1.e-14, "H2O D3(BJ)/pbe gradient vs. s-dftd3");
    }
    {
        // the atoms are 6 bohr apart along z and attract, so the correction pulls
        // atom 0 towards +z and atom 1 towards -z
        DispersionCorrection d("d3bj", "pbe", "", false);
        const Tensor<double> g = d.gradient(world, ne2());
        t.checkpoint(g[2], -4.3742513392993383e-05, 1.e-14, "Ne2 dE/dz(0) vs. s-dftd3");
        t.checkpoint(g[5], 4.3742513392993383e-05, 1.e-14, "Ne2 dE/dz(1) vs. s-dftd3");
    }

    // finite differences: independent of the reference numbers above, this
    // catches an analytic gradient that is transposed, misordered or scaled
    for (const std::string spec : {"d3bj", "d3zero"}) {
        DispersionCorrection d(spec, "pbe", "", false);
        Molecule mol = h2o();
        const Tensor<double> g = d.gradient(world, mol);
        const Tensor<double> x0 = mol.get_all_coords();
        const double h = 1.e-4;
        double err = 0.0;
        for (long iatom = 0; iatom < long(mol.natom()); ++iatom) {
            for (int axis = 0; axis < 3; ++axis) {
                Tensor<double> x = copy(x0);
                x(iatom, axis) += h;
                mol.set_all_coords(x);
                const double ep = d.energy(world, mol);
                x(iatom, axis) -= 2.0 * h;
                mol.set_all_coords(x);
                const double em = d.energy(world, mol);
                mol.set_all_coords(x0);
                err = std::max(err, std::fabs((ep - em) / (2.0 * h) - g[3 * iatom + axis]));
            }
        }
        t.checkpoint(err, 1.e-11, spec + ": analytic vs. central-difference gradient");
    }

    return t.end();
}


/// the Hessian: symmetric, translationally invariant, and consistent with the
/// gradient it is differenced from
int test_hessian(World& world) {
    test_output t("testing the D3 Hessian");

    DispersionCorrection d("d3bj", "pbe", "", false);
    const Molecule mol = h2o();
    const long n = 9;
    const Tensor<double> h = d.hessian(world, mol);

    t.checkpoint(h.ndim() == 2 and h.dim(0) == n and h.dim(1) == n,
                 "shape is 3*natom x 3*natom");

    // exact symmetry: hessian() averages the two off-diagonal estimates
    double asym = 0.0;
    for (long i = 0; i < n; ++i)
        for (long j = 0; j < n; ++j)
            asym = std::max(asym, std::fabs(h(i, j) - h(j, i)));
    t.checkpoint(asym, 1.e-16, "exactly symmetric");

    // translational invariance: displacing every atom along one axis cannot
    // change the energy, so each row must sum to zero over that axis
    double tsum = 0.0;
    for (long i = 0; i < n; ++i)
        for (int axis = 0; axis < 3; ++axis) {
            double s = 0.0;
            for (long jatom = 0; jatom < 3; ++jatom) s += h(i, 3 * jatom + axis);
            tsum = std::max(tsum, std::fabs(s));
        }
    t.checkpoint(tsum, 1.e-9, "rows obey the translational sum rule");

    // independent check of the whole matrix: second differences of the ENERGY,
    // which never touches the analytic gradient hessian() differences
    const Tensor<double> x0 = mol.get_all_coords();
    const double hh = 1.e-3;
    Molecule m = mol;
    double err = 0.0;
    for (long i = 0; i < n; ++i) {
        for (long j = 0; j <= i; ++j) {
            auto energy_at = [&](double di, double dj) {
                Tensor<double> x = copy(x0);
                x(i / 3, i % 3) += di;
                x(j / 3, j % 3) += dj;
                m.set_all_coords(x);
                return d.energy(world, m);
            };
            const double fd = (energy_at(hh, hh) - energy_at(hh, -hh)
                               - energy_at(-hh, hh) + energy_at(-hh, -hh))
                              / (4.0 * hh * hh);
            err = std::max(err, std::fabs(fd - h(i, j)));
        }
    }
    m.set_all_coords(x0);
    t.checkpoint(err, 1.e-6, "agrees with second differences of the energy");

    // against an outside implementation: the same second derivative obtained by
    // central-differencing the s-dftd3 CLI's own analytic gradient, LiH at
    // r = 3 bohr with D3(BJ)/hf --
    //   s-dftd3 --bj hf --grad g.txt --json o.json coord   at z = 3 +/- 1e-4
    //   d2E/dz(H)^2 = (g_plus[5] - g_minus[5]) / 2e-4
    {
        DispersionCorrection dl("d3bj", "hf", "", false);
        const Tensor<double> hl = dl.hessian(world, lih());
        t.checkpoint(hl(5, 5), -6.204984e-05, 1.e-10, "LiH zz(H,H) vs. the s-dftd3 CLI");
        t.checkpoint(hl(5, 2), 6.204984e-05, 1.e-10, "LiH zz(Li,H) vs. the s-dftd3 CLI");
    }

    // inactive -> exactly zero, same shape
    DispersionCorrection off("none", "pbe", "", false);
    const Tensor<double> hz = off.hessian(world, mol);
    t.checkpoint(hz.dim(0) == n and hz.absmax() == 0.0,
                 "inactive contributes a zero Hessian of the right shape");

    return t.end();
}


/// an inactive correction contributes exactly nothing
int test_inactive(World& world) {
    test_output t("testing the inactive dispersion correction");

    for (const std::string spec : {"", "none"}) {
        DispersionCorrection d(spec, "pbe", "", false);
        t.checkpoint(not d.active(), "'" + spec + "' is inactive");
        t.checkpoint(d.description() == "none", "'" + spec + "' has no description");
        t.checkpoint(d.energy(world, h2o()) == 0.0, "'" + spec + "' contributes no energy");
        const Tensor<double> g = d.gradient(world, h2o());
        t.checkpoint(g.size() == 9, "'" + spec + "' returns a 3*natom gradient");
        t.checkpoint(g.absmax() == 0.0, "'" + spec + "' contributes no gradient");
    }

    // the default constructed object is the one SCF holds when `dispersion` is
    // not set at all
    DispersionCorrection d;
    t.checkpoint(not d.active(), "the default constructed correction is inactive");

    // engines with no dispersion term in their energy expression call reject();
    // an inactive correction must pass straight through
    bool threw = false;
    try { d.reject("znemo"); } catch (...) { threw = true; }
    t.checkpoint(not threw, "reject() passes an inactive correction through");

    return t.end();
}


/// bad input is rejected loudly rather than silently returning zero
int test_bad_input(World& world) {
    test_output t("testing rejection of bad dispersion input");

    auto throws = [](auto&& f) {
        try { f(); } catch (...) { return true; }
        return false;
    };

    t.checkpoint(throws([] { (void)DispersionCorrection("d4", "pbe", "", false); }),
                 "an unknown correction is rejected");
    t.checkpoint(throws([] { (void)DispersionCorrection("d3bj", "", "", false); }),
                 "a missing functional name is rejected");
    t.checkpoint(throws([] { (void)DispersionCorrection("d3bj", "not_a_functional", "", false); }),
                 "an unparameterized functional is rejected");
    // a raw libxc spec is exactly the case the explicit override exists for
    t.checkpoint(throws([] { (void)DispersionCorrection("d3bj", "GGA_X_PBE 1.0 GGA_C_PBE 1.0", "", false); }),
                 "a raw libxc xc line is rejected");
    t.checkpoint(not throws([] { (void)DispersionCorrection("d3bj", "GGA_X_PBE 1.0", "pbe", false); }),
                 "... and dispersion_functional overrides it");
    // "none" is the parameter default for dispersion_functional and must fall
    // back to the xc line rather than be looked up as a functional name
    t.checkpoint(not throws([] { (void)DispersionCorrection("d3bj", "pbe", "none", false); }),
                 "dispersion_functional=none falls back to xc");

    // an active correction on an engine that cannot apply it must abort
    DispersionCorrection active("d3bj", "pbe", "", false);
    t.checkpoint(throws([&] { active.reject("znemo"); }),
                 "reject() aborts an active correction");

    return t.end();
}


int main(int argc, char** argv) {
    World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);

    int error = 0;
    if (not DispersionCorrection::available()) {
        print("MADNESS was built without simple-dftd3 -- nothing to test");
    } else {
        print("linked against simple-dftd3", DispersionCorrection::library_version());
        error += test_energies(world);
        error += test_gradient(world);
        error += test_hessian(world);
        error += test_inactive(world);
        error += test_bad_input(world);
    }

    world.gop.fence();
    finalize();
    return error;
}
