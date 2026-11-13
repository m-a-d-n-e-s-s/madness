// Molecule::nuclear_charge_density is the sum of the densities of all atoms
// in reach of the point, whatever the order of the atoms.

#include <madness/chem/molecule.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <vector>

using namespace madness;

namespace {

/// compares the density of `mol` with the sum of the densities of `parts`
/// at points along z, at distances `fracs * reach` from `z0`
bool is_sum_of(const Molecule& mol, const std::vector<Molecule>& parts,
               double z0, double reach, const std::vector<double>& fracs) {
  bool success = true;
  for (double frac : fracs) {
    const double z = z0 + frac * reach;
    const double total = mol.nuclear_charge_density(0.0, 0.0, z);
    double expected = 0.0;
    for (const auto& part : parts)
      expected += part.nuclear_charge_density(0.0, 0.0, z);
    const bool ok = expected != 0.0 &&
                    std::abs(total - expected) <= 1e-12 * std::abs(expected);
    print("  d/reach", frac, "density", total, "expected", expected,
          ok ? "ok" : "FAILED");
    success = success && ok;
  }
  return success;
}

}  // namespace

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  std::cout.precision(12);
  const double eprec = 1e-4;
  bool success = true;

  const Atom H(0.0, 0.0, 1.0, 1.0, 1);
  const Molecule mol_H({H}, eprec);
  const double reach_H = mol_H.nuclear_density_reach(0);
  const std::vector<double> fracs = {0.0, 0.1, 0.25, 0.5, 0.9};

  // an uncharged (ghost) atom listed first, 1 bohr from H: within 6 of its
  // default smoothing radii of 1 bohr, the old cutoff of its density, which
  // made it hide H. It has no density, hence no reach, and adds nothing
  {
    const Molecule mol({Atom(0.0, 0.0, 0.0, 0.0, 0), H}, eprec);
    const double rc_ghost = 1.0 / mol.get_rcut()[0];
    MADNESS_CHECK(Molecule::nuclear_density_reach_in_radii * rc_ghost >
                  1.0 + reach_H);
    const bool ok = mol.nuclear_density_reach(0) == 0.0 &&
                    mol.max_nuclear_density_reach() == reach_H;
    print("ghost atom listed before H: reach of ghost",
          mol.nuclear_density_reach(0), "max reach",
          mol.max_nuclear_density_reach(), ok ? "ok" : "FAILED");
    success = ok && success;
    success = is_sum_of(mol, {mol_H}, 1.0, reach_H, fracs) && success;
  }

  // two charged centers in reach of each other: the densities add
  {
    const Atom H2(0.0, 0.0, 1.0 + 0.5 * reach_H, 1.0, 1);
    const Molecule mol({H, H2}, eprec);
    const Molecule mol_H2({H2}, eprec);
    print("two H half a reach apart");
    success = is_sum_of(mol, {mol_H, mol_H2}, 1.0, reach_H, fracs) && success;
  }

  print("test_nuclear_charge_density", success ? "passed" : "FAILED");
  finalize();
  return success ? 0 : 1;
}
