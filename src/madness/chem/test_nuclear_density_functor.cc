// NuclearDensityFunctor with periodic boundaries: the density of a nucleus
// on a face of the cell must be the same on both sides of the face.

#include <madness/chem/molecule.h>
#include <madness/chem/potentialmanager.h>
#include <madness/mra/mra.h>

#include <cmath>
#include <vector>

using namespace madness;

namespace {

/// compares the functor on the two sides of the upper z face, at distance
/// `frac * reach` from a nucleus on that face
bool same_on_both_sides(const NuclearDensityFunctor& rho, double L,
                        double reach, const std::vector<double>& fracs) {
  bool success = true;
  for (double frac : fracs) {
    const double d = frac * reach;
    const double inside = rho(coord_3d{0.0, 0.0, 0.5 * L - d});
    const double wrapped = rho(coord_3d{0.0, 0.0, -0.5 * L + d});
    const bool ok = inside != 0.0 &&
                    std::abs(wrapped - inside) <= 1e-12 * std::abs(inside);
    print("  d/reach", frac, "inside", inside, "wrapped", wrapped,
          ok ? "ok" : "FAILED");
    success = success && ok;
  }
  return success;
}

}  // namespace

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  startup(world, argc, argv);
  std::cout.precision(12);

  const double L = 20.0;
  FunctionDefaults<3>::set_cubic_cell(-0.5 * L, 0.5 * L);
  const BoundaryConditions<3> bc(BC_PERIODIC);
  const Tensor<double> cell = FunctionDefaults<3>::get_cell();
  const double eprec = 1e-4;

  bool success = true;

  // a light nucleus on a face next to a heavier one: the light nucleus has
  // the larger smoothing radius, hence the longer reach across the face
  {
    const Molecule mol({Atom(0.0, 0.0, 0.0, 3.0, 3),
                        Atom(0.0, 0.0, 0.5 * L, 1.0, 1)},
                       eprec);
    const double reach_Li = mol.nuclear_density_reach(0);
    const double reach_H = mol.nuclear_density_reach(1);
    print("mixed-Z cell: reach of Li", reach_Li, "reach of H", reach_H);
    MADNESS_CHECK(reach_Li < reach_H && 4.0 * reach_H < L);
    const NuclearDensityFunctor rho(mol, bc, cell);
    // 0.25 and 0.4 are beyond the reach of Li if it is less than 1/4 of H's
    success =
        same_on_both_sides(rho, L, reach_H, {0.1, 0.25, 0.4}) && success;
  }

  // a nucleus enlarged by rscale reaches rscale times farther
  {
    const Molecule mol({Atom(0.0, 0.0, 0.5 * L, 1.0, 1)}, eprec);
    const double rscale = 2.0;
    const double reach = mol.nuclear_density_reach(0, rscale);
    print("rscale", rscale, ": reach of H", reach);
    MADNESS_CHECK(4.0 * reach < L);
    const NuclearDensityFunctor rho(mol, bc, cell, 15, rscale);
    // 0.55 is beyond the reach of the unscaled nucleus
    success = same_on_both_sides(rho, L, reach, {0.1, 0.25, 0.55}) && success;
  }

  // the projected periodic density of the mixed-Z cell holds all of the
  // nuclear charge, including the part of H that wraps across the face
  {
    FunctionDefaults<3>::set_bc(bc);
    FunctionDefaults<3>::set_k(8);
    FunctionDefaults<3>::set_thresh(1e-6);
    const Molecule mol({Atom(0.0, 0.0, 0.0, 3.0, 3),
                        Atom(0.0, 0.0, 0.5 * L, 1.0, 1)},
                       eprec);
    const real_function_3d rho =
        real_factory_3d(world)
            .functor(real_functor_3d(new NuclearDensityFunctor(mol, bc, cell)))
            .truncate_on_project();
    const double charge = rho.trace();
    const bool ok = std::abs(charge - 4.0) < 1e-3;
    print("projected charge", charge, "expected", 4.0, ok ? "ok" : "FAILED");
    success = ok && success;
  }

  print("test_nuclear_density_functor", success ? "passed" : "FAILED");
  finalize();
  return success ? 0 : 1;
}
