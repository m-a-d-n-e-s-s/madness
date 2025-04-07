#pragma once
#include "GroundStateData.hpp"
#include "ResponseManager.hpp"
#include "ResponseState.hpp"
#include "ResponseVector.hpp"
#include <madness/mra/nonlinsol.h>

namespace ResponseSolverUtils {

inline std::vector<poperatorT> make_bsh_operators_response(
    World &world, const double shift, const double omega,
    const Tensor<double> &ground_energies, const double &lo) {

  double tol = FunctionDefaults<3>::get_thresh();
  // Sizes inferred from ground and omega
  size_t num_orbitals = ground_energies.size(); // number of orbitals
  std::vector<poperatorT> ops(num_orbitals);
  // Run over occupied components
  int p = 0;
  std::for_each(ops.begin(), ops.end(), [&](auto &operator_p) {
    double mu = sqrt(-2.0 * (ground_energies(p++) + omega + shift));
    operator_p = poperatorT(BSHOperatorPtr3D(world, mu, lo, tol));
  });
  return ops;
  // End timer
}

inline double do_step_restriction(World &world, const vecfuncT &x,
                                  vecfuncT &x_new, const std::string &spin,
                                  const double &maxrotn) {

  double anorm = norm2(world, sub(world, x, x_new));
  
  int nres = 0;
  if (anorm > maxrotn) {
    double s = maxrotn / anorm;
    gaxpy(s, x_new, 1.0 - s, x, true);
  }

  world.gop.fence();
  if (world.rank() == 0)
    print("Norm of vector changes", spin, ": ", anorm);
  return anorm;
}
} // namespace ResponseSolverUtils
