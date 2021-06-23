
#include <madness/constants.h>
#include <madness/mra/mra.h>
#include <madness/mra/nonlinsol.h>  // The kain solver
#include <madness/mra/operator.h>
#include <math.h>
#include <stdio.h>

#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "../chem/xcfunctional.h"
#include "molresponse/TDDFT.h"
#include "molresponse/basic_operators.h"
#include "molresponse/ground_parameters.h"
#include "molresponse/load_balance.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/response_parameters.h"
#include "molresponse/response_potential.h"
#include "molresponse/timer.h"
#include "molresponse/x_space.h"

X_space TDDFT::Compute_Lambda_X(World& world,
                                X_space& Chi,
                                XCOperator<double, 3> xc,
                                std::string calc_type) {
  // compute
  X_space gamma;

  bool compute_Y = calc_type.compare("full") == 0;

  // compute
  if (calc_type.compare("full") == 0) {
    gamma = compute_gamma_full(world, Chi, xc);
  } else if (calc_type.compare("static") == 0) {
    gamma = compute_gamma_static(world, Chi, xc);
  } else {
    gamma = compute_gamma_TDA(world, Chi, xc);
  }

  X_space Lambda_X = X_space(world, size_states(Chi), size_orbitals(Chi));
  // Compute (V0-ham)X
  // V0 appliedo x response function
  response_space v0_X =
      CreatePotential(world, Chi.X, xc, r_params.print_level(), "x");
  response_space F0_X =
      CreateFock(world, v0_X, Chi.X, r_params.print_level(), "x");

  if (r_params.print_level() == 3) {
    print("norms of v0x");
    print(v0_X.norm2());
  }
  v0_X.truncate_rf();
  response_space HX =
      Chi.X * hamiltonian;  // scale_2d(world, x, ham_no_diagonal);
  if (r_params.print_level() == 3) {
    print("norms of x scaled by ham no diag");
    print(HX.norm2());
  }
  // Assemble Lambda_X and truncate
  Lambda_X.X = F0_X - HX + gamma.X;
  Lambda_X.X.truncate_rf();
  if (compute_Y) {
    // Compute (V0-ham_no_diag)X

    response_space v0_Y =
        CreatePotential(world, Chi.Y, xc, r_params.print_level(), "y");
    response_space F0_Y =
        CreateFock(world, v0_Y, Chi.Y, r_params.print_level(), "y");
    if (r_params.print_level() == 3) {
      print("norms of v0x");
      print(v0_Y.norm2());
    }
    v0_Y.truncate_rf();
    response_space HY =
        Chi.Y * hamiltonian;  // scale_2d(world, x, ham_no_diagonal);
    if (r_params.print_level() == 3) {
      print("norms of x scaled by ham no diag");
      print(HY.norm2());
    }
    Lambda_X.Y = F0_Y - HY + gamma.Y;
    Lambda_X.Y.truncate_rf();
  }
  return Lambda_X;
}
