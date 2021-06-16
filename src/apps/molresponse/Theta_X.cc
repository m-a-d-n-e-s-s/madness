
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

X_space TDDFT::Compute_Theta_X(World& world, X_space& Chi, XCOperator<double, 3> xc, bool compute_Y) {
  print("-------------------Compute Theta X-------------------");
  print("x_norms in Theta X ");
  print(Chi.X.norm2());
  print("y_norms in Theta X ");
  print(Chi.Y.norm2());

  X_space gamma;
  // compute
  if (compute_Y) {
    gamma = compute_gamma_full(world, Chi, xc);
  } else {
    gamma = compute_gamma_static(world, Chi, xc);
  }

  X_space Theta_X = X_space(world, size_states(Chi), size_orbitals(Chi));
  // Compute (V0-ham_no_diag)X
  // V0 appliedo x response function
  response_space v0_X = CreatePotential(world, Chi.X, xc, r_params.print_level(), "x");

  if (r_params.print_level() >= 2) {
    PrintRFExpectation(world, Chi.X, v0_X, "x", "V0X");
  }
  v0_X.truncate_rf();

  if (r_params.print_level() == 3) {
    print("norms of v0x");
    print(v0_X.norm2());
  }
  response_space ham_no_diag_X = Chi.X * ham_no_diag;
  // scale_2d(world, x, ham_no_diagonal);
  if (r_params.print_level() == 3) {
    print("norms of x scaled by ham no diag");
    print(ham_no_diag_X.norm2());
  }
  // Assemble Theta_X and truncate
  Theta_X.X = v0_X - ham_no_diag_X + gamma.X;

  if (compute_Y) {
    // Compute (V0-ham_no_diag)X

    response_space v0_Y = CreatePotential(world, Chi.Y, xc, r_params.print_level(), "y");
    if (r_params.print_level() == 3) {
      print("norms of v0x");
      print(v0_Y.norm2());
    }
    v0_Y.truncate_rf();
    if (r_params.print_level() >= 2) {
      PrintRFExpectation(world, Chi.Y, v0_Y, "y", "V0Y");
    }
    // scale_2d(world, x, ham_no_diagonal);
    response_space ham_no_diag_Y = Chi.Y * ham_no_diag;
    if (r_params.print_level() == 3) {
      print("norms of x scaled by ham no diag");
      print(ham_no_diag_Y.norm2());
    }
    Theta_X.Y = v0_Y - ham_no_diag_Y + gamma.Y;
  }

  return Theta_X;
}
