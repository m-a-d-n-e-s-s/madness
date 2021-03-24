
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

X_space TDDFT::Compute_Theta_X(World& world,
                                X_space& Chi,
                                XCOperator xc,
                                const GroundParameters& Gparams,
                                const ResponseParameters& Rparams,
                                Tensor<double> ham_no_diag,
                                bool compute_Y) {
  // compute
  GammaResponseFunctions gamma = ComputeGammaFunctions(
      world, Chi.X, Chi.Y, xc, Gparams, Rparams, compute_Y);

  X_space Theta_X;
  // Compute (V0-ham_no_diag)X
  // V0 appliedo x response function
  response_space v0_X =
      CreatePotential(world, Chi.X, xc, Rparams.print_level, "x");

  if (Rparams.print_level == 3) {
    print("norms of v0x");
    print(v0_X.norm2());
  }
  v0_X.truncate_rf();
  response_space HX =
      Chi.X * ham_no_diag;  // scale_2d(world, x, ham_no_diagonal);
  if (Rparams.print_level == 3) {
    print("norms of x scaled by ham no diag");
    print(HX.norm2());
  }
  // Assemble Theta_X and truncate
  Theta_X.X = v0_X - HX + gamma.gamma;
  Theta_X.X.truncate_rf();
  if (compute_Y) {
    // Compute (V0-ham_no_diag)X

    response_space v0_Y =
        CreatePotential(world, Chi.Y, xc, Rparams.print_level, "y");
    if (Rparams.print_level == 3) {
      print("norms of v0x");
      print(v0_Y.norm2());
    }
    v0_Y.truncate_rf();
    response_space HY =
        Chi.Y * ham_no_diag;  // scale_2d(world, x, ham_no_diagonal);
    if (Rparams.print_level == 3) {
      print("norms of x scaled by ham no diag");
      print(HY.norm2());
    }
    Theta_X.Y = v0_Y - HY + gamma.gamma_conjugate;
    Theta_X.Y.truncate_rf();
  }
  return Theta_X;
}