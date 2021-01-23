#include <math.h>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "NWChem.h"  // For nwchem interface
#include "Plot_VTK.h"
#include "TDDFT.h"
#include "TDHF_Basic_Operators2.h"
#include "adrian/ResponseFunction2.h"
#include "adrian/density.h"
#include "adrian/global_functions.h"
#include "adrian/property.h"
#include "adrian/timer.h"
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "load_balance.h"
#include "madness/mra/funcdefaults.h"

void TDHF::IterateXY(
    World& world,
    std::vector<real_function_3d> rho_omega,
    response_space orbital_products,
    response_space& x,
    response_space& y,
    response_space rhs_x,
    response_space rhs_y,
    XCOperator xc,
    double x_shift,
    const GroundParameters& Gparams,
    const ResponseParameters& Rparams,
    std::vector<std::shared_ptr<real_convolution_3d>> bsh_x_operators,
    std::vector<std::shared_ptr<real_convolution_3d>> bsh_y_operators,
    Tensor<double> ham_no_diagonal,
    int iteration) {
  // compute
  int m = Rparams.states;
  int n = Gparams.num_orbitals;
  double small = Rparams.small;
  double thresh = FunctionDefaults<3>::get_thresh();
  response_space bsh_x_resp(world, m, n);  // Holds wave function corrections
  response_space bsh_y_resp(world, m, n);  // Holds wave function corrections

  real_convolution_3d op = CoulombOperator(world, small, thresh);
  Zfunctions Z;
  /* We first compute the 1 electron potentials
   * We compute the the pieces that depend on x response functions
   */
  // x functions
  // V0 applied to x response function
  print("norms of x");
  print(x.norm2());
  Z.v0_x = CreatePotential(world, x, xc, Rparams.print_level, "x");
  if (Rparams.print_level == 3) {
    print("norms of v0x");
    print(Z.v0_x.norm2());
    // + \Delta xp
    print("Value of xshift", x_shift);
  }

  Z.v0_x += (x * x_shift);  // scale(x, x_shifts);
  Z.v0_x.truncate_rf();

  if (Rparams.print_level == 3) {
    print("norms of v0x after scaling");
    print(Z.v0_x.norm2());

    print("norm of psi");
    print(Gparams.orbitals[0].norm2());
    print("trace of psi");
    print(Gparams.orbitals[0].trace());
  }

  // x response scaled by off diagonal ham
  Z.x_f_no_diag = x * ham_no_diagonal;  // scale_2d(world, x, ham_no_diagonal);
  if (Rparams.print_level == 3) {
    print("norms of x scaled by ham no diag");
    print(Z.x_f_no_diag.norm2());
  }
  // If not static we compute the y components
  if (Rparams.omega != 0.0) {
    Z.v0_y = CreatePotential(world, y, xc, Rparams.print_level, "y");
    Z.y_f_no_diag = y * ham_no_diagonal;  // scale_2d(world, y,
                                          // ham_no_diagonal);
  }
  // Some printing for debugging
  if (Rparams.print_level >= 2) {
    { PrintRFExpectation(world, x, Z.v0_x, "x", "V0X"); }

    if (Rparams.omega != 0.0) {
      PrintRFExpectation(world, y, Z.v0_y, "y", "V0Y");
    }
  }
  // Last we compute the 2-electron peices
  //
  //
  if (Rparams.old_two_electron) {
    Z.Hx = ComputeHf(
        world, x, Gparams.orbitals, small, thresh, Rparams.print_level, "x");

    Z.Gy = ComputeGf(
        world, y, Gparams.orbitals, small, thresh, Rparams.print_level, "y");
    if (Rparams.omega != 0.0) {
      Z.Hy = ComputeHf(
          world, y, Gparams.orbitals, small, thresh, Rparams.print_level, "x");

      Z.Gx = ComputeGf(
          world, x, Gparams.orbitals, small, thresh, Rparams.print_level, "y");
    }

    // Z.Z_x = (Z.v0_x - Z.x_f_no_diag + rhs_x) * -2;
    Z.Z_x = (Z.v0_x - Z.x_f_no_diag + Z.Hx + Z.Gy + rhs_x) * -2;
    if (Rparams.omega != 0.0) {
      Z.Z_y = (Z.v0_y - Z.y_f_no_diag + Z.Hy + Z.Gx + rhs_y) * -2;
      // Z.Z_y = (Z.v0_y - Z.y_f_no_diag + rhs_y) * -2;
    }
  } else {
    GammaResponseFunctions gamma = ComputeGammaFunctions(
        world, rho_omega, orbital_products, x, y, xc, Gparams, Rparams);
    // We can use the old algorithm here for testings
    // we then assemble the right hand side vectors
    Z.Z_x = (Z.v0_x - Z.x_f_no_diag + gamma.gamma + rhs_x) * -2;
    // Z.Z_x = Z.v0_x - Z.x_f_no_diag + rhs_x;
    if (Rparams.omega != 0.0) {
      Z.Z_y = (Z.v0_y - Z.y_f_no_diag + gamma.gamma_conjugate + rhs_y) * -2;
      // Z.Z_y = Z.v0_y - Z.y_f_no_diag + rhs_y;
    }
  }
  Z.Z_x.truncate_rf();
  Z.Z_y.truncate_rf();

  if (Rparams.print_level >= 2) {
    if (world.rank() == 0)
      print("   Norms of RHS x components before application of BSH:");
    print_norms(world, Z.Z_x);

    if (Rparams.omega != 0.0) {
      if (world.rank() == 0)
        print("   Norms of RHS y components before application BSH:");
      print_norms(world, Z.Z_x);
    }
  }
  // Load Balancing
  if (world.size() > 1 && ((iteration < 2) or (iteration % 5 == 0))) {
    // Start a timer
    if (Rparams.print_level >= 1) start_timer(world);
    if (world.rank() == 0) print("");  // Makes it more legible
    // (TODO Ask Robert about load balancing)
    LoadBalanceDeux<3> lb(world);
    for (int j = 0; j < n; j++) {
      for (int k = 0; k < Rparams.states; k++) {
        lb.add_tree(x_response[k][j], lbcost<double, 3>(1.0, 8.0), true);
        lb.add_tree(Z.v0_x[k][j], lbcost<double, 3>(1.0, 8.0), true);
        lb.add_tree(Z.Hx[k][j], lbcost<double, 3>(1.0, 8.0), true);
        // lb.add_tree(Z.Gy[k][j], lbcost<double, 3>(1.0, 8.0), true);
      }
    }
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

    if (Rparams.print_level >= 1) end_timer(world, "Load balancing:");
  }
  if (Rparams.print_level >= 1) start_timer(world);

  bsh_x_resp = apply(world, bsh_x_operators, Z.Z_x);
  if (Rparams.omega != 0.0) bsh_y_resp = apply(world, bsh_y_operators, Z.Z_y);
  if (Rparams.print_level >= 1) end_timer(world, "Apply BSH:");

  // Scale by -2.0 (coefficient in eq. 37 of reference paper)
  /*
  for (int i = 0; i < m; i++)
    bsh_x_resp[i] =
        -2 * bsh_x_resp[i];  // * (std::max(1.0, x_norms[i]) * -2.0);
  if (Rparams.omega != 0.0) {
    for (int i = 0; i < m; i++)
      bsh_y_resp[i] =
          -2 * bsh_y_resp[i];  // * (std::max(1.0, x_norms[i]) * -2.0);
  }
  */
  // Debugging output
  if (Rparams.print_level >= 2) {
    if (world.rank() == 0)
      print("   Norms after application of BSH to x components:");
    print_norms(world, bsh_x_resp);

    if (Rparams.omega != 0.0) {
      if (world.rank() == 0)
        print("   Norms after application of BSH to y components:");
      print_norms(world, bsh_y_resp);
    }
  }
  x_response = bsh_x_resp;
  if (Rparams.omega != 0.0) y_response = bsh_y_resp;
  x_response.truncate_rf();
  y_response.truncate_rf();
}
