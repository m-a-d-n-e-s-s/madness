// Copyright 2021 Adrian Hurtado
#include <math.h>

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "../chem/NWCHEM.h"  // For nwchem interface
#include "Plot_VTK.h"
#include "TDDFT.h"
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "madness/mra/funcdefaults.h"
#include "molresponse/basic_operators.h"
#include "molresponse/density.h"
#include "molresponse/global_functions.h"
#include "molresponse/load_balance.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/timer.h"
#include "molresponse/x_space.h"

/**
 * @brief Computes two electron interaction functions
 *
 * @param world
 * @param f
 * @param g
 * @param phi
 * @param small
 * @param thresh
 * @param print_level
 * @param xy
 * @return response_space
 */
response_space TDDFT::CreateGamma(World& world,
                                  response_space& f,
                                  response_space& g,
                                  std::vector<real_function_3d>& phi,
                                  double small,
                                  double thresh,
                                  size_t print_level,
                                  std::string xy) {
  // Start timer
  if (print_level >= 1) molresponse::start_timer(world);

  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();

  // The gamma function to be returned, intialized to zero
  response_space gamma(world, m, n);

  // Intermediaries, initialized to zero
  response_space deriv_K(world, m, n);
  response_space deriv_XC(world, m, n);

  // Perturbed coulomb piece
  response_space deriv_J =
      CreateCoulombDerivativeRF(world, f, phi, small, thresh);
  // ResponseFunction deriv_XC=CreateXCDerivative

  // If including any HF exchange:
  if (xcf.hf_exchange_coefficient()) {
    deriv_K = CreateExchangeDerivativeRF(world, f, phi, small, thresh);
  }
  // CreateXcDerivativeOnF(world,f,phi,small,thresh);
  // Get the DFT contribution
  if (xcf.hf_exchange_coefficient() != 1.0) {
    // Get v_xc
    std::vector<real_function_3d> vxc = create_fxc(world, phi, f, g);

    // Apply xc kernel to ground state orbitals
    for (size_t i = 0; i < m; i++) {
      deriv_XC[i] = mul_sparse(world, vxc[i], phi, thresh, false);
    }
    world.gop.fence();
  }

  // Now assemble pieces together to get gamma
  // Spin integration gives 2.0
  gamma = deriv_J * 2.0 + deriv_XC - deriv_K * xcf.hf_exchange_coefficient();

  // Project out groundstate
  QProjector<double, 3> projector(world, Gparams.orbitals);
  for (size_t i = 0; i < m; i++) gamma[i] = projector(gamma[i]);

  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0)
      printf("   Coulomb Deriv matrix for %s components:\n", xy.c_str());
    response_space t = deriv_J * 2.0;
    Tensor<double> temp = expectation(world, f, t);
    if (world.rank() == 0) print(temp);
    if (Rparams.xc == "hf") {
      if (world.rank() == 0)
        printf("   Exchange Deriv matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, deriv_K);
    } else {
      if (world.rank() == 0)
        printf("   XC Deriv matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, deriv_XC);
    }
    if (world.rank() == 0) print(temp);
    if (world.rank() == 0)
      printf("   Gamma matrix for %s components:\n", xy.c_str());
    temp = expectation(world, f, gamma);
    if (world.rank() == 0) print(temp);
  }

  // Basic output
  if (print_level >= 1) molresponse::end_timer(world, "Creating gamma:");

  truncate(world, gamma);

  // Done
  return gamma;
}
/**
 * @brief Compute two electron interaction functions
 *
 * @param world
 * @param f
 * @param phi
 * @param small
 * @param thresh
 * @param print_level
 * @param xy
 * @return response_space
 */
response_space TDDFT::ComputeHf(World& world,
                                const response_space& f,
                                const std::vector<real_function_3d>& phi,
                                double small,
                                double thresh,
                                size_t print_level,
                                std::string xy) {
  // Start timer
  if (print_level >= 1) molresponse::start_timer(world);
  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();

  // The gamma function to be returned, intialized to zero
  response_space H(world, m, n);
  // Intermediaries, initialized to zero
  //
  //
  response_space deriv_J(world, m, n);
  response_space deriv_K(world, m, n);
  response_space deriv_XC(world, m, n);

  // Perturbed coulomb piece
  deriv_J = CreateCoulombDerivativeRF(world, f, phi, small, thresh);
  // Spin integration gives 2.0

  // If including any HF exchange:
  if (xcf.hf_exchange_coefficient()) {
    deriv_K = CreateExchangeDerivativeRF(world, f, phi, small, thresh);
  }
  deriv_K = deriv_K * xcf.hf_exchange_coefficient();
  // Get the DFT contribution
  if (xcf.hf_exchange_coefficient() != 1.0) {
    // Get v_xc
    deriv_XC = CreateXCDerivativeRF(world, f, phi, small, thresh);
  }
  // Now assemble pieces together to get gamma
  // Spin integration gives 2.0
  // J+K+W
  world.gop.fence();
  H = (deriv_J * 2) - deriv_K * xcf.hf_exchange_coefficient() + deriv_XC;

  // Project out groundstate
  QProjector<double, 3> projector(world, Gparams.orbitals);
  for (size_t i = 0; i < m; i++) {
    H[i] = projector(H[i]);
  }

  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0)
      printf("   Coulomb Deriv matrix for %s components:\n", xy.c_str());
    response_space t = deriv_J * 2.0;
    Tensor<double> temp = expectation(world, f, t);
    if (world.rank() == 0) print(temp);
    if (Rparams.xc == "hf") {
      if (world.rank() == 0)
        printf("   Exchange Deriv matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, deriv_K);
    } else {
      if (world.rank() == 0)
        printf("   XC Deriv matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, deriv_XC);
    }
    if (world.rank() == 0) print(temp);
    if (world.rank() == 0)
      printf("   H matrix for %s components:\n", xy.c_str());
    temp = expectation(world, f, H);
    if (world.rank() == 0) print(temp);
  }

  // Basic output
  if (print_level >= 1) molresponse::end_timer(world, "Creating H:");

  truncate(world, H);

  // Done
  return H;
}
/**
 * @brief computes Gf function two electron conjugate interactions
 *
 * @param world
 * @param f
 * @param orbitals
 * @param small
 * @param thresh
 * @param print_level
 * @param xy
 * @return response_space
 */
response_space TDDFT::ComputeGf(World& world,
                                const response_space& f,
                                const std::vector<real_function_3d>& orbitals,
                                double small,
                                double thresh,
                                size_t print_level,
                                std::string xy) {
  // Start a timer
  if (print_level >= 1) molresponse::start_timer(world);

  // Get sizes
  size_t m = f.size();
  size_t n = f[0].size();

  response_space G(world, m, n);

  // Initialize function
  response_space Jdagger(world, m, n);
  response_space Kdagger(world, m, n);
  response_space XCdagger(world, m, n);

  Jdagger = CreateCoulombDerivativeRFDagger(world, f, orbitals, small, thresh);

  // Exchange
  // Determine if including HF exchange
  if (xcf.hf_exchange_coefficient()) {
    Kdagger =
        CreateExchangeDerivativeRFDagger(world, f, orbitals, small, thresh);
  }
  Kdagger = Kdagger * xcf.hf_exchange_coefficient();
  // Determine if DFT potential is needed
  if (xcf.hf_exchange_coefficient() != 1.0) {
    // Get v_xc
    XCdagger = CreateXCDerivativeRFDagger(world, f, orbitals, small, thresh);
  }
  world.gop.fence();

  // Take care of coeficients
  // G = Jdagger * 2.0 + XCdagger - Kdagger * xcf.hf_exchange_coefficient();
  // (J-K)+W
  G = (Jdagger * 2) - Kdagger * xcf.hf_exchange_coefficient() + XCdagger;
  // Project out groundstate
  QProjector<double, 3> projector(world, Gparams.orbitals);
  for (size_t i = 0; i < m; i++) {
    G[i] = projector(G[i]);
  }

  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0) printf("   G coulomb deriv matrix:\n");
    response_space t = Jdagger * 2.0;
    Tensor<double> temp = expectation(world, f, t);
    if (world.rank() == 0) print(temp);
    if (Rparams.xc == "hf") {
      if (world.rank() == 0) printf("   G exchange deriv matrix:\n");
      temp = expectation(world, f, Kdagger);
    } else {
      if (world.rank() == 0) printf("   G XC deriv matrix:\n");
      temp = expectation(world, f, XCdagger);
    }
    if (world.rank() == 0) print(temp);
    printf("   G matrix for %s components:\n", xy.c_str());
    temp = expectation(world, f, G);
    if (world.rank() == 0) print(temp);
  }

  // End timer
  if (print_level >= 1) molresponse::end_timer(world, "   Creating G:");

  // Done
  return G;
}
/**
 * @brief
 *
 * @param world
 * @param rho_omega
 * @param orbital_products
 * @param x
 * @param y
 * @param rhs_x
 * @param rhs_y
 * @param xc
 * @param x_shift
 * @param Gparams
 * @param Rparams
 * @param bsh_x_operators
 * @param bsh_y_operators
 * @param ham_no_diagonal
 * @param iteration
 */
void TDDFT::IterateXY(
    World& world,
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
    size_t iteration) {
  // compute
  size_t m = Rparams.states;
  size_t n = Gparams.num_orbitals;
  double small = Rparams.small;
  double thresh = FunctionDefaults<3>::get_thresh();

  response_space bsh_x_resp(world, m, n);  // Holds wave function corrections
  response_space bsh_y_resp(world, m, n);  // Holds wave function corrections
                                           // coulomb operator
  real_convolution_3d op = CoulombOperator(world, small, thresh);
  Zfunctions Z;
  /* We first compute the 1 electron potentials
   * We compute the the pieces that depend on x response functions
   */
  // x functions
  // V0 applied to x response function
  Z.v0_x = CreatePotential(world, x, xc, Rparams.print_level, "x");
  if (Rparams.print_level == 3) {
    print("norms of v0x");
    print(Z.v0_x.norm2());
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
  Z.v0_y = CreatePotential(world, y, xc, Rparams.print_level, "y");
  Z.y_f_no_diag = y * ham_no_diagonal;  // scale_2d(world, y,
                                        // ham_no_diagonal);
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
        world, x, y, xc, Gparams, Rparams, Rparams.omega != 0.0);
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
  if (world.size() > 1 && (iteration < 2 or iteration % 5 == 0)) {
    // Start a timer
    if (Rparams.print_level >= 1) molresponse::start_timer(world);
    if (world.rank() == 0) print("");  // Makes it more legible
    // (TODO Ask Robert about load balancing)
    LoadBalanceDeux<3> lb(world);
    for (size_t j = 0; j < n; j++) {
      for (size_t k = 0; k < Rparams.states; k++) {
        lb.add_tree(x_response[k][j], lbcost<double, 3>(1.0, 8.0), true);
        lb.add_tree(Z.v0_x[k][j], lbcost<double, 3>(1.0, 8.0), true);
        lb.add_tree(Z.Hx[k][j], lbcost<double, 3>(1.0, 8.0), true);
        // lb.add_tree(Z.Gy[k][j], lbcost<double, 3>(1.0, 8.0), true);
      }
    }
    FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

    if (Rparams.print_level >= 1)
      molresponse::end_timer(world, "Load balancing:");
  }
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  bsh_x_resp = apply(world, bsh_x_operators, Z.Z_x);
  if (Rparams.omega != 0.0) bsh_y_resp = apply(world, bsh_y_operators, Z.Z_y);
  if (Rparams.print_level >= 1) molresponse::end_timer(world, "Apply BSH:");

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
  x_response = bsh_x_resp.copy();
  if (Rparams.omega != 0.0) {
    y_response = bsh_y_resp.copy();
  } else {
    y_response = x_response.copy();
  }

  x_response.truncate_rf();
  y_response.truncate_rf();
}
