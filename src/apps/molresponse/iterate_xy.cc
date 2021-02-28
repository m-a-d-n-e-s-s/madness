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
#include "molresponse/response_functions.h"
#include "molresponse/x_space.h"
#include "molresponse/density.h"
#include "molresponse/global_functions.h"
#include "molresponse/property.h"
#include "molresponse/timer.h"
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "load_balance.h"
#include "madness/mra/funcdefaults.h"

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
response_space TDHF::CreateGamma(World& world,
                                 response_space& f,
                                 response_space& g,
                                 std::vector<real_function_3d>& phi,
                                 double small,
                                 double thresh,
                                 int print_level,
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
response_space TDHF::ComputeHf(World& world,
                               const response_space& f,
                               const std::vector<real_function_3d>& phi,
                               double small,
                               double thresh,
                               int print_level,
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
response_space TDHF::ComputeGf(World& world,
                               const response_space& f,
                               const std::vector<real_function_3d>& orbitals,
                               double small,
                               double thresh,
                               int print_level,
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
 * @brief  Computes two electron interaction resonse_spaces 
 * 
 * @param world 
 * @param rho_omega vector of response densities 
 * @param phi_phi  orbital orbital products for exchange
 * @param x x_response states
 * @param y y_response states
 * @param xc  xc functional
 * @param Gparams 
 * @param Rparams 
 * @return GammaResponseFunctions 
 */
GammaResponseFunctions TDHF::ComputeGammaFunctions(
    World& world,
    std::vector<real_function_3d> rho_omega,
    response_space phi_phi,
    response_space& x,
    response_space& y,
    XCOperator xc,
    const GroundParameters& Gparams,
    const ResponseParameters& Rparams) {
  // Start a timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  size_t m = Rparams.states;
  int n = Gparams.num_orbitals;
  double small = Rparams.small;
  double thresh = FunctionDefaults<3>::get_thresh();
  // x functions
  real_convolution_3d op = CoulombOperator(world, small, thresh);

  GammaResponseFunctions gamma;
  gamma.gamma = response_space(world, m, n);
  gamma.gamma_conjugate = response_space(world, m, n);
  response_space J(world, m, n);

  response_space Kx(world, m, n);
  response_space Kx_conjugate(world, m, n);

  response_space Ky(world, m, n);
  response_space Ky_conjugate(world, m, n);

  response_space W(world, m, n);
  std::vector<real_function_3d> Wphi;

  // apply the exchange kernel to rho if necessary
  if (xcf.hf_exchange_coefficient() != 1.0) {
    for (size_t i = 0; i < m; i++) {
      Wphi.push_back(xc.apply_xc_kernel(rho_omega[i]));
    }
  }

  std::vector<response_space> y_phi;
  std::vector<response_space> x_phi;
  for (int b = 0; b < m; b++) {
    y_phi.push_back(response_space(world, n, n));
    if (Rparams.omega != 0.0) {
      x_phi.push_back(response_space(world, n, n));
    }
  }

  //
  for (int b = 0; b < m; b++) {
    for (int k = 0; k < n; k++) {
      // multiply the kth orbital to vector of y[b] response funtions...apply op      // to each product
      // (TODO) //split apply and
      y_phi[b][k] = apply(world, op, mul(world, Gparams.orbitals[k], y[b]));
      if (Rparams.omega != 0.0) {
        x_phi[b][k] = apply(world, op, mul(world, Gparams.orbitals[k], x[b]));
      }
    }
  }
  for (int b = 0; b < m; b++) {
    if (Rparams.omega != 0.0) {
      x_phi[b].truncate_rf();
    }
    y_phi[b].truncate_rf();
  }
  real_function_3d temp_J;
  // for each response state we compute the Gamma response functions
  for (int b = 0; b < m; b++) {
    // apply - save and truncate - multiply
    //
    //
    temp_J = apply(op, rho_omega[b]);
    temp_J.truncate();
    J[b] = mul(world, temp_J,
               Gparams.orbitals);  // multiply by k
    // if xcf not zero
    if (xcf.hf_exchange_coefficient() != 1.0) {
      W[b] = mul(world, Wphi[b], Gparams.orbitals);
    }
    //  compute each of the k response parts...
    for (int k = 0; k < n; k++) {
      // Jcoulb*orbital[k]
      //
      Kx[b][k] = dot(world, x[b], phi_phi[k]);
      Ky[b][k] = dot(world, y_phi[b][k], Gparams.orbitals);

      if (Rparams.omega != 0.0) {
        Ky_conjugate[b][k] = dot(world, y[b], phi_phi[k]);
        Kx_conjugate[b][k] = dot(world, x_phi[b][k], Gparams.orbitals);
      }
    }
  }
  // trucate all response functions
  J.truncate_rf();
  Kx_conjugate.truncate_rf();
  Ky_conjugate.truncate_rf();
  Kx.truncate_rf();
  Ky.truncate_rf();
  W.truncate_rf();

  if (Rparams.print_level >= 2) {
    print("2-Electron Potential for Iteration of x");
    PrintResponseVectorNorms(world, J * 2, "J");
    PrintResponseVectorNorms(world, Kx, "Kx");
    PrintResponseVectorNorms(world, Ky, "Ky");
    PrintResponseVectorNorms(world, Kx + Ky, "Kx+Ky");
    if (Rparams.omega != 0.0) {
      print("2-Electron Potential for Iteration of y");
      PrintResponseVectorNorms(world, Kx_conjugate, "Kx_conjugate");
      PrintResponseVectorNorms(world, Ky_conjugate, "Ky_conjugate");
      PrintResponseVectorNorms(
          world, Kx_conjugate + Ky_conjugate, "Kx_conjugate+Ky_conjugate");
    }
  }
  // update gamma functions
  QProjector<double, 3> projector(world, Gparams.orbitals);
  gamma.gamma = (J * 2) - (Kx + Ky) * xcf.hf_exchange_coefficient() + W;

  if (Rparams.omega != 0.0) {
    gamma.gamma_conjugate =
        (J * 2) -
        (Kx_conjugate + Ky_conjugate) * xcf.hf_exchange_coefficient() + W;
  }
  // project out ground state
  for (size_t i = 0; i < m; i++) {
    gamma.gamma[i] = projector(gamma.gamma[i]);
    truncate(world, gamma.gamma[i]);

    if (Rparams.omega != 0.0) {
      gamma.gamma_conjugate[i] = projector(gamma.gamma_conjugate[i]);
      truncate(world, gamma.gamma_conjugate[i]);
    }
  }

  // put it all together
  // no 2-electron

  if (Rparams.print_level >= 2) {
    print("<X ,Gamma(X,Y) Phi>");
    PrintRFExpectation(world, x, gamma.gamma, "x", "Gamma)");
    if (Rparams.omega != 0.0) {
      print("<Y ,Gamma_Conjugate(X,Y) Phi>");
      PrintRFExpectation(world, y, gamma.gamma_conjugate, "x", "Gamma)");
    }
  }

  // End timer
  if (Rparams.print_level >= 1) molresponse::end_timer(world, "   Creating Gamma:");

  // Done
  world.gop.fence();
  return gamma;
  // Get sizes
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
  size_t m = Rparams.states;
  int n = Gparams.num_orbitals;
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

    if (Rparams.print_level >= 1) molresponse::end_timer(world, "Load balancing:");
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
  if (Rparams.omega != 0.0){
    y_response = bsh_y_resp.copy();
  }else{
    y_response = x_response.copy();
  }
  
  x_response.truncate_rf();
  y_response.truncate_rf();
}
