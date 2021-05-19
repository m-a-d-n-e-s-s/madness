// Copyright 2021 Adrian Hurtado
#include <math.h>

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/NWChem.h"  // For nwchem interface
#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
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
  response_space deriv_J = CreateCoulombDerivativeRF(world, f, phi, small, thresh);
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
  QProjector<double, 3> projector(world, ground_orbitals);
  for (size_t i = 0; i < m; i++) gamma[i] = projector(gamma[i]);

  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0) printf("   Coulomb Deriv matrix for %s components:\n", xy.c_str());
    response_space t = deriv_J * 2.0;
    Tensor<double> temp = expectation(world, f, t);
    if (world.rank() == 0) print(temp);
    if (r_params.xc() == "hf") {
      if (world.rank() == 0) printf("   Exchange Deriv matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, deriv_K);
    } else {
      if (world.rank() == 0) printf("   XC Deriv matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, deriv_XC);
    }
    if (world.rank() == 0) print(temp);
    if (world.rank() == 0) printf("   Gamma matrix for %s components:\n", xy.c_str());
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
  QProjector<double, 3> projector(world, ground_orbitals);
  for (size_t i = 0; i < m; i++) {
    H[i] = projector(H[i]);
  }

  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0) printf("   Coulomb Deriv matrix for %s components:\n", xy.c_str());
    response_space t = deriv_J * 2.0;
    Tensor<double> temp = expectation(world, f, t);
    if (world.rank() == 0) print(temp);
    if (r_params.xc() == "hf") {
      if (world.rank() == 0) printf("   Exchange Deriv matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, deriv_K);
    } else {
      if (world.rank() == 0) printf("   XC Deriv matrix for %s components:\n", xy.c_str());
      temp = expectation(world, f, deriv_XC);
    }
    if (world.rank() == 0) print(temp);
    if (world.rank() == 0) printf("   H matrix for %s components:\n", xy.c_str());
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
    Kdagger = CreateExchangeDerivativeRFDagger(world, f, orbitals, small, thresh);
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
  QProjector<double, 3> projector(world, ground_orbitals);
  for (size_t i = 0; i < m; i++) {
    G[i] = projector(G[i]);
  }

  // Debugging output
  if (print_level >= 2) {
    if (world.rank() == 0) printf("   G coulomb deriv matrix:\n");
    response_space t = Jdagger * 2.0;
    Tensor<double> temp = expectation(world, f, t);
    if (world.rank() == 0) print(temp);
    if (r_params.xc() == "hf") {
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
