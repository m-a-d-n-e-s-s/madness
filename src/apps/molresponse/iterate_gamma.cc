
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
#include "basic_operators.h"
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "load_balance.h"
#include "madness/mra/funcdefaults.h"
#include "molresponse/density.h"
#include "molresponse/global_functions.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/timer.h"
#include "molresponse/x_space.h"

GammaResponseFunctions TDDFT::ComputeGammaFunctions(
    World& world,
    response_space& x,
    response_space& y,
    XCOperator xc,
    const GroundParameters& Gparams,
    const ResponseParameters& Rparams,
    bool compute_Y) {
  // Start a timer
  if (Rparams.print_level >= 1) molresponse::start_timer(world);

  size_t m = Rparams.states;
  size_t n = Gparams.num_orbitals;
  double small = Rparams.small;
  double thresh = FunctionDefaults<3>::get_thresh();
  // x functions
  real_convolution_3d op = CoulombOperator(world, small, thresh);
  // Two ways single vector or vector vector style
  // here I create the orbital products for elctron interaction terms
  response_space phi_phi(world, n, n);

  for (size_t k = 0; k < n; k++) {
    // important to do orb[i]*all orbs
    phi_phi[k] =
        apply(world, op, mul(world, Gparams.orbitals[k], Gparams.orbitals));
  }
  phi_phi.truncate_rf();
  print("orbital_products norms");
  print(phi_phi.norm2());

  vector_real_function_3d rho_omega =
      transition_density(world, Gparams.orbitals, x_response, y_response);
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
  for (size_t b = 0; b < m; b++) {
    y_phi.push_back(response_space(world, n, n));
    if (compute_Y) {
      x_phi.push_back(response_space(world, n, n));
    }
  }

  //
  for (size_t b = 0; b < m; b++) {
    for (size_t k = 0; k < n; k++) {
      // multiply the kth orbital to vector of y[b] response funtions...apply op
      // to each product
      // (TODO) //split apply and
      y_phi[b][k] = apply(world, op, mul(world, Gparams.orbitals[k], y[b]));
      if (compute_Y) {
        x_phi[b][k] = apply(world, op, mul(world, Gparams.orbitals[k], x[b]));
      }
    }
  }
  for (size_t b = 0; b < m; b++) {
    if (compute_Y) {
      x_phi[b].truncate_rf();
    }
    y_phi[b].truncate_rf();
  }
  real_function_3d temp_J;
  // for each response state we compute the Gamma response functions
  for (size_t b = 0; b < m; b++) {
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
    for (size_t k = 0; k < n; k++) {
      // Jcoulb*orbital[k]
      //
      Kx[b][k] = dot(world, x[b], phi_phi[k]);
      Ky[b][k] = dot(world, y_phi[b][k], Gparams.orbitals);

      if (compute_Y) {
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
    if (compute_Y) {
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
  print("Gamma norms");
  print(gamma.gamma.norm2());

  if (compute_Y) {
    gamma.gamma_conjugate =
        (J * 2) -
        (Kx_conjugate + Ky_conjugate) * xcf.hf_exchange_coefficient() + W;
  }
  // project out ground state
  for (size_t i = 0; i < m; i++) {
    gamma.gamma[i] = projector(gamma.gamma[i]);
    truncate(world, gamma.gamma[i]);

    if (compute_Y) {
      gamma.gamma_conjugate[i] = projector(gamma.gamma_conjugate[i]);
      truncate(world, gamma.gamma_conjugate[i]);
    }
  }

  // put it all together
  // no 2-electron

  if (Rparams.print_level >= 2) {
    print("<X ,Gamma(X,Y) Phi>");
    PrintRFExpectation(world, x, gamma.gamma, "x", "Gamma)");
    if (compute_Y) {
      print("<Y ,Gamma_Conjugate(X,Y) Phi>");
      PrintRFExpectation(world, y, gamma.gamma_conjugate, "x", "Gamma)");
    }
  }

  // End timer
  if (Rparams.print_level >= 1)
    molresponse::end_timer(world, "   Creating Gamma:");

  // Done
  world.gop.fence();
  return gamma;
  // Get sizes
}