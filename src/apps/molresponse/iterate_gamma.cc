
#include <math.h>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/NWChem.h"  // For nwchem interface
#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
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

X_space TDDFT::compute_gamma_full(World& world,
                                  X_space& Chi,
                                  XCOperator<double, 3> xc) {
  size_t m = r_params.n_states();
  size_t n = r_params.num_orbitals();
  // shallow copy
  std::shared_ptr<WorldDCPmapInterface<Key<3>>> oldpmap =
      FunctionDefaults<3>::get_pmap();

  X_space Chi_copy = Chi;
  vecfuncT phi0_copy = ground_orbitals;
  orbital_load_balance(world, ground_orbitals, phi0_copy, Chi, Chi_copy);

  molresponse::start_timer(world);
  X_space gamma(world, m, n);
  // x functions
  // Two ways single vector or vector vector style
  // here I create the orbital products for elctron interaction terms
  vecfuncT phi_phi;
  vecfuncT x_phi;
  vecfuncT y_phi;
  functionT temp_J;

  response_space J(world, m, n);
  response_space k1_x(world, m, n);
  response_space k2_y(world, m, n);
  response_space k1_y(world, m, n);
  response_space k2_x(world, m, n);
  response_space W(world, m, n);
  molresponse::end_timer(world, "Create Zero functions for Gamma calc");

  // apply the exchange kernel to rho if necessary
  molresponse::start_timer(world);
  // Create Coulomb potential on ground_orbitals
  for (size_t b = 0; b < m; b++) {
    temp_J = apply(*coulop, rho_omega[b]);
    temp_J.truncate();
    J[b] = mul(world, temp_J, phi0_copy);
  }
  molresponse::end_timer(world, "J[omega] phi:");

  // Create Coulomb potential on ground_orbitals
  if (xcf.hf_exchange_coefficient() != 1.0) {
    molresponse::start_timer(world);
    std::vector<real_function_3d> Wphi;
    for (size_t b = 0; b < m; b++) {
      Wphi.push_back(xc.apply_xc_kernel(rho_omega[b]));
      W[b] = mul(world, Wphi[b], phi0_copy);
    }
    molresponse::end_timer(world, "XC[omega] phi:");
  }

  molresponse::start_timer(world);
  for (size_t b = 0; b < m; b++) {
    for (size_t p = 0; p < n; p++) {
      phi_phi = mul(world, phi0_copy[p], phi0_copy);
      truncate(world, phi_phi);
      phi_phi = apply(world, *coulop, phi_phi);
      truncate(world, phi_phi);
      k1_x[b][p] = dot(world, Chi_copy.X[b], phi_phi);
      k1_y[b][p] = dot(world, Chi_copy.Y[b], phi_phi);

      // K2
      y_phi = mul(world, phi0_copy[p], Chi_copy.Y[b]);
      truncate(world, y_phi);
      y_phi = apply(world, *coulop, y_phi);
      truncate(world, y_phi);
      k2_y[b][p] = dot(world, Chi_copy.Y[b], phi_phi);

      x_phi = mul(world, phi0_copy[p], Chi_copy.X[b]);
      truncate(world, x_phi);
      x_phi = apply(world, *coulop, x_phi);
      truncate(world, x_phi);
      k2_x[b][p] = dot(world, phi0_copy, x_phi);
    }
  }

  // for each response state we compute the Gamma response functions
  // trucate all response functions
  J.truncate_rf();
  k1_x.truncate_rf();
  k2_x.truncate_rf();
  k1_y.truncate_rf();
  k2_y.truncate_rf();
  W.truncate_rf();
  molresponse::end_timer(world, "K[omega] phi:");

  molresponse::start_timer(world);
  if (r_params.print_level() >= 2) {
    print("-------------------------Gamma Functions ------------------");
    print("2-Electron Potential for Iteration of x");
    PrintResponseVectorNorms(world, J * 2, "J");
    PrintResponseVectorNorms(world, k1_x, "k1_x");
    PrintResponseVectorNorms(world, k2_y, "k2_y");
    PrintResponseVectorNorms(world, k1_x + k2_y, "k1_x+k2_y");
    print("2-Electron Potential for Iteration of y");
    PrintResponseVectorNorms(world, k1_y, "k1_y");
    PrintResponseVectorNorms(world, k2_x, "k2_x");
    PrintResponseVectorNorms(world, k1_y + k2_x, "k1_x+k2_y");
  }
  molresponse::end_timer(world, "Print Response Vector Norms:");
  // update gamma functions
  molresponse::start_timer(world);
  QProjector<double, 3> projector(world, phi0_copy);
  gamma.X = (J * 2) - (k1_x + k2_y) * xcf.hf_exchange_coefficient() + W;
  gamma.Y = (J * 2) - (k1_y + k2_x) * xcf.hf_exchange_coefficient() + W;
  molresponse::end_timer(world, "Add Gamma parts J-K+W  :");

  // project out ground state
  molresponse::start_timer(world);
  for (size_t i = 0; i < m; i++) {
    gamma.X[i] = projector(gamma.X[i]);
    truncate(world, gamma.X[i]);
    gamma.Y[i] = projector(gamma.Y[i]);
    truncate(world, gamma.Y[i]);
  }
  molresponse::end_timer(world, "Project Gamma:");
  molresponse::start_timer(world);
  if (r_params.print_level() >= 2) {
    print("------------------------ Gamma Functions Norms  ------------------");
    print("Gamma X norms");
    print(gamma.X.norm2());
    print("Gamma Y norms");
    print(gamma.Y.norm2());
  }
  molresponse::end_timer(world, " Print Gamma:");

  // put it all together
  // no 2-electron
  molresponse::start_timer(world);
  if (r_params.print_level() >= 2) {
    print("<X ,Gamma(X,Y) Phi>");
    PrintRFExpectation(world, Chi_copy.X, gamma.X, "x", "Gamma)");
    print("<Y ,Gamma_Conjugate(X,Y) Phi>");
    PrintRFExpectation(world, Chi_copy.Y, gamma.Y, "y", "Gamma)");
  }
  // End timer
  molresponse::end_timer(world, "Print Expectation Creating Gamma:");

  molresponse::start_timer(world);
  J.clear();
  k1_x.clear();
  k2_x.clear();
  k1_y.clear();
  k2_y.clear();
  W.clear();
  Chi_copy.clear();

  if (world.size() > 1) {
    FunctionDefaults<3>::set_pmap(oldpmap);  // ! DON'T FORGET !
  }
  molresponse::end_timer(world, "Clear functions and set old pmap");
  // Done
  world.gop.fence();
  return gamma;
  // Get sizes
}

X_space TDDFT::compute_gamma_static(World& world,
                                    X_space& Chi,
                                    XCOperator<double, 3> xc) {
  size_t m = r_params.n_states();
  size_t n = r_params.num_orbitals();
  // shallow copy
  std::shared_ptr<WorldDCPmapInterface<Key<3>>> oldpmap =
      FunctionDefaults<3>::get_pmap();

  X_space Chi_copy = Chi;
  vecfuncT phi0_copy = ground_orbitals;

  orbital_load_balance(world, ground_orbitals, phi0_copy, Chi, Chi_copy);

  molresponse::start_timer(world);
  X_space gamma(world, m, n);
  // x functions
  // here I create the orbital products for elctron interaction terms
  vecfuncT phi_phi;
  vecfuncT x_phi;
  functionT temp_J;

  response_space W(world, m, n);
  response_space J(world, m, n);
  response_space k1_x(world, m, n);
  response_space k2_x(world, m, n);
  molresponse::end_timer(world, "Create Zero functions for Gamma calc");

  // apply the exchange kernel to rho if necessary
  molresponse::start_timer(world);
  // Create Coulomb potential on ground_orbitals
  for (size_t b = 0; b < m; b++) {
    temp_J = apply(*coulop, rho_omega[b]);
    temp_J.truncate();
    J[b] = mul(world, temp_J, phi0_copy);
  }
  molresponse::end_timer(world, "J[omega] phi:");

  // Create Coulomb potential on ground_orbitals
  if (xcf.hf_exchange_coefficient() != 1.0) {
    molresponse::start_timer(world);
    std::vector<real_function_3d> Wphi;
    for (size_t b = 0; b < m; b++) {
      Wphi.push_back(xc.apply_xc_kernel(rho_omega[b]));
      W[b] = mul(world, Wphi[b], phi0_copy);
    }
    molresponse::end_timer(world, "XC[omega] phi:");
  }

  molresponse::start_timer(world);
  for (size_t b = 0; b < m; b++) {
    for (size_t p = 0; p < n; p++) {
      phi_phi = mul(world, phi0_copy[p], phi0_copy);
      truncate(world, phi_phi);
      phi_phi = apply(world, *coulop, phi_phi);
      truncate(world, phi_phi);
      k1_x[b][p] = dot(world, Chi_copy.X[b], phi_phi);
      x_phi = mul(world, phi0_copy[p], Chi_copy.X[b]);
      truncate(world, x_phi);
      x_phi = apply(world, *coulop, x_phi);
      // TODO maybe do not truncate here
      truncate(world, x_phi);
      k2_x[b][p] = dot(world, phi0_copy, x_phi);
    }
  }

  // for each response state we compute the Gamma response functions
  // trucate all response functions
  J.truncate_rf();
  k1_x.truncate_rf();
  k2_x.truncate_rf();
  W.truncate_rf();
  molresponse::end_timer(world, "K[omega] phi:");

  molresponse::start_timer(world);
  if (r_params.print_level() >= 2) {
    print("-------------------------Gamma Functions ------------------");
    print("2-Electron Potential for Iteration of x");
    PrintResponseVectorNorms(world, J * 2, "J");
    PrintResponseVectorNorms(world, k1_x, "k1_x");
    PrintResponseVectorNorms(world, k1_x + k2_x, "k1_x+k2_x");
  }
  molresponse::end_timer(world, "Print Response Vector Norms:");
  molresponse::start_timer(world);
  // update gamma functions
  QProjector<double, 3> projector(world, phi0_copy);
  gamma.X = (J * 2) - (k1_x + k2_x) * xcf.hf_exchange_coefficient() + W;
  molresponse::end_timer(world, "Add Gamma parts J-K+W  :");

  // project out ground state
  molresponse::start_timer(world);
  for (size_t i = 0; i < m; i++) {
    gamma.X[i] = projector(gamma.X[i]);
    truncate(world, gamma.X[i]);
  }
  molresponse::end_timer(world, "Project Gamma:");
  molresponse::start_timer(world);
  if (r_params.print_level() >= 2) {
    print("------------------------ Gamma Functions Norms  ------------------");
    print("Gamma X norms");
    print(gamma.X.norm2());
  }
  molresponse::end_timer(world, " Print Gamma:");

  // put it all together
  // no 2-electron

  molresponse::start_timer(world);
  if (r_params.print_level() >= 2) {
    print("<X ,Gamma(X,Y) Phi>");
    PrintRFExpectation(world, Chi.X, gamma.X, "x", "Gamma)");
  }
  molresponse::end_timer(world, "Print Expectation Creating Gamma:");

  // End timer

  molresponse::start_timer(world);

  J.clear();
  k1_x.clear();
  k2_x.clear();
  W.clear();
  Chi_copy.clear();

  if (world.size() > 1) {
    FunctionDefaults<3>::set_pmap(oldpmap);  // ! DON'T FORGET !
  }

  molresponse::end_timer(world, "Clear functions and set old pmap");
  // Done
  world.gop.fence();
  return gamma;
  // Get sizes
}

X_space TDDFT::compute_gamma_TDA(World& world,
                                 X_space& Chi,
                                 XCOperator<double, 3> xc) {
  // Start a timer
  //
  if (r_params.print_level() >= 1) molresponse::start_timer(world);
  //
  print("-------------------Gamma Functions-------------------");
  print("x_norms in Gamma Functions ");
  print(Chi.X.norm2());
  print("y_norms in Gamma Functions ");
  print(Chi.Y.norm2());

  size_t m = Chi.X.size();
  size_t n = Chi.X.size_orbitals();
  double lo = r_params.lo();
  double thresh = FunctionDefaults<3>::get_thresh();

  X_space gamma(world, m, n);
  // x functions
  real_convolution_3d op = CoulombOperator(world, lo, thresh);
  // Two ways single vector or vector vector style
  // here I create the orbital products for elctron interaction terms

  vector_real_function_3d phi_phi;
  real_function_3d temp_J;

  response_space J(world, m, n);
  response_space k1_x(world, m, n);
  // xc functional
  response_space W(world, m, n);
  std::vector<real_function_3d> Wphi;
  if (xcf.hf_exchange_coefficient() != 1.0) {
    for (size_t i = 0; i < m; i++) {
      Wphi.push_back(xc.apply_xc_kernel(rho_omega[i]));
    }
  }

  for (size_t b = 0; b < m; b++) {
    if (xcf.hf_exchange_coefficient() != 1.0) {
      W[b] = mul(world, Wphi[b], ground_orbitals);
    }
    temp_J = apply(op, rho_omega[b]);
    temp_J.truncate();
    J[b] = mul(world, temp_J,
               ground_orbitals);  // multiply by k
    for (size_t p = 0; p < n; p++) {
      // multiply the kth orbital to vector of y[b] response funtions...apply
      // op
      phi_phi = mul(world, ground_orbitals[p], ground_orbitals);
      truncate(world, phi_phi);
      phi_phi = apply(world, op, phi_phi);
      truncate(world, phi_phi);
      k1_x[b][p] = dot(world, Chi.X[b], phi_phi);
    }
  }
  k1_x.truncate_rf();
  J.truncate_rf();
  W.truncate_rf();
  if (r_params.print_level() >= 2) {
    print("-------------------------Gamma Functions ------------------");
    print("2-Electron Potential for Iteration of x");
    PrintResponseVectorNorms(world, J * 2, "J");
    PrintResponseVectorNorms(world, k1_x, "k1_x");
  }

  QProjector<double, 3> projector(world, ground_orbitals);
  gamma.X = (J * 2) - k1_x * xcf.hf_exchange_coefficient() + W;
  for (size_t i = 0; i < m; i++) {
    gamma.X[i] = projector(gamma.X[i]);
    truncate(world, gamma.X[i]);
  }
  if (r_params.print_level() >= 2) {
    print("------------------------ Gamma Functions Norms  ------------------");
    print("Gamma X norms");
    print(gamma.X.norm2());
  }

  if (r_params.print_level() >= 2) {
    print("<X ,Gamma(X,Y) Phi>");
    PrintRFExpectation(world, Chi.X, gamma.X, "x", "Gamma)");
  }

  // End timer
  if (r_params.print_level() >= 1)
    molresponse::end_timer(world, "   Creating Gamma X:");

  // Done
  world.gop.fence();
  return gamma;
}
