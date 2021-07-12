
#include <math.h>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/NWChem.h"  // For nwchem interface
#include "../chem/molecule.h"
#include "Plot_VTK.h"
#include "TDDFT.h"
#include "basic_operators.h"
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "load_balance.h"
#include "madness/mra/funcdefaults.h"
#include "madness/mra/vmra.h"
#include "molresponse/density.h"
#include "molresponse/global_functions.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/timer.h"
#include "molresponse/x_space.h"

// compute exchange |i><i|J|p>
vecfuncT K(vecfuncT& ket, vecfuncT& bra, vecfuncT& vf) {
  World& world = ket[0].world();
  int n = bra.size();
  int nf = ket.size();
  double tol = FunctionDefaults<3>::get_thresh();  /// Important this is
  double mul_tol = 0.0;
  const double lo = 1.e-4;
  const double econv = FunctionDefaults<3>::get_thresh();

  std::shared_ptr<real_convolution_3d> poisson;
  poisson = std::shared_ptr<real_convolution_3d>(
      CoulombOperatorPtr(world, lo, econv));
  /// consistent with Coulomb
  vecfuncT Kf = zero_functions_compressed<double, 3>(world, nf);

  reconstruct(world, bra);
  reconstruct(world, ket);
  reconstruct(world, vf);

  // i-j sym
  for (int i = 0; i < n; ++i) {
    // for each |i> <i|phi>
    vecfuncT psif = mul_sparse(world, bra[i], vf, mul_tol);  /// was vtol
    truncate(world, psif);
    // apply to vector of products <i|phi>..<i|1> <i|2>...<i|N>
    psif = apply(world, *poisson.get(), psif);
    truncate(world, psif);
    // multiply by ket i  <i|phi>|i>: <i|1>|i> <i|2>|i> <i|2>|i>
    psif = mul_sparse(world, ket[i], psif, mul_tol);  /// was vtol
    /// Generalized A*X+Y for vectors of functions ---- a[i] = alpha*a[i] +
    // 1*Kf+occ[i]*psif
    gaxpy(world, double(1.0), Kf, double(1.0), psif);
  }
  truncate(world, Kf, tol);
  return Kf;
}
// sum_i |i><i|J|p> for each p

X_space TDDFT::compute_gamma_full(World& world,
                                  X_space& Chi,
                                  XCOperator<double, 3> xc) {
  size_t m = Chi.num_states();
  size_t n = Chi.num_orbitals();
  //  copy old pmap
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

  X_space J(world, m, n);
  X_space W(world, m, n);
  X_space KX(world, m, n);
  X_space KY(world, m, n);
  molresponse::end_timer(world, "Create Zero functions for Gamma calc");

  // apply the exchange kernel to rho if necessary
  molresponse::start_timer(world);
  // Create Coulomb potential on ground_orbitals
  for (size_t b = 0; b < m; b++) {
    temp_J = apply(*coulop, rho_omega[b]);
    temp_J.truncate();
    J.X[b] = mul(world, temp_J, phi0_copy);
  }
  J.Y = J.X.copy();
  molresponse::end_timer(world, "J[omega] phi:");

  // Create Coulomb potential on ground_orbitals
  if (xcf.hf_exchange_coefficient() != 1.0) {
    molresponse::start_timer(world);
    std::vector<real_function_3d> Wphi;
    for (size_t b = 0; b < m; b++) {
      Wphi.push_back(xc.apply_xc_kernel(rho_omega[b]));
      W.X[b] = mul(world, Wphi[b], phi0_copy);
    }
    molresponse::end_timer(world, "XC[omega] phi:");
  }
  W.Y = W.X.copy();

  molresponse::start_timer(world);
  for (size_t b = 0; b < m; b++) {
    vecfuncT x, y;
    x = Chi_copy.X[b];
    y = Chi_copy.Y[b];
    // |x><i|p>
    KX.X[b] = K(x, phi0_copy, phi0_copy);
    KY.X[b] = K(phi0_copy, y, phi0_copy);
    // |y><i|p>
    KY.Y[b] = K(y, phi0_copy, phi0_copy);
    KX.Y[b] = K(phi0_copy, x, phi0_copy);
    // |i><x|p>
  }
  molresponse::end_timer(world, "K[omega] phi:");

  // for each response state we compute the Gamma response functions
  // trucate all response functions
  J.truncate();
  W.truncate();
  KX.truncate();
  KY.truncate();

  molresponse::start_timer(world);
  if (r_params.print_level() >= 2) {
    print("J(rho1)phi0>");
    J.print_norm2();
    print("K(rho1X)phi0>");
    KX.print_norm2();
    print("K(rho1Y)phi0>");
    KY.print_norm2();
    print("W(rho1)phi0>");
    W.print_norm2();
  }
  // End timer
  molresponse::end_timer(world, "Print Expectation Creating Gamma:");

  // update gamma functions
  molresponse::start_timer(world);
  QProjector<double, 3> projector(world, phi0_copy);
  gamma = (2 * J) - (KX + KY) * xcf.hf_exchange_coefficient() + W;
  molresponse::end_timer(world, "Add Gamma parts J-K+W  :");

  // project out ground state
  molresponse::start_timer(world);
  for (size_t i = 0; i < m; i++) {
    gamma.X[i] = projector(gamma.X[i]);
    gamma.Y[i] = projector(gamma.Y[i]);
  }
  gamma.truncate();

  molresponse::end_timer(world, "Project Gamma:");

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
  KX.clear();
  KY.clear();
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

  X_space W(world, m, n);
  X_space J(world, m, n);
  X_space KX(world, m, n);
  X_space KY(world, m, n);
  molresponse::end_timer(world, "Create Zero functions for Gamma calc");

  // apply the exchange kernel to rho if necessary
  molresponse::start_timer(world);
  // Create Coulomb potential on ground_orbitals
  for (size_t b = 0; b < m; b++) {
    temp_J = apply(*coulop, rho_omega[b]);
    temp_J.truncate();
    J.X[b] = mul(world, temp_J, phi0_copy);
  }

  molresponse::end_timer(world, "J[omega] phi:");

  // Create Coulomb potential on ground_orbitals
  if (xcf.hf_exchange_coefficient() != 1.0) {
    molresponse::start_timer(world);
    std::vector<real_function_3d> Wphi;
    for (size_t b = 0; b < m; b++) {
      Wphi.push_back(xc.apply_xc_kernel(rho_omega[b]));
      W.X[b] = mul(world, Wphi[b], phi0_copy);
    }
    molresponse::end_timer(world, "XC[omega] phi:");
  }

  molresponse::start_timer(world);

  for (size_t b = 0; b < m; b++) {
    vecfuncT x, y;
    x = Chi_copy.X[b];
    y = Chi_copy.Y[b];
    // |x><i|p>
    KX.X[b] = K(x, phi0_copy, phi0_copy);
    // |i><x|p>
    KY.X[b] = K(phi0_copy, y, phi0_copy);
    // |y><i|p>
  }

  molresponse::end_timer(world, "K[omega] phi:");
  // for each response state we compute the Gamma response functions
  // trucate all response functions
  J.truncate();
  KX.truncate();
  KY.truncate();
  W.truncate();

  molresponse::start_timer(world);
  if (r_params.print_level() >= 2) {
    print("J(rho1)phi0>");
    J.print_norm2();
    print("K(rho1X)phi0>");
    KX.print_norm2();
    print("K(rho1Y)phi0>");
    KY.print_norm2();
    print("W(rho1)phi0>");
    W.print_norm2();
  }
  // End timer
  molresponse::end_timer(world, "Print Expectation Creating Gamma:");

  // update gamma functions
  QProjector<double, 3> projector(world, phi0_copy);
  gamma = (2 * J) - (KX + KY) * xcf.hf_exchange_coefficient() + W;
  molresponse::end_timer(world, "Add Gamma parts J-K+W  :");

  // project out ground state
  molresponse::start_timer(world);
  for (size_t i = 0; i < m; i++) {
    gamma.X[i] = projector(gamma.X[i]);
  }
  molresponse::end_timer(world, "Project Gamma:");
  if (r_params.print_level() >= 2) {
    print("------------------------ Gamma Functions Norms  ------------------");
    print("Gamma X norms");
    print(gamma.X.norm2());
  }

  // put it all together
  // no 2-electron

  molresponse::start_timer(world);
  if (r_params.print_level() >= 2) {
    print("<X ,Gamma(X,Y) Phi>");
    PrintRFExpectation(world, Chi_copy.X, gamma.X, "x", "Gamma)");
  }
  molresponse::end_timer(world, "Print Expectation Creating Gamma:");

  // End timer

  molresponse::start_timer(world);

  J.clear();
  KX.clear();
  KY.clear();
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

X_space TDDFT::compute_gamma_tda(World& world,
                                 X_space& Chi,
                                 XCOperator<double, 3> xc) {
  size_t m = Chi.num_states();
  size_t n = Chi.num_orbitals();

  std::shared_ptr<WorldDCPmapInterface<Key<3>>> oldpmap =
      FunctionDefaults<3>::get_pmap();

  X_space Chi_copy = Chi;
  vecfuncT phi0_copy = ground_orbitals;

  orbital_load_balance(world, ground_orbitals, phi0_copy, Chi, Chi_copy);

  molresponse::start_timer(world);
  X_space gamma(world, m, n);
  // x functions

  vector_real_function_3d phi_phi;
  real_function_3d temp_J;

  response_space J(world, m, n);
  response_space k1_x(world, m, n);
  response_space W(world, m, n);
  molresponse::end_timer(world, "Create Zero functions for Gamma calc");

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
    vecfuncT x;
    x = Chi_copy.X[b];
    k1_x[b] = K(x, phi0_copy, phi0_copy);
  }

  molresponse::end_timer(world, "K[omega] phi:");

  k1_x.truncate_rf();
  J.truncate_rf();
  W.truncate_rf();

  if (r_params.print_level() >= 2) {
    molresponse::start_timer(world);
    print("-------------------------Gamma Functions ------------------");
    print("2-Electron Potential for Iteration of x");
    PrintResponseVectorNorms(world, J * 2, "J");
    PrintResponseVectorNorms(world, k1_x, "k1_x");
    molresponse::end_timer(world, "Print Response Vector Norms:");
  }

  molresponse::start_timer(world);
  QProjector<double, 3> projector(world, ground_orbitals);
  gamma.X = (J * 2) - k1_x * xcf.hf_exchange_coefficient() + W;
  molresponse::end_timer(world, "Add Gamma parts J-K+W  :");

  molresponse::start_timer(world);
  for (size_t i = 0; i < m; i++) {
    gamma.X[i] = projector(gamma.X[i]);
    truncate(world, gamma.X[i]);
  }
  molresponse::end_timer(world, "Project Gamma:");

  if (r_params.print_level() >= 2) {
    print("------------------------ Gamma Functions Norms  ------------------");
    print("Gamma X norms");
    print(gamma.X.norm2());
  }

  if (r_params.print_level() >= 2) {
    print("<X ,Gamma(X,Y) Phi>");
    PrintRFExpectation(world, Chi_copy.X, gamma.X, "x", "Gamma)");
  }

  molresponse::start_timer(world);

  J.clear();
  k1_x.clear();
  W.clear();
  Chi_copy.clear();

  if (world.size() > 1) {
    FunctionDefaults<3>::set_pmap(oldpmap);  // ! DON'T FORGET !
  }

  molresponse::end_timer(world, "Clear functions and set old pmap");
  // Done
  world.gop.fence();
  return gamma;
}

// Returns the ground state potential applied to functions f
// (V0 f) V0=(Vnuc+J0-K0+W0)
// J0=J[rho0]
// K0=K[rho0]f
// EXC0=W[rho0]
X_space TDDFT::compute_V0X(World& world,
                           X_space& Chi,
                           XCOperator<double, 3> xc,
                           bool compute_Y) {
  // Start a timer

  size_t m = Chi.num_states();
  size_t n = Chi.num_orbitals();

  X_space V0 = X_space(world, m, n);
  X_space K0 = X_space(world, m, n);

  vecfuncT phi0 = ground_orbitals;
  // v_nuc first
  real_function_3d v_nuc, v_j0, v_k0, v_xc;

  molresponse::start_timer(world);
  if (not r_params.store_potential()) {
    v_nuc = potentialmanager->vnuclear();
    v_nuc.truncate();
  } else {  // Already pre-computed
    v_nuc = stored_v_nuc;
  }
  molresponse::end_timer(world, "Nuclear energy");
  // Coulomb Potential J0*f
  molresponse::start_timer(world);
  if (not r_params.store_potential()) {
    // "a" is the core type
    // scale rho by 2 TODO
    // J^0 x^alpha
    v_j0 = apply(*coulop, rho0);
    v_j0.scale(2.0);
  } else {  // Already pre-computed
    v_j0 = stored_v_coul;
  }
  molresponse::end_timer(world, "Coulomb Potential J[rho0]");

  if (xcf.hf_exchange_coefficient() != 1.0) {
    v_xc = xc.make_xc_potential();
  } else {
    v_xc = Function<double, 3>(
        FunctionFactory<double, 3>(world).fence(false).initial_level(1));
  }

  // Intermediaries

  molresponse::start_timer(world);

  // If including any exact HF exchange
  if (xcf.hf_exchange_coefficient()) {
    for (size_t b = 0; b < m; b++) {
      K0.X[b] = K(phi0, phi0, Chi.X[b]);
      if (compute_Y) {
        K0.Y[b] = K(phi0, phi0, Chi.Y[b]);
      }
    }
  }
  real_function_3d v0 = v_j0 + v_nuc + v_xc;

  V0.X = v0 * Chi.X;
  V0.X += (-1 * K0.X * xcf.hf_exchange_coefficient());

  if (compute_Y) {
    V0.Y = v0 * Chi.Y;
    V0.Y += (-1 * K0.Y * xcf.hf_exchange_coefficient());
  }

  // Basic output

  // Done
  return V0;
}
// kinetic energy operator on response vector
response_space T(World& world, response_space& f) {
  response_space T;  // Fock = (T + V) * orbitals
  real_derivative_3d Dx(world, 0);
  real_derivative_3d Dy(world, 1);
  real_derivative_3d Dz(world, 2);
  // Apply derivatives to orbitals
  f.reconstruct_rf();
  response_space dvx = apply(world, Dx, f);
  response_space dvy = apply(world, Dy, f);
  response_space dvz = apply(world, Dz, f);
  // Apply again for 2nd derivatives
  response_space dvx2 = apply(world, Dx, dvx);
  response_space dvy2 = apply(world, Dy, dvy);
  response_space dvz2 = apply(world, Dz, dvz);
  T = (dvx2 + dvy2 + dvz2) * (-0.5);
  return T;
}
// Returns the ground state fock operator applied to functions f
X_space TDDFT::compute_F0X(World& world,
                           X_space& Chi,
                           XCOperator<double, 3> xc,
                           bool compute_Y) {
  // Debugging output
  size_t m = Chi.num_states();
  size_t n = Chi.num_orbitals();

  X_space F0X = X_space(world, m, n);
  X_space T0X = X_space(world, m, n);
  T0X.X = T(world, Chi.X);
  if (compute_Y) {
    T0X.Y = T(world, Chi.Y);
  }

  X_space V0X = compute_V0X(world, Chi, xc, compute_Y);

  F0X = T0X + V0X;

  // Done
  return F0X;
}
