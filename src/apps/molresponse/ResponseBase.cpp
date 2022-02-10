//
// Created by adrianhurtado on 1/24/22.
//

#include "ResponseBase.hpp"

// Initializes calculation object for both excited state and frequency dependent
// Copies both the response and ground state
/// Constructs the Base Response
/// \param world
/// \param params
ResponseBase::ResponseBase(World &world, const CalcParams &params)
    : r_params(params.response_parameters), molecule(params.molecule),
      ground_calc(params.ground_calculation),
      ground_orbitals(ground_calc.orbitals()),
      ground_energies(ground_calc.get_energies()),
      Chi(world, r_params.num_states(), r_params.num_orbitals()) {

  // Broadcast to all the other nodes
  world.gop.broadcast_serializable(r_params, 0);
  world.gop.broadcast_serializable(molecule, 0);
  xcf.initialize(r_params.xc(), !r_params.spinrestricted(), world,
                 r_params.print_level() >= 3);
  r_params.to_json(j_molresponse);

  // Set the Box Size and Truncate Mode
  FunctionDefaults<3>::set_cubic_cell(-r_params.L(), r_params.L());
  FunctionDefaults<3>::set_truncate_mode(1);
}

/// Checks the polynomial of each function in the ResponseBase
/// First checks the ground state orbitals.  If the orbitals
/// have incorrect k we read from the archive project and truncate
/// We follow by computing the hamiltonian with new orbitals
/// \param world
/// \param thresh
/// \param k
void ResponseBase::check_k(World &world, double thresh, int k) {
  // Boolean to redo ground hamiltonian calculation if
  // ground state orbitals change
  bool redo = false;
  // Verify ground state orbitals have correct k
  if (FunctionDefaults<3>::get_k() != ground_orbitals[0].k()) {
    // Re-read orbitals from the archive (assuming
    // the archive has orbitals stored at a higher
    // k value than what was previously computed
    ground_calc.read(world);
    reconstruct(world, ground_orbitals);
    // Reset correct k (its set in g_params.read)
    FunctionDefaults<3>::set_k(k);
    // Project each ground state to correct k
    for (auto &orbital : ground_orbitals) {
      orbital = project(orbital, FunctionDefaults<3>::get_k(), thresh, false);
    }
    world.gop.fence();
    // Clean up a bit
    truncate(world, ground_orbitals);
    // Now that ground orbitals have correct k lets make the ground density
    // again
    ground_density = make_ground_density(world);
    // Ground state orbitals changed, clear old hamiltonian
    redo = true;
  }
  // Recalculate ground state hamiltonian here
  if (redo or !hamiltonian.has_data()) {
    auto [HAM, HAM_NO_DIAG] = ComputeHamiltonianPair(world);
    // TODO this doesn't seem right...
    hamiltonian = HAM;
    ham_no_diag = HAM_NO_DIAG;
  }

  // If we stored the potential, check that too
  if (r_params.store_potential()) {
    if (FunctionDefaults<3>::get_k() != stored_potential[0][0].k()) {
      // Project the potential into correct k
      for (auto &potential_vector : stored_potential) {
        reconstruct(world, potential_vector);
        for (auto &vi : potential_vector) {
          vi = project(vi, FunctionDefaults<3>::get_k(), thresh, false);
        }
        world.gop.fence();
      }
    }
    if (FunctionDefaults<3>::get_k() != stored_v_coul.k())
      stored_v_coul =
          project(stored_v_coul, FunctionDefaults<3>::get_k(), thresh, false);
    if (FunctionDefaults<3>::get_k() != stored_v_nuc.k())
      stored_v_nuc =
          project(stored_v_nuc, FunctionDefaults<3>::get_k(), thresh, false);
  }
  // Don't forget the mask function as well
  if (FunctionDefaults<3>::get_k() != mask.k()) {
    mask = project(mask, FunctionDefaults<3>::get_k(), thresh, false);
  }
  ::check_k(world, Chi, thresh, k);

  // Make sure everything is done before leaving
  world.gop.fence();
}

/// @brief Computes the Hamiltonian from the ground state  orbitals
///  A side effect of this function is that the stored potentials
///  stored_v_coul and stored V nuc get modified
//   Returns both the hamiltonian and the hamiltonian without diagonal
//
//
/// \param world
/// \return
std::pair<Tensor<double>, Tensor<double>>
ResponseBase::ComputeHamiltonianPair(World &world) const {
  // Basic output
  if (r_params.print_level() >= 1)
    molresponse::start_timer(world);
  auto phi = ground_orbitals;
  // Get sizes
  auto num_orbitals = phi.size();
  // Debugging
  if (r_params.print_level() > 2) {
    Tensor<double> S = matrix_inner(world, phi, phi);
    if (world.rank() == 0)
      print("   Ground state overlap:");
    if (world.rank() == 0)
      print(S);
  }
  // Calculate T
  // Make the derivative operators in each direction
  real_derivative_3d Dx(world, 0);
  real_derivative_3d Dy(world, 1);
  real_derivative_3d Dz(world, 2);

  // Apply derivatives once, and take inner products
  // according to this formula (faster / less noise):
  //  < f | \nabla^2 | f > = - < \nabla f | \nabla f >
  reconstruct(world, phi);
  std::vector<real_function_3d> fx = apply(world, Dx, phi);
  std::vector<real_function_3d> fy = apply(world, Dy, phi);
  std::vector<real_function_3d> fz = apply(world, Dz, phi);
  compress(world, fx, false);
  compress(world, fy, false);
  compress(world, fz, false);
  world.gop.fence();

  // Construct T according to above formula
  // Note: No negative as the formula above
  // has one as well, so they cancel
  Tensor<double> T =
      1.0 / 2.0 *
      (matrix_inner(world, fx, fx) + matrix_inner(world, fy, fy) +
       matrix_inner(world, fz, fz));

  // Construct phiVphi
  // v_nuc first
  // TODO Here I am computing the potential using the potential manager.  Should
  // I do
  // TODO earlier on in set_protocol since this maybe used in the response
  // portion as well?
  real_function_3d v_nuc = potential_manager->vnuclear();
  v_nuc.truncate();

  // V_coul next
  // This does not include final multiplication of each orbital
  // 2 is from integrating out spin
  real_function_3d v_coul = 2.0 * Coulomb(world);

  // Clear old stored potentials  Made this mutable
  stored_v_coul.clear();
  stored_v_nuc.clear();

  // If storing potentials, save them here
  if (r_params.store_potential()) {
    stored_v_nuc = copy(v_nuc);
    stored_v_coul = copy(v_coul);
  }

  // v_nuc comes out negative from potential manager, so add it
  real_function_3d v = v_coul + v_nuc;

  // Apply phiVphi to f functions
  std::vector<real_function_3d> vf = v * phi;
  // Clear stored_potential
  stored_potential.clear();
  // ALWAYS DO THIS FOR THE STORED POTENTIAL!!
  // exchange last
  // 'small memory' algorithm from SCF.cc
  auto op = coulop;

  auto Kphi = zero_functions_compressed<double, 3>(world, num_orbitals);

  for (const auto &phi_i : phi) {
    /// Multiplies a function against a vector of functions using sparsity of a
    /// and v[i] --- q[i] = a * v[i]
    auto psif =
        mul_sparse(world, phi_i, phi, FunctionDefaults<3>::get_thresh());
    truncate(world, psif);
    psif = apply(world, *op, psif);
    truncate(world, psif);
    // Save the potential here if we are saving it
    if (r_params.store_potential()) {
      stored_potential.push_back(psif);
    }
    psif = mul_sparse(world, phi_i, psif, FunctionDefaults<3>::get_thresh());
    gaxpy(world, 1.0, Kphi, 1.0, psif);
  }
  // Only use the exchange above if HF:
  Tensor<double> phiVphi;
  real_function_3d v_xc;
  if (r_params.xc() == "hf") {
    // Construct phiVphi
    phiVphi = matrix_inner(world, phi, vf) - matrix_inner(world, phi, Kphi);
  } else { // DFT

    XCOperator<double, 3> xcop = make_xc_operator(world);

    real_function_3d v_xc = xcop.make_xc_potential();
    v = v + v_xc;
    std::vector<real_function_3d> vf = v * phi;
    if ((*xcop.xc).hf_exchange_coefficient() > 0.0) {
      // XCOperator<double,3>  has member variable xc, which is an
      // xcfunctional which has the hf_exchange_coeff we need here
      gaxpy(world, 1.0, vf, -(*xcop.xc).hf_exchange_coefficient(), Kphi);
    }
    phiVphi = matrix_inner(world, phi, vf);
  }
  // Now create the new_hamiltonian
  auto new_hamiltonian = T + phiVphi;

  for (int64_t i = 0; i < new_hamiltonian.dim(0); i++) {
    for (int64_t j = i + 1; j < new_hamiltonian.dim(1); j++) {
      //      print(i, j);
      //      print(xAx(i, j));
      //     print(xAx(j, i));
      new_hamiltonian(j, i) = new_hamiltonian(i, j);
    }
  }
  double traceOfHamiltonian(0);
  for (int64_t i = 0; i < new_hamiltonian.dim(0); i++) {
    traceOfHamiltonian += new_hamiltonian(i, i);
  }
  print("Trace of Hamiltonian");
  print(traceOfHamiltonian);
  // Save a matrix that is
  // (T+phiVphi) - Lambda * eye
  // Copy new_hamiltonian and zero the diagonal
  auto new_hamiltonian_no_diag = copy(new_hamiltonian);
  for (size_t i = 0; i < num_orbitals; i++)
    new_hamiltonian_no_diag(i, i) = 0.0;

  // Debug output
  if (r_params.print_level() >= 2 and world.rank() == 0) {
    print("   Ground state new_hamiltonian:");
    print(new_hamiltonian);
  }
  // End timer
  if (r_params.print_level() >= 1)
    molresponse::end_timer(world, "   Create grnd ham:");
  return {new_hamiltonian, new_hamiltonian_no_diag};
}

functionT ResponseBase::make_ground_density(World &world) const {

  std::vector<real_function_3d> vsq = square(world, ground_orbitals);
  compress(world, vsq);
  functionT rho = factoryT(world);
  rho.compress();
  for (const auto &phi_squared : vsq)
    rho.gaxpy(1.0, phi_squared, 1.0, false);
  world.gop.fence();
  vsq.clear();
  return rho;
}
// @brief Calculates ground state coulomb potential function
//
//  The coulomb potential of the ground state is
//  \f$ f=\int \sum \frac{\phi_i(r)\phi(r)}{|r-r'|}dr
//
/// \param world
/// \return
real_function_3d ResponseBase::Coulomb(World &world) const {
  return apply(*coulop, ground_density).truncate();
}

// TODO  Create apply_operator<T>(f)  generalized function in place of coulomb

XCOperator<double, 3> ResponseBase::make_xc_operator(World &world) const {
  return XCOperator<double, 3>(world, r_params.xc(), false, ground_density,
                               ground_density);
}

// Save the current response calculation
void ResponseBase::save(World &world, const std::string &name) {
  // Archive to write everything to
  archive::ParallelOutputArchive ar(world, name.c_str(), 1);
  // Just going to enforce 1 io server

  // Saving, in this order;
  //  string           ground-state archive name (garch_name)
  //  bool             TDA flag
  // size_t                number of ground state orbitals (n)
  // size_t                number of excited state orbitals (m)
  //  Tensor<double>   energies of m x-components
  //  for i from 0 to m-1
  //     for j from 0 to n-1
  //        Function<double,3> x_response[i][j]
  //  (If TDA flag == True)
  //  (Tensor<double>  energies of m y-components    )
  //  (for i from 0 to m-1                       )
  //  (   for j from 0 to n-1                    )
  //  (      Function<double,3> y_response[i][j] )
  ar &r_params.archive();
  ar &r_params.tda();
  ar &r_params.num_orbitals();
  ar &r_params.num_states();
  ar &omega;

  for (size_t i = 0; i < r_params.num_states(); i++)
    for (size_t j = 0; j < r_params.num_orbitals(); j++)
      ar &Chi.X[i][j];
  if (not r_params.tda()) {
    for (size_t i = 0; i < r_params.num_states(); i++)
      for (size_t j = 0; j < r_params.num_orbitals(); j++)
        ar &Chi.Y[i][j];
  }
}

// Load a response calculation
void ResponseBase::load(World &world, const std::string &name) {
  // The archive to read from
  archive::ParallelInputArchive ar(world, name.c_str());

  // Reading in, in this order;
  //  string           ground-state archive name (garch_name)
  //  bool             TDA flag
  // size_t                number of ground state orbitals (n)
  // size_t                number of excited state orbitals (m)
  //  Tensor<double>   energies of m x-components
  //  for i from 0 to m-1
  //     for j from 0 to n-1
  //        Function<double,3> x_response[i][j]
  //  (If TDA flag == True)
  //  (Tensor<double>  energies of m y-components    )
  //  (for i from 0 to m-1                       )
  //  (   for j from 0 to n-1                    )
  //  (      Function<double,3> y_response[i][j] )

  ar &r_params.archive();
  ar &r_params.tda();
  ar &r_params.num_orbitals();
  ar &r_params.num_states();
  ar &omega;

  Chi = X_space(world, r_params.num_states(), r_params.num_orbitals());

  for (size_t i = 0; i < r_params.num_states(); i++)
    for (size_t j = 0; j < r_params.num_orbitals(); j++)
      ar &Chi.X[i][j];
  world.gop.fence();

  if (not r_params.tda()) {
    for (size_t i = 0; i < r_params.num_states(); i++)
      for (size_t j = 0; j < r_params.num_orbitals(); j++)
        ar &Chi.Y[i][j];
    world.gop.fence();
  }
}
vecfuncT ResponseBase::make_density(World &world, const X_space &chi) const {
  molresponse::start_timer(world);
  vecfuncT density;
  auto calc_type = r_params.calc_type();
  if (calc_type == "full") {
    density = transition_density(world, ground_orbitals, chi.X, chi.Y);
  } else if (calc_type == "static") {
    density = transition_density(world, ground_orbitals, chi.X, chi.X);
  } else {
    density = transition_densityTDA(world, ground_orbitals, chi.X);
  }
  molresponse::end_timer(world, "Make density omega");
  world.gop.fence();
  return density;
}

void ResponseBase::load_balance_chi(World &world) {
  molresponse::start_timer(world);
  if (world.size() == 1)
    return;

  LoadBalanceDeux<3> lb(world);
  real_function_3d v_nuclear;
  v_nuclear = potential_manager->vnuclear();
  for (size_t i = 0; i < Chi.X.size(); ++i) {
    for (size_t j = 0; j < Chi.X.size_orbitals(); ++j) {
      lb.add_tree(Chi.X[i][j], lbcost<double, 3>(1.0, 8.0), false);
    }
  }
  if (r_params.omega() != 0) {
    for (size_t i = 0; i < Chi.X.size(); ++i) {
      for (size_t j = 0; j < Chi.X.size_orbitals(); ++j) {
        lb.add_tree(Chi.Y[i][j], lbcost<double, 3>(1.0, 8.0), false);
      }
    }
  }

  world.gop.fence();

  FunctionDefaults<3>::redistribute(
      world,
      lb.load_balance(r_params.loadbalparts())); // 6.0 needs retuning after

  world.gop.fence();
  molresponse::end_timer(world, "Load balancing");
}

std::vector<poperatorT>
ResponseBase::make_bsh_operators_response(World &world, double &shift,
                                          double &omega) const {
  if (r_params.print_level() >= 1)
    molresponse::start_timer(world);
  double tol = FunctionDefaults<3>::get_thresh();
  // Sizes inferred from ground and omega
  size_t num_orbitals = ground_energies.size(); // number of orbitals
  std::vector<poperatorT> ops(num_orbitals);
  // Run over occupied components
  for (size_t p = 0; p < num_orbitals; p++) {
    double mu = sqrt(-2.0 * (ground_energies(p) + omega + shift));
    ops[p] = poperatorT(BSHOperatorPtr3D(world, mu, r_params.lo(), tol));
  }
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "make bsh operators response");
  }
  return ops;
  // End timer
}
X_space ResponseBase::compute_theta_X(World &world, const X_space &chi,
                                      XCOperator<double, 3> xc,
                                      std::string calc_type) const {
  bool compute_Y = calc_type.compare("full") == 0;
  X_space Theta_X = X_space(world, chi.num_states(), chi.num_orbitals());
  // compute
  X_space V0X = compute_V0X(world, chi, xc, compute_Y);

  V0X.truncate();
  if (r_params.print_level() >= 20) {
    print("---------------Theta ----------------");
    print("<X|V0|X>");
    print(inner(chi, V0X));
  }

  X_space E0X(world, chi.num_states(), chi.num_orbitals());
  if (r_params.localize().compare("canon") == 0) {
    E0X = chi.copy();
    E0X.truncate();
    E0X.X = E0X.X * ham_no_diag;
    if (compute_Y) {
      E0X.Y = E0X.Y * ham_no_diag;
    }

    E0X.truncate();
  }

  if (r_params.print_level() >= 20) {
    print("<X|(E0-diag(E0)|X>");
    print(inner(chi, E0X));
  }

  X_space gamma;
  // compute
  if (calc_type == "full") {
    gamma =
        compute_gamma_full(world, {chi, ground_orbitals, ground_orbitals}, xc);
  } else if (calc_type == "static") {
    gamma = compute_gamma_static(world, {chi, ground_orbitals, ground_orbitals},
                                 xc);
  } else {
    gamma =
        compute_gamma_tda(world, {chi, ground_orbitals, ground_orbitals}, xc);
  }

  Theta_X = (V0X - E0X) + gamma;
  Theta_X.truncate();

  if (r_params.print_level() >= 20) {
    print("<X|Theta|X>");
    print(inner(chi, Theta_X));
  }

  return Theta_X;
}


X_space
ResponseBase::compute_gamma_full(World &world, const gamma_orbitals &density,
                                 const XCOperator<double, 3> &xc) const {
  std::shared_ptr<WorldDCPmapInterface<Key<3>>> old_pmap =
      FunctionDefaults<3>::get_pmap();

  auto [d_alpha, phi0, vf] =
      orbital_load_balance(world, density, r_params.loadbalparts());

  size_t num_states = d_alpha.num_states();
  size_t num_orbitals = d_alpha.num_orbitals();

  size_t m = d_alpha.num_states();
  size_t n = d_alpha.num_orbitals();
  //  copy old pmap

  truncate(world, phi0);
  truncate(world, vf);
  d_alpha.truncate();

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
  response_space j_x(world, m, n);
  response_space j_y(world, m, n);

  X_space W(world, m, n);
  X_space KX(world, m, n);
  X_space KY(world, m, n);
  molresponse::end_timer(world, "Create Zero functions for Gamma calc");

  // apply the exchange kernel to rho if necessary
  molresponse::start_timer(world);
  // Create Coulomb potential on ground_orbitals
  functionT rho_x_b;
  functionT rho_y_b;

  // note that x can refer to x or y
  auto compute_j = [&, &phi0 = phi0](auto dx) {
    // compute density with response function dx and orbitals phi0
    auto rho_x_b = dot(world, dx, phi0);
    rho_x_b.truncate();
    // apply the coulomb operator to rho_b
    rho_x_b = apply(*coulop, rho_x_b);
    return mul(world, rho_x_b, phi0);
  };

  // compute j_x = op(rho_x)*phi0
  std::transform(d_alpha.X.begin(), d_alpha.X.end(), j_x.begin(), compute_j);
  // compute j_y = op(rho_y)*phi0
  std::transform(d_alpha.Y.begin(), d_alpha.Y.end(), j_y.begin(), compute_j);

  J.X = j_x + j_y;
  // TODO is copy better than adding? probably?
  // J.Y=j_x+j_y;
  J.Y = J.X.copy();

  molresponse::end_timer(world, "J[omega] phi:");

  // Create Coulomb potential on ground_orbitals
  if (xcf.hf_exchange_coefficient() != 1.0) {
    auto rho = transition_density(world, phi0, d_alpha.X, d_alpha.X);
    auto compute_wx = [&, &phi0 = phi0](auto rho_alpha) {
      auto xc_rho = xc.apply_xc_kernel(rho_alpha);
      return mul(world, xc_rho, phi0);
    };
    molresponse::start_timer(world);
    std::transform(rho.begin(), rho.end(), W.X.begin(), compute_wx);
    W.Y = W.X.copy();
    molresponse::end_timer(world, "XC[omega] phi:");
  }

  molresponse::start_timer(world);

  for (size_t b = 0; b < m; b++) {
    vecfuncT x, y;
    x = d_alpha.X[b];
    y = d_alpha.Y[b];
    // |x><i|p>
    KX.X[b] = K(x, phi0, vf);
    KY.X[b] = K(phi0, y, vf);
    // |y><i|p>
    KY.Y[b] = K(y, phi0, vf);
    KX.Y[b] = K(phi0, x, vf);
    // |i><x|p>
  }
  molresponse::end_timer(world, "K[omega] phi:");
  J.truncate();
  KX.truncate();
  KY.truncate();
  W.truncate();
  // for each response state we compute the Gamma response functions
  // trucate all response functions

  // update gamma functions
  molresponse::start_timer(world);
  gamma = (2 * J) - (KX + KY) * xcf.hf_exchange_coefficient() + W;
  molresponse::end_timer(world, "Add Gamma parts J-K+W  :");

  // project out ground state
  molresponse::start_timer(world);
  QProjector<double, 3> projector(world, phi0);
  for (size_t i = 0; i < m; i++) {
    gamma.X[i] = projector(gamma.X[i]);
    gamma.Y[i] = projector(gamma.Y[i]);
  }

  molresponse::end_timer(world, "Project Gamma:");
  if (r_params.print_level() >= 3) {
      molresponse::start_timer(world);
      print("inner <X|J|X>");
      print(inner(d_alpha, J));
      print("inner <X|KX|X>");
      print(inner(d_alpha, KX));
      print("inner <X|KY|X>");
      print(inner(d_alpha, KY));
      print("inner <X|K|X>");
      X_space K = KX + KY;
      print(inner(d_alpha, K));
      print("inner <X|W|X>");
      print(inner(d_alpha, W));
      print("inner <X|Gamma|X>");
      print(inner(d_alpha, gamma));

      molresponse::end_timer(world, "Print Expectation Creating Gamma:");
  }
  // put
  // put it all together
  // no 2-electron
  // End timer

  molresponse::start_timer(world);
  J.clear();
  j_x.clear();
  j_y.clear();
  KX.clear();
  KY.clear();
  W.clear();

  d_alpha.clear();
  phi0.clear();
  vf.clear();

  if (world.size() > 1) {
    FunctionDefaults<3>::set_pmap(old_pmap); // ! DON'T FORGET !
  }
  molresponse::end_timer(world, "Clear functions and set old pmap");
  // Done
  world.gop.fence();
  return gamma;
  // Get sizes
}

X_space
ResponseBase::compute_gamma_static(World &world, const gamma_orbitals &density,
                                   const XCOperator<double, 3> &xc) const {

  // X contains the response vector that makes up the response density at a
  // given order

  std::shared_ptr<WorldDCPmapInterface<Key<3>>> old_pmap =
      FunctionDefaults<3>::get_pmap();

  auto [d_alpha, phi0, vf] =
      orbital_load_balance(world, density, r_params.loadbalparts());

  size_t num_states = d_alpha.num_states();
  size_t num_orbitals = d_alpha.num_orbitals();
  // shallow copy

  molresponse::start_timer(world);
  X_space gamma(world, num_states, num_orbitals);
  // x functions
  // here I create the orbital products for elctron interaction terms
  vecfuncT phi_phi;
  vecfuncT x_phi;
  functionT temp_J;

  X_space W(world, num_states, num_orbitals);
  X_space J(world, num_states, num_orbitals);
  X_space KX(world, num_states, num_orbitals);
  X_space KY(world, num_states, num_orbitals);
  molresponse::end_timer(world, "Create Zero functions for Gamma calc");

  // apply the exchange kernel to rho if necessary
  molresponse::start_timer(world);

  auto rho = transition_density(world, phi0, d_alpha.X, d_alpha.X);
  // Create Coulomb potential on ground_orbitals

  auto compute_jx = [&, &phi0 = phi0](auto rho_alpha) {
    auto temp_J = apply(*coulop, rho_alpha);
    temp_J.truncate();
    return mul(world, temp_J, phi0);
  };

  std::transform(rho.begin(), rho.end(), J.X.begin(), compute_jx);

  molresponse::end_timer(world, "J[omega] phi:");

  // Create Coulomb potential on ground_orbitals
  if (xcf.hf_exchange_coefficient() != 1.0) {
    auto compute_wx = [&, &phi0 = phi0](auto rho_alpha) {
      auto xc_rho = xc.apply_xc_kernel(rho_alpha);
      return mul(world, xc_rho, phi0);
    };
    molresponse::start_timer(world);
    // for every transition density apply the exchange kernel and multiply the
    // vector of orbitals
    std::transform(rho.begin(), rho.end(), W.X.begin(), compute_wx);
    W.Y = W.X.copy();

    molresponse::end_timer(world, "XC[omega] phi:");
  }

  molresponse::start_timer(world);

  for (size_t b = 0; b < num_states; b++) {
    vecfuncT x, y;
    x = d_alpha.X[b];
    y = d_alpha.Y[b];
    // |x><i|p>
    KX.X[b] = K(x, phi0, vf);
    // |i><x|p>
    KY.X[b] = K(phi0, y, vf);
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
  // update gamma functions
  gamma = (2 * J) - (KX + KY) * xcf.hf_exchange_coefficient() + W;
  molresponse::end_timer(world, "Add Gamma parts J-K+W  :");

  // project out ground state
  molresponse::start_timer(world);
  QProjector<double, 3> projector(world, phi0);
  for (size_t i = 0; i < num_states; i++) {
    gamma.X[i] = projector(gamma.X[i]);
  }
  molresponse::end_timer(world, "Project Gamma:");

  molresponse::start_timer(world);

  J.clear();
  KX.clear();
  KY.clear();
  W.clear();
  d_alpha.clear();
  phi0.clear();
  vf.clear();

  if (world.size() > 1) {
    FunctionDefaults<3>::set_pmap(old_pmap); // ! DON'T FORGET !
  }

  molresponse::end_timer(world, "Clear functions and set old pmap");
  // Done
  world.gop.fence();
  return gamma;
  // Get sizes
}

X_space ResponseBase::compute_gamma_tda(World &world,
                                        const gamma_orbitals &density,
                                        const XCOperator<double, 3> &xc) const {

  auto [d_alpha, phi0, vf] =
      orbital_load_balance(world, density, r_params.loadbalparts());
  std::shared_ptr<WorldDCPmapInterface<Key<3>>> oldpmap =
      FunctionDefaults<3>::get_pmap();
  size_t num_states = d_alpha.num_states();
  size_t num_orbitals = d_alpha.num_orbitals();

  molresponse::start_timer(world);
  X_space gamma(world, num_states, num_orbitals);
  // x functions

  vector_real_function_3d phi_phi;
  real_function_3d temp_J;

  response_space J(world, num_states, num_orbitals);
  response_space k1_x(world, num_states, num_orbitals);
  response_space W(world, num_states, num_orbitals);
  molresponse::end_timer(world, "Create Zero functions for Gamma calc");

  molresponse::start_timer(world);

  auto rho = transition_densityTDA(world, phi0, d_alpha.X);

  auto compute_jx = [&, &phi0 = phi0](auto rho_alpha) {
    auto temp_J = apply(*coulop, rho_alpha);
    temp_J.truncate();
    return mul(world, temp_J, phi0);
  };

  std::transform(rho.begin(), rho.end(), J.begin(), compute_jx);
  molresponse::end_timer(world, "J[omega] phi:");

  // Create Coulomb potential on ground_orbitals
  if (xcf.hf_exchange_coefficient() != 1.0) {
    auto compute_wx = [&, &phi0 = phi0](auto rho_alpha) {
      auto xc_rho = xc.apply_xc_kernel(rho_alpha);
      return mul(world, xc_rho, phi0);
    };
    molresponse::start_timer(world);
    // for every transition density apply the exchange kernel and multiply the
    // vector of orbitals
    std::transform(rho.begin(), rho.end(), W.begin(), compute_wx);
    W = W.copy();

    molresponse::end_timer(world, "XC[omega] phi:");
  }

  molresponse::start_timer(world);

  for (size_t b = 0; b < num_states; b++) {
    vecfuncT x;
    x = d_alpha.X[b];
    k1_x[b] = K(x, phi0, vf);
  }

  molresponse::end_timer(world, "K[omega] phi:");

  k1_x.truncate_rf();
  J.truncate_rf();
  W.truncate_rf();

  molresponse::start_timer(world);
  QProjector<double, 3> projector(world, ground_orbitals);
  gamma.X = (J * 2) - k1_x * xcf.hf_exchange_coefficient() + W;
  molresponse::end_timer(world, "Add Gamma parts J-K+W  :");

  molresponse::start_timer(world);
  for (size_t i = 0; i < num_states; i++) {
    gamma.X[i] = projector(gamma.X[i]);
    truncate(world, gamma.X[i]);
  }
  molresponse::end_timer(world, "Project Gamma:");

  if (r_params.print_level() >= 2) {
    print("------------------------ Gamma Functions Norms  ------------------");
    print("Gamma X norms");
    print(gamma.X.norm2());
  }

  molresponse::start_timer(world);

  J.clear();
  k1_x.clear();
  W.clear();

  d_alpha.clear();
  phi0.clear();
  vf.clear();

  if (world.size() > 1) {
    FunctionDefaults<3>::set_pmap(oldpmap); // ! DON'T FORGET !
  }

  molresponse::end_timer(world, "Clear functions and set old pmap");
  // Done
  world.gop.fence();
  return gamma;
}
X_space ResponseBase::compute_lambda_X(World &world, const X_space &chi,
                                       XCOperator<double, 3> xc,
                                       const std::string &calc_type) const {
  // compute
  bool compute_Y = calc_type.compare("full") == 0;

  X_space Lambda_X; // = X_space(world, chi.num_states(), chi.num_orbitals());
  X_space F0X = compute_F0X(world, chi, xc, compute_Y);
  X_space Chi_truncated = chi.copy();
  Chi_truncated.truncate();
  if (r_params.print_level() >= 20) {
    print("---------------Lambda ----------------");
    print("<X|F0|X>");
    print(inner(Chi_truncated, F0X));
  }
  // put it all together

  X_space E0X = Chi_truncated.copy();
  E0X.truncate();
  E0X.X = E0X.X * hamiltonian;

  if (compute_Y) {
    E0X.Y = E0X.Y * hamiltonian;
  }
  if (r_params.print_level() >= 20) {
    print("<X|E0|X>");
    print(inner(Chi_truncated, E0X));
  }

  // put it all together
  X_space gamma;

  // compute
  if (calc_type == "full") {
    gamma =
        compute_gamma_full(world, {chi, ground_orbitals, ground_orbitals}, xc);
  } else if (calc_type == "static") {
    gamma = compute_gamma_static(world, {chi, ground_orbitals, ground_orbitals},
                                 xc);
  } else {
    gamma =
        compute_gamma_tda(world, {chi, ground_orbitals, ground_orbitals}, xc);
  }
  if (r_params.print_level() >= 20) {
    print("<X|Gamma|X>");
    print(inner(Chi_truncated, gamma));
  }

  Lambda_X = (F0X - E0X) + gamma;
  if (r_params.print_level() >= 20) {
    print("<X|Lambda not truncated|X>");
    print(inner(Chi_truncated, Lambda_X));
  }
  Lambda_X.truncate();

  if (r_params.print_level() >= 20) {
    print("<X|Lambda_truncated|X>");
    print(inner(Chi_truncated, Lambda_X));
  }

  return Lambda_X;
}
// Returns the ground state potential applied to functions f
// (V0 f) V0=(Vnuc+J0-K0+W0)
// J0=J[ground_density]
// K0=K[ground_density]f
// EXC0=W[ground_density]
X_space ResponseBase::compute_V0X(World &world, const X_space &X,
                                  const XCOperator<double, 3> &xc,
                                  bool compute_Y) const {
  // Start a timer

  size_t m = X.num_states();
  size_t n = X.num_orbitals();

  X_space V0 = X_space(world, m, n);
  X_space K0 = X_space(world, m, n);

  X_space Chi_copy = X;
  vecfuncT phi0_copy = ground_orbitals;
  Chi_copy.truncate();
  truncate(world, phi0_copy);
  // v_nuc first
  real_function_3d v_nuc, v_j0, v_k0, v_xc;

  molresponse::start_timer(world);
  if (not r_params.store_potential()) {
    v_nuc = potential_manager->vnuclear();
    v_nuc.truncate();
  } else { // Already pre-computed
    v_nuc = stored_v_nuc;
  }
  molresponse::end_timer(world, "Nuclear energy");
  // Coulomb Potential J0*f
  molresponse::start_timer(world);
  if (not r_params.store_potential()) {
    // "a" is the core type
    // scale rho by 2 TODO
    // J^0 x^alpha
    v_j0 = apply(*coulop, ground_density);
    v_j0.scale(2.0);
  } else { // Already pre-computed
    v_j0 = stored_v_coul;
  }
  molresponse::end_timer(world, "Coulomb Potential J[ground_density]");

  if (xcf.hf_exchange_coefficient() != 1.0) {
    v_xc = xc.make_xc_potential();
  } else {
    // make a zero function
    v_xc = Function<double, 3>(
        FunctionFactory<double, 3>(world).fence(false).initial_level(1));
  }

  // Intermediaries

  molresponse::start_timer(world);

  auto k = [&phi0_copy](vector_real_function_3d xi) {
    return K(phi0_copy, phi0_copy, xi);
  };
  // If including any exact HF exchange
  if (xcf.hf_exchange_coefficient() != 0.0) {
    std::transform(Chi_copy.X.begin(), Chi_copy.X.end(), K0.X.begin(), k);
    if (compute_Y) {
      std::transform(Chi_copy.Y.begin(), Chi_copy.Y.end(), K0.Y.begin(), k);
    }
  }
  if (r_params.print_level() >= 10) {
    print("inner <X|K0|X>");
    print(inner(Chi_copy, K0));
  }
  molresponse::end_timer(world, "K[ground_density]");
  // Vnuc+V0+VXC
  molresponse::start_timer(world);
  real_function_3d v0 = v_j0 + v_nuc + v_xc;

  V0.X = v0 * X.X;
  V0.X += (-1 * K0.X * xcf.hf_exchange_coefficient());

  if (compute_Y) {
    V0.Y = v0 * X.Y;
    V0.Y += (-1 * K0.Y * xcf.hf_exchange_coefficient());
  }
  if (r_params.print_level() >= 3) {
    print("inner <X|V0|X>");
    print(inner(Chi_copy, V0));
  }
  molresponse::end_timer(world, "V0X");

  // Basic output

  // Done
  return V0;
}

// Returns the ground state fock operator applied to functions f
X_space ResponseBase::compute_F0X(World &world, const X_space &X, const XCOperator<double, 3> &xc,
                                  bool compute_Y) const {
  // Debugging output

  molresponse::start_timer(world);
  size_t m = X.num_states();
  size_t n = X.num_orbitals();

  X_space chi_copy = X.copy();
  // chi_copy.truncate();
  X_space F0X = X_space(world, m, n);
  X_space T0X = X_space(world, m, n);
  T0X.X = T(world, chi_copy.X);
  if (compute_Y) {
    T0X.Y = T(world, chi_copy.Y);
  }
  if (r_params.print_level() >= 20) {
    print("_________________compute F0X _______________________");
    print("inner <X|T0|X>");
    print(inner(chi_copy, T0X));
  }

  X_space V0X = compute_V0X(world, chi_copy, xc, compute_Y);
  if (r_params.print_level() >= 20) {
    print("_________________compute F0X _______________________");
    print("inner <X|V0|X>");
    print(inner(chi_copy, V0X));
  }

  F0X = T0X + V0X;

  if (r_params.print_level() >= 20) {
    print("_________________compute F0X _______________________");
    print("inner <X|F0|X>");
    print(inner(chi_copy, F0X));
  }

  molresponse::end_timer(world, "F0X:");
  // Done
  return F0X;
}
X_space ResponseBase::compute_residual(World &world, X_space &old_Chi,
                                       X_space &temp,
                                       Tensor<double> &bsh_residualsX,
                                       Tensor<double> &bsh_residualsY,
                                       std::string calc_type) {
  size_t m = old_Chi.X.size();
  size_t n = old_Chi.X.size_orbitals();
  molresponse::start_timer(world);
  //	compute residual
  X_space res(world, m, n);
  res.X = old_Chi.X - temp.X;
  //
  if (calc_type.compare("full") == 0) {
    res.Y = old_Chi.Y - temp.Y;
  } else if (calc_type.compare("static") == 0) {
    res.Y = res.X.copy();
  } else {
  }
  //*************************
  Tensor<double> errX(m);
  Tensor<double> errY(m);
  // rmsX and maxvalX for each m response states
  std::vector<double> rmsX(m), maxvalX(m);
  std::vector<std::vector<double>> rnormsX;
  std::vector<std::vector<double>> rnormsY;
  // find the norms of each of the residual response vectors
  for (size_t i = 0; i < m; i++) {
    // the 2norms of each of the orbitals in response vector
    rnormsX.push_back(norm2s(world, res.X[i]));
    if (world.rank() == 0 and (r_params.print_level() > 1))
      print("residuals X: state ", i, " : ", rnormsX[i]);
    // maxabsval = std::max<double>(maxabsval, std::abs(v[i]));
    // maxvalX= largest abs(v[i])
    vector_stats(rnormsX[i], rmsX[i], maxvalX[i]);
    // errX[i] is the largest residual orbital value
    errX[i] = maxvalX[i];
    if (world.rank() == 0 and (r_params.print_level() > 1))
      print("BSH residual: rms", rmsX[i], "   max", maxvalX[i]);
  }
  if (calc_type.compare("full") == 0) {
    std::vector<double> rmsY(m), maxvalY(m);
    for (size_t i = 0; i < m; i++) {
      rnormsY.push_back(norm2s(world, res.Y[i]));
      if (world.rank() == 0 and (r_params.print_level() > 1))
        print("residuals Y: state ", i, " : ", rnormsY[i]);
      vector_stats(rnormsY[i], rmsY[i], maxvalY[i]);
      errY[i] = maxvalY[i];
      if (world.rank() == 0 and (r_params.print_level() > 1))
        print("BSH residual: rms", rmsY[i], "   max", maxvalY[i]);
    }
  } else if (calc_type.compare("static") == 0) {
    // copy by value?
    errY = errX;
  } else { // tda
    errY.fill(0);
  }
  molresponse::end_timer(world, "BSH residual");

  if (r_params.print_level() >= 1) {
    print("res.X norms in iteration after compute_residual function: ");
    print(res.X.norm2());

    print("res.Y norms in iteration after compute_residual function: ");
    print(res.Y.norm2());
  }
  // the max error in residual
  bsh_residualsX = errX;
  bsh_residualsY = errY;
  // Apply shifts and rhs
  // Next calculate 2-norm of these vectors of differences
  return res;
}
X_space ResponseBase::kain_x_space_update(World &world, const X_space &temp,
                                          const X_space &res,
                                          NonLinearXsolver &kain_x_space,
                                          std::vector<X_vector> &Xvector,
                                          std::vector<X_vector> &Xresidual) {
  size_t m = temp.num_states();
  size_t n = temp.num_orbitals();
  X_space kain_update(world, m, n);
  molresponse::start_timer(world);
  for (size_t b = 0; b < m; b++) {
    Xvector[b].X[0] = copy(world, temp.X[b]);
    Xvector[b].Y[0] = copy(world, temp.Y[b]);
    Xresidual[b].X[0] = copy(world, res.X[b]);
    Xresidual[b].Y[0] = copy(world, res.Y[b]);
  }

  for (size_t b = 0; b < m; b++) {
    // passing xvectors
    X_vector kain_X = kain_x_space[b].update(
        Xvector[b], Xresidual[b], FunctionDefaults<3>::get_thresh(), 3.0);
    // deep copy of vector of functions
    kain_update.X[b] = copy(world, kain_X.X[0]);
    kain_update.Y[b] = copy(world, kain_X.Y[0]);
  }
  molresponse::end_timer(world, " KAIN update:");
  return kain_update;
}

void ResponseBase::x_space_step_restriction(World &world, X_space &old_Chi,
                                            X_space &temp, bool restrict_y,
                                            Tensor<double> &maxrotn) {
  size_t m = old_Chi.num_states();
  molresponse::start_timer(world);
  print(maxrotn);

  for (size_t b = 0; b < m; b++) {
    if (true) {
      // do_step_restriction(world, old_Chi.X[b], temp.X[b], old_Chi.Y[b],
      // temp.Y[b], "x and y_response"); if the norm(new-old) > max
      do_step_restriction(world, old_Chi.X[b], temp.X[b], "x_response",
                          maxrotn[b]);

      do_step_restriction(world, old_Chi.Y[b], temp.Y[b], "y_response",
                          maxrotn[b]);
      // do_step_restriction(world, old_Chi.X[b], temp.X[b], "x_response");
      // do_step_restriction(world, old_Chi.Y[b], temp.Y[b], "y_response");
    } else {
      do_step_restriction(world, old_Chi.X[b], temp.X[b], "x_response",
                          maxrotn[b]);
    }
  }

  molresponse::end_timer(world, " Step Restriction:");
}

// compute rms and maxabsval of vector of doubles
void ResponseBase::vector_stats(const std::vector<double> &v, double &rms,
                                double &maxabsval) const {
  rms = 0.0;
  maxabsval = v[0];
  for (size_t i = 0; i < v.size(); ++i) {
    rms += v[i] * v[i];
    maxabsval = std::max<double>(maxabsval, std::abs(v[i]));
  }
  rms = sqrt(rms / v.size());
}

void ResponseBase::vector_stats_new(const Tensor<double> v, double &rms,
                                    double &maxabsval) const {
  rms = 0.0;
  for (size_t i = 0; i < v.size(); ++i) {
    rms += v[i] * v[i];
  }
  rms = sqrt(rms / v.size());
  maxabsval = v.max();
}

double ResponseBase::do_step_restriction(World &world, const vecfuncT &x,
                                         vecfuncT &x_new,
                                         std::string spin) const {
  std::vector<double> anorm = norm2s(world, sub(world, x, x_new));
  size_t nres = 0;
  for (unsigned int i = 0; i < x.size(); ++i) {
    print("anorm ", i, " : ", anorm[i]);
    if (anorm[i] > r_params.maxrotn()) {
      double s = r_params.maxrotn() / anorm[i];
      ++nres;
      if (world.rank() == 0) {
        if (nres == 1 and (r_params.print_level() > 1))
          printf("  restricting step for %s orbitals:", spin.c_str());
        printf(" %d", i);
      }
      x_new[i].gaxpy(s, x[i], 1.0 - s, false);
      x_new[i].truncate();
    }
  }
  if (nres > 0 && world.rank() == 0 and (r_params.print_level() > 1))
    printf("\n");

  world.gop.fence();
  double rms, maxval;
  vector_stats(anorm, rms, maxval);
  if (world.rank() == 0 and (r_params.print_level() > 1))
    print("Norm of vector changes", spin, ": rms", rms, "   max", maxval);
  return maxval;
}

double ResponseBase::do_step_restriction(World &world, const vecfuncT &x,
                                         vecfuncT &x_new, std::string spin,
                                         double maxrotn) const {
  Tensor<double> anorm = norm2s_T(world, sub(world, x, x_new));
  print("ANORM", anorm);
  print("maxrotn: ", maxrotn);
  for (unsigned int i = 0; i < x_new.size(); ++i) {
    print("anorm ", i, " : ", anorm[i]);
    if (anorm[i] > maxrotn) {
      double s = maxrotn / anorm[i];
      /*
if (world.rank() == 0) {
if (nres == 1 and (r_params.print_level() > 1)) printf("  restricting step for
%s orbitals:", spin.c_str()); printf(" %d", i);
}
*/
      x_new[i].gaxpy(s, x[i], 1.0 - s, false);
      // x_new[i].truncate();
    }
  }
  world.gop.fence();
  double rms, maxval;
  vector_stats_new(anorm, rms, maxval);
  if (world.rank() == 0 and (r_params.print_level() > 1))
    print("Norm of vector changes", spin, ": rms", rms, "   max", maxval);
  return maxval;
}

double ResponseBase::do_step_restriction(World &world, const vecfuncT &x,
                                         const vecfuncT &y, vecfuncT &x_new,
                                         vecfuncT &y_new,
                                         std::string spin) const {
  // sub(world, x, x_new)
  vecfuncT x_diff = sub(world, x, x_new);
  vecfuncT y_diff = sub(world, y, y_new);

  // sub(world, x, x_new)
  Tensor<double> anorm_x = norm2s_T(world, x_diff);
  Tensor<double> anorm_y = norm2s_T(world, y_diff);
  Tensor<double> anorm(x.size());
  for (unsigned int i = 0; i < x.size(); ++i) {
    anorm[i] = std::sqrt(anorm_x[i] * anorm_x[i] + anorm_y[i] * anorm_y[i]);
  }

  for (unsigned int i = 0; i < x.size(); ++i) {
    if (anorm[i] > r_params.maxrotn()) {
      double s = r_params.maxrotn() / anorm[i];
      size_t nres = 0;
      if (world.rank() == 0) {
        if (nres == 1 and (r_params.print_level() > 1))
          printf("  restricting step for %s orbitals:", spin.c_str());
        printf(" %d", i);
      }
      x_new[i].gaxpy(s, x[i], 1.0 - s, false);
      y_new[i].gaxpy(s, y[i], 1.0 - s, false);
    }
  }

  world.gop.fence();
  double rms, maxval;
  vector_stats_new(anorm, rms, maxval);
  if (world.rank() == 0 and (r_params.print_level() > 1))
    print("Norm of vector changes", spin, ": rms", rms, "   max", maxval);
  return maxval;
}
void ResponseBase::PlotGroundandResponseOrbitals(
    World &world, size_t iteration, response_space &x_response,
    response_space &y_response, ResponseParameters const &r_params,
    GroundStateCalculation const &g_params) {
  std::filesystem::create_directories("plots/densities");
  std::filesystem::create_directory("plots/orbitals");

  // TESTING
  // get transition density
  // num orbitals
  size_t n = x_response[0].size();
  size_t m = x_response.size();

  real_function_3d rho0 = dot(world, ground_orbitals, ground_orbitals);
  std::vector<real_function_3d> rho1 =
      transition_density(world, ground_orbitals, x_response, y_response);
  std::string dir("xyz");
  // for plotname size
  size_t buffSize = 500;
  char plotname[buffSize];
  double Lp = std::min(r_params.L(), 24.0);
  // Doing line plots along each axis
  for (int d = 0; d < 3; d++) {
    // print ground_state
    plotCoords plt(d, Lp);
    // plot ground density
    if (iteration == 1) {
      snprintf(plotname, buffSize, "plots/densities/rho0_%c_0.plt", dir[d]);
      plot_line(plotname, 5001, plt.lo, plt.hi, rho0);
    }
    for (int i = 0; i < static_cast<int>(n); i++) {
      // print ground_state
      // plot gound_orbitals
      snprintf(plotname, buffSize, "plots/orbitals/phi0_%c_0_%d.plt", dir[d],
               static_cast<int>(i));
      plot_line(plotname, 5001, plt.lo, plt.hi, ground_orbitals[i]);
    }

    for (int b = 0; b < static_cast<int>(m); b++) {
      // plot rho1 direction d state b
      snprintf(plotname, buffSize, "plots/densities/rho1_%c_%d.plt", dir[d],
               static_cast<int>(b));
      plot_line(plotname, 5001, plt.lo, plt.hi, rho1[b]);

      for (int i = 0; i < static_cast<int>(n); i++) {
        // print ground_state
        // plot x function  x_dir_b_i__k_iter
        snprintf(plotname, buffSize, "plots/orbitals/phix_%c_%d_%d.plt", dir[d],
                 static_cast<int>(b), static_cast<int>(i));
        plot_line(plotname, 5001, plt.lo, plt.hi, x_response[b][i]);

        // plot y functione  y_dir_b_i__k_iter
        snprintf(plotname, buffSize, "plots/orbitals/phiy_%c_%d_%d.plt", dir[d],
                 static_cast<int>(b), static_cast<int>(i));
        plot_line(plotname, 5001, plt.lo, plt.hi, y_response[b][i]);
      }
    }
  }
  world.gop.fence();

  // END TESTING
}

void PlotGroundDensityVTK(World &world, const ResponseBase &calc) {

  auto [ground_calc, molecule, r_params] = calc.get_parameter();
  auto ground_orbitals = calc.get_orbitals();

  if (r_params.plot_initial()) {
    if (world.rank() == 0)
      print("\n   Plotting ground state densities.\n");
    if (r_params.plot_l() > 0.0)
      do_vtk_plots(world, r_params.plot_pts(), r_params.plot_l(), 0,
                   r_params.num_orbitals(), molecule,
                   square(world, ground_orbitals), "ground");
    else
      do_vtk_plots(world, r_params.plot_pts(), r_params.L() / 2.0, 0,
                   r_params.num_orbitals(), molecule,
                   square(world, ground_orbitals), "ground");
  }
}

/// Push back empty json onto "protocol_data" field
/// Writes protocol to that new field
/// Sets up the iteration data json with an empty {}
/// TODO This can work for any madness iteration therefore maybe I should move
/// it \param j \param proto
void protocol_to_json(json &j, const double proto) {
  j["protocol_data"].push_back({});
  auto proto_index = j["protocol_data"].size() - 1;
  j["protocol_data"][proto_index]["proto"] = proto;
  j["protocol_data"][proto_index]["iter_data"] = {};
}

void ResponseBase::solve(World &world) {
  // Get start time
  molresponse::start_timer(world);

  // Plotting input orbitals
  if (r_params.plot_initial()) {
    PlotGroundDensityVTK(world, *this);
  }
  // TODO Why would I plot the ground state density here if the protocol or k is
  // not set yet? Warm and fuzzy
  /*
  //if (world.rank() == 0) {
      print("\n\n     Excited State Calculation");
      print("   ------------------------");
  }
   */
  // Here print the relevant parameters
  // Ready to iterate!
  const auto protocol = r_params.protocol();
  if (world.rank() == 0) {
    print("Response State Calculation for the following protocols");
    print("Protocol: ", protocol);
  }
  bool first_protocol = true;
  for (const auto &iter_thresh : protocol) {
    // We set the protocol and function defaults here for the given threshold of
    // protocol
    set_protocol(world, iter_thresh);
    check_k(world, iter_thresh, FunctionDefaults<3>::get_k());
    protocol_to_json(j_molresponse, iter_thresh);
    if (first_protocol) {
      if (r_params.restart()) {
        if (world.rank() == 0) {
          print("   Restarting from file:", r_params.restart_file());
        }
        load(world, r_params.restart_file());
        check_k(world, iter_thresh, FunctionDefaults<3>::get_k());
        first_protocol = false;
      } else {
        this->initialize(world);
      }
      first_protocol = false;
    }
    // Now actually ready to iterate...
    this->iterate(world);
  }
  // At this point we should know if calc converged maybe add a flag to response.json which states if it has
  converged_to_json(j_molresponse);



  // Plot the response function if desired
}

void check_k(World &world, X_space &Chi, double thresh, int k) {
  if (-1 != Chi.X.size()) {
    if (FunctionDefaults<3>::get_k() != Chi.X[0].at(0).k()) {
      // Project all x components into correct k

      for (auto &xi : Chi.X) {
        reconstruct(world, xi);
        for (auto &xij : xi) {
          xij = project(xij, FunctionDefaults<3>::get_k(), thresh, false);
        }
        world.gop.fence();
      }
      for (auto &yi : Chi.Y) {
        reconstruct(world, yi);
        for (auto &yij : yi) {
          yij = project(yij, FunctionDefaults<3>::get_k(), thresh, false);
        }
        world.gop.fence();
      }
      Chi.truncate();
    }
  }
}

///  @brief Adds random noise to functions in response space
///
///
///
/// \param world
/// \param f
/// \param magnitude
/// \return
response_space add_randomness(World &world, const response_space &f,
                              double magnitude) {
  // Copy input functions
  response_space f_copy = f.copy();

  // Lambda function to add in noise
  auto noise = [](const Key<3> &key, Tensor<double> &x) mutable {
    Tensor<double> y(x.size());
    y.fillrandom();
    // y.scale(magnitude);
    y.scale(1e3);
    x = x + y;
    // x(0,0,0) += y(0,0,0)-0.5;
  };
  // TODO
  // Go through each function in f_copy and add in random noise

  for (auto &fi : f_copy) {
    for (auto &fij : fi) {
      fij.unaryop(noise);
    }
  }
  // Done
  return f_copy;
}
/***
 * @brief normalize a single response space
 *
 * \note a single response space consists of only x or y states
 *
 *
 * @param world
 * @param f
 */
void normalize(World &world, response_space &f) {
  // Run over rows
  for (auto &fi : f) {
    double norm = inner(fi, fi);
    norm = sqrt(norm);
    // Doing this to deal with zero functions.
    // Maybe not smrt.
    if (norm == 0)
      continue;
    // And scale
    fi = fi * (1.0 / norm);
  }
}

void normalize(World &world, X_space &Chi) {
  // Run over rows

  for (size_t i = 0; i < Chi.num_states(); i++) {
    // Get the normalization constant
    // (Sum included inside inner)
    double norm_x = inner(Chi.X[i], Chi.X[i]);
    double norm_y = inner(Chi.Y[i], Chi.Y[i]);
    double norm = sqrt(norm_x - norm_y);
    // Doing this to deal with zero functions.
    // Maybe not smrt.
    if (norm == 0)
      continue;
    Chi.X[i] = Chi.X[i] * (1.0 / norm);
    Chi.Y[i] = Chi.Y[i] * (1.0 / norm);
  }
}
std::map<std::vector<int>, real_function_3d> solid_harmonics(World &world,
                                                             int n) {
  // Container to return
  std::map<std::vector<int>, real_function_3d> result;

  // Create the basic x, y, z, constant and zero
  real_function_3d x = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{1, 0, 0})));
  real_function_3d y = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 1, 0})));
  real_function_3d z = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 1})));
  real_function_3d c = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 0})));
  real_function_3d zero = real_factory_3d(world);

  // Add in first few, since they're simple
  // Assuming n >= 1
  result[std::vector<int>{0, 0}] = copy(c);
  result[std::vector<int>{0, -1}] = zero;
  result[std::vector<int>{0, 1}] = zero;
  result[std::vector<int>{-1, 0}] = zero;

  // Generate the solid harmonics recursively from here
  for (int l = 0; l < n; l++) {
    // Calculate ends of this row first
    result[std::vector<int>{l + 1, l + 1}] =
        sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
        (x * result[std::vector<int>{l, l}] -
         (1 - kronecker(l, 0) * y * result[std::vector<int>{l, -l}]));
    result[std::vector<int>{l + 1, -l - 1}] =
        sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
        (y * result[std::vector<int>{l, l}] +
         (1 - kronecker(l, 0) * x * result[std::vector<int>{l, -l}]));

    // Formula below calls for some functions that don't exist.
    // Need zeroes where that would occur
    result[std::vector<int>{l + 1, l + 2}] = zero;
    result[std::vector<int>{l + 1, -l - 2}] = zero;

    // Run over quantum number m
    for (int m = -l; m < l + 1; m++) {
      // Calculate remaining terms
      result[std::vector<int>{l + 1, m}] =
          1.0 / std::sqrt((l + m + 1) * (l - m + 1)) *
          ((2 * l + 1) * z * result[std::vector<int>{l, m}] -
           sqrt((l + m) * (l - m)) * (x * x + y * y + z * z) *
               result[std::vector<int>{l - 1, m}]);
    }
  }

  // Get rid of any zero functions we added
  for (auto it = result.begin(); it != result.end();) {
    if (it->second.norm2() == 0)
      it = result.erase(it);
    else
      ++it;
  }

  // Also get rid of the constant
  result.erase(std::vector<int>{0, 0});

  // Done
  return result;
}
vector_real_function_3d
transition_density(World &world, const vector_real_function_3d &orbitals,
                   const response_space &x, const response_space &y) {
  // Get sizes
  // Check sizes and then run the algorithm
  size_t m = x.size();
  auto xx = x.copy();
  auto yy = y.copy();
  auto phi0 = copy(world, orbitals);

  xx.truncate_rf();
  yy.truncate_rf();
  truncate(world, phi0);

  // Return container
  std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);

  std::transform(x.begin(), x.end(), y.begin(), densities.begin(),
                 [&world, &phi0](auto x_alpha, auto y_alpha) {
                   auto dx = dot(world, x_alpha, phi0);
                   world.gop.fence();
                   auto dy = dot(world, phi0, y_alpha);
                   return dx + dy;
                 });
  truncate(world, densities);
  world.gop.fence();
  // Done!
  return densities;
}
/***
 * @brief returns orbitals and response functions with new p map
 *
 * @param world
 * @param psi0
 * @param X
 * @param load_balance
 * @return
 */
gamma_orbitals ResponseBase::orbital_load_balance(World &world,
                                                  const gamma_orbitals &input,
                                                  const double load_balance) {

  auto X = std::get<0>(input);
  auto psi0 = std::get<1>(input);
  auto vf = std::get<2>(input);

  size_t m = X.num_states();
  size_t n = X.num_orbitals();

  if (world.size() > 1) {

    molresponse::start_timer(world);
    LoadBalanceDeux<3> lb(world);

    for (unsigned int i = 0; i < m; ++i) {
      lb.add_tree(psi0[i], lbcost<double, 3>(1.0, 8.0), false);
      lb.add_tree(vf[i], lbcost<double, 3>(1.0, 8.0), false);
      for (unsigned int j = 0; j < n; ++j) {
        // add a tree for orbitals
        lb.add_tree(X.X[i][j], lbcost<double, 3>(1.0, 8.0), false);
        lb.add_tree(X.Y[i][j], lbcost<double, 3>(1.0, 8.0), false);
      }
    }
    world.gop.fence();
    // newpamap is the new pmap just based on the orbitals
    std::shared_ptr<WorldDCPmapInterface<Key<3>>> new_process_map =
        lb.load_balance(load_balance);
    molresponse::end_timer(world, "Gamma compute load_balance_chi");
    // default process map
    // We set the new_process_map
    molresponse::start_timer(world);
    FunctionDefaults<3>::set_pmap(new_process_map); // set default to be new

    world.gop.fence();
    // copy orbitals using new pmap
    auto X_copy = X.copy(new_process_map, false);
    world.gop.fence(); // then fence

    auto psi0_copy = copy(world, psi0, new_process_map, false);
    auto vf_copy = copy(world, vf, new_process_map, false);
    world.gop.fence(); // then fence
    molresponse::end_timer(world, "Gamma redist");
    return {X_copy, psi0_copy, vf_copy};
  } else {
    // return a copy with the same process map since we only have one world
    return {X.copy(), copy(world, psi0), copy(world, vf)};
  }
}
void ResponseBase::analyze_vectors(World &world, const vecfuncT &x,
                                   const std::string &response_state) {
  molresponse::start_timer(world);
  AtomicBasisSet sto3g("sto-3g");
  vecfuncT ao = project_ao_basis(world, sto3g);

  tensorT C = matrix_inner(world, ao, x);
  int nmo1 = x.size();
  tensorT rsq, dip(3, nmo1);
  {
    functionT frsq = factoryT(world).f(rsquared).initial_level(4);
    // <x r^2 x>
    //<x[i] | r^2 | x[i]>
    rsq = inner(world, x, mul_sparse(world, frsq, x, vtol));
    for (int axis = 0; axis < 3; ++axis) {
      // x y z
      functionT fdip = factoryT(world)
                           .functor(functorT(new madness::DipoleFunctor(axis)))
                           .initial_level(4);
      dip(axis, _) = inner(world, x, mul_sparse(world, fdip, x, vtol));
      //<x r^2 x> - <x|x|x>^2-<x|y|x>^2-<x|z|x>^2
      for (int i = 0; i < nmo1; ++i)
        rsq(i) -= dip(axis, i) * dip(axis, i);
    }
  }
  molresponse::end_timer(world, "Analyze vectors");

  long nmo = x.size();
  size_t ncoeff = 0;
  for (long i = 0; i < nmo; ++i) {
    size_t ncoeffi = x[i].size();
    ncoeff += ncoeffi;
    if (world.rank() == 0 and (r_params.print_level() > 1)) {
      print(response_state + " orbital : ", i);

      printf("ncoeff=%.2e:", (double)ncoeffi);

      printf("center=(%.2f,%.2f,%.2f) : radius=%.2f\n", dip(0, i), dip(1, i),
             dip(2, i), sqrt(rsq(i)));
      sto3g.print_anal(molecule, C(i, _));
      printf("total number of coefficients = %.8e\n\n", double(ncoeff));
    }
  }
}
vecfuncT ResponseBase::project_ao_basis_only(World &world,
                                             const AtomicBasisSet &aobasis,
                                             const Molecule &molecule) {
  vecfuncT ao = vecfuncT(aobasis.nbf(molecule));
  for (int i = 0; i < aobasis.nbf(molecule); ++i) {
    functorT aofunc(
        new AtomicBasisFunctor(aobasis.get_atomic_basis_function(molecule, i)));
    ao[i] = factoryT(world)
                .functor(aofunc)
                .truncate_on_project()
                .nofence()
                .truncate_mode(1);
  }
  world.gop.fence();
  truncate(world, ao);
  madness::normalize(world, ao);
  return ao;
}
vecfuncT ResponseBase::project_ao_basis(World &world,
                                        const AtomicBasisSet &aobasis) {
  // Make at_to_bf, at_nbf ... map from atom to first bf on atom, and nbf/atom
  std::vector<int> at_to_bf, at_nbf;
  aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);

  return project_ao_basis_only(world, aobasis, molecule);
}
void ResponseBase::output_json() const {
  std::ofstream ofs("response_base.json");
  ofs << j_molresponse;
}
void ResponseBase::converged_to_json(json &j) {
    j["converged"]=converged;
}


vector_real_function_3d
transition_densityTDA(World &world, const vector_real_function_3d &orbitals,
                      const response_space &x) {

  // Get sizes
  size_t m = x.size();

  auto xx = x;
  auto phi0 = copy(world, orbitals);

  xx.truncate_rf();
  truncate(world, phi0);

  std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);

  // dot xi with phi0
  auto f = [&world, &phi0](auto xi) { return dot(world, xi, phi0); };

  // for each vector is response space x dot and
  std::transform(x.begin(), x.end(), densities.begin(), f);

  truncate(world, densities);
  world.gop.fence();
  return densities;
}

// Transforms the given matrix of functions according to the give
// transformation matrix. Used to update orbitals / potential
response_space transform(World &world, const response_space &f,
                         const Tensor<double> &U) {
  // Return container
  response_space result;

  // Go element by element
  for (unsigned int i = 0; i < f.size(); i++) {
    // Temp for the result of one row
    std::vector<real_function_3d> temp =
        zero_functions_compressed<double, 3>(world, f[0].size());

    for (unsigned int j = 0; j < f.size(); j++) {
      gaxpy(world, 1.0, temp, U(j, i), f[j]);
    }

    // Add to temp to result
    result.push_back(temp);
  }

  result.truncate_rf();

  // Done
  return result;
}
Tensor<double> expectation(World &world, const response_space &A,
                           const response_space &B) {
  // Get sizes
  MADNESS_ASSERT(A.size() > 0);
  MADNESS_ASSERT(A.size() == B.size());
  MADNESS_ASSERT(A[0].size() > 0);
  MADNESS_ASSERT(A[0].size() == B[0].size());

  size_t dim_1 = A.size();
  size_t dim_2 = A[0].size();
  // Need to take transpose of each input ResponseFunction
  response_space A_t(world, dim_2, dim_1);
  response_space B_t(world, dim_2, dim_1);
  for (size_t i = 0; i < dim_1; i++) {
    for (size_t j = 0; j < dim_2; j++) {
      A_t[j][i] = A[i][j];
      B_t[j][i] = B[i][j];
    }
  }
  // Container for result
  Tensor<double> result(dim_1, dim_1);
  /**
   * @brief
   * [x1 x2 x3]T[x1 x2 x3]
   *
   */
  // Run over dimension two
  // each vector in orbital has dim_1 response functoins associated
  for (size_t p = 0; p < dim_2; p++) {
    result += matrix_inner(world, A_t[p], B_t[p]);
  }

  // Done
  return result;
}
void print_norms(World &world, const response_space &f) {

  print(f[0].size());
  Tensor<double> norms(f.size() * f[0].size());
  // Calc the norms
  long i = 0;
  for (const auto &fi : f) {

    for (const auto &fij : fi) {
      norms(i++) = fij.norm2();
    }
  }
  norms = norms.reshape(f.size(), f[0].size());
  // Print em in a smart way
  if (world.rank() == 0)
    print(norms);
}
response_space select_functions(World &world, response_space f,
                                Tensor<double> &energies, size_t k,
                                size_t print_level) {
  // Container for result
  response_space answer;

  // Debugging output
  if (print_level >= 1) {
    if (world.rank() == 0)
      print("\n   Selecting the", k, "lowest excitation energy components.\n");
  }

  // Get rid of extra functions and save
  // the first k
  while (f.size() > k)
    f.pop_back();
  answer = f;
  answer.truncate_rf();

  // Get rid of extra energies and save
  // the first k
  energies = energies(Slice(0, k - 1));

  // Basic output
  if (print_level >= 1) {
    if (world.rank() == 0)
      print("   The selected components have excitation energies:");
    if (world.rank() == 0)
      print(energies);
  }

  // Done
  return answer;
}
void sort(World &world, Tensor<double> &vals, response_space &f) {
  // Get relevant sizes
  size_t k = vals.size();

  // Copy everything...
  response_space f_copy(f);
  Tensor<double> vals_copy = copy(vals);
  Tensor<double> vals_copy2 = copy(vals);

  // Now sort vals_copy
  std::sort(vals_copy.ptr(), vals_copy.ptr() + vals_copy.size());

  // Now sort the rest of the things, using the sorted energy list
  // to find the correct indices
  for (size_t i = 0; i < k; i++) {
    // Find matching index in sorted vals_copy
    size_t j = 0;
    while (fabs(vals_copy(i) - vals_copy2(j)) > 1e-8 && j < k)
      j++;

    // Put corresponding function, difference function, value residual and
    // value in the correct place
    f[i] = f_copy[j];
    vals(i) = vals_copy(i);

    // Change the value of vals_copy2[j] to help deal with duplicates?
    vals_copy2(j) = 10000.0;
  }
}
// Sorts the given tensor of eigenvalues and
// response functions
void sort(World &world, Tensor<double> &vals, X_space &f) {
  // Get relevant sizes
  size_t k = vals.size();

  // Copy everything...
  X_space f_copy(f);
  Tensor<double> vals_copy = copy(vals);
  Tensor<double> vals_copy2 = copy(vals);

  // Now sort vals_copy
  std::sort(vals_copy.ptr(), vals_copy.ptr() + vals_copy.size());

  // Now sort the rest of the things, using the sorted energy list
  // to find the correct indices
  for (size_t i = 0; i < k; i++) {
    // Find matching index in sorted vals_copy
    size_t j = 0;
    while (fabs(vals_copy(i) - vals_copy2(j)) > 1e-8 && j < k)
      j++;

    // Put corresponding function, difference function, value residual and
    // value in the correct place
    f.X[i] = f_copy.X[j];
    f.Y[i] = f_copy.Y[j];

    vals(i) = vals_copy(i);

    // Change the value of vals_copy2[j] to help deal with duplicates?
    vals_copy2(j) = 10000.0;
  }
}
response_space gram_schmidt(World &world, const response_space &f) {
  // Sizes inferred
  size_t m = f.size();

  // Return container
  response_space result = f.copy();

  // Orthogonalize
  for (size_t j = 0; j < m; j++) {
    // Need to normalize the row
    double norm = norm2(world, result[j]);

    // Now scale each entry
    scale(world, result[j], 1.0 / norm);

    // Project out from the rest of the vectors
    for (size_t k = j + 1; k < m; k++) {
      // Temp function to hold the sum
      // of inner products
      // vmra.h function, line 627
      double temp = inner(result[j], result[k]);

      // Now subtract
      gaxpy(world, 1.0, result[k], -temp, result[j]);
    }
  }
  result.truncate_rf();

  // Done
  return result;
}
vector_real_function_3d make_xyz_functions(World &world) {
  // Container to return

  // Create the basic x, y, z, constant and zero
  real_function_3d x = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{1, 0, 0})));
  real_function_3d y = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 1, 0})));
  real_function_3d z = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 1})));

  std::vector<real_function_3d> funcs = {x, y, z};
  return funcs;
}
void gram_schmidt(World &world, response_space &f, response_space &g) {
  // Sizes inferred
  size_t m = f.size();

  // Orthogonalize
  for (size_t j = 0; j < m; j++) {
    // Need to normalize the row
    double norm = inner(f[j], f[j]) - inner(g[j], g[j]);

    // Now scale each entry
    scale(world, f[j], 1.0 / sqrt(norm));
    scale(world, g[j], 1.0 / sqrt(norm));

    // Project out from the rest of the vectors
    for (size_t k = j + 1; k < m; k++) {
      // Temp function to hold the sum
      // of inner products
      // vmra.h function, line 627
      double temp = inner(f[j], f[k]) - inner(g[j], g[k]);

      // Now subtract
      gaxpy(world, 1.0, f[k], -temp, f[j]);
      gaxpy(world, 1.0, g[k], -temp, g[j]);
    }
  }

  f.truncate_rf();
  g.truncate_rf();
}
