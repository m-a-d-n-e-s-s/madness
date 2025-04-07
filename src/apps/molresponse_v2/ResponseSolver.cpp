#include "ResponseSolver.hpp"
#include "../molresponse/response_macrotask.hpp"
#include "ResponseManager.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseVector.hpp"
#include "functypedefs.h"
#include "projector.h"

vector_real_function_3d StaticRestrictedSolver::ComputeRSH(
    World &world, const GroundStateData &gs, const ResponseVector &vecs,
    const vector_real_function_3d &vp, const std::vector<poperatorT> &bsh_x,
    const ResponseManager &rm) {

  auto &x = std::get<StaticRestrictedResponse>(vecs).x_alpha;
  auto &all_x = std::get<StaticRestrictedResponse>(vecs).flat;

  auto num_orbitals = gs.orbitals.size();

  std::vector<int> state_index;
  std::vector<int> ii;

  int i = 0;
  int si = 0;
  for (int j = 0; j < num_orbitals; j++) {
    state_index.push_back(si);
    ii.push_back(i++);
  }

  ResponseComputeGroundExchange t0;
  MacroTask g0_task(world, t0);
  ResponseComputeGammaX tresponse;
  MacroTask gx_task(world, tresponse);

  auto c_xc = gs.xcf_.hf_exchange_coefficient();
  auto k0 = g0_task(ii, state_index, x, gs.orbitals, true);
  if (world.rank() == 0) {
    print("k0");
  }
  auto gx = gx_task(ii, state_index, x, gs.orbitals, true);
  if (world.rank() == 0) {
    print("gx");
  }
  auto v_local = gs.V_local * x;
  auto v0x = v_local - c_xc * k0;
  auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);

  auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
  if (world.rank() == 0) {
    print("thetax");
  }
  truncate(world, thetax);
  auto rsh = apply(world, bsh_x, thetax);
  QProjector<double, 3> projector(gs.orbitals);
  rsh =
      projector(rsh); // project out the ground state density from the response
  if (world.rank() == 0) {
    print("bsh_x");
  }
  return rsh;
}

vector_real_function_3d DynamicRestrictedSolver::ComputeRSH(
    World &world, const GroundStateData &gs, const ResponseVector &vecs,
    const vector_real_function_3d &vp, const std::vector<poperatorT> &bsh_x,
    const ResponseManager &rm) {

  auto &x = std::get<DynamicRestrictedResponse>(vecs).x_alpha;
  auto &y = std::get<DynamicRestrictedResponse>(vecs).y_alpha;
  auto &all_x = std::get<StaticRestrictedResponse>(vecs).flat;

  auto num_orbitals = gs.orbitals.size();

  std::vector<int> state_index;
  std::vector<int> ii;

  int i = 0;
  int si = 0;
  for (int j = 0; j < 2 * num_orbitals; j++) {
    state_index.push_back(si);
    ii.push_back(i++);
  }

  ResponseComputeGroundExchange t0;
  MacroTask g0_task(world, t0);
  ResponseComputeGammaX tresponse;
  MacroTask gx_task(world, tresponse);

  auto c_xc = gs.xcf_.hf_exchange_coefficient();
  auto k0 = g0_task(ii, state_index, all_x, gs.orbitals, true);
  auto gx = gx_task(ii, state_index, all_x, gs.orbitals, true);
  auto v_local = gs.V_local * all_x;
  auto v0x = v_local - c_xc * k0;
  auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);
  auto epsilony = transform(world, y, gs.Hamiltonian_no_diag, true);
  epsilonx.insert(epsilonx.end(), epsilony.begin(), epsilony.end());
  auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
  truncate(world, thetax);
  auto rsh = apply(world, bsh_x, thetax);
  truncate(world, thetax);
  QProjector<double, 3> projector(gs.orbitals);
  rsh =
      projector(rsh); // project out the ground state density from the response

  return rsh;
}

// (v0 + g0) * x + g[x] * phi0 + vp

bool StaticRestrictedSolver::iterate(World &world, const ResponseManager &rm,
                                     const GroundStateData &gs,
                                     const ResponseState &state,
                                     ResponseVector &response, size_t max_iter,
                                     double conv_thresh) {

  auto &rvec = std::get<StaticRestrictedResponse>(response);
  auto &x = std::get<StaticRestrictedResponse>(response).x_alpha;
  auto &all_x = std::get<StaticRestrictedResponse>(response).flat;

  auto thresh = FunctionDefaults<3>::get_thresh();
  const double dconv = std::max(FunctionDefaults<3>::get_thresh(), conv_thresh);
  auto density_target =
      dconv * static_cast<double>(std::max(size_t(5.0), gs.molecule.natom()));
  auto freq = state.current_frequency();
  const size_t n = x.size();

  vector_real_function_3d residual(n);
  vector_real_function_3d dx(n);

  // Set up RHS (perturbation vector)
  vector_real_function_3d Vp = state.perturbation_vector(world, gs);
  auto &phi0 = gs.orbitals;

  const auto &orbital_energies = gs.getEnergies();
  // set up the bsh operators
  auto bsh_x = std::vector<poperatorT>(n);
  double x_shifts = 0.0;
  if ((orbital_energies[long(n) - 1] + freq) >= 0.0) {
    // Calculate minimum shift needed such that \eps + \omega + shift < 0
    print("*** we are shifting just so you know!!!");
    x_shifts = -.05 - (freq + orbital_energies[long(n) - 1]);
  }
  bsh_x = ResponseSolverUtils::make_bsh_operators_response(
      world, x_shifts, freq, orbital_energies, rm.params().lo());

  response_solver solver(
      response_vector_allocator(world, static_cast<int>(all_x.size())), false);

  auto drho = StaticRestrictedSolver::compute_density(world, x, phi0);

  functionT drho_old;
  // Main iteration
  for (size_t iter = 0; iter < max_iter; ++iter) {
    drho_old = copy(drho);

    // 1. Update x: x_new = -2 * G(kx)[V0 x - εx + g(x) + Vp ]
    auto x_new =
        StaticRestrictedSolver::ComputeRSH(world, gs, rvec, Vp, bsh_x, rm);
    auto norm_xnew = norm2(world, x_new);
    if (world.rank() == 0)
      print("Norm of x_new:", norm_xnew);

    // 2. Form residual r = x_new - x
    auto rx = x_new - all_x;

    // 3. DIIS / KAIN extrapolation
    auto temp_vec = solver.update(all_x, rx);
    auto temp_norm = norm2(world, temp_vec);
    if (world.rank() == 0)
      print("Norm of temp_vec:", temp_norm);

    // 4. Apply step restriction
    auto res_norm = ResponseSolverUtils::do_step_restriction(world, all_x,
                                                             temp_vec, "a", 10);
    temp_norm = norm2(world, temp_vec);
    if (world.rank() == 0)
      print("Norm of temp_vec after restriction:", temp_norm);

    rvec.flat = copy(world, temp_vec);
    rvec.sync();
    // 5. Compute updated response density
    drho = StaticRestrictedSolver::compute_density(world, all_x, phi0);
    double drho_change = (drho - drho_old).norm2();

    rvec.flat = copy(world, temp_vec);
    rvec.sync();

    auto prop_i = -2.0 * dot(world, all_x, Vp).trace();

    if (world.rank() == 0)

      // 6. Convergence check
      if (world.rank() == 0) {
        print("Iter", iter, "  Δρ =", drho_change, "  xres =", res_norm,
              " prop_i =", prop_i);
      }

    if (drho_change < density_target && res_norm < conv_thresh) {
      if (world.rank() == 0)
        print("✓ Converged in", iter, "iterations.");

      // 🔄 Sync flat vector back into structured view
      auto &rvec = std::get<StaticRestrictedResponse>(response);
      rvec.sync();

      return true;
    }
  }

  if (world.rank() == 0)
    print("⚠️  Reached max iterations without convergence.");
  return false;
}

bool DynamicRestrictedSolver::iterate(World &world, const ResponseManager &rm,
                                      const GroundStateData &gs,
                                      const ResponseState &state,
                                      ResponseVector &response, size_t max_iter,
                                      double conv_thresh) {

  auto &rvec = std::get<DynamicRestrictedResponse>(response);
  auto &x = rvec.x_alpha;
  auto &y = rvec.y_alpha;
  auto &all_x = rvec.flat;

  const auto thresh = FunctionDefaults<3>::get_thresh();
  const double dconv = std::max(thresh, conv_thresh);
  const double density_target =
      dconv * static_cast<double>(std::max(size_t(5), gs.molecule.natom()));
  const double freq = state.current_frequency();
  const size_t n = x.size();

  // Set up perturbation
  auto Vp = state.perturbation_vector(world, gs);
  auto vp_copy = copy(world, Vp);
  Vp.insert(Vp.end(), vp_copy.begin(), vp_copy.end());

  const auto &phi0 = gs.orbitals;
  const auto &orbital_energies = gs.getEnergies();

  // Construct BSH operators for ω and -ω
  double x_shift = 0.0;
  if ((orbital_energies[long(n) - 1] + freq) >= 0.0) {
    print("*** Applying energy shift to BSH operator");
    x_shift = -.05 - (freq + orbital_energies[long(n) - 1]);
  }

  auto bsh_x = ResponseSolverUtils::make_bsh_operators_response(
      world, x_shift, freq, orbital_energies, rm.params().lo());
  auto bsh_y = ResponseSolverUtils::make_bsh_operators_response(
      world, 0.0, -freq, orbital_energies, rm.params().lo());

  std::vector<poperatorT> bsh_ops(2 * n);
  for (size_t i = 0; i < n; ++i) {
    bsh_ops[i] = bsh_x[i];
    bsh_ops[i + n] = bsh_y[i];
  }

  // Initial density and solver
  auto drho = DynamicRestrictedSolver::compute_density(world, all_x, phi0);
  auto drho_old = copy(drho);

  response_solver solver(
      response_vector_allocator(world, static_cast<int>(all_x.size())), false);

  // Main iteration loop
  for (size_t iter = 0; iter < max_iter; ++iter) {
    drho_old = copy(drho);

    // 1. Compute next x/y vector
    auto x_new = DynamicRestrictedSolver::ComputeRSH(world, gs, response, Vp,
                                                     bsh_ops, rm);

    auto norm_xnew = norm2(world, x_new);
    if (world.rank() == 0)
      print("Norm of x_new:", norm_xnew);

    // 2. Residual
    auto rx = x_new - all_x;

    // 3. DIIS / KAIN step
    auto temp_vec = solver.update(all_x, rx);
    auto temp_norm = norm2(world, temp_vec);
    if (world.rank() == 0)
      print("Norm of temp_vec:", temp_norm);

    // 4. Step restriction
    auto res_norm = ResponseSolverUtils::do_step_restriction(
        world, all_x, temp_vec, "a", rm.params().maxrotn());
    temp_norm = norm2(world, temp_vec);
    if (world.rank() == 0)
      print("Norm of temp_vec after restriction:", temp_norm);

    rvec.flat = copy(world, temp_vec);
    rvec.sync();

    // 5. Update and compute new response density
    drho = DynamicRestrictedSolver::compute_density(world, all_x, phi0);
    double drho_change = (drho - drho_old).norm2();

    // 6. Logging
    if (world.rank() == 0) {
      print("Iter", iter, "  Δρ =", drho_change, "  xres =", res_norm);
    }

    // 7. Convergence check
    if (drho_change < density_target && res_norm < conv_thresh) {
      if (world.rank() == 0)
        print("✓ Converged in", iter, "iterations.");

      // Sync flat → structured pointers
      rvec.sync();
      return true;
    }
  }

  if (world.rank() == 0)
    print("⚠️  Reached max iterations without convergence.");
  return false;
}
