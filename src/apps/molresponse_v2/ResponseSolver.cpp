#include "ResponseSolver.hpp"
#include "../molresponse/response_macrotask.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseDebugLoggerMacros.hpp"
#include "ResponseManager.hpp"
#include "ResponseSolverUtils.hpp"
#include "ResponseVector.hpp"
#include "functypedefs.h"
#include "projector.h"

// Compute density function
real_function_3d compute_density(World &world, const StaticRestrictedResponse &rvec,
                                 const vector_real_function_3d &phi0) {
  auto &x = rvec.x_alpha;
  auto xphi = mul(world, x, phi0, true);
  return 2.0 * sum(world, xphi, true);
}

real_function_3d compute_density(World &world, const DynamicRestrictedResponse &rvec,
                                 const vector_real_function_3d &phi0) {
  auto phi_phi = phi0;
  phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());
  auto &x = rvec.flat;

  auto xphi = mul(world, x, phi_phi, true);
  return 2.0 * sum(world, xphi, true);
}
real_function_3d compute_density(World &world, const StaticUnrestrictedResponse &rvec,
                                 const vector_real_function_3d &phi0) {
  auto phi_phi = phi0;
  phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());
  auto &x = rvec.flat;

  auto xphi = mul(world, x, phi_phi, true);
  return 2.0 * sum(world, xphi, true);
}

real_function_3d compute_density(World &world, const DynamicUnrestrictedResponse &rvec,
                                 const vector_real_function_3d &phi0) {
  auto phi_phi = phi0;
  phi_phi.insert(phi_phi.end(), phi0.begin(), phi0.end());
  auto &x = rvec.flat;

  auto xphi = mul(world, x, phi_phi, true);
  return 2.0 * sum(world, xphi, true);
}

std::vector<poperatorT> make_bsh_operators(World &world, const ResponseManager &rm, const double freq,
                                           const Tensor<double> &orbital_energies, const int n,
                                           ResponseDebugLogger &logger, const StaticRestrictedResponse & /* vecs */) {

  auto bsh_x = std::vector<poperatorT>(n);
  double x_shifts = 0.0;
  if ((orbital_energies[long(n) - 1] + freq) >= 0.0) {
    x_shifts = -.05 - (freq + orbital_energies[long(n) - 1]);
  }
  bsh_x = ResponseSolverUtils::make_bsh_operators_response(world, x_shifts, freq, orbital_energies, rm.params().lo());
  return bsh_x;
}

std::vector<poperatorT> make_bsh_operators(World &world, const ResponseManager &rm, const double freq,
                                           const Tensor<double> &orbital_energies, const int n,
                                           ResponseDebugLogger &logger, const DynamicRestrictedResponse & /* vecs */) {

  auto bsh_x = std::vector<poperatorT>(2 * n);
  double x_shifts = 0.0;
  if ((orbital_energies[long(n) - 1] + freq) >= 0.0) {
    x_shifts = -.05 - (freq + orbital_energies[long(n) - 1]);
  }
  bsh_x = ResponseSolverUtils::make_bsh_operators_response(world, x_shifts, freq, orbital_energies, rm.params().lo());
  auto bsh_y = ResponseSolverUtils::make_bsh_operators_response(world, 0.0, -freq, orbital_energies, rm.params().lo());

  bsh_x.insert(bsh_x.end(), bsh_y.begin(), bsh_y.end());
  return bsh_x;
}

vector_real_function_3d CoupledResponseEquations(World &world, const GroundStateData &gs,
                                                 const StaticRestrictedResponse &vecs,
                                                 const vector_real_function_3d &vp,
                                                 const std::vector<poperatorT> &bsh_x, const ResponseManager &rm,
                                                 ResponseDebugLogger &logger) {

  auto &x = vecs.x_alpha;
  auto &all_x = vecs.flat;

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
  vector_real_function_3d k0;
  vector_real_function_3d gx;

  DEBUG_TIMED_BLOCK(world, &logger, "g0_task", { k0 = g0_task(ii, state_index, all_x, gs.orbitals, true); });
  DEBUG_TIMED_BLOCK(world, &logger, "gx_task", { gx = gx_task(ii, state_index, all_x, gs.orbitals, true); });

  auto v_local = gs.V_local * x;
  auto v0x = v_local - c_xc * k0;
  auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);

  auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
  truncate(world, thetax);
  auto rsh = apply(world, bsh_x, thetax);
  rsh = gs.Qhat(rsh); // project out the ground state density from the response
  return rsh;
}

vector_real_function_3d CoupledResponseEquations(World &world, const GroundStateData &gs,
                                                 const DynamicRestrictedResponse &vecs,
                                                 const vector_real_function_3d &vp,
                                                 const std::vector<poperatorT> &bsh_x, const ResponseManager &rm,
                                                 ResponseDebugLogger &logger) {

  auto &x = vecs.x_alpha;
  auto &y = vecs.y_alpha;
  auto &all_x = vecs.flat;

  auto num_orbitals = gs.orbitals.size();

  std::vector<int> state_index;
  std::vector<int> ii;

  int i = 0;
  for (int j = 0; j < 2 * num_orbitals; j++) {
    state_index.push_back(0);
    ii.push_back(i++);
  }

  ResponseComputeGroundExchange t0;
  MacroTask g0_task(world, t0);
  ResponseComputeGammaX tresponse;
  MacroTask gx_task(world, tresponse);

  vector_real_function_3d k0;
  vector_real_function_3d gx;

  DEBUG_TIMED_BLOCK(world, &logger, "g0_task", { k0 = g0_task(ii, state_index, all_x, gs.orbitals, false); });
  DEBUG_TIMED_BLOCK(world, &logger, "gx_task", { gx = gx_task(ii, state_index, all_x, gs.orbitals, false); });

  auto c_xc = gs.xcf_.hf_exchange_coefficient();
  auto v_local = gs.V_local * all_x;
  auto v0x = v_local - c_xc * k0;
  auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);
  auto epsilony = transform(world, y, gs.Hamiltonian_no_diag, true);
  epsilonx.insert(epsilonx.end(), epsilony.begin(), epsilony.end());
  auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
  truncate(world, thetax);
  auto rsh = apply(world, bsh_x, thetax);
  truncate(world, thetax);
  rsh = gs.Qhat(rsh);
  return rsh;
}

// vector_real_function_3d ResponseSolverPolicy<StaticRestrictedResponse>::CoupledResponseEquations(
//     World &world, const GroundStateData &gs, const StaticRestrictedResponse &vecs, const vector_real_function_3d &vp,
//     const std::vector<poperatorT> &bsh_x, const ResponseManager &rm, ResponseDebugLogger &logger) {
//
//   auto &x = vecs.x_alpha;
//   auto &all_x = vecs.flat;
//
//   auto num_orbitals = gs.orbitals.size();
//
//   std::vector<int> state_index;
//   std::vector<int> ii;
//
//   int i = 0;
//   int si = 0;
//   for (int j = 0; j < num_orbitals; j++) {
//     state_index.push_back(si);
//     ii.push_back(i++);
//   }
//
//   ResponseComputeGroundExchange t0;
//   MacroTask g0_task(world, t0);
//   ResponseComputeGammaX tresponse;
//   MacroTask gx_task(world, tresponse);
//
//   auto c_xc = gs.xcf_.hf_exchange_coefficient();
//   vector_real_function_3d k0;
//   vector_real_function_3d gx;
//
//   DEBUG_TIMED_BLOCK(world, &logger, "g0_task", { k0 = g0_task(ii, state_index, all_x, gs.orbitals, true); });
//   DEBUG_TIMED_BLOCK(world, &logger, "gx_task", { gx = gx_task(ii, state_index, all_x, gs.orbitals, true); });
//
//   auto v_local = gs.V_local * x;
//   auto v0x = v_local - c_xc * k0;
//   auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);
//
//   auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
//   truncate(world, thetax);
//   auto rsh = apply(world, bsh_x, thetax);
//   rsh = gs.Qhat(rsh); // project out the ground state density from the response
//   return rsh;
// }
//
// vector_real_function_3d ResponseSolverPolicy<DynamicRestrictedResponse>::CoupledResponseEquations(
//     World &world, const GroundStateData &gs, const DynamicRestrictedResponse &vecs, const vector_real_function_3d &vp,
//     const std::vector<poperatorT> &bsh_x, const ResponseManager &rm, ResponseDebugLogger &logger) {
//
//   auto &x = vecs.x_alpha;
//   auto &y = vecs.y_alpha;
//   auto &all_x = vecs.flat;
//
//   auto num_orbitals = gs.orbitals.size();
//
//   std::vector<int> state_index;
//   std::vector<int> ii;
//
//   int i = 0;
//   for (int j = 0; j < 2 * num_orbitals; j++) {
//     state_index.push_back(0);
//     ii.push_back(i++);
//   }
//
//   ResponseComputeGroundExchange t0;
//   MacroTask g0_task(world, t0);
//   ResponseComputeGammaX tresponse;
//   MacroTask gx_task(world, tresponse);
//
//   vector_real_function_3d k0;
//   vector_real_function_3d gx;
//
//   DEBUG_TIMED_BLOCK(world, &logger, "g0_task", { k0 = g0_task(ii, state_index, all_x, gs.orbitals, false); });
//   DEBUG_TIMED_BLOCK(world, &logger, "gx_task", { gx = gx_task(ii, state_index, all_x, gs.orbitals, false); });
//
//   auto c_xc = gs.xcf_.hf_exchange_coefficient();
//   auto v_local = gs.V_local * all_x;
//   auto v0x = v_local - c_xc * k0;
//   auto epsilonx = transform(world, x, gs.Hamiltonian_no_diag, true);
//   auto epsilony = transform(world, y, gs.Hamiltonian_no_diag, true);
//   epsilonx.insert(epsilonx.end(), epsilony.begin(), epsilony.end());
//   auto thetax = -2.0 * (v0x - epsilonx + gx + vp);
//   truncate(world, thetax);
//   auto rsh = apply(world, bsh_x, thetax);
//   truncate(world, thetax);
//   rsh = gs.Qhat(rsh);
//   return rsh;
// }
