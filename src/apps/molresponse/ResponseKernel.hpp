#ifndef MOLRESPONSE_V3_RESPONSEKERNEL_HPP
#define MOLRESPONSE_V3_RESPONSEKERNEL_HPP

#include "GroundState.hpp"
#include "ResponseFunctions.hpp"

#include <madness/chem/SCFOperators.h>
#include <madness/mra/mra.h>
#include <madness/mra/operator.h>

namespace molresponse_v3 {

using namespace madness;

/// Response type determines density definition, exchange terms, BSH shifts.
/// The storage (ResponseState) is the same for all three — only the
/// algorithm differs.
enum class ResponseType { Static, Full, TDA };

// =========================================================================
// Response density
// =========================================================================

/// Compute response density from one spin channel (no spin factor).
///
/// Returns the bare spin-channel density:
///   Static: sum(x_i * phi_i)        [y implied = x]
///   Full:   sum((x_i + y_i) * phi_i) [x,y independent]
///   TDA:    sum(x_i * phi_i)          [y = 0]
///
/// The spin factor (2x for restricted Static/TDA) is applied by
/// compute_response_density, NOT here.
inline real_function_3d compute_spin_density(
    World &world, ResponseType type, const vector_real_function_3d &x,
    const vector_real_function_3d &y, const vector_real_function_3d &phi) {
  auto xphi = mul(world, x, phi, true);

  switch (type) {
  case ResponseType::Static:
    return sum(world, xphi, true);

  case ResponseType::Full: {
    auto yphi = mul(world, y, phi, true);
    return sum(world, xphi, true) + sum(world, yphi, true);
  }

  case ResponseType::TDA:
    return sum(world, xphi, true);
  }
  MADNESS_EXCEPTION("Unknown ResponseType", 0);
}

/// Compute the physical first-order response density.
///
/// Returns the full spin-summed response density δρ(r):
///
///   δρ = spin_factor × y_factor × Σ_channels [compute_spin_density]
///
///   spin_factor: 2 for restricted (each spatial orbital = 2 electrons),
///                1 for unrestricted (each spin-orbital = 1 electron)
///   y_factor:    2 for Static/TDA (y=x implied but not solved),
///                1 for Full (x and y both solved)
///
/// This is the PHYSICAL density — no halving convention. The Coulomb
/// coupling in compute_gamma uses J[ρ] directly (not 2*J).
inline real_function_3d compute_response_density(World &world,
                                                 ResponseType type,
                                                 const RealResponseState &state,
                                                 const GroundState &gs) {

  auto rho_alpha = compute_spin_density(world, type, state.x_alpha,
                                        state.y_alpha, gs.orbitals_alpha());

  real_function_3d rho;
  if (gs.is_spin_restricted()) {
    rho = rho_alpha;
  } else {
    auto rho_beta = compute_spin_density(world, type, state.x_beta,
                                         state.y_beta, gs.orbitals_beta());
    rho = rho_alpha + rho_beta;
  }

  // Physical occupancy factors
  double spin_factor = gs.is_spin_restricted() ? 2.0 : 1.0;
  double y_factor =
      (type == ResponseType::Static || type == ResponseType::TDA) ? 2.0 : 1.0;
  rho = (spin_factor * y_factor) * rho;
  return rho;
}

// =========================================================================
// Response coupling (gamma)
// =========================================================================

/// Compute response coupling potential for one spin channel.
///
/// gamma = Q * [2*J[rho]*phi - c_xc * exchange_terms]
///
/// Exchange terms depend on ResponseType:
///   Static: K[phi,x](phi) + K[x,phi](phi)
///   Full x: K[phi,x](phi) + K[y,phi](phi)
///   TDA:    K[phi,x](phi)  [single term]
///
/// For Full y-channel, call with x and y swapped.
///
/// Spin handling: phi, x, y must all be from the SAME spin channel.
/// Coulomb uses rho_total (sum over both spins). Exchange is same-spin
/// only — this matches SCF::apply_potential(ispin) which constructs
/// Exchange(world, calc, ispin) using only the ispin orbital set.
///
/// ispin parameter (0=alpha, 1=beta) is reserved for future XC kernel
/// (fxc) which needs spin-dependent second derivatives.
inline vector_real_function_3d compute_gamma(
    World &world, ResponseType type, const vector_real_function_3d &x,
    const vector_real_function_3d &y, const vector_real_function_3d &phi,
    const real_function_3d &rho_total, const poperatorT &coulop,
    const QProjector<double, 3> &Q, double c_xc, double lo, int ispin = 0) {

  long n = phi.size();

  // Coulomb: J[rho] * phi
  // rho_total is the physical response density (no halving convention),
  // so the coupling is simply J[rho]*phi — no factor of 2.
  auto J_rho = apply(*coulop, rho_total);
  auto Jphi = mul(world, J_rho, phi, true);
  J_rho.clear();
  auto gamma = std::move(Jphi);

  // Exchange terms (scaled by c_xc)
  if (c_xc > 0.0) {
    Exchange<double, 3> Kpx(world, lo);
    Kpx.set_bra_and_ket(phi, x);
    Kpx.set_algorithm(
        Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
    auto Kpx_phi = Kpx(phi); // K[phi,x](phi)

    switch (type) {
    case ResponseType::Static: {
      // K[phi,x](phi) + K[x,phi](phi)
      Exchange<double, 3> Kxp(world, lo);
      Kxp.set_bra_and_ket(x, phi);
      Kxp.set_algorithm(
          Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
      auto Kxp_phi = Kxp(phi); // K[x,phi](phi)
      gaxpy(world, 1.0, gamma, -c_xc, Kpx_phi);
      gaxpy(world, 1.0, gamma, -c_xc, Kxp_phi);
      break;
    }
    case ResponseType::Full: {
      // For x-channel: K[phi,x](phi) + K[y,phi](phi)
      Exchange<double, 3> Kyp(world, lo);
      Kyp.set_bra_and_ket(y, phi);
      Kyp.set_algorithm(
          Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
      auto Kyp_phi = Kyp(phi); // K[y,phi](phi)
      gaxpy(world, 1.0, gamma, -c_xc, Kpx_phi);
      gaxpy(world, 1.0, gamma, -c_xc, Kyp_phi);
      break;
    }
    case ResponseType::TDA: {
      // K[phi,x](phi) only
      gaxpy(world, 1.0, gamma, -c_xc, Kpx_phi);
      break;
    }
    }
  }

  // TODO: XC kernel fxc for DFT (hybrid and pure DFT functionals)
  // The remaining correlation beyond c_xc*K is captured by the XC kernel
  // second derivative fxc. For unrestricted, ispin selects spin-dependent
  // fxc (alpha or beta channel). See XCOperator(world, scf, ispin).
  // if (has_fxc) {
  //     XCOperator<double,3> xc_kernel(world, xc_string, spin_restricted,
  //                                     arho, brho, ispin);
  //     auto vxc_response = xc_kernel.apply_kernel(rho_response);
  //     gaxpy(world, 1.0, gamma, 1.0, mul(world, vxc_response, phi));
  // }

  // Q-project
  gamma = Q(gamma);
  truncate(world, gamma);
  return gamma;
}

// =========================================================================
// BSH operators
// =========================================================================

/// Compute BSH level-shift to keep mu^2 positive.
inline double compute_bsh_shift(const Tensor<double> &eps, double omega) {
  constexpr double guard = 0.05;
  double lumo_shifted = eps(eps.size() - 1) + omega;
  return (lumo_shifted >= 0.0) ? -guard - lumo_shifted : 0.0;
}

/// Create BSH operators for one spin channel.
///
/// mu_p = sqrt(-2 * (eps_p + omega + shift))
inline std::vector<poperatorT>
make_bsh_operators(World &world, const Tensor<double> &orbital_energies,
                   double omega, double lo) {

  double tol = FunctionDefaults<3>::get_thresh();
  double shift = compute_bsh_shift(orbital_energies, omega);
  long n = orbital_energies.size();

  std::vector<poperatorT> ops(n);
  for (long p = 0; p < n; p++) {
    double mu = std::sqrt(-2.0 * (orbital_energies(p) + omega + shift));
    ops[p] = poperatorT(BSHOperatorPtr3D(world, mu, lo, tol));
  }
  return ops;
}

/// Create BSH operators for y-channel (Full TDDFT only).
/// Uses -omega: mu_p = sqrt(-2 * (eps_p - omega))
inline std::vector<poperatorT>
make_bsh_operators_y(World &world, const Tensor<double> &orbital_energies,
                     double omega, double lo) {

  // y-channel: eps_p - omega < 0 guaranteed for bound states, no shift needed
  return make_bsh_operators(world, orbital_energies, -omega, lo);
}

// =========================================================================
// FD iteration step (one spin channel)
// =========================================================================

/// One BSH update step for frequency-dependent response (one spin channel).
///
/// Computes:
///   k0x = K[phi,phi](x)
///   v0x = V_local * x - c_xc * k0x
///   eps_x = transform(x, fock_no_diag)
///   theta = -2 * (v0x - eps_x + gamma + v_p)
///   x_new = Q * BSH(theta)
///
/// gamma must be precomputed via compute_gamma().
inline vector_real_function_3d fd_iteration_step(
    World &world, const vector_real_function_3d &x,
    const vector_real_function_3d &phi,
    const vector_real_function_3d &v_perturbation,
    const real_function_3d &V_local, const Tensor<double> &fock_no_diag,
    const vector_real_function_3d &gamma, const QProjector<double, 3> &Q,
    const std::vector<poperatorT> &bsh_ops, double c_xc, double lo) {

  double vtol = FunctionDefaults<3>::get_thresh() * 0.1;

  // Ground-state exchange: K[phi,phi](x)
  vector_real_function_3d k0x;
  if (c_xc > 0.0) {
    Exchange<double, 3> K0(world, lo);
    K0.set_bra_and_ket(phi, phi);
    K0.set_algorithm(
        Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
    k0x = K0(x);
  } else {
    k0x = zero_functions<double, 3>(world, x.size());
  }

  // V_local * x
  auto Vx = mul_sparse(world, V_local, x, vtol);

  // V0*x = V_local*x - c_xc * K0*x
  gaxpy(world, 1.0, Vx, -c_xc, k0x);
  k0x.clear();

  // Off-diagonal Fock coupling: sum_j F_ij * x_j
  auto eps_x = transform(world, x, fock_no_diag, vtol, true);

  // Theta = -2 * (V0*x - eps*x + gamma + v_p)
  // Assemble: theta = V0x - eps_x + gamma + v_p
  gaxpy(world, 1.0, Vx, -1.0, eps_x);         // Vx = V0x - eps_x
  gaxpy(world, 1.0, Vx, 1.0, gamma);          // Vx += gamma
  gaxpy(world, 1.0, Vx, 1.0, v_perturbation); // Vx += v_p
  scale(world, Vx, -2.0);                     // theta = -2 * (...)
  truncate(world, Vx);

  // Apply BSH
  auto x_new = apply(world, bsh_ops, Vx);

  // Q-project
  x_new = Q(x_new);
  truncate(world, x_new);

  return x_new;
}

/// Full FD iteration step handling both spins and x/y channels.
///
/// For Static/TDA: updates x only.
/// For Full: updates both x and y (with swapped exchange coupling for y).
/// For unrestricted: applies alpha and beta channels separately with
///   shared total response density but separate Fock/exchange/BSH.
///
/// @param debug  if >= 3, print per-channel norms after each step
inline RealResponseState
fd_iteration(World &world, ResponseType type, const RealResponseState &current,
             const RealResponseState &perturbation, const GroundState &gs,
             const poperatorT &coulop,
             const std::vector<poperatorT> &bsh_alpha_x,
             const std::vector<poperatorT> &bsh_alpha_y, // empty for Static/TDA
             const std::vector<poperatorT> &bsh_beta_x,  // empty if restricted
             const std::vector<poperatorT>
                 &bsh_beta_y, // empty if restricted or Static/TDA
             double omega, int debug = 0) {

  // Total response density (shared across spins for Coulomb)
  auto rho = compute_response_density(world, type, current, gs);

  double c_xc = gs.hf_exchange_coefficient();
  double lo = gs.params().lo();
  RealResponseState result;

  // --- Alpha x-channel (ispin=0) ---
  auto gamma_alpha_x = compute_gamma(
      world, type, current.x_alpha, current.y_alpha, gs.orbitals_alpha(), rho,
      coulop, gs.Q_alpha(), c_xc, lo, /*ispin=*/0);

  result.x_alpha =
      fd_iteration_step(world, current.x_alpha, gs.orbitals_alpha(),
                        perturbation.x_alpha, gs.V_local(), gs.focka_no_diag(),
                        gamma_alpha_x, gs.Q_alpha(), bsh_alpha_x, c_xc, lo);

  // --- Alpha y-channel (Full only, ispin=0) ---
  if (type == ResponseType::Full) {
    // For y-channel: exchange terms are K[phi,y]+K[x,phi] (swapped from x)
    auto gamma_alpha_y = compute_gamma(
        world, type, current.y_alpha, current.x_alpha, // swap x and y
        gs.orbitals_alpha(), rho, coulop, gs.Q_alpha(), c_xc, lo, /*ispin=*/0);

    result.y_alpha =
        fd_iteration_step(world, current.y_alpha, gs.orbitals_alpha(),
                          perturbation.y_alpha.empty() ? perturbation.x_alpha
                                                       : perturbation.y_alpha,
                          gs.V_local(), gs.focka_no_diag(), gamma_alpha_y,
                          gs.Q_alpha(), bsh_alpha_y, c_xc, lo);
  }

  // --- Beta channels (unrestricted only, ispin=1) ---
  if (!gs.is_spin_restricted()) {
    auto gamma_beta_x = compute_gamma(
        world, type, current.x_beta, current.y_beta, gs.orbitals_beta(), rho,
        coulop, gs.Q_beta(), c_xc, lo, /*ispin=*/1);

    result.x_beta =
        fd_iteration_step(world, current.x_beta, gs.orbitals_beta(),
                          perturbation.x_beta, gs.V_local(), gs.fockb_no_diag(),
                          gamma_beta_x, gs.Q_beta(), bsh_beta_x, c_xc, lo);

    if (type == ResponseType::Full) {
      auto gamma_beta_y = compute_gamma(
          world, type, current.y_beta, current.x_beta, gs.orbitals_beta(), rho,
          coulop, gs.Q_beta(), c_xc, lo, /*ispin=*/1);

      result.y_beta =
          fd_iteration_step(world, current.y_beta, gs.orbitals_beta(),
                            perturbation.y_beta.empty() ? perturbation.x_beta
                                                        : perturbation.y_beta,
                            gs.V_local(), gs.fockb_no_diag(), gamma_beta_y,
                            gs.Q_beta(), bsh_beta_y, c_xc, lo);
    }
  }

  // Debug: per-channel norms (collective — all ranks must participate in
  // norm2s)
  if (debug >= 3) {
    auto print_norm = [&](const char *name, const vector_real_function_3d &v) {
      if (!v.empty()) {
        auto norms = norm2s(world, v);
        double total = 0;
        for (auto n : norms)
          total += n * n;
        if (world.rank() == 0)
          print("  fd_iter", name, "norm=", std::sqrt(total));
      }
    };
    print_norm("x_alpha", result.x_alpha);
    print_norm("y_alpha", result.y_alpha);
    print_norm("x_beta", result.x_beta);
    print_norm("y_beta", result.y_beta);
  }

  return result;
}

// =========================================================================
// ES building blocks (excited-state specific)
// =========================================================================

/// Compute Lambda: potential operator matrix on response space.
///
/// Lambda = (V0*x - eps_diag*x + gamma)
///
/// Unlike FD theta which uses V_local (no kinetic), Lambda uses the full
/// Fock matrix for the diagonal term. This is used for subspace rotation.
inline vector_real_function_3d compute_lambda(
    World &world, ResponseType type, const vector_real_function_3d &x,
    const vector_real_function_3d &y, const vector_real_function_3d &phi,
    const real_function_3d &V_local, const Tensor<double> &fock_no_diag,
    const real_function_3d &rho_response, const poperatorT &coulop,
    const QProjector<double, 3> &Q, double c_xc, double lo) {

  double vtol = FunctionDefaults<3>::get_thresh() * 0.1;

  // Ground-state potential
  vector_real_function_3d k0x;
  if (c_xc > 0.0) {
    Exchange<double, 3> K0(world, lo);
    K0.set_bra_and_ket(phi, phi);
    K0.set_algorithm(
        Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
    k0x = K0(x);
  } else {
    k0x = zero_functions<double, 3>(world, x.size());
  }

  auto Vx = mul_sparse(world, V_local, x, vtol);
  gaxpy(world, 1.0, Vx, -c_xc, k0x); // V0*x

  // Off-diagonal Fock
  auto eps_x = transform(world, x, fock_no_diag, vtol, true);

  // Response coupling
  auto gamma =
      compute_gamma(world, type, x, y, phi, rho_response, coulop, Q, c_xc, lo);

  // Lambda = V0*x - eps*x + gamma
  gaxpy(world, 1.0, Vx, -1.0, eps_x);
  gaxpy(world, 1.0, Vx, 1.0, gamma);
  truncate(world, Vx);

  return Vx;
}

/// Compute subspace potential matrix: A_ij = <x_i | lambda_j>
inline Tensor<double>
compute_subspace_matrix(World &world,
                        const std::vector<vector_real_function_3d> &states,
                        const std::vector<vector_real_function_3d> &lambda) {

  long m = states.size();
  Tensor<double> A(m, m);
  for (long i = 0; i < m; i++) {
    for (long j = 0; j < m; j++) {
      A(i, j) = madness::inner(states[i], lambda[j]);
    }
  }
  return A;
}

/// Compute overlap matrix: S_ij = <x_i | x_j>
/// For Full TDDFT with symplectic metric: S_ij = <x_i|x_j> - <y_i|y_j>
inline Tensor<double>
compute_overlap_matrix(World &world, ResponseType type,
                       const std::vector<vector_real_function_3d> &states_x,
                       const std::vector<vector_real_function_3d> &states_y) {

  long m = states_x.size();
  Tensor<double> S(m, m);
  for (long i = 0; i < m; i++) {
    for (long j = 0; j < m; j++) {
      S(i, j) = madness::inner(states_x[i], states_x[j]);
      if (type == ResponseType::Full && !states_y.empty()) {
        S(i, j) -= madness::inner(states_y[i], states_y[j]);
      }
    }
  }
  return S;
}

/// Apply T = -1/2 ∇² to a vector of functions via three gradient operators.
/// Used to assemble the full F0 = T + V_local - c_xc·K0 acting on response
/// functions; needed for the SUBSPACE matrix in ES (NOT for the BSH driver
/// where T is absorbed into the Helmholtz Green's function).
/// Direct port of molresponse_v2/ResponseKernels.hpp:apply_kinetic_flat.
inline vector_real_function_3d apply_kinetic(World &world,
                                             const vector_real_function_3d &v) {
  if (v.empty())
    return {};
  int dirs[] = {0, 1, 2};
  vector_real_function_3d result;
  for (int d : dirs) {
    real_derivative_3d D(world, d);
    auto dv = apply(world, D, v);
    auto dv2 = apply(world, D, dv);

    result += dv2;
  }
  scale(world, result, -0.5);
  truncate(world, result);
  return result;

  // real_derivative_3d Dx(world, 0);
  // real_derivative_3d Dy(world, 1);
  // real_derivative_3d Dz(world, 2);
  // auto dvx  = apply(world, Dx, v);
  // auto dvy  = apply(world, Dy, v);
  // auto dvz  = apply(world, Dz, v);
  // auto dvx2 = apply(world, Dx, dvx);
  // auto dvy2 = apply(world, Dy, dvy);
  // auto dvz2 = apply(world, Dz, dvz);
  // auto result = dvx2;
  // gaxpy(world, 1.0, result, 1.0, dvy2);
  // gaxpy(world, 1.0, result, 1.0, dvz2);
}

/// Compute the FULL response action Λ for SUBSPACE matrix construction:
///   Λ·X = (T + V_local − c_xc·K0)·X  −  hamiltonian·X  +  γ
///
/// vs `compute_lambda` (which is actually Θ and is correct for the BSH
/// driver but missing T and the diagonal of the Fock matrix), this function
/// includes T·X via `apply_kinetic` and uses the *full* hamiltonian
/// (`gs.focka()` not `focka_no_diag()`). It mirrors the legacy
/// `Compute_Lambda_X` (Lambda_X.cc:33-93): Λ = (F0X − E0X) + γ.
///
/// Use `compute_lambda` for the BSH driver (theta = −2·Λ_theta), and use
/// `compute_lambda_subspace` for `compute_subspace_matrix` (so the
/// generalized eigenvalues from `diagonalize_subspace` are the true
/// excitation energies).
inline vector_real_function_3d compute_lambda_subspace(
    World &world, ResponseType type, const vector_real_function_3d &x,
    const vector_real_function_3d &y, const vector_real_function_3d &phi,
    const real_function_3d &V_local, const Tensor<double> &fock_full,
    const real_function_3d &rho_response, const poperatorT &coulop,
    const QProjector<double, 3> &Q, double c_xc, double lo) {

  double vtol = FunctionDefaults<3>::get_thresh() * 0.1;

  // T·X
  auto Tx = apply_kinetic(world, x);

  // Ground-state exchange K0·X
  vector_real_function_3d k0x;
  if (c_xc > 0.0) {
    Exchange<double, 3> K0(world, lo);
    K0.set_bra_and_ket(phi, phi);
    K0.set_algorithm(
        Exchange<double, 3>::ExchangeAlgorithm::multiworld_efficient_row);
    k0x = K0(x);
  } else {
    k0x = zero_functions<double, 3>(world, x.size());
  }

  // V_local·X − c_xc·K0·X = V0·X
  auto Vx = mul_sparse(world, V_local, x, vtol);
  gaxpy(world, 1.0, Vx, -c_xc, k0x);

  // F0·X = T·X + V0·X
  auto F0x = Tx;
  gaxpy(world, 1.0, F0x, 1.0, Vx);

  // E0·X = transform(X, hamiltonian_full)
  auto E0x = transform(world, x, fock_full, vtol, true);

  // Response coupling γ
  auto gamma =
      compute_gamma(world, type, x, y, phi, rho_response, coulop, Q, c_xc, lo);

  // Λ = F0·X − E0·X + γ
  gaxpy(world, 1.0, F0x, -1.0, E0x);
  gaxpy(world, 1.0, F0x, 1.0, gamma);
  truncate(world, F0x);

  return F0x;
}

/// Diagonalize generalized eigenvalue problem A*U = S*U*diag(omega).
/// Returns {eigenvalues, transformation_matrix}.
inline std::pair<Tensor<double>, Tensor<double>>
diagonalize_subspace(const Tensor<double> &A, const Tensor<double> &S) {

  Tensor<double> evals;
  Tensor<double> U;
  sygvp(World::get_default(), A, S, 1, U, evals);
  return {evals, U};
}

// =========================================================================
// ES bundle helpers — shared by all ES stages (TDA RHF, Full RHF,
// TDA UHF, Full UHF). See docs/11_excited_state_solver_design.md.
// =========================================================================

/// Inner product between two states under the type-appropriate metric:
///   TDA / Static : ⟨a.x|b.x⟩  (positive-definite, identity)
///   Full         : ⟨a.x|b.x⟩ - ⟨a.y|b.y⟩  (symplectic, indefinite)
/// Spin: alpha + beta channels are summed when present.
inline double bundle_inner(World &world, ResponseType type,
                           const RealResponseState &a,
                           const RealResponseState &b) {
  double s = 0.0;
  if (!a.x_alpha.empty() && !b.x_alpha.empty())
    s += madness::inner(a.x_alpha, b.x_alpha);
  if (!a.x_beta.empty() && !b.x_beta.empty())
    s += madness::inner(a.x_beta, b.x_beta);
  if (type == ResponseType::Full) {
    if (!a.y_alpha.empty() && !b.y_alpha.empty())
      s -= madness::inner(a.y_alpha, b.y_alpha);
    if (!a.y_beta.empty() && !b.y_beta.empty())
      s -= madness::inner(a.y_beta, b.y_beta);
  }
  return s;
}

/// Gram-Schmidt orthonormalization of a bundle of response states under
/// the type-appropriate metric.
///
/// For TDA the metric is positive-definite; standard Gram-Schmidt works.
/// For Full the metric is indefinite (S = ⟨x|x⟩ - ⟨y|y⟩) and exact
/// metric-orthonormalization needs hyperbolic Gram-Schmidt. For now we
/// fall back to taking |⟨a|a⟩|^{1/2} as the norm — adequate when the
/// guess is already roughly metric-orthonormal (the legacy code does
/// this same pre-iter pass at TDDFT.cc:1667-1672).
inline void orthonormalize_bundle(World &world, ResponseType type,
                                  std::vector<RealResponseState> &X) {
  long n = static_cast<long>(X.size());
  if (n == 0)
    return;

  auto subtract = [&](RealResponseState &dst, double c,
                      const RealResponseState &src) {
    if (!src.x_alpha.empty())
      gaxpy(world, 1.0, dst.x_alpha, -c, src.x_alpha);
    if (!src.x_beta.empty())
      gaxpy(world, 1.0, dst.x_beta, -c, src.x_beta);
    if (!src.y_alpha.empty())
      gaxpy(world, 1.0, dst.y_alpha, -c, src.y_alpha);
    if (!src.y_beta.empty())
      gaxpy(world, 1.0, dst.y_beta, -c, src.y_beta);
  };
  auto rescale = [&](RealResponseState &s, double f) {
    if (!s.x_alpha.empty())
      scale(world, s.x_alpha, f);
    if (!s.x_beta.empty())
      scale(world, s.x_beta, f);
    if (!s.y_alpha.empty())
      scale(world, s.y_alpha, f);
    if (!s.y_beta.empty())
      scale(world, s.y_beta, f);
  };

  for (long i = 0; i < n; i++) {
    for (long j = 0; j < i; j++) {
      double c = bundle_inner(world, type, X[j], X[i]);
      subtract(X[i], c, X[j]);
    }
    double norm2_val = bundle_inner(world, type, X[i], X[i]);
    double n_abs = std::sqrt(std::abs(norm2_val));
    if (n_abs > 1e-12)
      rescale(X[i], 1.0 / n_abs);
  }
}

/// Apply rotation U to the bundle: X_new[i] = sum_j U(j, i) * X[j].
/// In place. U is N×N where N = X.size().
inline void rotate_bundle(World &world, std::vector<RealResponseState> &X,
                          const Tensor<double> &U) {
  long n = static_cast<long>(X.size());
  if (n == 0)
    return;

  long na = X[0].num_alpha();
  long nb = X[0].is_restricted() ? 0 : X[0].num_beta();
  bool has_y = X[0].has_y();

  std::vector<RealResponseState> Y(n);
  for (long i = 0; i < n; i++) {
    Y[i] = RealResponseState::allocate(world, na, nb, has_y);
  }

  // Y[new] += U(old, new) * X[old]
  for (long new_idx = 0; new_idx < n; new_idx++) {
    for (long old_idx = 0; old_idx < n; old_idx++) {
      double c = U(old_idx, new_idx);
      gaxpy(world, 1.0, Y[new_idx].x_alpha, c, X[old_idx].x_alpha);
      if (has_y)
        gaxpy(world, 1.0, Y[new_idx].y_alpha, c, X[old_idx].y_alpha);
      if (nb > 0) {
        gaxpy(world, 1.0, Y[new_idx].x_beta, c, X[old_idx].x_beta);
        if (has_y)
          gaxpy(world, 1.0, Y[new_idx].y_beta, c, X[old_idx].y_beta);
      }
    }
  }

  double thresh = FunctionDefaults<3>::get_thresh();
  for (auto &s : Y) {
    truncate(world, s.x_alpha, thresh);
    if (has_y)
      truncate(world, s.y_alpha, thresh);
    if (nb > 0) {
      truncate(world, s.x_beta, thresh);
      if (has_y)
        truncate(world, s.y_beta, thresh);
    }
  }
  X = std::move(Y);
}

/// Apply the same N×N rotation U to a vector of vector_real_function_3d
/// (e.g., a bundle of single-channel state-vectors like Λ_X). Same
/// convention: out[new] = sum_old U(old, new) * in[old].
inline void rotate_vector_bundle(World &world,
                                 std::vector<vector_real_function_3d> &X,
                                 const Tensor<double> &U) {

  long n = static_cast<long>(X.size());
  if (n == 0)
    return;
  long n_orb = static_cast<long>(X[0].size());

  std::vector<vector_real_function_3d> Y(n);
  for (long i = 0; i < n; i++) {
    Y[i] = zero_functions<double, 3>(world, n_orb);
  }
  for (long new_idx = 0; new_idx < n; new_idx++) {
    for (long old_idx = 0; old_idx < n; old_idx++) {
      double c = U(old_idx, new_idx);
      gaxpy(world, 1.0, Y[new_idx], c, X[old_idx]);
    }
  }
  double thresh = FunctionDefaults<3>::get_thresh();
  for (auto &v : Y)
    truncate(world, v, thresh);
  X = std::move(Y);
}

/// Build the ES theta for one root, one spin channel:
///   theta_for_BSH = -2 * (Lambda + shift·x)
/// where `shift` matches the level-shift used inside `make_bsh_operators`
/// when ε+ω > 0 (so μ² stays negative).  When shift==0, this reduces to
///   theta = -2 * Lambda
/// Mirrors molresponse_legacy/TDDFT.cc:1830-1831 where `apply_shift`
/// adds `shift·X` to theta_X before the *(-2) and BSH apply.
inline vector_real_function_3d
build_es_theta(World &world, const vector_real_function_3d &Lambda,
               const vector_real_function_3d &x, double shift = 0.0) {
  auto theta = copy(world, Lambda, true);
  if (shift != 0.0) {
    gaxpy(world, 1.0, theta, shift, x);
  }
  scale(world, theta, -2.0);
  truncate(world, theta);
  return theta;
}

// =========================================================================
// Property computation helpers
// =========================================================================

/// Alpha factor for polarizability: alpha = factor * (ip_xa + ip_ya + ip_xb +
/// ip_yb)
///
/// The physical formula is:  α = -Σ_k n_k [⟨x_k|vp_k⟩ + ⟨y_k|vp_k⟩]
/// where n_k = 2 (restricted) or 1 (unrestricted per spin channel).
///
/// For Static/TDA, y is not solved (y=x implied) so ip_ya = ip_yb = 0
/// — the y_factor of 2 in af absorbs the missing y contribution.
/// Empirically validated <1% vs Dalton (d-aug-cc-pVQZ) on H2 / H2O /
/// C2H4 for Static ClosedShell, ~0.6% on Li UHF Static.
///
///   factor = -(spin_factor × y_factor)
///
///   spin_factor: 2 restricted (n_k=2), 1 unrestricted (n_k=1)
///   y_factor:    2 Static/TDA (y=x implied),
///                1 Full (y is solved separately)
///
///   Restricted Static:     -(2×2) = -4.0
///   Restricted Dynamic:    -(2×1) = -2.0
///   Unrestricted Static:   -(1×2) = -2.0
///   Unrestricted Dynamic:  -(1×1) = -1.0
///
/// Practical caveat: the density-scale factor inside `compute_density`
/// (e.g. 4 for Static ClosedShell, 2 for Full ClosedShell) enters the
/// Coulomb feedback inside compute_gamma and so is part of the self-
/// consistent equation. The (density-factor, alpha-factor) pair has to
/// be chosen consistently — if you halve af without compensating
/// somewhere in the kernel, α will shift (confirmed empirically:
/// F=4 + af=-2 halved α on H2).
inline double alpha_factor(ResponseType type, bool restricted) {
  double spin_factor = restricted ? 2.0 : 1.0;
  double y_factor = (type == ResponseType::Full) ? 1.0 : 2.0;
  return -(spin_factor * y_factor);
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_RESPONSEKERNEL_HPP
