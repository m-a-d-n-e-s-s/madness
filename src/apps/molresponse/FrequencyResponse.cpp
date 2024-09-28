//
// Created by adrianhurtado on 2/3/22.
//

#include "FrequencyResponse.hpp"
#include "funcdefaults.h"
#include "functypedefs.h"
#include "timer.h"
#include "x_space.h"
#include <filesystem>

namespace fs = std::filesystem;

void FrequencyResponse::initialize(World& world) {
  if (world.rank() == 0) {
    print("FrequencyResponse::initialize()");
  }
  Chi = PQ.copy();
}

void FrequencyResponse::iterate(World& world) {
  size_t iter;
  // Variables needed to iterate
  madness::QProjector<double, 3> projector(world, ground_orbitals);
  size_t n = r_params.num_orbitals();  // num orbitals
  size_t m = r_params.num_states();    // num response states

  real_function_3d v_xc;

  const double dconv =
      std::max(FunctionDefaults<3>::get_thresh(), r_params.dconv());
  auto thresh = FunctionDefaults<3>::get_thresh();
  auto density_target =
      dconv * static_cast<double>(std::max(size_t(5.0), molecule.natom()));
  const double x_residual_target = density_target * 5.0;
  Tensor<double> x_residual((int(m)));
  Tensor<double> delta_density((int(m)));

  bool static_res = (omega == 0.0);
  bool compute_y = not static_res;
  int r_vector_size;
  all_done = false;
  r_vector_size = (compute_y) ? static_cast<int>(2 * n) : static_cast<int>(n);

  vecfuncT rho_omega_old(m);
  // initialize DFT XC functional operator
  XCOperator<double, 3> xc = make_xc_operator(world);
  // create X space residuals
  X_space residuals = X_space::zero_functions(world, m, n);
  // create a std vector of XNONLinearsolvers
  response_solver kain_x_space;
  for (size_t b = 0; b < m; b++) {
    kain_x_space.emplace_back(response_matrix_allocator(world, r_vector_size),
                              false);
  }
  if (r_params.kain()) {
    for (auto& kain_space_b : kain_x_space) {
      kain_space_b.set_maxsub(static_cast<int>(r_params.maxsub()));
    }
  }
  // We compute with positive frequencies
  if (world.rank() == 0) {
    print("Warning input frequency is assumed to be positive");
    print("Computing at positive frequency omega = ", omega);
  }
  double x_shifts = 0.0;
  double y_shifts = 0.0;
  // if less negative orbital energy + frequency is positive or greater than 0
  if ((ground_energies[long(n) - 1] + omega) >= 0.0) {
    // Calculate minimum shift needed such that \eps + \omega + shift < 0
    print("*** we are shifting just so you know!!!");
    x_shifts = -.05 - (omega + ground_energies[long(n) - 1]);
  }
  auto bsh_x_ops = make_bsh_operators_response(world, x_shifts, omega);
  std::vector<poperatorT> bsh_y_ops;
  bsh_y_ops = (compute_y) ? make_bsh_operators_response(world, y_shifts, -omega)
                          : bsh_x_ops;
  auto max_rotation = .5 * x_residual_target + x_residual_target;
  PQ = generator(world, *this);
  PQ.truncate();

  vector<bool> converged(Chi.num_states(), false);
  Chi.reset_active();
  // make density for the first time
  auto rho_omega = response_context.compute_density(
      world, Chi, ground_orbitals, vector_real_function_3d(Chi.num_states()),
      false);

  for (iter = 0; iter < r_params.maxiter(); ++iter) {
    // if (world.rank() == 0) { print("At the start of iterate x", checkx); }
    iter_timing.clear();

    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
      if (world.rank() == 0)
        printf("\n   Iteration %d at time %.1fs\n", static_cast<int>(iter),
               wall_time());
      if (world.rank() == 0)
        print("-------------------------------------------");
    }
    if (iter < 2 || (iter % 5) == 0) {
      load_balance_chi(world);
    }
    if (iter > 0) {
      if (delta_density.max() > 20 && iter > 5) {
        if (world.rank() == 0) {
          print("d-residual > 20...break");
        }
        break;
      }

      auto chi_norms = (compute_y) ? Chi.norm2s() : Chi.x.norm2();
      auto rho_norms = madness::norm2s_T(world, rho_omega);

      // Todo add chi norm and chi_x
      if (world.rank() == 0) {
        function_data_to_json(j_molresponse, iter, chi_norms, x_residual,
                              rho_norms, delta_density);
        frequency_to_json(j_molresponse, iter, polar);
      }
      if (r_params.print_level() >= 1) {
        if (world.rank() == 0) {
          print("r_params.dconv(): ", r_params.dconv());
          print("thresh: ", FunctionDefaults<3>::get_thresh());
          print("k: ", FunctionDefaults<3>::get_k());
          print("Chi Norms at start of iteration: ", iter);
          print("||X||: ", chi_norms);
          print("||f(x)-x||: ", x_residual);
          print("<< XI | XJ >>(omega): \n", polar);
          print("targets : ||x||", x_residual_target, "    ||delta_rho||",
                density_target);
        }
      }
      auto check_convergence = [&](auto& ri, auto& di) {
        if (world.rank() == 0) {
          print("              ", ri, "    ", di);
        }
        return ((ri < x_residual_target) && (di < density_target));
      };

      for (const auto& b : Chi.active) {
        converged[b] = check_convergence(x_residual[static_cast<int>(b)],
                                         delta_density[static_cast<int>(b)]);
      }
      int b = 0;
      auto remove_converged = [&]() {
        Chi.reset_active();
        Chi.active.remove_if([&](auto x) { return converged[b++]; });
      };
      remove_converged();

      if (world.rank() == 0) {
        print("converged", converged);
        print("active", Chi.active);
      }
      b = 0;
      all_done = std::all_of(converged.begin(), converged.end(),
                             [](const auto& ci) { return ci; });
      if (all_done || iter == r_params.maxiter()) {
        // if converged print converged
        if (world.rank() == 0 && all_done and (r_params.print_level() > 1)) {
          print("\nConverged!\n");
        }
        if (r_params.save()) {
          molresponse::start_timer(world);
          save(world, r_params.save_file());
          if (r_params.print_level() >= 1)
            molresponse::end_timer(world, "Save:");
        }

        break;
      }
    }
    auto x_inner = ((compute_y) ? 2 : 1) * response_context.inner(Chi, Chi);

    auto rho_omega_norm = norm2s_T(world, rho_omega);
    auto [new_chi, new_res, new_rho] = update_response(
        world, Chi, xc, bsh_x_ops, bsh_y_ops, projector, x_shifts, omega,
        kain_x_space, iter, max_rotation, rho_omega, x_residual, residuals);

    auto old_rho = copy(world, rho_omega);
    rho_omega = copy(world, new_rho);

    for (const auto& b : Chi.active) {
      auto drho_b = rho_omega[b] - old_rho[b];
      auto drho_b_norm = drho_b.norm2();
      world.gop.fence();
      delta_density[static_cast<int>(b)] = drho_b_norm;
    }
    world.gop.fence();

    auto old_density_residual = copy(delta_density);

    Chi = new_chi.copy();

    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
    }
    x_residual = copy(new_res.residual_norms);
    if (r_params.print_level() >= 1) {
      molresponse::end_timer(world, "copy_response_data", "copy_response_data",
                             iter_timing);
    }

    auto dnorm = norm2s_T(world, rho_omega);

    polar = ((compute_y) ? -2 : -4) * response_context.inner(Chi, PQ);
    if (r_params.print_level() >= 1) {
      molresponse::end_timer(world, "Iteration Timing", "iter_total",
                             iter_timing);
    }
    time_data.add_data(iter_timing);
  }
  function_data.add_convergence_targets(FunctionDefaults<3>::get_thresh(),
                                        density_target, x_residual_target);

  Chi.reset_active();
  if (world.rank() == 0)
    print("\n");
  if (world.rank() == 0)
    print("   Finished Response Calculation ");
  if (world.rank() == 0)
    print("   ------------------------");
  if (world.rank() == 0)
    print("\n");

  // Did we converge?
  if (iter == r_params.maxiter() && not all_done) {
    if (world.rank() == 0)
      print("   Failed to converge. Reason:");
    if (world.rank() == 0)
      print("\n  ***  Ran out of iterations  ***\n");
  }

  if (world.rank() == 0) {
    print(" Final energy residuals X:");
    print(x_residual);
    print(" Final density change:");
    print(delta_density);
  }
  // compute_and_print_polarizability(world, Chi, PQ, "Converged");
}

auto FrequencyResponse::update_response(
    World& world, X_space& chi, XCOperator<double, 3>& xc,
    std::vector<poperatorT>& bsh_x_ops, std::vector<poperatorT>& bsh_y_ops,
    QProjector<double, 3>& projector, double& x_shifts, double& omega_n,
    response_solver& kain_x_space, size_t iteration, const double& max_rotation,
    const vector_real_function_3d& rho_old, const Tensor<double>& old_residuals,
    const X_space& xres_old)
    -> std::tuple<X_space, residuals, vector_real_function_3d> {

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }

  auto x = chi.copy();
  auto chi_norm = x.norm2s();
  // if (world.rank() == 0) { print("x before bsh update: ", chi_norm); }
  X_space theta_X =
      compute_theta_X(world, x, rho_old, xc, r_params.calc_type());
  X_space new_chi = bsh_update_response(world, theta_X, bsh_x_ops, bsh_y_ops,
                                        projector, x_shifts);

  chi_norm = new_chi.norm2s();
  // if (world.rank() == 0) { print("new_chi_norm after bsh update: ", chi_norm); }

  auto [new_res, bsh_norms] = update_residual(
      world, chi, new_chi, r_params.calc_type(), old_residuals, xres_old);
  if (iteration >= 0) {  // & (iteration % 3 == 0)) {
    new_chi = kain_x_space_update(world, chi, new_res, kain_x_space);
  }
  chi_norm = new_chi.norm2s();
  if (world.rank() == 0) {
    print("new_chi_norm after kain update: ", chi_norm);
  }

  // bool compute_y = r_params.calc_type() == "full";
  // x_space_step_restriction(world, chi, new_chi, compute_y, max_rotation);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "update response", "update", iter_timing);
  }

  auto new_rho = response_context.compute_density(
      world, new_chi, ground_orbitals, rho_old, true);

  return {new_chi, {new_res, bsh_norms}, new_rho};
}

auto FrequencyResponse::bsh_update_response(World& world, X_space& theta_X,
                                            std::vector<poperatorT>& bsh_x_ops,
                                            std::vector<poperatorT>& bsh_y_ops,
                                            QProjector<double, 3>& projector,
                                            double& x_shifts) -> X_space {
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  size_t m = theta_X.x.size();
  size_t n = theta_X.x.size_orbitals();
  bool compute_y = omega != 0.0;

  if (compute_y) {
    theta_X += theta_X * x_shifts;
    theta_X += PQ;
    theta_X = -2 * theta_X;
    theta_X.truncate();
  } else {
    theta_X.x += theta_X.x * x_shifts;
    theta_X.x += PQ.x;
    theta_X.x = theta_X.x * -2;
    theta_X.x.truncate_rf();
  }
  auto chi_norm = theta_X.norm2s();
  if (world.rank() == 0) {
    print("In bsh after theta_X+PQ: ", chi_norm);
  }
  // apply bsh
  X_space bsh_X = X_space::zero_functions(world, m, n);
  bsh_X.set_active(theta_X.active);
  bsh_X.x = apply(world, bsh_x_ops, theta_X.x);
  if (compute_y) {
    bsh_X.y = apply(world, bsh_y_ops, theta_X.y);
  }

  if (compute_y) {
    bsh_X.truncate();
  } else {
    bsh_X.x.truncate_rf();
  }

  chi_norm = bsh_X.norm2s();
  if (world.rank() == 0) {
    print("In bsh after apply: ", chi_norm);
  }
  auto apply_projector = [&](auto& xi) {
    return projector(xi);
  };
  if (compute_y) {
    bsh_X = oop_apply(bsh_X, apply_projector);
  } else {
    for (const auto& i : bsh_X.active)
      bsh_X.x[i] = projector(bsh_X.x[i]);
  }
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "bsh_update", "bsh_update", iter_timing);
  }
  if (compute_y) {
    bsh_X.truncate();
  } else {
    bsh_X.x.truncate_rf();
  }
  return bsh_X;
}

void FrequencyResponse::frequency_to_json(json& j_mol_in, size_t iter,
                                          const Tensor<double>& polar_ij) {
  json j = {};
  j["iter"] = iter;
  j["polar"] = tensor_to_json(polar_ij);
  auto index = j_mol_in["protocol_data"].size() - 1;
  j_mol_in["protocol_data"][index]["property_data"].push_back(j);
}

void FrequencyResponse::save(World& world, const std::string& name) {
  // Archive to write everything to
  archive::ParallelOutputArchive ar(world, name.c_str(), 1);

  ar & r_params.archive();
  ar & r_params.tda();
  ar & r_params.num_orbitals();
  ar & r_params.num_states();

  for (size_t i = 0; i < r_params.num_states(); i++)
    for (size_t j = 0; j < r_params.num_orbitals(); j++)
      ar & Chi.x[i][j];
  for (size_t i = 0; i < r_params.num_states(); i++)
    for (size_t j = 0; j < r_params.num_orbitals(); j++)
      ar & Chi.y[i][j];
}

// Load a response calculation
void FrequencyResponse::load(World& world, const std::string& name) {
  if (world.rank() == 0) {
    print("FrequencyResponse::load() -state");
  }
  // The archive to read from
  archive::ParallelInputArchive ar(world, name.c_str());
  ar & r_params.archive();
  ar & r_params.tda();
  ar & r_params.num_orbitals();
  ar & r_params.num_states();
  Chi = X_space(world, r_params.num_states(), r_params.num_orbitals());
  for (size_t i = 0; i < r_params.num_states(); i++)
    for (size_t j = 0; j < r_params.num_orbitals(); j++)
      ar & Chi.x[i][j];
  world.gop.fence();
  for (size_t i = 0; i < r_params.num_states(); i++)
    for (size_t j = 0; j < r_params.num_orbitals(); j++)
      ar & Chi.y[i][j];
  world.gop.fence();
}

auto nuclear_generator(World& world, ResponseBase& calc) -> X_space {
  auto [gc, molecule, r_params] = calc.get_parameter();
  X_space PQ(world, r_params.num_states(), r_params.num_orbitals());
  auto num_operators = size_t(molecule.natom() * 3);
  auto nuclear_vector = vecfuncT(num_operators);

  for (long atom = 0; atom < molecule.natom(); ++atom) {
    for (long axis = 0; axis < 3; ++axis) {
      functorT func(new madchem::MolecularDerivativeFunctor(
          molecule, static_cast<int>(atom), static_cast<int>(axis)));
      nuclear_vector.at(atom * 3 + axis) = functionT(factoryT(world)
                                                         .functor(func)
                                                         .nofence()
                                                         .truncate_on_project()
                                                         .truncate_mode(0));
    }
  }
  PQ.x = vector_to_PQ(world, nuclear_vector, calc.get_orbitals());
  PQ.y = PQ.x;
  return PQ;
}

auto dipole_generator(World& world, ResponseBase& calc) -> X_space {
  auto [gc, molecule, r_params] = calc.get_parameter();
  X_space PQ(world, r_params.num_states(), r_params.num_orbitals());
  vector_real_function_3d dipole_vectors(3);
  size_t i = 0;
  for (auto& d : dipole_vectors) {
    std::vector<int> f(3, 0);
    f[i++] = 1;
    d = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
  }
  // truncate(world, dipole_vectors, true);
  world.gop.fence();
  PQ.x = vector_to_PQ(world, dipole_vectors, calc.get_orbitals());
  PQ.y = PQ.x.copy();
  if (world.rank() == 0) {
    print("Made new PQ");
  }
  return PQ;
}

auto vector_to_PQ(World& world, const vector_real_function_3d& rhs_operators,
                  const vector_real_function_3d& ground_orbitals)
    -> response_space {
  response_space rhs(world, rhs_operators.size(), ground_orbitals.size());
  auto orbitals = copy(world, ground_orbitals);
  reconstruct(world, orbitals);
  truncate(world, orbitals);
  QProjector<double, 3> Qhat(world, orbitals);
  int b = 0;
  for (const functionT& pi : rhs_operators) {
    auto op_phi = mul(world, pi, ground_orbitals, true);
    rhs[b] = Qhat(op_phi);
    b++;
  }
  return rhs;
}

void QuadraticResponse::load(World& world, const std::string& name) {}
void QuadraticResponse::save(World& world, const std::string& name) {}

void QuadraticResponse::iterate(World& world) {}
void QuadraticResponse::initialize(World& world) {}

// To compute 2nd order we need rhs of xx,yy,zz and yz only.
auto QuadraticResponse::setup_XBC(World& world, const double& omega_b,
                                  const double& omega_c)
    -> std::pair<X_space, X_space> {

  this->index_B = {0, 0, 0, 1, 1, 1, 2, 2, 2};
  this->index_C = {0, 1, 2, 0, 1, 2, 0, 1, 2};

  auto B = x_data[1].first.copy();
  auto C = x_data[2].first.copy();

  auto num_states = index_B.size();
  auto num_orbitals = B.num_orbitals();

  auto new_B = X_space(world, num_states, num_orbitals);
  auto new_C = X_space(world, num_states, num_orbitals);

  for (int i = 0; i < num_states; i++) {
    new_B.x[i] = copy(world, B.x[index_B[i]]);
    new_B.y[i] = copy(world, B.y[index_B[i]]);

    new_C.x[i] = copy(world, C.x[index_C[i]]);
    new_C.y[i] = copy(world, C.y[index_C[i]]);
  }
  this->bc_directions.clear();
  for (int i = 0; i < num_states; i++) {
    auto bc_direction_i = std::string(xyz[index_B[i]] + xyz[index_C[i]]);
    this->bc_directions.push_back(bc_direction_i);
  }
  if (world.rank() == 0) {
    print("bc_directions: ", bc_directions);
  }

  return {new_B, new_C};
}

vector_real_function_3d QuadraticResponse::compute_vbc(
    World& world, const response_pair& B, const response_pair& C,
    const response_pair& BxCy, const response_pair& phiBC,
    const vector_real_function_3d& phi0, const real_function_3d& vb) {

  madness::QProjector<double, 3> Q(world, phi0);
  auto compute_g = [&](const vector_real_function_3d& x,
                       const vector_real_function_3d& y,
                       const vector_real_function_3d& phi) {
    const double lo = 1.e-10;
    Exchange<double, 3> k{world, lo};
    k.set_bra_and_ket(x, y);
    std::string algorithm_ = r_params.hfexalg();

    if (algorithm_ == "multiworld") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
    } else if (algorithm_ == "multiworld_row") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient_row);
    } else if (algorithm_ == "largemem") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
    } else if (algorithm_ == "smallmem") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
    }

    auto rho = sum(world, mul(world, x, y, true), true);
    auto tempJ = apply(*shared_coulomb_operator, rho);
    auto J = mul(world, tempJ, phi, true);
    auto K = k(phi);

    auto result = Q(2.0 * J - K);
    world.gop.fence();
    return result;
  };
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }

  auto gzeta = -1.0 * (compute_g(BxCy.x, BxCy.y, phi0) +
                       compute_g(phiBC.x, phiBC.y, phi0));

  auto FBX = -1.0 * (compute_g(B.x, phi0, C.x) + compute_g(phi0, B.y, C.x) +
                     Q(mul(world, vb, C.x, true)));

  // Terms that be added to VB
  auto FB = compute_g(B.x, phi0, phi0) + compute_g(phi0, B.y, phi0) +
            Q(mul(world, vb, phi0, true));
  auto matrix_fb = matrix_inner(world, phi0, FB);
  world.gop.fence();
  FB = transform(world, C.x, matrix_fb, true);

  return truncate(gzeta + FBX + FB, FunctionDefaults<3>::get_thresh(), true);
}
//
//
// <C;A,B> = <V(BC);X(A)> + <zeta(bc)_x| v(a) | zeta_(bc)_y> + < zeta(cb)_x| v(a) | zeta_(cb)_y >
//
//
// vc = x x x x y y y z z z
//
std::pair<Tensor<double>, std::vector<std::string>>
QuadraticResponse::compute_beta_tensor(World& world, const X_space& BC_left,
                                       const X_space& BC_right,
                                       const X_space& CB_left,
                                       const X_space& CB_right,
                                       const X_space& XA, const X_space& VBC) {
  auto create_dipole = [&]() {
    vector_real_function_3d dipole_vectors(3);
    size_t i = 0;
    // creates a vector of x y z dipole functions
    for (auto& d : dipole_vectors) {
      std::vector<int> f(3, 0);
      f[i++] = 1;
      d = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
    }
    return dipole_vectors;
  };

  auto dipole_vectors = create_dipole();  // x y z
  truncate(dipole_vectors, FunctionDefaults<3>::get_thresh(), true);

  int num_elements = static_cast<int>(XA.num_states() * BC_left.num_states());
  std::vector<std::string> beta_indices(num_elements);
  Tensor<double> beta(num_elements);

  int i = 0;
  for (int a = 0; a < XA.num_states(); a++) {
    for (int bc = 0; bc < BC_left.num_states(); bc++) {

      auto one =
          dot(world, BC_left.x[bc], BC_right.x[bc] * dipole_vectors[a], true);
      auto two =
          dot(world, BC_left.y[bc], BC_right.y[bc] * dipole_vectors[a], true);
      auto three =
          dot(world, CB_left.x[bc], CB_right.x[bc] * dipole_vectors[a], true);
      auto four =
          dot(world, CB_left.y[bc], CB_right.y[bc] * dipole_vectors[a], true);
      auto five = dot(world, XA.x[a], VBC.x[bc], true);
      auto six = dot(world, XA.y[a], VBC.y[bc], true);

      // Truncation here might be a bad idea, scheisse
      // one.truncate();
      // two.truncate();
      // three.truncate();
      // four.truncate();
      // five.truncate();
      // six.truncate();
      //
      auto one_trace = one.trace();
      auto two_trace = two.trace();
      auto three_trace = three.trace();
      auto four_trace = four.trace();
      auto five_trace = five.trace();
      auto six_trace = six.trace();

      beta[i] = one_trace + two_trace + three_trace + four_trace + five_trace +
                six_trace;

      beta_indices[i] = xyz[a] + bc_directions[bc];
      i++;
    }
  }

  return {-2.0 * beta, beta_indices};
}

std::pair<Tensor<double>, std::vector<std::string>>
QuadraticResponse::compute_beta_tensor_v2(World& world, const X_space& B,
                                          const X_space& C,
                                          const response_space& phiBC,
                                          const response_space& phiCB,
                                          const X_space& XA,
                                          const X_space& VBC) {
  auto create_dipole = [&]() {
    vector_real_function_3d dipole_vectors(3);
    size_t i = 0;
    // creates a vector of x y z dipole functions
    for (auto& d : dipole_vectors) {
      std::vector<int> f(3, 0);
      f[i++] = 1;
      d = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
    }
    return dipole_vectors;
  };

  auto dipole_vectors = create_dipole();  // x y z
  truncate(dipole_vectors, FunctionDefaults<3>::get_thresh(), true);

  int num_elements = static_cast<int>(XA.num_states() * phiBC.num_states);
  std::vector<std::string> beta_indices(num_elements);
  Tensor<double> beta(num_elements);

  int i = 0;
  for (int a = 0; a < XA.num_states(); a++) {
    int bc = 0;
    for (const auto& [b, c] : this->BC_index_pairs) {

      auto one = dot(world, B.x[b], C.y[c] * dipole_vectors[a]);
      auto three = dot(world, C.x[c], B.y[b] * dipole_vectors[a]);

      auto two = dot(world, phiBC.x[bc], ground_orbitals * dipole_vectors[a]);
      auto four = dot(world, phiCB.x[bc], ground_orbitals * dipole_vectors[a]);

      auto five = dot(world, XA.x[a], VBC.x[bc], true);
      auto six = dot(world, XA.y[a], VBC.y[bc], true);

      // Truncation here might be a bad idea, scheisse
      // one.truncate();
      // two.truncate();
      // three.truncate();
      // four.truncate();
      // five.truncate();
      // six.truncate();
      //
      auto one_trace = one.trace();
      auto two_trace = two.trace();
      auto three_trace = three.trace();
      auto four_trace = four.trace();
      auto five_trace = five.trace();
      auto six_trace = six.trace();

      beta[i] = one_trace + two_trace + three_trace + four_trace + five_trace +
                six_trace;

      beta_indices[i] = xyz[a] + (xyz[b] + xyz[c]);
      i++;
      bc++;
    }
  }

  return {-2.0 * beta, beta_indices};
}

std::pair<Tensor<double>, std::vector<std::string>>
QuadraticResponse::compute_beta_v2(World& world, const double& omega_b,
                                   const double& omega_c) {
  // step 0: construct all response vectors
  auto XA = -1.0 * x_data[0].first.copy();
  auto [XB, XC] = setup_XBC(world, omega_b, omega_c);
  X_space phi0 = X_space(world, XB.num_states(), XC.num_orbitals());
  for (auto i = 0; i < phi0.num_states(); i++) {
    phi0.x[i] = copy(world, ground_orbitals);
    phi0.y[i] = copy(world, ground_orbitals);
  }
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto [zeta_bc_left, zeta_bc_right, zeta_cb_left, zeta_cb_right] =
      compute_zeta_response_vectors(world, XB, XC);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "Zeta(BC) and Zeta(CB)");
  }

  auto B = x_data[1].first.copy();
  auto C = x_data[2].first.copy();

  path vbc_archive = "vbc_archive";

  // if vbc archive exists load it
  if (fs::exists(vbc_archive.replace_extension(".00000"))) {
    if (world.rank() == 0) {
      print("Loading VBC from archive");
    }
    this->VBC = load_x_space(world, vbc_archive.stem().string());
  } else {
    this->VBC = compute_second_order_perturbation_terms_v2(
        world, XB, XC, zeta_bc_left, zeta_bc_right, zeta_cb_left, zeta_cb_right,
        phi0);
    save_x_space(world, vbc_archive.string(), this->VBC);
  }
  // step 1: compute all exchange terms because they are the most expensive

  auto VBC_2 = compute_second_order_perturbation_terms_v3(
      world, B, C, zeta_bc_left.y, zeta_cb_left.y, ground_orbitals);

  auto rVBC = VBC_2 - VBC;
  auto rVBC_norm = rVBC.norm2s();
  if (world.rank() == 0) {
    print("rVBC_norm: ", rVBC_norm);
  }

  auto [beta0, beta0_dir] = compute_beta_tensor(
      world, zeta_bc_left, zeta_bc_right, zeta_cb_left, zeta_cb_right, XA, VBC);

  auto [beta2, beta2_dir] =
      compute_beta_tensor(world, zeta_bc_left, zeta_bc_right, zeta_cb_left,
                          zeta_cb_right, XA, VBC_2);

  if (world.rank() == 0) {
    print("beta0: ", beta0);
    print("beta2: ", beta2);
  }

  return {beta0, beta0_dir};
}

Tensor<double> QuadraticResponse::compute_beta(World& world) {

  // construct an X_space containing phi0 copies

  // bsh_X = oop_apply(bsh_X, apply_projector);
  QProjector<double, 3> projector(world, ground_orbitals);
  auto apply_projector = [&](auto& xi) {
    return projector(xi);
  };

  auto perturbation_A = generator(world, *this);
  auto XA = -1.0 * x_data[0].first.copy();

  // first step to compute beta is to construct the X_space representations of the virt/virt and occ/occ blocks of gamma

  auto [XB, XC] = setup_XBC(world, 0.0, 0.0);
  X_space phi0 = X_space(world, XB.num_states(), XC.num_orbitals());
  for (auto i = 0; i < phi0.num_states(); i++) {
    phi0.x[i] = copy(world, ground_orbitals);
    phi0.y[i] = copy(world, ground_orbitals);
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto [zeta_bc_left, zeta_bc_right, zeta_cb_left, zeta_cb_right] =
      compute_zeta_response_vectors(world, XB, XC);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "Zeta(BC) and Zeta(CB)");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto VBC = compute_second_order_perturbation_terms(
      world, XB, XC, zeta_bc_left, zeta_bc_right, zeta_cb_left, zeta_cb_right,
      phi0);

  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "V(BC) and V(CB)");
  }
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto beta = compute_beta_tensor(world, zeta_bc_left, zeta_bc_right,
                                  zeta_cb_left, zeta_cb_right, XA, VBC);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "compute_beta_tensor");
  }
  return beta.first;
}

std::pair<X_space, X_space>
QuadraticResponse::compute_first_order_fock_matrix_terms_v2(
    World& world, const X_space& B, const X_space& C, const X_space& g1b,
    const X_space& g1c, const X_space& VB, const X_space& VC,
    const X_space& phi0) const {

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto f1b = g1b + VB;
  auto f1c = g1c + VC;

  auto FBC = X_space(world, VB.num_states(), VB.num_orbitals());
  auto FCB = X_space(world, VB.num_states(), VB.num_orbitals());
  // Here contains the y components
  for (auto i = 0; i < VB.num_states(); i++) {

    auto fbx = matrix_inner(world, phi0.x[i], f1b.x[i]);
    auto fbx_dagger = matrix_inner(world, phi0.y[i], f1b.y[i]);
    auto fcx = matrix_inner(world, phi0.x[i], f1c.x[i]);
    auto fc_dagger = matrix_inner(world, phi0.y[i], f1c.y[i]);

    FBC.x[i] = copy(world, transform(world, C.x[i], fbx, true), true);
    FBC.y[i] = copy(world, transform(world, C.y[i], fbx_dagger, true), true);
    FCB.x[i] = copy(world, transform(world, B.x[i], fcx, true), true);
    FCB.y[i] = copy(world, transform(world, B.y[i], fc_dagger, true), true);
  }

  FBC.truncate();
  FCB.truncate();
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "Fock transformation terms");
  }
  return {FBC, FCB};
}
// computes <phi0|Fa|phi0> * XB
// where Fa=g1[xa]+va
std::pair<X_space, X_space>
QuadraticResponse::compute_first_order_fock_matrix_terms(
    World& world, const X_space& B, const X_space& phi0,
    const X_space& C) const {

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto g1a = compute_g1_term(world, B, phi0, phi0);
  auto g1b = compute_g1_term(world, C, phi0, phi0);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "VBC: first_order_terms: compute g1 terms");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto [VB, VC] = dipole_perturbation(world, phi0, phi0);

  auto f1a = g1a + VB;
  auto f1b = g1b + VC;

  auto FAXB = X_space(world, B.num_states(), B.num_orbitals());
  auto FBXA = X_space(world, B.num_states(), B.num_orbitals());
  // Here contains the y components
  for (auto i = 0; i < B.num_states(); i++) {

    auto fax = matrix_inner(world, phi0.x[i], f1a.x[i]);
    auto fax_dagger = matrix_inner(world, phi0.y[i], f1a.y[i]);
    auto fbx = matrix_inner(world, phi0.x[i], f1b.x[i]);
    auto fb_dagger = matrix_inner(world, phi0.y[i], f1b.y[i]);

    FAXB.x[i] = copy(world, transform(world, C.x[i], fax, true), true);
    FAXB.y[i] = copy(world, transform(world, C.y[i], fax_dagger, true), true);
    FBXA.x[i] = copy(world, transform(world, B.x[i], fbx, true), true);
    FBXA.y[i] = copy(world, transform(world, B.y[i], fb_dagger, true), true);
  }

  FAXB.truncate();
  FBXA.truncate();

  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "VBC: first_order_terms: compute FBC terms");
  }
  world.gop.fence();

  return {FAXB, FBXA};
}
std::tuple<X_space, X_space, X_space, X_space, X_space, X_space>
QuadraticResponse::compute_beta_coulomb(World& world, const X_space& B,
                                        const X_space& C,
                                        const X_space& zeta_bc_left,
                                        const X_space& zeta_bc_right,
                                        const X_space& zeta_cb_left,
                                        const X_space& zeta_cb_right,
                                        const X_space& phi0) {
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto zeta_bc = compute_coulomb_term(world, zeta_bc_left, zeta_bc_right, phi0);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "j: zeta_bc");
  }
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto zeta_cb = compute_coulomb_term(world, zeta_cb_left, zeta_cb_right, phi0);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "j: zeta_cb");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto bxc = compute_coulomb_term(world, B, phi0, C);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "j: bphi0c");
  }
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto cxb = compute_coulomb_term(world, C, phi0, B);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "j: cphi0B");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto bphi0 = compute_coulomb_term(world, B, phi0, phi0);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "j: bphi0_phi0");
  }
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto cphi0 = compute_coulomb_term(world, C, phi0, phi0);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "j: cphi0_phi0");
  }

  world.gop.fence();

  return {zeta_bc, zeta_cb, bxc, cxb, bphi0, cphi0};
}
std::tuple<X_space, X_space, X_space, X_space, X_space, X_space>
QuadraticResponse::compute_beta_exchange(World& world, const X_space& B,
                                         const X_space& C,
                                         const X_space& zeta_bc_left,
                                         const X_space& zeta_bc_right,
                                         const X_space& zeta_cb_left,
                                         const X_space& zeta_cb_right,
                                         const X_space& phi0) {
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto zeta_bc =
      compute_exchange_term(world, zeta_bc_left, zeta_bc_right, phi0);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "k: zeta_bc");
  }
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto zeta_cb =
      compute_exchange_term(world, zeta_cb_left, zeta_cb_right, phi0);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "k: zeta_cb");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto bxc = compute_exchange_term(world, B, phi0, C);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "k: bphi0c");
  }
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto cxb = compute_exchange_term(world, C, phi0, B);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "k: cphi0B");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto bphi0 = compute_exchange_term(world, B, phi0, phi0);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "k: bphi0_phi0");
  }
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto cphi0 = compute_exchange_term(world, C, phi0, phi0);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "k: cphi0_phi0");
  }

  world.gop.fence();
  return {zeta_bc, zeta_cb, bxc, cxb, bphi0, cphi0};
}

X_space QuadraticResponse::compute_second_order_perturbation_terms_v2(
    World& world, const X_space& B, const X_space& C,
    const X_space& zeta_bc_left, const X_space& zeta_bc_right,
    const X_space& zeta_cb_left, const X_space& zeta_cb_right,
    const X_space& phi0) {

  auto [j_zeta_bc, j_zeta_cb, j_bxc, j_cxb, j_bphi0, j_cphi0] =
      compute_beta_coulomb(world, B, C, zeta_bc_left, zeta_bc_right,
                           zeta_cb_left, zeta_cb_right, phi0);
  auto [k_zeta_bc, k_zeta_cb, k_bxc, k_cxb, k_bphi0, k_cphi0] =
      compute_beta_exchange(world, B, C, zeta_bc_left, zeta_bc_right,
                            zeta_cb_left, zeta_cb_right, phi0);

  // sum k and j terms
  //
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto g_zeta_bc = 2.0 * j_zeta_bc - k_zeta_bc;
  auto g_zeta_cb = 2.0 * j_zeta_cb - k_zeta_cb;

  auto g_bxc = 2.0 * j_bxc - k_bxc;
  auto g_cxb = 2.0 * j_cxb - k_cxb;

  auto g1b = 2.0 * j_bphi0 - k_bphi0;
  auto g1c = 2.0 * j_cphi0 - k_cphi0;

  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "summing k and j terms");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }

  auto [v_bxc, v_cxb] = dipole_perturbation(world, B, C);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "Create dipole terms");
  }

  auto f_bxc = v_bxc + g_bxc;
  auto f_cxb = v_cxb + g_cxb;
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  // The first term to compute is -Q g1[K^BC], -Q g1[K^BC_conjugate]
  QProjector<double, 3> projector(world, ground_orbitals);
  auto apply_projector = [&](auto& xi) {
    return projector(xi);
  };

  g_zeta_bc = -1.0 * oop_apply(g_zeta_bc, apply_projector, true);
  g_zeta_cb = -1.0 * oop_apply(g_zeta_cb, apply_projector, true);
  f_bxc = -1.0 * oop_apply(f_bxc, apply_projector, true);
  f_cxb = -1.0 * oop_apply(f_cxb, apply_projector, true);

  g_zeta_bc.truncate();
  g_zeta_cb.truncate();
  f_bxc.truncate();
  f_cxb.truncate();

  world.gop.fence();
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "projecting terms");
  }

  auto [VBphi0, VCphi0] = dipole_perturbation(world, phi0, phi0);
  auto [FBC, FCB] = compute_first_order_fock_matrix_terms_v2(
      world, B, C, g1b, g1c, VBphi0, VCphi0, phi0);

  auto VBC = g_zeta_bc + g_zeta_cb + f_bxc + f_cxb + FBC + FCB;
  VBC.truncate();
  return VBC;
  // the next term we need to compute are the first order fock matrix terms
}

X_space QuadraticResponse::compute_second_order_perturbation_terms_v3(
    World& world, const X_space& B, const X_space& C,
    const response_space& phiBC, const response_space& phiCB,
    const vector_real_function_3d& phi0) {

  auto create_dipole = [&]() {
    vector_real_function_3d dipole_vectors(3);
    size_t i = 0;
    // creates a vector of x y z dipole functions
    for (auto& d : dipole_vectors) {
      std::vector<int> f(3, 0);
      f[i++] = 1;
      d = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
    }
    return dipole_vectors;
  };

  auto dipole_vectors = create_dipole();  // x y z
  truncate(dipole_vectors, FunctionDefaults<3>::get_thresh(), true);

  auto VBC_compare = this->VBC;

  X_space VBC(world, BC_index_pairs.size(), B.num_orbitals());
  int i = 0;
  auto num_states = BC_index_pairs.size();
  for (int i = 0; i < num_states; i++) {

    auto b = this->index_B[i];
    auto c = this->index_C[i];

    const auto& bx = B.x[b];
    const auto& by = B.y[b];
    const auto& cx = C.x[c];
    const auto& cy = C.y[c];
    const auto& phibc = phiBC[i];
    const auto& phicb = phiCB[i];
    const auto& vb = dipole_vectors[b];
    const auto& vc = dipole_vectors[c];

    std::string bc = a_directions[b] + a_directions[c];
    world.gop.fence();

    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
    }
    VBC.x[i] = compute_vbc(world, {bx, by}, {cx, cy}, {bx, cy}, {phi0, phibc},
                           phi0, vb);
    if (r_params.print_level() >= 1) {
      std::string message = "VBC.x[" + std::to_string(i) + "] BC=" + bc;
      molresponse::end_timer(world, message.c_str());
    }
    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
    }

    VBC.x[i] += compute_vbc(world, {cx, cy}, {bx, by}, {cx, by}, {phi0, phicb},
                            phi0, vc);
    if (r_params.print_level() >= 1) {
      std::string message = "VBC.x[" + std::to_string(i) + "] BC=" + bc;
      molresponse::end_timer(world, message.c_str());
    }

    // VBC is conjugate, threefore we just swap the x and y components
    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
    }
    VBC.y[i] = compute_vbc(world, {by, bx}, {cy, cx}, {cy, bx}, {phibc, phi0},
                           phi0, vb);
    if (r_params.print_level() >= 1) {
      std::string message = "VBC.y[" + std::to_string(i) + "] BC=" + bc;
      molresponse::end_timer(world, message.c_str());
    }
    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
    }
    VBC.y[i] += compute_vbc(world, {cy, cx}, {by, bx}, {by, cx}, {phicb, phi0},
                            phi0, vc);
    if (r_params.print_level() >= 1) {
      std::string message = "VBC.y[" + std::to_string(i) + "] BC=" + bc;
      molresponse::end_timer(world, message.c_str());
    }

    // auto rVBV = VBC_compare - VBC;
    //
    //
    // compare VBX.x[i] and VBX.y[i] to VBC_compare.x[i] and VBC_compare.y[i]
    //
    //
    auto vbx_norm = norm2(world, VBC.x[i]);
    auto vby_norm = norm2(world, VBC.y[i]);

    auto compare_norm = norm2(world, VBC_compare.x[i]);
    auto compare_norm_y = norm2(world, VBC_compare.y[i]);

    auto rxi = VBC_compare.x[i] - VBC.x[i];
    auto ryi = VBC_compare.y[i] - VBC.y[i];

    auto rxi_norm = norm2(world, rxi);
    auto ryi_norm = norm2(world, ryi);

    if (world.rank() == 0) {
      print("VBC.x[", i, ",", b, ",", c, "] norm: ", vbx_norm,
            " compare norm: ", compare_norm, " rxi norm: ", rxi_norm);
      print("VBC.y[", i, ",", b, ",", c, "] norm: ", vby_norm,
            " compare norm: ", compare_norm_y, " ryi norm: ", ryi_norm);
    }
  }
  return VBC;
}

X_space QuadraticResponse::compute_second_order_perturbation_terms(
    World& world, const X_space& B, const X_space& C, const X_space& zeta_bc_x,
    const X_space& zeta_bc_y, const X_space& zeta_cb_x,
    const X_space& zeta_cb_y, const X_space& phi0) {
  // The first term to compute is -Q g1[K^BC], -Q g1[K^BC_conjugate]
  QProjector<double, 3> projector(world, ground_orbitals);
  auto apply_projector = [&](auto& xi) {
    return projector(xi);
  };

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }

  // We have 2 terms to compute because we need to compute the contributions of BC and CB terms
  auto g1_kbc =
      -1.0 * oop_apply(compute_g1_term(world, zeta_bc_x, zeta_bc_y, phi0),
                       apply_projector);
  auto g1_kcb =
      -1.0 * oop_apply(compute_g1_term(world, zeta_cb_x, zeta_cb_y, phi0),
                       apply_projector);

  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "g1[zeta_bc] and g1[zeta_cb]");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  // the next term we need to compute are the first order fock matrix terms
  // -Q FB*XC
  auto g1bxc =
      -1.0 * oop_apply(compute_g1_term(world, B, phi0, C), apply_projector);
  // -Q FC*XB
  auto g1cxb =
      -1.0 * oop_apply(compute_g1_term(world, C, phi0, B), apply_projector);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "g1[B]C and g1[C]B");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  // -Q ( FB ) * ( XC )
  auto [vbxc, vcxb] = dipole_perturbation(world, C, B);

  vbxc = -1.0 * oop_apply(vbxc, apply_projector, false);
  vcxb = -1.0 * oop_apply(vcxb, apply_projector, false);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "V(B)C and VA(C)B");
  }

  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  auto [zFBzC, zFCzB] =
      compute_first_order_fock_matrix_terms(world, B, phi0, C);
  if (r_params.print_level() >= 1) {
    molresponse::end_timer(world, "FB[i,j]*C and FC[i,j]*B");
  }

  return g1_kbc + g1_kcb + g1bxc + g1cxb + zFBzC + zFCzB + vbxc + vcxb;

  // the next term we need to compute are the first order fock matrix terms
}

/**
 * @ brief computes the occ/occ and virt/virt blocks of second order gamma which leads to 4 space objects
 * ABX,ABY,BAX, and BAY where ABX contains A.X and phitildeab and ABY contains B.Y and phi0
 *
 *
 * @param world
 * @param B
 * @param C
 * @return
 */

std::tuple<X_space, X_space, X_space, X_space>
QuadraticResponse::compute_zeta_response_vectors(World& world, const X_space& B,
                                                 const X_space& C) {

  // zeta_bc_x=[x(b),tilde_phi_dagger(bc)]
  // zeta_bc_y=[y_dagger(c),phi_0]
  // zeta_cb_x=[x(c),tilde_phi_dagger(cb)]
  // zeta_cb_y=[y_dagger(b),phi_0]

  X_space zeta_bc_left = X_space(world, B.num_states(), C.num_orbitals());
  X_space zeta_bc_right = X_space(world, B.num_states(), C.num_orbitals());

  X_space zeta_cb_left = X_space(world, C.num_states(), B.num_orbitals());
  X_space zeta_cb_right = X_space(world, C.num_states(), B.num_orbitals());

  // Here are all the x components
  for (auto i = 0; i < B.num_states(); i++) {
    zeta_bc_left.x[i] = copy(world, B.x[i], false);
    zeta_cb_left.x[i] = copy(world, C.x[i], false);

    zeta_bc_right.x[i] = copy(world, C.y[i], false);
    zeta_cb_right.x[i] = copy(world, B.y[i], false);

    zeta_cb_right.y[i] = copy(world, ground_orbitals, false);
    zeta_bc_right.y[i] = copy(world, ground_orbitals, false);
  }
  world.gop.fence();

  // Here contains the y components
  // construct tilde_phi_dagger(bc) and tilde_phi_dagger(cb)
  // where tilde_phi_dagger(bc)_i = \sum_j <y_i(b) | x_j(c)> * phi0_j
  for (auto i = 0; i < B.num_states(); i++) {

    auto matrix_bycx = matrix_inner(world, B.y[i], C.x[i]);
    auto tilde_phi_bc =
        -1.0 * transform(world, ground_orbitals, 1.0 * matrix_bycx, true);
    zeta_bc_left.y[i] = copy(world, tilde_phi_bc, true);

    auto matrix_cybx = matrix_inner(world, C.y[i], B.x[i]);
    auto tilde_phi_cb =
        -1.0 * transform(world, ground_orbitals, 1.0 * matrix_cybx, true);
    zeta_cb_left.y[i] = copy(world, tilde_phi_cb, true);

    // print the norms of each vector
  }

  zeta_bc_left.truncate();
  zeta_bc_right.truncate();
  zeta_cb_left.truncate();
  zeta_cb_right.truncate();

  return {zeta_bc_left, zeta_bc_right, zeta_cb_left, zeta_cb_right};
}

auto QuadraticResponse::dipole_perturbation(World& world, const X_space& left,
                                            const X_space& right) const
    -> std::pair<X_space, X_space> {

  auto num_states = left.num_states();
  MADNESS_ASSERT(num_states == right.num_states());

  auto VB = X_space(world, left.num_states(), right.num_orbitals());
  auto VC = X_space(world, left.num_states(), right.num_orbitals());

  vector_real_function_3d dipole_vectors(3);
  size_t i = 0;
  // creates a vector of x y z dipole functions
  for (auto& d : dipole_vectors) {
    std::vector<int> f(3, 0);
    f[i++] = 1;
    d = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
  }
  truncate(world, dipole_vectors, true);

  for (int i = 0; i < num_states; i++) {

    VB.x[i] = mul(world, dipole_vectors[index_B[i]], right.x[i], true);
    VB.y[i] = mul(world, dipole_vectors[index_B[i]], right.y[i], true);

    VC.x[i] = mul(world, dipole_vectors[index_C[i]], left.x[i], true);
    VC.y[i] = mul(world, dipole_vectors[index_C[i]], left.y[i], true);
  }
  world.gop.fence();

  VB.truncate();
  VC.truncate();

  return {VB, VC};
}

// compute the g1 term
// we have 3 x_space which hold x and y components of the density matrix.
//
X_space QuadraticResponse::compute_g1_term(World& world, const X_space& left,
                                           const X_space& right,
                                           const X_space& apply) const {

  auto JBC = compute_coulomb_term(world, left, right, apply);
  auto KBC = compute_exchange_term(world, left, right, apply);

  return 2 * JBC - KBC;
}
// compute rhoBC and apply to D
X_space QuadraticResponse::compute_coulomb_term(World& world, const X_space& B,
                                                const X_space& C,
                                                const X_space& x_apply) const {

  X_space J = X_space::zero_functions(world, B.num_states(), B.num_orbitals());
  vector_real_function_3d rhoX(B.num_states());

  vector_real_function_3d x_phi, y_phi;

  // create the density for each state in B
  //B.x[j] and B.y[j] are vectors fo functions
  for (const auto& j : B.active) {
    x_phi = mul(world, B.x[j], C.x[j], true);
    y_phi = mul(world, B.y[j], C.y[j], true);

    rhoX[j] = sum(world, x_phi, true);
    rhoX[j] += sum(world, y_phi, true);
  }
  auto temp_J = apply(world, *shared_coulomb_operator, rhoX);

  for (const auto& j : B.active) {
    J.x[j] = mul(world, temp_J[j], x_apply.x[j], true);
    J.y[j] = mul(world, temp_J[j], x_apply.y[j], true);
  }

  // truncate(world, rhoX);
  J.truncate();

  return J;
}

X_space QuadraticResponse::compute_exchange_term(World& world, const X_space& B,
                                                 const X_space& C,
                                                 const X_space& x_apply) const {

  auto make_operator = [&](const vecfuncT& ket, const vecfuncT& bra) {
    const double lo = 1.e-10;
    auto& world = ket[0].world();
    Exchange<double, 3> k{world, lo};
    k.set_bra_and_ket(bra, ket);

    std::string algorithm_ = r_params.hfexalg();

    if (algorithm_ == "multiworld") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient);
    } else if (algorithm_ == "multiworld_row") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::multiworld_efficient_row);
    } else if (algorithm_ == "largemem") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::large_memory);
    } else if (algorithm_ == "smallmem") {
      k.set_algorithm(Exchange<double, 3>::Algorithm::small_memory);
    }

    return k;
  };

  // if the frequecy of A is 0 we run the static case
  // else we run the dynamic case
  auto K = X_space::zero_functions(world, B.num_states(), B.num_orbitals());

  vector_real_function_3d xb;
  vector_real_function_3d yb;
  // compute_x
  for (int k = 0; k < B.num_states(); k++) {

    auto K1 = make_operator(B.x[k], C.x[k]);
    auto k1 = K1(x_apply.x[k]);
    auto K2 = make_operator(C.y[k], B.y[k]);
    auto k2 = K2(x_apply.x[k]);
    K.x[k] = gaxpy_oop(1.0, k1, 1.0, k2, true);
  }
  //compute_y
  for (int k = 0; k < B.num_states(); k++) {
    auto K1_conjugate = make_operator(B.y[k], C.y[k]);
    auto k1_c = K1_conjugate(x_apply.y[k]);
    auto K2_conjugate = make_operator(C.x[k], B.x[k]);
    auto k2_c = K2_conjugate(x_apply.y[k]);
    K.y[k] = gaxpy_oop(1.0, k1_c, 1.0, k2_c, true);
  }

  K.truncate();

  return K;
}

//
void PODResponse::compute_pod_modes(World& world) {
  //
  // 1. Convert the x_data into a response matrix
  // 2. Compute the SVD of the response matrix
  // 3. Compute the POD modes

  // num calculations
  auto num_states = x_data.size() * 3;
  auto m = x_data.size();
  auto n = x_data[0].first.num_orbitals();

  auto all_chi = X_space(world, num_states, n);
  if (world.rank() == 0)
    print("m: ", m);
  if (world.rank() == 0)
    print("n: ", n);

  for (auto i = 0; i < m; i++) {

    // first step is to copy x_data_i into response_matrix_i
    auto x_data_i = x_data[i].first;
    for (auto j = 0; j < 3; j++) {
      auto index = i * 3 + j;
      if (world.rank() == 0)
        print("index: ", index);
      all_chi.x[index] = copy(world, x_data_i.x[j]);
      all_chi.y[index] = copy(world, x_data_i.y[j]);
    }
  }

  auto A = inner(all_chi, all_chi);

  // eigen values and eigen vectors

  // literally the same thing as the svd...
  Tensor<double> e(static_cast<long>(num_states), 1);
  Tensor<double> u(static_cast<long>(num_states),
                   static_cast<long>(num_states));
  syev(A, u, e);
  // A = U * diag(s) * VT    for A real
  // A = U * diag(s) * VH    for A complex
  Tensor<double> VT(static_cast<long>(num_states),
                    static_cast<long>(num_states));
  Tensor<double> U(static_cast<long>(num_states),
                   static_cast<long>(num_states));
  Tensor<double> sigma(static_cast<long>(num_states));
  svd(A, U, sigma, VT);
  // create a json printing out the singular values
  json j = {};
  j["sigma"] = tensor_to_json(sigma);
  j["VT"] = tensor_to_json(VT);
  j["U"] = tensor_to_json(U);
  j["A"] = tensor_to_json(A);
  j["eigenvectors"] = tensor_to_json(u);
  j["eigenvalues"] = tensor_to_json(e);

  if (world.rank() == 0)
    print("sigma: ", sigma);
  if (world.rank() == 0)
    print("VT: ", VT);
  if (world.rank() == 0)
    print("U: ", U);
  if (world.rank() == 0)
    print("A: ", A);

  // compute the POD modes

  if (world.rank() == 0) {
    std::ofstream o("pod.json");
    o << std::setw(4) << j << std::endl;
  }
}
void PODResponse::compute_pod_modes_2(World& world) {
  //
  // 1. Convert the x_data into a response matrix
  // 2. Compute the SVD of the response matrix
  // 3. Compute the POD modes
  // In this strategy we will unpack all x and y response functions into a single vector.

  // num calculations
  auto num_states = x_data.size() * 3;
  auto m = x_data.size();
  auto n = x_data[0].first.num_orbitals();
  auto total_orbitals = num_states * n * 2;  // 2 for x and y

  if (world.rank() == 0)
    print("m: ", m);
  if (world.rank() == 0)
    print("n: ", n);

  auto all_chi = X_space(world, num_states, n);
  if (world.rank() == 0)
    print("m: ", m);
  if (world.rank() == 0)
    print("n: ", n);

  // Create one large x_space with each state
  for (auto i = 0; i < m; i++) {
    auto x_data_i = x_data[i].first;
    for (auto j = 0; j < 3; j++) {
      auto index = i * 3 + j;
      all_chi.x[index] = copy(world, x_data_i.x[j]);
      all_chi.y[index] = copy(world, x_data_i.y[j]);
    }
  }

  // Turn into response matrix
  auto all_x = to_response_matrix(all_chi);
  vector_real_function_3d all_response_functions(total_orbitals);
  // unpack all_x into a single vector
  auto index = 0;
  for (const auto& all_x_i : all_x) {
    for (const auto& all_x_ij : all_x_i) {
      all_response_functions[index++] = madness::copy(all_x_ij, false);
    }
  }
  world.gop.fence();

  // Compute the 2 pt correlation matrix between response functions  // contains both x and y components
  auto A =
      matrix_inner(world, all_response_functions, all_response_functions, true);

  // auto A = inner(all_chi, all_chi);

  // eigen values and eigen vectors

  // literally the same thing as the svd...
  Tensor<double> e(static_cast<long>(num_states), 1);
  Tensor<double> u(static_cast<long>(num_states),
                   static_cast<long>(num_states));
  syev(A, u, e);

  Tensor<double> UU(static_cast<long>(total_orbitals),
                    static_cast<long>(total_orbitals));
  Tensor<double> EE(static_cast<long>(total_orbitals), 1);

  for (auto i = 0; i < total_orbitals; i++) {
    auto index_i = total_orbitals - i - 1;
    for (auto j = 0; j < total_orbitals; j++) {
      UU(i, j) = u(static_cast<long>(index_i), j);
    }
    EE(i, 0) = e(static_cast<long>(index_i), 0);
  }
  // A = U * diag(s) * VT    for A real
  // A = U * diag(s) * VH    for A complex
  Tensor<double> VT(static_cast<long>(num_states),
                    static_cast<long>(num_states));
  Tensor<double> U(static_cast<long>(num_states),
                   static_cast<long>(num_states));
  Tensor<double> sigma(static_cast<long>(num_states));
  svd(A, U, sigma, VT);
  // create a json printing out the singular values
  if (world.rank() == 0)
    print("Eigenvlaues: \n", EE);

  auto pod_modes = transform(world, all_response_functions, UU, false);
  world.gop.fence();

  double pod_threshold = 1e-6;
  // find index of first sigma less than threshold
  auto pod_index = 0;
  for (auto i = 0; i < sigma.size(); i++) {
    if (sigma[i] < pod_threshold) {
      pod_index = i;
      break;
    }
  }
  if (world.rank() == 0)
    print("pod_index: ", pod_index);
  // select the last pod_index modes
  auto num_modes = pod_index + 1;
  auto pod_modes_selected = vector_real_function_3d(num_modes);

  for (auto i = 0; i < num_modes; i++) {
    pod_modes_selected[i] = copy(pod_modes[i], false) * (1 / EE[i]);
  }
  // normalize the pod modes
  normalize(world, pod_modes_selected, true);
  // auto etas = orthonormalize(pod_modes_selected);
  auto etas = pod_modes_selected;
  world.gop.fence();

  // now compute the orbital energies of the pod modes
  auto kinetic_energy = vector<double>(num_modes);

  vecfuncT detax = apply(world, *(gradop[0]), etas, false);
  vecfuncT detay = apply(world, *(gradop[1]), etas, false);
  vecfuncT detaz = apply(world, *(gradop[2]), etas, false);
  world.gop.fence();

  auto dx2 = mul(world, detax, detax, false);
  auto dy2 = mul(world, detay, detay, false);
  auto dz2 = mul(world, detaz, detay, false);
  world.gop.fence();
  for (auto i = 0; i < num_modes; i++) {
    kinetic_energy[i] = 0.5 * inner(detax[i], detax[i]) +
                        0.5 * inner(detay[i], detay[i]) +
                        0.5 * inner(detaz[i], detaz[i]);
  }
  world.gop.fence();
  // print the kinetic energy of each mode
  if (world.rank() == 0)
    print("kinetic energy: \n", kinetic_energy);

  // now get the nuclear energy for the system
  auto nuclear_energy = vector<double>(num_modes);
  auto vnuc = this->potential_manager->vnuclear();
  auto vk = vnuc * square(world, etas);
  for (auto i = 0; i < num_modes; i++) {
    nuclear_energy[i] = inner(etas[i], vk[i]);
  }
  // print the nuclear energy of each mode
  if (world.rank() == 0)
    print("nuclear energy: \n", nuclear_energy);
  auto ground_orbitals = this->get_orbitals();
  // compute ground density
  auto ground_density = 2.0 * sum(world, square(world, ground_orbitals), true);
  auto coulomb_j = apply(*shared_coulomb_operator, ground_density);
  auto coulomb_energy = vector<double>(num_modes);
  auto jk = coulomb_j * square(world, etas);
  for (auto i = 0; i < num_modes; i++) {
    coulomb_energy[i] = inner(etas[i], jk[i]);
  }
  // print the coulomb energy of each mode
  if (world.rank() == 0)
    print("coulomb energy: \n", coulomb_energy);

  // exchange
  auto exchange_energy = vector<double>(num_modes);
  for (int a = 0; a < num_modes; a++) {
    auto phi_aj = copy(world, ground_orbitals, true);
    phi_aj = mul(world, etas[a], phi_aj, true);
    auto exchange_j = apply(world, *shared_coulomb_operator, phi_aj);
    auto exchange_jk = sum(world, exchange_j, true);

    exchange_energy[a] = -2 * inner(etas[a], exchange_jk);
  }
  // print the exchange energy of each mode
  if (world.rank() == 0)
    print("exchange energy: \n", exchange_energy);

  // compute the total energy of each mode
  auto total_energy = vector<double>(num_modes);
  for (auto i = 0; i < num_modes; i++) {
    total_energy[i] = kinetic_energy[i] + nuclear_energy[i] +
                      coulomb_energy[i] + exchange_energy[i];
  }
  // print the total energy of each mode
  if (world.rank() == 0)
    print("total energy: \n", total_energy);
  // Now sort the modes by total energy and only grab the positive modes

  std::vector<int> indices(etas.size());
  std::iota(indices.begin(), indices.end(),
            0);  // Fill indices with 0, 1, ..., n

  std::sort(indices.begin(), indices.end(),
            [&](int i, int j) { return total_energy[i] > total_energy[j]; });
  auto sorted_modes = vector_real_function_3d(num_modes);
  auto sorted_energies = vector<double>(num_modes);

  for (auto i = 0; i < num_modes; i++) {
    sorted_modes[i] = copy(etas[indices[i]], false);
    sorted_energies[i] = total_energy[indices[i]];
  }
  world.gop.fence();
  // Take only the positive
  auto pos_modes = vector_real_function_3d();
  auto pos_energies = vector<double>();
  for (auto i = 0; i < num_modes; i++) {
    // if the sorted energy is positive push back
    if (sorted_energies[i] > 0) {
      pos_modes.push_back(copy(sorted_modes[i], false));
      pos_energies.push_back(sorted_energies[i]);
    }
  }
  world.gop.fence();
  auto num_pos_modes = pos_modes.size();
  if (world.rank() == 0)
    print("num_pos_modes: ", num_pos_modes);
  // print the positive energies
  if (world.rank() == 0)
    print("pos_energies: \n", pos_energies);

  // now remove the positive modes

  json j = {};
  j["sigma"] = tensor_to_json(sigma);
  j["VT"] = tensor_to_json(VT);
  j["U"] = tensor_to_json(U);
  j["A"] = tensor_to_json(A);
  j["eigenvectors"] = tensor_to_json(UU);
  j["eigenvalues"] = tensor_to_json(EE);
  j["kinetic_energy"] = kinetic_energy;
  j["nuclear_energy"] = nuclear_energy;
  j["coulomb_energy"] = coulomb_energy;
  j["exchange_energy"] = exchange_energy;
  j["total_energy"] = pos_energies;

  if (world.rank() == 0) {
    std::ofstream o("pod.json");
    o << std::setw(4) << j << std::endl;
  }

  world.gop.fence();
}
void PODResponse::load(World& world, const std::string& name) {}
void PODResponse::save(World& world, const std::string& name) {}

void PODResponse::iterate(World& world) {}
void PODResponse::initialize(World& world) {}
