//
// Created by adrianhurtado on 1/24/22.
//

#include "densityVector.hpp"
#include "../timer.h"

namespace madness{

}

response_space FrequencySolver::PropertyRHS(World &world, PropertyBase &p) const {
  if (r_params.print_level() >= 1) {
    molresponse::start_timer(world);
  }
  response_space rhs(world, p.num_operators, r_params.num_orbitals());

  reconstruct(world, ground_orbitals);
  QProjector<double, 3> Qhat(world, ground_orbitals);
  // Set the dipoles (ground orbitals are probably
  // more accurate now, so recalc the dipoles)
  // why is it called dipole guess.
  // This is just orbitals times dipole operator
  std::vector<real_function_3d> orbitals = ground_orbitals;

  print("num operators ", p.num_operators);
  for (size_t i = 0; i < p.num_operators; i++) {
    // question here....MolecularDerivativeFunctor takes derivative with
    // respect to axis atom and axis
    // here we save
    // need to project

    rhs[i] = mul(world, p.operator_vector.at(i), ground_orbitals, r_params.lo());

    truncate(world, rhs[i]);
    // rhs[i].truncate_vec();

    // project rhs vectors for state
    rhs[i] = Qhat(rhs[i]);
    // truncate(world, rhs[i], true);
    for (size_t j = 0; j < orbitals.size(); j++) {
      print("RHS norm for after orbital ", j, "Response state  ", i, ": ", rhs[i][j].norm2());
    }

    world.gop.fence();
    // core projector contribution
  }

  // if (world.rank() ==dipole 0) print("derivatives:\n", r, ru, rc, ra);
  molresponse::end_timer(world, "rhs vectors");
  return rhs;
}
void FrequencySolver::compute(World &world, FrequencyVector& rho) {
  molresponse::start_timer(world);
  if (r_params.plot_initial()) {
    if (world.rank() == 0) print("\n   Plotting ground state densities.\n");
    if (r_params.plot_l() > 0.0)
      do_vtk_plots(world,
                   r_params.plot_pts(),
                   r_params.plot_l(),
                   0,
                   r_params.num_orbitals(),
                   molecule,
                   square(world, ground_orbitals),
                   "ground");
    else
      do_vtk_plots(world,
                   r_params.plot_pts(),
                   r_params.L() / 2.0,
                   0,
                   r_params.num_orbitals(),
                   molecule,
                   square(world, ground_orbitals),
                   "ground");
  }
  // Warm and fuzzy
  if (world.rank() == 0) {
    print("\n\n    Response Calculation");
    print("   ------------------------");
  }
  auto protocols=r_params.protocol();
  print(protocols);
  for (unsigned int proto = 0; proto < protocols.size(); proto++) {
    // Set defaults inside here
    // default value of
    auto thresh=protocols[proto];
    print(thresh);
    set_protocol(world,thresh);

    check_k(world, thresh, FunctionDefaults<3>::get_k());
    // Do something to ensure all functions have same k value
    j_molresponse["protocol_data"].push_back({});
    j_molresponse["protocol_data"][proto]["proto"] = thresh;//r_params.protocol()[proto];
    j_molresponse["protocol_data"][proto]["iter_data"]={};

    if (r_params.dipole()) {
      if (world.rank() == 0) print("creating dipole property operator");
       p = DipoleVector(world);
    } else if (r_params.nuclear()) {
      if (world.rank() == 0) print("creating nuclear property operator");
       p = NuclearVector(world, molecule);
    }
    if (proto == 0) {
      if (r_params.restart()) {
        if (world.rank() == 0) print("   Initial guess from file:", r_params.restart_file());
        load(world, r_params.restart_file());
        check_k(world, thresh, FunctionDefaults<3>::get_k());

        if (r_params.dipole()) {
          // set states
          PQ.X = PropertyRHS(world, p);
          PQ.Y = PQ.X.copy();
          print("P: ", PQ.X.norm2());
          print("Q: ", PQ.Y.norm2());

          // set RHS_Vector
        } else if (r_params.nuclear()) {
          PQ.X = PropertyRHS(world, p);
          PQ.Y = PQ.X.copy();

        } else if (r_params.order2()) {
          //
        } else if (r_params.order3()) {
          //
        } else {
          MADNESS_EXCEPTION("not a valid response state ", 0);
        }
      } else {  // Dipole guesses

        if (r_params.dipole()) {
          PQ.X = PropertyRHS(world, p);
          PQ.Y = PQ.X.copy();
          // set states
          //
          print("okay this is not a good idea if it comes up more than once");
          // set RHS_Vector
        } else if (r_params.nuclear()) {
          PQ.X = PropertyRHS(world, p);
          PQ.Y = PQ.X.copy();
        } else if (r_params.order2()) {
          //
        } else if (r_params.order3()) {
          //
        } else {
          MADNESS_EXCEPTION("not a valid response state ", 0);
        }
      }
    }
    //
    // Here i should print some information about the calculation we are
    // about to do
    print("Pre iteration Information");
    print("Number of Response States: ", r_params.n_states());
    print("Number of Ground States: ", r_params.num_orbitals());
    print("k = ", FunctionDefaults<3>::get_k());
    print("protocol threshold = ", FunctionDefaults<3>::get_k());

    print("Property rhs func k = ", PQ.X[0][0].k());
    print("Property func k thresh= ", PQ.X[0][0].thresh());

    print("Property rhs func Q k = ", PQ.Y[0][0].k());
    print("Property func Q k thresh = ", PQ.Y[0][0].thresh());

    print("Property rhs func P norms", PQ.X.norm2());
    print("Property rhs func Q norms", PQ.Y.norm2());

    if (proto > 0) {
      r_params.set_derived_value<bool>("first_run", false);
    }
    // Now actually ready to iterate...
    // iterate_freq(world);
    iterate(world);
    // IterateFrequencyResponse(world, P, Q);
  }  // end for --finished reponse density

  // Plot the response function if desired
  if (r_params.plot()) {
    // Need to get densities first
    std::vector<real_function_3d> densities = make_density(world);

    // For the instance where we don't plot all the orbitals
    std::vector<real_function_3d> plot_densities;
    for (size_t i : r_params.plot_data()) {
      plot_densities.push_back(densities[i]);
    }

    // Now plot
    if (world.rank() == 0) print("\n   Plotting response state densities.\n");
    if (r_params.plot_l() > 0.0)
      do_vtk_plots(world,
                   r_params.plot_pts(),
                   r_params.plot_l(),
                   0,
                   r_params.plot_data().size(),
                   molecule,
                   plot_densities,
                   "response-state");
    else
      do_vtk_plots(world,
                   r_params.plot_pts(),
                   r_params.L(),
                   0,
                   r_params.plot_data().size(),
                   molecule,
                   plot_densities,
                   "response-state");
  }

  std::cout.precision(2);
  std::cout << std::fixed;
  // Get start time
  molresponse::end_timer(world, "total:");
}  // end compute frequency response

void FrequencySolver::iterate(World& world) {
  size_t iter;
  // Variables needed to iterate
  QProjector<double, 3> projector(world, ground_orbitals);
  size_t m = r_params.n_states();      // Number of excited states
  size_t n = r_params.num_orbitals();  // Number of ground state orbitals

  real_function_3d v_xc;   // For TDDFT
  bool converged = false;  // Converged flag
  const double dconv = std::max(FunctionDefaults<3>::get_thresh(), r_params.dconv());
  // m residuals for x and y
  Tensor<double> bsh_residualsX(m);
  Tensor<double> bsh_residualsY(m);
  Tensor<double> density_residuals(m);

  vecfuncT rho_omega_old(m);

  // initialize DFT XC functional operator
  XCOperator<double, 3> xc = create_XCOperator(world,  r_params.xc());

  // create X space residuals
  X_space residuals(world, m, n);
  X_space old_Chi(world, m, n);
  // Create the X space
  // vector of Xvectors
  std::vector<X_vector> Xvector;
  std::vector<X_vector> Xresidual;
  for (size_t b = 0; b < m; b++) {
    Xvector.emplace_back(Chi, b);
    Xresidual.emplace_back(residuals, b);
  }
  // If DFT, initialize the XCOperator<double,3>

  // create a std vector of XNONLinearsolvers
  NonLinearXsolver kain_x_space;
  for (size_t b = 0; b < m; b++) {
    kain_x_space.push_back(
        XNonlinearSolver<X_vector, double, X_space_allocator>(X_space_allocator(world, n), true));
  }
  // if kain set the max sub of each of the solvers
  if (r_params.kain()){
    for (auto&  solver:kain_x_space){
      solver.set_maxsub(r_params.maxsub());
    }
  }
  //
  double omega_n = r_params.omega();
  omega_n = abs(omega_n);
  omega[0] = omega_n;
  // We compute with positive frequencies
  print("Warning input frequency is assumed to be positive");
  print("Computing at positive frequency omega = ", omega_n);
  double x_shifts = 0.0;
  double y_shifts = 0.0;
  // if less negative orbital energy + frequency is positive or greater than 0
  if ((ground_energies[n - 1] + omega_n) >= 0.0) {
    // Calculate minimum shift needed such that \eps + \omega + shift < 0
    print("*** we are shifting just so you know!!!");
    x_shifts = -.05 - (omega_n + ground_energies[n - 1]);
  }
  std::vector<poperatorT> bsh_x_ops = make_bsh_operators_response(world, x_shifts, omega_n);
  std::vector<poperatorT> bsh_y_ops;

  bool static_res = (omega_n == 0.0);
  bool compute_y = not static_res;
  // Negate omega to make this next set of BSH operators \eps - omega
  if (compute_y) {
    omega_n = -omega_n;
    bsh_y_ops = make_bsh_operators_response(world, y_shifts, omega_n);
    omega_n = -omega_n;
  }

  Tensor<double> maxrotn(m);
  maxrotn.fill(dconv * 100);

  for (iter = 0; iter <= r_params.maxiter(); ++iter) {
    // Basic output
    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
      if (world.rank() == 0) printf("\n   Iteration %d at time %.1fs\n", static_cast<int>(iter), wall_time());
      if (world.rank() == 0) print("-------------------------------------------");
    }
    if (r_params.print_level() >= 1) {
      if (world.rank() == 0) {
        print("Chi.x norms at start of iteration: ", iter);
        print(Chi.X.norm2());
        print("Chi.y norms at start of iteration ", iter);
        print(Chi.Y.norm2());
      }
    }

    old_Chi = Chi.copy();
    rho_omega_old = rho_omega;
    // compute rho_omega
    rho_omega = transition_density(world, ground_orbitals, Chi.X, Chi.Y);
    // rho_omega = make_density(world, Chi, compute_y);

    if (iter < 2 || (iter % 10) == 0) {
      load_balance(world);
    }

    // compute density residuals
    if (iter > 0) {
      density_residuals = norm2s_T(world, (rho_omega - rho_omega_old));
      // Take the max between this an a minimum maxrotn step
      maxrotn = (bsh_residualsX + bsh_residualsY) / 4;
      print("maxrotn", maxrotn);
      for (size_t i = 0; i < Chi.num_states(); i++) {
        if (maxrotn[i] < r_params.maxrotn()) {
          maxrotn[i] = r_params.maxrotn();
          print("less than maxrotn....set to maxrotn");
        }
      }
      if (world.rank() == 0 and (r_params.print_level() > 1)) {
        print("Density residuals");
        print(density_residuals);
        print("BSH  residuals");
        print("x", bsh_residualsX);
        print("y", bsh_residualsY);
        print("maxrotn", maxrotn);
      }
    }

    if (iter > 0) {
      if (density_residuals.max() > 2) {
        break;
      }
      double d_residual = density_residuals.max();
      double d_conv = dconv * std::max(size_t(5), molecule.natom());
      // Test convergence and set to true
      if ((d_residual < d_conv) and
          ((std::max(bsh_residualsX.absmax(), bsh_residualsY.absmax()) < d_conv * 5.0) or
           r_params.get<bool>("conv_only_dens"))) {
        converged = true;
      }

      if (converged || iter == r_params.maxiter() - 1) {
        // if converged print converged
        if (world.rank() == 0 && converged and (r_params.print_level() > 1)) {
          print("\nConverged!\n");
        }

        if (r_params.save()) {
          molresponse::start_timer(world);
          save(world, r_params.save_file());
          if (r_params.print_level() >= 1) molresponse::end_timer(world, "Save:");
        }
        // Basic output
        if (r_params.print_level() >= 1) molresponse::end_timer(world, " This iteration:");
        // plot orbitals
        if (r_params.plot_all_orbitals()) {
          PlotGroundandResponseOrbitals(world, iter, Chi.X, Chi.Y, r_params, ground_calc);
        }
        rho0 = make_ground_density(world, ground_orbitals);
        if (r_params.plot()) {
          do_vtk_plots(world, 200, r_params.L(), molecule, rho0, rho_omega, ground_orbitals, Chi);
        }
        break;
      }
    }
    update_x_space_response(world,
                            Chi,
                            residuals,
                            xc,
                            bsh_x_ops,
                            bsh_y_ops,
                            projector,
                            x_shifts,
                            omega_n,
                            kain_x_space,
                            Xvector,
                            Xresidual,
                            bsh_residualsX,
                            bsh_residualsY,
                            iter,
                            maxrotn);

    Tensor<double> polar = -2 * inner(Chi, PQ);

    frequency_to_json(j_molresponse, iter, bsh_residualsX, bsh_residualsY, density_residuals, polar);
  }

  if (world.rank() == 0) print("\n");
  if (world.rank() == 0) print("   Finished Response Calculation ");
  if (world.rank() == 0) print("   ------------------------");
  if (world.rank() == 0) print("\n");

  // Did we converge?
  if (iter == r_params.maxiter() && not converged) {
    if (world.rank() == 0) print("   Failed to converge. Reason:");
    if (world.rank() == 0) print("\n  ***  Ran out of iterations  ***\n");
    if (world.rank() == 0) print("    Running analysis on current values.\n");
  }
  if (world.rank() == 0) {
    print(" Final energy residuals X:");
    print(bsh_residualsX);
    print(" Final energy residuals Y:");
    print(bsh_residualsY);
    print(" Final density residuals:");
    print(density_residuals);
    compute_and_print_polarizability(world, "Converged");
  }
}

void FrequencySolver::compute_and_print_polarizability(World &world, std::string message) {
  Tensor<double> G = -2 * inner(Chi, PQ);
  if (world.rank() == 0) {
    print("Polarizability", message);
    print(G);
  }
}
void FrequencySolver::frequency_to_json(json& j_mol_in,
                              size_t iter,
                              const Tensor<double>& res_X,
                              const Tensor<double>& res_Y,
                              const Tensor<double>& density_res,
                              const Tensor<double>& omega) {
  json j = {};
  j["iter"] = iter;
  j["res_X"] = tensor_to_json(res_X);
  j["res_Y"] = tensor_to_json(res_Y);
  j["density_residuals"] = tensor_to_json(density_res);
  j["polar"] = tensor_to_json(omega);
  auto index = j_mol_in["protocol_data"].size() - 1;
  j_mol_in["protocol_data"][index]["iter_data"].push_back(j);
}
void FrequencySolver::update_x_space_response(World &world,
                                    X_space &Chi,
                                    X_space &res,
                                    XCOperator<double, 3> &xc,
                                    std::vector<poperatorT> &bsh_x_ops,
                                    std::vector<poperatorT> &bsh_y_ops,
                                    QProjector<double, 3> &projector,
                                    double &x_shifts,
                                    double &omega_n,
                                    NonLinearXsolver &kain_x_space,
                                    std::vector<X_vector> Xvector,
                                    std::vector<X_vector> Xresidual,
                                    Tensor<double> &bsh_residualsX,
                                    Tensor<double> &bsh_residualsY,
                                    size_t iteration,
                                    Tensor<double> &maxrotn) {
  size_t m = Chi.num_states();
  bool compute_y = r_params.omega() != 0.0;
  // size_t n = Chi.num_orbitals();

  Tensor<double> errX(m);
  Tensor<double> errY(m);

  X_space theta_X = Compute_Theta_X(world, Chi, xc, r_params.calc_type());
  // compute residual X_space
  print("BSH update iter = ", iteration);

  X_space temp = bsh_update_response(world, theta_X, bsh_x_ops, bsh_y_ops, projector, x_shifts);

  res = compute_residual(world, Chi, temp, bsh_residualsX, bsh_residualsY, r_params.calc_type());

  // kain update with temp adjusts temp
  if (r_params.kain() && (iteration > 0)) {
    temp = kain_x_space_update(world, Chi, res, kain_x_space, Xvector, Xresidual);
  }

  if (iteration > 0) {
    x_space_step_restriction(world, Chi, temp, compute_y, maxrotn);
  }
  // truncate x
  temp.X.truncate_rf();
  // truncate y if compute y
  if (compute_y) temp.Y.truncate_rf();
  //	if not compute y then copy x in to y
  if (!compute_y) temp.Y = temp.X.copy();
  Chi = temp.copy();
  // print x norms
}
X_space FrequencySolver::bsh_update_response(World &world,
                                   X_space &theta_X,
                                   std::vector<poperatorT> &bsh_x_ops,
                                   std::vector<poperatorT> &bsh_y_ops,
                                   QProjector<double, 3> &projector,
                                   double &x_shifts) {
  size_t m = theta_X.X.size();
  size_t n = theta_X.X.size_orbitals();
  bool compute_y = r_params.omega() != 0.0;

  molresponse::start_timer(world);

  theta_X.X += Chi.X * x_shifts;
  theta_X.X += PQ.X;
  theta_X.X = theta_X.X * -2;
  theta_X.X.truncate_rf();

  if (compute_y) {
    theta_X.Y += PQ.Y;
    theta_X.Y = theta_X.Y * -2;
    theta_X.Y.truncate_rf();
  }
  molresponse::end_timer(world, "Compute residual stuff theta_X");

  // apply bsh
  molresponse::start_timer(world);
  X_space bsh_X(world, m, n);

  bsh_X.X = apply(world, bsh_x_ops, theta_X.X);
  if (compute_y) {
    bsh_X.Y = apply(world, bsh_y_ops, theta_X.Y);
  }
  molresponse::end_timer(world, "Apply BSH to theta_X");

  molresponse::start_timer(world);
  // Project out ground state
  for (size_t i = 0; i < m; i++) bsh_X.X[i] = projector(bsh_X.X[i]);

  if (compute_y) {
    for (size_t i = 0; i < m; i++) {
      bsh_X.Y[i] = projector(bsh_X.Y[i]);
    }
    bsh_X.truncate();
  } else {
    bsh_X.X.truncate_rf();
    bsh_X.Y = bsh_X.X.copy();
  }
  molresponse::end_timer(world, "Project and truncate BSH_X");

  return bsh_X;
}



