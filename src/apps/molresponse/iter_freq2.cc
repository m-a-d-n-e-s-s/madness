
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
#include "apps/external_headers/tensor_json.hpp"
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "madness/mra/funcdefaults.h"
#include "molresponse/basic_operators.h"
#include "molresponse/density.h"
#include "molresponse/global_functions.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/timer.h"

json freq_iteration_to_json(int iter,
                            const Tensor<double>& bsh_residualsX,
                            const Tensor<double>& bsh_residualsY,
                            const Tensor<double>& density_residuals,
                            const Tensor<double>& polar) {
  json j_iter = {};
  j_iter["iter"] = iter;
  j_iter["bsh_residualsX"] = tensor_to_json(bsh_residualsX);
  j_iter["bsh_residualsY"] = tensor_to_json(bsh_residualsY);
  j_iter["density_residuals"] = tensor_to_json(density_residuals);
  j_iter["polarizability"] = tensor_to_json(polar);
  return j_iter;
}

// Iterate Frequency Response
void TDDFT::iterate_freq2(World& world) {
  size_t iter;
  // Variables needed to iterate
  QProjector<double, 3> projector(world, ground_orbitals);
  size_t n = r_params.num_orbitals();  // Number of ground state orbitals
  size_t m = r_params.n_states();      // Number of excited states

  real_function_3d v_xc;   // For TDDFT
  bool converged = false;  // Converged flag
  const double dconv =
      std::max(FunctionDefaults<3>::get_thresh(), r_params.dconv());
  // m residuals for x and y
  Tensor<double> bsh_residualsX(m);
  Tensor<double> bsh_residualsY(m);
  Tensor<double> density_residuals(m);

  vecfuncT rho_omega_old(m);

  // initialize DFT XC functional operator
  XCOperator<double, 3> xc =
      create_XCOperator(world, ground_orbitals, r_params.xc());

  // create X space residuals
  X_space residuals(world, m, n);
  X_space old_Chi(world, m, n);
  // Create the X space
  // vector of Xvectors
  std::vector<X_vector> Xvector;
  std::vector<X_vector> Xresidual;
  for (size_t b = 0; b < m; b++) {
    Xvector.push_back(X_vector(Chi, b));
    Xresidual.push_back(X_vector(residuals, b));
  }
  // If DFT, initialize the XCOperator<double,3>

  // create a std vector of XNONLinearsolvers
  NonLinearXsolver kain_x_space;
  for (size_t b = 0; b < m; b++) {
    kain_x_space.push_back(
        XNonlinearSolver<X_vector, double, X_space_allocator>(
            X_space_allocator(world, n), true));
  }
  for (size_t b = 0; b < m; b++) {
    if (r_params.kain()) kain_x_space[b].set_maxsub(r_params.maxsub());
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
  std::vector<poperatorT> bsh_x_ops =
      make_bsh_operators_response(world, x_shifts, omega_n);
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
  json j_frequency = {};
  j_frequency["num_states"] = m;
  j_frequency["num_orbitals"] = n;
  j_frequency["omega"] = omega_n;
  j_frequency["iter_data"] = json{};

  for (iter = 0; iter <= r_params.maxiter(); ++iter) {
    // Basic output
    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
      if (world.rank() == 0)
        printf("\n   Iteration %d at time %.1fs\n",
               static_cast<int>(iter),
               wall_time());
      if (world.rank() == 0)
        print("-------------------------------------------");
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
      load_balance(world, rho_omega, Chi, old_Chi);
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
          ((std::max(bsh_residualsX.absmax(), bsh_residualsY.absmax()) <
            d_conv * 5.0) or
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
          if (r_params.print_level() >= 1)
            molresponse::end_timer(world, "Save:");
        }
        // Basic output
        if (r_params.print_level() >= 1)
          molresponse::end_timer(world, " This iteration:");
        // plot orbitals
        if (r_params.plot_all_orbitals()) {
          PlotGroundandResponseOrbitals(
              world, iter, Chi.X, Chi.Y, r_params, g_params);
        }
        rho0 = make_ground_density(world, ground_orbitals);
        if (r_params.plot()) {
          do_vtk_plots(world,
                       200,
                       r_params.L(),
                       molecule,
                       rho0,
                       rho_omega,
                       ground_orbitals,
                       Chi);
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

  j_frequency["iter_data"].push_back(freq_iteration_to_json(
      iter, bsh_residualsX, bsh_residualsY, density_residuals, polar));

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
    compute_and_print_polarizability(world, Chi, PQ, "Converged");
  }
  std::ofstream ofs;  // open json file in append mode
  ofs.open("j_frequency.json");
  ofs << j_frequency;
}
