// Copyright 2021 Adrian Hurtado
#include <madness/world/worldmem.h>
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
#include "chem/potentialmanager.h"
#include "chem/projector.h"  // For easy calculation of (1 - \hat{\rho}^0)
#include "madness/mra/funcdefaults.h"
#include "molresponse/basic_operators.h"
#include "molresponse/density.h"
#include "molresponse/global_functions.h"
#include "molresponse/load_balance.h"
#include "molresponse/property.h"
#include "molresponse/response_functions.h"
#include "molresponse/timer.h"
#include "molresponse/x_space.h"

void TDDFT::iterate_excited(World& world, X_space& Chi) {
  size_t iter;
  QProjector<double, 3> projector(world, ground_orbitals);
  size_t m = r_params.n_states();      // Number of excited states
  size_t n = r_params.num_orbitals();  // Number of ground state orbitals

  bool all_converged = false;                   // For convergence
  bool relax = false;                           // For convergence
  size_t relax_start = r_params.maxiter() + 1;  // For convergence
  size_t num_conv = 0;                          // For convergence
  std::vector<bool> converged(m, false);        // For convergence
  const double dconv =
      std::max(FunctionDefaults<3>::get_thresh(), r_params.dconv());

  // m residuals for x and y
  Tensor<double> bsh_residualsX(m);
  Tensor<double> bsh_residualsY(m);
  // saved response densities
  vecfuncT rho_omega_old(m);
  // initialize DFT XC functional operator
  XCOperator<double, 3> xc =
      create_XCOperator(world, ground_orbitals, r_params.xc());

  // create X space residuals
  X_space residuals(world, m, n);
  X_space old_Chi(world, m, n);
  X_space old_Lambda_X(world, m, n);

  // Create the X space
  // vector of Xvectors
  std::vector<X_vector> Xvector;
  std::vector<X_vector> Xresidual;
  for (size_t b = 0; b < m; b++) {
    Xvector.push_back(X_vector(Chi, b));
    Xresidual.push_back(X_vector(residuals, b));
  }
  // If DFT, initialize the XCOperator<double,3>

  NonLinearXsolver kain_x_space;
  size_t nkain = m;  // (r_params.omega() != 0.0) ? 2 * m : m;
  for (size_t b = 0; b < nkain; b++) {
    kain_x_space.push_back(
        XNonlinearSolver<X_vector, double, X_space_allocator>(
            X_space_allocator(world, n), false));
    if (r_params.kain()) kain_x_space[b].set_maxsub(r_params.maxsub());
  }

  response_space bsh_x_resp(world, m, n);  // Holds wave function corrections
  response_space bsh_y_resp(world, m, n);  // Holds wave function corrections

  Tensor<double> old_energy(m);        // Holds previous iteration's energy
  Tensor<double> energy_residuals(m);  // Holds energy residuals
  // Holds the norms of y function residuals (for convergence)
  Tensor<double> x_norms(m);
  Tensor<double> y_norms(m);

  response_space x_differences(world, m, n);
  response_space y_differences(world, m, n);

  response_space x_residuals(world, m, n);
  response_space y_residuals(world, m, n);

  Tensor<double> x_shifts;  // Holds the shifted energy values
  Tensor<double> y_shifts;  // Holds the shifted energy values

  Tensor<double> S;       // Overlap matrix of response components for x states
  real_function_3d v_xc;  // For TDDFT
  Tensor<double> old_A;
  Tensor<double> A;
  Tensor<double> old_S;
  bool compute_y = !r_params.tda();

  // Now to iterate
  for (iter = 0; iter < r_params.maxiter(); ++iter) {
    // Start a timer for this iteration
    // Basic output
    if (r_params.print_level() >= 1) {
      molresponse::start_timer(world);
      if (world.rank() == 0)
        printf("\n   Iteration %d at time %.1fs\n",
               static_cast<int>(iter),
               wall_time());
      if (world.rank() == 0) print(" -------------------------------");
    }

    // compute rho_omega
    rho_omega = make_density(world, Chi, compute_y);

    // Normalize after projection
    if (r_params.tda()) {
      normalize(world, Chi.X);
    } else {
      normalize(world, Chi);
    }

    if (iter < 2 || (iter % 10) == 0) {
      loadbal(world, rho_omega, Chi, old_Chi);
    }

    // compute density residuals
    vector<double> density_residuals;
    if (iter > 0) {
      density_residuals = norm2s(world, (rho_omega - rho_omega_old));
      if (world.rank() == 0 and (r_params.print_level() > 2))
        print("delta rho",
              density_residuals,
              "residuals",
              bsh_residualsX,
              bsh_residualsY);
    }
    //
    old_Chi = Chi.copy();
    rho_omega_old = rho_omega;

    if (iter > 0) {
      // Only checking on X components even for full as Y are so small
      if (not relax) {
        for (size_t i = 0; i < m; i++) {
          // bsh_residual max orbital change after bsh apply
          if (not converged[i] && fabs(bsh_residualsX[i]) < dconv) {
            converged[i] = true;
            num_conv++;
            if (world.rank() == 0)
              print("   Response function", i, " has converged. Freezing it.");
          }
        }
        // Check if relaxing needs to start
        if (num_conv == m) {
          relax_start = iter;
          relax = true;
          if (world.rank() == 0)
            print(
                "   All components converged. Unfreezing all states for "
                "final "
                "relaxation.");

          num_conv = 0;
          for (size_t i = 0; i < m; i++) {
            converged[i] = false;
          }
        }
      } else {
        // Relaxing
        // Run at least 2 iterations
        if (iter >= relax_start + 2) {
          // Check each state again
          for (size_t i = 0; i < m; i++) {
            if (not converged[i] && fabs(bsh_residualsX[i]) < dconv) {
              converged[i] = true;
              num_conv++;
            }
          }
          if (num_conv == m || iter == r_params.maxiter() - 1) {
            all_converged = true;
            // Save
            if (r_params.save()) {
              molresponse::start_timer(world);
              save(world, r_params.save_file());
              molresponse::end_timer(world, "Saving:");
            }

            if (r_params.plot_all_orbitals()) {
              plot_excited_states(
                  world, iter, Chi.X, Chi.Y, r_params, g_params);
            }
            break;
          }
        }
      }
    }
    update_x_space_excited(world,
                           old_Chi,
                           Chi,
                           old_Lambda_X,
                           residuals,
                           xc,
                           projector,
                           omega,
                           kain_x_space,
                           Xvector,
                           Xresidual,
                           energy_residuals,
                           old_energy,
                           bsh_residualsX,
                           bsh_residualsY,
                           S,
                           old_S,
                           A,
                           old_A,
                           converged,
                           iter);

    // Basic output
    if (r_params.print_level() >= 1)
      molresponse::end_timer(world, " This iteration:");
  }

  if (world.rank() == 0) print("\n");
  if (world.rank() == 0) print("   Finished Excited State Calculation ");
  if (world.rank() == 0) print("   ------------------------");
  if (world.rank() == 0) print("\n");

  // Did we converge?
  if (iter == r_params.maxiter() && not all_converged) {
    if (world.rank() == 0) print("   Failed to converge. Reason:");
    if (world.rank() == 0) print("\n  ***  Ran out of iterations  ***\n");
    if (world.rank() == 0) print("    Running analysis on current values.\n");
  }

  // Sorstatict
  sort(world, omega, Chi.X);

  // Print final things
  if (world.rank() == 0) {
    print(" Final excitation energies:");
    print(omega);
    print(" Final energy residuals:");
    print(energy_residuals);
    print(" Final x-state response function residuals:");
    print(x_norms);

    if (not r_params.tda()) {
      if (world.rank() == 0)
        print(" Final y-state response function residuals:");
      if (world.rank() == 0) print(y_norms);
    }
  }

  analysis(world, Chi);
}

void TDDFT::analysis(World& world, X_space& Chi) {
  // Sizes get used a lot here, so lets get a local copy
  size_t n = Chi.X[0].size();
  size_t m = Chi.X.size();

  // Per response function, want to print the contributions from each
  // ground state So print the norm of each function?
  Tensor<double> x_norms(m, n);
  Tensor<double> y_norms(m, n);

  // Calculate the inner products
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      x_norms(i, j) = Chi.X[i][j].norm2();

      if (not r_params.tda()) y_norms(i, j) = Chi.Y[i][j].norm2();
    }
  }

  // 'sort' these inner products within in each row
  Tensor<double> cpy = copy(x_norms);
  Tensor<int> x_order(m, n);
  Tensor<int> y_order(m, n);
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      double x = cpy(i, _).max();
      size_t z = 0;
      while (x != cpy(i, z)) z++;
      cpy(i, z) = -100.0;
      x_order(i, j) = z;

      // Also sort y if full response
      if (not r_params.tda()) {
        y_order(i, j) = z;
      }
    }
  }

  // Need these to calculate dipole/quadrapole
  real_function_3d x = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{1, 0, 0})));
  real_function_3d y = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 1, 0})));
  real_function_3d z = real_factory_3d(world).functor(
      real_functor_3d(new BS_MomentFunctor(std::vector<int>{0, 0, 1})));

  // Calculate transition dipole moments for each response function
  Tensor<double> dipoles(m, 3);

  // Run over each excited state
  for (size_t i = 0; i < m; i++) {
    // Add in contribution from each ground state
    for (size_t j = 0; j < n; j++) {
      dipoles(i, 0) += inner(ground_orbitals[j], x * Chi.X[i][j]);
      dipoles(i, 1) += inner(ground_orbitals[j], y * Chi.X[i][j]);
      dipoles(i, 2) += inner(ground_orbitals[j], z * Chi.X[i][j]);

      if (not r_params.tda()) {
        dipoles(i, 0) += inner(ground_orbitals[j], x * Chi.Y[i][j]);
        dipoles(i, 1) += inner(ground_orbitals[j], y * Chi.Y[i][j]);
        dipoles(i, 2) += inner(ground_orbitals[j], z * Chi.Y[i][j]);
      }
    }

    // Normalization (negative?)
    dipoles(i, 0) *= -sqrt(2.0);
    dipoles(i, 1) *= -sqrt(2.0);
    dipoles(i, 2) *= -sqrt(2.0);
  }

  // Calculate oscillator strength
  Tensor<double> oscillator(m);
  for (size_t i = 0; i < m; i++) {
    oscillator(i) =
        2.0 / 3.0 *
        (dipoles(i, 0) * dipoles(i, 0) + dipoles(i, 1) * dipoles(i, 1) +
         dipoles(i, 2) * dipoles(i, 2)) *
        omega(i);
  }

  // Calculate transition quadrapole moments
  Tensor<double> quadrupoles(m, 3, 3);

  // Run over each excited state
  for (size_t i = 0; i < m; i++) {
    // Add in contribution from each ground state
    for (size_t j = 0; j < n; j++) {
      quadrupoles(i, 0, 0) += inner(ground_orbitals[j], x * x * Chi.X[i][j]);
      quadrupoles(i, 0, 1) += inner(ground_orbitals[j], x * y * Chi.X[i][j]);
      quadrupoles(i, 0, 2) += inner(ground_orbitals[j], x * z * Chi.X[i][j]);
      quadrupoles(i, 1, 0) += inner(ground_orbitals[j], y * x * Chi.X[i][j]);
      quadrupoles(i, 1, 1) += inner(ground_orbitals[j], y * y * Chi.X[i][j]);
      quadrupoles(i, 1, 2) += inner(ground_orbitals[j], y * z * Chi.X[i][j]);
      quadrupoles(i, 2, 0) += inner(ground_orbitals[j], z * x * Chi.X[i][j]);
      quadrupoles(i, 2, 1) += inner(ground_orbitals[j], z * y * Chi.X[i][j]);
      quadrupoles(i, 2, 2) += inner(ground_orbitals[j], z * z * Chi.X[i][j]);

      if (not r_params.tda()) {
        quadrupoles(i, 0, 0) += inner(ground_orbitals[j], x * x * Chi.Y[i][j]);
        quadrupoles(i, 0, 1) += inner(ground_orbitals[j], x * y * Chi.Y[i][j]);
        quadrupoles(i, 0, 2) += inner(ground_orbitals[j], x * z * Chi.Y[i][j]);
        quadrupoles(i, 1, 0) += inner(ground_orbitals[j], y * x * Chi.Y[i][j]);
        quadrupoles(i, 1, 1) += inner(ground_orbitals[j], y * y * Chi.Y[i][j]);
        quadrupoles(i, 1, 2) += inner(ground_orbitals[j], y * z * Chi.Y[i][j]);
        quadrupoles(i, 2, 0) += inner(ground_orbitals[j], z * x * Chi.Y[i][j]);
        quadrupoles(i, 2, 1) += inner(ground_orbitals[j], z * y * Chi.Y[i][j]);
        quadrupoles(i, 2, 2) += inner(ground_orbitals[j], z * z * Chi.Y[i][j]);
      }
    }
    // Normalization
    quadrupoles(i, 0, 0) *= sqrt(2.0);
    quadrupoles(i, 0, 1) *= sqrt(2.0);
    quadrupoles(i, 0, 2) *= sqrt(2.0);
    quadrupoles(i, 1, 0) *= sqrt(2.0);
    quadrupoles(i, 1, 1) *= sqrt(2.0);
    quadrupoles(i, 1, 2) *= sqrt(2.0);
    quadrupoles(i, 2, 0) *= sqrt(2.0);
    quadrupoles(i, 2, 1) *= sqrt(2.0);
    quadrupoles(i, 2, 2) *= sqrt(2.0);
  }

  // Now print?
  if (world.rank() == 0) {
    for (size_t i = 0; i < m; i++) {
      printf("   Response Function %d\t\t%7.8f a.u.",
             static_cast<int>(i),
             omega(i));
      print("\n   --------------------------------------------");

      print("\n   Transition Dipole Moments");
      printf("   X: %7.8f   Y: %7.8f   Z: %7.8f\n",
             dipoles(i, 0),
             dipoles(i, 1),
             dipoles(i, 2));

      printf("\n   Dipole Oscillator Strength: %7.8f\n", oscillator(i));

      print("\n   Transition Quadrupole Moments");
      printf("   %16s %16s %16s\n", "X", "Y", "Z");
      printf("   X %16.8f %16.8f %16.8f\n",
             quadrupoles(i, 0, 0),
             quadrupoles(i, 0, 1),
             quadrupoles(i, 0, 2));
      printf("   Y %16.8f %16.8f %16.8f\n",
             quadrupoles(i, 1, 0),
             quadrupoles(i, 1, 1),
             quadrupoles(i, 1, 2));
      printf("   Z %16.8f %16.8f %16.8f\n",
             quadrupoles(i, 2, 0),
             quadrupoles(i, 2, 1),
             quadrupoles(i, 2, 2));

      // Print contributions
      // Only print the top 5?
      if (r_params.tda()) {
        print("\n   Dominant Contributions:");
        for (size_t j = 0; j < std::min(size_t(5), n); j++) {
          printf("   Occupied %d   %7.8f\n",
                 x_order(i, j),
                 x_norms(i, x_order(i, j)));
        }

        print("\n");
      } else {
        print("\n   Dominant Contributions:");
        print("                  x          y");
        for (size_t j = 0; j < std::min(size_t(5), n); j++) {
          printf("   Occupied %d   %7.8f %7.8f\n",
                 x_order(i, j),
                 x_norms(i, x_order(i, j)),
                 y_norms(i, y_order(i, j)));
        }

        print("\n");
      }
    }
  }
}
