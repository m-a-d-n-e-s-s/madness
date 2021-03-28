// Copyright 2021 Adrian Hurtado
#include <math.h>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "../chem/SCFOperators.h"
#include "../chem/molecule.h"
#include "NWChem.h"  // For nwchem interface
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

void TDDFT::Iterate(World& world, X_space& Chi) {
  // Variables needed to iterate
  size_t iteration = 0;  // Iteration counter
  QProjector<double, 3> projector(world, Gparams.orbitals);
  size_t m = Rparams.states;        // Number of excited states
  size_t n = Gparams.num_orbitals;  // Number of ground state orbitals

  bool all_converged = false;                 // For convergence
  bool relax = false;                         // For convergence
  size_t relax_start = Rparams.max_iter + 1;  // For convergence
  size_t num_conv = 0;                        // For convergence
  std::vector<bool> converged(m, false);      // For convergence

  response_space bsh_x_resp(world, m, n);  // Holds wave function corrections
  response_space bsh_y_resp(world, m, n);  // Holds wave function corrections

  Tensor<double> old_energy(m);     // Holds previous iteration's energy
  Tensor<double> energy_residuals;  // Holds energy residuals
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
  Tensor<double> old_S;
  // initialize DFT XC functional operator
  XCOperator xc = create_xcoperator(world, Gparams.orbitals, Rparams.xc);

  /*
   * X space refers to X and Y vector spaces |X,Y>
   * X vector is a single |X_b,Y_b> b is a single response state
   * For kain we need a vector of X_vectors
   */

  /**
   * IF TDA then y_response =0
   */

  X_space residuals(world, m, n);
  X_space old_Chi(world, m, n);
  X_space old_Lambda_X(world, m, n);
  // create X space residuals
  // vector of Xvectors
  std::vector<X_vector> Xvector;
  std::vector<X_vector> Xresidual;

  for (size_t b = 0; b < m; b++) {
    Xvector.push_back(X_vector(Chi, b));
    Xresidual.push_back(X_vector(residuals, b));
  }
  // If DFT, initialize the XCOperator
  std::vector<XNonlinearSolver<X_vector, double, X_space_allocator>>
      kain_x_space;
  size_t nkain = m;  // (Rparams.omega != 0.0) ? 2 * m : m;
  for (size_t b = 0; b < nkain; b++) {
    kain_x_space.push_back(
        XNonlinearSolver<X_vector, double, X_space_allocator>(
            X_space_allocator(world, n), false));
    if (Rparams.kain) kain_x_space[b].set_maxsub(Rparams.maxsub);
  }

  // Here we create the kain solver for response
  // functions...TDHF_allocator returns a Response function of correct
  // size With correct number of occ and virtual orbitals all set to zero.

  // Set y things if not doing TDA
  if (not Rparams.tda) old_Chi.Y = response_space(world, m, n);

  // Now to iterate
  while (iteration < Rparams.max_iter && !all_converged) {
    // Start a timer for this iteration
    molresponse::start_timer(world);
    // Basic output
    if (Rparams.print_level >= 1) {
      if (world.rank() == 0)
        printf("\n   Iteration %d at time %.1fs\n",
               static_cast<int>(iteration),
               wall_time());
      if (world.rank() == 0) print(" -------------------------------");
    }
    Chi.X.truncate_rf();
    if (not Rparams.tda) Chi.Y.truncate_rf();

    // Normalize after projection
    if (Rparams.tda) {
      normalize(world, Chi.X);
    } else {
      normalize(world, Chi);
    }

    X_space Lambda_X = Compute_Lambda_X(world, Chi, xc, not Rparams.tda);
    // Load balance
    // Only balancing on x-components. Smart?
    if ((world.size() > 1 && (iteration < 2 or iteration % 5 == 0))) {
      // Start a timer
      if (Rparams.print_level >= 1) molresponse::start_timer(world);
      if (world.rank() == 0) print("");  // Makes it more legible

      // TODO Ask about load balancing
      LoadBalanceDeux<3> lb(world);
      for (size_t j = 0; j < n; j++) {
        for (size_t k = 0; k < Rparams.states; k++) {
          lb.add_tree(Chi.X[k][j], lbcost<double, 3>(1.0, 8.0), true);
          lb.add_tree(Lambda_X.X[k][j], lbcost<double, 3>(1.0, 8.0), true);
        }
      }
      FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

      if (Rparams.print_level >= 1)
        molresponse::end_timer(world, "Load balancing:");
    }

    if (Rparams.print_level >= 1 and world.rank() == 0) {
      print("Before Deflate");
      print("\n   Excitation Energies:");
      print("i=", iteration, " roots: ", iteration, omega);
    }
    // TDA approximation
    if (Rparams.tda) {
      deflateTDA(world,
                 Chi,
                 old_Chi,
                 Lambda_X,
                 old_Lambda_X,
                 S,
                 old_S,
                 old_A,
                 omega,
                 iteration,
                 m);
      // Constructing S
      // Full TDHF
    } else {
      // Constructing S
      deflateFull(world,
                  Chi,
                  old_Chi,
                  Lambda_X,
                  old_Lambda_X,
                  S,
                  old_S,
                  old_A,
                  omega,
                  iteration,
                  m);
    }

    // Basic output
    if (Rparams.print_level >= 1 and world.rank() == 0) {
      print("After Deflate");
      print("\n   Excitation Energies:");
      print("i=", iteration, " roots: ", iteration, omega);
      print("\n Hamiltonian: ");
      print(hamiltonian);
      print("\n Hamiltonian NoDiag: ");
      print(ham_no_diag);
      print("norm");
      print(hamiltonian.normf());
    }

    // Calculate energy residual and update old_energy
    energy_residuals = abs(omega - old_energy);

    old_energy = copy(omega);

    // Basic output
    if (Rparams.print_level >= 1) {
      if (world.rank() == 0) print("   Energy residuals:");
      if (world.rank() == 0) print("er: ", iteration, " ", energy_residuals);
    }

    // Analysis gets messed up if BSH is last thing applied
    // so exit early if last iteration
    if (iteration == Rparams.max_iter - 1) {
      molresponse::end_timer(world, " This iteration:");
      break;
    }

    //  Calculates shifts needed for potential / energies
    //  If none needed, the zero tensor is returned
    x_shifts =
        create_shift(world, Gparams.energies, omega, Rparams.print_level, "x");
    if (not Rparams.tda) {
      omega = -omega;  // Negative here is so that these Greens functions are
      // (eps - omega)
      y_shifts = create_shift_target(world,
                                     Gparams.energies,
                                     omega,
                                     Gparams.energies[n - 1],
                                     Rparams.print_level,
                                     "y");
      omega = -omega;
    }
    print(
        "----------------Before Compute_Theta_X After Deflate "
        "-----------------");
    if (Rparams.print_level == 3) {
      print("x norms in iteration after copy  : ", iteration);
      print(Chi.X.norm2());
      print("y norms in iteration after copy: ", iteration);
      print(Chi.Y.norm2());
    }

    // Compute Theta X
    X_space theta_X = Compute_Theta_X(world, Chi, xc, not Rparams.tda);
    // Apply the shifts
    theta_X.X = apply_shift(world, x_shifts, theta_X.X, Chi.X);
    theta_X.X = theta_X.X * -2;
    theta_X.X.truncate_rf();
    if (not Rparams.tda) {
      y_shifts = -y_shifts;
      theta_X.Y = apply_shift(world, y_shifts, theta_X.Y, Chi.Y);
      theta_X.Y = theta_X.Y * -2;
      theta_X.Y.truncate_rf();
    }
    if (not Rparams.tda) {
      // Debugging output
      if (Rparams.print_level >= 2) {
        if (world.rank() == 0) print("   Norms of RHS of main equation:");
        if (world.rank() == 0) print("   x components:");
        print_norms(world, theta_X.X);

        if (not Rparams.tda) {
          if (world.rank() == 0) print("   y components:");
          print_norms(world, theta_X.Y);
        }
      }

      // Construct BSH operators
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
          bsh_x_operators =
              create_bsh_operators(world,
                                   x_shifts,
                                   Gparams.energies,
                                   omega,
                                   Rparams.small,
                                   FunctionDefaults<3>::get_thresh());
      std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
          bsh_y_operators;
      if (not Rparams.tda) {
        omega = -omega;
        bsh_y_operators =
            create_bsh_operators(world,
                                 y_shifts,
                                 Gparams.energies,
                                 omega,
                                 Rparams.small,
                                 FunctionDefaults<3>::get_thresh());
        omega = -omega;
      }

      // Save current into old
      old_Chi = Chi.copy();
      X_space temp(world, m, n);

      // Apply BSH and get updated response components
      if (Rparams.print_level >= 1) molresponse::start_timer(world);
      bsh_x_resp = apply(world, bsh_x_operators, Chi.X);
      if (not Rparams.tda) bsh_y_resp = apply(world, bsh_y_operators, Chi.Y);
      if (Rparams.print_level >= 1) molresponse::end_timer(world, "Apply BSH:");

      // Debugging output
      if (Rparams.print_level >= 2) {
        if (world.rank() == 0) print("   Norms after application of BSH");
        if (world.rank() == 0) print("   x-components:");
        print_norms(world, bsh_x_resp);

        if (not Rparams.tda) {
          if (world.rank() == 0) print("   y-components:");
          print_norms(world, bsh_y_resp);
        }
      }

      // Project out ground state
      for (size_t i = 0; i < m; i++) bsh_x_resp[i] = projector(bsh_x_resp[i]);
      if (not Rparams.tda) {
        for (size_t i = 0; i < m; i++) bsh_y_resp[i] = projector(bsh_y_resp[i]);
      }

      // Only update non-converged components
      for (size_t i = 0; i < m; i++) {
        if (not converged[i]) {
          x_response[i] = bsh_x_resp[i];
          if (not Rparams.tda) y_response[i] = bsh_y_resp[i];
        }
      }
      // Scale by -2.0 (coefficient in eq. 37 of reference paper)

      temp.X = bsh_x_resp.copy();
      if (Rparams.omega != 0.0) {
        temp.Y = bsh_y_resp.copy();
      } else {
        temp.Y = temp.X.copy();
      }
      temp.X.truncate_rf();
      temp.Y.truncate_rf();
      // Get the difference between old and new
      x_differences = old_Chi.X - temp.X;
      if (not Rparams.tda) y_differences = old_Chi.Y - temp.Y;

      // Next calculate 2-norm of these vectors of differences
      // Remember: the entire vector is one state
      for (size_t i = 0; i < m; i++)
        x_norms(i) = norm2(world, x_differences[i]);
      if (not Rparams.tda) {
        for (size_t i = 0; i < m; i++)
          y_norms(i) = norm2(world, y_differences[i]);
      }

      // Basic output
      if (Rparams.print_level >= 1) {
        if (world.rank() == 0)
          print("\n   2-norm of response function residuals:");
        if (world.rank() == 0) print("   x components:");
        if (world.rank() == 0) print(x_norms);

        if (not Rparams.tda) {
          if (world.rank() == 0) print("   y components:");
          if (world.rank() == 0) print(y_norms);
        }
      }
      if (Rparams.kain) {
        residuals = X_space(x_differences, y_differences);

        // seperate X_space vectors into individual vectors
        for (size_t b = 0; b < m; b++) {
          Xvector[b] = (X_vector(temp, b));
          Xresidual[b] = (X_vector(residuals, b));
        }

        molresponse::start_timer(world);
        for (size_t b = 0; b < nkain; b++) {
          X_vector kain_X = kain_x_space[b].update(
              Xvector[b], Xresidual[b], FunctionDefaults<3>::get_thresh(), 3.0);
          x_response[b].assign(kain_X.X[0].begin(), kain_X.X[0].end());
          y_response[b].assign(kain_X.Y[0].begin(), kain_X.Y[0].end());
        }
        molresponse::end_timer(world, " KAIN update:");
      }  // end kain
      // KAIN solver update
      // Returns next set of components
      // If not kain, save the new components
      // do step restriction
      if (iteration > 0) {
        for (size_t b = 0; b < m; b++) {
          do_step_restriction(world, old_Chi.X[b], temp.X[b], "x_response");
          if (not Rparams.tda) {
            do_step_restriction(world, old_Chi.Y[b], temp.X[b], "y_response");
          }
        }
      }
      // Apply mask
      for (size_t i = 0; i < m; i++) temp.X[i] = mask * temp.X[i];
      if (not Rparams.tda) {
        for (size_t i = 0; i < m; i++) temp.Y[i] = mask * temp.Y[i];
      }
      temp.X.truncate_rf();
      if (not Rparams.tda) temp.Y.truncate_rf();
      // temp-> Chi
      Chi = temp.copy();
      if (Rparams.print_level >= 1) {
        print("Chi.x norms in iteration after truncate: ", iteration);
        print(Chi.X.norm2());

        print("Chi.y norms in iteration after truncate: ", iteration);
        print(Chi.Y.norm2());
      }

      // Only checking on X components even for full as Y are so small
      if (not relax) {
        for (size_t i = 0; i < m; i++) {
          if (iteration >= 1 && not converged[i] &&
              fabs(x_norms[i]) < Rparams.dconv) {
            converged[i] = true;
            num_conv++;
            if (world.rank() == 0)
              print("   Response function", i, " has converged. Freezing it.");
          }
        }
        // Check if relaxing needs to start
        if (num_conv == m) {
          relax_start = iteration;
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
        if (iteration >= relax_start + 2) {
          // Check each state again
          for (size_t i = 0; i < m; i++) {
            if (not converged[i] && fabs(x_norms[i]) < Rparams.dconv) {
              converged[i] = true;
              num_conv++;
            }
          }
          if (num_conv == m) all_converged = true;
        }
      }

      // Update counter
      iteration += 1;

      // Save
      if (Rparams.save) {
        molresponse::start_timer(world);
        save(world, Rparams.save_file);
        molresponse::end_timer(world, "Saving:");
      }

      // Basic output
      if (Rparams.print_level >= 1)
        molresponse::end_timer(world, " This iteration:");

      if (Rparams.plot_all_orbitals) {
        PlotGroundandResponseOrbitals(
            world, iteration, Chi.X, Chi.Y, Rparams, Gparams);
      }
    }

    if (world.rank() == 0) print("\n");
    if (world.rank() == 0) print("   Finished TDHF Calculation ");
    if (world.rank() == 0) print("   ------------------------");
    if (world.rank() == 0) print("\n");

    // Did we converge?
    if (iteration == Rparams.max_iter && not all_converged) {
      if (world.rank() == 0) print("   Failed to converge. Reason:");
      if (world.rank() == 0) print("\n  ***  Ran out of iterations  ***\n");
      if (world.rank() == 0) print("    Running analysis on current values.\n");
    }

    // Sort
    sort(world, omega, x_response);

    // Print final things
    if (world.rank() == 0) {
      print(" Final excitation energies:");
      print(omega);
      print(" Final energy residuals:");
      print(energy_residuals);
      print(" Final x-state response function residuals:");
      print(x_norms);

      if (not Rparams.tda) {
        if (world.rank() == 0)
          print(" Final y-state response function residuals:");
        if (world.rank() == 0) print(y_norms);
      }
    }

    // A little more detailed analysis
    analysis(world);

    // TEST
    // Doing line plots along z axis
    // if(world.rank() == 0) print("\n\nStarting plots");
    // coord_3d lo,hi;
    // char plotname[500];
    //// z axis
    // lo[0] = 0.0; lo[1] = 0.0; lo[2] = -Gparams.L;
    // hi[0] = 0.0; hi[1] = 0.0; hi[2] =  Gparams.L;
    // for(size_t    i = 0; i < Rparams.states; i++) {
    //  for(unsigned int j = 0; j < Gparams.num_orbitals; j++) {
    //    sprintf(plotname, "plot_exX_%d_%d.plt", i, j);
    //    plot_line(plotname, 500001, lo, hi, x_response[i][j]);
    //    sprintf(plotname, "plot_exY_%d_%d.plt", i, j);
    //    plot_line(plotname, 500001, lo, hi, y_response[i][j]);
    //  }
    //}
    // END TEST
  }  // Done with iterate.

  // More detailed analysis of the response functions
  // Uses member variables
}
