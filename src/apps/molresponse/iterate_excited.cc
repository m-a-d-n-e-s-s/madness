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

void TDDFT::Iterate(World& world) {
  // Variables needed to iterate
  size_t iteration = 0;  // Iteration counter
                         // Projector to project out ground state
  QProjector<double, 3> projector(world, Gparams.orbitals);
  size_t m = Rparams.states;                  // Number of excited states
  size_t n = Gparams.num_orbitals;            // Number of ground state orbitals
  bool all_converged = false;                 // For convergence
  bool relax = false;                         // For convergence
  size_t relax_start = Rparams.max_iter + 1;  // For convergence
  size_t num_conv = 0;                        // For convergence
  std::vector<bool> converged(m, false);      // For convergence
  Tensor<double> old_energy(m);     // Holds previous iteration's energy
  Tensor<double> energy_residuals;  // Holds energy residuals
  // Holds the norms of x function residuals (for convergence)
  Tensor<double> x_norms(m);
  // Holds the norms of y function residuals (for convergence)
  Tensor<double> y_norms(m);
  Tensor<double> x_shifts;       // Holds the shifted energy values
  Tensor<double> y_shifts;       // Holds the shifted energy values
  response_space rhs_x;          // Holds wave function corrections
  response_space rhs_y;          // Holds wave function corrections
  response_space bsh_x_resp;     // Holds wave function corrections
  response_space bsh_y_resp;     // Holds wave function corrections
  response_space x_differences;  // Holds wave function corrections
  response_space y_differences;  // Holds wave function corrections
  // Holds the shifted V^0 applied to response functions
  response_space shifted_V_x_response(world, m, n);
  response_space shifted_V_y_response(world, m, n);
  // Holds the old x_response vector of vectors
  response_space old_x_response(world, m, n);
  response_space old_y_response(world, m, n);
  Tensor<double> S;       // Overlap matrix of response components for x states
  real_function_3d v_xc;  // For TDDFT

  ElectronResponseFunctions ElectronResponses;
  ElectronResponseFunctions OldElectronResponses;
  // Versions from previous iteration that need to be stored
  // in order to diagonalize in a larger subspace
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

  // create X space residuals
  X_space residuals(world, m, n);
  // Create the X space
  X_space X(x_response, y_response);
  // vector of Xvectors
  std::vector<X_vector> Xvector;
  std::vector<X_vector> Xresidual;

  for (size_t b = 0; b < m; b++) {
    Xvector.push_back(X_vector(X, b));
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
  if (not Rparams.tda) old_y_response = response_space(world, m, n);

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
    // Truncate before doing expensive things
    // Truncate x response and y response meaning??? remove coefficients
    // based(TODO) on a global paramater??
    x_response.truncate_rf();
    if (not Rparams.tda) y_response.truncate_rf();

    // Normalize after projection
    if (Rparams.tda) {
      normalize(world, x_response);
    } else {
      normalize(world, x_response, y_response);
    }

    computeElectronResponse(world,
                            ElectronResponses,
                            x_response,
                            y_response,
                            Gparams.orbitals,
                            xc,
                            hamiltonian,
                            ham_no_diag,
                            Rparams.small,
                            FunctionDefaults<3>::get_thresh(),
                            Rparams.print_level,
                            "x");
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
          lb.add_tree(x_response[k][j], lbcost<double, 3>(1.0, 8.0), true);
          lb.add_tree(
              ElectronResponses.Vx[k][j], lbcost<double, 3>(1.0, 8.0), true);
          lb.add_tree(
              ElectronResponses.Hx[k][j], lbcost<double, 3>(1.0, 8.0), true);
          lb.add_tree(
              ElectronResponses.Hx[k][j], lbcost<double, 3>(1.0, 8.0), true);
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
                 S,
                 old_S,
                 old_A,
                 x_response,
                 old_x_response,
                 ElectronResponses,
                 OldElectronResponses,
                 omega,
                 iteration,
                 m);
      // Constructing S
      // Full TDHF
    } else {
      // Constructing S
      deflateFull(world,
                  S,
                  old_S,
                  old_A,
                  x_response,
                  y_response,
                  old_x_response,
                  old_y_response,
                  ElectronResponses,
                  OldElectronResponses,
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

    // Apply the shifts
    shifted_V_x_response =
        apply_shift(world, x_shifts, ElectronResponses.Vx, x_response);
    if (not Rparams.tda) {
      y_shifts = -y_shifts;
      shifted_V_y_response =
          apply_shift(world, y_shifts, ElectronResponses.Vy, y_response);
    }

    if (not Rparams.tda) {
      rhs_x = ElectronResponses.Hx + ElectronResponses.Gy;
      rhs_y = ElectronResponses.Hy + ElectronResponses.Gx;
      for (size_t i = 0; i < m; i++) {
        rhs_x[i] = projector(rhs_x[i]);
        rhs_y[i] = projector(rhs_y[i]);
      }
      rhs_x = shifted_V_x_response + rhs_x - ElectronResponses.EpsilonXNoDiag;
      rhs_y = shifted_V_y_response + rhs_y - ElectronResponses.EpsilonYNoDiag;
    } else {
      rhs_x = ElectronResponses.Hx;
      for (size_t i = 0; i < m; i++) {
        rhs_x[i] = projector(rhs_x[i]);
      }
      rhs_x = shifted_V_x_response + rhs_x - ElectronResponses.EpsilonXNoDiag;
    }

    // Debugging output
    if (Rparams.print_level >= 2) {
      if (world.rank() == 0) print("   Norms of RHS of main equation:");
      if (world.rank() == 0) print("   x components:");
      print_norms(world, rhs_x);

      if (not Rparams.tda) {
        if (world.rank() == 0) print("   y components:");
        print_norms(world, rhs_y);
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
      bsh_y_operators = create_bsh_operators(world,
                                             y_shifts,
                                             Gparams.energies,
                                             omega,
                                             Rparams.small,
                                             FunctionDefaults<3>::get_thresh());
      omega = -omega;
    }

    // Save current into old
    old_x_response = x_response.copy();
    if (not Rparams.tda) old_y_response = y_response.copy();

    // Apply BSH and get updated response components
    if (Rparams.print_level >= 1) molresponse::start_timer(world);
    bsh_x_resp = apply(world, bsh_x_operators, rhs_x);
    if (not Rparams.tda) bsh_y_resp = apply(world, bsh_y_operators, rhs_y);
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
    x_response = x_response * -2.0;  // scale(x_response, -2.0);
    if (not Rparams.tda)
      y_response = y_response * -2;  // scale(x_response, -2.0);

    // Get the difference between old and new
    x_differences = old_x_response - x_response;
    if (not Rparams.tda) y_differences = old_y_response - y_response;

    // Next calculate 2-norm of these vectors of differences
    // Remember: the entire vector is one state
    for (size_t i = 0; i < m; i++) x_norms(i) = norm2(world, x_differences[i]);
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
      X = X_space(x_response, y_response);

      residuals = X_space(x_differences, y_differences);

      // seperate X_space vectors into individual vectors
      for (size_t b = 0; b < m; b++) {
        Xvector[b] = (X_vector(X, b));
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
        do_step_restriction(
            world, old_x_response[b], x_response[b], "x_response");
        if (not Rparams.tda) {
          do_step_restriction(
              world, old_y_response[b], y_response[b], "y_response");
        }
      }
    }
    // Apply mask
    for (size_t i = 0; i < m; i++) x_response[i] = mask * x_response[i];
    if (not Rparams.tda) {
      for (size_t i = 0; i < m; i++) y_response[i] = mask * y_response[i];
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

    // Done with the iteration.. truncate
    x_response.truncate_rf();
    if (not Rparams.tda) y_response.truncate_rf();

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
          world, iteration, x_response, y_response, Rparams, Gparams);
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
void TDDFT::analysis(World& world) {
  // Sizes get used a lot here, so lets get a local copy
  size_t n = x_response[0].size();
  size_t m = x_response.size();

  // Per response function, want to print the contributions from each
  // ground state So print the norm of each function?
  Tensor<double> x_norms(m, n);
  Tensor<double> y_norms(m, n);

  // Calculate the inner products
  for (size_t i = 0; i < m; i++) {
    for (size_t j = 0; j < n; j++) {
      x_norms(i, j) = x_response[i][j].norm2();

      if (not Rparams.tda) y_norms(i, j) = y_response[i][j].norm2();
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
      if (not Rparams.tda) {
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
      dipoles(i, 0) += inner(Gparams.orbitals[j], x * x_response[i][j]);
      dipoles(i, 1) += inner(Gparams.orbitals[j], y * x_response[i][j]);
      dipoles(i, 2) += inner(Gparams.orbitals[j], z * x_response[i][j]);

      if (not Rparams.tda) {
        dipoles(i, 0) += inner(Gparams.orbitals[j], x * y_response[i][j]);
        dipoles(i, 1) += inner(Gparams.orbitals[j], y * y_response[i][j]);
        dipoles(i, 2) += inner(Gparams.orbitals[j], z * y_response[i][j]);
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
      quadrupoles(i, 0, 0) +=
          inner(Gparams.orbitals[j], x * x * x_response[i][j]);
      quadrupoles(i, 0, 1) +=
          inner(Gparams.orbitals[j], x * y * x_response[i][j]);
      quadrupoles(i, 0, 2) +=
          inner(Gparams.orbitals[j], x * z * x_response[i][j]);
      quadrupoles(i, 1, 0) +=
          inner(Gparams.orbitals[j], y * x * x_response[i][j]);
      quadrupoles(i, 1, 1) +=
          inner(Gparams.orbitals[j], y * y * x_response[i][j]);
      quadrupoles(i, 1, 2) +=
          inner(Gparams.orbitals[j], y * z * x_response[i][j]);
      quadrupoles(i, 2, 0) +=
          inner(Gparams.orbitals[j], z * x * x_response[i][j]);
      quadrupoles(i, 2, 1) +=
          inner(Gparams.orbitals[j], z * y * x_response[i][j]);
      quadrupoles(i, 2, 2) +=
          inner(Gparams.orbitals[j], z * z * x_response[i][j]);

      if (not Rparams.tda) {
        quadrupoles(i, 0, 0) +=
            inner(Gparams.orbitals[j], x * x * y_response[i][j]);
        quadrupoles(i, 0, 1) +=
            inner(Gparams.orbitals[j], x * y * y_response[i][j]);
        quadrupoles(i, 0, 2) +=
            inner(Gparams.orbitals[j], x * z * y_response[i][j]);
        quadrupoles(i, 1, 0) +=
            inner(Gparams.orbitals[j], y * x * y_response[i][j]);
        quadrupoles(i, 1, 1) +=
            inner(Gparams.orbitals[j], y * y * y_response[i][j]);
        quadrupoles(i, 1, 2) +=
            inner(Gparams.orbitals[j], y * z * y_response[i][j]);
        quadrupoles(i, 2, 0) +=
            inner(Gparams.orbitals[j], z * x * y_response[i][j]);
        quadrupoles(i, 2, 1) +=
            inner(Gparams.orbitals[j], z * y * y_response[i][j]);
        quadrupoles(i, 2, 2) +=
            inner(Gparams.orbitals[j], z * z * y_response[i][j]);
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
      if (Rparams.tda) {
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

// Main function, makes sure everything happens in correcct order
// Solves for response components
void TDDFT::solve_excited_states(World& world) {
  // Get start time
  molresponse::start_timer(world);

  // Plotting input orbitals
  if (Rparams.plot_initial) {
    if (world.rank() == 0) print("\n   Plotting ground state densities.\n");
    if (Rparams.plot_L > 0.0)
      do_vtk_plots(world,
                   Rparams.plot_pts,
                   Rparams.plot_L,
                   0,
                   Gparams.num_orbitals,
                   Gparams.molecule,
                   square(world, Gparams.orbitals),
                   "ground");
    else
      do_vtk_plots(world,
                   Rparams.plot_pts,
                   Gparams.L / 2.0,
                   0,
                   Gparams.num_orbitals,
                   Gparams.molecule,
                   square(world, Gparams.orbitals),
                   "ground");
  }

  // Warm and fuzzy
  if (world.rank() == 0) {
    print("\n\n     Response Calculation");
    print("   ------------------------");
  }
  // Here print the relevant parameters

  // Ready to iterate!
  for (unsigned int proto = 0; proto < Rparams.protocol_data.size(); proto++) {
    // Set defaults inside here
    set_protocol<3>(world, Rparams.protocol_data[proto]);

    // Do something to ensure all functions have same k value
    check_k(world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());

    // Create the active subspace (select which ground state orbitals to
    // calculate excitations from)
    // if(Rparams.e_window) select_active_subspace(world);

    if (proto == 0) {
      if (Rparams.restart) {
        if (world.rank() == 0) {
          print("   Restarting from file:", Rparams.restart_file);
        }
        load(world, Rparams.restart_file);
        check_k(
            world, Rparams.protocol_data[proto], FunctionDefaults<3>::get_k());
      } else {
        // Create trial functions by...
        // (Always creating (at least) twice the amount requested for
        // initial diagonalization)
        if (world.rank() == 0) print("\n   Creating trial functions.\n");
        if (Rparams.random) {
          // Random guess
          x_response = create_random_guess(world,
                                           2 * Rparams.states,
                                           Gparams.num_orbitals,
                                           Gparams.orbitals,
                                           Gparams.molecule);
        } else if (Rparams.nwchem != "") {
          // Virtual orbitals from NWChem
          x_response = create_nwchem_guess(world, 2 * Rparams.states);
        } else if (Rparams.guess_xyz) {
          // Use a symmetry adapted operator on ground state functions
          x_response = create_trial_functions2(
              world, Gparams.orbitals, Rparams.print_level);
        } else {
          x_response = create_trial_functions(
              world, 2 * Rparams.states, Gparams.orbitals, Rparams.print_level);
        }

        // Load balance
        // Only balancing on x-components. Smart?
        if (world.size() > 1) {
          // Start a timer
          if (Rparams.print_level >= 1) molresponse::start_timer(world);
          if (world.rank() == 0) print("");  // Makes it more legible

          LoadBalanceDeux<3> lb(world);
          for (size_t j = 0; j < Rparams.states; j++) {
            for (size_t k = 0; k < Gparams.num_orbitals; k++) {
              lb.add_tree(x_response[j][k], lbcost<double, 3>(1.0, 8.0), true);
            }
          }
          for (size_t j = 0; j < Gparams.num_orbitals; j++) {
            lb.add_tree(Gparams.orbitals[j], lbcost<double, 3>(1.0, 8.0), true);
          }
          FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

          if (Rparams.print_level >= 1)
            molresponse::end_timer(world, "Load balancing:");
        }

        // Project out groundstate from guesses
        QProjector<double, 3> projector(world, Gparams.orbitals);
        for (unsigned int i = 0; i < x_response.size(); i++)
          x_response[i] = projector(x_response[i]);

        // Ensure orthogonal guesses
        for (size_t i = 0; i < 2; i++) {
          molresponse::start_timer(world);
          // Orthog
          x_response = gram_schmidt(world, x_response);
          molresponse::end_timer(world, "orthog");

          molresponse::start_timer(world);
          // Normalize
          normalize(world, x_response);
          molresponse::end_timer(world, "normalize");
        }

        // Diagonalize guess
        if (world.rank() == 0)
          print(
              "\n   Iterating trial functions for an improved initial "
              "guess.\n");
        IterateGuess(world, x_response);
        // Sort
        sort(world, omega, x_response);
        // Basic output
        if (Rparams.print_level >= 1 and world.rank() == 0) {
          print("\n   Final initial guess excitation energies:");
          print(omega);
        }

        // Select lowest energy functions from guess
        x_response = select_functions(
            world, x_response, omega, Rparams.states, Rparams.print_level);

        // Initial guess for y are zero functions
        y_response =
            response_space(world, Rparams.states, Gparams.num_orbitals);
      }
    }

    // Now actually ready to iterate...
    Iterate(world);
  }

  // Plot the response function if desired
  if (Rparams.plot) {
    // Need to get densities first
    std::vector<real_function_3d> densities =
        transition_density(world, Gparams.orbitals, x_response, y_response);

    // For the instance where we don't plot all the orbitals
    std::vector<real_function_3d> plot_densities;
    for (size_t i : Rparams.plot_data) {
      plot_densities.push_back(densities[i]);
    }

    // Now plot
    if (world.rank() == 0) print("\n   Plotting response state densities.\n");
    if (Rparams.plot_L > 0.0)
      do_vtk_plots(world,
                   Rparams.plot_pts,
                   Rparams.plot_L,
                   0,
                   Rparams.plot_data.size(),
                   Gparams.molecule,
                   plot_densities,
                   "response-state");
    else
      do_vtk_plots(world,
                   Rparams.plot_pts,
                   Gparams.L,
                   0,
                   Rparams.plot_data.size(),
                   Gparams.molecule,
                   plot_densities,
                   "response-state");
  }

  // Print total time
  // Precision is set to 10 coming in, drop it to 2
  std::cout.precision(2);
  std::cout << std::fixed;

  // Get start time
  molresponse::end_timer(world, "total:");
}
