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
#include "../chem/NWChem.h"  // For nwchem interface
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
          Chi.X = create_random_guess(world,
                                           2 * Rparams.states,
                                           Gparams.num_orbitals,
                                           Gparams.orbitals,
                                           Gparams.molecule);
        } else if (Rparams.nwchem != "") {
          // Virtual orbitals from NWChem
          Chi.X = create_nwchem_guess(world, 2 * Rparams.states);
        } else if (Rparams.guess_xyz) {
          // Use a symmetry adapted operator on ground state functions
          Chi.X = create_trial_functions2(
              world, Gparams.orbitals, Rparams.print_level);
        } else {
          Chi.X = create_trial_functions(
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
              lb.add_tree(Chi.X[j][k], lbcost<double, 3>(1.0, 8.0), true);
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
        for (unsigned int i = 0; i < Chi.X.size(); i++)
          Chi.X[i] = projector(Chi.X[i]);

        // Ensure orthogonal guesses
        for (size_t i = 0; i < 2; i++) {
          molresponse::start_timer(world);
          // Orthog
          Chi.X = gram_schmidt(world, Chi.X);
          molresponse::end_timer(world, "orthog");

          molresponse::start_timer(world);
          // Normalize
          normalize(world, Chi.X);
          molresponse::end_timer(world, "normalize");
        }

        // Diagonalize guess
        if (world.rank() == 0)
          print(
              "\n   Iterating trial functions for an improved initial "
              "guess.\n");
        IterateGuess(world, Chi.X);
        // Sort
        sort(world, omega, Chi.X);
        // Basic output
        if (Rparams.print_level >= 1 and world.rank() == 0) {
          print("\n   Final initial guess excitation energies:");
          print(omega);
        }
        //Chi = X_space(world, Rparams.states, Gparams.num_orbitals);
        // Select lowest energy functions from guess
        Chi.X= select_functions(
            world, Chi.X, omega, Rparams.states, Rparams.print_level);
		Chi.Y=response_space(world,Rparams.states,Gparams.num_orbitals);
        // Initial guess for y are zero functions
      }
    }

    // Now actually ready to iterate...
    Iterate(world, Chi);
  }

  // Plot the response function if desired
  if (Rparams.plot) {
    // Need to get densities first
    std::vector<real_function_3d> densities =
        transition_density(world, Gparams.orbitals, Chi.X, Chi.Y);

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
