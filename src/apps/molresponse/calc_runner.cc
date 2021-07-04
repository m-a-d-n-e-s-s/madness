// Copyright 2021 Adrian Hurtado
#include <math.h>

#include <cstdint>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <utility>

#include "Plot_VTK.h"
#include "TDDFT.h"
#include "chem/NWChem.h"  // For nwchem interface
#include "chem/SCFOperators.h"
#include "chem/molecule.h"
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

// Masking function to switch from 0 to 1 smoothly at boundary
// Pulled from SCF.h
inline double mask1(double x) {
  /* Iterated first beta function to switch smoothly
     from 0->1 in [0,1].  n iterations produce 2*n-1
     zero derivatives at the end points. Order of polyn
     is 3^n.

     Currently use one iteration so that first deriv.
     is zero at interior boundary and is exactly representable
     by low order multiwavelet without refinement */

  x = (x * x * (3. - 2. * x));
  return x;
}

static double mask3(const coord_3d& ruser) {
  coord_3d rsim;
  user_to_sim(ruser, rsim);
  double x = rsim[0], y = rsim[1], z = rsim[2];
  double lo = 0.0625, hi = 1.0 - lo, result = 1.0;
  double rlo = 1.0 / lo;

  if (x < lo)
    result *= mask1(x * rlo);
  else if (x > hi)
    result *= mask1((1.0 - x) * rlo);
  if (y < lo)
    result *= mask1(y * rlo);
  else if (y > hi)
    result *= mask1((1.0 - y) * rlo);
  if (z < lo)
    result *= mask1(z * rlo);
  else if (z > hi)
    result *= mask1((1.0 - z) * rlo);

  return result;
}

template <std::size_t NDIM>
void TDDFT::set_protocol(World& world, double thresh) {
  size_t k;
  // Allow for imprecise conversion of threshold
  if (thresh >= 0.9e-2)
    k = 4;
  else if (thresh >= 0.9e-4)
    k = 6;
  else if (thresh >= 0.9e-6)
    k = 8;
  else if (thresh >= 0.9e-8)
    k = 10;
  else
    k = 12;

  // k defaults to make sense with thresh, override by providing k in
  // input file
  if (r_params.k() == -1) {
    FunctionDefaults<NDIM>::set_k(k);
  } else {
    FunctionDefaults<NDIM>::set_k(r_params.k());
  }

  // MolDFT sets all these, so copying
  FunctionDefaults<NDIM>::set_thresh(thresh);
  FunctionDefaults<NDIM>::set_refine(true);
  FunctionDefaults<NDIM>::set_initial_level(2);

  FunctionDefaults<NDIM>::set_autorefine(false);
  FunctionDefaults<NDIM>::set_apply_randomize(false);
  FunctionDefaults<NDIM>::set_project_randomize(false);
  GaussianConvolution1DCache<double>::map.clear();
  double safety = 0.1;
  vtol = FunctionDefaults<NDIM>::get_thresh() * safety;
  coulop = poperatorT(CoulombOperatorPtr(world, r_params.lo(), thresh));
  gradop = gradient_operator<double, 3>(world);
  // GaussianConvolution1DCache<double>::map.clear();//(TODO:molresponse-What
  // is this? Do i need it?)
  rho0 = make_ground_density(world, ground_orbitals);

  // Create the masking function
  mask = real_function_3d(real_factory_3d(world).f(mask3).initial_level(4).norefine());
  // dconv defaults to thresh*100, overrirde by providing dconv in input
  // file
  if (r_params.dconv_set() == false) {
    r_params.set_derived_value<double>("dconv", thresh * 100);
  }

  // Basic print
  if (world.rank() == 0) {
    print("\nSolving NDIM=",
          NDIM,
          " with thresh",
          thresh,
          "    k",
          FunctionDefaults<NDIM>::get_k(),
          "  dconv",
          std::max(thresh, r_params.dconv()),
          "\n");
  }
}

void TDDFT::check_k(World& world, double thresh, size_t k) {
  // Boolean to redo ground hamiltonian calculation if
  // ground state orbitals change
  bool redo = false;

  // Verify ground state orbitals have correct k
  if (FunctionDefaults<3>::get_k() != ground_orbitals[0].k()) {
    // Re-read orbitals from the archive (assuming
    // the archive has orbitals stored at a higher
    // k value than what was previously computed
    // with)
    g_params.read(world, r_params.archive());
    reconstruct(world, ground_orbitals);

    // Reset correct k (its set in g_params.read)
    FunctionDefaults<3>::set_k(k);

    // Project each ground state to correct k
    for (unsigned int i = 0; i < ground_orbitals.size(); i++)
      ground_orbitals[i] = project(ground_orbitals[i], FunctionDefaults<3>::get_k(), thresh, false);
    world.gop.fence();

    // Clean up a bit
    truncate(world, ground_orbitals);

    // Ground state orbitals changed, clear old hamiltonian
    redo = true;
  }

  // Recalculate ground state hamiltonian here
  if (redo or !hamiltonian.has_data()) {
    hamiltonian = CreateGroundHamiltonian(world, ground_orbitals, r_params.print_level());
  }

  // If we stored the potential, check that too
  if (r_params.store_potential()) {
    if (FunctionDefaults<3>::get_k() != stored_potential[0][0].k()) {
      // Project the potential into correct k
      for (unsigned int i = 0; i < stored_potential.size(); i++) {
        reconstruct(world, stored_potential[i]);
        for (unsigned int j = 0; j < stored_potential[0].size(); j++)
          stored_potential[i][j] = project(stored_potential[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
    }
    if (FunctionDefaults<3>::get_k() != stored_v_coul.k())
      stored_v_coul = project(stored_v_coul, FunctionDefaults<3>::get_k(), thresh, false);
    if (FunctionDefaults<3>::get_k() != stored_v_nuc.k())
      stored_v_nuc = project(stored_v_nuc, FunctionDefaults<3>::get_k(), thresh, false);
  }

  // Verify response functions have correct k
  if (Chi.X.size() != 0) {
    if (FunctionDefaults<3>::get_k() != Chi.X[0][0].k()) {
      // Project all x components into correct k
      for (unsigned int i = 0; i < Chi.X.size(); i++) {
        reconstruct(world, Chi.X[i]);
        for (unsigned int j = 0; j < Chi.X[0].size(); j++)
          Chi.X[i][j] = project(Chi.X[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      Chi.X.truncate_rf();

      // Do same for y components if applicable
      // (Always do this, as y will be zero
      //  and still used in doing DFT and TDA)
      // Project all y components into correct k
      for (unsigned int i = 0; i < Chi.Y.size(); i++) {
        reconstruct(world, Chi.Y[i]);
        for (unsigned int j = 0; j < Chi.Y[0].size(); j++)
          Chi.Y[i][j] = project(Chi.Y[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      Chi.Y.truncate_rf();
    }
  }
  // Verify response functions have correct k
  if (Chi.X.size() != 0) {
    if (FunctionDefaults<3>::get_k() != Chi.X[0][0].k()) {
      // Project all x components into correct k
      for (unsigned int i = 0; i < Chi.X.size(); i++) {
        reconstruct(world, Chi.X[i]);
        for (unsigned int j = 0; j < Chi.X[0].size(); j++)
          Chi.X[i][j] = project(Chi.X[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      Chi.X.truncate_rf();

      // Do same for y components if applicable
      // (Always do this, as y will be zero
      //  and still used in doing DFT and TDA)
      // Project all y components into correct k
      for (unsigned int i = 0; i < Chi.Y.size(); i++) {
        reconstruct(world, Chi.Y[i]);
        for (unsigned int j = 0; j < Chi.Y[0].size(); j++)
          Chi.Y[i][j] = project(Chi.Y[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      Chi.Y.truncate_rf();
    }
  }

  // Don't forget the mask function as well
  if (FunctionDefaults<3>::get_k() != mask.k()) {
    mask = project(mask, FunctionDefaults<3>::get_k(), thresh, false);
  }
  // Don't forget right hand side

  // Verify response functions have correct k
  if (PQ.X.size() != 0) {
    if (FunctionDefaults<3>::get_k() != PQ.X[0][0].k()) {
      // Project all x components into correct k
      for (unsigned int i = 0; i < Chi.X.size(); i++) {
        reconstruct(world, PQ.X[i]);
        for (unsigned int j = 0; j < PQ.X[0].size(); j++)
          PQ.X[i][j] = project(PQ.X[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      Chi.X.truncate_rf();

      // Do same for y components if applicable
      // (Always do this, as y will be zero
      //  and still used in doing DFT and TDA)
      // Project all y components into correct k
      for (unsigned int i = 0; i < PQ.Y.size(); i++) {
        reconstruct(world, PQ.Y[i]);
        for (unsigned int j = 0; j < PQ.Y[0].size(); j++)
          PQ.Y[i][j] = project(PQ.Y[i][j], FunctionDefaults<3>::get_k(), thresh, false);
        world.gop.fence();
      }
      Chi.Y.truncate_rf();
    }
  }
  // Make sure everything is done before leaving
  world.gop.fence();
}

// Creates random guess functions semi-intelligently(?)
X_space TDDFT::create_random_guess(World& world,
                                   size_t m,                                // m response states
                                   size_t n,                                // n ground states
                                   std::vector<real_function_3d>& grounds,  // guess should have size n
                                   Molecule& molecule) {
  // Basic output
  if (world.rank() == 0) print("   Using a random guess for initial response functions.\n");

  // Create empty container and add in randomness
  X_space f(world, m, n);
  f.X = add_randomness(world, f.X, 1e3);  // noise all over the world

  // Create and apply a centered gaussian on each atom so that the
  // randomness is localized around the atoms
  real_function_3d gaus = real_factory_3d(world);
  for (auto atom : molecule.get_atoms()) {
    real_function_3d x = real_factory_3d(world).functor(
        real_functor_3d(new GaussianGuess<3>(atom.get_coords(), 0.01, std::vector<int>{0, 0, 0})));
    gaus = gaus + x;
  }
  f = f * gaus;

  // Project out groundstate from guesses
  QProjector<double, 3> projector(world, grounds);
  for (unsigned int i = 0; i < f.num_states(); i++) f.X[i] = projector(f.X[i]);

  // Normalize
  normalize(world, f.X);

  return f;
}

// Creates random guess functions semi-intelligently(?)
std::vector<real_function_3d> TDDFT::create_random_guess(World& world,
                                                         size_t m,
                                                         std::vector<real_function_3d>& grounds,
                                                         Molecule& molecule) {
  // Basic output
  if (world.rank() == 0) print("   Using a random guess for initial response functions.");

  // Create empty container and add in randomness
  std::vector<real_function_3d> f = zero_functions_compressed<double, 3>(world, m);
  // create vector of m functions:w

  // Create and apply a centered gaussian on each atom so that the
  // randomness is localized around the atoms
  real_function_3d gaus = real_factory_3d(world);
  for (auto atom : molecule.get_atoms()) {
    real_function_3d x = real_factory_3d(world).functor(
        real_functor_3d(new GaussianGuess<3>(atom.get_coords(), 0.01, std::vector<int>{0, 0, 0})));
    gaus = gaus + x;
  }

  // Lambda function to add in noise
  auto lambda = [](const Key<3>& key, Tensor<double>& x) mutable {
    Tensor<double> y(x.size());
    y.fillrandom();
    y.scale(1e3);
    x = x + y;
  };

  // Go through each function in f_copy and add in random noise
  for (unsigned int i = 0; i < f.size(); i++) {
    // Add in random noise using rng and a the defined lambda function
    f[i].unaryop(lambda);

    // Apply mask to get boundary condition right
    f[i] = mask * f[i] * gaus;
  }

  // Project out groundstate from guesses
  QProjector<double, 3> projector(world, grounds);
  for (unsigned int i = 0; i < f.size(); i++) f[i] = projector(f[i]);

  // Normalize
  for (unsigned int i = 0; i < f.size(); i++) {
    double norm = f[i].norm2();
    f[i].scale(1.0 / norm);
  }

  return f;
}

// Creates an initial guess function from nwchem output files
X_space TDDFT::create_nwchem_guess(World& world, size_t m) {
  // Basic output
  if (world.rank() == 0) print("   Creating an initial guess from NWChem file", r_params.nwchem_dir());

  // Create empty containers
  response_space f;

  // Create the nwchem reader
  slymer::NWChem_Interface nwchem(r_params.nwchem_dir(), std::cout);

  // For parallel runs, silencing all but 1 slymer instance
  if (world.rank() != 0) {
    std::ostream dev_null(nullptr);
    nwchem.err = dev_null;
  }

  // Read in basis set
  nwchem.read(slymer::Properties::Basis);

  // Read in the molecular orbital coefficients, energies,
  // and occupancies
  nwchem.read(slymer::Properties::MOs | slymer::Properties::Energies | slymer::Properties::Occupancies);

  // Create the nwchem orbitals as madness functions
  std::vector<real_function_3d> temp1;
  for (auto basis : slymer::cast_basis<slymer::GaussianFunction>(nwchem.basis_set)) {
    // Get the center of gaussian as its special point
    std::vector<coord_3d> centers;
    coord_3d r;
    r[0] = basis.get().center[0];
    r[1] = basis.get().center[1];
    r[2] = basis.get().center[2];
    centers.push_back(r);

    // Now make the function
    temp1.push_back(FunctionFactory<double, 3>(world).functor(
        std::shared_ptr<FunctionFunctorInterface<double, 3>>(new slymer::Gaussian_Functor(basis.get(), centers))));

    // Let user know something is going on
    if (temp1.size() % 10 == 0 and world.rank() == 0) print("Created", temp1.size(), "functions.");
  }
  if (world.rank() == 0) print("Finished creating", temp1.size(), "functions.");

  // Normalize ao's
  madness::normalize(world, temp1);

  // Transform ao's now
  std::vector<real_function_3d> temp =
      madness::transform(world, temp1, nwchem.MOs, FunctionDefaults<3>::get_thresh(), true);

  // Now save the unoccupied orbitals
  std::vector<real_function_3d> temp2;
  size_t num_virt = 0;
  for (size_t i = 0; i < temp1.size(); i++) {
    if (nwchem.occupancies[i] == 0) {
      temp2.push_back(temp[i]);
      num_virt++;
    }
  }

  // Create as many vectors of functions as we can from these nwchem
  // virtual orbitals putting 1 virtual orbital from nwchem per vector
  for (size_t i = 0; i < std::max(m, num_virt); i++) {
    // Create the vector to add the new function to
    std::vector<real_function_3d> v1 = zero_functions_compressed<double, 3>(world, g_params.orbitals().size());

    // Put the "new" function into the vector
    v1[i % v1.size()] = temp2[i];

    // Add vector to return container
    f.push_back(v1);

    // See if we've made enough functions
    if (f.size() >= m) break;
  }
  if (world.rank() == 0) print("Created", f.size(), "guess functions from provided NWChem data.");

  // If not enough functions have been made, start adding symmetry adapted
  // functions
  size_t n = f.size();

  // If still not enough functions have been made, add in random guesses
  if (n < m) {
    // Tell user the bad news
    if (world.rank() == 0)
      print("\n   Only",
            n,
            "guess functions were provided by augmenting NWChem "
            "functions.\n   "
            "Augmenting with random functions.");

    // Create the random guess
    Molecule mol = g_params.molecule();
    X_space rand = create_random_guess(world, m - n, n, ground_orbitals, mol);

    // Add to vector of functions
    for (unsigned int i = 0; i < rand.num_states(); i++) f.push_back(rand.X[i]);
  }

  // Project out groundstate from guesses
  QProjector<double, 3> projector(world, ground_orbitals);
  for (unsigned int i = 0; i < f.size(); i++) f[i] = projector(f[i]);

  // Truncate and normalize
  f.truncate_rf();
  normalize(world, f);

  X_space trial(world, f.size(), n);
  trial.X = f;

  return trial;
}

// Creates potentials using the ResponsePotential object
// Potentials are modified in place
void TDDFT::create_all_potentials(World& world,
                                  response_space& x,
                                  response_space& x_gamma,
                                  response_space& x_V0,
                                  ResponsePotential& potentials,
                                  size_t print_level) {
  // Intermediaries
  response_space gammaJ, gammaK, groundJ, groundK;

  // Calc. coulomb like terms
  if (print_level >= 1) molresponse::start_timer(world);
  potentials.coulomb_terms(x, gammaK, groundJ);
  if (print_level >= 1) molresponse::end_timer(world, "Coulomb terms:");

  // Calc. exchange like terms
  if (print_level >= 1) molresponse::start_timer(world);
  potentials.exchange_terms(x, gammaJ, groundK);
  if (print_level >= 1) molresponse::end_timer(world, "Exchange terms:");

  // Assemble pieces together
  x_gamma = gammaJ - gammaK;
  x_V0 = groundJ - groundK;

  // Debugging output
  if (print_level >= 2) {
    // Coulomb
    if (world.rank() == 0) printf("   Coulomb Deriv matrix:\n");
    Tensor<double> temp = expectation(world, x, gammaJ);
    if (world.rank() == 0) print(temp);

    // Exchange or VXC
    if (r_params.xc() == "hf" and world.rank() == 0) printf("   Exchange Deriv matrix:\n");
    if (r_params.xc() != "hf" and world.rank() == 0) printf("   Negative of XC Deriv matrix:\n");
    temp = expectation(world, x, gammaK);
    if (world.rank() == 0) print(temp);

    // Total Gamma
    if (world.rank() == 0) printf("   Gamma matrix:\n");
    temp = expectation(world, x, x_gamma);
    if (world.rank() == 0) print(temp);

    // Coulomb (ground)
    if (world.rank() == 0) printf("   Coulomb + Nuclear potential matrix:\n");
    temp = expectation(world, x, groundJ);
    if (world.rank() == 0) print(temp);

    // Exchange or VXC (ground)
    if (r_params.xc() == "hf" and world.rank() == 0) printf("   Exchange potential matrix:\n");
    if (r_params.xc() != "hf" and world.rank() == 0) printf("   XC potential matrix:\n");
    temp = expectation(world, x, groundK);
    if (world.rank() == 0) print(temp);
    if (world.rank() == 0) printf("   Total Potential Energy matrix:\n");
    temp = expectation(world, x, x_V0);
    if (world.rank() == 0) print(temp);
  }
}

// Main function, makes sure everything happens in correcct order
// Solves for response components
void TDDFT::solve_excited_states(World& world) {
  // Get start time
  molresponse::start_timer(world);

  // Plotting input orbitals
  if (r_params.plot_initial()) {
    if (world.rank() == 0) print("\n   Plotting ground state densities.\n");
    if (r_params.plot_L() > 0.0)
      do_vtk_plots(world,
                   r_params.plot_pts(),
                   r_params.plot_L(),
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
    print("\n\n     Excited State Calculation");
    print("   ------------------------");
  }
  // Here print the relevant parameters

  // Ready to iterate!
  for (unsigned int proto = 0; proto < r_params.protocol().size(); proto++) {
    // Set defaults inside here
    set_protocol<3>(world, r_params.protocol()[proto]);

    // Do something to ensure all functions have same k value
    check_k(world, r_params.protocol()[proto], FunctionDefaults<3>::get_k());

    if (proto == 0) {
      if (r_params.restart()) {
        if (world.rank() == 0) {
          print("   Restarting from file:", r_params.restart_file());
        }
        load(world, r_params.restart_file());
        check_k(world, r_params.protocol()[proto], FunctionDefaults<3>::get_k());
      } else {
        X_space trial(world, 2 * r_params.n_states(), r_params.num_orbitals());
        // Create trial functions by...
        // (Always creating (at least) twice the amount requested for
        // initial diagonalization)
        if (world.rank() == 0) print("\n   Creating trial functions.\n");
        if (r_params.random()) {
          // Random guess
          trial =
              create_random_guess(world, 2 * r_params.n_states(), r_params.num_orbitals(), ground_orbitals, molecule);
        } else if (r_params.nwchem()) {
          // Virtual orbitals from NWChem
          trial = create_nwchem_guess(world, 2 * r_params.n_states());
        } else if (r_params.guess_xyz()) {
          // Use a symmetry adapted operator on ground state functions
          trial = create_trial_functions2(world, ground_orbitals, r_params.num_orbitals());
        } else {
          trial = create_trial_functions(world, 2 * r_params.n_states(), ground_orbitals, r_params.num_orbitals());
        }

        // Load balance
        // Only balancing on x-components. Smart?
        if (world.size() > 1) {
          // Start a timer
          if (r_params.num_orbitals() >= 1) molresponse::start_timer(world);
          if (world.rank() == 0) print("");  // Makes it more legible

          LoadBalanceDeux<3> lb(world);
          for (size_t j = 0; j < r_params.n_states(); j++) {
            for (size_t k = 0; k < r_params.num_orbitals(); k++) {
              lb.add_tree(trial.X[j][k], lbcost<double, 3>(1.0, 8.0), true);
            }
          }
          for (size_t j = 0; j < r_params.num_orbitals(); j++) {
            lb.add_tree(ground_orbitals[j], lbcost<double, 3>(1.0, 8.0), true);
          }
          FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

          if (r_params.num_orbitals() >= 1) molresponse::end_timer(world, "Load balancing:");
        }

        // Project out groundstate from guesses
        QProjector<double, 3> projector(world, ground_orbitals);
        for (unsigned int i = 0; i < trial.X.size(); i++) trial.X[i] = projector(trial.X[i]);

        // Ensure orthogonal guesses
        for (size_t i = 0; i < 2; i++) {
          molresponse::start_timer(world);
          // Orthog
          trial.X = gram_schmidt(world, trial.X);
          molresponse::end_timer(world, "orthog");

          molresponse::start_timer(world);
          // Normalize
          normalize(world, trial.X);
          molresponse::end_timer(world, "normalize");
        }

        // Diagonalize guess
        if (world.rank() == 0)
          print(
              "\n   Iterating trial functions for an improved initial "
              "guess.\n");
        iterate_guess(world, trial);
        // Sort
        sort(world, omega, trial.X);
        // Basic output
        if (r_params.num_orbitals() >= 1 and world.rank() == 0) {
          print("\n   Final initial guess excitation energies:");
          print(omega);
        }
        // Chi = X_space(world, r_params.n_states(), r_params.num_orbitals());
        // Select lowest energy functions from guess
        Chi.X = select_functions(world, trial.X, omega, r_params.n_states(), r_params.num_orbitals());
        Chi.Y = response_space(world, r_params.n_states(), r_params.num_orbitals());

        trial.clear();
        // Initial guess for y are zero functions
      }
    }

    // Now actually ready to iterate...
    iterate_excited(world, Chi);
    if (r_params.save()) save(world, r_params.save_file());
  }

  // Plot the response function if desired
  if (r_params.plot()) {
    // Need to get densities first
    std::vector<real_function_3d> densities = transition_density(world, ground_orbitals, Chi.X, Chi.Y);

    // For the instance where we don't plot all the orbitals
    std::vector<real_function_3d> plot_densities;
    for (size_t i : r_params.plot_data()) {
      plot_densities.push_back(densities[i]);
    }

    // Now plot
    if (world.rank() == 0) print("\n   Plotting response state densities.\n");
    if (r_params.plot_L() > 0.0)
      do_vtk_plots(world,
                   r_params.plot_pts(),
                   r_params.plot_L(),
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

  // Print total time
  // Precision is set to 10 coming in, drop it to 2
  std::cout.precision(2);
  std::cout << std::fixed;

  // Get start time
  molresponse::end_timer(world, "total:");
}

void TDDFT::solve_response_states(World& world) {
  molresponse::start_timer(world);
  // Warm and fuzzy
  if (world.rank() == 0) {
    print("\n\n    Response Calculation");
    print("   ------------------------");
  }
  for (unsigned int proto = 0; proto < r_params.protocol().size(); proto++) {
    // Set defaults inside here
    // default value of
    set_protocol<3>(world, r_params.protocol()[proto]);

    check_k(world, r_params.protocol()[proto], FunctionDefaults<3>::get_k());
    // Do something to ensure all functions have same k value

    if (r_params.dipole()) {
      if (world.rank() == 0) print("creating dipole property operator");
      p = DipoleVector(world);
    } else if (r_params.nuclear()) {
      Molecule molecule = g_params.molecule();
      if (world.rank() == 0) print("creating nuclear property operator");
      p = NuclearVector(world, molecule);
    }
    if (proto == 0) {
      if (r_params.restart()) {
        if (world.rank() == 0) print("   Initial guess from file:", r_params.restart_file());
        load(world, r_params.restart_file());
        check_k(world, r_params.protocol()[proto], FunctionDefaults<3>::get_k());

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
    print("Preiteration Information");
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

    // Now actually ready to iterate...
    // iterate_freq(world);
    iterate_freq2(world);
    // IterateFrequencyResponse(world, P, Q);
  }  // end for --finished reponse density

  // Have response function, now calculate polarizability for this axis
  // polarizability(world, polar_tensor);
  // PrintPolarizabilityAnalysis(world, polar_tensor, omega);

  // Print total time
  // Precision is set to 10 coming in, drop it to 2
  std::cout.precision(2);
  std::cout << std::fixed;

  // Get start time
  molresponse::end_timer(world, "total:");
}  // end compute frequency response
