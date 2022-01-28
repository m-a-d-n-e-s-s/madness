//
// Created by adrianhurtado on 1/28/22.
//

#include "ExcitedResponse.hpp"


void ExcitedResponse::initialize(World &world) {

    X_space trial(world, 2 * r_params.num_states(), r_params.num_orbitals());
    // make the trial functions
    if (world.rank() == 0) print("\n   Creating trial functions.\n");
    if (r_params.random()) {
        trial = make_random_trial(world, 2 * r_params.num_states());
    } else if (r_params.nwchem()) {
        // Virtual orbitals from NWChem
        trial = make_nwchem_trial(world, 2 * r_params.num_states());
    } else if (r_params.guess_xyz()) {
        // Use a symmetry adapted operator on ground state functions
        trial = create_trial_functions2(world, ground_orbitals, r_params.num_orbitals());
    } else {
        trial = create_trial_functions(world, 2 * r_params.num_states(), ground_orbitals,
                                       r_params.num_orbitals());
    }

    if (world.size() > 1) {
        // Start a timer
        if (r_params.num_orbitals() >= 1) molresponse::start_timer(world);
        if (world.rank() == 0) print("");// Makes it more legible

        LoadBalanceDeux<3> lb(world);
        for (size_t j = 0; j < r_params.num_states(); j++) {
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

    // Project out ground state from guesses
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
        print("\n   Iterating trial functions for an improved initial "
              "guess.\n");
    iterate_guess(world, trial);
    // Sort
    sort(world, omega, trial.X);
    // Basic output
    if (r_params.num_orbitals() >= 1 and world.rank() == 0) {
        print("\n   Final initial guess excitation energies:");
        print(omega);
    }
    // Chi = X_space(world, r_params.num_states(), r_params.num_orbitals());
    // Select lowest energy functions from guess
    Chi.X = select_functions(world, trial.X, omega, r_params.num_states(), r_params.num_orbitals());
    Chi.Y = response_space(world, r_params.num_states(), r_params.num_orbitals());

    trial.clear();
    // Initial guess for y are zero functions
}
// Creates random guess functions semi-intelligently(?)
/// & creates random guesses for nu
/// \param world
/// \param num_states
/// \return
X_space ExcitedResponse::make_random_trial(World &world, size_t m) const {
    // Basic output
    if (world.rank() == 0) print("   Using a random guess for initial response functions.\n");
    size_t n = r_params.num_orbitals();
    // Create empty container and add in randomness
    X_space f(world, m, n);
    f.X = add_randomness(world, f.X, 1e3);// noise all over the world
    f.X = mask * f.X;                     // make sure you mask after you add noise to the world

    // Create and apply a centered gaussian on each atom so that the
    // randomness is localized around the atoms
    real_function_3d gaus = real_factory_3d(world);
    for (auto atom: molecule.get_atoms()) {
        real_function_3d x = real_factory_3d(world).functor(real_functor_3d(
                new GaussianGuess<3>(atom.get_coords(), 0.01, std::vector<int>{0, 0, 0})));
        gaus = gaus + x;
    }
    f = f * gaus;

    // Project out groundstate from guesses
    QProjector<double, 3> projector(world, ground_orbitals);
    for (unsigned int i = 0; i < f.num_states(); i++) f.X[i] = projector(f.X[i]);

    // Normalize
    normalize(world, f.X);

    return f;
}


// Creates an initial guess function from nwchem output files
X_space ExcitedResponse::make_nwchem_trial(World &world, size_t m) const {
    // Basic output
    if (world.rank() == 0)
        print("   Creating an initial guess from NWChem file", r_params.nwchem_dir());

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
    nwchem.read(slymer::Properties::MOs | slymer::Properties::Energies |
                slymer::Properties::Occupancies);

    // Create the nwchem orbitals as madness functions
    std::vector<real_function_3d> temp1;
    for (auto basis: slymer::cast_basis<slymer::GaussianFunction>(nwchem.basis_set)) {
        // Get the center of gaussian as its special point
        std::vector<coord_3d> centers;
        coord_3d r;
        r[0] = basis.get().center[0];
        r[1] = basis.get().center[1];
        r[2] = basis.get().center[2];
        centers.push_back(r);

        // Now make the function
        temp1.push_back(FunctionFactory<double, 3>(world).functor(
                std::shared_ptr<FunctionFunctorInterface<double, 3>>(
                        new slymer::Gaussian_Functor(basis.get(), centers))));

        // Let user know something is going on
        if (temp1.size() % 10 == 0 and world.rank() == 0)
            print("Created", temp1.size(), "functions.");
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
        std::vector<real_function_3d> v1 =
                zero_functions_compressed<double, 3>(world, ground_orbitals.size());

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
            print("\n   Only", n,
                  "guess functions were provided by augmenting NWChem "
                  "functions.\n   "
                  "Augmenting with random functions.");

        // Create the random guess
        Molecule mol = molecule;
        X_space rand = make_random_trial(world, m - n);

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
// Returns initial guess functions as
// ground MO * solid harmonics
X_space ExcitedResponse::create_trial_functions(World &world, size_t k) {
    // Get size
    print("In create trial functions");
    auto n = r_params.num_orbitals();

    // Create solid harmonics such that num. solids * num. orbitals > k.
    // The total number of solid harmonics that exist up to level n is
    // (n+1)^2 (because we count from zero)
    // Always do at least 8 (through the d orbital angular momentum functions,
    // minus )
    std::map<std::vector<int>, real_function_3d> solids =
            solid_harmonics(world, std::max(2.0, ceil(sqrt(k / n) - 1)));

    // Useful info.
    if (world.rank() == 0) print("   Created", solids.size(), "solid harmonics.\n");

    // Container to return
    response_space trials_X;

    // Counter for number of trials created
    size_t count = 0;

    // Multiply each solid harmonic onto a ground state orbital
    for (size_t i = 0; i < n; i++) {
        // For each solid harmonic
        for (const auto &key: solids) {
            // Temp zero functions
            std::vector<real_function_3d> temp =
                    zero_functions_compressed<double, 3>(world, static_cast<int>(n));

            // Create one non-zero function and add to trials
            temp[count % n] = key.second * ground_orbitals[n - count % n - 1];
            trials_X.push_back(temp);
            count++;
        }

        // Stop when we first get beyond k components
        if (count >= k) break;
    }

    // Debugging output
    if (r_params.print_level() >= 2) {
        if (world.rank() == 0) print("   Norms of guess functions:");
        print_norms(world, trials_X);
    }

    // Truncate
    madness::truncate(world, trials_X, madness::FunctionDefaults<3>::get_thresh(), true);

    X_space trials(world, count, n);
    trials.X = trials_X.copy();
    trials_X.clear();

    // Done
    return trials;
}

// Returns initial guess functions as
// ground MO * <x,y,z>
X_space ExcitedResponse::create_trial_functions2(World &world) {
    // Get size
    size_t n = ground_orbitals.size();
    size_t directions = 3;
    // (n+1)^2 (because we count from zero)

    // make a vector of xyz functions {x,y,z}
    auto xyz = make_xyz_functions(world);
    // create 3 x n orbital functions
    std::vector<vector<real_function_3d>> functions;

    //

    for (const auto &d: xyz) {
        vector_real_function_3d temp;
        for (const auto &orb: ground_orbitals) { temp.push_back(d * orb); }
        functions.push_back(temp);
    }

    // Container to return
    size_t count = 0;

    X_space trials(world, 3 * n * n, n);
    for (size_t i = 0; i < n; i++) {
        for (size_t d = 0; d < directions; d++) {
            for (size_t o = 0; o < n; o++) {
                //        trials[i + j + o][o] = functions[i][j];
                trials.X[count][o] = copy(functions.at(d).at(o));
                count++;
            }
        }
    }
    // The above generates response function as follows
    // all functions start off as zeros
    // 1  [x1 0 0 ]
    // 2  [0 x1 0 ]
    // 3  [0 0 x1 ]
    // 4  [y1 0 0 ]
    // 5  [0 y1 0 ]
    // 6  [0 0 y1 ]
    // 7  [z1 0 0 ]
    // 8  [0 z1 0 ]
    // 9  [0 0 z1 ]
    // 10 [x2 0 0 ]
    // 11 [0 x2 0 ]
    // 12 [0 0 x2 ]
    // 13 [y2 0 0 ]
    // 14 [0 y2 0 ]
    // 15 [0 0 y2 ]
    // 16 [z2 0 0 ]
    // 17 [0 z2 0 ]
    // 18 [0 0 z2 ]
    // 19 [x3 0 0 ]
    // 20 [0 x3 0 ]
    // 21 [0 0 x3 ]
    // 22 [y3 0 0 ]
    // 23 [0 y3 0 ]
    // 24 [0 0 y3 ]
    // 25 [z3 0 0 ]
    // 26 [0 z3 0 ]
    // 27 [0 0 z3 ]
    // for each orbital for each direction
    // Counter for number of trials created
    // Multiply each solid harmonic onto a ground state orbital
    //  for each orbital we

    // For each solid harmonic
    // Temp zero functions

    // Debugging output
    if (r_params.print_level() >= 2) {
        if (world.rank() == 0) print("   Norms of guess functions:");
        print_norms(world, trials.X);
    }

    // Truncate
    madness::truncate(world, trials.X);

    // Done
    return trials;
}
// Simplified iterate scheme for guesses
void ExcitedResponse::iterate_trial(World &world, X_space &guesses) const {
    // Variables needed to iterate
    size_t iteration = 0;// Iteration counter
    QProjector<double, 3> projector(world,
                                    ground_orbitals);// Projector to project out ground state
    size_t m = r_params.num_states();                // Number of excited states
    size_t n = r_params.num_orbitals();              // Number of ground state orbitals
    Tensor<double> x_shifts;                         // Holds the shifted energy values
    response_space bsh_resp(world, m, n);            // Holds wave function corrections
    response_space V;                                // Holds V^0 applied to response functions
    response_space shifted_V;// Holds the shifted V^0 applied to response functions
    Tensor<double> S;        // Overlap matrix of response components for x states
    real_function_3d v_xc;   // For TDDFT

    auto xc = make_xc_operator(world);

    // Useful to have
    response_space zeros(world, m, n);

    // Now to iterate
    while (iteration < r_params.guess_max_iter()) {
        // Start a timer for this iteration
        molresponse::start_timer(world);
        //
        size_t N0 = guesses.X.size();
        size_t Ni = N0;


        // Basic output
        if (r_params.print_level() >= 1) {
            if (world.rank() == 0)
                printf("\n   Guess Iteration %d at time %.1fs\n", static_cast<int>(iteration),
                       wall_time());
            if (world.rank() == 0) print(" -------------------------------------");
        }

        // Load balance
        // Only balancing on x-components. Smart?
        if (world.size() > 1 && ((iteration < 2) or (iteration % 5 == 0)) and

            iteration != 0) {
            // Start a timer
            if (r_params.print_level() >= 1) molresponse::start_timer(world);
            if (world.rank() == 0) print("");// Makes it more legible

            LoadBalanceDeux<3> lb(world);
            for (size_t j = 0; j < n; j++) {
                for (size_t k = 0; k < r_params.num_states(); k++) {
                    lb.add_tree(guesses.X[k][j], lbcost<double, 3>(1.0, 8.0), true);
                }
            }
            FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

            if (r_params.print_level() >= 1) molresponse::end_timer(world, "Load balancing:");
        }

        // compute rho_omega
        rho_omega = transition_densityTDA(world,  guesses.X);
        // Project out ground state
        for (size_t i = 0; i < Ni; i++) guesses.X[i] = projector(guesses.X[i]);

        // Truncate before doing expensive things
        guesses.X.truncate_rf();

        // Normalize after projection
        if (r_params.tda()) normalize(world, guesses.X);

        // (TODO why not normalize if not tda)
        // compute Y = false
        X_space Lambda_X = Compute_Lambda_X(world, guesses, xc, "tda");

        deflateGuesses(world, guesses, Lambda_X, S, omega, iteration, m);
        // Debugging output

        // Ensure right number of omegas
        if (size_t(omega.dim(0)) != Ni) {
            if (world.rank() == 0)
                print("\n   Adding", Ni - omega.dim(0),
                      "eigenvalue(s) (counters subspace size "
                      "reduction in "
                      "diagonalizatoin).");
            Tensor<double> temp(Ni);
            temp(Slice(0, omega.dim(0) - 1)) = omega;
            for (size_t i = omega.dim(0); i < Ni; i++) temp[i] = 2.5 * i;
            omega = copy(temp);
        }

        // Basic output
        if (r_params.print_level() >= 1 and world.rank() == 0) {
            print("\n   Excitation Energies:");
            print("gi=", iteration, " roots: ", omega);
        }

        // Only do BSH if not the last iteration
        if (iteration + 1 < r_params.guess_max_iter()) {
            //  Calculates shifts needed for potential / energies
            //  If none needed, the zero tensor is returned
            x_shifts = create_shift(world, ground_energies, omega, r_params.print_level(), "x");

            X_space theta_X = Compute_Theta_X(world, guesses, xc, "tda");
            theta_X.X = apply_shift(world, x_shifts, theta_X.X, guesses.X);
            theta_X.X = theta_X.X * -2;
            theta_X.X.truncate_rf();

            // Construct BSH operators
            std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_x_operators =
                    create_bsh_operators(world, x_shifts, ground_energies, omega, r_params.lo(),
                                         FunctionDefaults<3>::get_thresh());

            // Apply BSH and get updated components
            if (r_params.print_level() >= 1) molresponse::start_timer(world);
            bsh_resp = apply(world, bsh_x_operators, theta_X.X);
            if (r_params.print_level() >= 1) molresponse::end_timer(world, "Apply BSH:");

            // Project out ground state
            for (size_t i = 0; i < Ni; i++) bsh_resp[i] = projector(bsh_resp[i]);
            // Save new components
            guesses.X = bsh_resp;
            // Apply mask
            for (size_t i = 0; i < Ni; i++) guesses.X[i] = mask * guesses.X[i];
        }

        // Ensure orthogonal guesses
        for (size_t i = 0; i < 2; i++) {
            molresponse::start_timer(world);
            // Orthog
            guesses.X = gram_schmidt(world, guesses.X);
            molresponse::end_timer(world, "orthog");

            molresponse::start_timer(world);
            // Normalize
            normalize(world, guesses.X);
            molresponse::end_timer(world, "normalize");
        }

        // Update counter
        iteration += 1;
        // Done with the iteration.. truncate
        guesses.X.truncate_rf();

        // Basic output
        if (r_params.print_level() >= 1) {//
            molresponse::end_timer(world, " This iteration:");
        }
    }
}// Done with iterate gues
 // Simplified iterate scheme for guesses