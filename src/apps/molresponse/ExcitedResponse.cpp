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
        trial = create_trial_functions2(world);
    } else {
        trial = create_trial_functions(world, 2 * r_params.num_states());
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
    iterate_trial(world, trial);
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
X_space ExcitedResponse::create_trial_functions(World &world, size_t k) const {
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
X_space ExcitedResponse::create_trial_functions2(World &world) const {
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
void ExcitedResponse::iterate_trial(World &world, X_space &guesses) {
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
        auto rho_omega = transition_densityTDA(world, ground_orbitals, guesses.X);
        // Project out ground state
        for (size_t i = 0; i < Ni; i++) guesses.X[i] = projector(guesses.X[i]);

        // Truncate before doing expensive things
        guesses.X.truncate_rf();

        // Normalize after projection
        if (r_params.tda()) normalize(world, guesses.X);

        // (TODO why not normalize if not tda)
        // compute Y = false
        auto xc = make_xc_operator(world);
        X_space Lambda_X = compute_lambda_X(world, guesses, xc, "tda");

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
            x_shifts = create_shift(world, ground_energies, omega, "x");

            X_space theta_X = compute_theta_X(world, guesses, xc, "tda");
            theta_X.X = apply_shift(world, x_shifts, theta_X.X, guesses.X);
            theta_X.X = theta_X.X * -2;
            theta_X.X.truncate_rf();

            // Construct BSH operators
            auto bsh_x_operators =
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

/**
 * @brief Diagonalize AX=SX omega
 *
 * @param world
 * @param S
 * @param old_S
 * @param old_A
 * @param x_response
 * @param old_x_response
 * @param ElectronResponses
 * @param OldElectronResponses
 * @param omega
 * @param iteration
 * @param m
 */

void ExcitedResponse::deflateGuesses(World &world, X_space &Chi, X_space &Lambda_X,
                                     Tensor<double> &S, Tensor<double> &omega, size_t &iteration,
                                     size_t &m) const {
    // XX =Omega XAX
    S = response_space_inner(Chi.X, Chi.X);
    Tensor<double> XAX = response_space_inner(Chi.X, Lambda_X.X);

    // Debugging output
    if (r_params.print_level() >= 2 and world.rank() == 0) {
        print(" Guess  Overlap matrix:");
        print(S);
        print(" Guess  XAX matrix:");
        print(XAX);
    }
    // Just to be sure dimensions work out, clear omega
    omega.clear();
    diagonalizeFockMatrix(world, Chi, Lambda_X, omega, XAX, S, FunctionDefaults<3>::get_thresh());
}

void ExcitedResponse::deflateTDA(World &world, X_space &Chi, X_space &old_Chi, X_space &Lambda_X,
                                 X_space &old_Lambda_X, Tensor<double> &S, Tensor<double> old_S,
                                 Tensor<double> old_A, Tensor<double> &omega, size_t &iteration,
                                 size_t &m) {
    S = response_space_inner(Chi.X, Chi.X);
    Tensor<double> XAX = response_space_inner(Chi.X, Lambda_X.X);

    // Debugging output
    if (r_params.print_level() >= 2 and world.rank() == 0) {
        print("   Overlap matrix:");
        print(S);
    }

    // Augment S_x, A_x, x_gamma, x_response, V_x_response and x_gamma
    // if using a larger subspace and not iteration zero (TODO ---Gotta
    // look at this and make sure it uses my new functions molresponse )
    // by default r_params.larger_subspace() = 0 therefore never uses this
    if (iteration < r_params.larger_subspace() and iteration > 0) {
        print("Using augmented subspace");
        augment(world, Chi, old_Chi, Lambda_X, old_Lambda_X, S, XAX, old_S, old_A,
                r_params.print_level());
    }

    // Solve Ax = Sxw
    // Just to be sure dimensions work out, clear omega
    omega.clear();
    diagonalizeFockMatrix(world, Chi, Lambda_X, omega, XAX, S, FunctionDefaults<3>::get_thresh());

    // If larger subspace, need to "un-augment" everything
    if (iteration < r_params.larger_subspace()) {
        print("Unaugmenting subspace");
        unaugment(world, Chi, old_Chi, Lambda_X, old_Lambda_X, omega, S, XAX, old_S, old_A,
                  r_params.num_states(), iteration, r_params.print_level());
    }
}

void ExcitedResponse::deflateFull(World &world, X_space &Chi, X_space &old_Chi, X_space &Lambda_X,
                                  X_space &old_Lambda_X, Tensor<double> &S, Tensor<double> old_S,
                                  Tensor<double> old_A, Tensor<double> &omega, size_t &iteration,
                                  size_t &m) {
    // Debugging output
    Tensor<double> A;

    if (iteration < r_params.larger_subspace() and iteration > 0) {
        print("Entering Augment Full");
        augment_full(world, Chi, old_Chi, Lambda_X, old_Lambda_X, S, A, old_S, old_A,
                     r_params.print_level());
        // computes Augments A and S

    } else {
        S = response_space_inner(Chi.X, Chi.X) - response_space_inner(Chi.Y, Chi.Y);
        if (world.rank() == 0 && (r_params.print_level() >= 10)) {
            print("\n   Overlap Matrix:");
            print(S);
        }
        X_space Chi_copy = Chi.copy();
        Chi_copy.truncate();
        A = inner(Chi_copy, Lambda_X);
        if (world.rank() == 0 && (r_params.print_level() >= 10)) {
            print("\n   Lambda Matrix:");
            print(A);
        }
    }

    omega.clear();

    Tensor<double> U = diagonalizeFullResponseMatrix(world, Chi, Lambda_X, omega, S, A,
                                                     FunctionDefaults<3>::get_thresh(),
                                                     r_params.print_level());

    if (iteration < r_params.larger_subspace() and iteration > 0) {
        print("Entering Unaugment Full");
        unaugment_full(world, Chi, old_Chi, Lambda_X, old_Lambda_X, omega, S, A, old_S, old_A,
                       r_params.num_states(), iteration, r_params.print_level());
    } else {
        old_Chi = Chi.copy();
        old_Lambda_X = Lambda_X.copy();
    }
}
void ExcitedResponse::augment(World &world, X_space &Chi, X_space &old_Chi, X_space &Lambda_X,
                              X_space &last_Lambda_X, Tensor<double> &S, Tensor<double> &A,
                              Tensor<double> &old_S, Tensor<double> &old_A, size_t print_level) {
    // Basic output
    if (print_level >= 1) molresponse::start_timer(world);

    // Get sizes
    size_t m = Chi.X.size();
    // Create work space, will overwrite S and A in the end
    Tensor<double> temp_S(2 * m, 2 * m);
    Tensor<double> temp_A(2 * m, 2 * m);
    /**
   * @brief Need to create off diagonal blocks of A
   *  A=
   *  [xAx      xAx_old ]
   *  [xoldAx  xoldAxold]
   *
   */
    // Calculate correct inner products of upper off diagonal

    Tensor<double> off = response_space_inner(Chi.X, last_Lambda_X.X);
    temp_A(Slice(0, m - 1), Slice(m, 2 * m - 1)) = copy(off);// top right
    // Now for lower off diagonal
    off = response_space_inner(old_Chi.X, Lambda_X.X);
    temp_A(Slice(m, 2 * m - 1), Slice(0, m - 1)) = copy(off);      // bottom left
    temp_A(Slice(0, m - 1), Slice(0, m - 1)) = copy(A);            // xAx top left
    temp_A(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) = copy(old_A);// xoldAxold bottom right
    // Debugging output
    if (print_level >= 2 and world.rank() == 0) {
        print("   Before symmeterizing A:");
        print(temp_A);
    }
    // Save temp_A as A_x
    // Need to symmeterize A as well (?)
    A = 0.5 * (temp_A + transpose(temp_A));
    /**
   * @brief Creating S
   * S= [<x|x>    <x|xold>   ]
   *    [<xold|x> <xold|xold>]
   */
    // Now create upper off diagonal block of S
    off = expectation(world, Chi.X, old_Chi.X);
    // Use slicing to put in correct spot
    temp_S(Slice(0, m - 1), Slice(m, 2 * m - 1)) = copy(off);// top right <x|xold>
    // Now the lower off diagonal block
    // (Go ahead and cheat and use the transpose...)
    off = transpose(off);// just transpose <xold|x>
    // Use slicing to put in correct spot
    temp_S(Slice(m, 2 * m - 1), Slice(0, m - 1)) = copy(off);// bottom right <xold|x>
    // Put together the rest of S
    temp_S(Slice(0, m - 1), Slice(0, m - 1)) = copy(S);            // top left <x|x>
    temp_S(Slice(m, 2 * m - 1), Slice(m, 2 * m - 1)) = copy(old_S);//<xold|xold>
    // Save temp_S as S_x
    S = copy(temp_S);
    // Add in old vectors to current vectors for the appropriate ones
    // Augment the vectors step
    for (size_t i = 0; i < m; i++) {
        Chi.X.push_back(old_Chi.X[i]);
        Lambda_X.X.push_back(last_Lambda_X.X[i]);
    }

    // End the timer
    if (print_level >= 1) molresponse::end_timer(world, "Aug. resp. matrix:");

    // Debugging output
    if (print_level >= 2 and world.rank() == 0) {
        print("\n   Augmented response matrix:");
        print(A);
    }

    // Debugging output
    if (print_level >= 2 and world.rank() == 0) {
        print("   Augmented overlap matrix:");
        print(S);
    }

    // SUPER debugging
    if (print_level >= 3) {
        if (world.rank() == 0) print("   Calculating condition number of aug. response matrix");
    }
}

// If using a larger subspace to diagonalize in, this will put everything in
// the right spot
void ExcitedResponse::augment_full(World &world, X_space &Chi, X_space &old_Chi, X_space &Lambda_X,
                                   X_space &last_Lambda_X, Tensor<double> &S, Tensor<double> &A,
                                   Tensor<double> &old_S, Tensor<double> &old_A,
                                   size_t print_level) {
    // Basic output
    if (print_level >= 1) molresponse::start_timer(world);

    size_t m = Chi.num_states();
    for (size_t i = 0; i < m; i++) {
        Chi.push_back(copy(world, old_Chi.X[i]), copy(world, old_Chi.Y[i]));
        Lambda_X.push_back(copy(world, last_Lambda_X.X[i]), copy(world, last_Lambda_X.Y[i]));
    }
    Tensor<double> temp_A = inner(Chi, Lambda_X);
    A = 0.5 * (temp_A + transpose(temp_A));
    S = response_space_inner(Chi.X, Chi.X) - response_space_inner(Chi.Y, Chi.Y);

    // End the timer
    if (print_level >= 1) molresponse::end_timer(world, "Aug. resp. matrix:");

    // Debugging output
    if (print_level >= 2 and world.rank() == 0) {
        print("\n   Augmented response matrix:");
        print(A);
    }

    // Debugging output
    if (print_level >= 2 and world.rank() == 0) {
        print("   Augmented overlap matrix:");
        print(S);
    }

    // SUPER debugging
    if (print_level >= 3) {
        if (world.rank() == 0) print("   Calculating condition number of aug. response matrix");
    }
}
// If using a larger subspace to diagonalize in, after diagonalization this
// will put everything in the right spot

void ExcitedResponse::unaugment(World &world, X_space &Chi, X_space &old_Chi, X_space &Lambda_X,
                                X_space &last_Lambda_X, Tensor<double> &omega, Tensor<double> &S_x,
                                Tensor<double> &A_x, Tensor<double> &old_S, Tensor<double> &old_A,
                                size_t num_states, size_t iter, size_t print_level) {
    // Basic output
    if (print_level >= 1) molresponse::start_timer(world);

    // Note: the eigenvalues and vectors were sorted after diagonalization
    // and hence all the functions are sorted in ascending order of energy

    // Quick copy of m lowest eigenvalues
    omega = omega(Slice(0, num_states - 1));
    // Pop off the "m" vectors off the back end of appropriate vectors
    // (only after first iteration)
    if (iter > 0) {
        for (size_t i = 0; i < num_states; i++) {
            Chi.X.pop_back();
            Lambda_X.X.pop_back();
        }
    }
    old_Chi.X = Chi.X.copy();
    last_Lambda_X.X = Lambda_X.X.copy();

    old_S = response_space_inner(Chi.X, Chi.X);
    old_A = Tensor<double>(num_states, num_states);
    for (size_t i = 0; i < num_states; i++) old_A(i, i) = omega(i);
    // End the timer
    if (print_level >= 1) molresponse::end_timer(world, "Unaug. resp. mat.:");
}
// If using a larger subspace to diagonalize in, after diagonalization this
// will put everything in the right spot

void ExcitedResponse::unaugment_full(World &world, X_space &Chi, X_space &old_Chi,
                                     X_space &Lambda_X, X_space &last_Lambda_X,
                                     Tensor<double> &omega, Tensor<double> &S_x,
                                     Tensor<double> &A_x, Tensor<double> &old_S,
                                     Tensor<double> &old_A, size_t num_states, size_t iter,
                                     size_t print_level) {
    // Basic output
    if (print_level >= 1) molresponse::start_timer(world);

    // Note: the eigenvalues and vectors were sorted after diagonalization
    // and hence all the functions are sorted in ascending order of energy

    // Quick copy of m lowest eigenvalues
    omega = omega(Slice(0, num_states - 1));

    // Pop off the "m" vectors off the back end of appropriate vectors
    // (only after first iteration)
    print("Entering Loop to pop_back Chi and LambdaX");
    if (iter > 0) {
        for (size_t i = 0; i < num_states; i++) {
            print("pop back Chi and LambdaX");
            Chi.pop_back();
            Lambda_X.pop_back();
        }
    }

    old_Chi = Chi.copy();
    last_Lambda_X = Lambda_X.copy();

    old_S = response_space_inner(Chi.X, Chi.X) - response_space_inner(Chi.Y, Chi.Y);
    old_A = Tensor<double>(num_states, num_states);
    for (size_t i = 0; i < num_states; i++) old_A(i, i) = omega(i);
    // End the timer
    if (print_level >= 1) molresponse::end_timer(world, "Unaug. resp. mat.:");
}
// Diagonalize the full response matrix, taking care of degenerate
// components Why diagonalization and then transform the x_fe vectors

Tensor<double> ExcitedResponse::diagonalizeFullResponseMatrix(
        World &world, X_space &Chi, X_space &Lambda_X, Tensor<double> &omega, Tensor<double> &S,
        Tensor<double> &A, const double thresh, size_t print_level) {
    // compute the unitary transformation matrix U that diagonalizes
    // the response matrix
    Tensor<double> U = GetFullResponseTransformation(world, S, A, omega, thresh);

    // Sort into ascending order
    // Tensor<int> selected = sort_eigenvalues(world, omega, U);

    // Start timer
    if (r_params.print_level() >= 1) molresponse::start_timer(world);

    Chi.X = transform(world, Chi.X, U);
    Chi.Y = transform(world, Chi.Y, U);
    Tensor<double> Sxa, Sya, Sa;

    Sxa = response_space_inner(Chi.X, Chi.X);
    Sya = response_space_inner(Chi.Y, Chi.Y);
    Sa = Sxa - Sya;

    if (world.rank() == 0 and r_params.print_level() >= 10) {
        print("\n  After apply transform Overlap Matrix:");
        print(Sxa);
        print(Sya);
        print(Sa);
    }

    Lambda_X.X = transform(world, Lambda_X.X, U);
    Lambda_X.Y = transform(world, Lambda_X.Y, U);
    // Transform the vectors of functions
    // Truncate happens in here
    // we do transform here
    // End timer
    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Transform orbs.:");

    // Normalize x and y
    normalize(world, Chi);

    // Debugging output
    if (world.rank() == 0 and print_level >= 2) {
        print("   Eigenvector coefficients from diagonalization:");
        print(U);
    }

    // Return the selected functions
    return U;
}

// Similar to what robert did above in "get_fock_transformation"
Tensor<double> ExcitedResponse::GetFullResponseTransformation(World &world, Tensor<double> &S,
                                                              Tensor<double> &A,
                                                              Tensor<double> &evals,
                                                              const double thresh_degenerate) {
    // Start timer
    if (r_params.print_level() >= 1) molresponse::start_timer(world);

    // Get size
    size_t m = S.dim(0);

    // Run an SVD on the overlap matrix and ignore values
    // less than thresh_degenerate
    Tensor<double> r_vecs, s_vals, l_vecs;
    Tensor<double> S_copy = copy(S);
    /**
   * @brief SVD on overlap matrix S
   * S=UsVT
   * S, U, s , VT
   *
   */
    svd(S_copy, l_vecs, s_vals, r_vecs);

    // Debugging output
    if (r_params.print_level() >= 2 and world.rank() == 0) {
        print("\n   Singular values of overlap matrix:");
        print(s_vals);
        print("   Left singular vectors of overlap matrix:");
        print(l_vecs);
    }

    // Check how many singular values are less than 10*thresh_degen
    size_t num_zero = 0;
    for (int64_t i = 0; i < s_vals.dim(0); i++) {
        if (s_vals(i) < 10 * thresh_degenerate) {
            if (world.rank() == 0 and num_zero == 0) print("");
            if (world.rank() == 0)
                printf("   Detected singular value (%.8f) below threshold (%.8f). "
                       "Reducing subspace size.\n",
                       s_vals(i), 10 * thresh_degenerate);
            num_zero++;
        }
        if (world.rank() == 0 and i == s_vals.dim(0) - 1 and num_zero > 0) print("");
    }

    // Going to use these a lot here, so just calculate them
    size_t size_l = s_vals.dim(0);    // number of singular values
    size_t size_s = size_l - num_zero;// smaller subspace size
    /**
   * @brief l_vecs_s(m,1)
   *
   * @return Tensor<double>
   */
    Tensor<double> l_vecs_s(size_l,
                            num_zero);// number of sv by number smaller than thress
    Tensor<double> copyA = copy(A);   // we copy xAx

    // Transform into this smaller space if necessary
    if (num_zero > 0) {
        print("num_zero = ", num_zero);
        // Cut out the singular values that are small
        // (singular values come out in descending order)

        // S(m-sl,m-sl)
        S = Tensor<double>(size_s, size_s);// create size of new size
        for (size_t i = 0; i < size_s; i++) S(i, i) = s_vals(i);
        // Copy the active vectors to a smaller container
        // left vectors [m,m-sl]
        l_vecs_s = copy(l_vecs(_, Slice(0, size_s - 1)));

        // Debugging output
        if (r_params.print_level() >= 2 and world.rank() == 0) {
            print("   Reduced size left singular vectors of overlap matrix:");
            print(l_vecs_s);
        }

        // Transform
        // Work(m,m-sl)
        Tensor<double> work(size_l, size_s);
        /*
  c(i,j) = c(i,j) + sum(k) a(i,k)*b(k,j)

  where it is assumed that the last index in each array is has unit
  stride and the dimensions are as provided.

  4-way unrolled k loop ... empirically fastest on PIII
  compared to 2/3 way unrolling (though not by much).
*/
        // dimi,dimj,dimk,c,a,b
        mxm(size_l, size_s, size_l, work.ptr(), A.ptr(), l_vecs_s.ptr());
        // A*left
        copyA = Tensor<double>(size_s, size_s);
        Tensor<double> l_vecs_t = transpose(l_vecs);
        // s s l, copyA=lvect_t*A*left
        mxm(size_s, size_s, size_l, copyA.ptr(), l_vecs_t.ptr(), work.ptr());

        // Debugging output
        if (r_params.print_level() >= 2 and world.rank() == 0) {
            print("   Reduced response matrix:");
            print(copyA);
            print("   Reduced overlap matrix:");
            print(S);
        }
    }
    // Diagonalize (NOT A SYMMETRIC DIAGONALIZATION!!!!)
    // Potentially complex eigenvalues come out of this
    Tensor<std::complex<double>> omega(size_s);
    Tensor<double> U(size_s, size_s);
    ggevp(world, copyA, S, U, omega);

    // Eigenvectors come out oddly packaged if there are
    // complex eigenvalues.
    // Currently only supporting real valued eigenvalues
    // so throw an error if any imaginary components are
    // not zero enough
    double max_imag = abs(imag(omega)).max();
    if (world.rank() == 0 and r_params.print_level() >= 2)
        print("\n   Max imaginary component of eigenvalues:", max_imag, "\n");
    MADNESS_ASSERT(max_imag <= r_params.dconv());// MUST BE REAL!
    evals = real(omega);

    // Easier to just resize here
    m = evals.dim(0);

    bool switched = true;
    while (switched) {
        switched = false;
        for (size_t i = 0; i < m; i++) {
            for (size_t j = i + 1; j < m; j++) {
                double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
                double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
                if (snew > sold) {
                    Tensor<double> tmp = copy(U(_, i));
                    U(_, i) = U(_, j);
                    U(_, j) = tmp;
                    std::swap(evals[i], evals[j]);
                    switched = true;
                }
            }
        }
    }

    // Fix phases.
    for (size_t i = 0; i < m; ++i)
        if (U(i, i) < 0.0) U(_, i).scale(-1.0);

    // Rotations between effectively degenerate components confound
    // the non-linear equation solver ... undo these rotations
    size_t ilo = 0;// first element of cluster
    while (ilo < m - 1) {
        size_t ihi = ilo;
        while (fabs(evals[ilo] - evals[ihi + 1]) <
               thresh_degenerate * 10.0 * std::max(fabs(evals[ilo]), 1.0)) {
            ++ihi;
            if (ihi == m - 1) break;
        }
        int64_t nclus = ihi - ilo + 1;
        if (nclus > 1) {
            Tensor<double> q = copy(U(Slice(ilo, ihi), Slice(ilo, ihi)));

            // Polar Decomposition
            Tensor<double> VH(nclus, nclus);
            Tensor<double> W(nclus, nclus);
            Tensor<double> sigma(nclus);

            svd(q, W, sigma, VH);
            q = transpose(inner(W, VH));// Should be conj. tranpose if complex
            U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), q);
        }
        ilo = ihi + 1;
    }

    // If we transformed into the smaller subspace, time to transform back
    if (num_zero > 0) {
        // Temp. storage
        Tensor<double> temp_U(size_l, size_l);
        Tensor<double> U2(size_l, size_l);

        // Copy U back to larger size
        temp_U(Slice(0, size_s - 1), Slice(0, size_s - 1)) = copy(U);
        for (size_t i = size_s; i < size_l; i++) temp_U(i, i) = 1.0;

        // Transform U back
        mxm(size_l, size_l, size_l, U2.ptr(), l_vecs.ptr(), temp_U.ptr());
        U = copy(U2);
    }

    // Sort into ascending order
    Tensor<int> selected = sort_eigenvalues(world, evals, U);

    // End timer
    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Diag. resp. mat.");

    return U;
}
Tensor<double> ExcitedResponse::diagonalizeFockMatrix(World &world, X_space &Chi, X_space &Lambda_X,
                                                      Tensor<double> &evals, Tensor<double> &A,
                                                      Tensor<double> &S,
                                                      const double thresh) const {
    // compute the unitary transformation matrix U that diagonalizes
    // the fock matrix
    Tensor<double> U = get_fock_transformation(world, S, A, evals, thresh);

    // Sort into ascending order
    Tensor<int> selected = sort_eigenvalues(world, evals, U);

    // Debugging output
    if (r_params.print_level() >= 2 and world.rank() == 0) {
        print("   U:");
        print(U);
    }

    // Start timer
    if (r_params.print_level() >= 1) molresponse::start_timer(world);

    // transform the orbitals and the potential
    // Truncate happens inside here
    Chi.X = transform(world, Chi.X, U);
    Lambda_X.X = transform(world, Lambda_X.X, U);

    // End timer
    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Transform orbs.:");

    // Normalize x
    normalize(world, Chi.X);

    // Debugging output
    if (r_params.print_level() >= 2 and world.rank() == 0) {
        print("   Eigenvector coefficients from diagonalization:");
        print(U);
    }

    return U;
}
Tensor<double> ExcitedResponse::get_fock_transformation(World &world, Tensor<double> &overlap,
                                                        Tensor<double> &fock, Tensor<double> &evals,
                                                        const double thresh_degenerate) const {
    // Run an SVD on the overlap matrix and ignore values
    // less than thresh_degenerate
    Tensor<double> r_vecs;
    Tensor<double> s_vals;
    Tensor<double> l_vecs;
    Tensor<double> overlap_copy = copy(overlap);
    svd(overlap_copy, l_vecs, s_vals, r_vecs);

    // Debugging output
    if (r_params.print_level() >= 2 and world.rank() == 0) {
        print("\n   Singular values of overlap matrix:");
        print(s_vals);
        print("   Left singular vectors of overlap matrix:");
        print(l_vecs);
    }

    // Check how many singular values are less than 10*thresh_degen
    size_t num_sv = 0;
    for (int64_t i = 0; i < s_vals.dim(0); i++) {
        if (s_vals(i) < 10 * thresh_degenerate) {
            if (world.rank() == 0 and num_sv == 0) print("");
            if (world.rank() == 0)
                printf("   Detected singular value (%.8f) below threshold (%.8f). "
                       "Reducing subspace size.\n",
                       s_vals(i), 10 * thresh_degenerate);
            num_sv++;
        }
        if (world.rank() == 0 and i == s_vals.dim(0) - 1 and num_sv > 0) print("");
    }

    // Going to use these a lot here, so just calculate them
    size_t size_l = s_vals.dim(0);
    size_t size_s = size_l - num_sv;
    Tensor<double> l_vecs_s(size_l, num_sv);

    // Transform into this smaller space if necessary
    if (num_sv > 0) {
        // Cut out the singular values that are small
        // (singular values come out in descending order)
        overlap = Tensor<double>(size_s, size_s);
        for (size_t i = 0; i < size_s; i++) overlap(i, i) = s_vals(i);

        // Copy the active vectors to a smaller container
        l_vecs_s = copy(l_vecs(_, Slice(0, size_s - 1)));

        // Debugging output
        if (r_params.print_level() >= 2 and world.rank() == 0) {
            print("   Reduced size left singular vectors of overlap matrix:");
            print(l_vecs_s);
        }

        // Transform
        Tensor<double> work(size_l, size_s);
        mxm(size_l, size_s, size_l, work.ptr(), fock.ptr(), l_vecs_s.ptr());
        fock = Tensor<double>(size_s, size_s);
        Tensor<double> l_vecs_t = transpose(l_vecs);
        mxm(size_s, size_s, size_l, fock.ptr(), l_vecs_t.ptr(), work.ptr());
    }

    // Diagonalize using lapack
    Tensor<double> U;
    sygv(fock, overlap, 1, U, evals);

    int64_t nmo = fock.dim(0);// NOLINT

    bool switched = true;
    while (switched) {
        switched = false;
        for (int64_t i = 0; i < nmo; i++) {
            for (int64_t j = i + 1; j < nmo; j++) {
                double sold = U(i, i) * U(i, i) + U(j, j) * U(j, j);
                double snew = U(i, j) * U(i, j) + U(j, i) * U(j, i);
                if (snew > sold) {
                    Tensor<double> tmp = copy(U(_, i));
                    U(_, i) = U(_, j);
                    U(_, j) = tmp;
                    std::swap(evals[i], evals[j]);
                    switched = true;
                }
            }
        }
    }

    // Fix phases.
    for (int64_t i = 0; i < nmo; ++i)// NOLINT
        if (U(i, i) < 0.0) U(_, i).scale(-1.0);

    // Rotations between effectively degenerate components confound
    // the non-linear equation solver ... undo these rotations
    int64_t ilo = 0;// first element of cluster NOLINT
    while (ilo < nmo - 1) {
        int64_t ihi = ilo;// NOLINT
        while (fabs(evals[ilo] - evals[ihi + 1]) <
               thresh_degenerate * 100.0 * std::max(fabs(evals[ilo]), 1.0)) {
            ++ihi;
            if (ihi == nmo - 1) break;
        }
        int64_t nclus = ihi - ilo + 1;// NOLINT
        if (nclus > 1) {
            Tensor<double> q = copy(U(Slice(ilo, ihi), Slice(ilo, ihi)));

            // Polar Decomposition
            Tensor<double> VH(nclus, nclus);
            Tensor<double> W(nclus, nclus);
            Tensor<double> sigma(nclus);

            svd(q, W, sigma, VH);
            q = transpose(inner(W, VH));// Should be conj. tranpose if complex
            U(_, Slice(ilo, ihi)) = inner(U(_, Slice(ilo, ihi)), q);
        }
        ilo = ihi + 1;
    }

    fock = 0;
    for (unsigned int i = 0; i < nmo; ++i) fock(i, i) = evals(i);

    // If we transformed into the smaller subspace, time to transform back
    if (num_sv > 0) {
        // Temp. storage
        Tensor<double> temp_U(size_l, size_l);
        Tensor<double> U2(size_l, size_l);

        // Copy U back to larger size
        temp_U(Slice(0, size_s - 1), Slice(0, size_s - 1)) = copy(U);
        for (size_t i = size_s; i < size_l; i++) temp_U(i, i) = 1.0;

        // Transform back
        mxm(size_l, size_l, size_l, U2.ptr(), l_vecs.ptr(), temp_U.ptr());

        U = copy(U2);
    }

    return U;
}
Tensor<int> ExcitedResponse::sort_eigenvalues(World &world, Tensor<double> &vals,
                                              Tensor<double> &vecs) const {
    // Get relevant sizes
    size_t k = vals.size();

    // Tensor to hold selection order
    Tensor<int> selected(k);

    // Copy everything...
    std::vector<double> vals_copy;
    for (size_t i = 0; i < k; i++) vals_copy.push_back(vals[i]);
    Tensor<double> vals_copy2 = copy(vals);
    Tensor<double> vecs_copy = copy(vecs);

    // Now sort vals_copy
    std::sort(vals_copy.begin(), vals_copy.end());

    // Now sort the rest of the things, using the sorted energy list
    // to find the correct indices
    for (size_t i = 0; i < k; i++) {
        // Find matching index in sorted vals_copy
        size_t j = 0;
        while (fabs(vals_copy[i] - vals_copy2[j]) > 1e-8 && j < k) j++;

        // Add in to list which one we're taking
        selected(i) = j;

        // Put corresponding things in the correct place
        vals(i) = vals_copy[i];
        vecs(_, i) = vecs_copy(_, j);

        // Change the value of vals_copy2[j] to help deal with duplicates?
        vals_copy2[j] = 10000.0;
    }

    // Done
    return selected;
}
Tensor<double> ExcitedResponse::create_shift(World &world, const Tensor<double> &ground,
                                             const Tensor<double> &omega, std::string xy) const {
    // Start a timer
    if (r_params.print_level() >= 1) molresponse::start_timer(world);

    // Get sizes
    size_t m = omega.size();
    size_t n = ground.size();

    // Container to hold shift
    Tensor<double> result(m, n);

    // Run over excited components
    for (size_t k = 0; k < m; k++) {
        // Run over ground components
        for (size_t p = 0; p < n; p++) {
            if (ground(p) + omega(k) > 0) {
                // Calculate the shift needed to get energy to -0.05,
                // which was arbitrary (same as moldft)
                result(k, p) = -(ground(p) + omega(k) + 0.05);

                // Basic output
                if (r_params.print_level() >= 2) {
                    if (world.rank() == 0)
                        printf("   Shift needed for transition from ground orbital %d to "
                               "response %s state %d\n",
                               static_cast<int>(p), xy.c_str(), static_cast<int>(k));
                    if (world.rank() == 0) print("   Ground energy =", ground(p));
                    if (world.rank() == 0) print("   Excited energy =", omega(k));
                    if (world.rank() == 0) print("   Shifting by", result(k, p));
                    if (world.rank() == 0) print("");
                }
            }
        }
    }

    // End timer
    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Create shift:");

    // Done
    return result;
}
response_space ExcitedResponse::apply_shift(World &world, const Tensor<double> &shifts,
                                            const response_space &V, const response_space &f) {
    // Start timer
    if (r_params.print_level() >= 1) molresponse::start_timer(world);

    // Sizes inferred from V
    size_t n = V[0].size();
    size_t m = V.size();

    // Container to return
    response_space shifted_V(world, m, n);

    // Run over occupied
    for (size_t k = 0; k < m; k++) {
        // Run over virtual
        for (size_t p = 0; p < n; p++) { shifted_V[k][p] = V[k][p] + shifts(k, p) * f[k][p]; }
    }

    shifted_V.truncate_rf();

    // End timer
    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Apply shift:");

    // Done
    return shifted_V;
}
std::vector<std::vector<std::shared_ptr<real_convolution_3d>>>
ExcitedResponse::create_bsh_operators(World &world, const Tensor<double> &shift,
                                      const Tensor<double> &ground, const Tensor<double> &omega,
                                      const double lo, const double thresh) const {
    // Start timer
    if (r_params.print_level() >= 1) molresponse::start_timer(world);

    // Sizes inferred from ground and omega
    size_t n = ground.size();
    size_t m = omega.size();

    // Make the vector
    std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> operators;

    // Make a BSH operator for each response function
    // Run over excited components
    for (size_t k = 0; k < m; k++) {
        // Container for intermediary
        std::vector<std::shared_ptr<real_convolution_3d>> temp(n);

        // Run over occupied components
        for (size_t p = 0; p < n; p++) {
            double mu = sqrt(-2.0 * (ground(p) + omega(k) + shift(k, p)));
            print("res state ", k, " orb ", p, " bsh exponent mu :", mu);
            temp[p] = std::shared_ptr<SeparatedConvolution<double, 3>>(
                    BSHOperatorPtr3D(world, mu, lo, thresh));
        }

        // Add intermediary to return container
        operators.push_back(temp);
    }

    // End timer
    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Creating BSH ops:");

    // Done
    return operators;
}
