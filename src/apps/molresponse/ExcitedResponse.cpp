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
        trial = create_trial_functions2(world);
        // Use a symmetry adapted operator on ground state functions
    } else {
        auto temp_trial = create_virtual_ao_guess(world);
        if (temp_trial.num_states() >= 2 * r_params.num_states()) {
            std::copy(temp_trial.x.begin(), temp_trial.x.begin() + 2 * r_params.num_states(),
                      trial.x.begin());

        } else if (temp_trial.num_states() >= r_params.num_states()) {
            trial = X_space(world, temp_trial.num_states(), r_params.num_orbitals());
            std::copy(temp_trial.x.begin(), temp_trial.x.begin() + temp_trial.num_states(),
                      trial.x.begin());

        } else {
            MADNESS_EXCEPTION("guess virtual ao did not produce enough states for calculation", 1);
        }
    }


    if (world.size() > 1) {
        // Start a timer
        if (r_params.num_orbitals() >= 1) molresponse::start_timer(world);
        if (world.rank() == 0) print("");// Makes it more legible

        LoadBalanceDeux<3> lb(world);
        for (size_t j = 0; j < r_params.num_states(); j++) {
            for (size_t k = 0; k < r_params.num_orbitals(); k++) {
                lb.add_tree(trial.x[j][k], lbcost<double, 3>(1.0, 8.0), true);
            }
        }
        for (size_t j = 0; j < r_params.num_orbitals(); j++) {
            lb.add_tree(ground_orbitals[j], lbcost<double, 3>(1.0, 8.0), true);
        }
        FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

        if (r_params.num_orbitals() >= 1) molresponse::end_timer(world, "Load balancing:");
    }

    // Project out ground state from guesses
    QProjector<double, 3> projector(ground_orbitals);
    for (unsigned int i = 0; i < trial.x.size(); i++) trial.x[i] = projector(trial.x[i]);

    // Ensure orthogonal guesses
    for (size_t i = 0; i < 2; i++) {
        molresponse::start_timer(world);
        // Orthog
        trial.x = gram_schmidt(world, trial.x);
        molresponse::end_timer(world, "orthog");

        molresponse::start_timer(world);
        // Normalize
        normalize(world, trial.x);
        molresponse::end_timer(world, "normalize");
    }

    // Diagonalize guess
    if (world.rank() == 0)
        print("\n   Iterating trial functions for an improved initial "
              "guess.\n");
    iterate_trial(world, trial);
    // Sort
    sort(world, omega, trial.x);
    // Basic output
    if (r_params.num_orbitals() >= 1 and world.rank() == 0) {
        print("\n   Final initial guess excitation energies:");
        print(omega);
    }
    // Chi = X_space(world, r_params.num_states(), r_params.num_orbitals());
    // Select lowest energy functions from guess
    Chi.x = select_functions(world, trial.x, omega, r_params.num_states(), r_params.num_orbitals());
    Chi.y = response_space(world, r_params.num_states(), r_params.num_orbitals());
    // save the guesses at the very least
    world.gop.fence();
    save(world, "guess_restart");

    trial.clear();
    // Initial guess for y are zero functions
}
/// an N-dimensional real-valued Gaussian function

/// the function looks like
/// \[
/// f(r) = x^i y^j .. z^k exp(-alpha r^2)
/// \]
template<std::size_t NDIM>
class GaussianGuess : public FunctionFunctorInterface<double, NDIM> {
    typedef Vector<double, NDIM> coordT;

public:
    /// ctor

    /// @param[in]  origin  the origin of the Gauss function
    /// @param[in]  alpha   the exponent exp(-alpha r^2)
    /// @param[in]  ijk     the monomial x^i y^j z^k exp(-alpha r^2) (for NDIM)
    GaussianGuess(const coordT &origin, const double alpha,
                  const std::vector<int> ijk = std::vector<int>(NDIM))
        : origin(origin), exponent(alpha), ijk(ijk) {}

    coordT origin;
    double exponent;     ///< exponent of the guess
    std::vector<int> ijk;///< cartesian exponents

    double operator()(const coordT &xyz) const {
        double arg = 0.0, prefac = 1.0;
        for (std::size_t i = 0; i < NDIM; ++i) {
            arg += (xyz[i] - origin[i]) * (xyz[i] - origin[i]);
            prefac *= pow(xyz[i], ijk[i]);
        }
        const double e = exponent * arg;
        return prefac * exp(-e);
    }
};

X_space ExcitedResponse::make_random_trial(World &world, size_t m) const {
    // Basic output
    if (world.rank() == 0) print("   Using a random guess for initial response functions.\n");
    size_t n = r_params.num_orbitals();
    // Create empty container and add in randomness
    X_space f(world, m, n);
    f.x = add_randomness(world, f.x, 1e3);// noise all over the world
    f.x = mask * f.x;                     // make sure you mask after you add noise to the world

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
    QProjector<double, 3> projector(ground_orbitals);
    for (unsigned int i = 0; i < f.num_states(); i++) f.x[i] = projector(f.x[i]);

    // Normalize
    normalize(world, f.x);

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
        for (unsigned int i = 0; i < rand.num_states(); i++) f.push_back(rand.x[i]);
    }

    // Project out groundstate from guesses
    QProjector<double, 3> projector(ground_orbitals);
    for (unsigned int i = 0; i < f.size(); i++) f[i] = projector(f[i]);

    // Truncate and normalize
    f.truncate_rf();
    normalize(world, f);

    X_space trial(world, f.size(), n);
    trial.x = f;

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

    // for every


    // Multiply each solid harmonic onto a ground state orbital
    for (size_t i = 0; i < n; i++) {
        // For each solid harmonic
        for (const auto &key: solids) {
            // Temp zero functions
            // // n temp zero functions
            std::vector<real_function_3d> temp =
                    zero_functions_compressed<double, 3>(world, static_cast<int>(n));

            // Create one non-zero function and add to trials
            print("index ", n - count % n - 1);
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
    trials.x = trials_X.copy();
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
                trials.x[count][o] = copy(functions.at(d).at(o));
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
        print_norms(world, trials.x);
    }

    // Truncate
    madness::truncate(world, trials.x);

    // Done
    return trials;
}

// Simplified iterate scheme for guesses
void ExcitedResponse::iterate_trial(World &world, X_space &guesses) {
    // Variables needed to iterate
    size_t iteration = 0;// Iteration counter
    QProjector<double, 3> projector(
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
        size_t N0 = guesses.x.size();
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
                    lb.add_tree(guesses.x[k][j], lbcost<double, 3>(1.0, 8.0), true);
                }
            }
            FunctionDefaults<3>::redistribute(world, lb.load_balance(2));

            if (r_params.print_level() >= 1) molresponse::end_timer(world, "Load balancing:");
        }

        // compute rho_omega
        auto rho_omega = transition_densityTDA(world, ground_orbitals, guesses.x);
        // Project out ground state
        for (size_t i = 0; i < Ni; i++) guesses.x[i] = projector(guesses.x[i]);

        // Truncate before doing expensive things
        guesses.x.truncate_rf();

        // Normalize after projection
        if (r_params.tda()) normalize(world, guesses.x);

        // (TODO why not normalize if not tda)
        // compute y = false
        auto xc = make_xc_operator(world);
        auto [temp_Lambda_X, temp_V0X, temp_gamma] =
                compute_response_potentials(world, guesses, xc, "tda");

        // Debugging output
        auto [new_omega, rotated_chi, rotated_lambda, rotated_v_x, rotated_gamma_x] =
                rotate_excited_space(world, guesses, temp_Lambda_X, temp_V0X, temp_gamma);

        omega = copy(new_omega);

        // Ensure right number of omegas
        if (size_t(omega.dim(0)) != Ni) {
            auto Ni_d = Ni - omega.dim(0);
            if (world.rank() == 0) {
                print("\n   Adding", Ni_d,
                      "eigenvalue(s) (counters subspace size "
                      "reduction in "
                      "diagonalizatoin).");
            }
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


            if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
            X_space E0X(world, rotated_chi.num_states(), rotated_chi.num_orbitals());
            if (r_params.localize() != "canon") {
                E0X = rotated_chi.copy();
                E0X.x = E0X.x * ham_no_diag;
            }
            world.gop.fence();

            if (r_params.print_level() >= 1) {
                molresponse::end_timer(world, "E0mDX", "E0mDX", iter_timing);
            }
            X_space theta_X = X_space(world, rotated_chi.num_states(), rotated_chi.num_orbitals());

            theta_X = rotated_v_x - E0X + rotated_gamma_x;


            theta_X.x = apply_shift(world, x_shifts, theta_X.x, guesses.x);
            theta_X.x = theta_X.x * -2;
            theta_X.x.truncate_rf();

            print("BARRIER before create BSH");
            world.mpi.Barrier();
            // Construct BSH operators
            auto bsh_x_operators =
                    create_bsh_operators(world, x_shifts, ground_energies, omega, r_params.lo(),
                                         FunctionDefaults<3>::get_thresh());

            print("BARRIER before apply BSH");
            world.mpi.Barrier();
            // Apply BSH and get updated components
            if (r_params.print_level() >= 1) molresponse::start_timer(world);
            bsh_resp = apply(world, bsh_x_operators, theta_X.x);
            if (r_params.print_level() >= 1) molresponse::end_timer(world, "Apply BSH:");

            // Project out ground state
            //for (size_t i = 0; i < Ni; i++) bsh_resp[i] = projector(bsh_resp[i]);

            for (auto &bsh_i: bsh_resp.x) { bsh_i = projector(bsh_i); }
            // Save new components
            guesses.x = bsh_resp;
            // Apply mask
            for (size_t i = 0; i < Ni; i++) guesses.x[i] = mask * guesses.x[i];
        }

        // Ensure orthogonal guesses
        for (size_t i = 0; i < 2; i++) {
            molresponse::start_timer(world);
            // Orthog
            guesses.x = gram_schmidt(world, guesses.x);
            molresponse::end_timer(world, "orthog");

            molresponse::start_timer(world);
            // Normalize
            normalize(world, guesses.x);
            molresponse::end_timer(world, "normalize");
        }

        // Update counter
        iteration += 1;
        // Done with the iteration.. truncate
        guesses.x.truncate_rf();

        // Basic output
        if (r_params.print_level() >= 1) {//
            molresponse::end_timer(world, " This iteration:");
        }
    }
}// Done with iterate gues
// Simplified iterate scheme for guesses

/**
 * @brief Diagonalize AX=SX frequencies
 *
 * @param world
 * @param S
 * @param old_S
 * @param old_A
 * @param x_response
 * @param old_x_response
 * @param ElectronResponses
 * @param OldElectronResponses
 * @param frequencies
 * @param iteration
 * @param m
 */

void ExcitedResponse::deflateGuesses(World &world, X_space &Chi, X_space &Lambda_X,
                                     Tensor<double> &S, Tensor<double> &frequencies,
                                     size_t &iteration, size_t &m) const {
    // XX =Omega XAX
    S = response_space_inner(Chi.x, Chi.x);
    Tensor<double> XAX = response_space_inner(Chi.x, Lambda_X.x);

    // Debugging output
    if (r_params.print_level() >= 2 and world.rank() == 0) {
        print(" Guess  Overlap matrix:");
        print(S);
        print(" Guess  XAX matrix:");
        print(XAX);
    }
    // Just to be sure dimensions work out, clear frequencies
    frequencies.clear();
    diagonalizeFockMatrix(world, Chi, Lambda_X, frequencies, XAX, S,
                          FunctionDefaults<3>::get_thresh());
}

void ExcitedResponse::deflateTDA(World &world, X_space &Chi, X_space &old_Chi, X_space &Lambda_X,
                                 X_space &old_Lambda_X, Tensor<double> &S, Tensor<double> old_S,
                                 Tensor<double> old_A, Tensor<double> &omega, size_t &iteration,
                                 size_t &m) {
    S = response_space_inner(Chi.x, Chi.x);
    Tensor<double> XAX = response_space_inner(Chi.x, Lambda_X.x);

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
        S = response_space_inner(Chi.x, Chi.x) - response_space_inner(Chi.y, Chi.y);
        if (world.rank() == 0 && (r_params.print_level() >= 10)) {
            print("\n   Overlap Matrix:");
            print(S);
            X_space Chi_copy = Chi.copy();
            Chi_copy.truncate();
            Lambda_X.truncate();
            A = inner(Chi_copy, Lambda_X);
            A = 0.5 * (A + transpose(A));
        }


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

std::tuple<Tensor<double>, X_space, X_space, X_space, X_space>
ExcitedResponse::rotate_excited_space(World &world, X_space &chi, X_space &lchi, X_space &v_chi,
                                      X_space &gamma_chi) {
    // Debugging output
    Tensor<double> A;
    X_space chi_copy = chi.copy();
    X_space l_copy = lchi.copy();

    Tensor<double> S = response_space_inner(chi_copy.x, chi_copy.x) -
                       response_space_inner(chi_copy.y, chi_copy.y);

    if (world.rank() == 0) {
        auto sm = S - transpose(S);
        print(sm.max());
    }
    S = 0.5 * (S + transpose(S));


    if (world.rank() == 0 && (r_params.print_level() >= 10)) {
        print("\n   Overlap Matrix:");
        print(S);
    }
    //
    A = inner(chi_copy, l_copy);
    if (world.rank() == 0) {
        auto am = A - transpose(A);
        print("largest non-symmetric :", am.max());
    }
    A = 0.5 * (A + transpose(A));
    if (world.rank() == 0 && (r_params.print_level() >= 10)) {
        print("\n   Lambda Matrix:");
        print(A);
    }

    auto [new_omega, U] = excited_eig(world, S, A, FunctionDefaults<3>::get_thresh());

    if (world.rank() == 0 and r_params.print_level() >= 2) {
        print("   Eigenvector coefficients from diagonalization:");
        print(U);
        print(new_omega);
    }
    auto [rotated_chi, rotated_l_chi, rotated_v_chi, rotated_gamma_chi] =
            rotate_excited_vectors(world, U, chi, lchi, v_chi, gamma_chi);

    return {new_omega, rotated_chi, rotated_l_chi, rotated_v_chi, rotated_gamma_chi};
}

std::tuple<Tensor<double>, Tensor<double>, Tensor<double>>
ExcitedResponse::reduce_subspace(World &world, Tensor<double> &S, Tensor<double> &A,
                                 const double thresh_degenerate) {


    // Get size
    size_t m = S.dim(0);
    // Run an SVD on the overlap matrix and ignore values
    // less than thresh_degenerate
    Tensor<double> r_vecs, s_vals, l_vecs;
    Tensor<double> S_copy = copy(S);
    // Step 1 find the svd of S
    svd(S_copy, l_vecs, s_vals, r_vecs);

    // Debugging output
    if (r_params.print_level() >= 2 and world.rank() == 0) {
        print("\n   Singular values of overlap matrix:");
        print(s_vals);
        print("   Left singular vectors of overlap matrix:");
        print(l_vecs);
    }

    // Step 2 find the number of singular values below threshold
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

    // in the overlap matrix
    size_t size_l = s_vals.dim(0);    // number of singular values
    size_t size_s = size_l - num_zero;// smaller subspace size

    Tensor<double> l_vecs_s(size_l, size_s);

    Tensor<double> copyA = copy(A);// we copy xAx
    // Transform into this smaller space if necessary
    if (num_zero > 0) {
        print("num_zero = ", num_zero);
        // Cut out the singular values that are small
        // (singular values come out in descending order)
        // S(m-sl,m-sl)
        S = Tensor<double>(size_s, size_s);// create size of new size

        // copy the singular values into the diagonal of smaller space
        for (size_t i = 0; i < size_s; i++) { S(i, i) = s_vals(i); }

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

    return {l_vecs, S, copyA};
}

std::pair<Tensor<double>, Tensor<double>>
ExcitedResponse::excited_eig(World &world, Tensor<double> &S, Tensor<double> &A,
                             const double thresh_degenerate) {
    // Start timer
    if (r_params.print_level() >= 1) molresponse::start_timer(world);
    auto size_l = S.dim(0);
    /*
    auto [l_vecs, copyS, copyA] = reduce_subspace(world, S, A, thresh_degenerate);
     */

    auto copyA = copy(A);
    auto copyS = copy(S);

    auto size_s = copyS.dim(0);
    auto num_zero = size_l - size_s;
    print("size_l: ", size_l);
    print("size_s: ", size_s);
    print("NUMZERO: ", num_zero);
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "reduce subspace", "subspace_reduce", iter_timing);
    }
    if (r_params.print_level() >= 1) molresponse::start_timer(world);
    // Diagonalize (NOT A SYMMETRIC DIAGONALIZATION!!!!)
    // Potentially complex eigenvalues come out of this
    Tensor<double> omega(size_s);
    Tensor<double> U(size_s, size_s);
    sygvp(world, copyA, copyS, 1, U, omega);
    //sygvp(world, fock, overlap, 1, c, e); from SCF.cc

    // not zero enough
    /*
    double max_imag = abs(imag(omega)).max();
    if (world.rank() == 0 and r_params.print_level() >= 2)
        print("\n   Max imaginary component of eigenvalues:", max_imag, "\n");
    if (max_imag > r_params.dconv()) {
        MADNESS_EXCEPTION("max imaginary component of eigenvalues > dconv", 0);
    }
     */
    Tensor<double> new_omega = real(omega);
    // Easier to just resize here
    auto m = new_omega.dim(0);
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
                    std::swap(new_omega[i], new_omega[j]);
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
        while (fabs(new_omega[ilo] - new_omega[ihi + 1]) <
               thresh_degenerate * 10.0 * std::max(fabs(new_omega[ilo]), 1.0)) {
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
    /*
    if (num_zero > 0) {
        // Temp. storage
        Tensor<double> temp_U(size_l, size_l);
        Tensor<double> temp_U2(size_l, size_l);
        Tensor<double> U2(size_l, size_l);

        // Copy U back to larger size
        temp_U(Slice(0, size_s - 1), Slice(0, size_s - 1)) = copy(U);
        for (size_t i = size_s; i < size_l; i++) temp_U(i, i) = 1.0;


        // Transform U back
        mxm(size_l, size_l, size_l, U2.ptr(), l_vecs.ptr(), temp_U.ptr());
        Tensor<double> l_vecs_t = transpose(l_vecs);
        mxm(size_l, size_l, size_l, temp_U2.ptr(), U2.ptr(), l_vecs_t.ptr());

        U = copy(temp_U2);
        if (world.rank() == 0 && r_params.print_level() >= 1) {
            print("Increasing subspace size and returning U");
            print(U);
        }
    }
     */

    // Sort into ascending order
    Tensor<int> selected = sort_eigenvalues(world, new_omega, U);

    // End timer
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "diagonalize response matrix", "diagonalize_response_matrix",
                               iter_timing);
    }

    return {new_omega, U};
}

void ExcitedResponse::augment(World &world, X_space &Chi, X_space &old_Chi, X_space &Lambda_X,
                              X_space &last_Lambda_X, Tensor<double> &S, Tensor<double> &A,
                              Tensor<double> &old_S, Tensor<double> &old_A, size_t print_level) {
    // Basic output
    if (print_level >= 1) molresponse::start_timer(world);

    // Get sizes
    size_t m = Chi.x.size();
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

    Tensor<double> off = response_space_inner(Chi.x, last_Lambda_X.x);
    temp_A(Slice(0, m - 1), Slice(m, 2 * m - 1)) = copy(off);// top right
    // Now for lower off diagonal
    off = response_space_inner(old_Chi.x, Lambda_X.x);
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
    off = expectation(world, Chi.x, old_Chi.x);
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
        Chi.x.push_back(old_Chi.x[i]);
        Lambda_X.x.push_back(last_Lambda_X.x[i]);
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
        Chi.push_back(copy(world, old_Chi.x[i]), copy(world, old_Chi.y[i]));
        Lambda_X.push_back(copy(world, last_Lambda_X.x[i]), copy(world, last_Lambda_X.y[i]));
    }
    Tensor<double> temp_A = inner(Chi, Lambda_X);
    A = 0.5 * (temp_A + transpose(temp_A));
    S = response_space_inner(Chi.x, Chi.x) - response_space_inner(Chi.y, Chi.y);

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
            Chi.x.pop_back();
            Lambda_X.x.pop_back();
        }
    }
    old_Chi.x = Chi.x.copy();
    last_Lambda_X.x = Lambda_X.x.copy();

    old_S = response_space_inner(Chi.x, Chi.x);
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

    old_S = response_space_inner(Chi.x, Chi.x) - response_space_inner(Chi.y, Chi.y);
    old_A = Tensor<double>(num_states, num_states);
    for (size_t i = 0; i < num_states; i++) old_A(i, i) = omega(i);
    // End the timer
    if (print_level >= 1) molresponse::end_timer(world, "Unaug. resp. mat.:");
}
// Diagonalize the full response matrix, taking care of degenerate
// components Why diagonalization and then transform the x_fe vectors


std::tuple<X_space, X_space, X_space, X_space>
ExcitedResponse::rotate_excited_vectors(World &world, const Tensor<double> &U, const X_space &chi,
                                        const X_space &l_chi, const X_space &v0_chi,
                                        const X_space &gamma_chi) {
    // compute the unitary transformation matrix U that diagonalizes
    // the response matrix

    // Start timer
    if (r_params.print_level() >= 1) molresponse::start_timer(world);


    auto rotated_chi = transform(world, chi, U);
    auto rotated_l_chi = transform(world, l_chi, U);
    auto rotated_v_chi = transform(world, v0_chi, U);
    auto rotated_gamma_chi = transform(world, gamma_chi, U);


    if (r_params.print_level() >= 10) {
        Tensor<double> S;
        S = response_space_inner(rotated_chi.x, rotated_chi.x) -
            response_space_inner(rotated_chi.y, rotated_chi.y);
        if (world.rank() == 0) {
            print("\n  After apply transform Overlap Matrix:");
            print(S);
        }
    }


    // Transform the vectors of functions
    // Truncate happens in here
    // we do transform here
    // End timer
    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Transform orbs.:");

    return {rotated_chi, rotated_l_chi, rotated_v_chi, rotated_gamma_chi};
}

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

    Chi.x = transform(world, Chi.x, U);
    Chi.y = transform(world, Chi.y, U);
    Tensor<double> Sxa, Sya, Sa;

    Sxa = response_space_inner(Chi.x, Chi.x);
    Sya = response_space_inner(Chi.y, Chi.y);
    Sa = Sxa - Sya;

    if (world.rank() == 0 and r_params.print_level() >= 10) {
        print("\n  After apply transform Overlap Matrix:");
        print(Sxa);
        print(Sya);
        print(Sa);
    }

    Lambda_X.x = transform(world, Lambda_X.x, U);
    Lambda_X.y = transform(world, Lambda_X.y, U);
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
    if (max_imag > r_params.dconv()) {
        MADNESS_EXCEPTION("max imaginary component of eigenvalues > dconv", 0);
    }
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
    Chi.x = transform(world, Chi.x, U);
    Lambda_X.x = transform(world, Lambda_X.x, U);

    // End timer
    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Transform orbs.:");

    // Normalize x
    normalize(world, Chi.x);

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
                if (r_params.print_level() >= 3) {
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

void ExcitedResponse::excited_to_json(json &j_mol_in, size_t iter, const Tensor<double> &omega) {
    json j = {};

    j["iter"] = iter;

    j["omega"] = tensor_to_json(omega);

    auto index = j_mol_in["protocol_data"].size() - 1;
    j_mol_in["protocol_data"][index]["property_data"].push_back(j);
}

void ExcitedResponse::iterate(World &world) {
    size_t iter;
    QProjector<double, 3> projector(ground_orbitals);
    size_t m = r_params.num_states();  // Number of excited states
    size_t n = r_params.num_orbitals();// Number of ground state orbitals

    const double conv_den = std::max(100 * FunctionDefaults<3>::get_thresh(), r_params.dconv());
    const double relative_max_target =
            std::max(50 * FunctionDefaults<3>::get_thresh(), .5 * r_params.dconv());

    auto thresh = FunctionDefaults<3>::get_thresh();
    auto max_rotation = .5;
    if (thresh >= 1e-2) {
        max_rotation = 2;
    } else if (thresh >= 1e-4) {
        max_rotation = .25;
    } else if (thresh >= 1e-6) {
        max_rotation = .1;
    } else if (thresh >= 1e-8) {
        max_rotation = .05;
    }

    // m residuals for x and y
    Tensor<double> bsh_residualsX(m);
    Tensor<double> density_residuals(m);
    Tensor<double> bsh_residualsY(m);

    Tensor<double> xij_norms(m, 2 * n);
    Tensor<double> xij_res_norms(m, 2 * n);
    // saved response densities
    vecfuncT rho_omega_old(m);
    // initialize DFT XC functional operator
    XCOperator<double, 3> xc = make_xc_operator(world);

    // create X space residuals
    X_space residuals(world, m, n);
    X_space old_Chi(world, m, n);
    X_space old_Lambda_X(world, m, n);

    // vector of Xvectors
    response_matrix x_vectors;
    response_matrix x_residuals;
    x_vectors = to_response_matrix(Chi);
    x_residuals = to_response_matrix(residuals);
    // If DFT, initialize the XCOperator<double,3>

    response_solver kain_x_space;
    size_t nkain = m;// (r_params.omega() != 0.0) ? 2 * m : m;
    for (size_t b = 0; b < m; b++) {
        kain_x_space.push_back(
                XNonlinearSolver<vector_real_function_3d, double, response_matrix_allocator>(
                        response_matrix_allocator(world, n), true));
    }
    if (r_params.kain()) {
        for (auto &kain_space_b: kain_x_space) { kain_space_b.set_maxsub(r_params.maxsub()); }
    }

    response_space bsh_x_resp(world, m, n);// Holds wave function corrections
    response_space bsh_y_resp(world, m, n);// Holds wave function corrections

    Tensor<double> old_energy(m);      // Holds previous iteration's energy
    Tensor<double> energy_residuals(m);// Holds energy residuals
    // Holds the norms of y function residuals (for convergence)
    Tensor<double> x_norms(m);
    Tensor<double> y_norms(m);

    response_space x_differences(world, m, n);
    response_space y_differences(world, m, n);


    Tensor<double> x_shifts;// Holds the shifted energy values
    Tensor<double> y_shifts;// Holds the shifted energy values

    Tensor<double> S;     // Overlap matrix of response components for x states
    real_function_3d v_xc;// For TDDFT
    Tensor<double> old_A;
    Tensor<double> A;
    Tensor<double> old_S;

    vector_real_function_3d rho_omega = make_density(world, Chi);

    // Create the X space
    all_done = false;// Converged flag
    // Now to iterate
    for (iter = 0; iter < r_params.maxiter(); ++iter) {

        iter_timing.clear();
        // Start a timer for this iteration
        // Basic output
        if (r_params.print_level() >= 1) {
            molresponse::start_timer(world);
            if (world.rank() == 0)
                printf("\n   Iteration %d at time %.1fs\n", static_cast<int>(iter), wall_time());
            if (world.rank() == 0) print(" -------------------------------");
        }

        print("Excited State Frequencies ");
        print(omega);


        // Normalize after projection

        if (iter < 2 || (iter % 10) == 0) { load_balance_chi(world); }

        if (iter > 0) {
            // Only checking on X components even for full as y are so small
            if (density_residuals.max() > 2) { break; }

            if (density_residuals.max() > 2) { break; }
            double d_residual = density_residuals.max();
            // Test convergence and set to true
            auto chi_norms = Chi.norm2s();
            auto rho_norms = norm2s_T(world, rho_omega);
            auto relative_bsh = copy(bsh_residualsX);

            std::transform(bsh_residualsX.ptr(), bsh_residualsX.ptr() + bsh_residualsX.size(),
                           chi_norms.ptr(), relative_bsh.ptr(),
                           [](auto bsh, auto norm_chi) { return bsh / norm_chi; });

            auto max_bsh = bsh_residualsX.absmax();
            auto relative_max_bsh = relative_bsh.absmax();


            function_data_to_json(j_molresponse, iter, chi_norms, bsh_residualsX, rho_norms,
                                  density_residuals);

            excited_to_json(j_molresponse, iter, omega);

            if (r_params.print_level() >= 1) {

                if (world.rank() == 0) {
                    print("thresh: ", FunctionDefaults<3>::get_thresh());
                    print("k: ", FunctionDefaults<3>::get_k());
                    print("Chi Norms at start of iteration: ", iter);
                    print("xij norms\n: ", xij_norms);
                    print("xij residual norms\n: ", xij_res_norms);
                    print("Chi_X: ", chi_norms);
                    print("bsh_residuals : ", bsh_residualsX);
                    print("relative_bsh : ", relative_bsh);
                    print("r_params.dconv(): ", r_params.dconv());
                    print("max rotation: ", max_rotation);
                    print("d_residual_max : ", d_residual);
                    print("d_residual_max target : ", conv_den);
                    print("bsh_residual_max : ", max_bsh);
                    print("relative_bsh_residual_max : ", relative_max_bsh);
                    print("relative_bsh_residual_max target : ", relative_max_target);
                }
            }
            if ((d_residual < conv_den) and ((relative_max_bsh < relative_max_target) or
                                             r_params.get<bool>("conv_only_dens"))) {
                all_done = true;
            }


            if (all_done || iter == r_params.maxiter() - 1) {
                // if converged print converged
                if (world.rank() == 0 && all_done and (r_params.print_level() > 1)) {
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
                    //plotResponseOrbitals(world, iter, Chi.x, Chi.y, r_params, ground_calc);
                }
                auto rho0 = make_ground_density(world);
//                if (r_params.plot()) {
//                    do_vtk_plots(world, 200, r_params.L(), molecule, rho0, rho_omega,
//                                 ground_orbitals, Chi);
//                }
                break;
            }
        }


        // We first rotate chi by diagonalizing AX=omegaBX
        // This provides new omegas
        // We then apply bsh on the rotated vector and compute the residual
        // The residual is then used to update KAIN
        // Followed by step restriction
        // residual is computed as new_chi-old_chi where both have been previously rotated.
        auto [new_omega, old_chi, new_chi, new_res] =
                update_response(world, Chi, xc, projector, kain_x_space, x_vectors, x_residuals,
                                iter, max_rotation, Tensor<double>(), residuals);

        residuals = new_res.residual.copy();

        if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
        rho_omega_old = make_density(world, old_chi);
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "make_density_old", "make_density_old", iter_timing);
        }
        if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
        rho_omega = make_density(world, new_chi);
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "make_density_new", "make_density_new", iter_timing);
        }

        if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
        bsh_residualsX = copy(new_res.residual_norms);
        bsh_residualsY = copy(new_res.residual_norms);
        omega = copy(new_omega);
        Chi = new_chi.copy();
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "copy_response_data", "copy_response_data", iter_timing);
        }

        xij_res_norms = new_res.residual.component_norm2s();
        xij_norms = Chi.component_norm2s();

        density_residuals = norm2s_T(world, (rho_omega - rho_omega_old));
        /*
        for (size_t i = 0; i < Chi.num_states(); i++) {
            if (maxrotn[i] < r_params.maxrotn()) {
                maxrotn[i] = r_params.maxrotn();
                print("less than maxrotn....set to maxrotn");
            }
        }
         */


        if (world.rank() == 0 and (r_params.print_level() > 2)) {
            print("Density residuals");
            print("dres", density_residuals);
            print("BSH  residuals");
            print("xres", bsh_residualsX);
            print("yres", bsh_residualsY);
            print("maxrotn", max_rotation);
        }


        // Basic output
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "Iteration Timing", "iter_total", iter_timing);
        }
        time_data.add_data(iter_timing);
    }

    if (world.rank() == 0) print("\n");
    if (world.rank() == 0) print("   Finished Excited State Calculation ");
    if (world.rank() == 0) print("   ------------------------");
    if (world.rank() == 0) print("\n");

    // Did we converge?
    if (iter == r_params.maxiter() && not all_done) {
        if (world.rank() == 0) print("   Failed to converge. Reason:");
        if (world.rank() == 0) print("\n  ***  Ran out of iterations  ***\n");
        if (world.rank() == 0) print("    Running analysis on current values.\n");
    }

    /*
    if (!r_params.tda()) {
        sort(world, omega, Chi);
    } else {
        sort(world, omega, Chi.X);
    }
     */

    // Print final things
    if (world.rank() == 0) {
        print(" Final excitation energies:");
        print(omega);
        print(" Final energy residuals X:");
        print(bsh_residualsX);
        print(" Final energy residuals y:");
        print(bsh_residualsY);
        print(" Final density residuals:");
        print(density_residuals);
    }

    /*
    analysis(world, Chi);
    print("--------------------------------------------------------");
    for (size_t i = 0; i < m; i++) {
        std::string x_state = "x_" + std::to_string(i) + "_";
        analyze_vectors(world, Chi.X[i], x_state);
        print("--------------------------------------------------------");
    }
    if (not r_params.tda()) {
        for (size_t i = 0; i < m; i++) {
            std::string y_state = "y_" + std::to_string(i) + "_";
            analyze_vectors(world, Chi.y[i], y_state);
            print("--------------------------------------------------------");
        }
    }
     */
}

auto ExcitedResponse::update_response(World &world, X_space &Chi, XCOperator<double, 3> &xc,
                                      QProjector<double, 3> &projector,
                                      response_solver &kain_x_space, response_matrix &Xvector,
                                      response_matrix &Xresidual, size_t iter,
                                      const double &maxrotn, const Tensor<double> old_residuals,
                                      const X_space &xres_old)
        -> std::tuple<Tensor<double>, X_space, X_space, residuals> {
    size_t m = Chi.num_states();
    bool compute_y = not r_params.tda();

    Tensor<double> x_shifts(m);
    Tensor<double> y_shifts(m);
    print("Entering Compute Lambda");

    /*
    if (compute_y) {
        gram_schmidt(world, Chi.X, Chi.y);
        normalize(world, Chi);
    } else {
        gram_schmidt(world, Chi.X);
        normalize(world, Chi.X);
    }
     */
    //
    // X_space Lambda_X = compute_lambda_X(world, Chi, xc, r_params.calc_type());

    auto [temp_Lambda_X, temp_V0X, temp_gamma] =
            compute_response_potentials(world, Chi, xc, r_params.calc_type());

    auto [new_omega, rotated_chi, rotated_lambda, rotated_v_x, rotated_gamma_x] =
            rotate_excited_space(world, Chi, temp_Lambda_X, temp_V0X, temp_gamma);

    print("omega_n before transform");
    print(omega);
    print("omega_n after transform");
    print(new_omega);
    // Analysis gets messed up if BSH is last thing applied
    // so exit early if last iteration


    //X_space theta_X = compute_theta_X(world, rotated_chi, xc, r_params.calc_type());
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    X_space rotated_EOX(world, rotated_chi.num_states(), rotated_chi.num_orbitals());
    if (r_params.localize() != "canon") {
        rotated_EOX = rotated_chi.copy();
        rotated_EOX.x = rotated_EOX.x * ham_no_diag;
        if (compute_y) { rotated_EOX.y = rotated_EOX.y * ham_no_diag; }
        if (r_params.print_level() >= 10) {
            print("<X|(E0-diag(E0)|X>");
            print(inner(rotated_chi, rotated_EOX));
        }
    }
    world.gop.fence();
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "E0mDX", "E0mDX", iter_timing);
    }

    X_space theta_X = X_space(world, rotated_chi.num_states(), rotated_chi.num_orbitals());

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    theta_X = rotated_v_x - rotated_EOX + rotated_gamma_x;
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "compute_ThetaX_add", "compute_ThetaX_add", iter_timing);
    }
    print("BSH update iter = ", iter);
    X_space new_chi = bsh_update_excited(world, new_omega, theta_X, projector);
    //res = Chi - new_chi;
    auto [new_res, bsh] = update_residual(world, rotated_chi, new_chi, r_params.calc_type(),
                                          old_residuals, xres_old);
    // kain if iteration >0 or first run where there should not be a problem
    // computed new_chi and res
    if (r_params.kain() && (iter > 0) && true) {
        new_chi = kain_x_space_update(world, rotated_chi, new_res, kain_x_space);
    }
    if (false) { x_space_step_restriction(world, rotated_chi, new_chi, compute_y, maxrotn); }


    if (compute_y) normalize(world, new_chi);
    else
        normalize(world, new_chi.x);


    new_chi.x.truncate_rf();
    if (compute_y) new_chi.y.truncate_rf();

    return {new_omega, rotated_chi, new_chi, {new_res, bsh}};
}

auto ExcitedResponse::bsh_update_excited(World &world, const Tensor<double> &omega,
                                         X_space &theta_X, QProjector<double, 3> &projector)
        -> X_space {
    size_t m = theta_X.num_states();
    size_t n = theta_X.num_orbitals();
    bool compute_y = !r_params.tda();
    Tensor<double> x_shifts(m);
    Tensor<double> y_shifts(m);
    print("omega before shifts");
    Tensor<double> omega_plus = omega;
    print(omega);
    x_shifts = create_shift(world, ground_energies, omega_plus, "x");
    // Compute Theta X
    // Apply the shifts
    theta_X.x = apply_shift(world, x_shifts, theta_X.x, Chi.x);
    theta_X.x = theta_X.x * -2;
    theta_X.x.truncate_rf();

    if (compute_y) {
        //   theta_X.y = apply_shift(world, y_shifts, theta_X.y, Chi.y);
        theta_X.y = theta_X.y * -2;
        theta_X.y.truncate_rf();
    }
    // Construct BSH operators
    std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_x_ops =
            create_bsh_operators(world, x_shifts, ground_energies, omega_plus, r_params.lo(),
                                 FunctionDefaults<3>::get_thresh());

    std::vector<std::vector<std::shared_ptr<real_convolution_3d>>> bsh_y_ops;
    if (compute_y) {
        Tensor<double> omega_minus = -omega;
        bsh_y_ops = create_bsh_operators(world, y_shifts, ground_energies, omega_minus,
                                         r_params.lo(), FunctionDefaults<3>::get_thresh());
    }
    X_space bsh_X(world, m, n);
    // Apply BSH and get updated response components
    bsh_X.x = apply(world, bsh_x_ops, theta_X.x);
    if (compute_y) bsh_X.y = apply(world, bsh_y_ops, theta_X.y);

    // Project out ground state
    for (size_t i = 0; i < m; i++) bsh_X.x[i] = projector(bsh_X.x[i]);
    if (compute_y) {
        for (size_t i = 0; i < m; i++) bsh_X.y[i] = projector(bsh_X.y[i]);
    }

    // Only update non-converged components
    /*
    for (size_t i = 0; i < m; i++) {
        bsh_X.X[i] = bsh_X.X[i];
        bsh_X.X[i] = mask * bsh_X.X[i];
        if (compute_y) {
            bsh_X.y[i] = bsh_X.y[i];
            bsh_X.y[i] = mask * bsh_X.y[i];
        }
    }
     */

    if (compute_y) normalize(world, bsh_X);
    else { normalize(world, bsh_X.x); }
    // Ensure orthogonal rguesses

    //bsh_X.truncate();

    return bsh_X;
}

void ExcitedResponse::analysis(World &world, const X_space &chi) {
    // Sizes get used a lot here, so lets get a local copy
    size_t n = chi.x[0].size();
    size_t m = chi.x.size();

    // Per response function, want to print the contributions from each
    // ground state So print the norm of each function?
    Tensor<double> x_norms(m, n);
    Tensor<double> y_norms(m, n);

    // Calculate the inner products
    for (long i = 0; i < m; i++) {
        for (long j = 0; j < n; j++) {
            x_norms(i, j) = chi.x[i][j].norm2();

            if (not r_params.tda()) y_norms(i, j) = chi.y[i][j].norm2();
        }
    }

    // 'sort' these inner products within in each row
    Tensor<double> cpy = copy(x_norms);
    Tensor<int> x_order(m, n);
    Tensor<int> y_order(m, n);
    for (long i = 0; i < m; i++) {
        for (long j = 0; j < n; j++) {
            double x = cpy(i, _).max();
            size_t z = 0;
            while (x != cpy(i, z)) z++;
            cpy(i, z) = -100.0;
            x_order(i, j) = z;

            // Also sort y if full response
            if (not r_params.tda()) { y_order(i, j) = z; }
        }
    }

    // Need these to calculate dipole/quadrapole
    real_function_3d x = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{1, 0, 0})));
    real_function_3d y = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{0, 1, 0})));
    real_function_3d z = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{0, 0, 1})));

    // Calculate transition dipole moments for each response function
    Tensor<double> dipoles(m, 3);

    // Run over each excited state
    for (size_t i = 0; i < m; i++) {
        // Add in contribution from each ground state
        for (size_t j = 0; j < n; j++) {
            dipoles(i, 0) += inner(ground_orbitals[j], x * chi.x[i][j]);
            dipoles(i, 1) += inner(ground_orbitals[j], y * chi.x[i][j]);
            dipoles(i, 2) += inner(ground_orbitals[j], z * chi.x[i][j]);

            if (not r_params.tda()) {
                dipoles(i, 0) += inner(ground_orbitals[j], x * chi.y[i][j]);
                dipoles(i, 1) += inner(ground_orbitals[j], y * chi.y[i][j]);
                dipoles(i, 2) += inner(ground_orbitals[j], z * chi.y[i][j]);
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
        oscillator(i) = 2.0 / 3.0 *
                        (dipoles(i, 0) * dipoles(i, 0) + dipoles(i, 1) * dipoles(i, 1) +
                         dipoles(i, 2) * dipoles(i, 2)) *
                        omega(i);
    }

    // Calculate transition quadrapole moments
    Tensor<double> quadrupoles(m, 3, 3);

    // Run over each excited state
    for (long i = 0; i < m; i++) {
        // Add in contribution from each ground state
        for (long j = 0; j < n; j++) {
            quadrupoles(i, 0, 0) += inner(ground_orbitals[j], x * x * chi.x[i][j]);
            quadrupoles(i, 0, 1) += inner(ground_orbitals[j], x * y * chi.x[i][j]);
            quadrupoles(i, 0, 2) += inner(ground_orbitals[j], x * z * chi.x[i][j]);
            quadrupoles(i, 1, 0) += inner(ground_orbitals[j], y * x * chi.x[i][j]);
            quadrupoles(i, 1, 1) += inner(ground_orbitals[j], y * y * chi.x[i][j]);
            quadrupoles(i, 1, 2) += inner(ground_orbitals[j], y * z * chi.x[i][j]);
            quadrupoles(i, 2, 0) += inner(ground_orbitals[j], z * x * chi.x[i][j]);
            quadrupoles(i, 2, 1) += inner(ground_orbitals[j], z * y * chi.x[i][j]);
            quadrupoles(i, 2, 2) += inner(ground_orbitals[j], z * z * chi.x[i][j]);

            if (not r_params.tda()) {
                quadrupoles(i, 0, 0) += inner(ground_orbitals[j], x * x * chi.y[i][j]);
                quadrupoles(i, 0, 1) += inner(ground_orbitals[j], x * y * chi.y[i][j]);
                quadrupoles(i, 0, 2) += inner(ground_orbitals[j], x * z * chi.y[i][j]);
                quadrupoles(i, 1, 0) += inner(ground_orbitals[j], y * x * chi.y[i][j]);
                quadrupoles(i, 1, 1) += inner(ground_orbitals[j], y * y * chi.y[i][j]);
                quadrupoles(i, 1, 2) += inner(ground_orbitals[j], y * z * chi.y[i][j]);
                quadrupoles(i, 2, 0) += inner(ground_orbitals[j], z * x * chi.y[i][j]);
                quadrupoles(i, 2, 1) += inner(ground_orbitals[j], z * y * chi.y[i][j]);
                quadrupoles(i, 2, 2) += inner(ground_orbitals[j], z * z * chi.y[i][j]);
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
        for (long i = 0; i < m; i++) {
            printf("   Response Function %d\t\t%7.8f a.u.", static_cast<int>(i), omega(i));
            print("\n   --------------------------------------------");
            printf("   Response Function %d\t\t%7.8f eV", static_cast<int>(i), omega(i) * 27.2114);
            print("\n   --------------------------------------------");

            print("\n   Transition Dipole Moments");
            printf("   X: %7.8f   y: %7.8f   Z: %7.8f\n", dipoles(i, 0), dipoles(i, 1),
                   dipoles(i, 2));

            printf("\n   Dipole Oscillator Strength: %7.8f\n", oscillator(i));

            print("\n   Transition Quadrupole Moments");
            printf("   %16s %16s %16s\n", "X", "y", "Z");
            printf("   X %16.8f %16.8f %16.8f\n", quadrupoles(i, 0, 0), quadrupoles(i, 0, 1),
                   quadrupoles(i, 0, 2));
            printf("   y %16.8f %16.8f %16.8f\n", quadrupoles(i, 1, 0), quadrupoles(i, 1, 1),
                   quadrupoles(i, 1, 2));
            printf("   Z %16.8f %16.8f %16.8f\n", quadrupoles(i, 2, 0), quadrupoles(i, 2, 1),
                   quadrupoles(i, 2, 2));

            // Print contributions
            // Only print the top 5?
            if (r_params.tda()) {
                print("\n   Dominant Contributions:");
                for (long j = 0; j < std::min(size_t(5), n); j++) {
                    printf("   Occupied %d   %7.8f\n", x_order(i, j), x_norms(i, x_order(i, j)));
                }

                print("\n");
            } else {
                print("\n   Dominant Contributions:");
                print("                  x          y");
                for (long j = 0; j < std::min(size_t(5), n); j++) {
                    printf("   Occupied %d   %7.8f %7.8f\n", x_order(i, j),
                           x_norms(i, x_order(i, j)), y_norms(i, y_order(i, j)));
                }

                print("\n");
            }
        }
    }
}
// Save the current response calculation

void ExcitedResponse::save(World &world, const std::string &name) {


    // Archive to write everything to
    archive::ParallelOutputArchive ar(world, name.c_str(), 1);
    // Just going to enforce 1 io server

    // Saving, in this order;
    //  string           ground-state archive name (garch_name)
    //  bool             TDA flag
    // size_t                number of ground state orbitals (n)
    // size_t                number of excited state orbitals (m)
    //  Tensor<double>   energies of m x-components
    //  for i from 0 to m-1
    //     for j from 0 to n-1
    //        Function<double,3> x_response[i][j]
    //  (If TDA flag == True)
    //  (Tensor<double>  energies of m y-components    )
    //  (for i from 0 to m-1                       )
    //  (   for j from 0 to n-1                    )
    //  (      Function<double,3> y_response[i][j] )
    ar & r_params.archive();
    ar & r_params.tda();
    ar & r_params.num_orbitals();
    ar & r_params.num_states();
    ar & omega;

    for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++) ar & Chi.x[i][j];
    if (not r_params.tda()) {
        for (size_t i = 0; i < r_params.num_states(); i++)
            for (size_t j = 0; j < r_params.num_orbitals(); j++) ar & Chi.y[i][j];
    }
}

// Load a response calculation
void ExcitedResponse::load(World &world, const std::string &name) {
    // The archive to read from
    archive::ParallelInputArchive ar(world, name.c_str());

    // Reading in, in this order;
    //  string           ground-state archive name (garch_name)
    //  bool             TDA flag
    // size_t                number of ground state orbitals (n)
    // size_t                number of excited state orbitals (m)
    //  Tensor<double>   energies of m x-components
    //  for i from 0 to m-1
    //     for j from 0 to n-1
    //        Function<double,3> x_response[i][j]
    //  (If TDA flag == True)
    //  (Tensor<double>  energies of m y-components    )
    //  (for i from 0 to m-1                       )
    //  (   for j from 0 to n-1                    )
    //  (      Function<double,3> y_response[i][j] )

    ar & r_params.archive();
    ar & r_params.tda();
    ar & r_params.num_orbitals();
    ar & r_params.num_states();
    ar & omega;

    Chi = X_space(world, r_params.num_states(), r_params.num_orbitals());

    for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++) ar & Chi.x[i][j];
    world.gop.fence();

    if (not r_params.tda()) {
        for (size_t i = 0; i < r_params.num_states(); i++)
            for (size_t j = 0; j < r_params.num_orbitals(); j++) ar & Chi.y[i][j];
        world.gop.fence();
    }
}

/**
 * Create Virtual Space Excited State Guess from atomic orbitals
 * @param world
 * @return
 */
X_space ExcitedResponse::create_virtual_ao_guess(World &world) const {

    auto phi_0 = copy(world, ground_orbitals);
    auto ao_basis_set = AtomicBasisSet{"aug-cc-pvdz"};
    if (world.rank() == 0) { ao_basis_set.print(molecule); }

    vecfuncT ao_vec = vecfuncT(ao_basis_set.nbf(molecule));

    for (int i = 0; i < ao_basis_set.nbf(molecule); ++i) {
        functorT aofunc(new madchem::AtomicBasisFunctor(
                ao_basis_set.get_atomic_basis_function(molecule, i)));
        ao_vec[i] = factoryT(world).functor(aofunc);
    }
    world.gop.fence();

    QProjector<double, 3> projector(phi_0);

    // project ground_state from each atomic orbital basis
    std::transform(ao_vec.begin(), ao_vec.end(), ao_vec.begin(),
                   [&](auto &f) { return projector(f); });
    auto overlap_S = matrix_inner(world, ao_vec, ao_vec);
    Tensor<double> U, sigma, VT;
    //S=U*sigma*VT
    svd(overlap_S, U, sigma, VT);

    auto first_small = std::find_if(sigma.ptr(), sigma.ptr() + sigma.size(),
                                    [](auto num) { return num < .05; });

    auto idx_small = std::distance(sigma.ptr(), first_small);
    auto xao = transform(world, ao_vec, U);
    xao.erase(xao.begin() + idx_small, xao.end());
    // remove small singular vectors here

    auto overlap_SU = matrix_inner(world, xao, xao);

    // copy pao to xao response space vector

    auto xc = make_xc_operator(world);
    // now compute the fock potential
    // kinetic energy
    auto F = [&](vector_real_function_3d x) {
        size_t num_orbs = x.size();
        real_derivative_3d Dx(world, 0);
        real_derivative_3d Dy(world, 1);
        real_derivative_3d Dz(world, 2);
        vector_real_function_3d dvx = apply(world, Dx, x);
        vector_real_function_3d dvy = apply(world, Dy, x);
        vector_real_function_3d dvz = apply(world, Dz, x);
        // Apply again for 2nd derivatives
        vector_real_function_3d dvx2 = apply(world, Dx, dvx);
        vector_real_function_3d dvy2 = apply(world, Dy, dvy);
        vector_real_function_3d dvz2 = apply(world, Dz, dvz);

        auto T = (dvx2 + dvy2 + dvz2) * (-0.5);
        real_function_3d v_nuc, v_j0, v_k0, v_xc;
        v_nuc = copy(potential_manager->vnuclear());
        v_nuc.truncate();
        auto N_elec = phi_0.size();
        v_nuc.scale((N_elec + 0.25) / N_elec);
        // J^0 x^alpha
        v_j0 = apply(*shared_coulomb_operator, ground_density);
        v_j0.scale(2.0);

        if (xcf.hf_exchange_coefficient() != 1.0) {
            v_xc = xc.make_xc_potential();
        } else {
            // make a zero function
            v_xc = Function<double, 3>(
                    FunctionFactory<double, 3>(world).fence(false).initial_level(1));
        }

        auto hf_exchange_X = zero_functions<double, 3>(world, num_orbs);
        // hf exchange
        std::transform(x.begin(), x.end(), hf_exchange_X.begin(), [&](auto &xi) {
            auto f = Function<double, 3>(FunctionFactory<double, 3>(world));
            std::accumulate(phi_0.begin(), phi_0.end(), f, [&](auto phi0, auto total) {
                auto sum = apply(*shared_coulomb_operator, xi * phi0) * phi0;
                return total + sum;
            });
            return f;
        });
        vector_real_function_3d V0 = zero_functions<double, 3>(world, num_orbs);
        real_function_3d v0 = v_j0 + v_nuc + v_xc;
        V0 = v0 * x;
        V0 += (-1 * hf_exchange_X * xcf.hf_exchange_coefficient());

        return T + V0;
    };

    auto phi_a = copy(world, xao);
    Tensor<double> e_a;
    for (int i = 0; i < 1; i++) {
        auto Fx = F(phi_a);
        auto xFx = matrix_inner(world, phi_a, Fx);
        auto [e, C] = syev(xFx);
        print("eigs  : \n", e);
        phi_a = transform(world, phi_a, C);
        projector(phi_a);
        e_a = e;
    }

    auto e_homo = ground_energies[phi_0.size() - 1];

    if (world.rank() == 0) {
        print("ground_orb energies  : \n", ground_energies);
        print("homo energy  : \n", e_homo);
        print("virtual orbital energies  : \n", e_a);
    }

    e_a = e_a.flat();
    auto is_positive = [&](auto num) { return num > e_homo; };
    // This seems dumb but i'm removing the negative vectors
    auto first_positive = std::find_if(e_a.ptr(), e_a.ptr() + e_a.size(), is_positive);
    auto num_negative = std::distance(e_a.ptr(), first_positive);
    auto virtual_phi = vector_real_function_3d(phi_a.size() - num_negative);
    Tensor<double> virtual_e(e_a.size() - num_negative);

    std::copy(first_positive, e_a.ptr() + e_a.size(), virtual_e.ptr());
    std::copy(phi_a.begin() + num_negative, phi_a.end(), virtual_phi.begin());

    phi_a = virtual_phi;
    e_a = virtual_e;


    auto t = phi_a.size() * r_params.num_orbitals();
    auto no = r_params.num_orbitals();

    X_space x_guess(world, t, no);

    int k = 0;
    // for each orbital
    // for each ground orbital j
    // add the virtual orbital in location j
    // therefore there should be a total of num_virt*num_ground_orbitals
    std::for_each(phi_a.begin(), phi_a.end(), [&](const auto virt) {
        for (int j = 0; j < no; j++) { x_guess.x[k++][j] = copy(virt); }
    });
    world.gop.fence();
    return x_guess;
}

/**
 * Create Virtual Space Excited State Guess from atomic orbitals
 * @param world
 * @return
 */
X_space ExcitedResponse::create_response_guess(World &world) const {

    print("thresh : ", FunctionDefaults<3>::get_thresh());
    print("k : ", FunctionDefaults<3>::get_k());
    auto phi_0 = copy(world, ground_orbitals);
    print("Ground Orbital norms: ", norm2s_T(world, phi_0));
    auto ao_basis_set = AtomicBasisSet{"6-31g"};
    ao_basis_set.print(molecule);

    vecfuncT ao_vec = vecfuncT(ao_basis_set.nbf(molecule));

    for (int i = 0; i < ao_basis_set.nbf(molecule); ++i) {
        functorT aofunc(new madchem::AtomicBasisFunctor(
                ao_basis_set.get_atomic_basis_function(molecule, i)));
        ao_vec[i] = factoryT(world).functor(aofunc);
    }
    world.gop.fence();
    print("number of ao_basis_functions: ", ao_vec.size());

    madness::print("norm of ao basis functions: ", norm2s_T(world, ao_vec));

    QProjector<double, 3> projector(phi_0);

    // project ground_state from each atomic orbital basis
    std::transform(ao_vec.begin(), ao_vec.end(), ao_vec.begin(),
                   [&](auto &f) { return projector(f); });
    print("overlap between phi0: ", matrix_inner(world, phi_0, ao_vec));
    madness::print("norm of ao basis after projection: ", norm2s_T(world, ao_vec));
    print("number of ao_basis_functions: ", ao_vec.size());
    auto overlap_S = matrix_inner(world, ao_vec, ao_vec);
    print("Overlap S : \n", overlap_S);
    Tensor<double> U, sigma, VT;
    //S=U*sigma*VT
    svd(overlap_S, U, sigma, VT);
    print("singular values of overlap: \n", sigma);
    print("left singular vectors of overlap:\n", U);
    print("right singular vectors of overlap:\n", VT);
    auto xao = transform(world, ao_vec, U);
    auto overlap_SU = matrix_inner(world, xao, xao);
    print("Overlap S after transform : \n", overlap_SU);

    // copy pao to xao response space vector

    auto xc = make_xc_operator(world);
    // now compute the fock potential
    // kinetic energy
    auto F = [&](vector_real_function_3d x) {
        size_t num_orbs = x.size();
        real_derivative_3d Dx(world, 0);
        real_derivative_3d Dy(world, 1);
        real_derivative_3d Dz(world, 2);
        vector_real_function_3d dvx = apply(world, Dx, x);
        vector_real_function_3d dvy = apply(world, Dy, x);
        vector_real_function_3d dvz = apply(world, Dz, x);
        // Apply again for 2nd derivatives
        vector_real_function_3d dvx2 = apply(world, Dx, dvx);
        vector_real_function_3d dvy2 = apply(world, Dy, dvy);
        vector_real_function_3d dvz2 = apply(world, Dz, dvz);

        auto T = (dvx2 + dvy2 + dvz2) * (-0.5);
        real_function_3d v_nuc, v_j0, v_k0, v_xc;
        v_nuc = potential_manager->vnuclear();
        v_nuc.truncate();
        // J^0 x^alpha
        v_j0 = apply(*shared_coulomb_operator, ground_density);
        v_j0.scale(2.0);

        if (xcf.hf_exchange_coefficient() != 1.0) {
            v_xc = xc.make_xc_potential();
        } else {
            // make a zero function
            v_xc = Function<double, 3>(
                    FunctionFactory<double, 3>(world).fence(false).initial_level(1));
        }

        auto hf_exchange_X = zero_functions<double, 3>(world, num_orbs);
        // hf exchange
        std::transform(x.begin(), x.end(), hf_exchange_X.begin(), [&](auto &xi) {
            auto f = Function<double, 3>(FunctionFactory<double, 3>(world));
            std::accumulate(phi_0.begin(), phi_0.end(), f, [&](auto phi0, auto total) {
                auto sum = apply(*shared_coulomb_operator, xi * phi0) * phi0;
                return total + sum;
            });
            return f;
        });
        vector_real_function_3d V0 = zero_functions<double, 3>(world, num_orbs);
        real_function_3d v0 = v_j0 + v_nuc + v_xc;
        V0 = v0 * x;
        V0 += (-1 * hf_exchange_X * xcf.hf_exchange_coefficient());

        return T + V0;
    };

    auto phi_a = copy(world, xao);
    Tensor<double> e_a;
    for (int i = 0; i < 1; i++) {
        auto Fx = F(phi_a);
        auto xFx = matrix_inner(world, phi_a, Fx);
        auto [e, C] = syev(xFx);
        print("eigs  : \n", e);
        phi_a = transform(world, phi_a, C);
        projector(phi_a);
        e_a = e;
    }

    print("ground_orb energies  : \n", ground_energies);
    print("virtual orbital energies  : \n", e_a);
    print("norm of virtual ", norm2s_T(world, phi_a));
    e_a = e_a.flat();
    auto homo_e = ground_energies[phi_0.size() - 1];
    print("homo_e  :", homo_e);
    auto is_positive = [&](auto num) { return num > homo_e; };
    // This seems dumb but i'm removing the negative vectors
    auto first_positive = std::find_if(e_a.ptr(), e_a.ptr() + e_a.size(), is_positive);
    auto num_negative = std::distance(e_a.ptr(), first_positive);
    print("num_negative :", num_negative);
    auto virtual_phi = vector_real_function_3d(phi_a.size() - num_negative);
    Tensor<double> virtual_e(e_a.size() - num_negative);
    std::copy(first_positive, e_a.ptr() + e_a.size(), virtual_e.ptr());
    std::copy(phi_a.begin() + num_negative, phi_a.end(), virtual_phi.begin());
    phi_a = virtual_phi;
    e_a = virtual_e;
    auto SOA = matrix_inner(world, phi_0, phi_a);
    print("Do I remove the vectors that are not completely 0 I hope it's zero\n", SOA);
    print("Sanity CHECK I hope it's diagonal\n", matrix_inner(world, phi_a, F(phi_a)));
    // create a vector of pairs
    std::vector<std::pair<int, int>> ia_indicies;
    print("forming response space pairs");
    print("( i , a )");
    for (int i = 0; i < phi_0.size(); i++) {
        for (int a = 0; a < phi_a.size(); a++) {
            ia_indicies.emplace_back(std::pair<int, int>{i, a});
        }
    }
    int kk = 0;
    Tensor<double> A(ia_indicies.size(), ia_indicies.size());
    for (const auto ia: ia_indicies) {
        std::cout << "( " << ia.first << " , " << ia.second << " )" << std::endl;
        A(kk, kk) = (e_a[ia.second] - ground_energies[ia.first]);
        kk++;
    }
    print("( ij , ab )");
    int ii = 0;
    int jj = 0;
    // create pairs ia, ii,ab
    auto two_int = [&](World &world, const vector_real_function_3d &phi0,
                       const vector_real_function_3d &phia) {
        // The easy case make n*m matrix
        auto tol = FunctionDefaults<3>::get_thresh();
        reconstruct(world, phi0);
        reconstruct(world, phia);

        vector_real_function_3d pairs_ij;
        vector_real_function_3d pairs_ia;
        vector_real_function_3d pairs_ab;

        for (auto &phi_ii: phi0) {
            for (auto &phi_aa: phia) { pairs_ia.push_back(mul_sparse(phi_ii, phi_aa, tol, false)); }
        }
        world.gop.fence();
        truncate(world, pairs_ia);
        vecfuncT Vpairs_ia = apply(world, *shared_coulomb_operator, pairs_ia);
        auto A_ia_jb = matrix_inner(world, pairs_ia, Vpairs_ia);
        for (auto &phi_ii: phi0) {
            for (auto &phi_jj: phi0) { pairs_ij.push_back(mul_sparse(phi_ii, phi_ii, tol, false)); }
        }
        for (auto &phi_aa: phia) {
            for (auto &phi_bb: phia) { pairs_ab.push_back(mul_sparse(phi_aa, phi_aa, tol, false)); }
        }
        world.gop.fence();
        truncate(world, pairs_ij);
        truncate(world, pairs_ab);
        vecfuncT Vpairs_ij = apply(world, *shared_coulomb_operator, pairs_ab);
        auto A_ij_ab = matrix_inner(world, pairs_ia, Vpairs_ij);
        // reshape A_ij_ab  n^2 x m^2
        auto a_len = phi_a.size() * phi_0.size();
        Tensor<double> A(a_len, a_len);
        print(A);
        int kk = 0;
        //[i,j]= will have a phia by phia block.  each row gets rearranged into ij
        for (int ii = 0; ii < phi0.size(); ii++) {
            for (int jj = 0; jj < phi0.size(); jj++) {
                auto ii_ab = A_ij_ab(kk++, _);
                print(ii_ab);
                auto b_i = ii * phia.size();
                auto b_ip1 = (ii * phia.size() + phia.size()) - 1;
                auto b_j = jj * phia.size();
                auto b_jp1 = (jj * phia.size() + phia.size()) - 1;
                print("( ", b_i, " , ", b_ip1, " ) ");
                print("( ", b_j, " , ", b_jp1, " ) ");
                print(Slice(b_i, b_ip1));
                //ii_ab = ii_ab.reshape(phi_a.size(), phi_a.size());
                A(Slice(b_i, b_ip1, 1), Slice(b_j, b_jp1, 1)) =
                        ii_ab.reshape(phi_a.size(), phi_a.size());
                //ii_ab = ii_ab.reshape(phi_a.size(), phi_a.size());
                // I think there is a bug in tensor where reshape doesn't automatically change dim
                //
                print("after add\n", A);
            }
        }
        return A_ia_jb - A;
    };
    auto Atwo = two_int(world, phi_0, phi_a);
    A = A + Atwo;
    print(A);
    auto [omega, X] = syev(A);


    print(omega);
    print(X);


    auto t = xao.size() * r_params.num_orbitals();
    auto no = r_params.num_orbitals();

    X_space x_guess(world, t, no);


    for (int i = 0; i < t; i++) {
        auto xt = copy(X(_, i));
        auto mt = xt.reshape(xao.size(), no);
        x_guess.x[i] = transform(world, phi_a, mt);
        // new size is xt column size
    }
    return x_guess;
}


X_space ExcitedTester::test_ao_guess(World &world, ExcitedResponse &calc) {

    print("thresh : ", FunctionDefaults<3>::get_thresh());
    print("k : ", FunctionDefaults<3>::get_k());

    auto phi_0 = copy(world, calc.ground_orbitals);
    print("Ground Orbital norms: ", norm2s_T(world, phi_0));
    auto ao_basis_set = AtomicBasisSet{"6-31g"};
    auto molecule = calc.molecule;
    ao_basis_set.print(molecule);

    vecfuncT ao_vec = vecfuncT(ao_basis_set.nbf(molecule));

    for (int i = 0; i < ao_basis_set.nbf(molecule); ++i) {
        functorT aofunc(new madchem::AtomicBasisFunctor(
                ao_basis_set.get_atomic_basis_function(calc.molecule, i)));
        ao_vec[i] = factoryT(world).functor(aofunc);
    }
    world.gop.fence();
    print("number of ao_basis_functions: ", ao_vec.size());

    madness::print("norm of ao basis functions: ", norm2s_T(world, ao_vec));

    QProjector<double, 3> projector(phi_0);

    // project ground_state from each atomic orbital basis
    std::transform(ao_vec.begin(), ao_vec.end(), ao_vec.begin(),
                   [&](auto &f) { return projector(f); });


    print("overlap between phi0: ", matrix_inner(world, phi_0, ao_vec));
    madness::print("norm of ao basis after projection: ", norm2s_T(world, ao_vec));
    print("number of ao_basis_functions: ", ao_vec.size());


    auto overlap_S = matrix_inner(world, ao_vec, ao_vec);
    print("Overlap S : \n", overlap_S);

    Tensor<double> U, sigma, VT;
    //S=U*sigma*VT
    svd(overlap_S, U, sigma, VT);

    print("singular values of overlap: \n", sigma);
    print("left singular vectors of overlap:\n", U);
    print("right singular vectors of overlap:\n", VT);
    auto xao = transform(world, ao_vec, U);
    auto overlap_SU = matrix_inner(world, xao, xao);
    print("Overlap S after transform : \n", overlap_SU);

    // copy pao to xao response space vector

    auto xc = calc.make_xc_operator(world);
    // now compute the fock potential

    // kinetic energy
    auto F = [&](vector_real_function_3d x) {
        size_t num_orbs = x.size();

        real_derivative_3d Dx(world, 0);
        real_derivative_3d Dy(world, 1);
        real_derivative_3d Dz(world, 2);

        vector_real_function_3d dvx = apply(world, Dx, x);
        vector_real_function_3d dvy = apply(world, Dy, x);
        vector_real_function_3d dvz = apply(world, Dz, x);
        // Apply again for 2nd derivatives
        vector_real_function_3d dvx2 = apply(world, Dx, dvx);
        vector_real_function_3d dvy2 = apply(world, Dy, dvy);
        vector_real_function_3d dvz2 = apply(world, Dz, dvz);

        auto T = (dvx2 + dvy2 + dvz2) * (-0.5);
        real_function_3d v_nuc, v_j0, v_k0, v_xc;
        v_nuc = calc.potential_manager->vnuclear();
        v_nuc.truncate();
        // J^0 x^alpha
        v_j0 = apply(*calc.shared_coulomb_operator, calc.ground_density);
        v_j0.scale(2.0);

        if (calc.xcf.hf_exchange_coefficient() != 1.0) {
            v_xc = xc.make_xc_potential();
        } else {
            // make a zero function
            v_xc = Function<double, 3>(
                    FunctionFactory<double, 3>(world).fence(false).initial_level(1));
        }

        auto hf_exchange_X = zero_functions<double, 3>(world, num_orbs);
        // hf exchange
        std::transform(x.begin(), x.end(), hf_exchange_X.begin(), [&](auto &xi) {
            auto f = Function<double, 3>(FunctionFactory<double, 3>(world));
            std::accumulate(phi_0.begin(), phi_0.end(), f, [&](auto phi0, auto total) {
                auto sum = apply(*calc.shared_coulomb_operator, xi * phi0) * phi0;
                return total + sum;
            });
            return f;
        });
        vector_real_function_3d V0 = zero_functions<double, 3>(world, num_orbs);
        real_function_3d v0 = v_j0 + v_nuc + v_xc;
        V0 = v0 * x;
        V0 += (-1 * hf_exchange_X * calc.xcf.hf_exchange_coefficient());

        return T + V0;
    };

    auto phi_a = copy(world, xao);
    Tensor<double> e_a;
    for (int i = 0; i < 1; i++) {
        auto Fx = F(phi_a);
        auto xFx = matrix_inner(world, phi_a, Fx);
        auto [e, C] = syev(xFx);
        print("eigs  : \n", e);
        phi_a = transform(world, phi_a, C);
        projector(phi_a);
        e_a = e;
    }

    print("ground_orb energies  : \n", calc.ground_energies);
    print("virtual orbital energies  : \n", e_a);
    print("norm of virtual ", norm2s_T(world, phi_a));
    e_a = e_a.flat();
    auto is_positive = [](auto num) { return num > 0; };
    // This seems dumb but i'm removing the negative vectors
    auto first_positive = std::find_if(e_a.ptr(), e_a.ptr() + e_a.size(), is_positive);
    auto num_negative = std::distance(e_a.ptr(), first_positive);
    auto virtual_phi = vector_real_function_3d(phi_a.size() - num_negative);
    Tensor<double> virtual_e(e_a.size() - num_negative);
    std::copy(first_positive, e_a.ptr() + e_a.size(), virtual_e.ptr());
    std::copy(phi_a.begin() + num_negative, phi_a.end(), virtual_phi.begin());
    phi_a = virtual_phi;
    e_a = virtual_e;
    auto SOA = matrix_inner(world, phi_0, phi_a);
    print("Do I remove the vectors that are not completely 0 I hope it's zero\n", SOA);
    print("Sanity CHECK I hope it's diagonal\n", matrix_inner(world, phi_a, F(phi_a)));
    // create a vector of pairs
    std::vector<std::pair<int, int>> ia_indicies;
    print("forming response space pairs");
    print("( i , a )");
    for (int i = 0; i < phi_0.size(); i++) {
        for (int a = 0; a < phi_a.size(); a++) {
            ia_indicies.emplace_back(std::pair<int, int>{i, a});
        }
    }
    int kk = 0;
    Tensor<double> A(ia_indicies.size(), ia_indicies.size());
    for (const auto ia: ia_indicies) {
        std::cout << "( " << ia.first << " , " << ia.second << " )" << std::endl;
        A(kk, kk) = (e_a[ia.second] - calc.ground_energies[ia.first]);
        kk++;
    }
    print("( ij , ab )");
    int ii = 0;
    int jj = 0;
    // create pairs ia, ii,ab
    auto two_int = [&](World &world, const vector_real_function_3d &phi0,
                       const vector_real_function_3d &phia) {
        // The easy case make n*m matrix
        auto tol = FunctionDefaults<3>::get_thresh();
        reconstruct(world, phi0);
        reconstruct(world, phia);

        vector_real_function_3d pairs_ij;
        vector_real_function_3d pairs_ia;
        vector_real_function_3d pairs_ab;

        for (auto &phi_ii: phi0) {
            for (auto &phi_aa: phia) { pairs_ia.push_back(mul_sparse(phi_ii, phi_aa, tol, false)); }
        }
        world.gop.fence();
        truncate(world, pairs_ia);
        vecfuncT Vpairs_ia = apply(world, *calc.shared_coulomb_operator, pairs_ia);
        auto A_ia_jb = matrix_inner(world, pairs_ia, Vpairs_ia);
        for (auto &phi_ii: phi0) {
            for (auto &phi_jj: phi0) { pairs_ij.push_back(mul_sparse(phi_ii, phi_ii, tol, false)); }
        }
        for (auto &phi_aa: phia) {
            for (auto &phi_bb: phia) { pairs_ab.push_back(mul_sparse(phi_aa, phi_aa, tol, false)); }
        }
        world.gop.fence();
        truncate(world, pairs_ij);
        truncate(world, pairs_ab);
        vecfuncT Vpairs_ij = apply(world, *calc.shared_coulomb_operator, pairs_ab);
        auto A_ij_ab = matrix_inner(world, pairs_ia, Vpairs_ij);
        // reshape A_ij_ab  n^2 x m^2
        auto a_len = phi_a.size() * phi_0.size();
        Tensor<double> A(a_len, a_len);
        print(A);
        int kk = 0;
        //[i,j]= will have a phia by phia block.  each row gets rearranged into ij
        for (int ii = 0; ii < phi0.size(); ii++) {
            for (int jj = 0; jj < phi0.size(); jj++) {
                auto ii_ab = A_ij_ab(kk++, _);
                print(ii_ab);
                auto b_i = ii * phia.size();
                auto b_ip1 = (ii * phia.size() + phia.size()) - 1;
                auto b_j = jj * phia.size();
                auto b_jp1 = (jj * phia.size() + phia.size()) - 1;
                print("( ", b_i, " , ", b_ip1, " ) ");
                print("( ", b_j, " , ", b_jp1, " ) ");
                print(Slice(b_i, b_ip1));
                //ii_ab = ii_ab.reshape(phi_a.size(), phi_a.size());
                A(Slice(b_i, b_ip1, 1), Slice(b_j, b_jp1, 1)) =
                        ii_ab.reshape(phi_a.size(), phi_a.size());
                //ii_ab = ii_ab.reshape(phi_a.size(), phi_a.size());
                // I think there is a bug in tensor where reshape doesn't automatically change dim
                //
                print("after add\n", A);
            }
        }
        return A_ia_jb - A;
    };
    auto Atwo = two_int(world, phi_0, phi_a);
    A = A + Atwo;
    print(A);
    auto [omega, X] = syev(A);


    print(omega);
    print(X);


    auto t = xao.size() * calc.r_params.num_orbitals();
    auto no = calc.r_params.num_orbitals();

    X_space x_guess(world, t, no);


    for (int i = 0; i < t; i++) {
        auto xt = copy(X(_, i));
        auto mt = xt.reshape(xao.size(), no);
        x_guess.x[i] = transform(world, phi_a, mt);
        // new size is xt column size
    }


    return x_guess;
}
