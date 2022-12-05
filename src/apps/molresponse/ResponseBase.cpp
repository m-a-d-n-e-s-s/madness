//
// Created by adrianhurtado on 1/24/22.
//

#include "ResponseBase.hpp"


#if defined(__has_include)
#if __has_include(<filesystem>)
#define MADCHEM_HAS_STD_FILESYSTEM
// <filesystem> is not reliably usable on Linux with gcc < 9
#if defined(__GNUC__)
#if __GNUC__ >= 7 && __GNUC__ < 9
#undef MADCHEM_HAS_STD_FILESYSTEM
#endif
#endif
#if defined(MADCHEM_HAS_STD_FILESYSTEM)

#include <filesystem>

#endif
#endif
#endif

// Initializes calculation object for both excited state and frequency dependent
// Copies both the response and ground state
/// Constructs the Base Response
/// \param world
/// \param params
ResponseBase::ResponseBase(World &world, const CalcParams &params)
    : r_params(params.response_parameters), molecule(params.molecule),
      ground_calc(params.ground_calculation), ground_orbitals(ground_calc.orbitals()),
      ground_energies(ground_calc.get_energies()),
      Chi(world, r_params.num_states(), r_params.num_orbitals()) {

    // Broadcast to all the other nodes
    world.gop.broadcast_serializable(r_params, 0);
    world.gop.broadcast_serializable(ground_energies, 0);
    world.gop.broadcast_serializable(molecule, 0);
    xcf.initialize(r_params.xc(), !r_params.spinrestricted(), world, r_params.print_level() >= 3);
    r_params.to_json(j_molresponse);

    // Set the Box Size and Truncate Mode
    FunctionDefaults<3>::set_cubic_cell(-r_params.L(), r_params.L());
    FunctionDefaults<3>::set_truncate_mode(1);
}

/// Checks the polynomial of each function in the ResponseBase
/// First checks the ground state orbitals.  If the orbitals
/// have incorrect k we read from the archive project and truncate
/// We follow by computing the hamiltonian with new orbitals
/// \param world
/// \param thresh
/// \param k
void ResponseBase::check_k(World &world, double thresh, int k) {
    // Boolean to redo ground hamiltonian calculation if
    // ground state orbitals change
    bool redo = false;
    // Verify ground state orbitals have correct k
    if (FunctionDefaults<3>::get_k() != ground_orbitals[0].k()) {
        // Re-read orbitals from the archive (assuming
        // the archive has orbitals stored at a higher
        if (world.rank() == 0) { print("check k: ground orbitals"); }
        // k value than what was previously computed
        ground_calc.read(world);
        if (world.rank() == 0) { print("check k: read ground orbitals"); }
        // k value than what was previously computed
        reconstruct(world, ground_orbitals);
        if (world.rank() == 0) { print("check k: reconstruct ground orbitals"); }
        // Reset correct k (its set in g_params.read)
        FunctionDefaults<3>::set_k(k);
        // Project each ground state to correct k
        for (auto &orbital: ground_orbitals) {
            orbital = project(orbital, FunctionDefaults<3>::get_k(), thresh, false);
        }
        world.gop.fence();
        if (world.rank() == 0) { print("check k: project ground orbitals"); }
        // Clean up a bit
        truncate(world, ground_orbitals);
        if (world.rank() == 0) { print("check k: truncate ground orbitals"); }
        // Now that ground orbitals have correct k lets make the ground density
        // again
        ground_density = make_ground_density(world);
        if (world.rank() == 0) { print("check k: make ground density"); }
        // Ground state orbitals changed, clear old hamiltonian
        redo = true;
    }
    // Recalculate ground state hamiltonian here
    if (redo or !hamiltonian.has_data()) {
        if (world.rank() == 0) { print("check k: re-do hamiltonian"); }
        auto [HAM, HAM_NO_DIAG] = ComputeHamiltonianPair(world);
        if (world.rank() == 0) { print("check k: output hamiltonian"); }
        // TODO this doesn't seem right...
        hamiltonian = HAM;
        ham_no_diag = HAM_NO_DIAG;
    }

    // If we stored the potential, check that too
    if (r_params.store_potential()) {
        if (FunctionDefaults<3>::get_k() != stored_potential[0][0].k()) {
            // Project the potential into correct k
            for (auto &potential_vector: stored_potential) {
                reconstruct(world, potential_vector);
                for (auto &vi: potential_vector) {
                    vi = project(vi, FunctionDefaults<3>::get_k(), thresh, false);
                }
                world.gop.fence();
            }
        }
        if (FunctionDefaults<3>::get_k() != stored_v_coul.k())
            stored_v_coul = project(stored_v_coul, FunctionDefaults<3>::get_k(), thresh, false);
        if (FunctionDefaults<3>::get_k() != stored_v_nuc.k())
            stored_v_nuc = project(stored_v_nuc, FunctionDefaults<3>::get_k(), thresh, false);
    }
    // Don't forget the mask function as well
    if (FunctionDefaults<3>::get_k() != mask.k()) {
        mask = project(mask, FunctionDefaults<3>::get_k(), thresh, false);
        if (world.rank() == 0) { print("check k: project mask"); }
    }
    ::check_k(world, Chi, thresh, k);
    if (world.rank() == 0) { print("check k: project Chi"); }

    // Make sure everything is done before leaving
    world.gop.fence();
}

/// @brief Computes the Hamiltonian from the ground state  orbitals
///  A side effect of this function is that the stored potentials
///  stored_v_coul and stored V nuc get modified
//   Returns both the hamiltonian and the hamiltonian without diagonal
//
//
/// \param world
/// \return
auto ResponseBase::ComputeHamiltonianPair(World &world) const
        -> std::pair<Tensor<double>, Tensor<double>> {
    // Basic output
    if (r_params.print_level() >= 1) molresponse::start_timer(world);
    auto phi = ground_orbitals;
    // Get sizes
    auto num_orbitals = phi.size();
    // Debugging
    // Make the derivative operators in each direction
    real_derivative_3d Dx(world, 0);
    real_derivative_3d Dy(world, 1);
    real_derivative_3d Dz(world, 2);

    // Apply derivatives once, and take inner products
    // according to this formula (faster / less noise):
    //  < f | \nabla^2 | f > = - < \nabla f | \nabla f >
    reconstruct(world, phi);
    std::vector<real_function_3d> fx = apply(world, Dx, phi);
    std::vector<real_function_3d> fy = apply(world, Dy, phi);
    std::vector<real_function_3d> fz = apply(world, Dz, phi);
    compress(world, fx, false);
    compress(world, fy, false);
    compress(world, fz, false);
    world.gop.fence();

    // Construct T according to above formula
    // Note: No negative as the formula above
    // has one as well, so they cancel
    Tensor<double> T = 1.0 / 2.0 *
                       (matrix_inner(world, fx, fx) + matrix_inner(world, fy, fy) +
                        matrix_inner(world, fz, fz));

    // Construct phiVphi
    // v_nuc first
    // TODO Here I am computing the potential using the potential manager.  Should
    // I do
    // TODO earlier on in set_protocol since this maybe used in the response
    // portion as well?
    real_function_3d v_nuc = potential_manager->vnuclear();
    v_nuc.truncate();

    // V_coul next
    // This does not include final multiplication of each orbital
    // 2 is from integrating out spin
    real_function_3d v_coul = 2.0 * Coulomb(world);

    // Clear old stored potentials  Made this mutable
    stored_v_coul.clear();
    stored_v_nuc.clear();

    // If storing potentials, save them here
    if (r_params.store_potential()) {
        stored_v_nuc = copy(v_nuc);
        stored_v_coul = copy(v_coul);
    }

    // v_nuc comes out negative from potential manager, so add it
    real_function_3d v = v_coul + v_nuc;

    // Apply phiVphi to f functions
    auto v_phi0 = v * phi;
    // Clear stored_potential
    stored_potential.clear();
    // ALWAYS DO THIS FOR THE STORED POTENTIAL!!
    // exchange last
    // 'small memory' algorithm from SCF.cc
    auto op = shared_coulomb_operator;

    auto Kphi = zero_functions_compressed<double, 3>(world, int(num_orbitals));

    for (const auto &phi_i: phi) {
        /// Multiplies a function against a vector of functions using sparsity of a
        /// and v[i] --- q[i] = a * v[i]
        auto psif = mul_sparse(world, phi_i, phi, FunctionDefaults<3>::get_thresh());
        truncate(world, psif);
        psif = apply(world, *op, psif);
        truncate(world, psif);
        // Save the potential here if we are saving it
        if (r_params.store_potential()) { stored_potential.push_back(psif); }
        psif = mul_sparse(world, phi_i, psif, FunctionDefaults<3>::get_thresh());
        gaxpy(world, 1.0, Kphi, 1.0, psif);
    }
    // Only use the exchange above if HF:
    Tensor<double> phiVphi;
    if (r_params.xc() == "hf") {
        // Construct phiVphi
        phiVphi = matrix_inner(world, phi, v_phi0) - matrix_inner(world, phi, Kphi);
    } else {// DFT

        XCOperator<double, 3> xcop = make_xc_operator(world);

        real_function_3d v_xc = xcop.make_xc_potential();
        v = v + v_xc;
        auto vf = v * phi;
        if ((*xcop.xc).hf_exchange_coefficient() > 0.0) {
            // XCOperator<double,3>  has member variable xc, which is an
            // xcfunctional which has the hf_exchange_coeff we need here
            gaxpy(world, 1.0, vf, -(*xcop.xc).hf_exchange_coefficient(), Kphi);
        }
        phiVphi = matrix_inner(world, phi, vf);
    }
    // Now create the new_hamiltonian
    auto new_hamiltonian = T + phiVphi;

    for (int64_t i = 0; i < new_hamiltonian.dim(0); i++) {
        for (int64_t j = i + 1; j < new_hamiltonian.dim(1); j++) {
            new_hamiltonian(j, i) = new_hamiltonian(i, j);
        }
    }
    double traceOfHamiltonian(0);
    for (int64_t i = 0; i < new_hamiltonian.dim(0); i++) {
        traceOfHamiltonian += new_hamiltonian(i, i);
    }
    if (world.rank() == 0) {
        print("Trace of Hamiltonian");
        print(traceOfHamiltonian);
    }
    // Save a matrix that is
    // (T+phiVphi) - Lambda * eye
    // Copy new_hamiltonian and zero the diagonal
    auto new_hamiltonian_no_diag = copy(new_hamiltonian);
    for (size_t i = 0; i < num_orbitals; i++) new_hamiltonian_no_diag(long(i), long(i)) = 0.0;

    // End timer
    if (r_params.print_level() >= 1) molresponse::end_timer(world, "   Create grnd ham:");
    return {new_hamiltonian, new_hamiltonian_no_diag};
}

auto ResponseBase::make_ground_density(World &world) const -> functionT {

    auto vsq = square(world, ground_orbitals);
    compress(world, vsq);
    functionT rho = factoryT(world);
    rho.compress();
    for (const auto &vsq_i: vsq) { rho.gaxpy(1.0, vsq_i, 1.0, false); }
    //for (const auto &phi_squared: vsq) rho.gaxpy(2.0, phi_squared, 1.0, false);
    world.gop.fence();
    vsq.clear();
    return rho;
}

// @brief Calculates ground state coulomb potential function
//
//  The coulomb potential of the ground state is
//  \f$ f=\int \sum \frac{\phi_i(r)\phi(r)}{|r-r'|}dr
//
/// \param world
/// \return
auto ResponseBase::Coulomb(World &world) const -> real_function_3d {
    return apply(*shared_coulomb_operator, ground_density).truncate();
}

// TODO  Create apply_operator<T>(f)  generalized function in place of coulomb

auto ResponseBase::make_xc_operator(World &world) const -> XCOperator<double, 3> {
    return {world, r_params.xc(), false, ground_density, ground_density};
}


auto ResponseBase::make_density(World &world, const X_space &chi) const -> vecfuncT {
    auto density = vector_real_function_3d(chi.num_states());
    auto calc_type = r_params.calc_type();
    if (calc_type == "full" || "static") {
        auto r_matrix = to_response_matrix(chi);
        if (world.rank() == 0) { print("make density: to response matrix"); }
        auto r_phi0 = to_response_vector(ground_orbitals);
        if (world.rank() == 0) { print("make density: to response vector"); }
        int b = 0;
        for (auto &rho_b: density) {
            rho_b = dot(world, r_matrix[b], r_phi0);
            b++;
        }

    } else {
        density = transition_densityTDA(world, ground_orbitals, chi.X);
    }
    if (world.rank() == 0) { print("make density: made density"); }
    truncate(world, density);
    if (world.rank() == 0) { print("make density: truncate"); }
    return density;
}

void ResponseBase::load_balance_chi(World &world) {
    molresponse::start_timer(world);
    if (world.size() == 1) return;

    LoadBalanceDeux<3> lb(world);
    real_function_3d v_nuclear;
    v_nuclear = potential_manager->vnuclear();
    for (auto &xi: Chi.X) {
        for (auto &xij: xi) { lb.add_tree(xij, lbcost<double, 3>(1.0, 8.0), false); }
    }
    if (r_params.omega() != 0) {
        for (auto &yi: Chi.Y) {
            for (auto &yij: yi) { lb.add_tree(yij, lbcost<double, 3>(1.0, 8.0), false); }
        }
    }
    world.gop.fence();
    FunctionDefaults<3>::redistribute(
            world,
            lb.load_balance(r_params.loadbalparts()));// 6.0 needs retuning after
    world.gop.fence();
    molresponse::end_timer(world, "Load balancing");
}

auto ResponseBase::make_bsh_operators_response(World &world, double &shift, double &omega) const
        -> std::vector<poperatorT> {
    if (r_params.print_level() >= 1) molresponse::start_timer(world);

    double tol = FunctionDefaults<3>::get_thresh();
    // Sizes inferred from ground and omega
    size_t num_orbitals = ground_energies.size();// number of orbitals
    std::vector<poperatorT> ops(num_orbitals);
    // Run over occupied components

    int p = 0;
    std::for_each(ops.begin(), ops.end(), [&](auto &operator_p) {
        double mu = sqrt(-2.0 * (ground_energies(p++) + omega + shift));
        operator_p = poperatorT(BSHOperatorPtr3D(world, mu, r_params.lo(), tol));
    });
    /*
    for (size_t p = 0; p < num_orbitals; p++) {
        double mu = sqrt(-2.0 * (ground_energies(p) + omega + shift));
        ops[p] = poperatorT(BSHOperatorPtr3D(world, mu, r_params.lo(), tol));
    }
     */
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "make bsh operators response");
    }
    return ops;
    // End timer
}

auto ResponseBase::compute_theta_X(World &world, const X_space &chi,
                                   const XCOperator<double, 3> &xc,
                                   const std::string &calc_type) const -> X_space {

    if (r_params.print_level() >= 1) {
        molresponse::start_timer(world);
        if (world.rank() == 0) { print("------------compute theta x_________"); }
    }
    //     std::cout << "MPI BARRIER 3 " << std::endl;
    //     world.mpi.Barrier();
    bool compute_Y = calc_type == "full";
    X_space Theta_X = X_space(world, chi.num_states(), chi.num_orbitals());
    world.gop.fence();

    //     std::cout << "MPI BARRIER 4 " << std::endl;
    //     world.mpi.Barrier();
    // compute
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    X_space V0X = compute_V0X(world, chi, xc, compute_Y);
    //V0X.truncate();
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "compute_V0X", "compute_V0X", iter_timing);
    }

    if (r_params.print_level() >= 20) { print_inner(world, "xV0x", chi, V0X); }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    X_space E0X(world, chi.num_states(), chi.num_orbitals());
    if (r_params.localize() != "canon") {
        E0X = chi.copy();
        E0X.X = E0X.X * ham_no_diag;
        if (compute_Y) { E0X.Y = E0X.Y * ham_no_diag; }
        if (r_params.print_level() >= 20) { print_inner(world, "xE0x", chi, E0X); }
    }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "compute_E0X", "compute_E0X", iter_timing);
    }

    X_space gamma;
    // compute

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    if (calc_type == "full") {
        gamma = compute_gamma_full(world, {chi, ground_orbitals}, xc);
    } else if (calc_type == "static") {
        gamma = compute_gamma_static(world, {chi, ground_orbitals}, xc);
    } else {
        gamma = compute_gamma_tda(world, {chi, ground_orbitals}, xc);
    }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_compute", "gamma_compute", iter_timing);
    }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    Theta_X = (V0X - E0X) + gamma;
    world.gop.fence();
    Theta_X.truncate();
    //    Theta_X.truncate();
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "compute_ThetaX_add", "compute_ThetaX_add", iter_timing);
    }
    if (r_params.print_level() >= 20) { print_inner(world, "xThetax", chi, Theta_X); }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "compute_ThetaX", "compute_ThetaX", iter_timing);
    }

    return Theta_X;
}


auto ResponseBase::compute_gamma_full(World &world, const gamma_orbitals &density,
                                      const XCOperator<double, 3> &xc) const -> X_space {
    std::shared_ptr<WorldDCPmapInterface<Key<3>>> old_pmap = FunctionDefaults<3>::get_pmap();

    auto [chi_alpha, phi0] = orbital_load_balance(world, density, r_params.loadbalparts());

    size_t num_states = chi_alpha.num_states();
    size_t num_orbitals = chi_alpha.num_orbitals();


    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    // x functions
    // here I create the orbital products for elctron interaction terms
    vecfuncT phi_phi;
    vecfuncT x_phi;
    vecfuncT y_phi;
    functionT temp_J;

    X_space J(world, num_states, num_orbitals);
    response_space j_x(world, num_states, num_orbitals);
    response_space j_y(world, num_states, num_orbitals);

    X_space W = X_space::zero_functions(world, num_states, num_orbitals);
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_zero_functions", "gamma_zero_functions", iter_timing);
    }

    // apply the exchange kernel to rho if necessary
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    // Create Coulomb potential on ground_orbitals
    functionT rho_x_b;
    functionT rho_y_b;

    auto mul_tol = FunctionDefaults<3>::get_thresh();
    // note that x can refer to x or y
    auto compute_j = [&, &phi0 = phi0](const auto &dx) {
        // compute density with response function dx and orbitals phi0
        auto rho_x_b = dot(world, dx, phi0);
        rho_x_b.truncate();
        // apply the coulomb operator to rho_b
        rho_x_b = apply(*shared_coulomb_operator, rho_x_b);
        return mul_sparse(world, rho_x_b, phi0, mul_tol, true);
    };

    // compute j_x = op(rho_x)*phi0

    std::transform(chi_alpha.X.begin(), chi_alpha.X.end(), j_x.begin(), compute_j);
    // compute j_y = op(rho_y)*phi0

    std::transform(chi_alpha.Y.begin(), chi_alpha.Y.end(), j_y.begin(), compute_j);

    J.X = j_x + j_y;
    // TODO is copy better than adding? probably?
    // J.Y=j_x+j_y;
    J.Y = J.X.copy();
    world.gop.fence();

    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "J[omega]", "J[omega]", iter_timing);
    }

    // Create Coulomb potential on ground_orbitals
    if (xcf.hf_exchange_coefficient() != 1.0) {
        auto rho = transition_density(world, phi0, chi_alpha.X, chi_alpha.X);
        auto compute_wx = [&, &phi0 = phi0](auto rho_alpha) {
            auto xc_rho = xc.apply_xc_kernel(rho_alpha);
            return mul(world, xc_rho, phi0);
        };
        if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
        std::transform(rho.begin(), rho.end(), W.X.begin(), compute_wx);
        W.Y = W.X.copy();
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "XC[omega]", "XC[omega]", iter_timing);
        }
    }
    /*
    X_space K1 = X_space::zero_functions(world, num_states, num_orbitals);
    X_space K2 = X_space::zero_functions(world, num_states, num_orbitals);
    auto phi0_c = madness::copy(world, phi0);
    vecfuncT x, y;
    for (size_t b = 0; b < num_states; b++) {
        x = chi_alpha.X[b];
        y = chi_alpha.Y[b];
        K1.X[b] = newK(x, phi0, phi0_c);
        world.gop.fence();
        K1.Y[b] = newK(y, phi0, phi0_c);
        world.gop.fence();
        K2.X[b] = newK(phi0, y, phi0_c);
        world.gop.fence();
        K2.Y[b] = newK(phi0, x, phi0_c);
        world.gop.fence();
    }
    auto K = K1 + K2;
    world.gop.fence();
    if (r_params.print_level() >= 20) { print_inner(world, "old xK1x", chi_alpha, K1); }
    if (r_params.print_level() >= 20) { print_inner(world, "old xK2x", chi_alpha, K2); }
    if (r_params.print_level() >= 20) { print_inner(world, "old xKx", chi_alpha, K); }
     */
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    auto K = response_exchange_multiworld(phi0, chi_alpha, true);
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "K[omega]", "K[omega]", iter_timing);
    }
    if (r_params.print_level() >= 20) { print_inner(world, "old xKx", chi_alpha, K); }
    /*
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    if (r_params.print_level() >= 1) {
        K = response_exchange(phi0, chi_alpha, true);
        molresponse::end_timer(world, "K[omega] multiworld");
        print_inner(world, "new xKx", chi_alpha, K);
    }
     */
    molresponse::start_timer(world);
    X_space gamma(world, num_states, num_orbitals);
    auto c_xc = xcf.hf_exchange_coefficient();
    gamma = 2 * J;
    if (world.rank() == 0) { print("gamma: 2 * J"); }
    gamma += -c_xc * K;
    if (world.rank() == 0) { print("gamma: += -c_xc * K"); }

    if (xcf.hf_exchange_coefficient() != 1.0) {
        gamma += W;
        if (world.rank() == 0) { print("gamma: += W"); }
    }
    //gamma.truncate();
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma add", "gamma_truncate_add", iter_timing);
    }
    // project out ground state
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    QProjector<double, 3> projector(world, phi0);
    for (size_t i = 0; i < num_states; i++) { gamma.X[i] = projector(gamma.X[i]); }
    world.gop.fence();
    for (size_t i = 0; i < num_states; i++) { gamma.Y[i] = projector(gamma.Y[i]); }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_project", "gamma_project", iter_timing);
    }

    if (r_params.print_level() >= 20) {
        molresponse::start_timer(world);
        print_inner(world, "xJx", chi_alpha, J);
        print_inner(world, "xKx", chi_alpha, K);
        print_inner(world, "xWx", chi_alpha, W);
        print_inner(world, "xGammax", chi_alpha, gamma);
        molresponse::end_timer(world, "Print Expectation Creating Gamma:");
    }

    molresponse::start_timer(world);
    J.clear();
    j_x.clear();
    j_y.clear();
    K.clear();
    W.clear();
    chi_alpha.clear();
    phi0.clear();
    molresponse::end_timer(world, "Clear functions and set old pmap");
    if (world.size() > 1) {
        FunctionDefaults<3>::set_pmap(old_pmap);// ! DON'T FORGET !
    }
    //gamma.truncate();
    return gamma;
    // Get sizes
}

auto ResponseBase::compute_gamma_static(World &world, const gamma_orbitals &gammaOrbitals,
                                        const XCOperator<double, 3> &xc) const -> X_space {

    // X contains the response vector that makes up the response gammaOrbitals at a
    // given order

    auto old_pmap = FunctionDefaults<3>::get_pmap();

    auto [xy, phi0] = orbital_load_balance(world, gammaOrbitals, r_params.loadbalparts());


    size_t num_states = xy.num_states();
    size_t num_orbitals = xy.num_orbitals();
    // shallow copy

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    X_space gamma(world, num_states, num_orbitals);
    // here I create the orbital products for elctron interaction terms
    vecfuncT phi_phi;
    vecfuncT x_phi;
    functionT temp_J;

    X_space W = X_space::zero_functions(world, num_states, num_orbitals);
    X_space J(world, num_states, num_orbitals);
    X_space K = X_space::zero_functions(world, num_states, num_orbitals);
    X_space KX = X_space::zero_functions(world, num_states, num_orbitals);
    X_space KY = X_space::zero_functions(world, num_states, num_orbitals);

    //     std::cout << "MPI BARRIER After create Zero functions gamma " << std::endl;
    //     world.mpi.Barrier();

    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_zero_functions", "gamma_zero_functions", iter_timing);
    }

    // apply the exchange kernel to rho if necessary
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    auto rho = make_density(world, xy);
    //auto rho = transition_density(world, phi0, xy.X, xy.X);

    if (r_params.print_level() >= 1) { molresponse::end_timer(world, "compute density J[omega]"); }

    // Create Coulomb potential on ground_orbitals

    /*
    auto compute_jx = [&, &phi0 = phi0](auto rho_alpha) {
        auto temp_J = apply(*shared_coulomb_operator, rho_alpha);
        return mul(world, temp_J, phi0);
    };
     */

    int b = 0;
    for (const auto rho_b: rho) {
        auto temp_J = apply(*shared_coulomb_operator, rho_b);
        J.X[b++] = mul(world, temp_J, phi0);
    }
    //std::transform(rho.begin(), rho.end(), J.X.begin(), compute_jx);
    J.Y = J.X.copy();

    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "J[omega]", "J[omega]", iter_timing);
    }

    if (xcf.hf_exchange_coefficient() != 1.0) {
        auto compute_wx = [&, &phi0 = phi0](auto rho_alpha) {
            auto xc_rho = xc.apply_xc_kernel(rho_alpha);
            return mul(world, xc_rho, phi0);
        };
        if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
        std::transform(rho.begin(), rho.end(), W.X.begin(), compute_wx);
        W.Y = W.X.copy();
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "XC[omega]", "XC[omega]", iter_timing);
        }
    }

    /*
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    std::transform(xy.X.begin(), xy.X.end(), KX.X.begin(),
                   [&](const auto &xi) { return newK(xi, phi0, phi0); });

    std::transform(xy.Y.begin(), xy.Y.end(), KY.X.begin(),
                   [&](const auto &yi) { return newK(phi0, yi, phi0); });


    K = KX + KY;
    world.gop.fence();

    if (r_params.print_level() >= 20) { print_inner(world, "old xK1x", xy, KX); }
    if (r_params.print_level() >= 20) { print_inner(world, "old xK2x", xy, KY); }
    if (r_params.print_level() >= 20) { print_inner(world, "old xKx", xy, K); }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "old K[omega]", "K[omega]", iter_timing);
    }
     */
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    K = response_exchange_multiworld(phi0, xy, false);
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "K[omega]", "K[omega]", iter_timing);
    }
    if (r_params.print_level() >= 20) { print_inner(world, "new static KX", xy, K); }

    /*
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    K = response_exchange(phi0, xy, false);
    if (r_params.print_level() >= 20) { print_inner(world, "new static KX", xy, K); }
    if (r_params.print_level() >= 1) { molresponse::end_timer(world, "new K[omega]"); }
     */
    // for each response state we compute the Gamma response functions
    // trucate all response functions
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    /*
    J.truncate();
    KX.truncate();
    KY.truncate();
    W.truncate();
     */

    // update gamma functions
    gamma = 2 * J - K * xcf.hf_exchange_coefficient() + W;
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_truncate_add", "gamma_truncate_add", iter_timing);
    }

    // project out ground state
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    QProjector<double, 3> projector(world, phi0);
    for (size_t i = 0; i < num_states; i++) { gamma.X[i] = projector(gamma.X[i]); }
    gamma.Y = gamma.X.copy();
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_project", "gamma_project", iter_timing);
    }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    J.clear();
    K.clear();
    W.clear();
    xy.clear();
    phi0.clear();

    if (world.size() > 1) {
        FunctionDefaults<3>::set_pmap(old_pmap);// ! DON'T FORGET !
    }

    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_clear_functions", "gamma_clear_functions",
                               iter_timing);
    }
    // Done
    // gamma.truncate();
    return gamma;
    // Get sizes
}

auto ResponseBase::compute_gamma_tda(World &world, const gamma_orbitals &density,
                                     const XCOperator<double, 3> &xc) const -> X_space {

    auto [d_alpha, phi0] = orbital_load_balance(world, density, r_params.loadbalparts());
    std::shared_ptr<WorldDCPmapInterface<Key<3>>> oldpmap = FunctionDefaults<3>::get_pmap();
    size_t num_states = d_alpha.num_states();
    size_t num_orbitals = d_alpha.num_orbitals();

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    X_space gamma(world, num_states, num_orbitals);
    // x functions
    vector_real_function_3d phi_phi;
    real_function_3d temp_J;
    response_space J(world, num_states, num_orbitals);
    response_space k1_x(world, num_states, num_orbitals);
    response_space W(world, num_states, num_orbitals);
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_zero_functions", "gamma_zero_functions", iter_timing);
    }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    auto rho = transition_densityTDA(world, phi0, d_alpha.X);

    auto compute_jx = [&, &phi0 = phi0](auto rho_alpha) {
        auto temp_J = apply(*shared_coulomb_operator, rho_alpha);
        temp_J.truncate();
        return mul(world, temp_J, phi0);
    };

    std::transform(rho.begin(), rho.end(), J.begin(), compute_jx);
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "J[omega]", "J[omega]", iter_timing);
    }

    // Create Coulomb potential on ground_orbitals
    if (xcf.hf_exchange_coefficient() != 1.0) {
        auto compute_wx = [&, &phi0 = phi0](auto rho_alpha) {
            auto xc_rho = xc.apply_xc_kernel(rho_alpha);
            return mul(world, xc_rho, phi0);
        };
        if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
        // for every transition density apply the exchange kernel and multiply the
        // vector of orbitals
        std::transform(rho.begin(), rho.end(), W.begin(), compute_wx);
        W = W.copy();

        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "XC[omega]", "XC[omega]", iter_timing);
        }
    }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    for (size_t b = 0; b < num_states; b++) {
        vecfuncT x;
        x = d_alpha.X[b];
        k1_x[b] = newK(x, phi0, phi0);
    }

    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "K[omega]", "K[omega]", iter_timing);
    }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    k1_x.truncate_rf();
    J.truncate_rf();
    W.truncate_rf();

    gamma.X = (J * 2) - k1_x * xcf.hf_exchange_coefficient() + W;
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_truncate_add", "gamma_truncate_add", iter_timing);
    }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    QProjector<double, 3> projector(world, ground_orbitals);
    for (size_t i = 0; i < num_states; i++) {
        gamma.X[i] = projector(gamma.X[i]);
        truncate(world, gamma.X[i]);
    }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_project", "gamma_project", iter_timing);
    }

    if (r_params.print_level() >= 20) {
        print("------------------------ Gamma Functions Norms  ------------------");
        print("Gamma X norms");
        print(gamma.X.norm2());
    }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    J.clear();
    k1_x.clear();
    W.clear();

    d_alpha.clear();
    phi0.clear();

    if (world.size() > 1) {
        FunctionDefaults<3>::set_pmap(oldpmap);// ! DON'T FORGET !
    }

    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "gamma_clear_functions", "gamma_clear_functions",
                               iter_timing);
    }
    // Done
    world.gop.fence();
    return gamma;
}

auto ResponseBase::compute_lambda_X(World &world, const X_space &chi, XCOperator<double, 3> &xc,
                                    const std::string &calc_type) const -> X_space {
    // compute
    bool compute_Y = calc_type == "full";

    X_space Lambda_X;// = X_space(world, chi.num_states(), chi.num_orbitals());

    X_space F0X = compute_F0X(world, chi, xc, compute_Y);
    X_space Chi_truncated = chi.copy();
    Chi_truncated.truncate();
    if (r_params.print_level() >= 5) {
        print("---------------Lambda ----------------");
        print("<X|F0|X>");
        print(inner(Chi_truncated, F0X));
    }
    // put it all together

    X_space E0X = Chi_truncated.copy();
    E0X.truncate();
    E0X.X = E0X.X * hamiltonian;

    if (compute_Y) { E0X.Y = E0X.Y * hamiltonian; }
    if (r_params.print_level() >= 20) {
        print("<X|E0|X>");
        print(inner(Chi_truncated, E0X));
    }

    // put it all together
    X_space gamma;

    // compute
    if (calc_type == "full") {
        gamma = compute_gamma_full(world, {chi, ground_orbitals}, xc);
    } else if (calc_type == "static") {
        gamma = compute_gamma_static(world, {chi, ground_orbitals}, xc);
    } else {
        gamma = compute_gamma_tda(world, {chi, ground_orbitals}, xc);
    }
    if (r_params.print_level() >= 5) {
        print("<X|Gamma|X>");
        print(inner(Chi_truncated, gamma));
    }

    Lambda_X = (F0X - E0X) + gamma;
    Lambda_X.truncate();

    if (r_params.print_level() >= 5) {
        print("<X|Lambda_truncated|X>");
        print(inner(Chi_truncated, Lambda_X));
    }

    return Lambda_X;
}

auto ResponseBase::compute_response_potentials(World &world, const X_space &chi,
                                               XCOperator<double, 3> &xc,
                                               const std::string &calc_type) const
        -> std::tuple<X_space, X_space, X_space> {
    // compute
    bool compute_Y = calc_type == "full";

    // first compute kinetic energy piece

    size_t m = chi.num_states();
    size_t n = chi.num_orbitals();
    X_space chi_copy = chi.copy();

    molresponse::start_timer(world);
    X_space T0X = X_space(world, m, n);
    T0X.X = T(world, chi_copy.X);
    if (compute_Y) { T0X.Y = T(world, chi_copy.Y); }
    if (r_params.print_level() >= 20) {
        print("inner <X|T0|X>");
        print(inner(chi_copy, T0X));
    }
    molresponse::end_timer(world, "TX", "TX", iter_timing);


    molresponse::start_timer(world);
    X_space E0X = chi_copy.copy();
    E0X.X = E0X.X * hamiltonian;
    if (compute_Y) { E0X.Y = E0X.Y * hamiltonian; }
    molresponse::end_timer(world, "E0X", "E0X", iter_timing);


    X_space V0X = compute_V0X(world, chi_copy, xc, compute_Y);

    // put it all together
    X_space gamma;
    // compute
    if (calc_type == "full") {
        gamma = compute_gamma_full(world, {chi, ground_orbitals}, xc);
    } else if (calc_type == "static") {
        gamma = compute_gamma_static(world, {chi, ground_orbitals}, xc);
    } else {
        gamma = compute_gamma_tda(world, {chi, ground_orbitals}, xc);
    }

    X_space Lambda_X(world, m, n);// = X_space(world, chi.num_states(), chi.num_orbitals());

    Lambda_X = (T0X + V0X - E0X) + gamma;


    return {Lambda_X, V0X, gamma};
}

// Returns the ground state potential applied to functions f
// (V0 f) V0=(Vnuc+J0-K0+W0)
// J0=J[ground_density]
// K0=K[ground_density]f
// EXC0=W[ground_density]
auto ResponseBase::compute_V0X(World &world, const X_space &X, const XCOperator<double, 3> &xc,
                               bool compute_Y) const -> X_space {
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    //     // std::cout << "MPI BARRIER V0X " << std::endl;
    //     // world.mpi.Barrier();
    // Start a timer
    size_t m = X.num_states();
    size_t n = X.num_orbitals();

    X_space V0 = X_space(world, m, n);
    X_space K0 = X_space(world, m, n);

    X_space Chi_copy = X;
    vecfuncT phi0_copy = madness::copy(world, ground_orbitals);
    world.gop.fence();
    Chi_copy.truncate();
    //Chi_copy.truncate();
    truncate(world, phi0_copy);
    // v_nuc first
    real_function_3d v_nuc, v_j0, v_k0, v_xc;

    world.gop.fence();

    if (not r_params.store_potential()) {
        v_nuc = potential_manager->vnuclear();
        //v_nuc.truncate();
    } else {// Already pre-computed
        v_nuc = stored_v_nuc;
    }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "V0_nuc", "V0_nuc", iter_timing);
    }
    // Coulomb Potential J0*f
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    if (not r_params.store_potential()) {
        // "a" is the core type
        // scale rho by 2 TODO
        // J^0 x^alpha
        v_j0 = apply(*shared_coulomb_operator, ground_density);
        v_j0.scale(2.0);
    } else {// Already pre-computed
        v_j0 = stored_v_coul;
    }
    if (r_params.print_level() >= 1) { molresponse::end_timer(world, "J[0]", "J[0]", iter_timing); }

    if (xcf.hf_exchange_coefficient() != 1.0) {
        if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
        v_xc = xc.make_xc_potential();
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "XC[0]", "XC[0]", iter_timing);
        }
    } else {
        // make a zero function
        v_xc = Function<double, 3>(FunctionFactory<double, 3>(world).fence(false).initial_level(1));
    }

    // Intermediaries
    world.gop.fence();

    /*
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    auto phi0_c = copy(world, phi0_copy);
    world.gop.fence();
    int b = 0;
    for (auto &k0x: K0.X) { k0x = newK(phi0_copy, phi0_c, Chi_copy.X[b++]); }
    if (compute_Y) {
        b = 0;
        for (auto &k0x: K0.Y) { k0x = newK(phi0_copy, phi0_c, Chi_copy.Y[b++]); }

    } else {
        K0.Y = K0.X.copy();
    }
    if (r_params.print_level() >= 20) { print_inner(world, "old xK0x", Chi_copy, K0); }
     */
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    /*
    K0 = ground_exchange(phi0_copy, X, compute_Y);
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "old K[0]");
        print_inner(world, "new xK0x", Chi_copy, K0);
    }
     */

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    K0 = ground_exchange_multiworld(phi0_copy, X, compute_Y);
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "new K[0]", "K[0]", iter_timing);
    }
    if (r_params.print_level() >= 20) { print_inner(world, "new xK0x", Chi_copy, K0); }

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    real_function_3d v0 = v_j0 + v_nuc + v_xc;
    auto c_xc = xcf.hf_exchange_coefficient();
    if (compute_Y) {

        V0 = v0 * X;
        if (world.rank() == 0) { print("vox: v0=v0*X"); }
        V0 += -c_xc * K0;
        if (world.rank() == 0) { print("vox: v0+=c_xc*K0"); }

    } else {
        V0.X = v0 * X.X;
        if (world.rank() == 0) { print("vox: v0=v0*X"); }
        V0.X += -c_xc * K0.X;
        if (world.rank() == 0) { print("vox: v0+=c_xc*K0"); }
        V0.Y = V0.X.copy();
    }
    if (r_params.print_level() >= 20) { print_inner(world, "xV0x", Chi_copy, V0); }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "V0_add", "V0_add", iter_timing);
    }

    return V0;
}

// Returns the ground state fock operator applied to functions f
auto ResponseBase::compute_F0X(World &world, const X_space &X, const XCOperator<double, 3> &xc,
                               bool compute_Y) const -> X_space {
    // Debugging output

    molresponse::start_timer(world);
    size_t m = X.num_states();
    size_t n = X.num_orbitals();

    X_space chi_copy = X.copy();
    // chi_copy.truncate();
    X_space F0X = X_space(world, m, n);
    X_space T0X = X_space(world, m, n);
    T0X.X = T(world, chi_copy.X);
    if (compute_Y) { T0X.Y = T(world, chi_copy.Y); }
    if (r_params.print_level() >= 20) {
        print("_________________compute F0X _______________________");
        print("inner <X|T0|X>");
        print(inner(chi_copy, T0X));
    }

    X_space V0X = compute_V0X(world, chi_copy, xc, compute_Y);
    if (r_params.print_level() >= 20) {
        print("_________________compute F0X _______________________");
        print("inner <X|V0|X>");
        print(inner(chi_copy, V0X));
    }

    F0X = T0X + V0X;

    if (r_params.print_level() >= 20) {
        print("_________________compute F0X _______________________");
        print("inner <X|F0|X>");
        print(inner(chi_copy, F0X));
    }

    molresponse::end_timer(world, "F0X:");
    // Done
    return F0X;
}

auto ResponseBase::compute_residual(World &world, const X_space &chi, const X_space &g_chi,
                                    const std::string &calc_type) -> residuals {
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    size_t m = chi.X.size();
    size_t n = chi.X.size_orbitals();
    //	compute residual
    X_space res(world, m, n);
    res = g_chi - chi;
    auto residual_norms = res.norm2s();
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "compute_bsh_residual", "compute_bsh_residual", iter_timing);
    }
    // Next calculate 2-norm of these vectors of differences
    return {res, residual_norms};
}

auto ResponseBase::kain_x_space_update(World &world, const X_space &chi,
                                       const X_space &residual_chi, response_solver &kain_x_space)
        -> X_space {
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    size_t m = chi.num_states();
    size_t n = chi.num_orbitals();
    X_space kain_update(world, m, n);
    response_matrix update(m);
    auto x_vectors = to_response_matrix(chi);
    auto x_residuals = to_response_matrix(residual_chi);
    if (world.rank() == 0) { print("----------------Start Kain Update -----------------"); }
    int b = 0;
    for (auto &kain_xb: kain_x_space) {
        update[b] = kain_xb.update(x_vectors[b], x_residuals[b]);
        b++;
    }
    world.gop.fence();
    kain_update = to_X_space(update);
    if (world.rank() == 0) { print("----------------End Kain Update -----------------"); }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "kain_x_update", "kain_x_update", iter_timing);
    }
    return kain_update;
}

void ResponseBase::x_space_step_restriction(World &world, const X_space &old_Chi, X_space &temp,
                                            bool restrict_y, const double &maxrotn) {
    size_t m = old_Chi.num_states();
    size_t n = old_Chi.num_orbitals();

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    print(maxrotn);
    auto diff = temp - old_Chi;
    auto m_old = to_response_matrix(old_Chi);
    auto m_new = to_response_matrix(temp);
    auto m_diff = to_response_matrix(diff);

    if (world.rank() == 0) { print("----------------Start Step Restriction -----------------"); }
    for (size_t b = 0; b < m; b++) {
        auto step_size = norm2(world, m_diff[b]);
        auto norm_xb = norm2(world, m_old[b]);
        auto max_step = maxrotn * norm_xb;
        if (world.rank() == 0) {
            print("---------------- step restriction :", b, " ------------------");
            if (world.rank() == 0) { print("X[b]: ", norm_xb); }
            if (world.rank() == 0) { print("deltaX[b]: ", step_size); }
            if (world.rank() == 0) { print("max_step = max_rotation*norm_X: ", max_step); }
        }
        if (step_size > max_step) {
            double s = max_step / step_size;
            if (world.rank() == 0) {
                if (r_params.print_level() > 1)
                    print("  restricting step for response-state: ", b, " step size", s);
            }
            gaxpy(world, s, m_new[b], (1.0 - s), m_old[b], false);
        }
    }
    if (world.rank() == 0) { print("----------------End Step Restriction -----------------"); }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "x_space_restriction", "x_space_restriction", iter_timing);
    }
}


void ResponseBase::plotResponseOrbitals(World &world, size_t iteration,
                                        const response_space &x_response,
                                        const response_space &y_response,
                                        ResponseParameters const &responseParameters,
                                        GroundStateCalculation const &g_params) {
    std::filesystem::create_directories("plots/densities");
    std::filesystem::create_directory("plots/orbitals");

    // TESTING
    // get transition density
    // num orbitals
    size_t n = x_response[0].size();
    size_t m = x_response.size();

    real_function_3d rho0 = dot(world, ground_orbitals, ground_orbitals);
    std::vector<real_function_3d> rho1 =
            transition_density(world, ground_orbitals, x_response, y_response);
    std::string dir("xyz");
    // for plot_name size
    size_t buffSize = 500;
    char plot_name[buffSize];
    double Lp = std::min(responseParameters.L(), 24.0);
    // Doing line plots along each axis
    for (int d = 0; d < 3; d++) {
        // print ground_state
        plotCoords plt(d, Lp);
        // plot ground density
        if (iteration == 1) {
            snprintf(plot_name, buffSize, "plots/densities/rho0_%c_0.plt", dir[d]);
            plot_line(plot_name, 5001, plt.lo, plt.hi, rho0);
        }
        for (int i = 0; i < static_cast<int>(n); i++) {
            // print ground_state
            // plot gound_orbitals
            snprintf(plot_name, buffSize, "plots/orbitals/phi0_%c_0_%d.plt", dir[d],
                     static_cast<int>(i));
            plot_line(plot_name, 5001, plt.lo, plt.hi, ground_orbitals[i]);
        }

        for (int b = 0; b < static_cast<int>(m); b++) {
            // plot rho1 direction d state b
            snprintf(plot_name, buffSize, "plots/densities/rho1_%c_%d.plt", dir[d],
                     static_cast<int>(b));
            plot_line(plot_name, 5001, plt.lo, plt.hi, rho1[b]);

            for (int i = 0; i < static_cast<int>(n); i++) {
                // print ground_state
                // plot x function  x_dir_b_i__k_iter
                snprintf(plot_name, buffSize, "plots/orbitals/phix_%c_%d_%d.plt", dir[d],
                         static_cast<int>(b), static_cast<int>(i));
                plot_line(plot_name, 5001, plt.lo, plt.hi, x_response[b][i]);

                // plot y functione  y_dir_b_i__k_iter
                snprintf(plot_name, buffSize, "plots/orbitals/phiy_%c_%d_%d.plt", dir[d],
                         static_cast<int>(b), static_cast<int>(i));
                plot_line(plot_name, 5001, plt.lo, plt.hi, y_response[b][i]);
            }
        }
    }
    world.gop.fence();

    // END TESTING
}


void PlotGroundDensityVTK(World &world, const ResponseBase &calc) {

    auto [ground_calc, molecule, r_params] = calc.get_parameter();
    auto ground_orbitals = calc.get_orbitals();

    if (r_params.plot_initial()) {
        if (world.rank() == 0) print("\n   Plotting ground state densities.\n");
        if (r_params.plot_l() > 0.0)
            do_vtk_plots(world, int(r_params.plot_pts()), r_params.plot_l(), 0,
                         int(r_params.num_orbitals()), molecule, square(world, ground_orbitals),
                         "ground");
        else
            do_vtk_plots(world, int(r_params.plot_pts()), r_params.L() / 2.0, 0,
                         int(r_params.num_orbitals()), molecule, square(world, ground_orbitals),
                         "ground");
    }
}

/// Push back empty json onto "protocol_data" field
/// Writes protocol to that new field
/// Sets up the iteration data json with an empty {}
/// TODO This can work for any madness iteration therefore maybe I should move
/// it \param j \param proto
void protocol_to_json(json &j, const double proto) {
    j["protocol_data"].push_back({});
    auto proto_index = j["protocol_data"].size() - 1;
    j["protocol_data"][proto_index]["proto"] = proto;
    j["protocol_data"][proto_index]["k"] = FunctionDefaults<3>::get_k();
    j["protocol_data"][proto_index]["iter_data"] = {};
}

void ResponseBase::function_data_to_json(json &j_mol_in, size_t iter, const Tensor<double> &x_norms,
                                         const Tensor<double> &x_abs_norms,
                                         const Tensor<double> &x_rel_norms,
                                         const Tensor<double> &xij_norms,
                                         const Tensor<double> &xij_abs_norms,
                                         const Tensor<double> &rho_norms,
                                         const Tensor<double> &rho_abs_norms) {
    json j = {};

    j["iter"] = iter;

    j["x_norms"] = tensor_to_json(x_norms);
    j["x_abs_error"] = tensor_to_json(x_abs_norms);
    j["x_rel_error"] = tensor_to_json(x_rel_norms);

    j["xij_norms"] = tensor_to_json(xij_norms);
    j["xij_abs_error"] = tensor_to_json(xij_abs_norms);

    j["rho_norms"] = tensor_to_json(rho_norms);
    j["rho_abs_error"] = tensor_to_json(rho_abs_norms);


    auto index = j_mol_in["protocol_data"].size() - 1;
    j_mol_in["protocol_data"][index]["iter_data"].push_back(j);
}

void ResponseBase::solve(World &world) {
    // Get start time
    molresponse::start_timer(world);

    // Plotting input orbitals
    if (r_params.plot_initial()) { PlotGroundDensityVTK(world, *this); }
    const auto protocol = r_params.protocol();
    if (world.rank() == 0) {
        print("Response State Calculation for the following protocols");
        print("Protocol: ", protocol);
    }
    bool first_protocol = true;
    for (const auto &iter_thresh: protocol) {
        // We set the protocol and function defaults here for the given threshold of
        set_protocol(world, iter_thresh);
        if (world.rank() == 0) { print("Succesfully set protocol"); }
        // protocol
        if (first_protocol) {
            if (r_params.restart()) {
                if (world.rank() == 0) {
                    print("   Restarting from file:", r_params.restart_file());
                }
                load(world, r_params.restart_file());
                first_protocol = false;
            } else {
                this->initialize(world);
                if (world.rank() == 0) { print("Succesfully initialized "); }
            }
            check_k(world, iter_thresh, FunctionDefaults<3>::get_k());
            if (world.rank() == 0) { print("Succesfully check K first initialization "); }
            first_protocol = false;
        } else {
            check_k(world, iter_thresh, FunctionDefaults<3>::get_k());
            if (world.rank() == 0) { print("Succesfully check K not first initialization "); }
        }
        protocol_to_json(j_molresponse, iter_thresh);
        // Now actually ready to iterate...
        this->iterate(world);
    }
    // At this point we should know if calc converged maybe add a flag to response.json which states if it has
    converged_to_json(j_molresponse);
    if (r_params.plot()) {
        auto r_matrix = to_response_matrix(Chi);
        do_response_orbital_vtk_plots(world, r_params.plot_pts(), r_params.L(), molecule,
                                      ground_orbitals, r_matrix);
        auto response_densities = make_density(world, Chi);
        do_response_density_vtk_plots(world, r_params.plot_pts(), r_params.L(), molecule,
                                      ground_density, response_densities);
    }


    // Plot the response function if desired
}

void check_k(World &world, X_space &Chi, double thresh = FunctionDefaults<3>::get_thresh(),
             int k = FunctionDefaults<3>::get_k()) {
    if (0 != Chi.X.size()) {
        if (FunctionDefaults<3>::get_k() != Chi.X[0].at(0).k()) {
            // Project all x components into correct k

            for (auto &xi: Chi.X) {
                reconstruct(world, xi);
                for (auto &xij: xi) {
                    xij = project(xij, FunctionDefaults<3>::get_k(), thresh, false);
                }
                world.gop.fence();
            }
            for (auto &yi: Chi.Y) {
                reconstruct(world, yi);
                for (auto &yij: yi) {
                    yij = project(yij, FunctionDefaults<3>::get_k(), thresh, false);
                }
                world.gop.fence();
            }
            Chi.truncate();
        }
    }
}

///  @brief Adds random noise to functions in response space
///
///
///
/// \param world
/// \param f
/// \param magnitude
/// \return
auto add_randomness(World &world, const response_space &f, double magnitude) -> response_space {
    // Copy input functions
    response_space f_copy = f.copy();

    // Lambda function to add in noise
    auto noise = [](const Key<3> &key, Tensor<double> &x) mutable {
        Tensor<double> y(x.size());
        y.fillrandom();
        // y.scale(magnitude);
        y.scale(1e3);
        x = x + y;
        // x(0,0,0) += y(0,0,0)-0.5;
    };
    // TODO
    // Go through each function in f_copy and add in random noise

    for (auto &fi: f_copy) {
        for (auto &fij: fi) { fij.unaryop(noise); }
    }
    // Done
    return f_copy;
}

/***
 * @brief normalize a single response space
 *
 * \note a single response space consists of only x or y states
 *
 *
 * @param world
 * @param f
 */
void normalize(World &world, response_space &f) {
    // Run over rows
    for (auto &fi: f) {
        double norm = inner(fi, fi);
        norm = sqrt(norm);
        // Doing this to deal with zero functions.
        // Maybe not smrt.
        if (norm == 0) continue;
        // And scale
        fi = fi * (1.0 / norm);
    }
}

void normalize(World &world, X_space &Chi) {
    // Run over rows

    for (size_t i = 0; i < Chi.num_states(); i++) {
        // Get the normalization constant
        // (Sum included inside inner)
        double norm_x = inner(Chi.X[i], Chi.X[i]);
        double norm_y = inner(Chi.Y[i], Chi.Y[i]);
        double norm = sqrt(norm_x - norm_y);
        // Doing this to deal with zero functions.
        // Maybe not smrt.
        if (norm == 0) continue;
        Chi.X[i] = Chi.X[i] * (1.0 / norm);
        Chi.Y[i] = Chi.Y[i] * (1.0 / norm);
    }
}

auto solid_harmonics(World &world, int n) -> std::map<std::vector<int>, real_function_3d> {
    // Container to return
    std::map<std::vector<int>, real_function_3d> result;

    // Create the basic x, y, z, constant and zero
    real_function_3d x = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{1, 0, 0})));
    real_function_3d y = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{0, 1, 0})));
    real_function_3d z = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{0, 0, 1})));
    real_function_3d c = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{0, 0, 0})));
    real_function_3d zero = real_factory_3d(world);

    // Add in first few, since they're simple
    // Assuming n >= 1
    result[std::vector<int>{0, 0}] = copy(c);
    result[std::vector<int>{0, -1}] = zero;
    result[std::vector<int>{0, 1}] = zero;
    result[std::vector<int>{-1, 0}] = zero;

    // Generate the solid harmonics recursively from here
    for (int l = 0; l < n; l++) {
        // Calculate ends of this row first
        result[std::vector<int>{l + 1, l + 1}] =
                sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
                (x * result[std::vector<int>{l, l}] -
                 (1 - kronecker(l, 0) * y * result[std::vector<int>{l, -l}]));
        result[std::vector<int>{l + 1, -l - 1}] =
                sqrt(pow(2, kronecker(l, 0) * (2 * l) / (2 * l + 1))) *
                (y * result[std::vector<int>{l, l}] +
                 (1 - kronecker(l, 0) * x * result[std::vector<int>{l, -l}]));

        // Formula below calls for some functions that don't exist.
        // Need zeroes where that would occur
        result[std::vector<int>{l + 1, l + 2}] = zero;
        result[std::vector<int>{l + 1, -l - 2}] = zero;

        // Run over quantum number m
        for (int m = -l; m < l + 1; m++) {
            // Calculate remaining terms
            result[std::vector<int>{l + 1, m}] =
                    1.0 / std::sqrt((l + m + 1) * (l - m + 1)) *
                    ((2 * l + 1) * z * result[std::vector<int>{l, m}] -
                     sqrt((l + m) * (l - m)) * (x * x + y * y + z * z) *
                             result[std::vector<int>{l - 1, m}]);
        }
    }

    // Get rid of any zero functions we added
    for (auto it = result.begin(); it != result.end();) {
        if (it->second.norm2() == 0) it = result.erase(it);
        else
            ++it;
    }

    // Also get rid of the constant
    result.erase(std::vector<int>{0, 0});

    // Done
    return result;
}


vector_real_function_3d transition_density(World &world, const vector_real_function_3d &orbitals,
                                           const response_space &x, const response_space &y) {
    // Get sizes
    // Check sizes and then run the algorithm
    //size_t m = x.size();
    // auto xx = x.copy();
    //  auto yy = y.copy();
    // auto phi0 = copy(world, orbitals);
    //world.gop.fence();

    //xx.truncate_rf();
    // yy.truncate_rf();
    //truncate(world, phi0);
    std::vector<real_function_3d> densities = zero_functions<double, 3>(world, x.size(), true);


    // Return container
    auto compute_density = [&world, &orbitals](const auto &x_alpha, const auto &y_alpha) {
        /*
        for (const auto &xij: x_alpha) {
            print("xij !!", xij.max_depth(), " ", (void *) xij.get_impl().get());
        }
         */
        auto dx = dot(world, x_alpha, orbitals, true);
        /*
        for (const auto &yij: y_alpha) {
            print("yij !!", yij.max_depth(), " ", (void *) yij.get_impl().get());
        }*/
        auto dy = dot(world, orbitals, y_alpha, true);
        return dx + dy;
    };

    /*
    for (const auto &phi_i: orbitals) {
        print("phi_i !!", phi_i.max_depth(), " ", (void *) phi_i.get_impl().get());
    }
     */
    std::transform(x.begin(), x.end(), y.begin(), densities.begin(), compute_density);
    world.gop.fence();
    //truncate(world, densities, FunctionDefaults<3>::get_thresh(), true);
    return densities;
}

/***
 * @brief returns orbitals and response functions with new p map
 *
 * @param world
 * @param psi0
 * @param X
 * @param load_balance
 * @return
 */
auto ResponseBase::orbital_load_balance(World &world, const gamma_orbitals &gammaOrbitals,
                                        const double load_balance) -> gamma_orbitals {

    auto X = std::get<0>(gammaOrbitals);
    auto psi0 = std::get<1>(gammaOrbitals);

    size_t m = X.num_states();
    size_t n = X.num_orbitals();

    if (world.size() > 1) {
        molresponse::start_timer(world);
        LoadBalanceDeux<3> lb(world);

        for (const auto &phi0_i: psi0) { lb.add_tree(phi0_i, lbcost<double, 3>(1.0, 8.0), false); }
        for (const auto &xi: X.X) {
            for (const auto &xij: xi) { lb.add_tree(xij, lbcost<double, 3>(1.0, 8.0), false); }
        }
        for (const auto &yi: X.Y) {
            for (const auto &yij: yi) { lb.add_tree(yij, lbcost<double, 3>(1.0, 8.0), false); }
        }
        world.gop.fence();
        // newpamap is the new pmap just based on the orbitals
        auto new_process_map = lb.load_balance(load_balance);
        // default process map
        // We set the new_process_map
        FunctionDefaults<3>::set_pmap(new_process_map);// set default to be new

        world.gop.fence();
        // copy orbitals using new pmap
        auto X_copy = X.copy(new_process_map, true);
        auto psi0_copy = copy(world, psi0, new_process_map, true);

        molresponse::end_timer(world, "Gamma Orbital Load Balance");
        return {X_copy, psi0_copy};
    } else {
        // return a copy with the same process map since we only have one world
        return {X.copy(), copy(world, psi0)};
    }
}

void ResponseBase::analyze_vectors(World &world, const vecfuncT &x,
                                   const std::string &response_state) {
    molresponse::start_timer(world);
    AtomicBasisSet sto3g("sto-3g");
    vecfuncT ao = project_ao_basis(world, sto3g);

    tensorT C = matrix_inner(world, ao, x);
    int nmo1 = x.size();
    tensorT rsq, dip(3, nmo1);
    {
        functionT frsq = factoryT(world).f(rsquared).initial_level(4);
        // <x r^2 x>
        //<x[i] | r^2 | x[i]>
        rsq = inner(world, x, mul_sparse(world, frsq, x, vtol));
        for (int axis = 0; axis < 3; ++axis) {
            // x y z
            functionT fdip = factoryT(world)
                                     .functor(functorT(new madness::DipoleFunctor(axis)))
                                     .initial_level(4);
            dip(axis, _) = inner(world, x, mul_sparse(world, fdip, x, vtol));
            //<x r^2 x> - <x|x|x>^2-<x|y|x>^2-<x|z|x>^2
            for (int i = 0; i < nmo1; ++i) rsq(i) -= dip(axis, i) * dip(axis, i);
        }
    }
    molresponse::end_timer(world, "Analyze vectors");

    long nmo = x.size();
    size_t ncoeff = 0;
    for (long i = 0; i < nmo; ++i) {
        size_t ncoeffi = x[i].size();
        ncoeff += ncoeffi;
        if (world.rank() == 0 and (r_params.print_level() > 1)) {
            print(response_state + " orbital : ", i);

            printf("ncoeff=%.2e:", (double) ncoeffi);

            printf("center=(%.2f,%.2f,%.2f) : radius=%.2f\n", dip(0, i), dip(1, i), dip(2, i),
                   sqrt(rsq(i)));
            sto3g.print_anal(molecule, C(i, _));
            printf("total number of coefficients = %.8e\n\n", double(ncoeff));
        }
    }
}

auto ResponseBase::project_ao_basis_only(World &world, const AtomicBasisSet &aobasis,
                                         const Molecule &mol) -> vecfuncT {
    vecfuncT ao = vecfuncT(aobasis.nbf(mol));
    for (int i = 0; i < aobasis.nbf(mol); ++i) {
        functorT aofunc(new madchem::AtomicBasisFunctor(aobasis.get_atomic_basis_function(mol, i)));
        ao[i] = factoryT(world).functor(aofunc).truncate_on_project().nofence().truncate_mode(1);
    }
    world.gop.fence();
    truncate(world, ao);
    madness::normalize(world, ao);
    return ao;
}

auto ResponseBase::project_ao_basis(World &world, const AtomicBasisSet &aobasis) -> vecfuncT {
    // Make at_to_bf, at_nbf ... map from atom to first bf on atom, and nbf/atom
    std::vector<int> at_to_bf, at_nbf;
    aobasis.atoms_to_bfn(molecule, at_to_bf, at_nbf);

    return project_ao_basis_only(world, aobasis, molecule);
}

void ResponseBase::output_json() {
    time_data.to_json(j_molresponse);
    auto print_time = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(print_time);
    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y-%m-%d %X");
    print(ss.str());

    nlohmann::json calc_precision = {};
    calc_precision["dconv"] = r_params.dconv();
    calc_precision["thresh"] = FunctionDefaults<3>::get_thresh();
    calc_precision["k"] = FunctionDefaults<3>::get_k();
    j_molresponse["precision"] = calc_precision;
    nlohmann::json timing = {};
    timing["datetime"] = ss.str();
    timing["wall_time"] = wall_time();
    timing["cpu_time"] = cpu_time();
    j_molresponse["time_data"] = timing;
    std::ofstream ofs("response_base.json");
    ofs << std::setw(4) << j_molresponse;
}

void ResponseBase::converged_to_json(json &j) { j["converged"] = converged; }

void ResponseBase::print_inner(World &world, const std::string &name, const X_space &left,
                               const X_space &right) {
    auto m_val = inner(left, right);
    world.gop.fence();
    if (world.rank() == 0) {
        print(name);
        print(m_val);
    }
}


auto transition_densityTDA(World &world, const vector_real_function_3d &orbitals,
                           const response_space &x) -> vector_real_function_3d {

    // Get sizes
    size_t m = x.size();

    auto xx = x;
    auto phi0 = copy(world, orbitals);

    xx.truncate_rf();
    truncate(world, phi0);

    std::vector<real_function_3d> densities = zero_functions<double, 3>(world, m);

    // dot xi with phi0
    auto f = [&world, &phi0](auto xi) { return dot(world, xi, phi0); };

    // for each vector is response space x dot and
    std::transform(x.begin(), x.end(), densities.begin(), f);

    truncate(world, densities);
    world.gop.fence();
    return densities;
}

// Transforms the given matrix of functions according to the give
// transformation matrix. Used to update orbitals / potential
response_space transform(World &world, const response_space &f, const Tensor<double> &U) {
    // Return container
    response_space result;

    // Go element by element
    for (unsigned int i = 0; i < f.size(); i++) {
        // Temp for the result of one row
        std::vector<real_function_3d> temp =
                zero_functions_compressed<double, 3>(world, f[0].size());

        for (unsigned int j = 0; j < f.size(); j++) { gaxpy(world, 1.0, temp, U(j, i), f[j]); }

        // Add to temp to result
        result.push_back(temp);
    }

    result.truncate_rf();

    // Done
    return result;
}

auto transform(World &world, const X_space &x, const Tensor<double> &U) -> X_space {
    // Return container
    X_space result(world, x.num_states(), x.num_orbitals());

    result.X = transform(world, x.X, U);
    result.Y = transform(world, x.Y, U);
    // Done
    return result;
}

auto expectation(World &world, const response_space &A, const response_space &B) -> Tensor<double> {
    // Get sizes
    MADNESS_ASSERT(!A[0].empty());
    MADNESS_ASSERT(A[0].size() == B[0].size());

    size_t dim_1 = A.size();
    size_t dim_2 = A[0].size();
    // Need to take transpose of each input ResponseFunction
    response_space A_t(world, dim_2, dim_1);
    response_space B_t(world, dim_2, dim_1);
    for (size_t i = 0; i < dim_1; i++) {
        for (size_t j = 0; j < dim_2; j++) {
            A_t[j][i] = A[i][j];
            B_t[j][i] = B[i][j];
        }
    }
    // Container for result
    Tensor<double> result(dim_1, dim_1);
    /**
   * @brief
   * [x1 x2 x3]T[x1 x2 x3]
   *
   */
    // Run over dimension two
    // each vector in orbital has dim_1 response functoins associated
    for (size_t p = 0; p < dim_2; p++) { result += matrix_inner(world, A_t[p], B_t[p]); }

    // Done
    return result;
}

void print_norms(World &world, const response_space &f) {

    print(f[0].size());
    Tensor<double> norms(f.size() * f[0].size());
    // Calc the norms
    long i = 0;
    for (const auto &fi: f) {

        for (const auto &fij: fi) { norms(i++) = fij.norm2(); }
    }
    norms = norms.reshape(f.size(), f[0].size());
    // Print em in a smart way
    if (world.rank() == 0) print(norms);
}

response_space select_functions(World &world, response_space f, Tensor<double> &energies, size_t k,
                                size_t print_level) {
    // Container for result
    response_space answer;

    // Debugging output
    if (print_level >= 1) {
        if (world.rank() == 0)
            print("\n   Selecting the", k, "lowest excitation energy components.\n");
    }

    // Get rid of extra functions and save
    // the first k
    while (f.size() > k) f.pop_back();
    answer = f;
    answer.truncate_rf();

    // Get rid of extra energies and save
    // the first k
    energies = energies(Slice(0, k - 1));

    // Basic output
    if (print_level >= 1) {
        if (world.rank() == 0) print("   The selected components have excitation energies:");
        if (world.rank() == 0) print(energies);
    }

    // Done
    return answer;
}

void sort(World &world, Tensor<double> &vals, response_space &f) {
    // Get relevant sizes
    size_t k = vals.size();

    // Copy everything...
    response_space f_copy(f);
    Tensor<double> vals_copy = copy(vals);
    Tensor<double> vals_copy2 = copy(vals);

    // Now sort vals_copy
    std::sort(vals_copy.ptr(), vals_copy.ptr() + vals_copy.size());

    // Now sort the rest of the things, using the sorted energy list
    // to find the correct indices
    for (size_t i = 0; i < k; i++) {
        // Find matching index in sorted vals_copy
        size_t j = 0;
        while (fabs(vals_copy(i) - vals_copy2(j)) > 1e-8 && j < k) j++;

        // Put corresponding function, difference function, value residual and
        // value in the correct place
        f[i] = f_copy[j];
        vals(i) = vals_copy(i);

        // Change the value of vals_copy2[j] to help deal with duplicates?
        vals_copy2(j) = 10000.0;
    }
}

// Sorts the given tensor of eigenvalues and
// response functions
void sort(World &world, Tensor<double> &vals, X_space &f) {
    // Get relevant sizes
    size_t k = vals.size();

    // Copy everything...
    X_space f_copy(f);
    Tensor<double> vals_copy = copy(vals);
    Tensor<double> vals_copy2 = copy(vals);

    // Now sort vals_copy
    std::sort(vals_copy.ptr(), vals_copy.ptr() + vals_copy.size());

    // Now sort the rest of the things, using the sorted energy list
    // to find the correct indices
    for (size_t i = 0; i < k; i++) {
        // Find matching index in sorted vals_copy
        size_t j = 0;
        while (fabs(vals_copy(i) - vals_copy2(j)) > 1e-8 && j < k) j++;

        // Put corresponding function, difference function, value residual and
        // value in the correct place
        f.X[i] = f_copy.X[j];
        f.Y[i] = f_copy.Y[j];

        vals(i) = vals_copy(i);

        // Change the value of vals_copy2[j] to help deal with duplicates?
        vals_copy2(j) = 10000.0;
    }
}

auto gram_schmidt(World &world, const response_space &f) -> response_space {
    // Sizes inferred
    size_t m = f.size();

    // Return container
    response_space result = f.copy();

    // Orthogonalize
    for (size_t j = 0; j < m; j++) {
        // Need to normalize the row
        double norm = norm2(world, result[j]);

        // Now scale each entry
        scale(world, result[j], 1.0 / norm);

        // Project out from the rest of the vectors
        for (size_t k = j + 1; k < m; k++) {
            // Temp function to hold the sum
            // of inner products
            // vmra.h function, line 627
            double temp = inner(result[j], result[k]);

            // Now subtract
            gaxpy(world, 1.0, result[k], -temp, result[j]);
        }
    }
    result.truncate_rf();

    // Done
    return result;
}

auto make_xyz_functions(World &world) -> vector_real_function_3d {
    // Container to return

    // Create the basic x, y, z, constant and zero
    real_function_3d x = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{1, 0, 0})));
    real_function_3d y = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{0, 1, 0})));
    real_function_3d z = real_factory_3d(world).functor(
            real_functor_3d(new MomentFunctor(std::vector<int>{0, 0, 1})));

    std::vector<real_function_3d> funcs = {x, y, z};
    return funcs;
}

// Here i should print some information about the calculation we are
// about to do
response_timing::response_timing() : iter(0) {

    wall_time_data.insert({"iter_total", std::vector<double>(0)});
    wall_time_data.insert({"update", std::vector<double>(0)});
    wall_time_data.insert({"compute_V0X", std::vector<double>(0)});
    wall_time_data.insert({"compute_E0X", std::vector<double>(0)});
    wall_time_data.insert({"compute_ThetaX_add", std::vector<double>(0)});
    wall_time_data.insert({"compute_ThetaX", std::vector<double>(0)});
    wall_time_data.insert({"gamma_compute", std::vector<double>(0)});
    wall_time_data.insert({"gamma_zero_functions", std::vector<double>(0)});
    wall_time_data.insert({"gamma_truncate_add", std::vector<double>(0)});
    wall_time_data.insert({"gamma_project", std::vector<double>(0)});
    wall_time_data.insert({"gamma_clear_functions", std::vector<double>(0)});
    wall_time_data.insert({"J[omega]", std::vector<double>(0)});
    wall_time_data.insert({"XC[omega]", std::vector<double>(0)});
    wall_time_data.insert({"K[omega]", std::vector<double>(0)});
    wall_time_data.insert({"bsh_update", std::vector<double>(0)});
    wall_time_data.insert({"compute_bsh_residual", std::vector<double>(0)});
    wall_time_data.insert({"kain_x_update", std::vector<double>(0)});
    wall_time_data.insert({"x_space_restriction", std::vector<double>(0)});
    wall_time_data.insert({"V0_nuc", std::vector<double>(0)});
    wall_time_data.insert({"J[0]", std::vector<double>(0)});
    wall_time_data.insert({"XC[0]", std::vector<double>(0)});
    wall_time_data.insert({"K[0]", std::vector<double>(0)});
    wall_time_data.insert({"V0_add", std::vector<double>(0)});
    wall_time_data.insert({"make_density_old", std::vector<double>(0)});
    wall_time_data.insert({"make_density_new", std::vector<double>(0)});
    wall_time_data.insert({"copy_response_data", std::vector<double>(0)});
    wall_time_data.insert({"TX", std::vector<double>(0)});
    wall_time_data.insert({"E0X", std::vector<double>(0)});
    wall_time_data.insert({"E0mDX", std::vector<double>(0)});
    wall_time_data.insert({"subspace_reduce", std::vector<double>(0)});
    wall_time_data.insert({"diagonalize_response_matrix", std::vector<double>(0)});

    cpu_time_data.insert({"iter_total", std::vector<double>(0)});
    cpu_time_data.insert({"update", std::vector<double>(0)});
    cpu_time_data.insert({"compute_V0X", std::vector<double>(0)});
    cpu_time_data.insert({"compute_E0X", std::vector<double>(0)});
    cpu_time_data.insert({"compute_ThetaX_add", std::vector<double>(0)});
    cpu_time_data.insert({"compute_ThetaX", std::vector<double>(0)});
    cpu_time_data.insert({"gamma_compute", std::vector<double>(0)});
    cpu_time_data.insert({"gamma_zero_functions", std::vector<double>(0)});
    cpu_time_data.insert({"gamma_truncate_add", std::vector<double>(0)});
    cpu_time_data.insert({"gamma_project", std::vector<double>(0)});
    cpu_time_data.insert({"gamma_clear_functions", std::vector<double>(0)});
    cpu_time_data.insert({"J[omega]", std::vector<double>(0)});
    cpu_time_data.insert({"XC[omega]", std::vector<double>(0)});
    cpu_time_data.insert({"K[omega]", std::vector<double>(0)});
    cpu_time_data.insert({"bsh_update", std::vector<double>(0)});
    cpu_time_data.insert({"compute_bsh_residual", std::vector<double>(0)});
    cpu_time_data.insert({"kain_x_update", std::vector<double>(0)});
    cpu_time_data.insert({"x_space_restriction", std::vector<double>(0)});
    cpu_time_data.insert({"V0_nuc", std::vector<double>(0)});
    cpu_time_data.insert({"J[0]", std::vector<double>(0)});
    cpu_time_data.insert({"XC[0]", std::vector<double>(0)});
    cpu_time_data.insert({"K[0]", std::vector<double>(0)});
    cpu_time_data.insert({"V0_add", std::vector<double>(0)});
    cpu_time_data.insert({"make_density_old", std::vector<double>(0)});
    cpu_time_data.insert({"make_density_new", std::vector<double>(0)});
    cpu_time_data.insert({"copy_response_data", std::vector<double>(0)});
    cpu_time_data.insert({"TX", std::vector<double>(0)});
    cpu_time_data.insert({"E0X", std::vector<double>(0)});
    cpu_time_data.insert({"E0mDX", std::vector<double>(0)});
    cpu_time_data.insert({"subspace_reduce", std::vector<double>(0)});
    cpu_time_data.insert({"diagonalize_response_matrix", std::vector<double>(0)});
}

void response_timing::print_data() {

    for (const auto &[key, value]: wall_time_data) { print(key, " : ", value); }
    for (const auto &[key, value]: cpu_time_data) { print(key, " : ", value); }
}

/**
 * add the pair of s wall_time and cpu_time to the time_data and wall_data maps
 *
 * values.first=wall_time
 * values.second=cpu_time
 * @param values
 */
void response_timing::add_data(std::map<std::string, std::pair<double, double>> values) {
    //   print("ADDING DATA");
    iter++;
    std::for_each(wall_time_data.begin(), wall_time_data.end(), [&values](auto &v) {
        // print(v.first, " : ", values[v.first]);
        v.second.push_back(values[v.first].first);// .first to get first value of pair wall_time
    });

    std::for_each(cpu_time_data.begin(), cpu_time_data.end(), [&values](auto &v) {
        //print(v.first, " : ", values[v.first]);
        v.second.push_back(values[v.first].second);// .first to get first value of pair wall_time
    });
}

void response_timing::to_json(json &j) {

    //::print("FREQUENCY TIME DATA TO JSON");

    j["time_data"] = json();
    j["time_data"]["iterations"] = iter;


    j["time_data"]["wall_time"] = json();
    for (const auto &e: wall_time_data) { j["time_data"]["wall_time"][e.first] = e.second; }

    j["time_data"]["cpu_time"] = json();
    for (const auto &e: cpu_time_data) { j["time_data"]["cpu_time"][e.first] = e.second; }
}
