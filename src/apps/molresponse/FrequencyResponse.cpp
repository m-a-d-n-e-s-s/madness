//
// Created by adrianhurtado on 2/3/22.
//

#include "FrequencyResponse.hpp"

#include "property.h"


void FrequencyResponse::initialize(World &world) {
    if (world.rank() == 0) { print("FrequencyResponse::initialize()"); }
    Chi = PQ.copy();
}

void FrequencyResponse::iterate(World &world) {
    size_t iter;
    // Variables needed to iterate
    QProjector<double, 3> projector(world, ground_orbitals);
    size_t n = r_params.num_orbitals();// Number of ground state orbitals
    size_t m = r_params.num_states();  // Number of excited states

    real_function_3d v_xc;// For TDDFT
    // the Final protocol should be equal to dconv at the minimum
    const double conv_den = std::max(100 * FunctionDefaults<3>::get_thresh(), r_params.dconv());
    const double relative_max_target =
            std::max(1500 * FunctionDefaults<3>::get_thresh(), 5 * r_params.dconv());
    // m residuals for x and y
    Tensor<double> bsh_residualsX((int(m)));
    Tensor<double> bsh_residualsY((int(m)));
    Tensor<double> density_residuals((int(m)));

    Tensor<double> xij_norms(m, 2 * n);
    Tensor<double> xij_res_norms(m, 2 * n);

    vecfuncT rho_omega_old(m);

    // initialize DFT XC functional operator
    XCOperator<double, 3> xc = make_xc_operator(world);

    // create X space residuals
    X_space residuals = X_space::zero_functions(world, m, n);

    /*
    for (size_t b = 0; b < m; b++) {
        x_vectors.emplace_back(Chi, b);
        x_residuals.emplace_back(residuals, b);
    }
     */
    // create a std vector of XNONLinearsolvers
    response_solver kain_x_space;
    for (size_t b = 0; b < m; b++) {
        kain_x_space.push_back(
                XNonlinearSolver<vector_real_function_3d, double, response_matrix_allocator>(
                        response_matrix_allocator(world, 2 * n), false));
    }
    if (r_params.kain()) {
        for (auto &kain_space_b: kain_x_space) { kain_space_b.set_maxsub(r_params.maxsub()); }
    }
    //
    // We compute with positive frequencies
    if (world.rank() == 0) {
        print("Warning input frequency is assumed to be positive");
        print("Computing at positive frequency omega = ", omega);
    }
    double x_shifts = 0.0;
    double y_shifts = 0.0;
    // if less negative orbital energy + frequency is positive or greater than 0
    if ((ground_energies[long(n) - 1] + omega) >= 0.0) {
        // Calculate minimum shift needed such that \eps + \omega + shift < 0
        print("*** we are shifting just so you know!!!");
        x_shifts = -.05 - (omega + ground_energies[long(n) - 1]);
    }
    auto bsh_x_ops = make_bsh_operators_response(world, x_shifts, omega);
    std::vector<poperatorT> bsh_y_ops;

    bool static_res = (omega == 0.0);
    bool compute_y = not static_res;
    // Negate omega to make this next set of BSH operators \eps - omega
    if (compute_y) {
        omega = -omega;
        bsh_y_ops = make_bsh_operators_response(world, y_shifts, omega);
        omega = -omega;
    }
    vector_real_function_3d rho_omega = make_density(world, Chi);
    converged = false;// Converged flag
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
    PQ = generator(world, *this);
    for (iter = 0; iter < r_params.maxiter(); ++iter) {
        iter_timing.clear();
        // Basic output
        if (r_params.print_level() >= 1) {
            molresponse::start_timer(world);
            if (world.rank() == 0)
                printf("\n   Iteration %d at time %.1fs\n", static_cast<int>(iter), wall_time());
            if (world.rank() == 0) print("-------------------------------------------");
        }
        if (iter < 2 || (iter % 10) == 0) { load_balance_chi(world); }
        if (iter > 0) {
            if (density_residuals.max() > 2) { break; }
            double d_residual = density_residuals.max();
            // Test convergence and set to true
            auto chi_norms = Chi.norm2s();
            auto relative_bsh = copy(bsh_residualsX);
            auto rho_norms = norm2s_T(world, rho_omega);
            std::transform(bsh_residualsX.ptr(), bsh_residualsX.ptr() + bsh_residualsX.size(),
                           chi_norms.ptr(), relative_bsh.ptr(),
                           [](auto bsh, auto norm_chi) { return bsh / norm_chi; });
            auto max_bsh = bsh_residualsX.absmax();
            auto relative_max_bsh = relative_bsh.absmax();
            Tensor<double> polar = -2 * inner(Chi, PQ);
            world.gop.fence();
            // Todo add chi norm and chi_x
            function_data_to_json(j_molresponse, iter, chi_norms, bsh_residualsX, relative_bsh,
                                  xij_norms, xij_res_norms, rho_norms, density_residuals);
            frequency_to_json(j_molresponse, iter, polar);
            world.gop.fence();

            if (r_params.print_level() >= 1) {

                if (world.rank() == 0) {
                    print("thresh: ", FunctionDefaults<3>::get_thresh());
                    print("k: ", FunctionDefaults<3>::get_k());
                    print("Chi Norms at start of iteration: ", iter);
                    print("xij norms: \n", xij_norms);
                    print("xij residual norms: \n", xij_res_norms);
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
                    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Save:");
                }
                if (r_params.plot_all_orbitals()) {
                    plotResponseOrbitals(world, iter, Chi.X, Chi.Y, r_params, ground_calc);
                }
                break;
            }
        }
        auto [new_chi, new_res] = update(world, Chi, xc, bsh_x_ops, bsh_y_ops, projector, x_shifts,
                                         omega, kain_x_space, iter, max_rotation);
        if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
        rho_omega_old = make_density(world, Chi);
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
        if (world.rank() == 0) { print("copy tensors: bshX"); }
        bsh_residualsY = copy(new_res.residual_norms);
        if (world.rank() == 0) { print("copy tensors: bshY"); }
        Chi = new_chi.copy();
        if (world.rank() == 0) { print("copy chi:"); }
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "copy_response_data", "copy_response_data", iter_timing);
        }
        xij_res_norms = new_res.residual.component_norm2s();
        if (world.rank() == 0) { print("computing residuals: xij residuals"); }
        xij_norms = Chi.component_norm2s();
        if (world.rank() == 0) { print("computing chi norms: xij residuals"); }
        density_residuals = norm2s_T(world, (rho_omega - rho_omega_old));
        if (world.rank() == 0) { print("computing residuals: density residuals"); }
        Tensor<double> polar = -2 * inner(Chi, PQ);
        if (world.rank() == 0) { print("computing polarizability:"); }

        if (r_params.print_level() >= 20) {
            auto [eval, evec] = syev(polar);
            if (world.rank() == 0) {
                printf("\n--------Response Properties after %d-------------\n",
                       static_cast<int>(iter));
                print("polarizability");
                print(polar);
                print("polarizability eigenvalues");
                print(eval);
                print("polarizability eigenvectors");
                print(evec);
            }
        }

        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "Iteration Timing", "iter_total", iter_timing);
        }
        time_data.add_data(iter_timing);
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
    }
    compute_and_print_polarizability(world, Chi, PQ, "Converged");
}

auto FrequencyResponse::update(World &world, X_space &chi, XCOperator<double, 3> &xc,
                               std::vector<poperatorT> &bsh_x_ops,
                               std::vector<poperatorT> &bsh_y_ops, QProjector<double, 3> &projector,
                               double &x_shifts, double &omega_n, response_solver &kain_x_space,
                               size_t iteration, const double &maxrotn)
        -> std::tuple<X_space, residuals> {

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    size_t m = chi.num_states();
    bool compute_y = omega_n != 0.0;
    // size_t n = Chi.num_orbitals();
    X_space theta_X = compute_theta_X(world, chi, xc, r_params.calc_type());

    X_space new_chi =
            bsh_update_response(world, theta_X, bsh_x_ops, bsh_y_ops, projector, x_shifts);

    auto [new_res, bsh] = compute_residual(world, chi, new_chi, r_params.calc_type());

    // kain update with temp adjusts temp
    if (r_params.kain() && (iteration > 2)) {
        new_chi = kain_x_space_update(world, chi, new_res, kain_x_space);
    }

    if (false) { x_space_step_restriction(world, chi, new_chi, compute_y, maxrotn); }
    // truncate x

    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "update response", "update", iter_timing);
    }
    //	if not compute y then copy x in to y
    return {new_chi, {new_res, bsh}};

    // print x norms
}

auto FrequencyResponse::bsh_update_response(World &world, X_space &theta_X,
                                            std::vector<poperatorT> &bsh_x_ops,
                                            std::vector<poperatorT> &bsh_y_ops,
                                            QProjector<double, 3> &projector, double &x_shifts)
        -> X_space {
    if (r_params.print_level() >= 1) {
        molresponse::start_timer(world);
        if (world.rank() == 0) { print("--------------- BSH UPDATE RESPONSE------------------"); }
    }

    size_t m = theta_X.X.size();
    size_t n = theta_X.X.size_orbitals();
    bool compute_y = omega != 0.0;
    theta_X.X += theta_X.X * x_shifts;
    theta_X.X += PQ.X;
    theta_X.X = theta_X.X * -2;
    world.gop.fence();
    if (world.rank() == 0) { print("--------------- ThetaX.X------------------"); }

    if (compute_y) {
        theta_X.Y += PQ.Y;
        theta_X.Y = theta_X.Y * -2;
    }
    world.gop.fence();

    if (world.rank() == 0) { print("--------------- ThetaX.Y------------------"); }
    // apply bsh
    X_space bsh_X(world, m, n);
    /*
    bsh_x_ops.insert(bsh_x_ops.end(), std::make_move_iterator(bsh_y_ops.begin()),
                     std::make_move_iterator(bsh_y_ops.end()));
                     */
    bsh_X.X = apply(world, bsh_x_ops, theta_X.X);
    if (world.rank() == 0) { print("--------------- Apply BSH X ------------------"); }
    if (compute_y) { bsh_X.Y = apply(world, bsh_y_ops, theta_X.Y); }

    if (world.rank() == 0) { print("--------------- Apply BSH------------------"); }
    // Project out ground state
    for (size_t i = 0; i < m; i++) bsh_X.X[i] = projector(bsh_X.X[i]);
    if (compute_y) {
        for (size_t i = 0; i < m; i++) { bsh_X.Y[i] = projector(bsh_X.Y[i]); }
    } else {
        bsh_X.Y = bsh_X.X.copy();
    }
    if (world.rank() == 0) { print("--------------- Project BSH------------------"); }

    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "bsh_update", "bsh_update", iter_timing);
    }

    return bsh_X;
}

void FrequencyResponse::frequency_to_json(json &j_mol_in, size_t iter,
                                          const Tensor<double> &polar_ij) {
    json j = {};
    j["iter"] = iter;
    j["polar"] = tensor_to_json(polar_ij);
    auto index = j_mol_in["protocol_data"].size() - 1;
    j_mol_in["protocol_data"][index]["property_data"].push_back(j);
}

void FrequencyResponse::compute_and_print_polarizability(World &world, X_space &Chi, X_space &pq,
                                                         std::string message) {
    Tensor<double> G = -2 * inner(Chi, pq);
    if (world.rank() == 0) {
        print("Polarizability", message);
        print(G);
    }
}

void FrequencyResponse::save(World &world, const std::string &name) {
    // Archive to write everything to
    archive::ParallelOutputArchive ar(world, name.c_str(), 1);

    ar &r_params.archive();
    ar &r_params.tda();
    ar &r_params.num_orbitals();
    ar &r_params.num_states();

    for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++) ar &Chi.X[i][j];
    if (not r_params.tda()) {
        for (size_t i = 0; i < r_params.num_states(); i++)
            for (size_t j = 0; j < r_params.num_orbitals(); j++) ar &Chi.Y[i][j];
    }
}

// Load a response calculation
void FrequencyResponse::load(World &world, const std::string &name) {
    if (world.rank() == 0) { print("FrequencyResponse::load() -state"); }
    // The archive to read from
    archive::ParallelInputArchive ar(world, name.c_str());
    ar &r_params.archive();
    ar &r_params.tda();
    ar &r_params.num_orbitals();
    ar &r_params.num_states();
    Chi = X_space(world, r_params.num_states(), r_params.num_orbitals());
    for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++) ar &Chi.X[i][j];
    world.gop.fence();
    if (not r_params.tda()) {
        for (size_t i = 0; i < r_params.num_states(); i++)
            for (size_t j = 0; j < r_params.num_orbitals(); j++) ar &Chi.Y[i][j];
        world.gop.fence();
    }
}

auto nuclear_generator(World &world, FrequencyResponse &calc) -> X_space {
    auto [gc, molecule, r_params] = calc.get_parameter();
    X_space PQ(world, r_params.num_states(), r_params.num_orbitals());
    auto num_operators = size_t(molecule.natom() * 3);
    auto nuclear_vector = vecfuncT(num_operators);

    for (long atom = 0; atom < molecule.natom(); ++atom) {
        for (long axis = 0; axis < 3; ++axis) {
            functorT func(new madchem::MolecularDerivativeFunctor(molecule, atom, axis));
            nuclear_vector.at(atom * 3 + axis) = functionT(
                    factoryT(world).functor(func).nofence().truncate_on_project().truncate_mode(0));
        }
    }
    PQ.X = vector_to_PQ(world, nuclear_vector, calc.get_orbitals());
    PQ.Y = PQ.X;
    return PQ;
}

auto dipole_generator(World &world, FrequencyResponse &calc) -> X_space {
    auto [gc, molecule, r_params] = calc.get_parameter();
    X_space PQ(world, r_params.num_states(), r_params.num_orbitals());
    vector_real_function_3d dipole_vectors(3);
    size_t i = 0;
    for (auto &d: dipole_vectors) {
        std::vector<int> f(3, 0);
        f[i++] = 1;
        d = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
    }
    //truncate(world, dipole_vectors, true);
    world.gop.fence();
    PQ.X = vector_to_PQ(world, dipole_vectors, calc.get_orbitals());
    PQ.Y = PQ.X;
    if (world.rank() == 0) { print("Made new PQ"); }
    return PQ;
}

auto vector_to_PQ(World &world, const vector_real_function_3d &rhs_operators,
                  const vector_real_function_3d &ground_orbitals) -> response_space {
    response_space rhs(world, rhs_operators.size(), ground_orbitals.size());
    auto orbitals = copy(world, ground_orbitals);
    reconstruct(world, orbitals);
    truncate(world, orbitals);
    QProjector<double, 3> Qhat(world, orbitals);
    int b = 0;
    for (const functionT &pi: rhs_operators) {
        auto op_phi = mul(world, pi, ground_orbitals, true);
        rhs[b] = Qhat(op_phi);
        b++;
    }
    return rhs;
}
//
