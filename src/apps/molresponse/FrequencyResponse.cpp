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

    real_function_3d v_xc;
    // the Final protocol should be equal to dconv at the minimum
    const double dconv =
            std::max(FunctionDefaults<3>::get_thresh() * 100, r_params.dconv());//.01 .0001 .1e-5
    auto thresh = FunctionDefaults<3>::get_thresh();
    auto density_target = dconv * std::max(size_t(5.0), molecule.natom());
    const double a_pow {0.655};
    const double b_pow {0.519};
// Last attempt 1.035 2.121

    const double bsh_abs_target = pow(thresh, a_pow) * pow(10, b_pow);//thresh^a*10^b
    // m residuals for x and y
    Tensor<double> bsh_residualsX((int(m)));
    Tensor<double> bsh_relative_residualsX((int(m)));
    Tensor<double> density_residuals((int(m)));

    bool static_res = (omega == 0.0);
    bool compute_y = not static_res;
    int r_vector_size;

    r_vector_size = (compute_y) ? 2 * n : n;

    Tensor<double> v_polar(m, m);
    Tensor<double> polar;
    Tensor<double> res_polar;

    vecfuncT rho_omega_old(m);
    // initialize DFT XC functional operator
    XCOperator<double, 3> xc = make_xc_operator(world);
    // create X space residuals
    X_space residuals = X_space::zero_functions(world, m, n);
    // create a std vector of XNONLinearsolvers
    response_solver kain_x_space;
    for (size_t b = 0; b < m; b++) {
        kain_x_space.emplace_back(response_matrix_allocator(world, r_vector_size), false);
    }
    if (r_params.kain()) {
        for (auto &kain_space_b: kain_x_space) { kain_space_b.set_maxsub(r_params.maxsub()); }
    }
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
    bsh_y_ops = (compute_y) ? make_bsh_operators_response(world, y_shifts, -omega) : bsh_x_ops;
    vector_real_function_3d rho_omega = make_density(world, Chi);
    converged = false;// Converged flag
    auto max_rotation = .5;
    if (thresh >= 1e-2) {
        max_rotation = 2;
    } else if (thresh >= 1e-4) {
        max_rotation = 2 * bsh_abs_target;
    } else if (thresh >= 1e-6) {
        max_rotation = 2 * bsh_abs_target;
    } else if (thresh >= 1e-7) {
        max_rotation = .01;
    }
    functionT mask;
    mask = real_function_3d(real_factory_3d(world).f(mask3).initial_level(4).norefine());
    PQ = generator(world, *this);
    //PQ = PQ * mask;
    //PQ.truncate();
    for (iter = 0; iter < r_params.maxiter(); ++iter) {
        iter_timing.clear();
        if (r_params.print_level() >= 1) {
            molresponse::start_timer(world);
            if (world.rank() == 0)
                printf("\n   Iteration %d at time %.1fs\n", static_cast<int>(iter), wall_time());
            if (world.rank() == 0) print("-------------------------------------------");
        }
        if (iter < 2 || (iter % 10) == 0) { load_balance_chi(world); }
        if (iter > 0) {
            if (density_residuals.max() > 20 && iter > 5) {
                if (world.rank() == 0) { print("d-residual > 20...break"); }
                break;
            }
            double d_residual = density_residuals.max();
            auto chi_norms = (compute_y) ? Chi.norm2s() : Chi.X.norm2();
            auto rho_norms = norm2s_T(world, rho_omega);
            std::transform(bsh_residualsX.ptr(), bsh_residualsX.ptr() + bsh_residualsX.size(),
                           chi_norms.ptr(), bsh_relative_residualsX.ptr(),
                           [](auto bsh, auto norm_chi) { return bsh / norm_chi; });
            auto max_bsh = bsh_residualsX.absmax();
            auto max_relative_bsh = bsh_relative_residualsX.absmax();
            max_rotation = 1.0 * max_bsh;
            world.gop.fence();
            // Todo add chi norm and chi_x
            if (world.rank() == 0) {
                function_data_to_json(j_molresponse, iter, chi_norms, bsh_residualsX, rho_norms,
                                      density_residuals);
                frequency_to_json(j_molresponse, iter, polar, res_polar);
            }
            if (r_params.print_level() >= 1) {
                if (world.rank() == 0) {
                    print("r_params.dconv(): ", r_params.dconv());
                    print("thresh: ", FunctionDefaults<3>::get_thresh());
                    print("k: ", FunctionDefaults<3>::get_k());
                    print("Chi Norms at start of iteration: ", iter);
                    print("Chi_X: ", chi_norms);
                    print("max rotation: ", max_rotation);
                    print("d_residuals : ", density_residuals);
                    print("bsh_residuals : ", bsh_residualsX);
                    print("d_residual_max : ", d_residual);
                    print("density target : ", density_target);
                    print("d_residual_max target : ", density_target);
                    print("bsh_residual_max : ", max_bsh);
                    print("bsh_relative_residual_max : ", max_relative_bsh);
                    print("bsh abs target : ", bsh_abs_target);
                }
            }
            if ((d_residual < density_target) and
                ((max_bsh < bsh_abs_target) or r_params.get<bool>("conv_only_dens"))) {
                converged = true;
            }
            if (converged || iter == r_params.maxiter()) {
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
        if (compute_y) {
            Chi = new_chi.copy();
        } else {
            Chi.X = new_chi.X.copy();
        }
        if (world.rank() == 0) { print("copy chi:"); }
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "copy_response_data", "copy_response_data", iter_timing);
        }
        if (world.rank() == 0) { print("computing chi norms: xij residuals"); }

        auto density_change = madness::sub(world, rho_omega, rho_omega_old, true);
        density_residuals = norm2s_T(world, density_change);
        if (world.rank() == 0) { print("computing residuals: density residuals"); }

        if (world.rank() == 0) { print("computing polarizability:"); }
        if (compute_y) {
            polar = -2 * inner(Chi, PQ);
            res_polar = -2 * inner(new_res.residual, PQ);
        } else {
            polar = -4 * response_space_inner(Chi.X, PQ.X);
            res_polar = -2 * response_space_inner(new_res.residual.X, PQ.X);
        }
        if (r_params.print_level() >= 20) {
            if (world.rank() == 0) {
                printf("\n--------Response Properties after %d-------------\n",
                       static_cast<int>(iter));
                print("polarizability");
                print(polar);
                print("res polarizability");
                print(res_polar);
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
    }
    if (world.rank() == 0) {
        print(" Final energy residuals X:");
        print(bsh_residualsX);
        print(" Final density residuals:");
        print(density_residuals);
    }
    compute_and_print_polarizability(world, Chi, PQ, "Converged");
}

auto FrequencyResponse::update(World &world, X_space &chi, XCOperator<double, 3> &xc,
                               std::vector<poperatorT> &bsh_x_ops,
                               std::vector<poperatorT> &bsh_y_ops, QProjector<double, 3> &projector,
                               double &x_shifts, double &omega_n, response_solver &kain_x_space,
                               size_t iteration, const double &max_rotation)
        -> std::tuple<X_space, residuals> {

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    bool compute_y = omega_n != 0.0;
    auto x = chi.copy();// copy chi
    X_space theta_X = compute_theta_X(world, chi, xc, r_params.calc_type());
    X_space new_chi =
            bsh_update_response(world, theta_X, bsh_x_ops, bsh_y_ops, projector, x_shifts);
    auto [new_res, bsh] = compute_residual(world, chi, new_chi, r_params.calc_type());
    //&& iteration < 7
    if (iteration > 0) {// & (iteration % 3 == 0)) {
        new_chi = kain_x_space_update(world, chi, new_res, kain_x_space);
    }
    if (false) { x_space_step_restriction(world, chi, new_chi, compute_y, max_rotation); }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "update response", "update", iter_timing);
    }
    //	if not compute y then copy x in to y
    return {new_chi, {new_res, bsh}};

    // print x norms
}

auto FrequencyResponse::new_kain_x_space_update(World &world, const X_space &x, const X_space &fx,
                                                response_function_solver &rf_solver) -> X_space {
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    bool compute_y = omega != 0.0;

    size_t m = x.num_states();
    size_t n = x.num_orbitals();
    size_t p = compute_y ? 2 : 1;


    X_space kain_update(world, m, n);
    // step 1 is to place all functions into a single vector

    if (world.rank() == 0) { print("----------------Start Kain Update -----------------"); }

    vector_real_function_3d vect_x(m * n * p);
    vector_real_function_3d vect_fx(m * n * p);
    vector_real_function_3d vect_rx(m * n * p);

    int orb_x;
    int orb_y;
    for (int i = 0; i < m; i++) {
        orb_x = i * p * n;
        for (int j = 0; j < n; j++) {
            vect_x[orb_x + j] = x.X[i][j];
            vect_fx[orb_x + j] = fx.X[i][j];
        }
        if (compute_y) {
            orb_y = orb_x + n;
            for (int j = 0; j < n; j++) {
                vect_x[orb_y + j] = x.Y[i][j];
                vect_fx[orb_y + j] = fx.Y[i][j];
            }
        }
    }
    vect_rx = sub(world, vect_fx, vect_x);
    truncate(world, vect_rx);

    for (int i = 0; i < m; i++) {
        orb_x = i * p * n;
        for (int j = 0; j < n; j++) {
            kain_update.X[i][j] =
                    rf_solver[orb_x + j].update(vect_x[orb_x + j], vect_rx[orb_x + j]);
        }
        if (compute_y) {
            orb_y = orb_x + n;
            for (int j = 0; j < n; j++) {
                kain_update.Y[i][j] =
                        rf_solver[orb_y + j].update(vect_x[orb_y + j], vect_rx[orb_y + j]);
            }
        }
    }
    if (world.rank() == 0) { print("----------------End Kain Update -----------------"); }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "kain_x_update", "kain_x_update", iter_timing);
    }
    return kain_update;
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
    // construct lhs for 2nd order property

    theta_X.X += theta_X.X * x_shifts;
    theta_X.X += PQ.X;
    theta_X.X = theta_X.X * -2;
    world.gop.fence();
    if (compute_y) {
        theta_X.Y += PQ.Y;
        theta_X.Y = theta_X.Y * -2;
    }
    world.gop.fence();
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
    }
    if (world.rank() == 0) { print("--------------- Project BSH------------------"); }
    if (r_params.print_level() >= 1) {
        molresponse::end_timer(world, "bsh_update", "bsh_update", iter_timing);
    }
    if (compute_y) {
        bsh_X.truncate();
    } else {
        bsh_X.X.truncate_rf();
    }
    return bsh_X;
}

void FrequencyResponse::frequency_to_json(json &j_mol_in, size_t iter,
                                          const Tensor<double> &polar_ij,
                                          const Tensor<double> &v_polar_ij) {
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
    PQ.Y = PQ.X.copy();
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
