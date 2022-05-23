//
// Created by adrianhurtado on 2/3/22.
//

#include "FrequencyResponse.hpp"

#include "property.h"


void FrequencyResponse::initialize(World &world) {}
void FrequencyResponse::iterate(World &world) {
    size_t iter;
    // Variables needed to iterate
    QProjector<double, 3> projector(world, ground_orbitals);
    size_t n = r_params.num_orbitals();// Number of ground state orbitals
    size_t m = r_params.num_states();  // Number of excited states

    real_function_3d v_xc;// For TDDFT
    const double dconv = std::max(FunctionDefaults<3>::get_thresh(), r_params.dconv());
    // m residuals for x and y
    Tensor<double> bsh_residualsX(m);
    Tensor<double> bsh_residualsY(m);
    Tensor<double> density_residuals(m);

    vecfuncT rho_omega_old(m);

    // initialize DFT XC functional operator
    XCOperator<double, 3> xc = make_xc_operator(world);

    // create X space residuals
    X_space residuals(world, m, n);
    X_space old_Chi(world, m, n);
    // Create the X space
    // vector of Xvectors
    std::vector<X_vector> Xvector;
    std::vector<X_vector> Xresidual;
    for (size_t b = 0; b < m; b++) {
        Xvector.push_back(X_vector(Chi, b));
        Xresidual.push_back(X_vector(residuals, b));
    }
    // If DFT, initialize the XCOperator<double,3>

    // create a std vector of XNONLinearsolvers
    NonLinearXsolver kain_x_space;
    for (size_t b = 0; b < m; b++) {
        kain_x_space.push_back(XNonlinearSolver<X_vector, double, X_space_allocator>(
                X_space_allocator(world, n), true));
    }
    for (size_t b = 0; b < m; b++) {
        if (r_params.kain()) kain_x_space[b].set_maxsub(r_params.maxsub());
    }
    //
    // We compute with positive frequencies
    print("Warning input frequency is assumed to be positive");
    print("Computing at positive frequency omega = ", omega);
    double x_shifts = 0.0;
    double y_shifts = 0.0;
    // if less negative orbital energy + frequency is positive or greater than 0
    if ((ground_energies[n - 1] + omega) >= 0.0) {
        // Calculate minimum shift needed such that \eps + \omega + shift < 0
        print("*** we are shifting just so you know!!!");
        x_shifts = -.05 - (omega + ground_energies[n - 1]);
    }
    std::vector<poperatorT> bsh_x_ops = make_bsh_operators_response(world, x_shifts, omega);
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

    Tensor<double> maxrotn(m);
    maxrotn.fill(dconv * 100);

    for (iter = 0; iter <= r_params.maxiter(); ++iter) {
        // Basic output
        if (r_params.print_level() >= 1) {
            molresponse::start_timer(world);
            if (world.rank() == 0)
                printf("\n   Iteration %d at time %.1fs\n", static_cast<int>(iter), wall_time());
            if (world.rank() == 0) print("-------------------------------------------");
        }
        if (r_params.print_level() >= 1) {
            if (world.rank() == 0) {
                print("Chi.x norms at start of iteration: ", iter);
                print(Chi.X.norm2());
                print("Chi.y norms at start of iteration ", iter);
                print(Chi.Y.norm2());
            }
        }

        // rho_omega = make_density(world, Chi, compute_y);

        if (iter < 2 || (iter % 10) == 0) { load_balance_chi(world); }

        if (iter > 0) {
            if (density_residuals.max() > 2) { break; }
            double d_residual = density_residuals.max();
            double d_conv = dconv * std::max(size_t(5), molecule.natom());
            // Test convergence and set to true
            print("dconv: ", dconv);
            if ((d_residual < d_conv) and

                ((std::max(bsh_residualsX.absmax(), bsh_residualsY.absmax()) < d_conv * 5.0) or
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
                    PlotGroundandResponseOrbitals(world, iter, Chi.X, Chi.Y, r_params, ground_calc);
                }
                auto rho0 = make_ground_density(world);
                if (r_params.plot()) {
                    do_vtk_plots(world, 200, r_params.L(), molecule, rho0, rho_omega,
                                 ground_orbitals, Chi);
                }
                break;
            }
        }

        auto [new_chi, old_chi, new_res] =
                update(world, Chi, xc, bsh_x_ops, bsh_y_ops, projector, x_shifts, omega,
                       kain_x_space, Xvector,Xresidual,iter,maxrotn);


        rho_omega_old = make_density(world, old_chi);
        rho_omega = make_density(world, new_chi);

        bsh_residualsX = copy(new_res.x);
        bsh_residualsY = copy(new_res.y);

        Chi=new_chi.copy();

        density_residuals = norm2s_T(world, (rho_omega - rho_omega_old));
        maxrotn = (bsh_residualsX + bsh_residualsY) / 4;
        for (size_t i = 0; i < Chi.num_states(); i++) {
            if (maxrotn[i] < r_params.maxrotn()) {
                maxrotn[i] = r_params.maxrotn();
                print("less than maxrotn....set to maxrotn");
            }
        }


        if (world.rank() == 0 and (r_params.print_level() > 2)) {
            print("Density residuals");
            print("dres", density_residuals);
            print("BSH  residuals");
            print("xres", bsh_residualsX);
            print("yres", bsh_residualsY);
            print("maxrotn", maxrotn);
        }


        Tensor<double> polar = -2 * inner(Chi, PQ);

        auto [p, pval] = syev(polar);
        print(p);
        print(pval);


        frequency_to_json(j_molresponse, iter, bsh_residualsX, bsh_residualsY, density_residuals,
                          polar);
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
        compute_and_print_polarizability(world, Chi, PQ, "Converged");
    }
}
std::tuple<X_space, X_space, residuals> FrequencyResponse::update(
        World &world, X_space &Chi, XCOperator<double, 3> &xc, std::vector<poperatorT> &bsh_x_ops,
        std::vector<poperatorT> &bsh_y_ops, QProjector<double, 3> &projector, double &x_shifts,
        double &omega_n, NonLinearXsolver &kain_x_space, vector<X_vector> &Xvector,
        vector<X_vector> &Xresidual,
        size_t iteration, Tensor<double> &maxrotn) {
    size_t m = Chi.num_states();
    bool compute_y = omega_n != 0.0;
    // size_t n = Chi.num_orbitals();

    Tensor<double> errX(m);
    Tensor<double> errY(m);

    X_space theta_X = compute_theta_X(world, Chi, xc, r_params.calc_type());
    // compute residual X_space
    print("BSH update iter = ", iteration);

    X_space new_chi =
            bsh_update_response(world, theta_X, bsh_x_ops, bsh_y_ops, projector, x_shifts);

    auto [new_res, bsh_x, bsh_y] = compute_residual(world, Chi, new_chi, r_params.calc_type());

    // kain update with temp adjusts temp
    if (r_params.kain() && (iteration > 0)) {
        new_chi = kain_x_space_update(world, Chi, new_res, kain_x_space, Xvector, Xresidual);
        if (r_params.print_level() >= 1) {
            compute_and_print_polarizability(world, new_chi, PQ, "<KAIN|PQ>");
        }
    }

    if (iteration > 0) {
        x_space_step_restriction(world, Chi, new_chi, compute_y, maxrotn);
        if (r_params.print_level() >= 1) {
            compute_and_print_polarizability(world, new_chi, PQ, "<STEP_RESTRICTED|PQ>");
        }
    }


    // truncate x
    new_chi.X.truncate_rf();
    // truncate y if compute y
    if (compute_y) new_chi.Y.truncate_rf();
    //	if not compute y then copy x in to y
    return {Chi, new_chi, {new_res, bsh_x, bsh_y}};

    // print x norms
}
void FrequencyResponse::update(World &world, X_space &Chi, X_space &res, XCOperator<double, 3> &xc,
                               std::vector<poperatorT> &bsh_x_ops,
                               std::vector<poperatorT> &bsh_y_ops, QProjector<double, 3> &projector,
                               double &x_shifts, double &omega_n, NonLinearXsolver &kain_x_space,
                               vector<X_vector> &Xvector, vector<X_vector> &Xresidual,
                               Tensor<double> &bsh_residualsX, Tensor<double> &bsh_residualsY,
                               size_t iteration, Tensor<double> &maxrotn) {
    size_t m = Chi.num_states();
    bool compute_y = omega_n != 0.0;
    // size_t n = Chi.num_orbitals();

    Tensor<double> errX(m);
    Tensor<double> errY(m);

    X_space theta_X = compute_theta_X(world, Chi, xc, r_params.calc_type());
    // compute residual X_space
    print("BSH update iter = ", iteration);

    X_space temp = bsh_update_response(world, theta_X, bsh_x_ops, bsh_y_ops, projector, x_shifts);

    res = compute_residual(world, Chi, temp, bsh_residualsX, bsh_residualsY, r_params.calc_type());

    // kain update with temp adjusts temp
    if (r_params.kain() && (iteration > 0)) {
        temp = kain_x_space_update(world, Chi, res, kain_x_space, Xvector, Xresidual);
        if (r_params.print_level() >= 1) {
            compute_and_print_polarizability(world, temp, PQ, "<KAIN|PQ>");
        }
    }

    if (iteration > 0) {
        x_space_step_restriction(world, Chi, temp, compute_y, maxrotn);
        if (r_params.print_level() >= 1) {
            compute_and_print_polarizability(world, temp, PQ, "<STEP_RESTRICTED|PQ>");
        }
    }

    if (r_params.print_level() >= 1) {
        compute_and_print_polarizability(world, temp, PQ, "<BSHX|PQ>");
    }

    // truncate x
    temp.X.truncate_rf();
    // truncate y if compute y
    if (compute_y) temp.Y.truncate_rf();
    //	if not compute y then copy x in to y
    if (!compute_y) temp.Y = temp.X.copy();

    Chi = temp.copy();
    if (r_params.print_level() >= 1) {
        compute_and_print_polarizability(world, Chi, PQ, "<ChiNew|PQ>");
    }
    // print x norms
}
X_space FrequencyResponse::bsh_update_response(World &world, X_space &theta_X,
                                               std::vector<poperatorT> &bsh_x_ops,
                                               std::vector<poperatorT> &bsh_y_ops,
                                               QProjector<double, 3> &projector, double &x_shifts) {
    size_t m = theta_X.X.size();
    size_t n = theta_X.X.size_orbitals();
    bool compute_y = omega != 0.0;

    molresponse::start_timer(world);

    theta_X.X += Chi.X * x_shifts;
    theta_X.X += PQ.X;
    theta_X.X = theta_X.X * -2;
    theta_X.X.truncate_rf();

    if (compute_y) {
        theta_X.Y += PQ.Y;
        theta_X.Y = theta_X.Y * -2;
        theta_X.Y.truncate_rf();
    }
    molresponse::end_timer(world, "Compute residual stuff theta_X");

    // apply bsh
    molresponse::start_timer(world);
    X_space bsh_X(world, m, n);

    bsh_X.X = apply(world, bsh_x_ops, theta_X.X);
    if (compute_y) { bsh_X.Y = apply(world, bsh_y_ops, theta_X.Y); }
    molresponse::end_timer(world, "Apply BSH to theta_X");

    molresponse::start_timer(world);
    // Project out ground state
    for (size_t i = 0; i < m; i++) bsh_X.X[i] = projector(bsh_X.X[i]);

    if (compute_y) {
        for (size_t i = 0; i < m; i++) { bsh_X.Y[i] = projector(bsh_X.Y[i]); }
        bsh_X.truncate();
    } else {
        bsh_X.X.truncate_rf();
        bsh_X.Y = bsh_X.X.copy();
    }
    molresponse::end_timer(world, "Project and truncate BSH_X");

    return bsh_X;
}
void FrequencyResponse::frequency_to_json(json &j_mol_in, size_t iter, const Tensor<double> &res_X,
                                          const Tensor<double> &res_Y,
                                          const Tensor<double> &density_res,
                                          const Tensor<double> &omega) {
    json j = {};
    j["iter"] = iter;
    j["res_X"] = tensor_to_json(res_X);
    j["res_Y"] = tensor_to_json(res_Y);
    j["density_residuals"] = tensor_to_json(density_res);
    j["polar"] = tensor_to_json(omega);
    auto index = j_mol_in["protocol_data"].size() - 1;
    j_mol_in["protocol_data"][index]["iter_data"].push_back(j);
}

void FrequencyResponse::compute_and_print_polarizability(World &world, X_space &Chi, X_space &PQ,
                                                         std::string message) {
    Tensor<double> G = -2 * inner(Chi, PQ);
    if (world.rank() == 0) {
        print("Polarizability", message);
        print(G);
    }
}
void FrequencyResponse::save(World &world, const std::string &name) {


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
X_space nuclear_generator(World &world, FrequencyResponse &calc) {
    auto [gc, molecule, r_params] = calc.get_parameter();
    X_space PQ(world, r_params.num_states(), r_params.num_orbitals());
    auto num_operators = size_t(molecule.natom() * 3);
    auto nuclear_vector = vecfuncT(num_operators);

    for (size_t atom = 0; atom < molecule.natom(); ++atom) {
        for (size_t axis = 0; axis < 3; ++axis) {
            FunctorT func(new MolecularDerivativeFunctor(molecule, atom, axis));
            nuclear_vector.at(atom * 3 + axis) = FunctionT(
                    FactoryT(world).functor(func).nofence().truncate_on_project().truncate_mode(0));
        }
    }
    PQ.X = vector_to_PQ(world, nuclear_vector, calc.get_orbitals(), r_params.lo());
    PQ.Y = PQ.X;
    return PQ;
}
X_space dipole_generator(World &world, FrequencyResponse &calc) {
    auto [gc, molecule, r_params] = calc.get_parameter();
    X_space PQ(world, r_params.num_states(), r_params.num_orbitals());
    vector_real_function_3d dipole_vectors(3);
    size_t i = 0;
    for (auto &d: dipole_vectors) {

        std::vector<int> f(3, 0);
        f[i++] = 1;
        d = real_factory_3d(world).functor(real_functor_3d(new BS_MomentFunctor(f)));
    }
    truncate(world, dipole_vectors, true);
    PQ.X = vector_to_PQ(world, dipole_vectors, calc.get_orbitals(), r_params.lo());
    PQ.Y = PQ.X;
    return PQ;
}
response_space vector_to_PQ(World &world, const vector_real_function_3d &p,
                            const vector_real_function_3d &ground_orbitals, double lo) {

    response_space rhs(world, p.size(), ground_orbitals.size());

    reconstruct(world, ground_orbitals);

    QProjector<double, 3> Qhat(world, ground_orbitals);

    std::vector<real_function_3d> orbitals = ground_orbitals;

    auto f = [&](auto property) {
        auto phat_phi = mul(world, property, ground_orbitals, lo);
        truncate(world, phat_phi);
        // rhs[i].truncate_vec();

        // project rhs vectors for state
        phat_phi = Qhat(phat_phi);
        world.gop.fence();
        return phat_phi;
    };
    std::transform(p.begin(), p.end(), rhs.begin(), f);
    return rhs;
}
//
// Here i should print some information about the calculation we are
// about to do
