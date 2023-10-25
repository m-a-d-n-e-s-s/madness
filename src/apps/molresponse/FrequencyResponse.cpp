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
    madness::QProjector<double, 3> projector(world, ground_orbitals);
    size_t n = r_params.num_orbitals();// Number of ground state orbitals
    size_t m = r_params.num_states();  // Number of excited states

    real_function_3d v_xc;
    const double dconv = std::max(FunctionDefaults<3>::get_thresh() * 10,
                                  r_params.dconv());//.01 .0001 .1e-5
    auto thresh = FunctionDefaults<3>::get_thresh();
    auto density_target = dconv * std::max(size_t(5.0), molecule.natom());
    const double a_pow{0.70254};
    const double b_pow{.73735};
    const double x_residual_target = pow(thresh, a_pow) * pow(10, b_pow);//thresh^a*10^b
    Tensor<double> x_residual((int(m)));
    Tensor<double> density_residuals((int(m)));

    bool static_res = (omega == 0.0);
    bool compute_y = not static_res;
    int r_vector_size;
    all_done = false;
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
    auto max_rotation = .5;
    if (thresh >= 1e-2) {
        max_rotation = 2;
    } else if (thresh >= 1e-4) {
        max_rotation = 2 * x_residual_target;
    } else if (thresh >= 1e-6) {
        max_rotation = 2 * x_residual_target;
    } else if (thresh >= 1e-7) {
        max_rotation = .01;
    }
    PQ = generator(world, *this);
    PQ.truncate();

    vector<bool> converged(Chi.num_states(), false);
    Chi.reset_active();
    // make density for the first time
    auto rho_omega = response_context.compute_density(world, Chi, ground_orbitals,
                                                      vector_real_function_3d(Chi.num_states()), false);

    for (iter = 0; iter < r_params.maxiter(); ++iter) {
        //if (world.rank() == 0) { print("At the start of iterate x", checkx); }
        iter_timing.clear();
        iter_function_data.clear();

        if (r_params.print_level() >= 1) {
            molresponse::start_timer(world);
            if (world.rank() == 0) printf("\n   Iteration %d at time %.1fs\n", static_cast<int>(iter), wall_time());
            if (world.rank() == 0) print("-------------------------------------------");
        }
        if (iter < 2 || (iter % 5) == 0) { load_balance_chi(world); }
        if (iter > 0) {
            if (density_residuals.max() > 20 && iter > 5) {
                if (world.rank() == 0) { print("d-residual > 20...break"); }
                break;
            }

            auto chi_norms = (compute_y) ? Chi.norm2s() : Chi.x.norm2();
            auto rho_norms = madness::norm2s_T(world, rho_omega);

            // Todo add chi norm and chi_x
            if (world.rank() == 0) {
                function_data_to_json(j_molresponse, iter, chi_norms, x_residual, rho_norms, density_residuals);
                frequency_to_json(j_molresponse, iter, polar, res_polar);
            }
            if (r_params.print_level() >= 1) {
                if (world.rank() == 0) {
                    print("r_params.dconv(): ", r_params.dconv());
                    print("thresh: ", FunctionDefaults<3>::get_thresh());
                    print("k: ", FunctionDefaults<3>::get_k());
                    print("Chi Norms at start of iteration: ", iter);
                    print("||X||: ", chi_norms);
                    print("<< XI | XJ >>(omega): \n", polar);
                    print("targets : ||x||", x_residual_target, "    ||delta_rho||", density_target);
                }
            }
            auto check_convergence = [&](auto &ri, auto &di) {
                if (world.rank() == 0) { print("              ", ri, "    ", di); }
                return ((ri < x_residual_target) && (di < density_target));
            };

            for (const auto &b: Chi.active) { converged[b] = check_convergence(x_residual[b], density_residuals[b]); }
            int b = 0;
            auto remove_converged = [&]() {
                Chi.reset_active();
                Chi.active.remove_if([&](auto x) { return converged[b++]; });
            };
            remove_converged();

            if (world.rank() == 0) {
                print("converged", converged);
                print("active", Chi.active);
            }
            b = 0;
            all_done = std::all_of(converged.begin(), converged.end(), [](const auto &ci) { return ci; });
            if (all_done || iter == r_params.maxiter()) {
                // if converged print converged
                if (world.rank() == 0 && all_done and (r_params.print_level() > 1)) { print("\nConverged!\n"); }
                if (r_params.save()) {
                    molresponse::start_timer(world);
                    save(world, r_params.save_file());
                    if (r_params.print_level() >= 1) molresponse::end_timer(world, "Save:");
                }


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
                if (r_params.plot_all_orbitals()) {
                    //plotResponseOrbitals(world, iter, Chi.x, Chi.y, r_params,
                     //                    ground_calc);
                }
#endif
#endif
#endif

                break;
            }
        }
        auto x_inner = ((compute_y) ? 2 : 1) * response_context.inner(Chi, Chi);
        inner_to_json(world, "x", x_inner, iter_function_data);

        auto rho_omega_norm = norm2s_T(world, rho_omega);
        inner_to_json(world, "density_norms", rho_omega_norm, iter_function_data);
        auto [new_chi, new_res, new_rho] =
                update_response(world, Chi, xc, bsh_x_ops, bsh_y_ops, projector, x_shifts, omega, kain_x_space, iter,
                                max_rotation, rho_omega, x_residual, residuals);

        auto old_rho = copy(world, rho_omega);
        rho_omega = copy(world, new_rho);
        // first thing we should do is update the density residuals
        // drho = rho(x)-rho(g(x))
        // new_rho= rho(g(x))

        for (const auto &b: Chi.active) {
            auto drho_b = rho_omega[b] - old_rho[b];
            auto drho_b_norm = drho_b.norm2();
            world.gop.fence();
            density_residuals[b] = drho_b_norm;
        }
        world.gop.fence();

        auto old_density_residual = copy(density_residuals);
        iter_function_data["r_d"] = old_density_residual;


        if (compute_y) {
            Chi = new_chi.copy();
        } else {
            Chi.x = new_chi.x.copy();
        }

        if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
        x_residual = copy(new_res.residual_norms);
        iter_function_data["x_residuals"] = x_residual;
        residuals = new_res.residual.copy();
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "copy_response_data", "copy_response_data", iter_timing);
        }
        inner_to_json(world, "x_residual", x_residual, iter_function_data);

        inner_to_json(world, "density_residuals", old_density_residual, iter_function_data);

        auto dnorm = norm2s_T(world, rho_omega);
        iter_function_data["d"] = dnorm;

        polar = ((compute_y) ? -2 : -4) * response_context.inner(Chi, PQ);
        res_polar = ((compute_y) ? -2 : -4) * response_context.inner(new_res.residual, PQ);
        inner_to_json(world, "alpha", polar, iter_function_data);
        inner_to_json(world, "r_alpha", res_polar, iter_function_data);
        if (r_params.print_level() >= 20) {
            if (world.rank() == 0) {
                printf("\n--------Response Properties after %d-------------\n", static_cast<int>(iter));
                print("<<X,P>> at omega =", omega);
                print(polar);
                print("res r<<X,P>> at omega =", omega);
                print(res_polar);
            }
        }
        if (r_params.print_level() >= 1) {
            molresponse::end_timer(world, "Iteration Timing", "iter_total", iter_timing);
        }
        time_data.add_data(iter_timing);
        function_data.add_data(iter_function_data);
    }
    function_data.add_convergence_targets(FunctionDefaults<3>::get_thresh(), density_target, x_residual_target);
    Chi.reset_active();
    if (world.rank() == 0) print("\n");
    if (world.rank() == 0) print("   Finished Response Calculation ");
    if (world.rank() == 0) print("   ------------------------");
    if (world.rank() == 0) print("\n");

    // Did we converge?
    if (iter == r_params.maxiter() && not all_done) {
        if (world.rank() == 0) print("   Failed to converge. Reason:");
        if (world.rank() == 0) print("\n  ***  Ran out of iterations  ***\n");
    }
    if (world.rank() == 0) {
        print(" Final energy residuals X:");
        print(x_residual);
        print(" Final density residuals:");
        print(density_residuals);
    }
    //compute_and_print_polarizability(world, Chi, PQ, "Converged");
}

auto FrequencyResponse::update_response(World &world, X_space &chi, XCOperator<double, 3> &xc,
                                        std::vector<poperatorT> &bsh_x_ops, std::vector<poperatorT> &bsh_y_ops,
                                        QProjector<double, 3> &projector, double &x_shifts, double &omega_n,
                                        response_solver &kain_x_space, size_t iteration, const double &max_rotation,
                                        const vector_real_function_3d &rho_old, const Tensor<double> &old_residuals,
                                        const X_space &xres_old)
        -> std::tuple<X_space, residuals, vector_real_function_3d> {

    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }

    auto x = chi.copy();
    X_space theta_X = compute_theta_X(world, x, rho_old, xc, r_params.calc_type());
    X_space new_chi = bsh_update_response(world, theta_X, bsh_x_ops, bsh_y_ops, projector, x_shifts);


    inner_to_json(world, "x_new", response_context.inner(new_chi, new_chi), iter_function_data);

    auto [new_res, bsh] = update_residual(world, chi, new_chi, r_params.calc_type(), old_residuals, xres_old);
    inner_to_json(world, "r_x", response_context.inner(new_res, new_res), iter_function_data);
    if (iteration >= 0) {// & (iteration % 3 == 0)) {
        new_chi = kain_x_space_update(world, chi, new_res, kain_x_space);
    }
    inner_to_json(world, "x_update", response_context.inner(new_chi, new_chi), iter_function_data);

    //bool compute_y = r_params.calc_type() == "full";
    //x_space_step_restriction(world, chi, new_chi, compute_y, max_rotation);
    if (r_params.print_level() >= 1) { molresponse::end_timer(world, "update response", "update", iter_timing); }

    auto new_rho = response_context.compute_density(world, new_chi, ground_orbitals, rho_old, true);

    return {new_chi, {new_res, bsh}, new_rho};
}

auto FrequencyResponse::bsh_update_response(World &world, X_space &theta_X, std::vector<poperatorT> &bsh_x_ops,
                                            std::vector<poperatorT> &bsh_y_ops, QProjector<double, 3> &projector,
                                            double &x_shifts) -> X_space {
    if (r_params.print_level() >= 1) { molresponse::start_timer(world); }
    size_t m = theta_X.x.size();
    size_t n = theta_X.x.size_orbitals();
    bool compute_y = omega != 0.0;

    if (compute_y) {
        theta_X += theta_X * x_shifts;
        theta_X += PQ;
        theta_X = -2 * theta_X;
        theta_X.truncate();
    } else {
        theta_X.x += theta_X.x * x_shifts;
        theta_X.x += PQ.x;
        theta_X.x = theta_X.x * -2;
        theta_X.x.truncate_rf();
    }
    // apply bsh
    X_space bsh_X(world, m, n);
    bsh_X.set_active(theta_X.active);
    bsh_X.x = apply(world, bsh_x_ops, theta_X.x);
    if (compute_y) { bsh_X.y = apply(world, bsh_y_ops, theta_X.y); }

    if (compute_y) {
        bsh_X.truncate();
    } else {
        bsh_X.x.truncate_rf();
    }

    auto apply_projector = [&](auto &xi) { return projector(xi); };
    if (compute_y) {
        bsh_X = oop_apply(bsh_X, apply_projector);
    } else {
        for (const auto &i: bsh_X.active) bsh_X.x[i] = projector(bsh_X.x[i]);
    }
    if (r_params.print_level() >= 1) { molresponse::end_timer(world, "bsh_update", "bsh_update", iter_timing); }
    if (compute_y) {
        bsh_X.truncate();
    } else {
        bsh_X.x.truncate_rf();
    }
    return bsh_X;
}

void FrequencyResponse::frequency_to_json(json &j_mol_in, size_t iter, const Tensor<double> &polar_ij,
                                          const Tensor<double> &res_polar_ij) {
    json j = {};
    j["iter"] = iter;
    j["polar"] = tensor_to_json(polar_ij);
    j["res_polar"] = tensor_to_json(res_polar_ij);
    auto index = j_mol_in["protocol_data"].size() - 1;
    j_mol_in["protocol_data"][index]["property_data"].push_back(j);
}

void FrequencyResponse::compute_and_print_polarizability(World &world, X_space &Chi, X_space &pq, std::string message) {
    Tensor<double> G = -2 * inner(Chi, pq);
    if (world.rank() == 0) {
        print("Polarizability", message);
        print(G);
    }
}

void FrequencyResponse::save(World &world, const std::string &name) {
    // Archive to write everything to
    archive::ParallelOutputArchive ar(world, name.c_str(), 1);

    ar & r_params.archive();
    ar & r_params.tda();
    ar & r_params.num_orbitals();
    ar & r_params.num_states();

    for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++) ar & Chi.x[i][j];
    for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++) ar & Chi.y[i][j];
}

// Load a response calculation
void FrequencyResponse::load(World &world, const std::string &name) {
    if (world.rank() == 0) { print("FrequencyResponse::load() -state"); }
    // The archive to read from
    archive::ParallelInputArchive ar(world, name.c_str());
    ar & r_params.archive();
    ar & r_params.tda();
    ar & r_params.num_orbitals();
    ar & r_params.num_states();
    Chi = X_space(world, r_params.num_states(), r_params.num_orbitals());
    for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++) ar & Chi.x[i][j];
    world.gop.fence();
    for (size_t i = 0; i < r_params.num_states(); i++)
        for (size_t j = 0; j < r_params.num_orbitals(); j++) ar & Chi.y[i][j];
    world.gop.fence();
}

auto nuclear_generator(World &world, ResponseBase &calc) -> X_space {
    auto [gc, molecule, r_params] = calc.get_parameter();
    X_space PQ(world, r_params.num_states(), r_params.num_orbitals());
    auto num_operators = size_t(molecule.natom() * 3);
    auto nuclear_vector = vecfuncT(num_operators);

    for (long atom = 0; atom < molecule.natom(); ++atom) {
        for (long axis = 0; axis < 3; ++axis) {
            functorT func(new madchem::MolecularDerivativeFunctor(molecule, atom, axis));
            nuclear_vector.at(atom * 3 + axis) =
                    functionT(factoryT(world).functor(func).nofence().truncate_on_project().truncate_mode(0));
        }
    }
    PQ.x = vector_to_PQ(world, nuclear_vector, calc.get_orbitals());
    PQ.y = PQ.x;
    return PQ;
}

auto dipole_generator(World &world, ResponseBase &calc) -> X_space {
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
    PQ.x = vector_to_PQ(world, dipole_vectors, calc.get_orbitals());
    PQ.y = PQ.x.copy();
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

void QuadraticResponse::load(World &world, const std::string &name) {}
void QuadraticResponse::save(World &world, const std::string &name) {}

void QuadraticResponse::iterate(World &world) {}
void QuadraticResponse::initialize(World &world) {}

auto QuadraticResponse::setup_rhs_BC(World &world, const X_space &PQ) -> std::pair<X_space, X_space> {

    auto pq = to_response_matrix(x_data[1].first.copy());

    auto PQB = create_response_matrix(6, r_params.num_orbitals());
    auto PQC = create_response_matrix(6, r_params.num_orbitals());


    vector<double> times_push_back{3, 2, 1};

    auto new_index = 0;
    auto xb_index = 0;
    for (const auto &xbi: pq) {
        auto copy_val = times_push_back[xb_index];
        for (auto j = 0; j < copy_val; j++) {
            // copy xbi into XB[xbi]
            PQB[new_index] = copy(world, xbi, false);
            new_index++;
        };
        xb_index++;
    }


    world.gop.fence();


    auto xc_index = 0;
    for (int j = 0; j < 3; j++) {
        for (int k = j; k < 3; k++) {
            PQC[xc_index] = copy(world, pq[k], false);
            xc_index++;
        }
    }


    world.gop.fence();
    auto r_PB = to_X_space(PQB);
    auto r_PC = to_X_space(PQC);
    world.gop.fence();
    return {r_PB, r_PC};
}

// Take the second and third x_data and organize them by copying into a new X_space
// xb xxxyyz
// xc xyzyzz

auto QuadraticResponse::setup_XBC(World &world) -> std::pair<X_space, X_space> {

    // copy x_data[1] and x_data[2] into new response matrix holding 6 vectors

    auto xb = to_response_matrix(x_data[1].first.copy());
    auto xc = to_response_matrix(x_data[2].first.copy());

    print("x_data norms x", x_data[1].first.x.norm2());
    print("x_data norms y", x_data[1].first.y.norm2());

    // create new response matrices to hold the organized data
    auto XB = create_response_matrix(6, r_params.num_orbitals());
    auto XC = create_response_matrix(6, r_params.num_orbitals());


    vector<double> times_push_back{3, 2, 1};

    auto new_index = 0;
    auto xb_index = 0;
    // for each vector in xb
    for (const auto &xbi: xb) {
        auto copy_j_times = times_push_back[xb_index];
        for (auto j = 0; j < copy_j_times; j++) {
            // copy xbi into XB[xbi]
            XB[new_index] = copy(world, xbi, false);
            new_index++;
        };
        xb_index++;
    }


    world.gop.fence();


    auto xc_index = 0;

    for (int j = 0; j < 3; j++) {
        for (int k = j; k < 3; k++) {
            XC[xc_index] = copy(world, xc[k], false);
            xc_index++;
        }
    }


    world.gop.fence();
    auto r_XB = to_X_space(XB);
    auto r_XC = to_X_space(XC);
    world.gop.fence();
    return {r_XB, r_XC};
}

X_space QuadraticResponse::compute_gamma_bc_virt(World &world, const X_space &XB, const X_space &XC) {}

Tensor<double> QuadraticResponse::compute_beta_unrelaxed(World &world, const X_space &AB_left, const X_space &AB_right,
                                                         X_space &BA_left, X_space &BA_right) {

    Tensor<double> beta(3, 6);

    auto create_dipole = [&]() {
        vector_real_function_3d dipole_vectors(3);
        size_t i = 0;
        // creates a vector of x y z dipole functions
        for (auto &d: dipole_vectors) {
            std::vector<int> f(3, 0);
            f[i++] = 1;
            d = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
        }
        return dipole_vectors;
    };

    auto dipole_vectors = create_dipole();
    truncate(world, dipole_vectors, true);
    std::vector<std::string> names_i = {"x", "y", "z"};
    std::vector<std::string> names_j = {"xx", "xy", "xz", "yy", "yz", "zz"};
    // for each vector in dipole vectors
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 6; j++) {
            auto beta_AB_x = (dot(world, AB_left.x[j], AB_right.x[j]) * dipole_vectors[i]).trace();
            auto beta_BA_x = (dot(world, BA_left.x[j], BA_right.x[j]) * dipole_vectors[i]).trace();
            auto beta_AB_y = (dot(world, AB_left.y[j], AB_right.y[j]) * dipole_vectors[i]).trace();
            auto beta_BA_y = (dot(world, BA_left.y[j], BA_right.y[j]) * dipole_vectors[i]).trace();
            if (r_params.print_level() > 0) {
                print(names_i[i], names_j[j], " beta_AB_x ", beta_AB_x, " beta_BA_x ", beta_BA_x, " beta_AB_y ",
                      beta_AB_y, " beta_BA_y ", beta_BA_y, "\n");
            }
            beta(i, j) = beta_AB_x + beta_AB_y + beta_BA_x + beta_BA_y;
        }
    }
    if (world.rank() == 0) {
        cout << "beta unrelaxed " << endl;
        cout << beta << endl;
    }

    return beta;
}

Tensor<double> QuadraticResponse::compute_beta(World &world) {

    // construct an X_space containing phi0 copies

    // bsh_X = oop_apply(bsh_X, apply_projector);
    QProjector<double, 3> projector(world, ground_orbitals);
    auto apply_projector = [&](auto &xi) { return projector(xi); };

    auto perturbation_A = generator(world, *this);
    auto XA = -1.0 * x_data[0].first.copy();

    // first step to compute beta is to construct the X_space representations of the virt/virt and occ/occ blocks of gamma

    auto [XB, XC] = setup_XBC(world);
    X_space phi0 = X_space(world, XB.num_states(), XC.num_orbitals());
    for (auto i = 0; i < phi0.num_states(); i++) {
        phi0.x[i] = copy(world, ground_orbitals);
        phi0.y[i] = copy(world, ground_orbitals);
    }


    auto [AB_left, AB_right, BA_left, BA_right] = compute_gamma_unrelaxed_block(world, XB, XC);

    auto beta_unrelaxed = compute_beta_unrelaxed(world, AB_left, AB_right, BA_left, BA_right);


    auto second_order_x =
            compute_second_order_perturbation_terms(world, XB, XC, AB_left, AB_right, BA_left, BA_right, phi0);
    second_order_x.truncate();

    auto beta_relaxed = inner(XA, second_order_x);

    if (world.rank() == 0) {
        cout << "beta unrelaxed " << endl;
        cout << beta_unrelaxed << endl;
    }

    if (world.rank() == 0) {
        cout << "beta relaxed " << endl;
        cout << beta_relaxed << endl;
    }

    auto beta = beta_relaxed + beta_unrelaxed;

    if (world.rank() == 0) {
        cout << "beta" << endl;
        cout << beta << endl;
    }

    return -2.0 * beta;
}

std::pair<X_space, X_space> QuadraticResponse::compute_first_order_fock_matrix_terms(World &world, const X_space &A,
                                                                                     const X_space &phi0,
                                                                                     const X_space &B) const {


    auto g1a = compute_g1_term(world, A, phi0, phi0);
    auto g1b = compute_g1_term(world, B, phi0, phi0);

    auto norm_g1ax = g1a.x.norm2();
    auto norm_g1bx = g1b.x.norm2();
    auto norm_g1ay = g1a.y.norm2();
    auto norm_g1by = g1b.y.norm2();

    if (world.rank() == 0) {
        cout << "norm_g1ax " << norm_g1ax << " norm_g1bx " << norm_g1bx << " norm_g1ay " << norm_g1ay << " norm_g1by "
             << norm_g1by << endl;
    }

    auto [VA, VB] = dipole_perturbation(world, phi0, phi0);
    auto normVAx = VA.x.norm2();
    auto normVAy = VA.y.norm2();
    auto normVBx = VB.x.norm2();
    auto normVBy = VB.y.norm2();
    if (world.rank() == 0) {
        cout << "normVAx " << normVAx << " normVAy " << normVAy << " normVBx " << normVBx << " normVBy " << normVBy
             << endl;
    }

    auto f1a = g1a + VA;
    auto f1b = g1b + VB;

    auto normf1ax = f1a.x.norm2();
    auto normf1bx = f1b.x.norm2();
    auto normf1ay = f1a.y.norm2();
    auto normf1by = f1b.y.norm2();
    if (world.rank() == 0) {
        cout << "normf1ax " << normf1ax << " normf1bx " << normf1bx << " normf1ay " << normf1ay << " normf1by "
             << normf1by << endl;
    }

    auto FAXB = X_space(world, A.num_states(), A.num_orbitals());
    auto FBXA = X_space(world, A.num_states(), A.num_orbitals());

    std::vector<std::string> names_1 = {"x", "x", "x", "y", "y", "z"};
    std::vector<std::string> names_2 = {"x", "y", "z", "y", "z", "z"};

    // Here contains the y components
    for (auto i = 0; i < A.num_states(); i++) {
        print("i", i);


        auto fax = matrix_inner(world, phi0.x[i], f1a.x[i]);
        auto fax_dagger = matrix_inner(world, phi0.y[i], f1a.y[i]);
        auto fbx = matrix_inner(world, phi0.x[i], f1b.x[i]);
        auto fb_dagger = matrix_inner(world, phi0.y[i], f1b.y[i]);


        FAXB.x[i] = copy(world, transform(world, B.x[i], fax, true), true);
        FAXB.y[i] = copy(world, transform(world, B.y[i], fax_dagger, true), true);
        FBXA.x[i] = copy(world, transform(world, A.x[i], fbx, true), true);
        FBXA.y[i] = copy(world, transform(world, A.y[i], fb_dagger, true), true);

        print("<0|F[", names_1[i], "]0>\n", fax);
        auto fxab_norm = norm2s_T(world, (FAXB.x[i]));
        print("norm of FAXB.x[i]", names_1[i], names_2[i], fxab_norm);

        print("<0|F_dagger[", names_1[i], "]0>\n", fax_dagger);
        auto fxab_dagger_norm = norm2s_T(world, (FAXB.y[i]));
        print("norm of FAXB.y[i]", names_1[i], names_2[i], fxab_dagger_norm);


        print("<0|F[", names_2[i], "]0>\n", fbx);
        auto fbx_norm = norm2s_T(world, (FBXA.x[i]));
        print("norm of FBXA.x[i]", names_2[i], names_1[i], fbx_norm);

        print("<0|F_dagger[", names_2[i], "]0>\n", fb_dagger);
        auto fbx_dagger_norm = norm2s_T(world, (FBXA.y[i]));
        print("norm of FBXA.x[i]", names_2[i], names_1[i], fbx_norm);
    }

    // print norms of FAXB and FBXA
    auto norm_FAXB_x = FAXB.x.norm2();
    auto norm_FAXB_y = FAXB.y.norm2();
    auto norm_FBXA_x = FBXA.x.norm2();
    auto norm_FBXA_y = FBXA.y.norm2();
    if (world.rank() == 0) {
        cout << "norm FAXB:\n x " << norm_FAXB_x << " y " << norm_FAXB_y << endl;
        cout << "norm FBXA:\n x " << norm_FBXA_x << " y " << norm_FBXA_y << endl;
    }
    FAXB.truncate();
    FBXA.truncate();


    world.gop.fence();

    return {FAXB, FBXA};
}

X_space QuadraticResponse::compute_second_order_perturbation_terms(World &world, const X_space &B, const X_space &C,
                                                                   const X_space &BC_left, const X_space &BC_right,
                                                                   const X_space &CB_left, const X_space &CB_right,
                                                                   const X_space &phi0) {
    // The first term to compute is -Q g1[K^BC], -Q g1[K^BC_conjugate]
    QProjector<double, 3> projector(world, ground_orbitals);
    auto apply_projector = [&](auto &xi) { return projector(xi); };

    // We have 2 terms to compute because we need to compute the contributions of BC and CB terms
    auto g1_kbc = -1.0 * oop_apply(compute_g1_term(world, BC_left, BC_right, phi0), apply_projector);
    auto g1_kcb = -1.0 * oop_apply(compute_g1_term(world, CB_left, CB_right, phi0), apply_projector);

    auto norm_g1_kbc_x = g1_kbc.x.norm2();
    auto norm_g1_kbc_y = g1_kbc.y.norm2();
    auto norm_g1_kcb_x = g1_kcb.x.norm2();
    auto norm_g1_kcb_y = g1_kcb.y.norm2();
    if (r_params.print_level() > 15) {

        if (world.rank() == 0) {
            cout << "norm_g1_kbc_x " << norm_g1_kbc_x << " norm_g1_kbc_y " << norm_g1_kbc_y << " norm_g1_kcb_x "
                 << norm_g1_kcb_x << " norm_g1_kcb_y " << norm_g1_kcb_y << endl;
        }
    }


    // the next term we need to compute are the first order fock matrix terms
    // -Q FB*XC
    auto g1bxc = -1.0 * oop_apply(compute_g1_term(world, B, phi0, C), apply_projector);
    // -Q FC*XB
    auto g1cxb = -1.0 * oop_apply(compute_g1_term(world, C, phi0, B), apply_projector);

    // -Q ( FB ) * ( XC )
    auto [vbxc, vcxb] = dipole_perturbation(world, C, B);

    vbxc = -1.0 * oop_apply(vbxc, apply_projector);
    vcxb = -1.0 * oop_apply(vcxb, apply_projector);

    auto [zFBzC, zFCzB] = compute_first_order_fock_matrix_terms(world, B, phi0, C);

    auto XA = -1.0 * x_data[0].first.copy();


    if (r_params.print_level() > 15) {

        // compute inner of each with XA
        auto g1_kbc_inner = inner(XA, g1_kbc);
        auto g1_kcb_inner = inner(XA, g1_kcb);
        auto g1cxb_inner = inner(XA, g1cxb);
        auto g1bxc_inner = inner(XA, g1bxc);
        auto zFBzC_inner = inner(XA, zFBzC);
        auto zFCzB_inner = inner(XA, zFCzB);
        auto VBXC_inner = inner(XA, vbxc);
        auto VCXB_inner = inner(XA, vcxb);
        // print the inner of each

        if (world.rank() == 0) {
            print("Inner products of XA with each term in rhs vectors\n");
            print("g1_kbc\n", g1_kbc_inner);
            print("g1_kcb\n", g1_kcb_inner);
            print("g1bxc\n", g1bxc_inner);
            print("g1cxb\n", g1cxb_inner);
            print("zFBzC\n", zFBzC_inner);
            print("zFCzB\n", zFCzB_inner);
            print("vbxc\n", VBXC_inner);
            print("vcxb\n", VCXB_inner);
        }
    }

    return g1_kbc + g1_kcb + g1bxc + g1cxb + zFBzC + zFCzB + vbxc + vcxb;

    // the next term we need to compute are the first order fock matrix terms
}


/**
 * @ brief computes the occ/occ and virt/virt blocks of second order gamma which leads to 4 space objects
 * ABX,ABY,BAX, and BAY where ABX contains A.X and phitildeab and ABY contains B.Y and phi0
 *
 *
 * @param world
 * @param A
 * @param B
 * @return
 */

std::tuple<X_space, X_space, X_space, X_space>
QuadraticResponse::compute_gamma_unrelaxed_block(World &world, const X_space &A, const X_space &B) {


    X_space AB_left = X_space(world, A.num_states(), B.num_orbitals());
    X_space AB_right = X_space(world, A.num_states(), B.num_orbitals());
    X_space BA_left = X_space(world, B.num_states(), A.num_orbitals());
    X_space BA_right = X_space(world, B.num_states(), A.num_orbitals());

    // Here are all the x components
    for (auto i = 0; i < A.num_states(); i++) {
        AB_left.x[i] = copy(world, A.x[i], false);
        AB_right.x[i] = copy(world, B.y[i], false);
        BA_left.x[i] = copy(world, B.x[i], false);
        BA_right.x[i] = copy(world, A.y[i], false);
    }

    std::vector<std::string> names = {"xx", "xy", "xz", "yy", "yz", "zz"};
    // Here contains the y components
    for (auto i = 0; i < A.num_states(); i++) {

        auto AB = matrix_inner(world, A.y[i], B.x[i]);
        auto BA = matrix_inner(world, B.y[i], A.x[i]);
        if (world.rank() == 0) { print(names[i], AB, BA); }
        // although we know there is a negative sign in the density matrix we do not include it here for _left.y terms
        auto temp = -1.0 * transform(world, ground_orbitals, 1.0 * AB, true);
        AB_left.y[i] = copy(world, temp, true);

        temp = -1.0 * transform(world, ground_orbitals, 1.0 * BA, true);
        BA_left.y[i] = copy(world, temp, true);

        auto norm_temp = norm2s_T(world, temp);
        auto norm_ground = norm2s_T(world, ground_orbitals);
        if (world.rank() == 0) { print("norm of temp: ", norm_temp); }
        if (world.rank() == 0) { print("norm of ground ", norm_ground); }

        AB_right.y[i] = copy(world, ground_orbitals, true);
        BA_right.y[i] = copy(world, ground_orbitals, true);


        //print the norms of each vector
    }


    if (r_params.print_level() > 0) {
        auto AB_left_x = AB_left.x.norm2();
        auto AB_right_x = AB_right.x.norm2();
        auto AB_left_y = AB_left.y.norm2();
        auto AB_right_y = AB_left.y.norm2();

        auto BA_left_x = BA_left.x.norm2();
        auto BA_right_x = BA_right.x.norm2();
        auto BA_left_y = BA_left.y.norm2();
        auto BA_right_y = BA_left.y.norm2();

        if (world.rank() == 0) {
            print("NORMS of UNRELAXED BLOCKS");
            print("AB_left_x");
            print(AB_left_x);
            print("AB_right_x");
            print(AB_right_x);
            print("BA_left_x");
            print(BA_left_x);
            print("BA_right_x");
            print(BA_right_x);
            print("AB_left_y");
            print(AB_left_y);
            print("AB_right_y");
            print(AB_right_y);
            print("BA_left_y");
            print(BA_left_y);
            print("BA_right_y");
            print(BA_right_y);
        }
    }

    AB_left.truncate();
    AB_right.truncate();
    BA_left.truncate();
    BA_right.truncate();


    return {AB_left, AB_right, BA_left, BA_right};
}


auto QuadraticResponse::dipole_perturbation(World &world, const X_space &left, const X_space &right) const
        -> std::pair<X_space, X_space> {
    vector_real_function_3d dipole_vectors(3);
    size_t i = 0;
    // creates a vector of x y z dipole functions
    for (auto &d: dipole_vectors) {
        std::vector<int> f(3, 0);
        f[i++] = 1;
        d = real_factory_3d(world).functor(real_functor_3d(new MomentFunctor(f)));
    }
    truncate(world, dipole_vectors, true);

    std::array<int, 6> case_1_indices = {0, 0, 0, 1, 1, 2};
    std::array<int, 6> case_2_indices = {0, 1, 2, 1, 2, 2};

    auto VB = X_space(world, left.num_states(), right.num_orbitals());

    for (int i = 0; i < 6; i++) {
        VB.x[i] = dipole_vectors[case_1_indices[i]] * left.x[i];
        VB.y[i] = dipole_vectors[case_1_indices[i]] * left.y[i];
    }

    auto VC = X_space(world, left.num_states(), right.num_orbitals());

    for (int i = 0; i < 6; i++) {
        VC.x[i] = dipole_vectors[case_2_indices[i]] * right.x[i];
        VC.y[i] = dipole_vectors[case_2_indices[i]] * right.y[i];
    }

    VB.truncate();
    VC.truncate();

    return {1.0 * VB, 1.0 * VC};
    // in the first case i need to mulitply x x x y y z to vectors in right
    // in the second case i need to multiply x y z y z a to vectors in left
}


// compute the g1 term
// we have 3 x_space which hold x and y components of the density matrix.
//
X_space QuadraticResponse::compute_g1_term(World &world, const X_space &left, const X_space &right,
                                           const X_space &apply) const {


    auto JBC = compute_coulomb_term(world, left, right, apply);
    auto KBC = compute_exchange_term(world, left, right, apply);

    return 2 * JBC - KBC;
}
// compute rhoBC and apply to D
X_space QuadraticResponse::compute_coulomb_term(World &world, const X_space &x_left, const X_space &x_right,
                                                const X_space &x_apply) const {

    X_space J = X_space::zero_functions(world, x_left.num_states(), x_left.num_orbitals());
    vector_real_function_3d rhoX(x_left.num_states());
    vector_real_function_3d temp_J(x_left.num_states());

    vector_real_function_3d x_phi, y_phi;

    // create the density for each state in B
    for (const auto &b: x_left.active) {
        x_phi = mul(world, x_left.x[b], x_right.x[b], false);
        y_phi = mul(world, x_left.y[b], x_right.y[b], false);
        world.gop.fence();
        rhoX[b] = sum(world, x_phi, true);
        rhoX[b] += sum(world, y_phi, true);
        world.gop.fence();
    }

    truncate(world, rhoX);

    for (const auto &k: x_left.active) {
        temp_J[k] = apply(*shared_coulomb_operator, rhoX[k]);
        J.x[k] = mul(world, temp_J[k], x_apply.x[k], false);
        J.y[k] = mul(world, temp_J[k], x_apply.y[k], false);
    }
    world.gop.fence();
    J.truncate();

    return J;
}

auto Koperator(const vecfuncT &ket, const vecfuncT &bra) {
    const double lo = 1.e-10;
    auto &world = ket[0].world();
    Exchange<double, 3> k{world, lo};
    k.set_bra_and_ket(bra, ket);
    k.set_algorithm(k.multiworld_efficient);
    return k;
};

X_space QuadraticResponse::compute_exchange_term(World &world, const X_space &x_left, const X_space &x_right,
                                                 const X_space &x_apply) const {


    // if the frequecy of x_left is 0 we run the static case
    // else we run the dynamic case
    auto K = X_space::zero_functions(world, x_left.num_states(), x_left.num_orbitals());

    vector_real_function_3d xb;
    vector_real_function_3d yb;

    auto k1 = create_response_matrix(x_left.num_states(), x_left.num_orbitals());
    auto k2 = create_response_matrix(x_left.num_states(), x_left.num_orbitals());
    auto k1_conjugate = create_response_matrix(x_left.num_states(), x_left.num_orbitals());
    auto k2_conjugate = create_response_matrix(x_left.num_states(), x_left.num_orbitals());


    for (int k = 0; k < x_left.num_states(); k++) {

        auto K1 = Koperator(x_left.x[k], x_right.x[k]);
        auto K2 = Koperator(x_right.y[k], x_left.y[k]);

        auto K1_conjugate = Koperator(x_left.y[k], x_right.y[k]);
        auto K2_conjugate = Koperator(x_right.x[k], x_left.x[k]);
        world.gop.fence();

        k1[k] = K1(x_apply.x[k]);
        k2[k] = K2(x_apply.x[k]);


        k1_conjugate[k] = K1_conjugate(x_apply.y[k]);
        k2_conjugate[k] = K2_conjugate(x_apply.y[k]);
        world.gop.fence();
        K.x[k] = gaxpy_oop(1.0, k1[k], 1.0, k2[k], false);
        K.y[k] = gaxpy_oop(1.0, k1_conjugate[k], 1.0, k2_conjugate[k], false);
    }
    world.gop.fence();
    K.truncate();

    return K;
}


//
