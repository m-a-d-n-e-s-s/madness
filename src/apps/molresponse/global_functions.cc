#include "global_functions.h"

#include <ResponseBase.hpp>
#include <memory>

#include "madness/chem/SCFOperators.h"
#include "response_parameters.h"


static auto set_poisson(World &world, const double lo, const double econv = FunctionDefaults<3>::get_thresh()) {
    return std::shared_ptr<real_convolution_3d>(CoulombOperatorPtr(world, lo, econv));
}


auto initialize_calc_params(World &world, const std::string &input_file) -> CalcParams {
    ResponseParameters r_params{};
    commandlineparser parser;
    parser.set_keyval("input", input_file);
    r_params.read_input_and_commandline_options(world, parser, "response");
    GroundStateCalculation ground_calculation{world, r_params.archive()};
    if (world.rank() == 0) { ground_calculation.print_params(); }
    Molecule molecule = ground_calculation.molecule();
    r_params.set_ground_state_calculation_data(ground_calculation);
    r_params.set_derived_values(world, molecule);
    if (world.rank() == 0) { r_params.print(); }
    world.gop.fence();
    return {ground_calculation, molecule, r_params};
}
// TODO some operator definitions that I will need to move to a separate file
auto T(World &world, response_space &f) -> response_space {
    response_space T;// Fock = (T + V) * orbitals
    real_derivative_3d Dx(world, 0);
    real_derivative_3d Dy(world, 1);
    real_derivative_3d Dz(world, 2);
    // Apply derivatives to orbitals
    f.reconstruct_rf();
    response_space dvx = apply(world, Dx, f);
    response_space dvy = apply(world, Dy, f);
    response_space dvz = apply(world, Dz, f);
    // Apply again for 2nd derivatives
    response_space dvx2 = apply(world, Dx, dvx);
    response_space dvy2 = apply(world, Dy, dvy);
    response_space dvz2 = apply(world, Dz, dvz);
    T = (dvx2 + dvy2 + dvz2) * (-0.5);
    return T;
}

// compute exchange |i><i|J|p>
/**
 * Computes ground density exchange on response vectors
 *  This algorithm places all functions in a single vector,
 *  computes exchange and returns the result into a response
 *  matrix
 *
 * @param v1
 * @param v2
 * @param f
 * @return
 */
auto ground_exchange(const vecfuncT &phi0, const X_space &x, const bool compute_y) -> X_space {
    World &world = phi0[0].world();
    molresponse::start_timer(world);
    auto num_states = x.num_states();
    auto num_orbitals = x.num_orbitals();
    X_space K0 = X_space::zero_functions(world, x.num_states(), x.num_orbitals());
    long n{};
    response_matrix xx;
    // place all x and y functions into a single response vector
    if (compute_y) {
        n = 2;
        xx = to_response_matrix(x);
        // place all x
    } else {
        n = 1;
        xx = x.x.x;
        // if not compute y we are only working with the x functions
    }
    auto n_exchange = n * num_states * num_orbitals * num_orbitals;
    // should have num_states * num_orbitals * n  if compute y  n=2 else n=1
    vecfuncT x_vector(n_exchange);
    long b_index = 0;
    long p_index = 0;
    long p = 0;
    for (long b = 0; b < num_states; b++) {
        b_index = b * num_orbitals * num_orbitals * n;
        p = 0;
        // for each function in a response vector copy num orbital times
        std::for_each(xx[b].begin(), xx[b].end(), [&](const auto &xb_p) {
            p_index = p * num_orbitals;
            for (long j = 0; j < num_orbitals; j++) { x_vector.at(b_index + p_index + j) = copy(xb_p, false); }
            p++;
        });
    }
    vecfuncT phi1(n_exchange);
    vecfuncT phi2(n_exchange);
    int orb_i = 0;
    // copy ground-state orbitals into a single long vector
    std::for_each(phi1.begin(), phi1.end(), [&](auto &phi_i) { phi_i = copy(phi0[orb_i++ % num_orbitals]); });
    world.gop.fence();
    phi2 = madness::copy(world, phi1);
    world.gop.fence();
    molresponse::end_timer(world, "ground exchange copy");
    auto K = molresponseExchange(world, phi1, phi2, x_vector, n, num_states, num_orbitals);
    if (world.rank() == 0) { print("made it out of molresponseExchange"); }
    return K;
}
// compute full response exchange |i><i|J|p>
/**
 * Computes ground density exchange on response vectors
 *  This algorithm places all functions in a single vector,
 *  computes exchange and returns the result into a response
 *  matrix
 *
 * @param v1
 * @param v2
 * @param f
 * @return
 */
auto response_exchange(const vecfuncT &phi0, const X_space &x, const bool compute_y) -> X_space {
    World &world = phi0[0].world();
    molresponse::start_timer(world);
    auto num_states = x.num_states();
    auto num_orbitals = x.num_orbitals();
    X_space K1 = X_space(world, num_states, num_orbitals);
    X_space K2 = X_space(world, num_states, num_orbitals);
    X_space K = X_space(world, num_states, num_orbitals);

    long n = compute_y ? 2 : 1;
    auto n_exchange = n * num_states * num_orbitals * num_orbitals;
    vecfuncT phi_right(n_exchange);
    long b_index = 0;
    long p_index = 0;
    long p = 0;
    // 111 222 333 111 222 333 111 222 333 111 222 333
    for (long b = 0; b < num_states; b++) {
        b_index = n * b * num_orbitals * num_orbitals;
        p = 0;
        std::for_each(phi0.begin(), phi0.end(), [&](const auto &phi_p) {
            p_index = num_orbitals * p;
            for (long j = 0; j < num_orbitals; j++) { phi_right.at(b_index + p_index + j) = copy(phi_p, false); }
            p++;
        });
        if (compute_y) {
            p = 0;
            std::for_each(phi0.begin(), phi0.end(), [&](const auto &phi_p) {
                p_index = num_orbitals * p;
                for (long j = 0; j < num_orbitals; j++) { phi_right.at(num_orbitals * num_orbitals + b_index + p_index + j) = copy(phi_p, false); }
                p++;
            });
        }
    }
    world.gop.fence();

    vecfuncT x_vector(n_exchange);
    vecfuncT x_vector_conjugate(n_exchange);
    for (long b = 0; b < num_states; b++) {// for each state copy the vector n times
        b_index = b * num_orbitals * num_orbitals * n;
        for (long j = 0; j < num_orbitals; j++) {
            p_index = j * num_orbitals;
            std::transform(x.x[b].begin(), x.x[b].end(), x_vector.begin() + b_index + p_index, [&](const auto &xbi) { return copy(xbi, false); });

            std::transform(x.y[b].begin(), x.y[b].end(), x_vector_conjugate.begin() + b_index + p_index, [&](const auto &xbi) { return copy(xbi, false); });
        }
        if (compute_y) {
            long y_shift = num_orbitals * num_orbitals;
            for (long j = 0; j < num_orbitals; j++) {
                p_index = j * num_orbitals;
                std::transform(x.y[b].begin(), x.y[b].end(), x_vector.begin() + b_index + p_index + y_shift, [&](const auto &ybi) { return copy(ybi, false); });
                std::transform(x.x[b].begin(), x.x[b].end(), x_vector_conjugate.begin() + b_index + p_index + y_shift,
                               [&](const auto &ybi) { return copy(ybi, false); });
            }
        }
    }
    world.gop.fence();

    vecfuncT phi_left(n_exchange);
    for (long b = 0; b < num_states; b++) {
        b_index = b * num_orbitals * num_orbitals * n;
        for (long p = 0; p < n * num_orbitals; p++) {
            p_index = p * num_orbitals;
            for (int i = 0; i < num_orbitals; i++) { phi_left[b_index + p_index + i] = copy(phi0[i], false); }
        }
    }
    world.gop.fence();
    molresponse::end_timer(world, "response exchange copy");

    molresponse::start_timer(world);
    K1 = molresponseExchange(world, x_vector, phi_left, phi_right, n, num_states, num_orbitals);
    world.gop.fence();
    K2 = molresponseExchange(world, phi_left, x_vector_conjugate, phi_right, n, num_states, num_orbitals);
    world.gop.fence();
    /*
    auto xk1 = inner(x, K1);
    if (world.rank() == 0) { print("new xk1\n", xk1); }
    auto xk2 = inner(x, K2);
    if (world.rank() == 0) { print("new xk2\n", xk2); }
     */
    K = K1 + K2;
    world.gop.fence();
    molresponse::end_timer(world, "response exchange K1+K2 ");
    return K;
}
// compute exchange |i><i|J|p>
auto newK(const vecfuncT &ket, const vecfuncT &bra, const vecfuncT &vf) -> vecfuncT {
    World &world = ket[0].world();
    const double lo = 1.e-10;

    Exchange<double, 3> op{world, lo};
    op.set_bra_and_ket(bra, ket);
    op.set_algorithm(op.multiworld_efficient);
    auto vk = op(vf);
    return vk;
}
// sum_i |i><i|J|p> for each p
auto molresponseExchange(World &world, const vecfuncT &ket_i, const vecfuncT &bra_i, const vecfuncT &fp, const int &n, const int &num_states,
                         const int &num_orbitals) -> X_space {
    molresponse::start_timer(world);
    reconstruct(world, ket_i, false);
    reconstruct(world, bra_i, false);
    reconstruct(world, fp, false);
    if (world.rank() == 0) { print("exchange reconstruct tree"); }
    world.gop.fence();
    /*
    norm_tree(world, ket_i, false);
    norm_tree(world, bra_i, false);
    norm_tree(world, fp, false);
     */
    if (world.rank() == 0) { print("exchange norm tree"); }
    const double lo = 1.0e-10;
    auto tol = FunctionDefaults<3>::get_thresh();
    auto poisson = set_poisson(world, lo);
    world.gop.fence();
    if (world.rank() == 0) { print("create poisson v23"); }
    auto v23 = mul(world, bra_i, fp, true);
    //mul(world, bra_i, fp, true);
    if (world.rank() == 0) { print("multiply v23"); }
    truncate(world, v23, tol, true);
    // truncate
    if (world.rank() == 0) { print("truncate v23"); }
    v23 = apply(world, *poisson, v23);
    if (world.rank() == 0) { print("apply v23"); }
    truncate(world, v23, tol, true);
    if (world.rank() == 0) { print("truncate after apply v23"); }
    auto v123 = mul(world, ket_i, v23, true);
    if (world.rank() == 0) { print("multiply  apply v123"); }
    const long n_exchange{num_states * n * num_orbitals};
    auto exchange_vector = vecfuncT(n_exchange);
    long b = 0;
    long b_shift = 0;
    compress(world, v123);
    for (auto &kij: exchange_vector) {
        b_shift = b * num_orbitals;
        auto phi_phiX_i = vecfuncT(num_orbitals);
        // this right here is a sketch of how sums with iterators could work
        kij = FunctionFactory<double, 3>(world).compressed();
        std::for_each(v123.begin() + b_shift, v123.begin() + b_shift + num_orbitals, [&](const auto &v123_i) { kij.gaxpy(1.0, v123_i, 1.0, false); });
        b++;
    }
    world.gop.fence();
    if (world.rank() == 0) print("exchange sum");
    truncate(world, exchange_vector, tol, true);
    molresponse::end_timer(world, "exchange apply");

    molresponse::start_timer(world);
    auto exchange_matrix = create_response_matrix(num_states, n * num_orbitals);
    b = 0;
    for (auto &xi: exchange_matrix) {
        for (auto &xij: xi) { xij = copy(exchange_vector[b++]); }
    }
    X_space K0 = X_space::zero_functions(world, num_states, num_orbitals);
    if (n == 2) {
        K0 = to_X_space(exchange_matrix);
    } else {
        K0.x = exchange_matrix;
    }
    world.gop.fence();
    molresponse::end_timer(world, "ground exchange reorganize");
    return K0;
}
auto make_k(const vecfuncT &ket, const vecfuncT &bra, const ResponseParameters &r_params) {
    const double lo = 1.e-10;
    auto &world = ket[0].world();
    Exchange<double, 3> k{world, lo};
    k.set_bra_and_ket(bra, ket);
    if (r_params.hfexalg() == "multiworld") {
        k.set_algorithm(madness::Exchange<double, 3>::multiworld_efficient);
    } else if (r_params.hfexalg() == "smallmem") {
        k.set_algorithm(madness::Exchange<double, 3>::small_memory);
    } else if (r_params.hfexalg() == "largemem") {
        k.set_algorithm(madness::Exchange<double, 3>::large_memory);
    };
    return k;
}
/**
 * Computes ground density exchange on response vectors
 *  This algorithm places all functions in a single vector,
 *  computes exchange and returns the result into a response
 *  matrix
 *
 * @param v1
 * @param v2
 * @param f
 * @return
 */
auto response_exchange(const vecfuncT &phi0, const X_space &chi, const bool &compute_y, const ResponseParameters &r_params) -> X_space {
    World &world = phi0[0].world();
    molresponse::start_timer(world);
    auto num_states = chi.num_states();
    auto num_orbitals = chi.num_orbitals();
    auto K = X_space::zero_functions(world, num_states, num_orbitals);
    vector_real_function_3d k1x, k1y, k2x, k2y;


    if (compute_y) {
        for (const auto &b: chi.active) {
            auto x = chi.x[b];
            auto y = chi.y[b];
            auto K1X = make_k(x, phi0, r_params);
            auto K2X = make_k(y, phi0, r_params);
            auto K1Y = make_k(phi0, y, r_params);
            auto K2Y = make_k(phi0, x, r_params);
            k1x = K1X(phi0);
            k1y = K1Y(phi0);
            k2x = K2X(phi0);
            k2y = K2Y(phi0);
            world.gop.fence();
            K.x[b] = gaxpy_oop(1.0, k1x, 1.0, k1y, false);
            K.y[b] = gaxpy_oop(1.0, k2x, 1.0, k2y, false);
            world.gop.fence();
        }
    } else {
        for (const auto &b: chi.active) {
            auto x = chi.x[b];
            auto y = chi.x[b];
            auto K1X = make_k(x, phi0, r_params);
            auto K1Y = make_k(phi0, y, r_params);
            k1x = K1X(phi0);
            k1y = K1Y(phi0);
            K.x[b] = gaxpy_oop(1.0, k1x, 1.0, k1y, true);
        }
    }
    return K;
}

// sum_i |i><i|J|p> for each p
/**
 * Computes ground density exchange on response vectors
 *  This algorithm places all functions in a single vector,
 *  computes exchange and returns the result into a response
 *  matrix
 *
 * @param v1
 * @param v2
 * @param f
 * @return
 */
auto ground_exchange(const vecfuncT &phi0, const X_space &chi, const bool &compute_y, const ResponseParameters &r_params) -> X_space {
    World &world = phi0[0].world();
    molresponse::start_timer(world);

    auto num_states = chi.num_states();
    auto num_orbitals = chi.num_orbitals();

    auto K0 = X_space::zero_functions(world, num_states, num_orbitals);
    K0.set_active(chi.active);
    // the question is copying pointers mpi safe
    world.gop.fence();
    auto k0 = make_k(phi0, phi0, r_params);
    if (compute_y) {
        for (const auto &i: chi.active) {
            K0.x[i] = k0(chi.x[i]);
            K0.y[i] = k0(chi.y[i]);
        }
    } else {
        for (const auto &i: chi.active) { K0.x[i] = k0(chi.x[i]); }
        K0.y = K0.x.copy();
    }
    K0.truncate();
    return K0;
}
;
