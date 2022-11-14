#include "global_functions.h"

#include <ResponseBase.hpp>
#include <memory>

#include "madness/chem/SCFOperators.h"
#include "response_parameters.h"


auto initialize_calc_params(World &world, const std::string &input_file) -> CalcParams {
    ResponseParameters r_params{};
    commandlineparser parser;
    parser.set_keyval("input", input_file);
    r_params.read_input_and_commandline_options(world, parser, "response");
    GroundStateCalculation ground_calculation{world, r_params.archive()};
    if (world.rank() == 0) {
        ground_calculation.print_params();
    }
    Molecule molecule = ground_calculation.molecule();
    r_params.set_ground_state_calculation_data(ground_calculation);
    r_params.set_derived_values(world, molecule);
    if (world.rank() == 0) {
        r_params.print();
    }
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
    X_space K0 = X_space(world, x.num_states(), x.num_orbitals());
    long n{};
    response_matrix xx;
    // place all x and y functions into a single response vector
    if (compute_y) {
        n = 2;
        xx = to_response_matrix(x);
        // place all x
    } else {
        n = 1;
        xx = x.X.x;
        // if not compute y we are only working with the x functions
    }
    // should have num_states * num_orbitals * n  if compute y  n=2 else n=1
    vecfuncT x_vector(x.num_states() * n * x.num_orbitals());
    long ij = 0;
    for (const auto &xi: xx) {// copy the response matrix into a single vector of functions
        std::for_each(xi.begin(), xi.end(), [&](const auto &xij) {
            x_vector.at(ij++) = copy(xij);
        });
    }
    vecfuncT phi_vector(x.num_states() * n * phi0.size());
    int orb_i = 0;
    // copy ground-state orbitals into a single long vector
    std::for_each(phi_vector.begin(), phi_vector.end(), [&](auto &phi_i) { phi_i = copy(phi0[orb_i++ % x.num_orbitals()]); });
    world.gop.fence();
    molresponse::end_timer(world, "ground exchange copy");
    molresponse::start_timer(world);
    const double lo = 1.e-10;
    Exchange<double, 3> op{};
    // Do exchange by creating operator with parameters
    op.set_parameters(phi_vector, madness::copy(world, phi_vector), lo);
    op.set_algorithm(op.multiworld_efficient);
    world.gop.fence();
    // apply exchange phi phi x
    auto exchange_vector = op(x_vector);
    molresponse::end_timer(world, "ground exchange apply");
    molresponse::start_timer(world);

    auto exchange_matrix = create_response_matrix(x.num_states(), n * x.num_orbitals());
    long b = 0;
    for (auto &xi: exchange_matrix) {
        for (auto &xij: xi) {
            xij = copy(exchange_vector[b++]);
        }
    }
    if (compute_y) {
        K0 = to_X_space(exchange_matrix);
    } else {
        K0.X = exchange_matrix;
        K0.Y = K0.X.copy();
    }

    world.gop.fence();
    molresponse::end_timer(world, "ground exchange reorganize");
    return K0;
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
    X_space K = X_space(world, x.num_states(), x.num_orbitals());
    X_space conjugateK = X_space(world, x.num_states(), x.num_orbitals());
    auto num_orbitals = phi0.size();
    long n{};
    response_matrix xx;
    response_matrix xx_dagger;
    if (compute_y) {
        n = 2;
        xx = to_response_matrix(x);
        xx_dagger = to_conjugate_response_matrix(x);
        // place all x
    } else {
        n = 1;
        xx = x.X.x;
        xx_dagger = x.X.x;
        // if not compute y we are only working with the x functions
    }
    vecfuncT x_vector(x.num_states() * n * x.num_orbitals());
    long ij = 0;
    for (const auto &xi: xx) {// copy the response matrix into a single vector of functions
        std::for_each(xi.begin(), xi.end(), [&](const auto &xij) {
            x_vector.at(ij++) = copy(xij);
        });
    }
    ij = 0;
    vecfuncT x_dagger_vector(x.num_states() * n * x.num_orbitals());
    for (const auto &xdi: xx_dagger) {// copy the response matrix into a single vector of functions
        std::for_each(xdi.begin(), xdi.end(), [&](const auto &xdij) {
            x_dagger_vector.at(ij++) = copy(xdij);
        });
    }
    vecfuncT phi_vector(x.num_states() * n * phi0.size());
    int orb_i = 0;
    // copy ground-state orbitals into a single long vector
    std::for_each(phi_vector.begin(), phi_vector.end(), [&](auto &phi_i) { phi_i = copy(phi0[orb_i++ % x.num_orbitals()]); });
    world.gop.fence();
    molresponse::end_timer(world, "response exchange copy");
    molresponse::start_timer(world);
    const double lo = 1.e-10;
    auto phi_copy1 = madness::copy(world, phi_vector, true);
    auto phi_copy2 = madness::copy(world, phi_vector, true);
    // We have 2 versions of exchange  k(x,phi) phi and k(phi,x) phi
    // K1.X = k[x ,phi] phi , K2.X = k[phi,y] phi // if static we only compute the top line
    // K1.Y = k[y',phi] phi , K2.Y = k[phi,x] phi
    // K=K1+K2
    Exchange<double, 3> K1{};
    K1.set_parameters(x_vector, phi_copy1, lo);
    K1.set_algorithm(K1.multiworld_efficient);
    Exchange<double, 3> K2{};
    K2.set_parameters(phi_copy2, x_dagger_vector, lo);
    K2.set_algorithm(K2.multiworld_efficient);
    world.gop.fence();

    auto K1_vect = K1(phi_vector);
    world.gop.fence();
    auto K2_vect = K2(phi_vector);
    world.gop.fence();
    molresponse::end_timer(world, "response exchange apply");
    molresponse::start_timer(world);
    vecfuncT K_vector = K1_vect + K2_vect;
    world.gop.fence();

    auto exchange_matrix = create_response_matrix(x.num_states(), n * x.num_orbitals());
    long b = 0;
    for (auto &xi: exchange_matrix) {
        for (auto &xij: xi) {
            xij = copy(K_vector[b++]);
        }
    }
    world.gop.fence();
    if (compute_y) {
        K = to_X_space(exchange_matrix);
    } else {
        K.X = exchange_matrix;
        K.Y = K.X.copy();
    }
    world.gop.fence();
    molresponse::end_timer(world, "response exchange reorganize");
    return K;
}
// compute exchange |i><i|J|p>
auto newK(const vecfuncT &ket, const vecfuncT &bra, const vecfuncT &vf) -> vecfuncT {
    World &world = ket[0].world();
    const double lo = 1.e-10;

    Exchange<double, 3> op{};
    op.set_parameters(bra, ket, lo);
    op.set_algorithm(op.multiworld_efficient);
    auto vk = op(vf);
    return vk;
}
// sum_i |i><i|J|p> for each p