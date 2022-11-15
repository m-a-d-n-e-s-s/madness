#include "global_functions.h"

#include <ResponseBase.hpp>
#include <memory>

#include "madness/chem/SCFOperators.h"
#include "response_parameters.h"


static auto set_poisson(World &world, const double lo,
                        const double econv = FunctionDefaults<3>::get_thresh()) {
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
    vecfuncT x_vector(x.num_states() * n * x.num_orbitals() * x.num_orbitals());
    long ij = 0;


    size_t b = 0;
    for (const auto &xi: xx) {// copy the response matrix into a single vector of functions
        std::for_each(xi.begin(), xi.end(), [&](const auto &xij) {
            auto vect_xij = copy(world, xij, x.num_orbitals(), true);
            std::copy(vect_xij.begin(), vect_xij.end(),
                      x_vector.begin() +
                              b * x.num_orbitals());// copy the xij vector n times into a block
            b++;
        });
    }
    vecfuncT phi1(x.num_orbitals() * x.num_orbitals() * n * x.num_states());
    vecfuncT phi2(x.num_orbitals() * x.num_orbitals() * n * x.num_states());
    int orb_i = 0;
    // copy ground-state orbitals into a single long vector
    std::for_each(phi1.begin(), phi1.end(),
                  [&](auto &phi_i) { phi_i = copy(phi0[orb_i++ % x.num_orbitals()]); });
    phi2 = madness::copy(world, phi1);
    world.gop.fence();
    molresponse::end_timer(world, "ground exchange copy");
    return molresponse_exchange(world, phi1, phi2, x_vector, n, x.num_states(), x.num_orbitals());
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
    auto n_exchange = n * num_states * num_orbitals * num_orbitals;
    vecfuncT x_vector(n_exchange);
    long ij = 0;
    long x_index = 0;
    for (const auto &xi: xx) {                      // num_states
        for (const auto &xij: xi) {                 //num_orbitals*n
            for (int z = 0; z < num_orbitals; z++) {//*num_orbitals
                x_index = z + num_orbitals * ij;
                x_vector.at(x_index) = copy(xij, false);
            }
            world.gop.fence();
            ij++;
        }
    };
    ij = 0;
    x_index = 0;
    vecfuncT x_vector_conjugate(n_exchange);
    for (const auto &xdi: xx_dagger) {// copy the response matrix into a single vector of functions
        for (const auto &xdij: xdi) {
            for (int z = 0; z < num_orbitals; z++) {
                x_index = z + num_orbitals * ij;
                x_vector_conjugate.at(x_index) = copy(xdij, false);
            }
            world.gop.fence();
            ij++;
        }
    }
    vecfuncT phi_right(n_exchange);
    // copy ground-state orbitals into a single long vector
    long ii = 0;
    long z_index = 0;
    long bshift = 0;
    auto orb_shift = 0;
    std::for_each(phi0.begin(), phi0.end(), [&](const auto &phi_i) {
        orb_shift = n * num_orbitals * ii;
        for (long b = 0; b < num_states; b++) {
            bshift = b * n * num_orbitals*num_orbitals;
            for (long z = 0; z < n * num_orbitals; z++) {
                z_index = bshift + orb_shift + z;
                phi_right.at(z_index) = copy(phi_i);
            }
        }
        ii++;
        world.gop.fence();
    });
    world.gop.fence();
    vecfuncT phi_left(n_exchange);
    long orb_i = 0;
    for (auto &phi_i: phi_left) { phi_i = copy(phi0[orb_i++ % num_orbitals]); }
    world.gop.fence();
    molresponse::end_timer(world, "response exchange copy");

    molresponse::start_timer(world);
    K1 = molresponse_exchange(world, x_vector, phi_left, phi_right, n, num_states, num_orbitals);
    K2 = molresponse_exchange(world, phi_left, x_vector_conjugate, phi_right, n, num_states,
                              num_orbitals);
    K = K1 + K2;
    world.gop.fence();
    molresponse::end_timer(world, "response exchange K1+K2 ");
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
auto molresponse_exchange(World &world, const vecfuncT &v1, const vecfuncT &v2, const vecfuncT &v3,
                          const int &n, const int &num_states, const int &num_orbitals) -> X_space {
    molresponse::start_timer(world);

    reconstruct(world, v1, false);
    reconstruct(world, v2, false);
    reconstruct(world, v3, false);
    world.gop.fence();
    norm_tree(world, v1, false);
    norm_tree(world, v2, false);
    norm_tree(world, v3, false);
    world.gop.fence();
    const double lo = 1.0e-10;
    auto poisson = set_poisson(world, lo);
    auto v23 = mul(world, v2, v3, true);
    truncate(world, v23);
    v23 = apply(world, *poisson, v23);
    truncate(world, v23);
    auto v123 = mul(world, v1, v23, true);
    const long n_exchange{num_states * n * num_orbitals};
    auto exchange_vector = vecfuncT(n_exchange);

    long b = 0;
    for (auto &kij: exchange_vector) {
        auto phi_phiX_i = vecfuncT(num_orbitals);
        std::copy(v123.begin() + (b * num_orbitals),
                  v123.begin() + (b * num_orbitals) + num_orbitals, phi_phiX_i.begin());
        world.gop.fence();
        kij = sum(world, phi_phiX_i, true);
        b++;
        // option to use inner product kij=std::inner_product(phi_phiX.begin()+(b*x.num_orbitals()),phi_phiX.begin()+(b*x.num_orbitals(),)
    }
    molresponse::end_timer(world, "ground exchange apply");
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
        K0.X = exchange_matrix;
    }
    world.gop.fence();
    molresponse::end_timer(world, "ground exchange reorganize");
    return K0;
}
// sum_i |i><i|J|p> for each p
