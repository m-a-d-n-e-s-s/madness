#include "global_functions.h"

#include <ResponseBase.hpp>
#include <memory>
#include <string>
#include <vector>

#include "madness/chem/SCFOperators.h"
#include "response_parameters.h"

void print_molecule(World &world, const GroundStateCalculation &g_params) {
    if (world.rank() == 0) {
        // Precision is set to 10 coming in, drop it to 5
        std::cout.precision(5);
        std::cout << std::fixed;

        // First get atom
        const std::vector<Atom> atoms = g_params.molecule().get_atoms();
        size_t num_atoms = atoms.size();

        // Now print
        print("\n   Geometry Information");
        print("   --------------------\n");
        print("   Units: a.u.\n");
        print(" Atom            x                 y                 z");
        print("----------------------------------------------------------------");
        for (size_t j = 0; j < num_atoms; j++) {
            Vector<double, 3> coords = atoms[j].get_coords();
            std::cout << std::setw(3) << atomic_number_to_symbol(atoms[j].get_atomic_number());
            std::cout << std::setw(18) << std::right << coords[0] << std::setw(18) << coords[1]
                      << std::setw(18) << coords[2] << endl;
        }
        print("");

        // Reset precision
        std::cout.precision(10);
        std::cout << std::scientific;
    }
}

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
        K0.X = K0.Y.copy();
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
auto response_exchange(const vecfuncT &phi0, const response_matrix &x, const response_matrix &x_dagger, const bool static_response) -> response_matrix {
    World &world = phi0[0].world();
    molresponse::start_timer(world);
    auto num_orbitals = phi0.size();
    long n{};
    n = (static_response) ? 1 : 2;
    vecfuncT phi_vect(x.size() * n * phi0.size());
    vecfuncT x_vect(x.size() * n * phi0.size());
    vecfuncT xd_vect(x.size() * n * phi0.size());
    int orb_i = 0;
    for (auto &phi_i: phi_vect) { phi_i = copy(phi0[orb_i++ % num_orbitals]); }
    long j = 0;
    for (const auto &xi: x) {
        for (const auto &xij: xi) {
            x_vect[j++] = copy(xij);
        }
    }
    world.gop.fence();
    j = 0;
    for (const auto &xi: x_dagger) {
        for (const auto &xij: xi) {
            xd_vect[j++] = copy(xij);
        }
    }
    world.gop.fence();
    molresponse::end_timer(world, "response exchange copy");
    molresponse::start_timer(world);
    const double lo = 1.e-10;
    auto phi_copy1 = madness::copy(world, phi_vect, true);
    auto phi_copy2 = madness::copy(world, phi_vect, true);
    Exchange<double, 3> x_phi_K{};
    x_phi_K.set_parameters(x_vect, phi_copy1, lo);
    x_phi_K.set_algorithm(x_phi_K.multiworld_efficient);
    Exchange<double, 3> phi_xd_K{};
    phi_xd_K.set_parameters(phi_copy2, xd_vect, lo);
    phi_xd_K.set_algorithm(phi_xd_K.multiworld_efficient);

    world.gop.fence();
    auto exchange_vector = x_phi_K(phi_vect);
    world.gop.fence();
    auto exchange_conjugate_vector = phi_xd_K(phi_vect);
    world.gop.fence();
    molresponse::end_timer(world, "response exchange apply");
    molresponse::start_timer(world);

    vecfuncT exchange = exchange_vector + exchange_conjugate_vector;


    auto exchange_matrix = create_response_matrix(x.size(), n * phi0.size());
    int b = 0;
    for (auto &xi: exchange_matrix) {
        for (auto &xij: xi) {
            xij = copy(exchange[b++]);
        }
    }
    world.gop.fence();
    molresponse::end_timer(world, "response exchange reorganize");
    return exchange_matrix;
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