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
    ground_calculation.print_params();
    Molecule molecule = ground_calculation.molecule();
    r_params.set_ground_state_calculation_data(ground_calculation);
    r_params.set_derived_values(world, molecule);
    r_params.print();
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
auto newK(const vecfuncT &ket, const vecfuncT &bra, const vecfuncT &vf) -> vecfuncT {
    World &world = ket[0].world();
    const double lo = 1.e-10;

    Exchange<double, 3> op{};
    op.set_parameters(bra, ket, lo);
    op.set_algorithm(op.multiworld_efficient);
    return op(vf);
}
// sum_i |i><i|J|p> for each p