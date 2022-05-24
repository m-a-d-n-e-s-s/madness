#include "global_functions.h"

#include <memory>
#include <string>
#include <vector>

#include "response_parameters.h"

void print_molecule(World &world, const GroundStateCalculation & g_params) {
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
      std::cout << std::setw(18) << std::right << coords[0] << std::setw(18) << coords[1] << std::setw(18) << coords[2]
                << endl;
    }
    print("");

    // Reset precision
    std::cout.precision(10);
    std::cout << std::scientific;
  }
}

CalcParams initialize_calc_params(World &world, std::string input_file) {
  ResponseParameters r_params{};
  commandlineparser parser;
  parser.set_keyval("input",input_file);
  r_params.read_input_and_commandline_options(world, parser, "response");
  GroundStateCalculation ground_calculation{world};
  ground_calculation.print_params();
  Molecule molecule = ground_calculation.molecule();
  r_params.set_ground_state_calculation_data(ground_calculation);
  r_params.set_derived_values(world,  molecule);
  r_params.print();
  return CalcParams{ground_calculation, molecule, r_params};
}
