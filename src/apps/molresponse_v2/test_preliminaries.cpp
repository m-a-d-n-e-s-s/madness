#include "../madqc/parameter_manager.hpp"
#include "ResponseManager.hpp"
#include "madness/world/world.h"
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>

using namespace madness;

// Utility function to compare values within tolerance
bool approximately_equal(double computed, double reference, double tol = 1e-6) {
  return std::abs(computed - reference) <= tol;
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  {

    startup(world, argc, argv, true);
    if (world.rank() == 0)
      print(info::print_revision_information());

    if (argc != 3) {
      if (world.rank() == 0)
        std::cerr
            << "Usage: test_preliminaries [input_file.json] [reference.txt]\n";
      finalize();
      return 1;
    }

    const std::string input_filename = argv[1];
    const std::string reference_filename = argv[2];
    commandlineparser parser(argc, argv);
    ParameterManager params;
    // Initialize the necessary components
    path input_file(argv[1]);
    // Extract relevant parameters
    const std::string ground_state_archive = "moldft.restartdata";
    // Initialize the ResponseManager with ground-state archive
    ResponseManager response_manager(world, ground_state_archive,
                                     params.get_molecule());

    auto protocol = params.get_moldft_params().protocol();

    for (auto thresh : protocol) {
      if (world.rank() == 0)
        print("Setting protocol with threshold: ", thresh);
      response_manager.setProtocol(thresh);
      auto k = FunctionDefaults<3>::get_k();
      if (world.rank() == 0) {

        print("K is now: ", k);
      }
      response_manager.prepareOrbitalsForAccuracyStep();
      response_manager.computePreliminaries();

      auto prelims = response_manager.currentPreliminaries();

      if (world.rank() == 0) {
        print("Computed: \n", prelims.Hamiltonian, "\n ",
              prelims.Hamiltonian_no_diag);
      }
    }
    world.gop.fence();
    world.gop.fence();
    print_stats(world);
  }
  finalize();
  return 0;
}

