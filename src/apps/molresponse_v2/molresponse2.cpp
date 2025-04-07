#include "../madqc/parameter_manager.hpp"
#include "FrequencyLoop.hpp" // Make sure this is included
#include "GroundStateData.hpp"
#include "MolecularProperty.hpp"
#include "ResponseManager.hpp"
#include "StateGenerator.hpp"
#include "funcdefaults.h"
#include "madness/world/world.h"
#include <iostream>
#include <memory>

using namespace madness;

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  {
    startup(world, argc, argv, true);

    if (argc != 2) {
      if (world.rank() == 0)
        std::cerr << "Usage: molresponse2 [input_file.json]\n";
      finalize();
      return 1;
    }

    // Load input parameters explicitly from JSON file
    ParameterManager params;
    params = ParameterManager(world, path(argv[1]));

    auto molresponse_params = params.get_molresponse_params();
    std::string fock_json_file = molresponse_params.fock_json_file();
    Molecule molecule = params.get_molecule();
    auto protocol = molresponse_params.protocol();
    molresponse_params.set_user_defined_value<std::string>(
        "archive", "moldft.restartdata");

    auto ground_state_archive = molresponse_params.archive();
    if (world.rank() == 0) {
      std::cout << "Ground state archive: " << ground_state_archive
                << std::endl;
    }

    std::vector<MolecularProperty> requested_properties = {
        MolecularProperty(MolecularPropertyType::Polarizability, {0.0, 0.0656},
                          {'z'}),
    };
    /*MolecularProperty(MolecularPropertyType::Raman, {0.0656})};*/

    // Initialize the ResponseManager with ground-state archive
    auto ground_state =
        GroundStateData(world, ground_state_archive, params.get_molecule());
    ResponseManager rm =
        ResponseManager(world, params.get_molresponse_params());

    std::vector<double> protocols = {1e-4, 1e-6};
    StateGenerator state_generator(molecule, requested_properties, protocols);
    auto all_states = state_generator.generateStates();

    bool all_states_converged = false;
    while (!all_states_converged) {
      all_states_converged = true;

      // Extract all unique thresholds needed for this round
      std::set<double> active_thresholds;
      for (const auto &state : all_states)
        if (!state.is_converged)
          active_thresholds.insert(state.current_threshold());

      for (double thresh : active_thresholds) {
        rm.setProtocol(world, ground_state.getL(), thresh);
        ground_state.prepareOrbitals(world, FunctionDefaults<3>::get_k(),
                                     thresh);
        ground_state.computePreliminaries(world, *rm.getCoulombOp(),
                                          rm.getVtol(), fock_json_file);

        for (auto &state : all_states) {
          if (state.is_converged || state.current_threshold() != thresh)
            continue;
          ResponseMetadata metadata(state.perturbationDescription());
          metadata.load();
          computeFrequencyLoop(world, rm, state, ground_state, metadata);
          metadata.save();
          // Check if we reached final protocol or should advance
          if (state.at_final_threshold()) {
            state.is_converged = true;
            if (world.rank() == 0) {
              madness::print("✓ Final convergence reached for",
                             state.description());
            }
          } else {
            state.advance_threshold();
            all_states_converged = false;
            if (world.rank() == 0) {

              madness::print("→ Converged at thresh", thresh,
                             "→ advancing to next protocol for",
                             state.description());
            }
          }
        }
      }
    }

    // Example placeholder for a future implementation
    // ResponseState response_state(freq, thresh, "x");
    // auto response_functions =
    // response_manager.computeResponse(response_state);
    // response_manager.saveResponse(response_state, response_functions);

    // Currently just logging the progress explicitly
    world.gop.fence();

    if (world.rank() == 0) {
      madness::print(
          "\nMolecular response calculation completed successfully.");
    }
  }

  finalize();
  return 0;
}
