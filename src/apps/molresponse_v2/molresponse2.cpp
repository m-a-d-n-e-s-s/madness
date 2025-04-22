#include "../madqc/parameter_manager.hpp"
#include "FrequencyLoop.hpp" // Make sure this is included
#include "GroundStateData.hpp"
#include "MolecularProperty.hpp"
#include "PropertyManager.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseManager.hpp"
#include "ResponseMetaData.hpp"
#include "StateGenerator.hpp"
#include "VBCMacrotask.hpp"
#include "funcdefaults.h"
#include "madness/world/world.h"
#include <iostream>

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

    auto input_freq = molresponse_params.freq_range();

    std::vector<MolecularProperty> requested_properties = {
        MolecularProperty(MolecularPropertyType::Polarizability, input_freq,
                          {'x', 'y', 'z'}),
        MolecularProperty(MolecularPropertyType::Hyperpolarizability,
                          input_freq, {'x', 'y', 'z'}),
    };
    PropertyManager prop_manager("responses/properties.json");
    initialize_property_structure(prop_manager, requested_properties);
    if (world.rank() == 0) {
      prop_manager.print_alpha_table();
      prop_manager.print_beta_table();
    }

    // Initialize the ResponseManager with ground-state archive
    auto ground_state =
        GroundStateData(world, ground_state_archive, params.get_molecule());
    ResponseManager rm =
        ResponseManager(world, params.get_molresponse_params());

    std::vector<double> protocols = {1e-4, 1e-6};
    StateGenerator state_generator(molecule, requested_properties, protocols,
                                   ground_state.isSpinRestricted());
    auto generated_states = state_generator.generateStates();

    if (world.rank() == 0) {
      GeneratedStateData::print_generated_state_map(generated_states.state_map);
    }
    world.gop.fence();
    bool all_states_converged = false;
    std::string response_metadata_file = "responses/response_metadata.json";
    ResponseMetadata metadata(world, response_metadata_file);
    metadata.initialize_states(generated_states.states);
    if (world.rank() == 0) {
      metadata.print_summary();
    }
    world.gop.fence();
    ResponseDebugLogger debug_logger(true);
    // Extract all unique thresholds needed for this round
    for (double thresh : protocol) {
      rm.setProtocol(world, ground_state.getL(), thresh);
      ground_state.prepareOrbitals(world, FunctionDefaults<3>::get_k(), thresh);
      ground_state.computePreliminaries(world, *rm.getCoulombOp(), rm.getVtol(),
                                        fock_json_file);
      if (world.rank() == 0) {
        print("hamiltonian:\n", ground_state.Hamiltonian);
      }

      for (auto &state : generated_states.states) {
        if (state.is_converged || state.current_threshold() != thresh)
          continue;

        const auto state_id = state.perturbationDescription();
        if (metadata.final_converged(state_id)) {
          continue;
        }

        computeFrequencyLoop(world, rm, state, ground_state, metadata,
                             debug_logger);
        if (debug_logger.enabled()) {
          debug_logger.write_to_disk("response_log.json");
        }
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

    compute_alpha(world, generated_states.state_map, ground_state, input_freq,
                  "xy", prop_manager);
    prop_manager.save();
    if (world.rank() == 0) {
      prop_manager.print_alpha_table();
    }
    auto t1 = VBC_save_load(ground_state, "xyz", input_freq,
                            ground_state.isSpinRestricted());

    MacroTask task_one(world, t1);
    task_one.set_debug(true);
    std::vector<int> freq = {0, 1, 2, 3, 4, 5};
    auto result_one = task_one(freq);

    // Example placeholder for a future implementation
    // ResponseState response_state(freq, thresh, "x");
    // auto response_functions =
    // response_manager.computeResponse(response_state);
    // response_manager.saveResponse(response_state, response_functions);

    if (world.rank() == 0) {
      madness::print(
          "\nMolecular response calculation completed successfully.");
    }
    // Currently just logging the progress explicitly
    world.gop.fence();
    world.gop.fence();
    print_stats(world);
  }
  finalize();
  return 0;
}
