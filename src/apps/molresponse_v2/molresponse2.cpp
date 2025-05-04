#include <iostream>

#include "FrequencyLoop.hpp"  // Make sure this is included
#include "GroundStateData.hpp"
#include "MolecularProperty.hpp"
#include "ParameterManager.hpp"
#include "PropertyManager.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseManager.hpp"
#include "ResponseMetaData.hpp"
#include "StateGenerator.hpp"
#include "funcdefaults.h"

using namespace madness;

int main(int argc, char **argv) {
  // Initialize MADNESS world (passes MPI args, etc.)
  World &world = initialize(argc, argv);
  {
    startup(world, argc, argv, true);
    if (argc != 2) {
      if (world.rank() == 0) std::cerr << "Usage: molresponse2 [input_file.json]\n";
      finalize();
      return 1;
    }
    // Load input parameters explicitly from JSON file
    // Define a concrete aliased ParameterManager type
    using MyParamMgr = ParameterManager<CalculationParameters, ResponseParameters, OptimizationParameters, Molecule>;
    commandlineparser parser(argc, argv);
    if (argc != 2) {
      std::cerr << "Usage: " << argv[0] << " <input_file>\n";
      return 1;
    }
    std::string input_file = argv[1];
    // Construct the manager, reading .inp or JSON as needed
    MyParamMgr pm(world, input_file);

    auto rp = pm.get<ResponseParameters>();
    std::string fock_json_file = rp.fock_json_file();
    Molecule molecule = pm.get<Molecule>();
    auto params = pm.get<CalculationParameters>();
    auto protocol = params.protocol();
    rp.set_user_defined_value<std::string>("archive", "moldft.restartdata");
    auto ground_state_archive = rp.archive();
    if (world.rank() == 0) {
      std::cout << "Ground state archive: " << ground_state_archive << std::endl;
    }
    // Initialize the ResponseManager with ground-state archive
    auto ground_state = GroundStateData(world, ground_state_archive, molecule);
    ResponseManager rm = ResponseManager(world, params);
    std::vector<double> protocols = {1e-4, 1e-6};
    StateGenerator state_generator(molecule, protocols, ground_state.isSpinRestricted(), rp);
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
      ground_state.computePreliminaries(world, *rm.getCoulombOp(), rm.getVtol(), fock_json_file);
      if (world.rank() == 0) {
        print("hamiltonian:\n", ground_state.Hamiltonian);
      }
      for (auto &state : generated_states.states) {
        if (state.is_converged || state.current_threshold() != thresh) continue;

        computeFrequencyLoop(world, rm, state, ground_state, metadata, debug_logger);

        if (debug_logger.enabled()) {
          debug_logger.write_to_disk("response_log.json");
        }
        // Check if we reached final protocol or should advance
        if (state.at_final_threshold()) {
          state.is_converged = true;
          if (world.rank() == 0) {
            madness::print("✓ Final convergence reached for", state.description());
          }
        } else {
          state.advance_threshold();
          all_states_converged = false;
          if (world.rank() == 0) {
            madness::print("→ Converged at thresh", thresh, "→ advancing to next protocol for", state.description());
          }
        }
      }
    }

    // After all_states have converged & been saved…

    // 1) Initialize PropertyManager
    PropertyManager properties(world, "responses/properties.json");
    initialize_property_structure(properties, rp);
    if (world.rank() == 0) {
      properties.print_alpha_table();
      properties.print_beta_table();
    }

    auto dipole_dirs = rp.dipole_directions();
    auto dipole_freqs = rp.dipole_frequencies();
    auto nuclear_dirs = rp.nuclear_directions();
    auto nuclear_atom_indices = rp.nuclear_atom_indices();
    auto nuclear_freqs = rp.nuclear_frequencies();

    // 2) Loop over requested properties and dispatch
    for (auto const &prop : rp.requested_properties()) {
      if (prop == "polarizability") {
        if (world.rank() == 0)
          print("▶️ Computing polarizability α for directions ", std::string(dipole_dirs.begin(), dipole_dirs.end()), " at freqs ", dipole_dirs);
        compute_alpha(world, generated_states.state_map, ground_state, dipole_freqs, std::string(dipole_dirs.begin(), dipole_dirs.end()), properties);
        properties.save();
        if (world.rank() == 0) properties.print_alpha_table();
      } else if (prop == "hyperpolarizability") {
        if (world.rank() == 0) print("▶️ Computing hyperpolarizability β for directions ", dipole_dirs, " at freqs ", dipole_freqs);
        compute_hyperpolarizability(world, ground_state, dipole_freqs, std::vector<char>{dipole_dirs.begin(), dipole_dirs.end()}, properties);
        properties.save();
        if (world.rank() == 0) properties.print_beta_table();
      } else if (prop == "raman") {
        if (world.rank() == 0)
          print("▶️ Computing Raman response for directions ", std::string(dipole_dirs.begin(), dipole_dirs.end()), " at freqs ", dipole_freqs);
        compute_Raman(world, ground_state, {dipole_freqs, dipole_freqs}, std::vector<char>{dipole_dirs.begin(), dipole_dirs.end()}, properties);
        properties.save();
        if (world.rank() == 0) properties.print_beta_table();
      }
    }

    // 3) Final message
    if (world.rank() == 0) madness::print("\n✅ Molecular response & property calculation complete.");

    // Currently just logging the progress explicitly
    world.gop.fence();
    world.gop.fence();
    print_stats(world);
  }
  finalize();
  return 0;
}
