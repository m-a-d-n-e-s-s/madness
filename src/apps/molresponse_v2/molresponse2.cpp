#include <iostream>

#include <madness/chem/oep.h>
#include <madness/chem/TDHF.h>
#include "FrequencyLoop.hpp"  // Make sure this is included
#include "GroundStateData.hpp"
#include <madness/chem/ParameterManager.hpp>
#include "PropertyManager.hpp"
#include "ResponseDebugLogger.hpp"
#include "ResponseManager.hpp"
#include "ResponseMetaData.hpp"
#include "StateGenerator.hpp"
#include <madness/mra/funcdefaults.h>

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
    // std::string input_file = argv[1];
    // Construct the manager, reading .inp or JSON as needed
    MyParamMgr pm(world, parser);

    const auto &tmpl = pm.getAllInputJson();
    if (world.rank() == 0) {
      std::cout << "Input JSON: " << parser.value("input") << "\n" << std::setw(2) << tmpl.dump(2) << std::endl;
    }

    auto rp = pm.get<ResponseParameters>();
    std::string fock_json_file = rp.fock_json_file();
    const auto &molecule = pm.get<Molecule>();
    const auto &params = pm.get<CalculationParameters>();
    auto protocol = params.protocol();
    rp.set_user_defined_value<std::string>("archive", "moldft.restartdata");
    auto ground_state_archive = rp.archive();
    if (world.rank() == 0) {
      std::cout << "Ground state archive: " << ground_state_archive << std::endl;
    }
    // Initialize the ResponseManager with ground-state archive
    auto ground_state = GroundStateData(world, ground_state_archive, molecule);
    ResponseManager rm = ResponseManager(world, params);
    StateGenerator state_generator(molecule, params.protocol(), ground_state.isSpinRestricted(), rp);
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
    ResponseDebugLogger debug_logger("responses/response_log.json", true);
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
          debug_logger.write_to_disk();
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

   // compute requested properties
    PropertyManager properties(world, "properties.json");
    std::string dip_dirs = rp.dipole_directions();
    std::string nuc_dirs = rp.nuclear_directions();

    for (auto const& prop : rp.requested_properties()) {
      if (prop == "polarizability") {
        if (world.rank() == 0)
          madness::print("▶️ Computing polarizability α...");
        compute_alpha(world, generated_states.state_map, ground_state,
                      rp.dipole_frequencies(), rp.dipole_directions(),
                      properties);
        properties.save();

      } else if (prop == "hyperpolarizability") {
        if (world.rank() == 0)
          madness::print("▶️ Computing hyperpolarizability β...");

        compute_hyperpolarizability(world, ground_state, rp.dipole_frequencies(),
                                    dip_dirs, properties);
        properties.save();

      } else if (prop == "raman") {
        auto nuclear_dirs = rp.nuclear_directions();
        auto dipole_dirs = rp.dipole_directions();
        if (world.rank() == 0) madness::print("▶️ Computing Raman response...");
        // compute_Raman(world, ground, rp.dipole_frequencies(), dipole_dirs,
        // properties);
        properties.save();
      }
      if (world.rank() == 0) {
        properties.print_table();
      }
    }

    // global inner-product contributions
    auto contribs = global_inner_contributions();
    if (world.rank() == 0 && !contribs.empty()) {
      std::ofstream out("all_inner_contributions.json");
      out << std::setw(2) << contribs << std::endl;
      madness::print("📂 Wrote all inner‐product contributions");
    }
    // Currently just logging the progress explicitly
    world.gop.fence();
    world.gop.fence();
    print_stats(world);
  }
  finalize();
  return 0;
}
