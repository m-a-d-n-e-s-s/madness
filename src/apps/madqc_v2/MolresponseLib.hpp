#pragma once

#include <filesystem>

#include "../molresponse_v2/FrequencyLoop.hpp"
#include "../molresponse_v2/GroundStateData.hpp"
#include "../molresponse_v2/PropertyManager.hpp"
#include "../molresponse_v2/ResponseDebugLogger.hpp"
#include "../molresponse_v2/ResponseManager.hpp"
#include "../molresponse_v2/ResponseMetaData.hpp"
#include "../molresponse_v2/StateGenerator.hpp"
#include "Utils.hpp"
#include "funcdefaults.h"

namespace molresponse_lib {

// -----------------------------------------------------------------------------
// Container for structured JSON fragments produced by the workflow
// -----------------------------------------------------------------------------
struct Results {
  nlohmann::json metadata;    // convergence metadata per state
  nlohmann::json properties;  // computed α, β, Raman property tables
};

/**
 * @brief Run the full molecular response & property workflow.
 *
 * @param world      The MADNESS world communicator
 * @param params     Unified parameters containing response and molecule info
 * @param indir      Path to ground-state calculation directory
 * @param outdir     Directory where all outputs will be written
 * @return Results   Structured JSON fragments: metadata + properties
 */
inline Results run_response(World& world, const Params& params, const std::filesystem::path& indir, const std::filesystem::path& outdir) {
  // create top-level output directory
  std::filesystem::create_directories(outdir);

  // --- configure the ground-state archive location ---
  auto& rp = params.get<ResponseParameters>();
  auto& gp = params.get<CalculationParameters>();
  auto& molecule = params.get<Molecule>();
  std::filesystem::path ground_state_archive = indir / gp.prefix() / ".00000";

  // initialize ground-state data and response manager
  GroundStateData ground(world, ground_state_archive.string(), molecule);
  ResponseManager rm(world, gp);

  // generate response states
  StateGenerator state_generator(molecule, gp.protocol(), ground.isSpinRestricted(), rp);
  auto generated_states = state_generator.generateStates();
  if (world.rank() == 0) GeneratedStateData::print_generated_state_map(generated_states.state_map);
  world.gop.fence();

  // initialize metadata
  auto responses_dir = outdir / "responses";
  std::filesystem::create_directories(responses_dir);
  std::string meta_file = (responses_dir / "response_metadata.json").string();

  ResponseMetadata metadata(world, meta_file);
  metadata.initialize_states(generated_states.states);
  if (world.rank() == 0) metadata.print_summary();
  world.gop.fence();

  // debug logger
  ResponseDebugLogger debug_logger((responses_dir / "response_log.json").string(), true);

  // loop over thresholds
  for (double thresh : gp.protocol()) {
    rm.setProtocol(world, ground.getL(), thresh);
    ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(), thresh);
    ground.computePreliminaries(world, *rm.getCoulombOp(), rm.getVtol(), rp.fock_json_file());
    if (world.rank() == 0) madness::print("hamiltonian:\n", ground.Hamiltonian);

    for (auto& state : generated_states.states) {
      if (state.is_converged || state.current_threshold() != thresh) continue;

      computeFrequencyLoop(world, rm, state, ground, metadata, debug_logger);

      if (debug_logger.enabled()) debug_logger.write_to_disk();

      if (state.at_final_threshold()) {
        state.is_converged = true;
        if (world.rank() == 0) madness::print("✓ Final convergence reached for", state.description());
      } else {
        state.advance_threshold();
        if (world.rank() == 0) madness::print("→ advancing to next protocol for", state.description());
      }
    }
  }

  // compute requested properties
  PropertyManager properties(world, (responses_dir / "properties.json").string());
  initialize_property_structure(properties, rp);
  const std::vector<char> dipole_dirs(rp.dipole_directions().begin(), rp.dipole_directions().end());
  const std::vector<char> nuclear_dirs(rp.nuclear_directions().begin(), rp.nuclear_directions().end());

  for (auto const& prop : rp.requested_properties()) {
    if (prop == "polarizability") {
      if (world.rank() == 0) madness::print("▶️ Computing polarizability α...");
      compute_alpha(world, generated_states.state_map, ground, rp.dipole_frequencies(), rp.dipole_directions(), properties);
      properties.save();
      if (world.rank() == 0) properties.print_alpha_table();

    } else if (prop == "hyperpolarizability") {
      if (world.rank() == 0) madness::print("▶️ Computing hyperpolarizability β...");
      compute_hyperpolarizability(world, ground, rp.dipole_frequencies(), dipole_dirs, properties);
      properties.save();
      if (world.rank() == 0) properties.print_beta_table();

    } else if (prop == "raman") {
      auto nuclear_dirs = rp.nuclear_directions();
      auto dipole_dirs = rp.dipole_directions();
      if (world.rank() == 0) madness::print("▶️ Computing Raman response...");
      // compute_Raman(world, ground, rp.dipole_frequencies(), dipole_dirs, properties);
      properties.save();
      if (world.rank() == 0) properties.print_beta_table();
    }
  }

  // global inner-product contributions
  auto contribs = global_inner_contributions();
  if (world.rank() == 0 && !contribs.empty()) {
    std::ofstream out((outdir / "all_inner_contributions.json").string());
    out << std::setw(2) << contribs << std::endl;
    madness::print("📂 Wrote all inner‐product contributions");
  }

  // finalize & stats
  if (world.rank() == 0) madness::print("\n✅ Molecular response & property calculation complete.");
  world.gop.fence();
  world.gop.fence();
  print_stats(world);

  // aggregate JSON results
  Results results;
  results.metadata = metadata.to_json();
  results.properties = properties.to_json();
  return results;
}

}  // namespace molresponse_lib
