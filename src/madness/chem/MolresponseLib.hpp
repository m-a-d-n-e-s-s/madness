/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA

  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov
  tel:   865-241-3937
  fax:   865-572-0680


  $Id$
 */
#pragma once

#include <apps/molresponse_v2/GroundStateData.hpp>
#include <apps/molresponse_v2/PropertyManager.hpp>
#include <apps/molresponse_v2/ResponseDebugLogger.hpp>
#include <apps/molresponse_v2/ResponseManager.hpp>
#include <apps/molresponse_v2/ResponseRecord.hpp>
#include <apps/molresponse_v2/StateGenerator.hpp>
#include <filesystem>
#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/Results.h>

struct molresponse_lib
{
  // -----------------------------------------------------------------------------
  // Container for structured JSON fragments produced by the workflow
  // -----------------------------------------------------------------------------
  struct Results
  {
    nlohmann::json metadata;             // convergence metadata per state
    nlohmann::json properties;           // computed α, β, Raman property tables
    nlohmann::json vibrational_analysis; // vibrational analysis results
    nlohmann::json debug_log;            // debug log of response calculations
  };
  static constexpr const char* label()
  {
    return "molresponse";
  }

  /**
   * @brief Run the full molecular response & property workflow.
   *
   * @param world      The MADNESS world communicator
   * @param params     Unified parameters containing response and molecule info
   * @param outdir     Directory where all outputs will be written
   * @return Results   Structured JSON fragments: metadata + properties
   */
  inline static Results run_response(World& world, const Params& params,
                                     const std::shared_ptr<SCF>& scf_calc,
                                     const std::filesystem::path& outdir)
  {
    // --- configure the ground-state archive location ---
    const auto& calc_params = params.get<CalculationParameters>();
    const auto& molecule = scf_calc->molecule;
    const auto& rp_copy = params.get<ResponseParameters>();

    if (world.rank() == 0)
    {
      json response_input_json = {};
      response_input_json["response"] = rp_copy.to_json_if_precedence("defined");
      print("response_input_json: ", response_input_json.dump(4));
      std::ofstream ofs("response.in");
      write_json_to_input_file(response_input_json, {"response"}, ofs);
      ofs.close();
    }
    world.gop.fence();
    commandlineparser parser;
    parser.set_keyval("input", "response.in");
    if (world.rank() == 0)
    {
      ::print("input filename: ", parser.value("input"));
    }

    auto response_params = ResponseParameters(world, parser);
    auto indir = scf_calc->work_dir;
    auto rel = std::filesystem::relative(indir, outdir);
    auto prox = std::filesystem::proximate(indir, outdir);
    if (world.rank() == 0)
    {
      print("Running MolresponseLib::run_response() in directory: ", outdir);
      print("Ground state archive: ", indir);
      print("Relative path: ", rel);
      print("Proximate path: ", prox);
    }

    std::string archive_name = "mad.restartdata";
    std::string fock_json_file = prox / "moldft.fock.json";
    auto relative_archive = prox / archive_name;

    // // initialize ground-state data and response manager
    // if (scf_calc){
    //   ground = GroundStateData(world, scf_calc);
    // }else{

    GroundStateData ground(world, relative_archive.string(), molecule);
    ResponseManager rm(world, calc_params);

    // generate response states
    StateGenerator state_generator(molecule, calc_params.protocol(), ground.isSpinRestricted(),
                                   response_params);
    auto generated_states = state_generator.generateStates();
    if (world.rank() == 0)
    {
      GeneratedStateData::print_generated_state_map(generated_states.state_map);
    }
    world.gop.fence();

    // initialize metadata
    std::string meta_file = "response_metadata.json";

    ResponseRecord response_record(world, meta_file);
    response_record.initialize_states(generated_states.states);
    if (world.rank() == 0)
    {

      response_record.print_summary();
    }
    world.gop.fence();

    // debug logger
    ResponseDebugLogger debug_logger("response_log.json", true);

    // loop over thresholds
    for (double thresh : calc_params.protocol())
    {
      rm.setProtocol(world, ground.getL(), thresh);
      ground.prepareOrbitals(world, FunctionDefaults<3>::get_k(), thresh);
      ground.computePreliminaries(world, *rm.getCoulombOp(), rm.getVtol(), fock_json_file);
      // if (world.rank() == 0)
      //   madness::print("hamiltonian:\n", ground.Hamiltonian);

      for (auto& state : generated_states.states)
      {
        //     if (state.is_converged || state.current_threshold() != thresh)
        //     continue;

        computeFrequencyLoop(world, rm, state, ground, response_record, debug_logger);

        if (debug_logger.enabled())
          debug_logger.write_to_disk();
        // Check if we reached final protocol or should advance
        if (state.at_final_threshold())
        {
          if (world.rank() == 0)
          {
            madness::print("✓ Final convergence reached for", state.description());
          }
        }
        else
        {
          state.advance_threshold();
          if (world.rank() == 0)
          {
            madness::print("→ Converged at thresh", thresh, "→ advancing to next protocol for",
                           state.description());
          }
        }
      }
    }

    // compute requested properties
    PropertyManager properties(world, "properties.json");
    std::string dip_dirs = response_params.dipole_directions();
    enum class PropertyType
    {
      Alpha,
      Beta,
      Raman
    };

    std::map<std::string, PropertyType> prop_map = {{"polarizability", PropertyType::Alpha},
                                                    {"hyperpolarizability", PropertyType::Beta},
                                                    {"raman", PropertyType::Raman}};

    VibrationalResults vibrational_results;
    for (const std::string& prop : response_params.requested_properties())
    {
      // get rid of first and last characters
      auto prop_type = prop_map[prop.substr(1, prop.size() - 2)];
      if (prop_type == PropertyType::Alpha)
      {
        if (world.rank() == 0)
        {
          madness::print("▶️ Computing polarizability α...");
        }
        compute_alpha(world, generated_states.state_map, ground,
                      response_params.dipole_frequencies(), response_params.dipole_directions(),
                      properties);
        properties.save();
      }
      else if (prop_type == PropertyType::Beta)
      {
        if (world.rank() == 0)
        {

          madness::print("▶️ Computing hyperpolarizability β...");
        }

        compute_hyperpolarizability(world, ground, response_params.dipole_frequencies(), dip_dirs,
                                    properties);
        properties.save();
      }
      else if (prop_type == PropertyType::Raman)
      {

        // vibrational analysis (Hessian + frequencies + intensities)
        if (world.rank() == 0)
        {

          madness::print("▶️ Computing Hessian...");
        }
        auto vib = compute_hessian(world, generated_states.state_map, ground,
                                   response_params.dipole_directions(), scf_calc);
        vibrational_results = vib;
        if (world.rank() == 0)
        {

          madness::print("▶️ Computing Raman response...");
        }
        compute_Raman(world, ground, response_params.dipole_frequencies(),
                      response_params.nuclear_atom_indices(), response_params.dipole_directions(),
                      response_params.nuclear_directions(), properties);
        properties.save();
      }
      if (world.rank() == 0)
      {
        properties.print_table();
      }
    }

    // // global inner-product contributions
    // auto contribs = global_inner_contributions();
    // if (world.rank() == 0 && !contribs.empty())
    // {
    //   std::ofstream out("all_inner_contributions.json");
    //   print(std::setw(2),contribs);
    //   out << std::setw(2) << contribs << std::endl;
    //   madness::print("📂 Wrote all inner‐product contributions");
    // }

    // finalize & stats
    if (world.rank() == 0)
    {

      madness::print("\n✅ Molecular response & property calculation complete.");
    }
    world.gop.fence();
    world.gop.fence();
    print_stats(world);

    // aggregate JSON results
    Results results;
    results.metadata = response_record.to_json();
    results.properties = properties.to_json();
    results.debug_log = debug_logger.to_json();
    results.vibrational_analysis = vibrational_results.to_json();
    return results;
  }

}; // namespace molresponse_lib
