#pragma once

#include <filesystem>
#include "InputWriter.hpp"
#include "ParameterManager.hpp"
#include <GroundStateData.hpp>

struct cc_lib {
  // -----------------------------------------------------------------------------
  // Container for structured JSON fragments produced by the workflow
  // -----------------------------------------------------------------------------
  struct Results {
    double energy;
    tensorT dipole;
    std::optional<Tensor<double>> gradient;

    nlohmann::json metadata;    // convergence metadata per state
    nlohmann::json properties;  // computed α, β, Raman property tables

    std::optional<nlohmann::json>
        debug_log;  // debug log of response calculations
  };
  static constexpr char const* label() { return "cc2"; }

  /**
   * @brief Run the full molecular response & property workflow.
   *
   * @param world      The MADNESS world communicator
   * @param params     Unified parameters containing response and molecule info
   * @param indir      Path to ground-state calculation directory
   * @param outdir     Directory where all outputs will be written
   * @return Results   Structured JSON fragments: metadata + properties
   */
  inline static Results run_cc2(World& world, const Params& params,
                                     const std::filesystem::path& indir,
                                     const std::filesystem::path& outdir) {
    // --- configure the ground-state archive location ---
    auto rp = params.get<CCParameters>();
    const auto& gp = params.get<CalculationParameters>();
    const auto& molecule = params.get<Molecule>();

    if (world.rank() == 0) {
      json response_input_json = {};
      response_input_json["response"] = rp.to_json_if_precedence("defined");
      print("response_input_json: ", response_input_json.dump(4));
      std::ofstream ofs("response.in");
      write_json_to_input_file(response_input_json, {"response"}, ofs);
      ofs.close();
    }
    world.gop.fence();
    commandlineparser parser;
    parser.set_keyval("input", "response.in");
    if (world.rank() == 0) ::print("input filename: ", parser.value("input"));


    auto rel = std::filesystem::relative(indir, outdir);
    auto prox = std::filesystem::proximate(indir, outdir);
    auto prefix = std::filesystem::path(indir).stem().string();
    if (world.rank() == 0) {
      std::cout << "Running cc2 calculation in: " << outdir << std::endl;
      std::cout << "Ground state archive: " << indir << std::endl;
      std::cout << "Relative path: " << rel << std::endl;
      std::cout << "Proximate path: " << prox << std::endl;
    }

    std::string archive_name = "moldft.restartdata";
    std::string archive_file = archive_name + ".00000";
    std::string fock_json_file = prox / "moldft.fock.json";
    auto relative_archive = prox / archive_name;

    // initialize ground-state data and response manager
    GroundStateData ground(world, relative_archive.string(), molecule);

    // aggregate JSON results
    Results results;
    return results;
  }

};  // namespace molresponse_lib
