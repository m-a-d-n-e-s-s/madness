#pragma once

#include <filesystem>
#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <apps/molresponse_v2/GroundStateData.hpp>

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
    const auto& ccparam = params.get<CCParameters>();
    const auto& scfparam = params.get<CalculationParameters>();
    const auto& molecule = params.get<Molecule>();

    auto rel = std::filesystem::relative(indir, outdir);
    if (world.rank() == 0) {
      std::cout << "Running cc2 calculation in: " << outdir << std::endl;
      std::cout << "Ground state archive: " << indir << std::endl;
      std::cout << "Relative path: " << rel << std::endl;
    }

    // aggregate JSON results
    Results results;
    return results;
  }

};  // namespace molresponse_lib
