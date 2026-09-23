#include <madness/chem/ResultsSummary.hpp>
#include <madness/mra/mra.h>
#include <madness/world/test_utilities.h>

#include <sstream>
#include <string>

using namespace madness;

static bool has_text(const std::string &s, const std::string &needle) {
  return s.find(needle) != std::string::npos;
}

int main(int argc, char **argv) {
  World &world = madness::initialize(argc, argv);
  int success = 0;

  if (world.rank() == 0) {
    test_output t("results summary citations");

    const nlohmann::json base_scf = {
        {"model", "scf"},
        {"properties", {{"energy", -1.0}}},
        {"citations", {{"dftd3", false}, {"pcm", false}, {"libxc", false}}}};

    {
      std::ostringstream os;
      qcapp::write_results_summary(os, {{"tasks", nlohmann::json::array({base_scf})}});
      const auto out = os.str();
      const bool ok = !has_text(out, "Citations for external modules");
      t.checkpoint(ok, "omits citations block when no module is requested");
      success += !ok;
    }

    {
      std::ostringstream os;
      nlohmann::json scf = base_scf;
      scf["citations"]["dftd3"] = true;
      scf["citations"]["pcm"] = true;
      scf["citations"]["libxc"] = true;
      qcapp::write_results_summary(os, {{"tasks", nlohmann::json::array({scf})}});
      const auto out = os.str();
      const bool ok = has_text(out, "Citations for external modules") &&
                      has_text(out, "DFT-D3 / simple-dftd3") &&
                      has_text(out, "PCMSolver (PCM)") &&
                      has_text(out, "Libxc");
      t.checkpoint(ok, "prints module citations only when requested");
      success += !ok;
    }
  }

  world.gop.fence();
  madness::finalize();
  return success;
}
