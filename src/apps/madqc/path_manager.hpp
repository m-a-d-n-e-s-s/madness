#include <apps/molresponse/response_parameters.h>
#include <madchem.h>
#include <madness/chem/CalculationParameters.h>
#include <madness/chem/SCF.h>
#include <madness/chem/molecule.h>
#include <filesystem>
#include <madness/external/nlohmann_json/json.hpp>
#include <utility>

using path = std::filesystem::path;
using json = nlohmann::json;
using commandlineparser = madness::commandlineparser;

using namespace madness;

class PathStrategy {
 public:
  // Generate the paths for the calculation, restarts, and output files in the form of a json object.
  virtual json generateCalcPaths(const path& root) = 0;
  virtual ~PathStrategy() = default;
};

// Allows us to combine multiple path strategies into a single strategy e.g. moldft, response, mp2, etc.
class CompositePathStrategy : public PathStrategy {
 private:
  std::vector<std::unique_ptr<PathStrategy>> strategies;

 public:
  void addStrategy(std::unique_ptr<PathStrategy> strategy) {
    strategies.push_back(std::move(strategy));
  }

  json generateCalcPaths(const path& root) override {
    json result ={};
    for (const auto& strategy : strategies) {
      json paths = strategy->generateCalcPaths(root);
      // I know there should only be one key in the json object.
      for (auto& [key, value] : paths.items()) {
        result[key] = value;
      }
    }
    return result;
  }
};

// The pass strategies to the PathManager to construct the paths for different types of calculations.
class PathManager {
 private:
  std::unique_ptr<CompositePathStrategy> strategy;

 public:
  PathManager() : strategy(std::make_unique<CompositePathStrategy>()) {}

  void addStrategy(std::unique_ptr<PathStrategy> newStrategy) {
    strategy->addStrategy(std::move(newStrategy));
  }

  json generateCalcPaths(const path& root) {
    return strategy->generateCalcPaths(root);
  }
};

class MoldftPathStrategy : public PathStrategy {

  // The name of the calculation in the case we want multiple
  // types of moldft calculations.
  //
  // moldft is the default name.
  std::string calc_name = "moldft";

 public:
  json generateCalcPaths(const path& root) override {
    json paths;

    paths[calc_name] = {};
    auto& moldft = paths[calc_name];

    auto base_root = root / calc_name;

    moldft["calculation"] = base_root;
    moldft["restart"] = base_root / "moldft.restartdata.00000";
    moldft["output"] = {};
    moldft["output"]["calc_info"] = base_root / "moldft.calc_info.json";
    moldft["output"]["scf_info"] = base_root / "moldft.scf_info.json";

    return paths;
  }

  explicit MoldftPathStrategy() = default;
  explicit MoldftPathStrategy(std::string calc_name)
      : calc_name(std::move(calc_name)) {}
};
