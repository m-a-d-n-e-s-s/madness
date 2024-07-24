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

using ResponseInput = std::tuple<std::string, std::string, std::vector<double>>;

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
    json result = {};
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

class ResponsePathStrategy : public PathStrategy {

  std::string calc_name = "response";
  std::string perturbation;
  std::string xc;
  std::vector<double> frequencies;
  // sets the current path to the save path
  /**
     * Generates the frequency save path with format
     * /frequency_run_path/restart_[frequency_run_filename].00000
     *
     * @param frequency_run_path
     * @return
     */
  static path restart_path(const std::filesystem::path& calc_path) {

    auto save_path = std::filesystem::path(calc_path);
    auto run_name = calc_path.filename();
    std::string save_string = "restart_" + run_name.string();
    save_path += "/";
    save_path += save_string;
    save_path += ".00000";

    return save_path;
  }

  /**
     * generates the frequency response path using the format
     * [property]_[xc]_[1-100]
     *
     * where 1-100 corresponds a frequency of 1.100
     *
     * @param moldft_path
     * @param property
     * @param frequency
     * @param xc
     * @return
     */
  [[nodiscard]] auto calc_path(const path& root, const double& frequency) const
      -> std::filesystem::path {
    std::string s_frequency = std::to_string(frequency);
    auto sp = s_frequency.find('.');
    s_frequency = s_frequency.replace(sp, sp, "-");
    std::string run_name = perturbation + "_" + xc + "_" + s_frequency;
    return root / std::filesystem::path(run_name);
  }

 public:
  json generateCalcPaths(const path& root) override {

    json paths;
    paths[calc_name] = {};
    auto& response = paths[calc_name];
    response["calculation"] = {};
    response["output"] = {};
    response["restart"] = {};

    auto base_path = root / calc_name;

    for (const auto& frequency : frequencies) {
      auto frequency_run_path = calc_path(base_path, frequency);
      restart_path(frequency_run_path);
      response["calculation"].push_back(frequency_run_path);
      response["restart"].push_back(restart_path(frequency_run_path));
      response["output"].push_back(frequency_run_path / "response_base.json");
    }
    return paths;
  }

  explicit ResponsePathStrategy() = delete;
  explicit ResponsePathStrategy(std::string calc_name, ResponseInput input)
      : calc_name(std::move(calc_name)),
        perturbation(std::get<0>(input)),
        xc(std::get<1>(input)),
        frequencies(std::get<2>(input)) {}
  explicit ResponsePathStrategy(ResponseInput input)
      : ResponsePathStrategy("response", std::move(input)) {}
};
