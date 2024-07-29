#ifndef MADCHEM_PATH_MANAGER_HPP
#define MADCHEM_PATH_MANAGER_HPP
//
//
// Purpose: Contains the classes for managing the paths for the calculations.
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

class ResponseConfig {

 public:
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
    std::string run_name =
        this->perturbation + "_" + this->xc + "_" + s_frequency;
    return root / std::filesystem::path(run_name);
  }
  explicit ResponseConfig(std::string calc_name, ResponseInput input)
      : calc_name(std::move(calc_name)),
        perturbation(std::get<0>(input)),
        xc(std::get<1>(input)),
        frequencies(std::get<2>(input)) {}
  explicit ResponseConfig(ResponseInput input)
      : perturbation(std::get<0>(input)),
        xc(std::get<1>(input)),
        frequencies(std::get<2>(input)) {}
};

class ResponsePathStrategy : public PathStrategy {

 private:
  std::string calc_name = "response";
  ResponseConfig config;

 public:
  json generateCalcPaths(const path& root) override {

    json paths;
    paths[calc_name] = {};
    auto& response = paths[calc_name];
    response["calculation"] = {};
    response["output"] = {};
    response["restart"] = {};

    auto base_path = root / calc_name;
    // TODO: (@ahurta92) Feels hacky, should I always be creating new directories?  If I am, should I just create all of them from the start?
    // I added this because later down the line I was getting an error that the directory didn't exist for /base/response when trying to create /base/response/frequency
    if (!std::filesystem::exists(base_path)) {
      std::filesystem::create_directory(base_path);
    }

    response["properties"] = {};
    response["properties"]["alpha"] = base_path / "alpha.json";
    response["properties"]["beta"] = base_path / "beta.json";

    for (const auto& frequency : config.frequencies) {
      auto frequency_run_path = config.calc_path(base_path, frequency);
      auto restart_path = ResponseConfig::restart_path(frequency_run_path);
      response["calculation"].push_back(frequency_run_path);
      response["restart"].push_back(restart_path);
      response["output"].push_back(frequency_run_path / "response_base.json");
    }
    return paths;
  }

  explicit ResponsePathStrategy(ResponseInput input)
      : config(std::move(input)) {}

  explicit ResponsePathStrategy(ResponseInput input,
                                const std::string& calc_name)
      : config(calc_name, std::move(input)) {}
};

class HyperPolarizabilityPathStrategy : public PathStrategy {

  ResponseConfig config;

 public:
  json generateCalcPaths(const path& root) override {

    json paths;
    paths[config.calc_name] = {};
    auto& response = paths[config.calc_name];
    response["calculation"] = {};
    response["output"] = {};
    response["restart"] = {};

    auto base_path = root / config.calc_name;
    // TODO: (@ahurta92) Feels hacky, should I always be creating new directories?  If I am, should I just create all of them from the start?
    // I added this because later down the line I was getting an error that the directory didn't exist for /base/response when trying to create /base/response/frequency
    if (!std::filesystem::exists(base_path)) {
      std::filesystem::create_directory(base_path);
    }

    // TODO: @ahurta92 This could include a few different stratgies dependent on what the user wants

    auto set_freqs = [&]() {
      vector<double> freqs_copy = config.frequencies;
      auto num_freqs = config.frequencies.size();

      auto compare_freqs = [](double x, double y) {
        return std::abs(x - y) < 1e-3;
      };

      for (int i = 0; i < num_freqs; i++) {  // for i=0:n-1
        for (int j = i; j < num_freqs;
             j++) {  // for j = i  omega_3=-(omega_1+omega_2)
          auto omega_1 = config.frequencies[i];
          auto omega_2 = config.frequencies[j];
          auto omega_3 = omega_1 + omega_2;

          if (std::find_if(freqs_copy.begin(), freqs_copy.end(), [&](double x) {
                return compare_freqs(x, omega_3);
              }) != freqs_copy.end()) {
            continue;
          }
          if (omega_2 == 0.0)
            continue;
          freqs_copy.push_back(omega_3);
        }
      }

      config.frequencies = freqs_copy;
      std::sort(config.frequencies.begin(), config.frequencies.end());
      // only unique frequencies
      config.frequencies.erase(
          std::unique(config.frequencies.begin(), config.frequencies.end()),
          config.frequencies.end());
      return config.frequencies;
    };
    config.frequencies = set_freqs();

    response["properties"] = {};
    response["properties"]["alpha"] = base_path / "alpha.json";
    response["properties"]["beta"] = base_path / "beta.json";

    for (const auto& frequency : config.frequencies) {
      auto frequency_run_path = config.calc_path(base_path, frequency);
      auto restart_path = ResponseConfig::restart_path(frequency_run_path);
      response["calculation"].push_back(frequency_run_path);
      response["restart"].push_back(restart_path);
      response["output"].push_back(frequency_run_path / "response_base.json");
    }
    return paths;
  }

  explicit HyperPolarizabilityPathStrategy(ResponseInput input)
      : config(std::move(input)){};
  explicit HyperPolarizabilityPathStrategy(ResponseInput input,
                                           const std::string& calc_name)
      : config(calc_name, std::move(input)) {}
};
class ExcitedStatePathStrategy : public PathStrategy {

  std::string calc_name = "excited-state";
  std::string xc;
  int num_states;

  [[nodiscard]] auto calc_path(const path& root) const
      -> std::filesystem::path {
    std::string s_num_states = std::to_string(num_states);
    std::string run_name = "excited-state_" + xc + "_" + s_num_states;
    return root / std::filesystem::path(run_name);
  }

  static path restart_path(const std::filesystem::path& calc_path) {

    auto save_path = std::filesystem::path(calc_path);
    auto run_name = calc_path.filename();
    std::string save_string = "restart_" + run_name.string();
    save_path += "/";
    save_path += save_string;
    save_path += ".00000";

    return save_path;
  }

 public:
  json generateCalcPaths(const path& root) override {

    json paths;
    paths[calc_name] = {};

    auto base_path = root / calc_name;
    auto excited_path = calc_path(base_path);

    auto& excited = paths[calc_name];
    excited["calculation"] = excited_path;
    excited["output"] = excited_path / "response_base.json";
    excited["restart"] = restart_path(excited_path);

    return paths;
  }

  explicit ExcitedStatePathStrategy() = delete;
  explicit ExcitedStatePathStrategy(std::string calc_name, std::string xc,
                                    int num_states)
      : calc_name(std::move(calc_name)),
        xc(std::move(xc)),
        num_states(num_states) {}
};

class MP2PathStrategy : public PathStrategy {

  std::string calc_name = "mp2";

 public:
  json generateCalcPaths(const path& root) override {

    json paths;
    paths[calc_name] = {};
    auto& mp2 = paths[calc_name];

    auto base_path = root / calc_name;

    mp2["calculation"] = base_path;
    //mp2["output"] = excited_path / "response_base.json";
    mp2["restart"] = base_path / "mp2_restartdata.00000";

    return paths;
  }

  explicit MP2PathStrategy() = default;
  explicit MP2PathStrategy(std::string calc_name)
      : calc_name(std::move(calc_name)){};
};
#endif  // MADCHEM_PATH_MANAGER_HPP
