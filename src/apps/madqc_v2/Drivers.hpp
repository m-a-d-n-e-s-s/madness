
#pragma once

// Workflow.hpp
// Defines the Driver interface, a generic SinglePointDriver, and the Workflow class

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <memory>
#include <vector>

#include "Applications.hpp"  // Interface for SCFApplication / ResponseApplication

namespace qcapp {

/**
 * @brief Abstract base class for all drivers that encapsulate one or more Applications.
 */
class Driver {
 public:
  virtual ~Driver() = default;

  /**
   * @brief Execute the driver, writing outputs under the given directory.
   * @param workdir Base directory for this driver's outputs.
   */
  virtual void execute(const std::filesystem::path& workdir) = 0;

  /**
   * @brief Return a JSON summary of results produced by this driver.
   */
  [[nodiscard]] virtual nlohmann::json summary() const = 0;
};

/**
 * @brief Runs a single Application (e.g. SCF or Response) in its own subdirectory.
 */
class SinglePointDriver : public Driver {
 public:
  explicit SinglePointDriver(std::unique_ptr<Application> app) : app_(std::move(app)) {}

  void execute(const std::filesystem::path& workdir) override {
    // Create workdir for this application
    std::filesystem::create_directories(workdir);

    // Delegate to the Application
    app_->run(workdir);
    result_ = app_->results();
  }

  nlohmann::json summary() const override { return result_; }

 private:
  std::unique_ptr<Application> app_;
  nlohmann::json result_;
};

/**
 * @brief Orchestrates multiple drivers in sequence and writes a global output.json.
 */
class Workflow {
 public:
  Workflow() = default;

  /**
   * @brief Add a driver to the workflow.
   * @param driver Unique pointer to a Driver instance.
   */
  void addDriver(std::unique_ptr<Driver> driver) { drivers_.push_back(std::move(driver)); }

  /**
   * @brief Run all added drivers under the top-level directory, then emit output.json.
   * @param topDir Root directory for the entire workflow.
   */
  void run(const std::filesystem::path& topDir) {
    std::filesystem::create_directories(topDir);
    nlohmann::json all;
    all["tasks"] = nlohmann::json::array();

    for (size_t i = 0; i < drivers_.size(); ++i) {
      auto taskDir = topDir / ("task_" + std::to_string(i));
      drivers_[i]->execute(taskDir);
      all["tasks"].push_back(drivers_[i]->summary());
    }

    // Write out aggregate results
    std::ofstream ofs(topDir / "output.json");
    ofs << std::setw(2) << all;
    ofs.close();
  }

 private:
  std::vector<std::unique_ptr<Driver>> drivers_;
};

}  // namespace qcapp
