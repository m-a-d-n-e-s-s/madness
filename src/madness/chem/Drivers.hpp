
#pragma once

// Workflow.hpp
// Defines the Driver interface, a generic SinglePointDriver, and the Workflow
// class

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <memory>
#include <vector>

#include <madness/chem/Applications.hpp>  // Interface for SCFApplication / ResponseApplication
#include <madness/chem/SCFTargetAdapter.hpp>  // SCFTarget

namespace qcapp {

/**
 * @brief Abstract base class for all drivers that encapsulate one or more
 * Applications.
 */
class Driver {
 public:
  virtual ~Driver() = default;

  virtual void print_parameters(World& world) const =0;

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
 * @brief Runs a single Application (e.g. SCF or Response) in its own
 * subdirectory.
 */
class SinglePointDriver : public Driver {
 public:
  explicit SinglePointDriver(std::shared_ptr<Application> app) : app_(app) {}

  void print_parameters(World& world) const override {
    app_->print_parameters(world);
  }

  void execute(const std::filesystem::path& workdir) override {
    // Create workdir for this application
    std::filesystem::create_directories(workdir);

    // Delegate to the Application
    app_->run(workdir);
    result_ = app_->results();
  }

  nlohmann::json summary() const override {
    return result_;
  }

 private:
  std::shared_ptr<Application> app_;
  nlohmann::json result_;
};

class OptimizeDriver : public Driver {
 public:
  OptimizeDriver(World& w,
                 std::function<std::unique_ptr<Application>(Params)> factory,
                 Params p)
      : world_(w), factory_(std::move(factory)), params_(std::move(p)) {}

    void print_parameters(World& world) const override {
      params_.print_all();
    }

  void execute(const std::filesystem::path& workdir) override {
    // 1) make our single "opt" folder
    std::filesystem::create_directories(workdir);
    PathManager pm(workdir, "opt");
    pm.create();

    // 2) switch into it
    ScopedCWD guard(pm.dir());
    if (world_.rank() == 0)
      std::cout << "Running geometry optimization in " << pm.dir() << "\n";

    // 3) build MolOpt from Params
    auto& op = params_.get<OptimizationParameters>();
    MolOpt optimizer(op.get_maxiter(), 0.1, op.get_value_precision(),
                     op.get_geometry_tolerence(), 1e-3, 1e-5,
                     op.get_gradient_precision(), 1, op.get_algopt());
    // seed the Hessian
    optimizer.initialize_hessian(params_.get<Molecule>());

    // 4) build our target adaptor
    SCFTarget target(world_, factory_, params_);

    // 5) run the optimization
    auto mol0 = params_.get<Molecule>();
    auto mol_opt = optimizer.optimize(mol0, target);

    // 6) update params (if you plan further drivers)
    params_.set(mol_opt);

    // 7) record final results
    summary_ = {
        {"type", "optimization"}, {"final_energy", target.last_energy}
        // you could add geometry, gradient norms, etc.
    };

    // 8) optionally dump optimized geometry
    if (world_.rank() == 0) {
      auto geom_j = mol_opt.to_json();
      std::ofstream f("optimized_geometry.json");
      f << std::setw(2) << geom_j << "\n";
    }
  }

  nlohmann::json summary() const override { return summary_; }

 private:
  World& world_;
  std::function<std::unique_ptr<Application>(Params)> factory_;
  Params params_;
  nlohmann::json summary_;
};

/**
 * @brief Orchestrates multiple drivers in sequence and writes a global
 * output.json.
 */
class Workflow {
 public:
  Workflow() = default;

  /**
   * @brief Add a driver to the workflow.
   * @param driver Unique pointer to a Driver instance.
   */
  void addDriver(std::unique_ptr<Driver> driver) {
    drivers_.push_back(std::move(driver));
  }

  void print_parameters(World& world) const {
      for (const auto& d : drivers_) d->print_parameters(world);
  }

  /**
   * @brief Run all added drivers under the top-level directory, then emit
   * output.json.
   * @param topDir Root directory for the entire workflow.
   */
  void run(const std::filesystem::path& topDir) {
    std::filesystem::create_directories(topDir);
    nlohmann::json all;
    all["tasks"] = nlohmann::json::array();

    for (size_t i = 0; i < drivers_.size(); ++i) {
      auto taskDir = topDir / ("task_" + std::to_string(i));
      drivers_[i]->execute(taskDir);
      auto current_output= drivers_[i]->summary();

      // write out the current output to a file
      {
        std::ofstream ofs(taskDir / "output.json");
        ofs << std::setw(4) << current_output;
        ofs.close();
      }
      /// append current output to all
      if (current_output.is_array()) {
        for (const auto& item : current_output) {
          all["tasks"].push_back(item);
        }
      } else {
        all["tasks"].push_back(current_output);
      }


      // Write out aggregate results
      {
        std::ofstream ofs(topDir / "output.json");
        ofs << std::setw(4) << all;
        ofs.close();
      }
    }

  }

 private:
  std::vector<std::unique_ptr<Driver>> drivers_;
};

}  // namespace qcapp
