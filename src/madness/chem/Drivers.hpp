
#pragma once

// Workflow.hpp
// Defines the Driver interface, a generic SinglePointDriver, and the Workflow
// class

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
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
  void addDriver(std::unique_ptr<Driver> driver) { drivers_.push_back(std::move(driver)); }

  void print_parameters(World &world) const {
    for (const auto &d : drivers_)
      d->print_parameters(world);
  }

  /**
   * @brief Run all added drivers under the top-level directory, then emit
   * output.json.
   * @param topDir Root directory for the entire workflow.
   * @param outputfile Name of the output file to write the aggregated results.
   */
  void run(const std::string prefix) {
    std::filesystem::path topDir = prefix;
    std::filesystem::create_directories(topDir);
    all_ = nlohmann::json::object();
    all_["tasks"] = nlohmann::json::array();

    for (size_t i = 0; i < drivers_.size(); ++i) {
      auto taskDir = topDir / ("task_" + std::to_string(i));
      try {
        drivers_[i]->execute(taskDir);
      } catch (...) {
        // A failing task must still leave a complete, diagnosable
        // calc_info.json: record the failure as a task entry, persist, and
        // rethrow so the app-level handler reports the error (madqc.cpp).
        nlohmann::json failed;
        failed["type"] = "task_failed";
        failed["task_index"] = i;
        try {
          throw;
        } catch (const std::exception &e) {
          failed["error"] = e.what();
        } catch (...) {
          failed["error"] = "unknown exception";
        }
        all_["tasks"].push_back(failed);
        write_calc_info(prefix);
        throw;
      }
      auto current_output = drivers_[i]->summary();

      /// append current output to all
      if (current_output.is_array()) {
        for (const auto &item : current_output) {
          all_["tasks"].push_back(item);
        }
      } else {
        all_["tasks"].push_back(current_output);
      }

      write_calc_info(prefix);
    }
  }

  /// Aggregated results of all drivers (the calc_info JSON); valid after run().
  const nlohmann::json &results() const { return all_; }

private:
  /// Rank-0-only, ATOMIC (tmp+rename) aggregate write. Every rank used to
  /// stream the file directly — N concurrent writers to one path on a shared
  /// FS, and a crash mid-write truncated it.
  void write_calc_info(const std::string &prefix) const {
    if (madness::World::get_default().rank() != 0) return;
    const std::string outputfile = prefix + ".calc_info.json";
    const std::string tmpfile = outputfile + ".tmp";
    {
      std::ofstream ofs(tmpfile);
      ofs << std::setw(4) << all_;
    }
    std::filesystem::rename(tmpfile, outputfile);
  }

  std::vector<std::unique_ptr<Driver>> drivers_;
  nlohmann::json all_;
};

} // namespace qcapp
// class OptimizeDriver : public Driver {
// public:
//   OptimizeDriver(World &w, std::function<std::unique_ptr<Application>(Params)> factory, Params p)
//       : world_(w), factory_(std::move(factory)), params_(std::move(p)) {}

//   void print_parameters(World &world) const override {
//     if (world.rank() == 0) {
//       params_.print_all();
//     }
//   }

//   void execute(const std::filesystem::path &workdir) override {
//     // 1) make our single "opt" folder
//     std::filesystem::create_directories(workdir);
//     PathManager pm(workdir, "opt");
//     pm.create();

//     // 2) switch into it
//     ScopedCWD guard(pm.dir());
//     if (world_.rank() == 0)
//       std::cout << "Running geometry optimization in " << pm.dir() << "\n";

//     // 3) build MolOpt from Params
//     auto &op = params_.get<OptimizationParameters>();
//     MolOpt optimizer(op.get_maxiter(), 0.1, op.get_value_precision(), op.get_geometry_tolerence(), 1e-3, 1e-5,
//                      op.get_gradient_precision(), 1, op.get_algopt());
//     // seed the Hessian
//     optimizer.initialize_hessian(params_.get<Molecule>());

//     // 4) build our target adaptor
//     SCFTarget target(world_, factory_, params_);

//     // 5) run the optimization
//     auto mol0 = params_.get<Molecule>();
//     Molecule optimized_mol = optimizer.optimize(mol0, target);

//     OptimizationResults results;
//     results.final_geometry = optimized_mol;
//     results.final_energy = target.last_energy;
//     // 6) update params (if you plan further drivers)
//     params_.set(optimized_mol);

//     // 7) record final results
//     summary_ = {
//         {"type", "optimization"}, {"final_energy", target.last_energy}
//         // you could add geometry, gradient norms, etc.
//     };

//     // 8) optionally dump optimized geometry
//     if (world_.rank() == 0) {
//       auto geom_j = optimized_mol.to_json();
//       std::ofstream f("optimized_geometry.json");
//       f << std::setw(2) << geom_j << "\n";
//     }
//   }

//   [[nodiscard]] nlohmann::json summary() const override { return summary_; }

// private:
//   World &world_;
//   std::function<std::unique_ptr<Application>(Params)> factory_;
//   Params params_;
//   nlohmann::json summary_;
// };
