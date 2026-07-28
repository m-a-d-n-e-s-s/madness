
#pragma once

// Workflow.hpp
// Defines the Driver interface, a generic SinglePointDriver, and the Workflow
// class

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <madness/external/nlohmann_json/json.hpp>
#include <memory>
#include <stdexcept>
#include <system_error>
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
   * @param ctx     Shared StepContext threaded task-to-task; a driver reads
   *                artifacts published by upstream steps and publishes its own.
   */
  virtual void execute(const std::filesystem::path& workdir,
                       madness::StepContext& ctx) = 0;

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

  void execute(const std::filesystem::path& workdir,
               madness::StepContext& ctx) override {
    // Create workdir for this application
    std::filesystem::create_directories(workdir);

    // Read upstream artifacts (no-op unless the app overrides), run, then
    // publish this step's artifacts for downstream steps.
    app_->consume_context(ctx);
    app_->run(workdir);
    result_ = app_->results();
    app_->publish_to_context(ctx);
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

  /// Run-level provenance stamped into EVERY calc_info write (e.g. the
  /// effective restart-I/O configuration — roadmap change 5: io provenance
  /// for all tasks, not just response). Recorded here and re-applied at the
  /// top of run() (which resets the aggregate), so it is present even in the
  /// partial calc_info a failed task leaves behind.
  void set_provenance(const std::string &key, nlohmann::json value) {
    provenance_[key] = std::move(value);
  }

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
    for (const auto &kv : provenance_.items()) all_[kv.key()] = kv.value();
    all_["tasks"] = nlohmann::json::array();

    // One StepContext threaded through the whole chain: each driver reads what
    // upstream steps published and publishes its own artifacts (roadmap 1).
    madness::StepContext ctx;

    for (size_t i = 0; i < drivers_.size(); ++i) {
      auto taskDir = topDir / ("task_" + std::to_string(i));
      try {
        drivers_[i]->execute(taskDir, ctx);
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
  ///
  /// Two additional hazards handled here (raman thread brief, defect 2):
  ///  (a) an unchecked ofstream means an ENOSPC-truncated tmp still SUCCEEDS
  ///      the rename and replaces a good calc_info.json — so we verify the
  ///      stream and the rename before committing, and drop the tmp on failure;
  ///  (b) a rank-0-only filesystem failure would unwind only rank 0 while the
  ///      other ranks march into the next driver's collectives → a hang that
  ///      burns the whole allocation. We broadcast an ok/fail flag so every
  ///      rank throws together (or none does).
  void write_calc_info(const std::string &prefix) const {
    auto &world = madness::World::get_default();
    int ok = 1;
    if (world.rank() == 0) {
      const std::string outputfile = prefix + ".calc_info.json";
      const std::string tmpfile = outputfile + ".tmp";
      {
        std::ofstream ofs(tmpfile);
        ofs << std::setw(4) << all_;
        ofs.flush();
        if (!ofs) ok = 0;
      }
      if (ok) {
        std::error_code ec;
        std::filesystem::rename(tmpfile, outputfile, ec);
        if (ec) ok = 0;
      }
      if (!ok) {
        std::error_code ec;
        std::filesystem::remove(tmpfile, ec);
      }
    }
    world.gop.broadcast(&ok, 1, 0);
    if (!ok)
      throw std::runtime_error(
          "write_calc_info: failed to persist " + prefix +
          ".calc_info.json (disk full or filesystem error)");
  }

  std::vector<std::unique_ptr<Driver>> drivers_;
  nlohmann::json all_;
  nlohmann::json provenance_ = nlohmann::json::object();
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
