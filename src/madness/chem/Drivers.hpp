
#pragma once

// Workflow.hpp
// Defines the Driver interface, a generic SinglePointDriver, and the Workflow
// class

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <nlohmann/json.hpp>
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

/**
 * @brief Geometry optimization as its own task (madqc ARCHITECTURE_ROADMAP change 2).
 *
 * Drives MolOpt over the reference engine the Library policy supplies, so the same
 * code optimizes on moldft (`Calc = SCF`) and nemo (`Calc = Nemo`). Selected by
 * `--optimize` together with `--wf=<scf|nemo>`, which names the reference method.
 *
 * A Driver rather than an Application, deliberately: numerical gradients are to be
 * computed from displaced sub-runs, each with its own directory and calc_info (see
 * GeometryTarget), and owning sub-runs is a Driver's job. The same seam serves
 * roadmap changes 3 and 4.
 *
 * Not restartable as a step: the optimizer keeps no checkpoint of its own (hessian
 * and step history), so an interrupted optimization restarts from the input
 * geometry. The SCF underneath still restarts per geometry as usual.
 */
template <typename Library> class OptimizeDriver : public Driver {
public:
  using Calc = typename Library::Calc;

  OptimizeDriver(World &world, const Params &params)
      : world_(world), params_(params) {}

  void print_parameters(World &world) const override {
    if (world.rank() != 0)
      return;
    params_.get<OptimizationParameters>().print(OptimizationParameters::tag,
                                                "end");
    if constexpr (std::is_same_v<Calc, madness::SCF>) {
      params_.get<CalculationParameters>().print(CalculationParameters::tag,
                                                 "end");
    } else {
      params_.get<CalculationParameters>().print(CalculationParameters::tag);
      params_.get<madness::Nemo::NemoCalculationParameters>().print();
      madness::print("end");
    }
  }

  void execute(const std::filesystem::path &workdir,
               madness::StepContext &ctx) override {
    // An upstream step may have moved the molecule; honour it before the engine is
    // built, since construction freezes the geometry.
    if (ctx.molecule && ctx.molecule->natom() > 0)
      params_.get<madness::Molecule>() = *ctx.molecule;

    madness::PathManager pm(workdir, Library::label());
    pm.create();
    world_.gop.fence();
    {
      madness::ScopedCWD scwd(pm.dir());
      if (world_.rank() == 0)
        madness::print("Running geometry optimization on", Library::label(), "in",
                       pm.dir().string());

      // Tie the geometry thresholds to the accuracy of the wavefunction, the same
      // Dalton-style rules the in-SCF gopt path uses -- but as *derived* values, so
      // anything the deck sets in the `optimization` group wins.
      auto &op = params_.get<OptimizationParameters>();
      const double wf_thresh =
          params_.get<CalculationParameters>().protocol().back();
      op.set_derived_value("etol", std::max(1.0e-7, 2.0 * wf_thresh));
      op.set_derived_value("gtol", std::max(1.0e-5, 2.0 * wf_thresh));
      op.set_derived_value("xtol", std::max(1.0e-5, 2.0 * wf_thresh));
      op.set_derived_value("value_precision", std::max(1.0e-7, 2.0 * wf_thresh));
      op.set_derived_value("gradient_precision", std::max(1.0e-6, wf_thresh));

      MADNESS_CHECK_THROW(
          !op.get_initial_hessian(),
          "optimization: initial_hessian is not implemented for the optimizer -- "
          "MolOpt starts from a mass-weighted diagonal guess. Remove the key "
          "rather than have it silently ignored");

      auto engine = lib_.calc(world_, params_);
      engine->work_dir = pm.dir();

      if constexpr (std::is_same_v<Calc, madness::SCF>) {
        // The same preparation moldft_lib::run does before handing the engine to
        // MolOpt. It is not optional: MolecularEnergy::value only calls
        // set_protocol when FunctionDefaults' thresh differs from protocol[0], so
        // when they already agree nothing would build the nuclear potential / data
        // map and the first geometry segfaults in the initial guess. Nemo needs no
        // equivalent -- Nemo::value sets its own protocol on every call.
        // `template` disambiguator: Calc is a template parameter here.
        if (world_.size() > 1) {
          engine->template set_protocol<3>(world_, 1e-4);
          engine->make_nuclear_potential(world_);
          engine->initial_load_bal(world_);
        }
        engine->template set_protocol<3>(world_,
                                         engine->param.protocol()[0]);
      } else {
        // Nemo carries work_dir both on itself and on its inner SCF; downstream
        // consumers read one or the other, so set both.
        engine->get_calc()->work_dir = pm.dir();
      }

      // The engine as an OptimizationTargetInterface: MolecularEnergy wraps an SCF,
      // Nemo is one itself. `scf_target` must outlive the optimization -- it holds
      // an SCF&.
      std::unique_ptr<madness::MolecularEnergy> scf_target;
      madness::OptimizationTargetInterface *engine_target = nullptr;
      if constexpr (std::is_same_v<Calc, madness::SCF>) {
        scf_target = std::make_unique<madness::MolecularEnergy>(world_, *engine);
        engine_target = scf_target.get();
      } else {
        engine_target = engine.get();
      }
      MADNESS_CHECK_THROW(engine_target->provides_gradient(),
                          "the reference engine provides no analytic gradient; "
                          "numerical gradients from displaced sub-runs are not "
                          "implemented yet (see GeometryTarget)");

      // THE SEAM: swap in a displaced-sub-run target here to get numerical
      // gradients; everything below is unchanged by that choice.
      madness::AnalyticTarget target(*engine_target);

      madness::MolOpt opt(op.get_maxiter(), op.get_maxstep(), op.get_etol(),
                          op.get_gtol(), op.get_xtol(), op.get_value_precision(),
                          op.get_gradient_precision(),
                          (world_.rank() == 0) ? 1 : 0, op.get_algopt());

      madness::OptimizationResults opt_res;
      if constexpr (std::is_same_v<Calc, madness::SCF>) {
        opt_res = opt.optimize_app(engine->molecule, target);
      } else {
        opt_res = opt.optimize_app(engine->molecule(), target);
      }

      // Leave the engine AT the optimized geometry and pick up the final energy and
      // gradient there (mirrors the in-SCF gopt path).
      double energy = 0.0;
      madness::Tensor<double> gradient;
      target.energy_and_gradient(opt_res.final_geometry, energy, gradient);
      opt_res.final_energy = energy;

      madness::PropertyResults prop_res;
      prop_res.energy = energy;
      prop_res.gradient = gradient;

      summary_["model"] = "optimize";
      summary_["optimization_results"] = opt_res.to_json();
      summary_["molecule"] = opt_res.final_geometry.to_json();
      summary_["properties"] = prop_res.to_json();
      summary_["metadata"] = {{"mpi_size", world_.size()},
                              {"method", Library::label()}};

      if (world_.rank() == 0) {
        const std::string geomfile =
            params_.get<CalculationParameters>().prefix() + "_opt.xyz";
        std::ofstream ofs(geomfile);
        opt_res.final_geometry.print(ofs);
        ofs.close();
        madness::print("optimized geometry written to", geomfile);
        opt_res.final_geometry.print();
      }
      final_geometry_ = opt_res.final_geometry;
    }

    // Hand the optimized geometry to whatever runs next -- the point of making this
    // a first-class step.
    if (final_geometry_.natom() > 0)
      ctx.molecule = final_geometry_;
    if constexpr (std::is_same_v<Calc, madness::SCF>)
      ctx.reference = lib_.calc(world_, params_);
    try {
      const auto &cp = params_.get<CalculationParameters>();
      ctx.archives["restartdata"] = pm.dir() / (cp.prefix() + ".restartdata");
    } catch (...) {
      // best effort, as in SCFApplication
    }
  }

  nlohmann::json summary() const override { return summary_; }

private:
  World &world_;
  Params params_;
  Library lib_;
  nlohmann::json summary_;
  madness::Molecule final_geometry_;
};

} // namespace qcapp
