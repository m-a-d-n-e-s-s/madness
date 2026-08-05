#pragma once

#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/PathManager.hpp>
#include <madness/chem/Results.h>
#include <madness/chem/SCFTargetAdapter.hpp>
#include <madness/chem/molopt.h>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <system_error>
#include <type_traits>

namespace madness {
enum class NextAction { Ok, ReloadOnly, Restart, Redo };

class SCF; // forward decl for StepContext::reference

/// Typed artifacts handed from one workflow step to the next, threaded through
/// Workflow::run and Driver::execute (madqc ARCHITECTURE_ROADMAP change 1).
///
/// This replaces the cwd/mad.in side channels for downstream geometry + archive
/// discovery. The build-time shared_ptr capture of the ground-state reference
/// remains a supported compatibility path — a producer that can expose a live
/// SCF publishes it here (`reference`), and a consumer prefers it when present.
///
/// Producers set fields in publish_to_context(); consumers read them in
/// consume_context(). Everything is optional: an unset field means "no upstream
/// value — fall back to the build-time / input value".
struct StepContext {
  /// The molecule the next step should use (possibly optimized/displaced).
  std::optional<Molecule> molecule;
  /// Live ground-state reference engine, if an upstream SCF step exposed one.
  std::shared_ptr<SCF> reference;
  /// Named archive/output paths (absolute), e.g. "restartdata" -> path.
  std::map<std::string, std::filesystem::path> archives;
  /// Free-form JSON for artifacts not yet first-class.
  nlohmann::json blob = nlohmann::json::object();
};

// Scoped CWD: changes the current directory to the given one, and restores when
// the object goes out of scope
struct ScopedCWD {
  std::filesystem::path old_cwd;

  explicit ScopedCWD(const std::filesystem::path &new_dir) {
    old_cwd = std::filesystem::current_path();
    std::filesystem::current_path(new_dir);
  }

  ~ScopedCWD() { std::filesystem::current_path(old_cwd); }
};

class Application {
public:
  explicit Application(const Params &p) : params_(p) {}

  virtual ~Application() = default;

  // run: write all outputs under the given directory
  virtual void run(const std::filesystem::path &workdir) = 0;

  /// Consume artifacts published by an upstream step. Called BEFORE run().
  /// Default no-op; override to read the shared StepContext (e.g. a response
  /// step preferring the upstream ground-state reference / geometry).
  virtual void consume_context(const StepContext & /*ctx*/) {}

  /// Publish this step's typed outputs into the shared StepContext for
  /// downstream steps. Called AFTER run(). Default no-op.
  virtual void publish_to_context(StepContext & /*ctx*/) {}

  // optional hook to return a JSON fragment of this app's main results
  [[nodiscard]] virtual nlohmann::json results() const = 0;

  virtual void print_parameters(World &world) const = 0;

  /// check if this calculation has a json with results
  [[nodiscard]] virtual bool has_results(const std::string &filename) const {
    // check if the results file exists
    // return std::filesystem::exists(workdir_ / filename);
    return std::filesystem::exists(filename);
  }

  [[nodiscard]] virtual bool verify_molecule(const nlohmann::json &j) const {
    // check if some key parameters of the calculation match:
    // molecule, box size, nmo_alpha, nmo_beta
    Molecule mol1 = params_.get<Molecule>();
    Molecule mol2;
    mol2.from_json(j["molecule"]);
    if (not(mol1 == mol2)) {
      print("molecule mismatch");
      mol1.print();
      mol2.print();
      return false;
    }
    return true;
  }

  /// read the results from a json file
  [[nodiscard]] virtual nlohmann::json
  read_results(const std::string &filename) const {
    if (has_results(filename)) {
      std::ifstream ifs(filename);
      nlohmann::json j;
      ifs >> j;
      ifs.close();
      // if (not verify_molecule(j))
      // {
      //   std::string msg =
      //       "Results file " + filename + " does not match the parameters of
      //       the calculation";
      //   print(msg);
      //   return nlohmann::json(); // return empty json
      // }
      return j;
    } else {
      std::string msg = "Results file " + filename + " does not exist in " +
                        std::filesystem::current_path().string();
      MADNESS_EXCEPTION(msg.c_str(), 1);
    }
    return nlohmann::json();
  }

protected:
  Params params_;
  nlohmann::json results_;
};

template <typename Library> class SCFApplication : public Application {
private:
public:
  using Calc = typename Library::Calc;

  explicit SCFApplication(World &w, const Params &p)
      : Application(p), world_(w) {}

  // Give downstream steps the live calc
  std::shared_ptr<Calc> calc() { return lib_.calc(world_, params_); }
  void set_calc_workdir(const std::filesystem::path &workdir) {
    calc()->work_dir = workdir;
  }

  // print parameters
  /// Print the *effective* parameters of this step (user-defined, derived and
  /// default values, as annotated by QCCalculationParametersBase::print), not
  /// the static template of available keys — Library::print_parameters() is the
  /// latter and is what `--print_parameters=<group>` is for.
  void print_parameters(World &world) const override {
    if (world.rank() != 0)
      return;
    if constexpr (std::is_same_v<Calc, SCF>) {
      params_.get<CalculationParameters>().print(
          CalculationParameters::tag, "end");
    } else {
      // Nemo-based engines carry an additional nemo block inside dft
      params_.get<CalculationParameters>().print(CalculationParameters::tag);
      params_.get<Nemo::NemoCalculationParameters>().print();
      print("end");
    }
  }

  // sets the calc working directory and runs the calculation
  void run(const std::filesystem::path &workdir) override {
    // 1) set up a namedspaced directory for this run
    std::string label = Library::label();
    PathManager pm(workdir, label);
    pm.create();
    {
      world_.gop.fence();
      ScopedCWD scwd(pm.dir());
      if (world_.rank() == 0) {
        std::cout << "Running SCF in " << pm.dir() << std::endl;
      }
      // 2) define the "checkpoint" file
      auto ckpt = label + ".calc_info.json";
      SCFResultsTuple empty_results;
      nlohmann::json j;
      NextAction action;
      if (has_results(ckpt)) {
        try {
          // Parse INSIDE the try: a truncated/corrupt checkpoint (e.g. a
          // process killed mid-write) makes read_results throw, and that
          // throw must degrade to Redo rather than escaping past this catch
          // and bricking the restart.
          j = read_results(ckpt); // which results are we reading
          auto &[scf_r, properties, convergence, optr] = scf_results;
          scf_r.from_json(j["scf"]);
          properties.from_json(j["properties"]);
          convergence.from_json(j["convergence"]);
          action = lib_.valid(world_, scf_results, params_);

        } catch (...) {
          print("Failed to parse checkpoint file: ", ckpt);
          j = nlohmann::json(); // drop any partially-parsed data
          scf_results = empty_results;
          action = madness::NextAction::Redo;
        }
      } else {
        scf_results = empty_results;
        action = madness::NextAction::Redo;
      }

      // Guard against reusing a checkpoint computed for a DIFFERENT molecule.
      if (action != madness::NextAction::Redo &&
          !checkpoint_geometry_matches(j)) {
        if (world_.rank() == 0)
          print("WARNING: checkpoint geometry does not match the requested "
                "molecule; ignoring checkpoint and recomputing.");
        scf_results = empty_results;
        action = madness::NextAction::Redo;
      }
      // If we are restarting from an existing archive, the SCF engine must be
      // told to load its MOs from <prefix>.restartdata. This flag has to be set
      // on params_ BEFORE the engine is constructed — set_calc_workdir() below
      // builds the SCF (via calc()->initialize_(), which freezes params into
      // mad.in), so a flag set afterwards (as the old moldft_lib::run did, on a
      // discarded local copy) never reaches the constructor and "Restart"
      // silently recomputed from scratch. (raman thread brief, defect 3)
      if (action == madness::NextAction::Restart) {
        params_.get<CalculationParameters>().set_user_defined_value("restart",
                                                                    true);
        if (world_.rank() == 0)
          print("Restart requested: loading MOs from restartdata archive");
      }
      world_.gop.fence();
      set_calc_workdir(pm.dir());
      auto params_copy = params_;

      // if okay and optimize we have to set Molecule to new geometry
      //
      if (world_.rank() == 0) {
        print("Molecule from params:");
        params_.get<Molecule>().print();
      }
      if (action == madness::NextAction::Ok ||
          action == madness::NextAction::ReloadOnly ||
          action == madness::NextAction::Restart) {
        if (world_.rank() == 0) {
          print("SCF results are valid, no need to rerun");
          print("Ensure we are running on the correct molecule geometry");
        }

        auto &mol = params_.get<Molecule>();
        mol.from_json(j["scf"]["molecule"]);
      }
      if (world_.rank() == 0) {
        print("Molecule from scf results :");
        params_.get<Molecule>().print();
      }

      world_.gop.fence();

      if (world_.rank() == 0)
        print("Next action is ", static_cast<int>(action),
              " (0=Ok,1=ReloadOnly,2=Restart,3=Redo)");

      if (action == madness::NextAction::Restart ||
          action == madness::NextAction::Redo) {
        // Restart vs Redo is now carried by the 'restart' flag set on params_
        // above (Restart => load restartdata; Redo => fresh). Pass the real
        // action through rather than a hardcoded Restart.
        scf_results = lib_.run(world_, params_, action);
      } else {
        lib_.calc(world_, params_); // just set up the calc without running
      }

      // // Need work (Restart or Redo) — both call run()
      // scf_results = lib_.run(world_, params_);

      results_["scf"] = std::get<0>(scf_results).to_json();
      results_["properties"] = std::get<1>(scf_results).to_json();
      results_["convergence"] = std::get<2>(scf_results).to_json();
      results_["molecule"] = std::get<0>(scf_results).scf_molecule.to_json();
      results_["optimization_results"] = std::get<3>(scf_results).to_json();
      // Backward-compatible top-level fields expected by existing scripted tests.
      // Keep these in sync with the nested "scf/properties/convergence" schema.
      results_["model"] = "scf";
      results_["scf_total_energy"] = results_["scf"]["scf_total_energy"];
      results_["scf_eigenvalues_a"] = results_["scf"]["scf_eigenvalues_a"];
      results_["scf_fock_a"] = results_["scf"]["scf_fock_a"];
      // Open-shell: mirror the beta channel too (review MED — SCFResults emits
      // scf_eigenvalues_b/scf_fock_b in the nested object, but only the alpha
      // channel was surfaced at top level, so the .out summary never showed
      // beta). Guarded on presence: closed-shell runs omit them.
      if (results_["scf"].contains("scf_eigenvalues_b"))
        results_["scf_eigenvalues_b"] = results_["scf"]["scf_eigenvalues_b"];
      if (results_["scf"].contains("scf_fock_b"))
        results_["scf_fock_b"] = results_["scf"]["scf_fock_b"];
      results_["convergence_info"] = results_["convergence"];
      results_["metadata"] = {{"mpi_size", world_.size()}};

      // write the checkpoint file atomically (tmp write + rename) so a crash
      // or ENOSPC mid-write cannot truncate a previously-good checkpoint and
      // brick the next restart. See defect-1 in the raman thread brief.
      if (world_.rank() == 0) {
        const std::string tmp = ckpt + ".tmp";
        bool ok = true;
        {
          std::ofstream ofs(tmp);
          ofs << results_.dump(4);
          ofs.flush();
          ok = static_cast<bool>(ofs);
        }
        if (ok) {
          std::error_code ec;
          std::filesystem::rename(tmp, ckpt, ec);
          if (ec) ok = false;
        }
        if (ok) {
          print("Written checkpoint file: ", ckpt);
        } else {
          std::error_code ec;
          std::filesystem::remove(tmp, ec);
          print("ERROR: failed to write checkpoint file (disk full?): ", ckpt);
        }
      }
    }
  }

  // std::shared_ptr<SCFApplicationT> scf_app =
  // std::dynamic_pointer_cast<SCFApplicationT>(reference_.shared_from_this());

  /// Publish this SCF step's typed outputs for downstream steps: the live
  /// engine (only when it is an SCF — nemo's Calc is Nemo and stays on the
  /// build-time capture path), the converged molecule, and the restartdata
  /// archive path. See StepContext / ARCHITECTURE_ROADMAP change 1.
  void publish_to_context(StepContext &ctx) override {
    if constexpr (std::is_same_v<Calc, SCF>) {
      ctx.reference = calc(); // shared_ptr<SCF>
    }
    if (results_.contains("molecule")) {
      try {
        Molecule m;
        m.from_json(results_["molecule"]);
        ctx.molecule = std::move(m);
      } catch (...) {
        // leave ctx.molecule unset -> downstream falls back to input geometry
      }
    }
    // Record the restartdata archive location (absolute) so a downstream step
    // can restart from this ground state without relying on the cwd.
    try {
      const auto &cp = params_.get<CalculationParameters>();
      const std::filesystem::path dir =
          calc()->work_dir.empty() ? std::filesystem::current_path()
                                    : std::filesystem::path(calc()->work_dir);
      ctx.archives["restartdata"] = dir / (cp.prefix() + ".restartdata");
    } catch (...) {
      // best-effort; absence just means downstream restart discovery falls back
    }
  }

  nlohmann::json results() const override { return results_; }

private:
  // Returns true iff the checkpoint's stored geometry matches the requested
  // molecule's nuclear framework. Compares only atoms (element + position via
  // Atom::operator==), NOT derived quantities (rcut/field/pointgroup/
  // parameters), which can differ on a JSON round-trip and would otherwise
  // force spurious recomputes.
  bool checkpoint_geometry_matches(const nlohmann::json &j) const {
    const nlohmann::json *molj = nullptr;
    if (j.contains("molecule"))
      molj = &j.at("molecule");
    else if (j.contains("scf") && j["scf"].contains("molecule"))
      molj = &j["scf"].at("molecule");
    if (!molj)
      return true; // nothing to compare against -> don't block reuse
    Molecule ckpt_mol;
    try {
      ckpt_mol.from_json(*molj);
    } catch (...) {
      return true; // can't parse -> let other validation decide
    }
    const Molecule &want = params_.get<Molecule>();
    if (ckpt_mol.natom() != want.natom())
      return false;
    for (unsigned int i = 0; i < want.natom(); ++i)
      if (!(want.get_atom(i) == ckpt_mol.get_atom(i)))
        return false;
    return true;
  }

  World &world_;
  Library lib_; // owns shared_ptr<Engine>
  SCFResultsTuple scf_results;
};

/// First-class geometry optimization step (madqc ARCHITECTURE_ROADMAP change 2).
///
/// Drives MolOpt over the reference engine the Library policy supplies, so the
/// same code optimizes with moldft (`Calc = SCF`) and nemo (`Calc = Nemo`); the
/// only engine-specific parts are which target adapter wraps it and how its
/// molecule is reached. Selected by `optimization { method moldft|nemo }`.
///
/// This is distinct from the in-SCF `dft gopt` path (moldft_lib::run), which
/// optimizes inside a single SCF task. Both drive the same MolOpt; the difference
/// is that this one is its own task, honours the `optimization` parameter group,
/// and publishes the optimized geometry into the StepContext for a later step.
///
/// Not yet a restartable step: an interrupted optimization restarts from the input
/// geometry, since no checkpoint of the optimizer state (hessian + step history) is
/// written. The SCF underneath still restarts per geometry as usual.
template <typename Library> class OptimizeApplication : public Application {
public:
  using Calc = typename Library::Calc;

  explicit OptimizeApplication(World &w, const Params &p)
      : Application(p), world_(w) {}

  std::shared_ptr<Calc> calc() { return lib_.calc(world_, params_); }

  void print_parameters(World &world) const override {
    if (world.rank() != 0)
      return;
    params_.get<OptimizationParameters>().print(OptimizationParameters::tag,
                                                "end");
    if constexpr (std::is_same_v<Calc, SCF>) {
      params_.get<CalculationParameters>().print(CalculationParameters::tag,
                                                 "end");
    } else {
      params_.get<CalculationParameters>().print(CalculationParameters::tag);
      params_.get<Nemo::NemoCalculationParameters>().print();
      print("end");
    }
  }

  void run(const std::filesystem::path &workdir) override {
    PathManager pm(workdir, Library::label());
    pm.create();
    world_.gop.fence();
    {
      ScopedCWD scwd(pm.dir());
      if (world_.rank() == 0)
        print("Running geometry optimization on", Library::label(), "in",
              pm.dir().string());

      // Tie the geometry thresholds to the accuracy of the wavefunction, the same
      // Dalton-style rules the in-SCF gopt path uses -- but as *derived* values, so
      // anything the deck sets in the `optimization` group wins. This is what makes
      // the group's knobs effective; they were previously read by nobody.
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
          "optimization: initial_hessian is not implemented for the `optimize` "
          "workflow -- MolOpt starts from a mass-weighted diagonal guess. Remove "
          "the key rather than have it silently ignored");

      auto engine = lib_.calc(world_, params_);
      engine->work_dir = pm.dir();

      if constexpr (std::is_same_v<Calc, SCF>) {
        // The same preparation moldft_lib::run does before handing the engine to
        // MolOpt. It is not optional: MolecularEnergy::value only calls
        // set_protocol when FunctionDefaults' thresh differs from protocol[0], so
        // when they already agree (the common case) nothing would build the
        // nuclear potential / data map and the first geometry segfaults in the
        // initial guess. Nemo needs no equivalent -- Nemo::value sets its own
        // protocol on every call.
        // `template` disambiguator: Calc is a template parameter here, so
        // set_protocol<3> is a dependent name.
        if (world_.size() > 1) {
          engine->template set_protocol<3>(world_, 1e-4);
          engine->make_nuclear_potential(world_);
          engine->initial_load_bal(world_);
        }
        engine->template set_protocol<3>(world_, engine->param.protocol()[0]);
      } else {
        // Nemo carries work_dir both on itself and on its inner SCF; downstream
        // consumers read one or the other (CC2 reads Nemo::work_dir, the response
        // adapter reads SCF::work_dir), so set both.
        engine->get_calc()->work_dir = pm.dir();
      }

      MolOpt opt(op.get_maxiter(), op.get_maxstep(), op.get_etol(), op.get_gtol(),
                 op.get_xtol(), op.get_value_precision(),
                 op.get_gradient_precision(), (world_.rank() == 0) ? 1 : 0,
                 op.get_algopt());

      OptimizationResults opt_res;
      double energy = 0.0;
      Tensor<double> gradient;
      if constexpr (std::is_same_v<Calc, SCF>) {
        MolecularEnergy target(world_, *engine);
        opt_res = opt.optimize_app(engine->molecule, target);
        // Leave the engine AT the optimized geometry and pick up the final energy
        // and gradient there (mirrors the in-SCF gopt path).
        target.energy_and_gradient(opt_res.final_geometry, energy, gradient);
      } else {
        EngineOptTarget<Calc> target(world_, *engine);
        opt_res = opt.optimize_app(engine->molecule(), target);
        target.energy_and_gradient(opt_res.final_geometry, energy, gradient);
      }
      opt_res.final_energy = energy;

      PropertyResults prop_res;
      prop_res.energy = energy;
      prop_res.gradient = gradient;

      results_["model"] = "optimize";
      results_["optimization_results"] = opt_res.to_json();
      results_["molecule"] = opt_res.final_geometry.to_json();
      results_["properties"] = prop_res.to_json();
      results_["metadata"] = {{"mpi_size", world_.size()},
                              {"method", Library::label()}};

      if (world_.rank() == 0) {
        const std::string geomfile =
            params_.get<CalculationParameters>().prefix() + "_opt.xyz";
        std::ofstream ofs(geomfile);
        opt_res.final_geometry.print(ofs);
        ofs.close();
        print("optimized geometry written to", geomfile);
        opt_res.final_geometry.print();
      }
      final_geometry_ = opt_res.final_geometry;
    }
  }

  /// Hand the optimized geometry (and, for moldft, the live engine) to whatever
  /// runs next -- the point of making this a first-class step.
  void publish_to_context(StepContext &ctx) override {
    if (final_geometry_.natom() > 0)
      ctx.molecule = final_geometry_;
    if constexpr (std::is_same_v<Calc, SCF>)
      ctx.reference = calc();
    try {
      const auto &cp = params_.get<CalculationParameters>();
      const std::filesystem::path dir = calc()->work_dir.empty()
                                            ? std::filesystem::current_path()
                                            : calc()->work_dir;
      ctx.archives["restartdata"] = dir / (cp.prefix() + ".restartdata");
    } catch (...) {
      // best effort, as in SCFApplication
    }
  }

  nlohmann::json results() const override { return results_; }

private:
  World &world_;
  Library lib_;
  Molecule final_geometry_;
};

/**
 * @brief Wrapper application to run the molresponse workflow
 *        via the molresponse_lib::run_response function.
 */
template <typename Library> class ResponseApplication : public Application {
public:
  /**
   * @param world   MADNESS world communicator
   * @param params  Unified Params containing ResponseParameters & Molecule
   * @param ref_dir   Directory of precomputed ground-state (SCF) outputs
   */
  ResponseApplication(World &world, Params params,
                      std::shared_ptr<SCF> reference)
      : Application(std::move(params)),
        world_(world),
        reference_(std::move(reference)) {}

  // print parameters
  void print_parameters(World &world) const override {
    if (world.rank() == 0)
      params_.get<ResponseParameters>().print(ResponseParameters::tag, "end");
  }

  /// Prefer the ground-state reference published by the upstream SCF step over
  /// the one captured at build time. Also adopt an upstream (optimized/
  /// displaced) geometry so response runs AT the geometry the chain computed.
  /// (ARCHITECTURE_ROADMAP change 1 acceptance.)
  void consume_context(const StepContext &ctx) override {
    if (ctx.reference)
      reference_ = ctx.reference;
    if (ctx.molecule)
      params_.get<Molecule>() = *ctx.molecule;
  }

  /**
   * @brief Execute response + property workflow, writing into workdir/response
   */
  void run(const std::filesystem::path &workdir) override {
    // create a namespaced subdirectory for response outputs
    PathManager pm(workdir, Library::label());
    pm.create();
    {
      ScopedCWD scwd(pm.dir());

      auto res = Library::run_response(world_, params_, reference_, pm.dir());

      metadata_ = std::move(res.metadata);
      properties_["response_properties"] = std::move(res.properties);
      properties_["vibrational_analysis"] = std::move(res.vibrational_analysis);
      properties_["raman_spectra"] = std::move(res.raman_spectra);
    }
  }

  /**
   * @brief Return a JSON fragment summarizing results
   */
  [[nodiscard]] nlohmann::json results() const override {
    return {{"type", "response"},
            {"metadata", metadata_},
            {"properties", properties_}};
  }

private:
  World &world_;
  nlohmann::json metadata_;
  nlohmann::json properties_;
  std::optional<nlohmann::json> vibrational_analysis_;
  std::shared_ptr<SCF> reference_;
};

class CC2Application : public Application, public CC2 {
public:
  explicit CC2Application(World &w, const Params &p,
                          const std::shared_ptr<Nemo> &reference)
      : Application(p),
        CC2(w, p.get<CCParameters>(), p.get<TDHFParameters>(), reference),
        world_(w), reference_(reference) {}

  // print_parameters
  void print_parameters(World &world) const override {
    if (world.rank() == 0)
      params_.get<CCParameters>().print(CCParameters::tag, "end");
  }

  void run(const std::filesystem::path &workdir) override {
    // 1) set up a namedspaced directory for this run
    std::string label = "cc2";
    PathManager pm(workdir, label);
    pm.create();
    world_.gop.fence();
    {
      ScopedCWD scwd(pm.dir());
      if (world_.rank() == 0) {
        std::cout << "Running CC2 in " << pm.dir() << std::endl;
      }

      // 2) define the "checkpoint" file
      auto ckpt = label + "_results.json";
      print("cc checkpoint file", ckpt);
      if (std::filesystem::exists(ckpt)) {
        if (world_.rank() == 0) {
          std::cout << "Found checkpoint file: " << ckpt << std::endl;
        }
        // read the checkpoint file
        std::ifstream ifs(ckpt);
        ifs >> results_;
        ifs.close();

        // bool ok = true;
        // bool needEnergy = true;
        // if (needEnergy && !results_.contains("energy"))
        //   ok = false;
      }

      auto rel = std::filesystem::relative(reference_->work_dir, pm.dir());
      if (world_.rank() == 0) {
        std::cout << "Running cc2 calculation in: " << pm.dir() << std::endl;
        std::cout << "Ground state archive: " << reference_->work_dir
                  << std::endl;
        std::cout << "Relative path: " << rel << std::endl;
      }

      results_ = this->solve();
    }
  }

  nlohmann::json results() const override { return results_; }

private:
  World &world_;
  const std::shared_ptr<Nemo> reference_;
};

class TDHFApplication : public Application, public TDHF {
public:
  explicit TDHFApplication(World &w, const Params &p,
                           const std::shared_ptr<Nemo> &reference)
      : Application(p), TDHF(w, p.get<TDHFParameters>(), reference), world_(w),
        reference_(reference) {}

  // print_parameters
  void print_parameters(World &world) const override {
    if (world.rank() == 0)
      params_.get<TDHFParameters>().print(TDHFParameters::tag, "end");
  }

  void run(const std::filesystem::path &workdir) override {
    // 1) set up a namedspaced directory for this run
    PathManager pm(workdir, "tdhf");
    pm.create();
    world_.gop.fence();
    {
      ScopedCWD scwd(pm.dir());
      if (world_.rank() == 0) {
        std::cout << "Running CIS in " << pm.dir() << std::endl;
      }

      // we could dump params_ to JSON and pass as argv if desired…
      try {
        const double time_scf_start = wall_time();
        this->prepare_calculation();
        const double time_scf_end = wall_time();
        if (world_.rank() == 0)
          printf(" at time %.1f\n", wall_time());

        const double time_cis_start = wall_time();
        std::vector<CC_vecfunction> roots = this->solve_cis();
        const double time_cis_end = wall_time();
        if (world_.rank() == 0)
          printf(" at time %.1f\n", wall_time());

        if (world_.rank() == 0) {
          std::cout << std::setfill(' ');
          std::cout << "\n\n\n";
          std::cout << "--------------------------------------------------\n";
          std::cout << "MRA-CIS ended \n";
          std::cout << "--------------------------------------------------\n";
          std::cout << std::setw(25) << "time scf" << " = "
                    << time_scf_end - time_scf_start << "\n";
          std::cout << std::setw(25) << "time cis" << " = "
                    << time_cis_end - time_cis_start << "\n";
          std::cout << "--------------------------------------------------\n";
        }
        auto j = this->analyze(roots);
        // funnel through CISResults to make sure we have the right format
        CISResults results(j);
        results_ = results.to_json();
      } catch (std::exception &e) {
        // Do not silently swallow: record the failure so the emitted
        // calc_info.json reflects it instead of looking like a clean run.
        if (world_.rank() == 0) {
          print("==================================================");
          print("CIS calculation FAILED with an exception:");
          print(e.what());
          print("==================================================");
        }
        results_["status"] = "failed";
        results_["error"] = e.what();
      }
    }
  }

  nlohmann::json results() const override { return results_; }

private:
  World &world_;
  std::shared_ptr<Nemo> reference_;
  std::filesystem::path ref_dir_;
};

class OEPApplication : public Application, public OEP {
public:
  explicit OEPApplication(World &w, const Params &p,
                          const std::shared_ptr<Nemo> &reference)
      : Application(p), OEP(w, p.get<OEP_Parameters>(), reference), world_(w),
        reference_(reference) {}

  // print_parameters
  void print_parameters(World &world) const override {
    if (world.rank() == 0)
      params_.get<OEP_Parameters>().print(OEP_Parameters::tag, "end");
  }

  void run(const std::filesystem::path &workdir) override {
    // 1) set up a namedspaced directory for this run
    PathManager pm(workdir, "oep");
    pm.create();
    world_.gop.fence();
    {
      ScopedCWD scwd(pm.dir());
      if (world_.rank() == 0) {
        std::cout << "Running OEP in " << pm.dir() << std::endl;
      }

      // 2) define the "checkpoint" file
      std::string label = "oep";
      auto ckpt = label + "_results.json";
      print("cc checkpoint file", ckpt);
      if (std::filesystem::exists(ckpt)) {
        if (world_.rank() == 0) {
          std::cout << "Found checkpoint file: " << ckpt << std::endl;
        }
        // read the checkpoint file
        std::ifstream ifs(ckpt);
        nlohmann::json j;
        ifs >> j;
        ifs.close();
      }

      // we could dump params_ to JSON and pass as argv if desired…
      try {
        const double time_scf_start = wall_time();
        this->value();
        const double time_scf_end = wall_time();
        if (world_.rank() == 0)
          printf(" at time %.1f\n", wall_time());

        if (world_.rank() == 0) {
          std::cout << std::setfill(' ');
          std::cout << "\n\n\n";
          std::cout << "--------------------------------------------------\n";
          std::cout << "MRA-OEP ended \n";
          std::cout << "--------------------------------------------------\n";
          std::cout << std::setw(25) << "time scf" << " = "
                    << time_scf_end - time_scf_start << "\n";
          std::cout << "--------------------------------------------------\n";
        }
        results_ = this->analyze();
      } catch (std::exception &e) {
        // Do not silently swallow: record the failure so the emitted
        // calc_info.json reflects it instead of looking like a clean run.
        if (world_.rank() == 0) {
          print("==================================================");
          print("OEP calculation FAILED with an exception:");
          print(e.what());
          print("==================================================");
        }
        results_["status"] = "failed";
        results_["error"] = e.what();
      }
    }
  }

  nlohmann::json results() const override { return results_; }

private:
  World &world_;
  std::shared_ptr<Nemo> reference_;

  // double energy_;
  // std::optional<Tensor<double>> dipole_;
  // std::optional<Tensor<double>> gradient_;
  // std::optional<real_function_3d> density_;
};

inline NextAction decide_next_action(bool at_protocol, bool archive_needed,
                                     bool archive_exists,
                                     bool all_properties_computed,
                                     bool restart_exists) {
  // We must recompute if any of these are true:
  const bool must_redo =
      !at_protocol // not at final protocol
      ||
      (archive_needed && !archive_exists) // user wants archive but it's missing
      || !all_properties_computed;        // a requested prop is missing

  if (!must_redo) {
    // We’re at final protocol, have all requested props, and either
    // the archive is present or not required.
    // If the archive isn't there but also not required, it's just a reload.
    if (!archive_exists && !archive_needed)
      return NextAction::ReloadOnly;
    return NextAction::Ok;
  }

  // We need work; decide between Restart vs Redo:
  return restart_exists ? NextAction::Restart : NextAction::Redo;
}

template <typename SCFParams>
NextAction valid(World &world, const SCFResultsTuple &results,
                 const SCFParams &params) {
  // Take a copy of the parameters
  auto [sr, pr, cr, optr] = results;

  // Required convergence for "final" protocol
  const auto vthresh = params.protocol().back(); // final protocol
  const auto vdconv = params.dconv();
  const bool archive_needed = params.save();

  // Requested outputs
  const bool need_energy = true;
  const bool need_dipole = params.dipole();
  const bool need_gradient = params.derivatives();

  if (world.rank() == 0) {
    print("Validating SCF results:");
    print(" Required protocol threshold: ", vthresh);
    print(" Required density convergence: ", vdconv);
    print(" Archive needed: ", archive_needed);
    print(" Need energy: ", need_energy);
    print(" Need dipole: ", need_dipole);
    print(" Need gradient: ", need_gradient);
  }

  // Files/paths
  const std::string archivename = params.prefix();
  const auto restart_path =
      std::filesystem::path(archivename + ".restartdata.00000");
  const bool archive_exists = std::filesystem::exists(restart_path);
  if (world.rank() == 0) {
    print("Restart file: ", restart_path.string());
  }

  // State in resultout the threshold refinement.
  //
  const bool at_protocol =
      (cr.converged_for_thresh == vthresh && cr.converged_for_dconv == vdconv);

  const auto pjson = sr.properties.to_json();
  const bool energy_ok = pjson.contains("energy");
  const bool dipole_ok = pjson.contains("dipole");
  const bool gradient_ok = pjson.contains("gradient");

  if (world.rank() == 0) {
    print("at_protocol: ", at_protocol);
    print("archive_needed: ", archive_needed);
    print("archive_exists: ", archive_exists);
    print("energy_ok: ", energy_ok);
    print("dipole_ok: ", dipole_ok);
    print("gradient_ok: ", gradient_ok);
  }

  // Only require props the user asked for
  const bool all_properties_computed = (need_energy ? energy_ok : true) &&
                                       (need_dipole ? dipole_ok : true) &&
                                       (need_gradient ? gradient_ok : true);

  bool is_gopt = false;
  bool gopt_ok = true;
  if (params.gopt()) {
    is_gopt = true;
  }
  if (is_gopt) {
    // check if the optimized geometry is available
    try {

      double gtol = params.gtol();
      // A geometry optimization is only valid if the max gradient has
      // dropped below gtol; otherwise it must be restarted/redone.
      // (Previously this set gopt_ok=true on non-convergence, which made an
      //  unconverged optimization look valid and suppressed the redo.)
      gopt_ok = (optr.max_gradient < gtol);
    } catch (...) {
      gopt_ok = false;
    }
  }

  // Decide action
  const bool must_redo = !at_protocol || (archive_needed && !archive_exists) ||
                         !all_properties_computed || !gopt_ok;

  // if we don't need to redo, we can either reload or return ok
  if (!must_redo)
    return (!archive_exists && !archive_needed) ? NextAction::ReloadOnly
                                                : NextAction::Ok;

  if (world.rank() == 0) {
    print("at_protocol: ", at_protocol);
    print("archive_needed: ", archive_needed);
    print("archive_exists: ", archive_exists);
    print("all_properties_computed: ", all_properties_computed);
    print("gopt_ok: ", gopt_ok);
  }
  // with we need to redo we can restart from the exisiting archive
  return archive_exists ? NextAction::Restart : NextAction::Redo;
}

struct moldft_lib {
  static constexpr const char *label() { return "moldft"; }

  vector<double> protocol;
  SCFResultsTuple last_results_;

  using Calc = SCF;

  NextAction valid(World &world, const SCFResultsTuple &results,
                   const Params &params) {
    last_results_ = results;
    return ::valid(world, results, params.get<CalculationParameters>());
  }

  // expose the live engine
  std::shared_ptr<Calc> calc(World &world, const Params &params) {
    if (!calc_)
      initialize_(world, params); // create once
    return calc_;
  }

  static void print_parameters() { Calc::print_parameters(); }
  // params get's changed by SCF constructor
  SCFResultsTuple run(World &world, const Params &params,
                      const NextAction next_action_) {
    const auto &molecule = params.get<Molecule>();
    const auto &params_copy = params;

    SCFResultsTuple results;
    auto &scf_res = std::get<0>(results);
    auto &opt_res = std::get<3>(results);
    auto &prop_res = std::get<1>(results);
    auto &conv_res = std::get<2>(results);

    if (next_action_ == NextAction::Ok ||
        next_action_ == NextAction::ReloadOnly) {
      // nothing to do
      return last_results_;
    }

    // NOTE: for NextAction::Restart, the SCF engine is already told to load
    // MOs from restartdata by the caller (SCFApplication::run sets the
    // 'restart' flag on params_ BEFORE the engine is constructed). Setting it
    // here on a local copy was dead code — the engine was already built.
    auto scf = calc(world, params_copy);
    // redirect any log files into outdir if needed…
    // Warm and fuzzy for the user
    if (world.rank() == 0) {
      print("\n\n");
      print(" MADNESS Hartree-Fock and Density Functional Theory Program");
      print(" ----------------------------------------------------------\n");
      print("\n");
      scf->molecule.print();
      print("\n");
      scf->param.print(CalculationParameters::tag);
    }
    // Come up with an initial OK data map
    if (world.size() > 1) {
      scf->set_protocol<3>(world, 1e-4);
      scf->make_nuclear_potential(world);
      scf->initial_load_bal(world);
    }
    // vama
    scf->set_protocol<3>(world, scf->param.protocol()[0]);
    double energy = 0.0;
    scf_res.is_opt = scf->param.gopt();

    if (scf_res.is_opt) {
      // Geometry convergence thresholds tied to final SCF protocol
      // following the Dalton-style rules:
      //
      //   wf_thresh = final SCF threshold (protocol().back())
      //
      //   |ΔE|   < max(1e-6, 2 * wf_thresh)
      //   ||g||  < max(1e-4, 2 * wf_thresh)
      //   ||dx|| < max(1e-4, 2 * wf_thresh)
      //
      // MolOpt uses:
      //   etol -> energy change tolerance
      //   gtol -> max gradient tolerance
      //   xtol -> max Cartesian step tolerance
      //
      const double wf_thresh = scf->param.protocol().back();
      const double etol = std::max(1.0e-7, 2.0 * wf_thresh);
      const double gxtol = std::max(1.0e-5, 2.0 * wf_thresh);
      const double gprec = std::max(1.0e-6, wf_thresh);

      MolOpt opt(scf->param.gmaxiter(), // maximum geometry iterations
                 0.1,   // maximum step in any Cartesian coordinate
                 etol,  // energy-change tolerance
                 gxtol, // gradient tolerance
                 gxtol, // step (Cartesian) tolerance
                 etol,  // assumed energy precision
                 gprec, // assumed gradient precision
                 (world.rank() == 0) ? 1 : 0, // print_level
                 scf->param.algopt());

      MolecularEnergy target(world, *scf);
      opt_res = opt.optimize_app(scf->molecule, target);
      auto new_mol = opt_res.final_geometry;

      Tensor<double> gradient;
      target.energy_and_gradient(new_mol, energy, gradient);
      scf_res.scf_molecule = new_mol;

      scf_res.properties.energy = energy;
      scf_res.properties.gradient = gradient;

      // write out the optimized geometry
      if (world.rank() == 0) {
        std::string geomfile = scf->param.prefix() + "_opt.xyz";
        std::ofstream ofs(geomfile);
        new_mol.print(ofs);
        ofs.close();
        print("optimized geometry written to ", geomfile);
        // write out mad.in with optimized geometry
      }

      // MolecularEnergy E(world, *scf);
      // energy = E.value(new_mol.get_all_coords().flat());
    } else {
      MolecularEnergy E(world, *scf);
      scf_res.scf_molecule = molecule;

      energy = E.value(scf->molecule.get_all_coords().flat());
      if (world.rank() == 0 && scf->param.print_level() > 0)
        E.output_calc_info_schema();
    }

    functionT rho = scf->make_density(world, scf->aocc, scf->amo);
    functionT brho = rho;
    if (scf->param.nbeta() != 0 && !scf->param.spin_restricted())
      brho = scf->make_density(world, scf->bocc, scf->bmo);
    rho.gaxpy(1.0, brho, 1.0);

    // optionally compute gradient, dipole, etc.
    Tensor<double> grad;
    if (scf->param.derivatives()) {
      grad = scf->derivatives(world, rho);
      scf->e_data.add_gradient(grad);
      scf_res.properties.gradient = grad;
    }

    tensorT dip;
    if (scf->param.dipole())
      dip = scf->dipole(world, scf->make_density(world, scf->aocc, scf->amo));

    scf->do_plots(world);

    conv_res.set_converged_thresh(FunctionDefaults<3>::get_thresh());
    conv_res.set_converged_dconv(scf->param.dconv());
    prop_res.energy = energy;
    prop_res.dipole = dip;
    prop_res.gradient = grad;

    scf_res.aeps = scf->aeps;
    scf_res.beps = scf->beps;
    scf_res.properties = prop_res;

    return results;
  }

private:
  void initialize_(World &world, const Params &params) {
    // write mad.in if missing
    const auto &cp = params.get<CalculationParameters>();
    const auto &mol = params.get<Molecule>();

    world.gop.fence();
    if (world.rank() == 0) {
      if (true) {
        // should always overwrite for now

        json in;
        in["dft"] = cp.to_json_if_precedence("defined");
        in["molecule"] = mol.to_json_if_precedence("defined");
        std::ofstream ofs("mad.in");
        write_json_to_input_file(in, {"dft"}, ofs);
        mol.print_defined_only(ofs);
      }
    }
    world.gop.fence();

    commandlineparser parser;
    parser.set_keyval("input", "mad.in");
    if (world.rank() == 0)
      ::print("input filename: ", parser.value("input"));

    FunctionDefaults<3>::set_pmap(pmapT(new LevelPmap<Key<3>>(world)));
    std::cout.precision(6);
    calc_ = std::make_shared<SCF>(world, parser);
  }

  std::shared_ptr<Calc> calc_;
}; // namespace moldft_lib

struct nemo_lib {
  using Calc = Nemo;
  static constexpr const char *label() { return "nemo"; }

  std::shared_ptr<Calc> calc(World &world, const Params &params) {
    if (!nemo_)
      initialize_(world, params);
    return nemo_;
  }

  static NextAction valid(World &world, const SCFResultsTuple &results,
                   const Params &params) {
    // Take a copy of the parameters
    return ::valid(world, results, params.get<CalculationParameters>());
  }

  static void print_parameters() { Calc::print_parameters(); }

  SCFResultsTuple run(World &world, const Params &params,
                      NextAction action = NextAction::Redo) {
    SCFResultsTuple results;
    auto nm = calc(world, params);
    nm->get_calc()->work_dir = std::filesystem::current_path();

    nm->value();
    PropertyResults pr = nm->analyze();
    // compute the hessian

    ConvergenceResults cr;
    cr.set_converged_thresh(nm->get_calc()->converged_for_thresh);
    cr.set_converged_dconv(nm->get_calc()->converged_for_dconv);

    SCFResults sr;
    sr.aeps = nm->get_calc()->aeps;
    sr.beps = nm->get_calc()->beps;
    sr.properties = pr;
    sr.scf_total_energy = nm->get_calc()->current_energy;

    if (nm->get_nemo_param().hessian())
      sr.properties.vibrations =
          nm->hessian(nm->get_calc()->molecule.get_all_coords());
    results = {sr, pr, cr, OptimizationResults()};

    return results;
  }

private:
  void initialize_(World &world, const Params &params) {
    nemo_ = std::make_shared<Nemo>(
        world, params.get<CalculationParameters>(),
        params.get<Nemo::NemoCalculationParameters>(), params.get<Molecule>());
  }

  std::shared_ptr<Calc> nemo_;
};
} // namespace madness
