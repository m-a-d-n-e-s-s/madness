#pragma once

#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/PathManager.hpp>
#include <madness/chem/Results.h>
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
  /// Live Nemo-based ground-state reference, if an upstream nemo step exposed
  /// one. Deliberately a SECOND field rather than a widening of `reference`:
  /// SCF (SCF.h) and Nemo (via NemoBase) share no base class, and aliasing a
  /// Nemo's inner SCF into `reference` would silently re-point consumers that
  /// expect a standalone SCF engine (ResponseApplication).
  std::shared_ptr<Nemo> nemo_reference;
  /// Named archive/output paths (absolute), e.g. "restartdata" -> path.
  std::map<std::string, std::filesystem::path> archives;
  /// Free-form JSON for artifacts not yet first-class.
  nlohmann::json blob = nlohmann::json::object();
};

/// True iff `a` and `b` describe the same nuclear framework: equal atom count
/// and, per atom, equal element/charge/mass and a displacement below the 1e-10
/// threshold Atom::operator== already applies.
///
/// Deliberately ignores derived state (rcut, pointgroup, field, molecular
/// parameters). Those differ across a JSON round trip -- from_json re-orients
/// the molecule -- which is why Molecule::operator== is not usable as a geometry
/// guard. Lives here rather than in molecule.h only to avoid making every
/// translation unit in the tree rebuild for it; molecule.h is its natural home
/// once something outside this header needs it.
inline bool same_nuclear_framework(const Molecule &a, const Molecule &b) {
  if (a.natom() != b.natom())
    return false;
  for (unsigned int i = 0; i < a.natom(); ++i)
    if (!(a.get_atom(i) == b.get_atom(i)))
      return false;
  return true;
}

/// Consumer-side precondition for the steps that cannot adopt an upstream
/// reference (madqc ARCHITECTURE_ROADMAP change 1).
///
/// CC2, TDHF and OEP build their engine from a reference captured when the
/// workflow is ASSEMBLED, not when it runs: CC2's ctor calls set_protocol on it
/// and builds its own TDHF + CCPotentials from it, TDHF freezes its derived
/// parameters off it, and OEP's own Nemo base is constructed FROM it. So a
/// reference or geometry published later cannot be adopted, and the partial
/// hooks that look like they would do it (TDHF::set_reference,
/// OEP::set_reference, CCPotentials::reset_nemo) each swap only part of the
/// state and would leave a silently inconsistent object. Until engine
/// construction moves into run(), the honest thing a consumer can do is verify
/// that the threaded context agrees with what it was built on -- so a
/// mis-chained workflow fails loudly instead of writing out a correlation
/// energy computed at the wrong geometry.
///
/// @param[in] ctx         the context Workflow::run threads task-to-task
/// @param[in] reference    the build-time reference this step's engine uses
/// @param[in] step_label  "cc2" / "cis" / "oep", for diagnostics only
/// @param[in] world       for rank-0 printing
inline void check_context_matches_reference(const StepContext &ctx,
                                            const Nemo &reference,
                                            const char *step_label,
                                            World &world) {
  // Messages passed to MADNESS_CHECK_THROW must be string literals:
  // MadnessException stores a bare const char* and does not copy it, so a
  // temporary std::string's c_str() would dangle by the time what() is called.
  // Runtime detail goes into the rank-0 print next to each check instead.
  if (ctx.nemo_reference && ctx.nemo_reference.get() != &reference) {
    if (world.rank() == 0)
      print("StepContext reference mismatch in step", step_label);
    MADNESS_CHECK_THROW(
        false,
        "StepContext carries a different ground-state reference than this step "
        "was constructed with; cc2/cis/oep build their engine when the workflow "
        "is assembled and cannot adopt a reference published later. See "
        "src/apps/madqc/STEP_CONTEXT.md");
  }

  // Only meaningful when the published molecule came from some producer other
  // than our own reference. When the upstream reference IS our engine,
  // ctx.molecule is by construction the geometry it ran at, and any difference
  // is JSON-round-trip / re-orientation noise against a 1e-10 threshold.
  if (!ctx.nemo_reference && ctx.molecule && ctx.molecule->natom() > 0) {
    if (!same_nuclear_framework(*ctx.molecule, reference.molecule())) {
      if (world.rank() == 0) {
        print("geometry published by an upstream step, in step", step_label);
        ctx.molecule->print();
        print("geometry this step's engine was constructed with:");
        reference.molecule().print();
      }
      MADNESS_CHECK_THROW(
          false,
          "An upstream step published a geometry that differs from the one this "
          "step's engine was constructed with; cc2/cis/oep freeze the nuclear "
          "framework at workflow-assembly time and cannot run at an upstream "
          "optimized or displaced geometry. Run the method as its own workflow "
          "at that geometry; see src/apps/madqc/STEP_CONTEXT.md");
    }
  }

  // A reference engine that was constructed but never solved or reloaded has no
  // orbitals: CC2 would call compute_fock_matrix on an empty amo, TDHF and OEP
  // would fail their own check_converged gate much deeper in. Reachable when the
  // upstream SCF results were reused (NextAction::ReloadOnly) and no restartdata
  // archive existed to reload the MOs from -- i.e. the deck set `save false`.
  if (reference.get_calc()->get_amo().empty()) {
    if (world.rank() == 0)
      print("empty reference orbitals in step", step_label);
    MADNESS_CHECK_THROW(
        false,
        "The upstream ground state has no orbitals: its results were reused from "
        "a checkpoint, but no restartdata archive was available to reload the "
        "MOs from, so the reference engine was only constructed. cc2/cis/oep "
        "need a live converged reference. Set `save 1` in the dft group, or "
        "delete the reference step's calc_info.json to recompute");
  }

  if (world.rank() == 0) {
    print("step", step_label,
          ": ground-state reference confirmed via StepContext");
    const auto it = ctx.archives.find("restartdata");
    if (it != ctx.archives.end())
      print("  upstream restartdata archive:", it->second.string());
  }
}

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
      // ... and against reusing one computed for a DIFFERENT Hamiltonian. The
      // engine makes this call for itself once it is built (plan_restart), but
      // nothing here builds it: valid() looks only at thresholds, properties and
      // the archive, all of which a changed `xc` leaves intact. So `madqc
      // --dft="xc=lda"` in a directory holding an `xc=hf` checkpoint returned the
      // HF energy as the LDA answer without ever constructing an SCF.
      if (action != madness::NextAction::Redo &&
          !checkpoint_hamiltonian_matches(j)) {
        if (world_.rank() == 0)
          print("WARNING: checkpoint was computed for a different Hamiltonian "
                "(or does not record which); ignoring checkpoint and recomputing.");
        scf_results = empty_results;
        action = madness::NextAction::Redo;
      }
      // NB: nothing needs to be pushed into params_ for a restart. The engine
      // decides for itself where its orbitals come from (plan_restart, called
      // from MolecularEnergy::value and Nemo::value), reading the archive's
      // header rather than a flag set from out here. What used to live at this
      // point -- set_user_defined_value("restart", true) before the engine was
      // constructed, because a flag set afterwards never reached the ctor -- is
      // therefore gone along with the boolean it set.
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
        // Both actions mean the same thing here -- run the engine. Restart vs
        // Redo used to select where the orbitals came from; that is now the
        // engine's own decision (plan_restart, from the archive header), so the
        // distinction survives only as a diagnostic. The action is still passed
        // through for the log rather than hardcoding one.
        scf_results = lib_.run(world_, params_, action);
      } else if (action == madness::NextAction::Ok) {
        // Results stand as checkpointed -- but the engine must not be left as a
        // bare construction: a downstream cc2/cis/oep step needs a reference with
        // orbitals. Ok implies the restartdata archive exists (see valid()), so
        // the MOs can be reloaded from it.
        //
        // NOTE: Ok is currently UNREACHABLE, for a different reason on each path,
        // so this branch is correct-but-dormant rather than exercised:
        //  - nemo: valid() tests at_protocol as
        //    `converged_for_thresh == protocol().back()`, but Nemo::value drives
        //    its SCFProtocol from econv, so converged_for_thresh never equals the
        //    final protocol rung -> must_redo stays true -> Restart.
        //  - moldft: valid() looks for `params.prefix() + ".restartdata.00000"`,
        //    while the engine is constructed from the mad.in written by
        //    moldft_lib::initialize_ and so writes `mad.restartdata.00000`
        //    -> archive_exists is always false -> Redo.
        // Both are recorded in ARCHITECTURE_ROADMAP.md. Reuse still works today
        // via Restart (which reloads the MOs), it just costs the archive read.
        lib_.reload(world_, params_);
      } else {
        // ReloadOnly: results stand, and there is no archive to reload from, so
        // the engine has no orbitals. Harmless for a standalone scf/nemo run;
        // check_context_matches_reference reports it if a chained step then
        // needs a live reference.
        lib_.calc(world_, params_); // just set up the calc without running
      }

      // // Need work (Restart or Redo) — both call run()
      // scf_results = lib_.run(world_, params_);

      results_["scf"] = std::get<0>(scf_results).to_json();
      results_["properties"] = std::get<1>(scf_results).to_json();
      results_["convergence"] = std::get<2>(scf_results).to_json();
      results_["molecule"] = std::get<0>(scf_results).scf_molecule.to_json();
      results_["optimization_results"] = std::get<3>(scf_results).to_json();
      // what these numbers are a solution of, so the next invocation can tell
      // whether they answer its question -- see checkpoint_hamiltonian_matches
      results_["hamiltonian"] = checkpoint_hamiltonian();
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
    } else if constexpr (std::is_base_of_v<Nemo, Calc>) {
      // The nemo-based ground state used by the cc2/cis/oep chains. NOT also
      // published as ctx.reference: that field means "a standalone SCF engine",
      // and re-pointing it at a Nemo's inner SCF would change what
      // ResponseApplication consumes.
      ctx.nemo_reference = calc(); // shared_ptr<Nemo>
    }
    // Publish the geometry this step actually ran at. The natom() guard matters:
    // a producer that leaves SCFResults::scf_molecule unset would otherwise
    // publish a DEFAULT-CONSTRUCTED empty Molecule, and a consumer that assigns
    // ctx.molecule straight into its params (as ResponseApplication does) would
    // wipe its geometry. params_'s Molecule is the correct fallback -- on the
    // Ok/ReloadOnly path run() has already refreshed it from the checkpoint.
    if (results_.contains("molecule")) {
      try {
        Molecule m;
        m.from_json(results_["molecule"]);
        if (m.natom() > 0)
          ctx.molecule = std::move(m);
        else
          ctx.molecule = params_.get<Molecule>();
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
    return same_nuclear_framework(params_.get<Molecule>(), ckpt_mol);
  }

  /// what operator this checkpoint's numbers are a solution OF
  ///
  /// The molecule alone does not identify a calculation. Changing `xc` in place
  /// and rerunning is an everyday thing to do, and the geometry, the thresholds
  /// and the archive are all still valid for it -- so without this the cached
  /// results passed every test valid() applies and the engine was never built.
  /// `localize` is included because it selects the orbitals the eigenvalues and
  /// the Fock matrix in this file describe, even though it leaves the energy
  /// alone.
  nlohmann::json checkpoint_hamiltonian() const {
    const auto &cp = params_.get<CalculationParameters>();
    nlohmann::json h;
    h["xc"] = cp.xc();
    h["localize"] = cp.localize_method();
    // an added/removed/re-parameterized dispersion correction shifts the total
    // energy without touching the orbitals, so a checkpoint written without it
    // would otherwise be reused and its energy reported as this run's answer
    h["dispersion"] = cp.dispersion();
    h["dispersion_functional"] = cp.dispersion_functional();
    h["dispersion_atm"] = cp.dispersion_atm();
    if constexpr (!std::is_same_v<Calc, SCF>) {
      // Same spelling SCF::restart_ncf uses, so the checkpoint and the
      // restartdata header agree on what "the same ncf" means.
      const auto ncf = params_.get<Nemo::NemoCalculationParameters>().ncf();
      h["ncf"] = ncf.first + ":" + std::to_string(ncf.second);
    }
    return h;
  }

  /// true if the checkpoint solves the operator this run is asking about
  ///
  /// A checkpoint with no `hamiltonian` block was written by a build that did
  /// not record one, so what it solved cannot be established. That is treated as
  /// a mismatch rather than waved through: the failure being guarded against is
  /// a wrong energy reported as this run's answer, and the cost of being wrong
  /// here is one recompute that then writes the block. This is deliberately
  /// stricter than checkpoint_geometry_matches(), which waves through what it
  /// cannot parse -- there a mismatch merely wastes work, here it is silent.
  bool checkpoint_hamiltonian_matches(const nlohmann::json &j) const {
    if (!j.contains("hamiltonian"))
      return false;
    return j.at("hamiltonian") == checkpoint_hamiltonian();
  }

  World &world_;
  Library lib_; // owns shared_ptr<Engine>
  SCFResultsTuple scf_results;
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

  /// Verify-only participation in the StepContext dataflow: this step's engine
  /// was built from a reference captured at workflow-assembly time and cannot be
  /// re-pointed. See check_context_matches_reference.
  void consume_context(const StepContext &ctx) override {
    check_context_matches_reference(ctx, *reference_, "cc2", world_);
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

  /// Verify-only participation in the StepContext dataflow: this step's engine
  /// was built from a reference captured at workflow-assembly time and cannot be
  /// re-pointed. See check_context_matches_reference.
  void consume_context(const StepContext &ctx) override {
    check_context_matches_reference(ctx, *reference_, "cis", world_);
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

  /// Verify-only participation in the StepContext dataflow: this step's engine
  /// was built from a reference captured at workflow-assembly time and cannot be
  /// re-pointed. See check_context_matches_reference.
  void consume_context(const StepContext &ctx) override {
    check_context_matches_reference(ctx, *reference_, "oep", world_);
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
  // Take a copy of the parameters. The 4th element (OptimizationResults) is
  // deliberately not bound: an SCF task no longer optimizes, so there is
  // nothing of it to validate here.
  auto [sr, pr, cr, optr_unused] = results;
  (void)optr_unused;

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
  // "at least as good as requested", not "exactly equal": these are thresholds,
  // and exact float equality on them made a run converged to 1e-6 look invalid
  // against a request for 1e-6 whenever the two were computed differently.
  const bool at_protocol =
      (cr.converged_for_thresh <= vthresh && cr.converged_for_dconv <= vdconv);

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

  // NB: no geometry-optimization validity test here any more. An SCF task no
  // longer optimizes anything -- that is qcapp::OptimizeDriver's job, and it
  // keeps no checkpoint of its own to validate (the optimizer's hessian and
  // step history are not persisted, so an interrupted optimization restarts
  // from the input geometry while the SCF underneath still restarts per
  // geometry as usual).

  // Decide action
  const bool must_redo = !at_protocol || (archive_needed && !archive_exists) ||
                         !all_properties_computed;

  // if we don't need to redo, we can either reload or return ok
  if (!must_redo)
    return (!archive_exists && !archive_needed) ? NextAction::ReloadOnly
                                                : NextAction::Ok;

  if (world.rank() == 0) {
    print("at_protocol: ", at_protocol);
    print("archive_needed: ", archive_needed);
    print("archive_exists: ", archive_exists);
    print("all_properties_computed: ", all_properties_computed);
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

  /// Rehydrate the engine when the checkpoint results are reused, so a
  /// downstream step gets a reference with orbitals rather than a bare
  /// freshly-constructed SCF. load_mos handles k-projection and the threshold
  /// itself; the protocol is set first to match the order Nemo::value uses.
  void reload(World &world, const Params &params) {
    auto scf = calc(world, params);
    scf->set_protocol<3>(world,
                         params.get<CalculationParameters>().protocol().back());
    scf->load_mos(world);
  }

  // params get's changed by SCF constructor
  SCFResultsTuple run(World &world, const Params &params,
                      const NextAction next_action_) {
    const auto &molecule = params.get<Molecule>();
    const auto &params_copy = params;

    SCFResultsTuple results;
    auto &scf_res = std::get<0>(results);
    auto &prop_res = std::get<1>(results);
    auto &conv_res = std::get<2>(results);

    if (next_action_ == NextAction::Ok ||
        next_action_ == NextAction::ReloadOnly) {
      // nothing to do
      return last_results_;
    }

    // NOTE: NextAction::Restart needs nothing done here. The engine reads the
    // restartdata header itself and decides whether to load, iterate or skip
    // (plan_restart); NextAction only says whether the engine has to run at all.
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
    scf->dispersion.print_citation(world);
    double energy = 0.0;
    // An SCF task computes an energy at one geometry. Geometry optimization is
    // its own workflow task now -- `madqc --optimize --wf=scf`,
    // qcapp::OptimizeDriver in chem/Drivers.hpp -- which drives the same MolOpt
    // over the same MolecularEnergy target, derives its thresholds from the
    // `optimization` group, and publishes the optimized geometry downstream.
    // scf_res.is_opt stays false; the field and the OptimizationResults slot in
    // SCFResultsTuple remain because the driver fills them.
    MolecularEnergy E(world, *scf);
    scf_res.scf_molecule = molecule;

    energy = E.value(scf->molecule.get_all_coords().flat());
    if (world.rank() == 0 && scf->param.print_level() > 0)
      E.output_calc_info_schema();

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

    // report what the SCF actually achieved, not what was requested -- taking
    // these from FunctionDefaults/param made an unconverged run checkpoint as
    // converged, so the next invocation skipped it.
    conv_res.set_converged_thresh(scf->converged_for_thresh);
    conv_res.set_converged_dconv(scf->converged_for_dconv);
    prop_res.energy = energy;
    prop_res.dipole = dip;
    prop_res.gradient = grad;

    scf_res.aeps = scf->aeps;
    scf_res.beps = scf->beps;
    scf_res.scf_dispersion_correction_energy =
        scf->dispersion.energy(world, scf->molecule);
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
        // `prefix` must be carried explicitly. It is the one parameter that is
        // DERIVED from information the engine cannot recompute -- the name of
        // the original input file (ParameterManager.hpp) -- and this round trip
        // keeps only user-defined values. Without it the engine falls back to
        // the "mad" default and writes mad.restartdata, while valid() looks for
        // <prefix>.restartdata.00000: archive_exists is then always false and
        // the restart can never fire. Everything else that set_derived_values()
        // computes is re-derived identically by the SCF ctor.
        in["dft"]["prefix"] = cp.prefix();
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

  /// Rehydrate the engine when the checkpoint results are reused, so downstream
  /// steps (cc2/cis/oep) get a reference with orbitals rather than a bare
  /// freshly-constructed Nemo. Goes through the engine's own no_compute path:
  /// Nemo::value then loads the MOs from the archive, builds the nuclear
  /// correlation factor via set_protocol, records coords_sum (so check_converged
  /// passes downstream) and skips the SCF iterations. Setting no_compute after
  /// construction is safe because it is read in value(), unlike `restart` which
  /// SCFApplication::run must set before the engine is built.
  void reload(World &world, const Params &params) {
    auto nm = calc(world, params);
    nm->get_calc()->param.set_user_defined_value("no_compute", true);
    nm->value(nm->molecule().get_all_coords());
  }

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
    sr.scf_dispersion_correction_energy = nm->get_calc()->dispersion.energy(
        world, nm->get_calc()->molecule);
    // The geometry this reference was solved at. Without it, results_["molecule"]
    // (and therefore the checkpoint and ctx.molecule) reports an empty molecule,
    // and checkpoint_geometry_matches compares 0 atoms against N and rejects
    // every nemo-path checkpoint -- so restarts always recomputed.
    sr.scf_molecule = nm->get_calc()->molecule;

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
