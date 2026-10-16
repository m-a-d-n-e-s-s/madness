#pragma once

#include <madness/chem/InputWriter.hpp>
#include <madness/chem/ParameterManager.hpp>
#include <madness/chem/PathManager.hpp>
#include <madness/chem/RestartPlan.h>
#include <madness/chem/Results.h>
#include <madness/chem/molopt.h>
#include <filesystem>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <system_error>
#include <type_traits>

namespace madness {
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
  /// The id of each archive in `archives` when it was published (see ArchiveId);
  /// a consumer that records it can later tell whether the archive changed.
  std::map<std::string, ArchiveId> archive_ids;
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
  // engine was only constructed, never solved or given orbitals from an archive.
  if (reference.get_calc()->get_amo().empty()) {
    if (world.rank() == 0)
      print("empty reference orbitals in step", step_label);
    MADNESS_CHECK_THROW(
        false,
        "The upstream ground state has no orbitals: the reference engine was "
        "only constructed. cc2/cis/oep need a live converged reference");
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

  /// Optional pre-run hook, invoked collectively INSIDE the SCF work directory
  /// (cwd = the task dir, before the engine constructs its restart plan) and
  /// only when the engine is about to run, i.e. the stored results are not reused. The app layer uses it
  /// to lay down a ground-state seed archive (e.g. madqc: `dalton.dir` ->
  /// <prefix>.restartdata projected from the DALTON molden), which `restart
  /// auto` then picks up like any other archive. chem/ stays ignorant of where
  /// the seed comes from.
  using PreRunHook = std::function<void(World &, const Params &,
                                        const std::filesystem::path &)>;
  void set_pre_run_hook(PreRunHook h) { pre_run_hook_ = std::move(h); }

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

  /// bump whenever a stored result would change for the same orbitals and inputs
  /// -- a corrected property formula -- so results written before the fix are
  /// recomputed rather than reused (see post_key)
  static constexpr int results_schema = 1;

  // sets the calc working directory and runs the calculation
  //
  // Whether the stored results can stand is decided in two parts. The restart
  // planner -- the same one the engine consults -- says which orbitals are on
  // disk and whether they answer this request without iterating (geometry,
  // Hamiltonian, localization, orbital count, convergence and the `restart`
  // mode). The results file is reused only if it was computed from that very
  // archive (its archive_id) with the same post-SCF requests (post_key).
  // Anything else runs the engine, which plans again and, for converged orbitals,
  // only recomputes the properties.
  void run(const std::filesystem::path &workdir) override {
    std::string label = Library::label();
    PathManager pm(workdir, label);
    pm.create();
    {
      world_.gop.fence();
      ScopedCWD scwd(pm.dir());
      if (world_.rank() == 0) {
        std::cout << "Running SCF in " << pm.dir() << std::endl;
      }
      set_calc_workdir(pm.dir());

      const std::string results_file = label + ".results.json";
      const RestartPlan plan = plan_restart_from_disk();
      const nlohmann::json key = post_key();
      const nlohmann::json stored = read_results_collective(results_file);

      std::string why;
      if (plan.source != RestartSource::restartdata or plan.iterate)
        why = "the orbitals on disk do not answer this request (" + plan.why + ")";
      else if (plan.archive_id == 0)
        why = "the archive records no id, so no results can be tied to it";
      else if (not stored.is_object())
        why = "there is no readable " + results_file;
      else if (archive_id_from_string(stored.value("archive_id", std::string())) !=
               plan.archive_id)
        why = results_file + " was computed from a different archive";
      else if (stored.value("post_key", nlohmann::json()) != key)
        why = results_file + " answers different post-SCF requests";

      if (why.empty() and reuse_stored_results(stored)) {
        if (world_.rank() == 0)
          print("reusing", results_file, "computed from archive",
                archive_id_to_string(plan.archive_id));
        lib_.reload(world_, params_);
        return;
      }
      if (world_.rank() == 0)
        print("running the SCF engine:", why.empty() ? "stored results unreadable" : why);

      if (pre_run_hook_) {
        pre_run_hook_(world_, params_, pm.dir());
        world_.gop.fence();
      }
      scf_results = lib_.run(world_, params_);

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

      // Task-entry envelope: every task states its type; precision
      // is mirrored so consumers need not know the nested layout.
      results_["type"] = results_["scf"].value("model", std::string("scf"));
      if (results_["scf"].contains("precision"))
        results_["precision"] = results_["scf"]["precision"];

      // the archive these results were computed from: the one the engine just
      // wrote, or loaded without iterating. 0 (save false) never matches a plan.
      const ArchiveId id = engine_scf()->archive_id;
      write_results_file(results_file, {{"schema", 1},
                                        {"archive_id", archive_id_to_string(id)},
                                        {"post_key", key},
                                        {"results", results_}});
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
    // wipe its geometry. params_'s Molecule is the correct fallback.
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
      ctx.archive_ids["restartdata"] = engine_scf()->archive_id;
    } catch (...) {
      // best-effort; absence just means downstream restart discovery falls back
    }
  }

  nlohmann::json results() const override { return results_; }

private:
  /// the SCF that owns the orbitals: the engine itself for moldft, nemo's inner one
  std::shared_ptr<SCF> engine_scf() {
    if constexpr (std::is_same_v<Calc, SCF>)
      return calc();
    else
      return calc()->get_calc();
  }

  /// what the engine's own restart planner will decide, without reading orbitals
  ///
  /// Uses the constructed engine's parameters and molecule, so the comparison is
  /// against exactly what MolecularEnergy::value / Nemo::value will compare.
  RestartPlan plan_restart_from_disk() {
    const auto scf = engine_scf();
    const RestartCapabilities can = std::is_same_v<Calc, SCF>
                                        ? RestartCapabilities::all()
                                        : RestartCapabilities::restartdata_only();
    return make_restart_plan(world_, restart_mode_from_string(scf->param.restart()),
                             scf->param, scf->molecule, scf->restart_representation,
                             can, scf->hamiltonian_key());
  }

  /// the requests that change the results but not the orbitals
  ///
  /// Geometry, Hamiltonian, localization and convergence are deliberately absent:
  /// the archive_id stands for them, because the planner has checked that archive
  /// against the request.
  nlohmann::json post_key() const {
    const auto &cp = params_.get<CalculationParameters>();
    nlohmann::json k;
    k["results_schema"] = results_schema;
    k["dipole"] = cp.dipole();
    k["derivatives"] = cp.derivatives();
    // an added/removed/re-parameterized dispersion correction shifts the total
    // energy without touching the orbitals
    k["dispersion"] = cp.dispersion();
    k["dispersion_functional"] = cp.dispersion_functional();
    k["dispersion_atm"] = cp.dispersion_atm();
    if constexpr (!std::is_same_v<Calc, SCF>)
      k["hessian"] = params_.get<Nemo::NemoCalculationParameters>().hessian();
    return k;
  }

  /// the results file, read on rank 0 and broadcast; null if absent or unparsable
  nlohmann::json read_results_collective(const std::string &filename) const {
    std::string text;
    if (world_.rank() == 0 and std::filesystem::exists(filename)) {
      try {
        // a truncated or corrupt file must degrade to "not there", not throw
        text = nlohmann::json::parse(std::ifstream(filename)).dump();
      } catch (...) {
        print("WARNING: could not parse", filename, "-- ignoring it");
      }
    }
    world_.gop.broadcast_serializable(text, 0);
    return text.empty() ? nlohmann::json() : nlohmann::json::parse(text);
  }

  /// adopt the stored results; false if they do not parse
  bool reuse_stored_results(const nlohmann::json &stored) {
    try {
      const nlohmann::json &r = stored.at("results");
      auto &[scf_r, properties, convergence, optr] = scf_results;
      scf_r.from_json(r.at("scf"));
      properties.from_json(r.at("properties"));
      convergence.from_json(r.at("convergence"));
      results_ = r;
      return true;
    } catch (...) {
      scf_results = SCFResultsTuple();
      return false;
    }
  }

  /// write the results file atomically (tmp + rename), rank 0
  void write_results_file(const std::string &filename,
                          const nlohmann::json &j) const {
    if (world_.rank() != 0) return;
    const std::string tmp = filename + ".tmp";
    bool ok = true;
    {
      std::ofstream ofs(tmp);
      ofs << j.dump(4);
      ofs.flush();
      ok = static_cast<bool>(ofs);
    }
    std::error_code ec;
    if (ok) {
      std::filesystem::rename(tmp, filename, ec);
      ok = not ec;
    }
    if (ok) {
      print("Written results file: ", filename);
    } else {
      std::filesystem::remove(tmp, ec);
      print("ERROR: failed to write results file (disk full?): ", filename);
    }
  }

  World &world_;

  PreRunHook pre_run_hook_;
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

struct moldft_lib {
  static constexpr const char *label() { return "moldft"; }

  vector<double> protocol;

  using Calc = SCF;

  // expose the live engine
  std::shared_ptr<Calc> calc(World &world, const Params &params) {
    if (!calc_)
      initialize_(world, params); // create once
    return calc_;
  }

  static void print_parameters() { Calc::print_parameters(); }

  /// Rehydrate the engine when the stored results are reused, so a downstream
  /// step gets a reference with orbitals rather than a bare freshly-constructed
  /// SCF. The restart plan has already checked the archive against the request;
  /// load_mos handles k-projection and the threshold itself, and the protocol
  /// is set first to match the order Nemo::value uses.
  void reload(World &world, const Params &params) {
    auto scf = calc(world, params);
    scf->set_protocol<3>(world,
                         params.get<CalculationParameters>().protocol().back());
    scf->load_mos(world);
  }

  // params get's changed by SCF constructor
  SCFResultsTuple run(World &world, const Params &params) {
    const auto &molecule = params.get<Molecule>();
    const auto &params_copy = params;

    SCFResultsTuple results;
    auto &scf_res = std::get<0>(results);
    auto &prop_res = std::get<1>(results);
    auto &conv_res = std::get<2>(results);

    // The engine reads the restartdata header itself and decides whether to
    // load, iterate or skip (plan_restart).
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

    // total density: 2*rho_alpha when restricted, rho_alpha + rho_beta otherwise
    // (no beta term when nbeta == 0, e.g. a fully spin-polarized reference)
    functionT rho = scf->make_density(world, scf->aocc, scf->amo);
    if (scf->param.spin_restricted())
      rho.scale(2.0);
    else if (scf->param.have_beta())
      rho.gaxpy(1.0, scf->make_density(world, scf->bocc, scf->bmo), 1.0);

    // optionally compute gradient, dipole, etc.
    Tensor<double> grad;
    if (scf->param.derivatives()) {
      grad = scf->derivatives(world, rho);
      scf->e_data.add_gradient(grad);
      scf_res.properties.gradient = grad;
    }

    tensorT dip;
    if (scf->param.dipole())
      dip = scf->dipole(world, rho);

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
    scf_res.uses_dftd3 = scf->dispersion.active();
#ifdef MADNESS_HAS_PCM
    scf_res.uses_pcm = (scf->pcm_param.solvent() != "none");
#else
    scf_res.uses_pcm = false;
#endif
#ifdef MADNESS_HAS_LIBXC
    scf_res.uses_libxc = scf->xc.uses_libxc_backend();
#else
    scf_res.uses_libxc = false;
#endif

    // SCF task record: what this SCF actually ran at and
    // what it produced, under QCSchema names where they exist. Pure bookkeeping
    // from quantities the solve already holds — nothing here changes a number.
    {
      const auto last = scf->e_data.last();
      auto get = [&](const char *k) {
        auto it = last.find(k);
        return it == last.end() ? 0.0 : it->second;
      };
      scf_res.scf_total_energy = energy;   // moldft never set this (0.0); only the nemo path did (Applications.hpp:1313)

      // Plan without iterations (converged archive, restart read_only): value()
      // returns without ever calling scf->e_data.add_data(), so last() is
      // empty here. Recording an all-zero energy decomposition and
      // scf_iterations = 0 next to the real (archive-derived) total energy
      // would misreport a converged reload as an unconverged, zero-energy
      // solve, so leave `energies` unset and iterations at their -1 defaults.
      if (!last.empty()) {
        nlohmann::json e;
        e["nuclear_repulsion_energy"]      = get("e_nrep");
        e["scf_kinetic_energy"]            = get("e_kinetic");
        e["scf_nuclear_attraction_energy"] = get("e_nuclear");
        e["scf_coulomb_energy"]            = get("e_coulomb");
        e["scf_pcm_energy"]                = get("e_pcm");
        e["scf_one_electron_energy"]       = get("e_kinetic") + get("e_nuclear") + get("e_local");
        e["scf_two_electron_energy"]       = get("e_coulomb") + get("e_xc");   // HF: e_xc is exact exchange
        if (scf->xc.is_dft()) e["scf_xc_energy"] = get("e_xc");
        scf_res.energies = e;
        scf_res.scf_iterations = scf->e_data.iterations();
        conv_res.iterations = scf->e_data.iterations();
      }
      scf_res.xc = scf->param.xc();

      // Collective: one reduction total instead of one per orbital
      // (Function::size() is itself a collective global sum).
      world.gop.fence();
      std::size_t ncoeff = 0;
      for (const auto &f : scf->amo) ncoeff += f.size_local();
      for (const auto &f : scf->bmo) ncoeff += f.size_local();
      world.gop.sum(ncoeff);
      scf_res.precision = {{"k", FunctionDefaults<3>::get_k()},
                           {"thresh", FunctionDefaults<3>::get_thresh()},
                           {"protocol", scf->param.protocol()},
                           {"econv", scf->param.econv()},
                           {"dconv", scf->param.dconv()},
                           {"L", scf->param.L()},
                           {"ncoeff", ncoeff}};

      // converged_for_thresh/dconv always come from the engine (on the reload
      // path, from the archive header -- SCF.cc:456) and describe the stored
      // wavefunction truthfully either way, so this status check runs
      // unconditionally. It deliberately checks against param.dconv() (what
      // the deck asked for), which is stricter than the engine's own reload
      // test max(protocol.back(), dconv) (SCF.h:625-629).
      const double finest = scf->param.protocol().empty()
                                ? FunctionDefaults<3>::get_thresh()
                                : scf->param.protocol().back();
      conv_res.status = (scf->converged_for_thresh <= finest &&
                         scf->converged_for_dconv <= scf->param.dconv())
                            ? "converged" : "unconverged";
    }

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
        // The `pcm` group has to make the same round trip as `dft`: the SCF is
        // rebuilt from this regenerated mad.in, so anything omitted here is
        // silently lost. Written unconditionally -- an empty `pcm/end` block is
        // harmless, and PCMParameters is inert unless dft's pcm_data is set.
        in["pcm"] = params.get<PCMParameters>().to_json_if_precedence("defined");
        // `prefix` must be carried explicitly. It is the one parameter that is
        // DERIVED from information the engine cannot recompute -- the name of
        // the original input file (ParameterManager.hpp) -- and this round trip
        // keeps only user-defined values. Without it the engine falls back to
        // the "mad" default and writes mad.restartdata, where the restart
        // planner (which reads <prefix>.restartdata) never finds it. Everything else that set_derived_values()
        // computes is re-derived identically by the SCF ctor.
        in["dft"]["prefix"] = cp.prefix();
        std::ofstream ofs("mad.in");
        write_json_to_input_file(in, {"dft", "pcm"}, ofs);
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

  SCFResultsTuple run(World &world, const Params &params) {
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
    sr.xc = nm->get_calc()->param.xc();
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
        params.get<Nemo::NemoCalculationParameters>(), params.get<Molecule>(),
        params.get<PCMParameters>());
  }

  std::shared_ptr<Calc> nemo_;
};
} // namespace madness
