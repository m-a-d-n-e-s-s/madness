#ifndef MOLRESPONSE_V3_MADQC_ADAPTER_HPP
#define MOLRESPONSE_V3_MADQC_ADAPTER_HPP

// -----------------------------------------------------------------------------
// molresponse_v3_lib — madqc adapter (doc 16 R3 / doc 05 seam).
//
// Satisfies the duck-typed interface ResponseApplication<Library> expects
// (Applications.hpp): `Library::label()` + `Library::run_response(world, params,
// scf, outdir) -> Results{metadata, properties, vibrational_analysis,
// raman_spectra}`. So `ResponseApplication<molresponse_v3_lib>` runs the v3
// pipeline through the same madqc workflow path as v2's molresponse_lib —
// enabling a SAME-INPUT calc_info.json parity check (engine = v2 vs v3).
//
// It builds a v3 GroundState from the moldft restart ARCHIVE (resolved from
// scf_calc->work_dir, like v2's make_ground_context — NOT from the in-memory
// SCF, whose MOs may be unloaded on a restart-in-place run), maps the
// ResponseParameters input deck to a v3 ResponsePlan, and calls
// run_response_with_ground.
//
// SCOPE: multi-property mapping (R3b) — polarizability + hyperpolarizability +
// single-component raman + resonant/excited, deduped via merge_plans.
// This header lives in the v3 app and is included by madqc.cpp (the app links
// both MADchem and v3's GroundState.cpp); MADchem/WorkflowBuilders must NOT
// include it (that would be a circular library dependency — the engine is
// selected in madqc.cpp instead).
// -----------------------------------------------------------------------------

#include <apps/molresponse/orchestrator/response_workflow.hpp>
#include <apps/molresponse/solvers/dalton_import.hpp>

#include <madness/chem/CalculationParameters.h>
#include <madness/chem/ParameterManager.hpp>   // Params
#include <madness/chem/ResponseParameters.hpp>
#include <madness/chem/SCF.h>
#include <nlohmann/json.hpp>
#include <madness/world/worldprofile.h>

#include <algorithm>
#include <filesystem>
#include <memory>
#include <vector>

/// Global namespace (mirrors `molresponse_lib`) so madqc can write
/// `ResponseApplication<molresponse_v3_lib>`.
struct molresponse_v3_lib {
  /// Output subdir name + the interface ResponseApplication reads.
  static const char *label() { return "molresponse"; }

  /// Structured result returned to ResponseApplication (→ calc_info.json).
  struct Results {
    nlohmann::json metadata;
    nlohmann::json properties;
    nlohmann::json vibrational_analysis;  // empty for alpha (R3a)
    nlohmann::json raman_spectra;         // empty for alpha (R3a)
  };

  inline static Results
  run_response(madness::World &world, const Params &params,
               const std::shared_ptr<madness::SCF> &scf_calc,
               const std::filesystem::path &outdir) {
    using namespace madness;
    using namespace molresponse_v3;

    const auto &cp = params.get<CalculationParameters>();
    const auto &rp = params.get<ResponseParameters>();
    std::vector<double> protocol = cp.protocol();
    MADNESS_CHECK(!protocol.empty());

    // Deck `seed.start_rung fine` — a dalton.dir-seeded run skips straight to
    // the FINEST rung (W6 between-pole finding: the coarse rung hits maxiter
    // unconverged and launders away the seed's head start). Must happen
    // BEFORE set_response_protocol / gs.prepare / run_dalton_import below so
    // the seed projection lands at the fine rung's (k, thresh). Default
    // 'coarse' = full ladder, unchanged. Shared helper with the standalone
    // driver (apply_seed_start_rung), so the two surfaces agree.
    if (apply_seed_start_rung(protocol, rp.seed_start_rung(),
                              !rp.dalton_dir().empty())) {
      if (world.rank() == 0)
        print("response: seed.start_rung=fine — dalton.dir seed starts the "
              "ladder at thresh", protocol.front());
    } else if (rp.seed_start_rung() == "fine" && rp.dalton_dir().empty() &&
               world.rank() == 0) {
      print("response: seed.start_rung=fine ignored — no dalton.dir seed "
            "configured (full ladder runs)");
    }

    // 1. Ground state from the moldft ARCHIVE (not the in-memory SCF). On a
    //    restart-in-place run madqc validates the SCF as "Ok" and never loads the
    //    MOs into memory (lib_.calc() just constructs the SCF), so scf_calc->amo is
    //    empty — building from the live SCF then segfaults in build_fock_matrices.
    //    Loading from the archive (exactly like v2's make_ground_context) avoids
    //    that AND lets prepare() reproject pristine MOs on each protocol climb.
    //
    //    Resolve the moldft work dir RELATIVE to the response `outdir`, exactly
    //    like v2's make_ground_context (MolresponseLib.hpp ~1154) and the
    //    CC2/TDHF/OEP applications. scf_calc->work_dir is stored relative to the
    //    top calc dir, but ResponseApplication::run has already chdir'd (ScopedCWD)
    //    into `outdir` (the response task dir), so using work_dir raw makes the
    //    archive lookup resolve against the wrong cwd → "could not find file:
    //    <work_dir>/<prefix>.restartdata" on multi-node madqc runs.
    namespace fs = std::filesystem;
    const double L = cp.L();
    set_response_protocol(world, L, protocol.front());
    const std::string prefix     = cp.prefix();
    const fs::path    moldft_dir = fs::proximate(scf_calc->work_dir, outdir);
    const std::string archive    = (moldft_dir / (prefix + ".restartdata")).string();
    const std::string fock_json  = (moldft_dir / (prefix + ".fock.json")).string();
    GroundState gs = GroundState::from_archive(world, archive, scf_calc->molecule);
    const double thresh = FunctionDefaults<3>::get_thresh();
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * thresh));
    gs.prepare(world, 0.001 * thresh, coulop, fock_json);

    // 2. Map the input deck → a Tier-A plan (R3b). requested_properties +
    //    beta.*/excited.* knobs select which ResponsePropertyRequests to build;
    //    merge_plans dedupes shared FD states across them.
    std::vector<char> axes;
    for (char c : rp.dipole_directions()) {
      const char l = static_cast<char>(std::tolower(c));
      if (l == 'x' || l == 'y' || l == 'z') axes.push_back(l);
    }
    if (axes.empty()) axes = {'x', 'y', 'z'};
    const std::vector<double> freqs = rp.dipole_frequencies();
    const auto props = rp.requested_properties();
    auto wants = [&](const char *p) {
      return std::find(props.begin(), props.end(), std::string(p)) != props.end();
    };

    std::vector<ResponsePlan> plans;
    auto add = [&](ResponsePropertyRequest r) {
      r.axes = axes;
      r.protocol_thresholds = protocol;
      plans.push_back(plan_one(r));
    };
    if (wants("polarizability")) {
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::Polarizability;
      r.frequencies = freqs;
      add(r);
    }
    // `quadratic true` is the legacy (v2-era) spelling of "compute beta" —
    // honor it alongside requested_properties so old decks keep their
    // hyperpolarizability instead of silently dropping it (M3 golden regen
    // caught this: the alpha+beta deck produced alpha-only output).
    if (wants("hyperpolarizability") || rp.quadratic()) {
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::Hyperpolarizability;
      r.beta_process = rp.beta_or() ? BetaProcess::OR : BetaProcess::SHG;
      r.frequencies = freqs;
      add(r);
    }
    if (wants("raman")) {
      // SINGLE-COMPONENT vibrational Raman (atom 0, z) — v3's full per-atom
      // tensor is deferred (post-state-parallel), so this won't match v2's full
      // Raman; it exercises the β(dipole;dipole,nuclear) path.
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::PolarizabilityGradient;
      r.gradient_mode = GradientMode::Nuclear;
      r.frequencies = freqs;
      r.raman_nuc_atom = 0;
      r.raman_nuc_axis = 2;
      add(r);
    }
    if (rp.excited_enable()) {
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::PolarizabilityGradient;
      r.gradient_mode = GradientMode::Resonant;
      r.n_roots = static_cast<int>(rp.excited_num_states());
      add(r);
    }
    if (plans.empty()) {  // default: polarizability
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::Polarizability;
      r.frequencies = freqs;
      add(r);
    }

    // 3. Build the workflow input + settings; run the core with the ground
    //    state already loaded. archive_file is still passed through: the GS
    //    fingerprint gate hashes it (restart safety), and the F2 subworld
    //    fan-out needs it to load per-subworld ground states.
    ResponseWorkflowInput in;
    in.archive_file = archive;
    in.protocols = protocol;
    in.plan = merge_plans(plans);
    // excited.tda=false → Full (X,Y) ES bundle (default TDA).
    if (rp.excited_enable() && !rp.excited_tda())
      for (auto &e : in.plan.es) e.tda = false;
    // ResponseApplication::run has already chdir'd (ScopedCWD) into `outdir`, and
    // `outdir` is RELATIVE — so the calc dir is the cwd ("."). Using outdir here
    // would double the path (outdir/outdir) and the metadata would be written/read
    // in different places, leaving Output.properties empty. (outdir is still used
    // above to resolve the ground archive relative to this cwd.)
    in.settings.calc_dir = ".";
    in.settings.max_iters = static_cast<int>(rp.maxiter());
    // Deck `subworlds N` -> the F2 state-parallel fan-out (same path as the
    // standalone --fd-subworlds flag; archive_file above makes it live).
    in.settings.fd_subworlds = std::max(0, rp.subworlds());
    if (world.rank() == 0 && in.settings.fd_subworlds > 0) {
      print("response: deck subworlds =", in.settings.fd_subworlds,
            "(F2 state-parallel fan-out requested)");
      // Review io HIGH (early warning, not a hard gate): on a MULTI-NODE run
      // the subworlds write native archives at subworld-rank-count, but
      // property assembly reloads them at universe scale — the native np-guard
      // then aborts AFTER the full solve. Warn up front so the user isn't
      // surprised late; the HDF5 backend gathers to one client and reloads at
      // any np. (Single-node subworlds, where writer==reader np, are fine.)
#ifdef MADNESS_HAS_HDF5
      const bool hdf5_on = hdf5_io_enabled();
#else
      const bool hdf5_on = false;
#endif
      if (!hdf5_on)
        print("response: WARNING — subworlds with the NATIVE backend will hit "
              "the cross-np archive guard at property assembly on a multi-node "
              "run (states saved per-subworld, reloaded at universe scale). Use "
              "`response { io { backend hdf5 } }` for multi-node subworld runs, "
              "or set MADRESPONSE_ALLOW_NP_MISMATCH=1 if the archives are known "
              "portable. Single-node subworld runs are unaffected.");
    }
    in.settings.policy.dconv_user = rp.dconv();
    in.settings.print_level =
        static_cast<PrintLevel>(std::max(0, std::min(3, rp.print_level())));
    // Deck-level HDF5 opt-in (response.hdf5 true) — the env var
    // MADRESPONSE_IO_HDF5 still works; the deck parameter wins when set.
    if (rp.hdf5()) {
#ifdef MADNESS_HAS_HDF5
      set_hdf5_io_enabled(true);
#else
      throw std::runtime_error(
          "response.hdf5 requested but this build has no HDF5 support — "
          "configure with -DMADNESS_ENABLE_HDF5=ON");
#endif
    }

    // Deck `dalton.dir <path>` — the seed-from-directory import contract
    // (showcase W3; import-only, madness never invokes DALTON). Runs BEFORE
    // the workflow so reconcile sees the seed bundles as restart sources.
    // The GS is prepared at protocol.front() above (the documented
    // precondition); geometry-fingerprint / frequency mismatches throw here,
    // failing the response task loudly instead of silently solving cold.
    if (!rp.dalton_dir().empty()) {
      run_dalton_import(world, gs, scf_calc->molecule, in.plan,
                        in.settings.calc_dir, rp.dalton_dir());
    }

    // Pass the RESOLVED fock path through (review finding: "" here made every
    // per-rung re-prepare and every subworld GS reload RECOMPUTE the Fock
    // instead of loading moldft's per-protocol entries — silently diverging
    // from the standalone CLI path this adapter is supposed to match).
    ResponseWorkflowOutput out =
        run_response_with_ground(world, gs, L, fock_json, in);

    // 4. Map Output → Results. Stash v3 timing/diagnostics under metadata so
    //    they ride into the workflow's calc_info.json.
    Results res;
    res.metadata = std::move(out.metadata);
    res.properties = std::move(out.properties);
    if (world.rank() == 0) {
      res.metadata["engine"]         = "molresponse_v3";
      res.metadata["v3_timing"]      = std::move(out.timing);
      res.metadata["v3_diagnostics"] = std::move(out.diagnostics);
    }
    return res;
  }
};

#endif // MOLRESPONSE_V3_MADQC_ADAPTER_HPP
