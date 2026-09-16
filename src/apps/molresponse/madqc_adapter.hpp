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
#include <apps/molresponse/solvers/dalton_gs_seed.hpp>   // GS seed from molden (dalton.dir)

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

namespace molresponse_v3 {

/// Run-wide seed directory: `io.dalton.dir` (ParameterManager IOParameters)
/// wins; the response-block `dalton.dir` is the alias. Same rule for HDF5:
/// `io.backend hdf5` or the alias `response.hdf5`.
template <typename Params>
inline std::string effective_dalton_dir(const Params &params) {
  const std::string io_dir = params.template get<IOParameters>().dalton_dir();
  if (!io_dir.empty()) return io_dir;
  return params.template get<ResponseParameters>().dalton_dir();
}
template <typename Params>
inline bool effective_hdf5(const Params &params) {
  return params.template get<IOParameters>().hdf5() ||
         params.template get<ResponseParameters>().hdf5();
}

// ---------------------------------------------------------------------------
// GROUND-STATE SEED HOOK (2026-09-09, seeding showcase). Installed by madqc on
// the SCF application (SCFApplication::set_pre_run_hook) when the deck carries
// `dalton.dir`. Runs collectively inside the SCF work dir right before the
// engine plans its restart:
//   * no-op when an archive <prefix>.restartdata* already exists (a real
//     restart always wins over a seed) or when `restart` is not auto/iterate;
//   * locates molden.inp in dalton.dir (same resolver as the FD/ES import),
//     fingerprints its geometry against the deck's Molecule (hard error);
//   * n_occ = (sum Z - charge)/2 (closed shell), L = deck `l`,
//     thresh = first protocol rung; writes <prefix>.restartdata via
//     write_gs_seed_from_molden, so RestartPlan(auto) -> restartdata/iterate.
// The projection happens at the deck's box; SCF re-sets its own
// FunctionDefaults afterwards, and load_mos re-projects k/thresh as needed.
// ---------------------------------------------------------------------------
inline void seed_gs_from_dalton_dir(World &world, const Params &params,
                                    const std::filesystem::path &workdir,
                                    const std::string &dalton_dir) {
  namespace fs = std::filesystem;
  const auto &cp  = params.get<CalculationParameters>();
  const auto &mol = params.get<Molecule>();
  const std::string prefix = cp.prefix();
  const std::string mode   = cp.restart();
  if (mode != "auto" && mode != "iterate") {
    if (world.rank() == 0)
      print("[DALTON-SEED] GS: restart =", mode, "-> not seeding the ground state");
    return;
  }
  // Existing archive wins (rank 0 decides, collective broadcast).
  int have = 0;
  if (world.rank() == 0) {
    for (const auto &e : fs::directory_iterator(workdir)) {
      const std::string n = e.path().filename().string();
      if (n.rfind(prefix + ".restartdata", 0) == 0) { have = 1; break; }
    }
  }
  world.gop.broadcast(have, 0);
  if (have) {
    if (world.rank() == 0)
      print("[DALTON-SEED] GS:", prefix + ".restartdata*", "already present in",
            workdir.string(), "-> restart from it, not from the DALTON seed");
    return;
  }
  // Locate + fingerprint (rank 0), broadcast the molden path or the error.
  std::string molden, err, report;
  if (world.rank() == 0) {
    try {
      auto m = locate_dalton_dir(dalton_dir, (workdir / "dalton_import").string(),
                                 "", "", "");
      auto check = fingerprint_dalton_geometry(m, mol, 1e-4);
      report = check.report;
      if (!check.ok)
        throw std::runtime_error(
            "dalton import (GS seed): GEOMETRY FINGERPRINT MISMATCH\n" + check.report);
      molden = m.molden_path;
    } catch (const std::exception &ex) { err = ex.what(); }
  }
  world.gop.broadcast_serializable(err, 0);
  if (!err.empty()) throw std::runtime_error(err);
  world.gop.broadcast_serializable(molden, 0);
  world.gop.broadcast_serializable(report, 0);

  const double Z     = mol.total_nuclear_charge();
  const double nelec = Z - cp.charge();
  const long   ne    = std::lround(nelec);
  if (std::abs(nelec - static_cast<double>(ne)) > 1e-6 || ne % 2 != 0)
    throw std::runtime_error("dalton import (GS seed): closed-shell seed needs an even "
                             "electron count, got " + std::to_string(nelec));
  GsSeedOptions opt;
  opt.L        = cp.L();
  opt.thresh   = cp.protocol().empty() ? 1e-4 : cp.protocol().front();
  opt.xc       = cp.get<std::string>("xc");
  opt.localize = cp.get<std::string>("localize");
  opt.extra_prefixes = {prefix + ".gs_seed"};   // preserved copy (save_mos overwrites <prefix>.restartdata)
  opt.active_molecule = &mol;   // RestartPlan matches geometry at 1e-8 and eprec exactly
  if (world.rank() == 0) {
    print("[DALTON-SEED] GS: metadata molecule (stamped) eprec =", mol.parameters.eprec());
    for (std::size_t i = 0; i < mol.natom(); ++i) {
      const auto at = mol.get_atom(i);
      printf("[DALTON-SEED] GS:   atom %zu  Z=%d  %.10f %.10f %.10f\n", i, at.atomic_number, at.x, at.y, at.z);
    }
  }
  if (world.rank() == 0) {
    print("[DALTON-SEED] GS: seeding", prefix + ".restartdata", "from", molden,
          " n_occ =", ne / 2, " L =", opt.L, " thresh =", opt.thresh);
    print(report);
  }
  auto rep = write_gs_seed_from_molden(world, molden, static_cast<int>(ne / 2), prefix, opt);
  if (world.rank() == 0) {
    nlohmann::json j;
    j["seed"] = "dalton_import"; j["stage"] = "ground_state";
    j["molden"] = molden; j["dalton_dir"] = dalton_dir;
    j["n_occ"] = rep.n_occ; j["n_ao"] = rep.n_ao; j["n_mo"] = rep.n_mo;
    j["L"] = opt.L; j["thresh"] = opt.thresh;
    j["max_offdiag_pre_loewdin"] = rep.max_offdiag_pre;
    j["archive"] = rep.archive;
    j["preserved_copy"] = prefix + ".gs_seed.restartdata";
    std::ofstream out((workdir / (prefix + ".gs_seed.json")).string());
    out << j.dump(2) << "\n";
  }
}
} // namespace molresponse_v3

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
    const std::string dalton_dir = effective_dalton_dir(params);
    if (apply_seed_start_rung(protocol, rp.seed_start_rung(),
                              !dalton_dir.empty())) {

      if (world.rank() == 0)
        print("response: seed.start_rung=fine — dalton.dir seed starts the "
              "ladder at thresh", protocol.front());
    } else if (rp.seed_start_rung() == "fine" && dalton_dir.empty() &&
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
      // Deck `beta.frequencies` restricts the driver frequencies (and thus the
      // 2w legs) to the ones actually wanted; empty = the whole dipole grid.
      const auto bf = rp.beta_frequencies();
      r.frequencies = bf.empty() ? freqs : bf;
      add(r);
    }
    if (wants("raman")) {
      // ONE nuclear coordinate per run (deck: raman.nuc_atom / raman.nuc_axis;
      // default atom 0, z as before). The run emits the whole dipole x dipole
      // block for that coordinate at every requested frequency, so a full 3N
      // Cartesian polarizability gradient is 3N runs over one calculation
      // directory: the dipole responses are computed once and reloaded. The
      // planner's full per-atom sentinel is still deferred (post-state-parallel).
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::PolarizabilityGradient;
      r.gradient_mode = GradientMode::Nuclear;
      r.frequencies = freqs;
      r.raman_nuc_atom = rp.raman_nuc_atom();
      r.raman_nuc_axis = rp.raman_nuc_axis();
      if (r.raman_nuc_atom < 0 || r.raman_nuc_axis < 0 || r.raman_nuc_axis > 2)
        throw std::runtime_error(
            "response: raman.nuc_atom must be >= 0 and raman.nuc_axis in 0..2 — "
            "the full per-atom Raman tensor is not implemented; drive the 3N "
            "coordinates as 3N runs over one calc dir");
      add(r);
    }
    if (rp.excited_enable()) {
      ResponsePropertyRequest r;
      r.kind = ResponsePropertyKind::PolarizabilityGradient;
      r.gradient_mode = GradientMode::Resonant;
      r.n_roots = static_cast<int>(rp.excited_num_states());
      // Roots only unless the deck also asks for the two-photon contraction
      // (excited.tpa) — the derived FD legs serve only the 2PA residue.
      r.tpa = rp.excited_tpa();
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
    // Deck `subworld_ranks R` -> universe-level split into subworlds of R ranks
    // that may span nodes (large-system regime); 0 keeps per-node packing.
    in.settings.fd_subworld_ranks = std::max(0, rp.subworld_ranks());
    // Deck `dalton.dir` + `seed.freq_tol` -> nearest-frequency DALTON guess for
    // the derived (two-photon) FD legs (calc_executor solve_fd seam).
    in.settings.dalton_dir    = dalton_dir;
    in.settings.seed_freq_tol = rp.seed_freq_tol();
    in.settings.es_seed_warmup = rp.seed_es_warmup();
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
        print("response: NOTE — subworlds with the NATIVE backend: per-subworld "
              "archives are treated as np-portable by default at property "
              "assembly (nio=1); set MADRESPONSE_STRICT_NP=1 to hard-block a "
              "cross-np reload instead. `response { io { backend hdf5 } }` "
              "sidesteps the question entirely on multi-node runs.");
    }
    in.settings.policy.dconv_user = rp.dconv();
    // Deck overrides for the iteration policy (2026-09-10): the response block's
    // `kain` / `maxrotn` reach the FD/ES solvers only when the deck SETS them.
    // The deck default kain=false predates the solvers' KAIN-on default;
    // honouring it unconditionally would switch KAIN off for every run.
    if (rp.is_user_defined("kain"))    in.settings.policy.kain    = rp.kain();
    if (rp.is_user_defined("maxrotn")) in.settings.policy.maxrotn = rp.maxrotn();
    if (rp.is_user_defined("kain.min_residual"))
      in.settings.policy.kain_min_residual = rp.kain_min_residual();
    // Plateau detector (deck defaults equal the policy defaults: window 6, ratio 0.1).
    if (rp.is_user_defined("stall.window")) in.settings.policy.stall_window = rp.stall_window();
    if (rp.is_user_defined("stall.ratio"))  in.settings.policy.stall_ratio  = rp.stall_ratio();
    in.settings.print_level =
        static_cast<PrintLevel>(std::max(0, std::min(3, rp.print_level())));
    // ES initial-guess knobs (deck: response { excited.guess virtual_ao,
    // excited.guess_basis aug-cc-pvtz }). Deck default (solid_harmonics /
    // aug-cc-pvdz) maps to the ExecutorSettings default — existing decks are
    // unchanged. virtual_ao is the energy-ordered AO-virtual guess; it is
    // required to reach totally-symmetric / radially-excited states on atoms,
    // which the angular-only solid-harmonic trials structurally cannot span.
    in.settings.es_guess       = parse_es_guess_mode(rp.excited_guess());
    in.settings.es_guess_basis = rp.excited_guess_basis();
    if (world.rank() == 0 && rp.excited_enable() &&
        in.settings.es_guess != ESGuessMode::SolidHarmonics)
      print("response: excited.guess =", to_string(in.settings.es_guess),
            (in.settings.es_guess == ESGuessMode::VirtualAO
                 ? std::string("(basis " + in.settings.es_guess_basis + ")")
                 : std::string()));
    // Deck-level HDF5 opt-in (response.hdf5 true) — the env var
    // MADRESPONSE_IO_HDF5 still works; the deck parameter wins when set.
    if (effective_hdf5(params)) {
#ifdef MADNESS_HAS_HDF5
      set_hdf5_io_enabled(true);
#else
      throw std::runtime_error(
          "response.hdf5 requested but this build has no HDF5 support — "
          "configure with -DMADNESS_ENABLE_HDF5=ON");
#endif
    }

#ifdef MADNESS_HAS_HDF5
    // Ground state in HDF5 (2026-09-09): with the HDF5 restart opt-in, mirror
    // moldft's native archive as <prefix>.restartdata.h5 (same stream, one
    // blob) so a calc dir can carry GS + response + ES in HDF5. Native stays
    // authoritative for moldft; GroundState::from_archive reads either.
    // Always (re)write: the DALTON GS seed hook writes <prefix>.restartdata.h5
    // BEFORE moldft runs (its HDF5 twin of the seed), so an existence check
    // here kept the unconverged seed as the .h5 GS (closeout attempt 9, 2026-09-10:
    // seed_h2o.restartdata.h5 == the 12:01:12 seed twin, native at 12:01:41).
    // Native was still preferred by from_archive, so nothing was mis-loaded,
    // but a native-less consumer would have read the seed. The preserved twin
    // lives under <prefix>.gs_seed.restartdata.h5.
    if (hdf5_io_enabled()) gs.save_archive_hdf5(world, archive + ".h5");
#endif

    // Deck `dalton.dir <path>` — the seed-from-directory import contract
    // (showcase W3; import-only, madness never invokes DALTON). Runs BEFORE
    // the workflow so reconcile sees the seed bundles as restart sources.
    // The GS is prepared at protocol.front() above (the documented
    // precondition); geometry-fingerprint / frequency mismatches throw here,
    // failing the response task loudly instead of silently solving cold.
    if (!dalton_dir.empty()) {
      run_dalton_import(world, gs, scf_calc->molecule, in.plan,
                        in.settings.calc_dir, dalton_dir, {}, {}, {}, 1e-4,
                        /*es_y_from_dalton=*/rp.seed_es_y() == "dalton");
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
