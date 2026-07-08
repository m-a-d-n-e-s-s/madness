// =========================================================================
// test_v3_es_skeleton.cpp — validates the new architecture skeleton
//
//   kernels/tags.hpp          (response-type + shell tags + traits)
//   kernels/tda.hpp           (Kernels<TDA, ClosedShell>)
//   kernels/response_space_ops.hpp  (rs::inner / rs::diagonalize /
//                                    rs::transform)
//   solvers/response_state.hpp      (ResponseStateX/XY × Shell)
//   solvers/iterate.hpp             (driver template)
//   solvers/es_solver.hpp           (ESSolver<Type, Shell>)
//   solvers/build_response_ground_state.hpp (build helpers)
//
// Closed-shell-TDA only. Reuses the H2 fixture and reference omegas
// from test_es_solver.cpp.
// =========================================================================

#include "../GroundState.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/tda.hpp"
#include "../solvers/build_response_ground_state.hpp"
#include "../solvers/convergence_policy.hpp"
#include "../solvers/es_save_load.hpp"
#include "../solvers/es_solver.hpp"
#include "../solvers/iterate.hpp"
#include "../solvers/iterate_protocol.hpp"
#include "../solvers/response_state.hpp"

#include <madness/chem/SCFOperators.h>
#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <madness/world/worldmem.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

// Reference omegas for sanity-checking v3 ESSolver<TDA, ClosedShell>.
//
// Two references per molecule when available:
//
//   legacy_tda  — molresponse_legacy at the same protocol, same basis
//                 (MRA). Tight: pass-fail comparison.
//   dalton_rpa  — Dalton RPA at the given Gaussian basis. Informative
//                 only (RPA != TDA but Y-amplitudes are small for the
//                 lowest singlets of these closed-shell systems, so the
//                 omegas agree to a few meV).
//
// Sourced from
//   madness_studies/es_benchmarks/data/<system>/run_legacy_tda/run.log
//   madness_studies/refs/dalton_tdhf.json
struct ESRef {
  const char* name;
  long n_alpha;
  std::vector<double> legacy_tda;      // empty if not available
  std::vector<double> dalton_rpa;      // empty if not available
  const char*         dalton_rpa_basis = "";
};

static const ESRef ES_REFS[] = {
  // H2 RHF (legacy from test_es_solver.cpp:42-43, k=8).
  // DALTON RPA at d-aug-cc-pVQZ (madness_studies/refs/dalton_tdhf.json:h2.rpa).
  {"H2",  1,
   {0.46535391, 0.47696774, 0.48065960, 0.48069615},
   {0.46537636, 0.47706849, 0.48088936, 0.48088936},
   "d-aug-cc-pVQZ"},
  // H2O RHF.
  // Legacy TDA: madness_studies/es_benchmarks/h2o_1942685/run_legacy_tda.
  // DALTON RPA at d-aug-cc-pVQZ (madness_studies/refs/dalton_tdhf.json:h2o.rpa).
  // NB: NOT aug-cc-pVDZ — the lowest roots include diffuse/Rydberg states
  // (esp. root 3) that need double augmentation to reach the basis-set
  // limit that complete-basis MRA is compared against. Against aug-cc-pVDZ
  // root 3 looks ~3.5e-2 off; against d-aug-cc-pVQZ it is ~8e-5.
  {"H2O", 5,
   {0.31953104, 0.38160118, 0.40059069, 0.42073252},
   {0.31778528, 0.37897397, 0.39965329, 0.40973190},
   "d-aug-cc-pVQZ"},
  // C2H4 RHF.
  // Legacy TDA: madness_studies/es_benchmarks/c2h4_1942683/run_legacy_tda.
  // DALTON RPA at d-aug-cc-pVQZ (refs/dalton_tdhf.json:c2h4.rpa).
  {"C2H4", 8,
   {0.26122188, 0.28272418, 0.28824251, 0.29641883},
   {0.26060565, 0.28217118, 0.28801612, 0.31409892},
   "d-aug-cc-pVQZ"},
  // ---- Li UHF (atomic Li doublet) ----
  // n_alpha = 2 (1s↑, 2s↑), n_beta = 1 (1s↓). Lowest excitation is
  // 2s↑→2p↑ around 0.068 au (~1.85 eV experimentally) — same energy
  // where Dalton's α diverged at ω=0.0656 in the FD benchmark.
  //
  // DALTON RPA at aug-cc-pVQZ for the Li doublet: PENDING — run
  //   python refs/compute_dalton_tdhf.py li --basis aug-cc-pVQZ --states 4
  // and paste omegas_au here once that script supports UHF.
  {"Li", 2, {}, {}, "aug-cc-pVQZ (pending)"},
};

static const ESRef* find_ref(long n_alpha) {
  for (const auto& r : ES_REFS) if (r.n_alpha == n_alpha) return &r;
  return nullptr;
}

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);

    if (world.rank() == 0) {
      print("\n=========================================================");
      print(" molresponse_v3 ES Skeleton Test  (TDA × ClosedShell)");
      print("=========================================================\n");
    }
    // Inner scope so every MADNESS-allocated object (solver, target,
    // state, GroundState, coulop, response functions...) destructs
    // BEFORE finalize(). Without this, destructors fire during world
    // teardown and hit destroyed mutexes inside the convolution /
    // function caches — produces a wall of "RecursiveMutex::lock()
    // failed acquiring mutex" lines after the run is otherwise done.
    {

    commandlineparser parser(argc, argv);
    if (!parser.key_exists("archive")) {
      if (world.rank() == 0) {
        print("Usage: test_v3_es_skeleton --archive=<path> "
              "[--num-roots=N] [--maxiter=N] [--dconv=X] "
              "[--thresh=X | --protocol=th1,th2,...] [--k=N] "
              "[--tol=X] [--print-level=0..3] "
              "[--log-convergence=PATH] "
              "[--stream-theta] [--es-batch] [--es-tensor] [--es-time] "
              "[--tda-warmup-iters=N] [--cluster-unmix-factor=F] "
              "[--warmup-oversample=K] "
              "[--save-roots=DIR] [--load-roots=DIR] "
              "[--load-roots-tda=DIR] "
              "[--load-best=CALC_DIR] "
              "[--type=tda|full] "
              "[--guess=random|solid_harmonics] "
              "[--no-kain] [--kain-maxsub=N] [--maxrotn=X] "
              "[--kain-cmax=X]");
      }
      finalize();
      return 1;
    }

    std::string archive_path = parser.value_raw("archive");
    const long num_roots = parser.key_exists("num-roots")
                               ? std::stol(parser.value("num-roots"))
                               : 4;
    const int max_iters = parser.key_exists("maxiter")
                              ? std::stoi(parser.value("maxiter"))
                              : 25;
    const double dconv =
        parser.key_exists("dconv") ? std::stod(parser.value("dconv")) : 1e-4;
    const double tol = parser.key_exists("tol")
                           ? std::stod(parser.value("tol"))
                           : 0.005;
    const int pl_int = parser.key_exists("print-level")
                           ? std::stoi(parser.value("print-level"))
                           : 1;
    const PrintLevel print_level = static_cast<PrintLevel>(
        std::max(0, std::min(3, pl_int)));
    const std::string log_path = parser.key_exists("log-convergence")
                                     ? parser.value_raw("log-convergence")
                                     : std::string();
    const bool stream_theta = parser.key_exists("stream-theta");
    // Stage-1 bundle batching (closed-shell TDA): V0x/T0x/BSH over the
    // flattened M·n_occ bundle in one pass. Default off → per-root reference.
    const bool es_batch = parser.key_exists("es-batch");
    // Inc-3c tensor-layer γ (closed-shell TDA): g0 cached per protocol +
    // batched Coulomb + per-root contractions. A/B-to-floor vs the reference
    // (~1e-3 REL) — a SEPARATE gate from the bit-identical --es-batch.
    const bool es_tensor = parser.key_exists("es-tensor");
    // Per-phase wall-time instrumentation: one ES_TIMING line per step.
    const bool es_time = parser.key_exists("es-time");
    const int tda_warmup_iters = parser.key_exists("tda-warmup-iters")
                                     ? std::stoi(parser.value("tda-warmup-iters"))
                                     : 3;
    const double cluster_unmix_factor =
        parser.key_exists("cluster-unmix-factor")
            ? std::stod(parser.value("cluster-unmix-factor"))
            : 100.0;
    const double warmup_oversample = parser.key_exists("warmup-oversample")
                                          ? std::stod(parser.value("warmup-oversample"))
                                          : 2.0;
    const bool no_kain = parser.key_exists("no-kain");
    const int kain_maxsub = parser.key_exists("kain-maxsub")
                                ? std::stoi(parser.value("kain-maxsub"))
                                : 10;
    const double maxrotn = parser.key_exists("maxrotn")
                               ? std::stod(parser.value("maxrotn"))
                               : 0.5;
    const double kain_cmax = parser.key_exists("kain-cmax")
                                 ? std::stod(parser.value("kain-cmax"))
                                 : 100.0;
    const std::string save_roots_dir = parser.key_exists("save-roots")
                                           ? parser.value_raw("save-roots")
                                           : std::string();
    const std::string load_best_dir = parser.key_exists("load-best")
                                          ? parser.value_raw("load-best")
                                          : std::string();
    const std::string load_roots_dir = parser.key_exists("load-roots")
                                           ? parser.value_raw("load-roots")
                                           : std::string();
    const std::string load_roots_tda_dir =
        parser.key_exists("load-roots-tda")
            ? parser.value_raw("load-roots-tda") : std::string();
    // ES iteration variant: "tda" (default) or "full" (RPA paired X,Y).
    // Full is ClosedShell-only today; open-shell + full is rejected
    // below before the dispatch.
    const std::string es_type = parser.key_exists("type")
                                    ? parser.value("type")
                                    : std::string("tda");
    if (es_type != "tda" && es_type != "full") {
      if (world.rank() == 0)
        print("ERROR: --type must be 'tda' or 'full', got '",
              es_type, "'");
      finalize();
      return 1;
    }
    const ESGuessMode guess_mode = parser.key_exists("guess")
                                       ? parse_es_guess_mode(parser.value("guess"))
                                       : ESGuessMode::SolidHarmonics;

    auto header = GroundState::read_archive_header(world, archive_path);
    int override_k = parser.key_exists("k") ? std::stoi(parser.value("k")) : -1;

    // Protocol ladder: --protocol=th1,th2,... wins; else fall back to
    // single-step --thresh=X; else default_thresh_for_k(header.k).
    std::vector<double> protocol_thresholds;
    if (parser.key_exists("protocol")) {
      std::string csv = parser.value("protocol");
      size_t pos = 0;
      while (pos < csv.size()) {
        size_t comma = csv.find(',', pos);
        std::string tok = csv.substr(pos, comma - pos);
        if (!tok.empty()) protocol_thresholds.push_back(std::stod(tok));
        if (comma == std::string::npos) break;
        pos = comma + 1;
      }
    } else {
      double single = parser.key_exists("thresh")
                          ? std::stod(parser.value("thresh"))
                          : default_thresh_for_k(header.k);
      protocol_thresholds.push_back(single);
    }
    // Set the first protocol step now so subsequent ground-state
    // prep + initial-guess construction sees the right FunctionDefaults.
    set_response_protocol(world, header.L, protocol_thresholds.front(),
                          override_k);

    Molecule molecule;
    auto archive_dir = std::filesystem::path(archive_path).parent_path();
    for (const auto &name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
      auto candidate = archive_dir / name;
      if (std::filesystem::exists(candidate)) {
        std::ifstream ifs(candidate);
        nlohmann::json j;
        ifs >> j;
        nlohmann::json mol_json;
        if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
          mol_json = j["tasks"][0]["molecule"];
        else if (j.contains("molecule"))
          mol_json = j["molecule"];
        if (!mol_json.is_null()) molecule.from_json(mol_json);
        break;
      }
    }

    auto gs = GroundState::from_archive(world, archive_path, molecule);
    if (world.rank() == 0) gs.print_info();

    std::string fock_json;
    for (const auto &name : {"moldft.fock.json", "mad.fock.json"}) {
      auto candidate = archive_dir / name;
      if (std::filesystem::exists(candidate)) {
        fock_json = candidate.string();
        break;
      }
    }
    const double cur_thresh = FunctionDefaults<3>::get_thresh();
    auto coulop = poperatorT(
        CoulombOperatorPtr(world, gs.params().lo(), 0.001 * cur_thresh));
    gs.prepare(world, 0.001 * cur_thresh, coulop, fock_json);

    const bool restricted = gs.is_spin_restricted();
    if (es_type == "full" && !restricted) {
      if (world.rank() == 0)
        print("ERROR: --type=", es_type, " requires a closed-shell "
              "(restricted) ground state. Open-shell Full ES is not yet "
              "implemented.");
      finalize();
      return 1;
    }
    const ESRef* ref = find_ref(gs.num_alpha());
    if (!ref) {
      if (world.rank() == 0) {
        print("\nNo reference omegas tabulated for num_alpha=",
              gs.num_alpha(),
              " — running anyway, comparisons will be skipped.");
      }
    } else if (world.rank() == 0) {
      print("\nReference data for ", ref->name,
            " (num_alpha=", ref->n_alpha, "):");
      if (!ref->legacy_tda.empty()) {
        print("  legacy molresponse TDA  (same protocol):");
        for (size_t i = 0; i < ref->legacy_tda.size(); ++i)
          print("    omega[", i, "] =", ref->legacy_tda[i]);
      }
      if (!ref->dalton_rpa.empty()) {
        print("  DALTON RPA  basis=", ref->dalton_rpa_basis, ":");
        for (size_t i = 0; i < ref->dalton_rpa.size(); ++i)
          print("    omega[", i, "] =", ref->dalton_rpa[i]);
      }
    }

    // -------------------------------------------------------------------
    // Drive the new ESSolver<TDA, *> via iterate_protocol. Both shells
    // produce a `final omegas + final residuals` summary that the
    // comparison block below reads.
    // -------------------------------------------------------------------
    madness::Tensor<double> final_omega;
    std::vector<double>     final_residuals;
    // Convergence gate inputs (captured per branch alongside the omegas).
    // A run is only a PASS if it actually CONVERGED — without this the
    // ω-tolerance check alone passes a run that stalled at maxiter with a
    // residual floor (exactly how the indefinite-metric bug hid as green).
    bool   final_diverged        = false;
    double final_residual_target = 0.0;  // solver's own bsh_residual target

    ConvergencePolicy policy;
    policy.dconv_user           = dconv;
    policy.stream_theta         = stream_theta;
    policy.tda_warmup_iters     = tda_warmup_iters;
    policy.cluster_unmix_factor = cluster_unmix_factor;
    policy.warmup_oversample_factor = warmup_oversample;
    policy.kain                 = !no_kain;
    policy.kain_maxsub          = kain_maxsub;
    policy.maxrotn              = maxrotn;
    policy.kain_cmax_cap        = kain_cmax;
    if (world.rank() == 0) {
      print("  stream_theta =", stream_theta ? "true" : "false",
            "  ←  step_", (stream_theta ? "recompute_pieces" : "rotate_pieces"));
      print("  tda_warmup_iters =", tda_warmup_iters,
            "  warmup_oversample =", warmup_oversample,
            "  cluster_unmix_factor =", cluster_unmix_factor,
            "  guess =", to_string(guess_mode));
      print("  kain =", policy.kain ? "on" : "off",
            "  kain_maxsub =", kain_maxsub,
            "  maxrotn =", maxrotn,
            "  kain_cmax =", kain_cmax);
    }
    const long n_roots_warmup = std::max<long>(
        num_roots,
        static_cast<long>(std::ceil(warmup_oversample *
                                     static_cast<double>(num_roots))));

    // After the explicit warmup phase, the main solver should NOT
    // re-run a warmup window — it's starting from a cleaned state
    // and KAIN should engage from iter 1.
    ConvergencePolicy main_policy = policy;
    main_policy.tda_warmup_iters  = 0;

    if (world.rank() == 0) {
      print("\n=== iterate_protocol(ESSolver<TDA, ",
            (restricted ? "ClosedShell" : "OpenShell"), ">, ...) ===");
      print("  num_roots  = ", num_roots);
      print("  dconv_user = ", dconv);
      print("  max_iters  = ", max_iters, "  per protocol step");
      print("  protocol   = ", protocol_thresholds);
    }

    solvers::IterateProtocolPolicy proto_policy;
    proto_policy.max_iters_per_step = max_iters;

    if (restricted && es_type == "full") {
      // --- Full (paired X,Y) RPA path, ClosedShell only ---
      using Solver = ESSolver<Full, ClosedShell>;
      auto problem = build_es_problem_full<ClosedShell>(
          world, gs, num_roots, /*c_xc=*/1.0, gs.params().lo());

      // Initial state: precedence is (a) same-type Full archive,
      // (b) TDA archive promoted to Full with Y=0, (c) fresh TDA
      // warmup promoted to Full with Y=0.
      Solver::State state0;
      if (!load_roots_dir.empty()) {
        state0 = load_es_roots<Full, ClosedShell>(world, load_roots_dir);
      } else if (!load_roots_tda_dir.empty()) {
        auto tda_state =
            load_es_roots<TDA, ClosedShell>(world, load_roots_tda_dir);
        state0 = promote_tda_to_full_closed_shell(world, tda_state);
      } else {
        auto tda_state = run_oversampled_tda_warmup<ClosedShell>(
            world, gs, num_roots, n_roots_warmup, tda_warmup_iters,
            policy, /*c_xc=*/1.0, gs.params().lo(), print_level,
            guess_mode);
        state0 = promote_tda_to_full_closed_shell(world, tda_state);
      }
      Solver solver(world, std::move(problem), main_policy, print_level);
      if (!log_path.empty()) solver.set_log_path(log_path);

      auto prepare = [&](double thresh, Solver &solv, Solver::State &st) {
        set_response_protocol(world, header.L, thresh, override_k);
        const int    new_k = FunctionDefaults<3>::get_k();
        const double new_t = FunctionDefaults<3>::get_thresh();
        if (world.rank() == 0)
          print("\n--- protocol step: thresh =", new_t,
                "  k =", new_k, " ---");
        auto coulop_new = poperatorT(
            CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
        gs.prepare(world, 0.001 * new_t, coulop_new, fock_json);
        solv.set_gs(build_response_ground_state_closed_shell(
            world, gs, /*c_xc=*/1.0, gs.params().lo()));
        for (auto &root : st.roots) {
          for (auto &fn : root.x_alpha)
            fn = madness::project(fn, new_k, new_t);
          for (auto &fn : root.y_alpha)
            fn = madness::project(fn, new_k, new_t);
        }
        for (auto &rho : st.rho_alpha_prev)
          rho = madness::project(rho, new_k, new_t);
      };
      auto sf = solvers::iterate_protocol(
          solver, state0, protocol_thresholds, prepare, proto_policy);
      final_omega     = sf.omega;
      final_residuals = sf.last_bsh_residual;
      final_diverged        = sf.diverged;
      final_residual_target = solver.targets().bsh_residual;
      if (!save_roots_dir.empty()) {
        save_es_roots<Full, ClosedShell>(world, sf, save_roots_dir,
                                          /*converged=*/!sf.diverged);
        if (world.rank() == 0)
          print("[SAVE] es_roots: wrote", num_roots,
                "root(s) +/", save_roots_dir, "/roots.json");
      }
    } else if (restricted) {
      using Solver = ESSolver<TDA, ClosedShell>;
      auto problem = build_es_problem_tda<ClosedShell>(
          world, gs, num_roots, /*c_xc=*/1.0, gs.params().lo());
      // Initial state: either reload from a previous --save-roots dir
      // (skips the warmup entirely), or run the oversampled-warmup
      // power iteration (no-op fast path when both knobs are off).
      Solver::State state0;
      bool seeded = false;
      // (1) --load-best: cross-protocol restart via response_metadata.json.
      // Non-exact matches are reprojected by the first prepare().
      if (!load_best_dir.empty()) {
        auto loaded =
            try_load_es_bundle<TDA, ClosedShell>(world, load_best_dir);
        if (loaded) { state0 = std::move(loaded->state); seeded = true; }
      }
      // (2) --load-roots: exact bundle dir (existing 13d/Inc-9 path).
      if (!seeded && !load_roots_dir.empty()) {
        state0  = load_es_roots<TDA, ClosedShell>(world, load_roots_dir);
        seeded  = true;
      }
      // (3) Fresh: oversampled-warmup power iteration.
      if (!seeded) {
        state0 = run_oversampled_tda_warmup<ClosedShell>(
            world, gs, num_roots, n_roots_warmup, tda_warmup_iters,
            policy, /*c_xc=*/1.0, gs.params().lo(), print_level,
            guess_mode);
      }
      Solver solver(world, std::move(problem), main_policy, print_level);
      solver.set_batched(es_batch);
      solver.set_gamma_tensor(es_tensor);
      solver.set_time_phases(es_time);
      if (!log_path.empty()) solver.set_log_path(log_path);

      auto prepare = [&](double thresh, Solver &solv, Solver::State &st) {
        set_response_protocol(world, header.L, thresh, override_k);
        const int    new_k = FunctionDefaults<3>::get_k();
        const double new_t = FunctionDefaults<3>::get_thresh();
        if (world.rank() == 0)
          print("\n--- protocol step: thresh =", new_t,
                "  k =", new_k, " ---");
        auto coulop_new = poperatorT(
            CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
        gs.prepare(world, 0.001 * new_t, coulop_new, fock_json);
        solv.set_gs(build_response_ground_state_closed_shell(
            world, gs, /*c_xc=*/1.0, gs.params().lo()));
        for (auto &root : st.roots)
          for (auto &fn : root.x_alpha)
            fn = madness::project(fn, new_k, new_t);
        for (auto &rho : st.rho_alpha_prev)
          rho = madness::project(rho, new_k, new_t);
      };
      auto sf = solvers::iterate_protocol(
          solver, state0, protocol_thresholds, prepare, proto_policy);
      final_omega     = sf.omega;
      final_residuals = sf.last_bsh_residual;
      final_diverged        = sf.diverged;
      final_residual_target = solver.targets().bsh_residual;
      if (!save_roots_dir.empty()) {
        save_es_roots<TDA, ClosedShell>(world, sf, save_roots_dir,
                                        /*converged=*/!sf.diverged);
        if (world.rank() == 0)
          print("[SAVE] es_roots: wrote", num_roots,
                "root(s) +/", save_roots_dir, "/roots.json");
      }
    } else {
      using Solver = ESSolver<TDA, OpenShell>;
      auto problem = build_es_problem_tda<OpenShell>(
          world, gs, num_roots, gs.hf_exchange_coefficient(),
          gs.params().lo());
      Solver::State state0;
      if (!load_roots_dir.empty()) {
        state0 = load_es_roots<TDA, OpenShell>(world, load_roots_dir);
      } else {
        state0 = run_oversampled_tda_warmup<OpenShell>(
            world, gs, num_roots, n_roots_warmup, tda_warmup_iters,
            policy, gs.hf_exchange_coefficient(), gs.params().lo(),
            print_level, guess_mode);
      }
      Solver solver(world, std::move(problem), main_policy, print_level);
      if (!log_path.empty()) solver.set_log_path(log_path);

      auto prepare = [&](double thresh, Solver &solv, Solver::State &st) {
        set_response_protocol(world, header.L, thresh, override_k);
        const int    new_k = FunctionDefaults<3>::get_k();
        const double new_t = FunctionDefaults<3>::get_thresh();
        if (world.rank() == 0)
          print("\n--- protocol step: thresh =", new_t,
                "  k =", new_k, " ---");
        auto coulop_new = poperatorT(
            CoulombOperatorPtr(world, gs.params().lo(), 0.001 * new_t));
        gs.prepare(world, 0.001 * new_t, coulop_new, fock_json);
        solv.set_gs(build_response_ground_state_open_shell(
            world, gs, gs.hf_exchange_coefficient(), gs.params().lo()));
        for (auto &root : st.roots) {
          for (auto &fn : root.x_alpha)
            fn = madness::project(fn, new_k, new_t);
          for (auto &fn : root.x_beta)
            fn = madness::project(fn, new_k, new_t);
        }
        for (auto &rho : st.rho_alpha_prev)
          rho = madness::project(rho, new_k, new_t);
      };
      auto sf = solvers::iterate_protocol(
          solver, state0, protocol_thresholds, prepare, proto_policy);
      final_omega     = sf.omega;
      final_residuals = sf.last_bsh_residual;
      final_diverged        = sf.diverged;
      final_residual_target = solver.targets().bsh_residual;
      if (!save_roots_dir.empty()) {
        save_es_roots<TDA, OpenShell>(world, sf, save_roots_dir,
                                      /*converged=*/!sf.diverged);
        if (world.rank() == 0)
          print("[SAVE] es_roots: wrote", num_roots,
                "root(s) +/", save_roots_dir, "/roots.json");
      }
    }

    // -------------------------------------------------------------------
    // Compare against legacy reference.
    // -------------------------------------------------------------------
    auto state_final_omega    = final_omega;     // alias for the rest
    auto state_final_residual = final_residuals;
    if (world.rank() == 0) {
      print("\n=== Final omegas ===");
      for (long s = 0; s < num_roots; ++s) {
        print("  omega[", s, "] = ", final_omega(s),
              "   res = ", final_residuals[s]);
      }
    }

    // Convergence gate — checked BEFORE the per-root ω comparison. A run
    // is only allowed to PASS if it actually converged: not diverged, and
    // every root's BSH/Davidson residual below the solver's own target.
    // Without this, the ω-tolerance check passes a run that stalled at
    // maxiter on a residual floor — which is exactly how the
    // indefinite-metric bug stayed green (ω within tol, residual ~7e-3,
    // "stopped at maxiter, not converged"). The reference number is only
    // meaningful for a converged state.
    double max_residual = 0.0;
    for (double r : final_residuals) max_residual = std::max(max_residual, r);
    const bool solve_converged =
        !final_diverged && !final_residuals.empty() &&
        (max_residual <= final_residual_target);
    if (world.rank() == 0) {
      print("\n=== convergence gate ===");
      print("  diverged       =", final_diverged);
      print("  max_residual   =", max_residual,
            "   target =", final_residual_target);
      print("  CONVERGED      =", solve_converged,
            solve_converged ? "" : "  <-- run did NOT converge → FAIL");
    }

    // Two comparisons:
    //   (a) PASS/FAIL vs legacy molresponse TDA at the same protocol.
    //       A root passes only if the solve CONVERGED and |Δω| < tol
    //       (default 5 mhartree, --tol overrides).
    //   (b) Informational vs DALTON RPA at the given basis. Not part
    //       of the pass/fail; difference scales with TDA-vs-RPA gap.
    int passed = 0, failed = 0;
    if (ref) {
      if (!ref->legacy_tda.empty()) {
        const long n_compare = std::min<long>(
            num_roots, static_cast<long>(ref->legacy_tda.size()));
        if (world.rank() == 0)
          print("\n=== PASS/FAIL  vs legacy molresponse TDA (",
                ref->name, " RHF) ===");
        for (long i = 0; i < n_compare; ++i) {
          double err = std::abs(final_omega(i) - ref->legacy_tda[i]);
          bool ok = solve_converged && (err < tol);
          if (world.rank() == 0) {
            print("  [", (ok ? "PASS" : "FAIL"),
                  "]  omega[", i, "]  ref=", ref->legacy_tda[i],
                  "  actual=", final_omega(i),
                  "  err=", err, "  tol=", tol,
                  (solve_converged ? "" : "  (not converged)"));
          }
          (ok ? passed : failed)++;
        }
        if (world.rank() == 0)
          print("\n  ", passed, " passed, ", failed, " failed\n");
      }

      if (!ref->dalton_rpa.empty()) {
        const long n_compare = std::min<long>(
            num_roots, static_cast<long>(ref->dalton_rpa.size()));
        if (world.rank() == 0)
          print("=== INFO  vs DALTON RPA (", ref->name, ", basis ",
                ref->dalton_rpa_basis,
                ") — informational, not pass/fail ===");
        for (long i = 0; i < n_compare; ++i) {
          double err = std::abs(final_omega(i) - ref->dalton_rpa[i]);
          if (world.rank() == 0) {
            print("  omega[", i, "]  ref=", ref->dalton_rpa[i],
                  "  actual=", final_omega(i),
                  "  err=", err);
          }
        }
        if (world.rank() == 0) print("");
      }
    }
    rc = (failed == 0) ? 0 : 1;
    world.gop.fence();
    }  // <-- end inner scope; solver/state/target/gs destruct here
    world.gop.fence();
    finalize();
    return rc;
  } catch (const std::exception &e) {
    if (world.rank() == 0) print("ERROR:", e.what());
    finalize();
    return 2;
  }
}
