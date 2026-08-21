// ===========================================================================
// test_calc_manager_run.cpp — end-to-end drive of the calc manager (doc 15,
// 15a). Sets up a ground state from a moldft archive (same recipe as the FD
// skeleton test), plans a polarizability request, then lets CalcManager +
// FdResponseExecutor schedule and solve every FD protocol step. Validates by
// reading back response_metadata.json and checking the expected FD states
// converged.
//
// This is an ALLOCATION test (it runs real MRA solves). Usage:
//   test_calc_manager_run --archive=<moldft restartdata>
//       [--omega=0.0,0.057] [--axes=xyz] [--protocol=1e-4,1e-6] [--es-roots=N]
//       [--maxiter=N] [--dconv=X] [--calc-dir=DIR] [--print-level=0..3]
//       [--conv-factor=F | --bsh-factor=F --density-factor=F]
//   --conv-factor loosens/tightens both convergence gates (target = F*dconv,
//   default F=5); --bsh-factor / --density-factor set each gate independently.
//   --accept-at-maxiter accepts a non-diverged FD that exhausts --maxiter
//   without hitting the target, so a stiff channel climbs the protocol ladder
//   and VBC proceeds (recorded with an `accepted` marker + its real residual).
// ===========================================================================

#include "../GroundState.hpp"
#include "../ResponsePropertyPlanner.hpp"
#include "../ResponseProtocol.hpp"
#include "../calc/calc_executor.hpp"
#include "../solvers/convergence_policy.hpp"
#include "../solvers/dalton_import.hpp"
#include "../solvers/response_metadata.hpp"

#include <nlohmann/json.hpp>
#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <madness/world/worldprofile.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

std::vector<double> parse_csv_doubles(const std::string &s) {
  return parse_protocol_csv(s); // reuses the protocol CSV parser
}

std::vector<char> parse_axes(const std::string &s) {
  std::vector<char> ax;
  for (char c : s) {
    char l = std::tolower(c);
    if (l == 'x' || l == 'y' || l == 'z')
      ax.push_back(l);
  }
  if (ax.empty())
    ax = {'x', 'y', 'z'};
  return ax;
}

} // namespace

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    if (world.rank() == 0) {
      print("\n=========================================================");
      print(" molresponse_v3 CalcManager run test  (15a, FD/alpha)");
      print("=========================================================\n");
    }
    {
      commandlineparser parser(argc, argv);
      if (!parser.key_exists("archive")) {
        if (world.rank() == 0)
          print("Usage: test_calc_manager_run --archive=<path> "
                "[--omega=0.0,0.057] [--axes=xyz] [--protocol=1e-4,1e-6] "
                "[--maxiter=N] [--dconv=X] [--calc-dir=DIR] "
                "[--print-level=0..3] "
                "[--dalton-dir=DIR [--dalton-molden=F] [--dalton-rspvec=F] "
                "[--dalton-out=F]]");
        finalize();
        return 1;
      }

      const std::string archive_path = parser.value_raw("archive");
      const std::vector<double> freqs =
          parser.key_exists("omega") ? parse_csv_doubles(parser.value("omega"))
                                     : std::vector<double>{0.0};
      const std::vector<char> axes =
          parse_axes(parser.key_exists("axes") ? parser.value("axes") : "xyz");
      const int max_iters = parser.key_exists("maxiter")
                                ? std::stoi(parser.value("maxiter"))
                                : 25;
      const double dconv =
          parser.key_exists("dconv") ? std::stod(parser.value("dconv")) : 1e-4;
      const int pl_int = parser.key_exists("print-level")
                             ? std::stoi(parser.value("print-level"))
                             : 1;
      const PrintLevel print_level =
          static_cast<PrintLevel>(std::max(0, std::min(3, pl_int)));

      auto header = GroundState::read_archive_header(world, archive_path);

      std::vector<double> protocol;
      if (parser.key_exists("protocol"))
        protocol = parse_protocol_csv(parser.value("protocol"));
      else
        protocol.push_back(default_thresh_for_k(header.k));
      // Analyze-only reports at the TOP protocol and never solves, so set that
      // protocol up front: the single gs.prepare() below loads + projects the
      // orbitals at the target (k,thresh) directly. This avoids a second
      // back-to-back gs.prepare() (which reruns exchange/MacroTaskQ and, with
      // no solve/fence between the two preparations, left the ground-state
      // orbitals bound to a torn-down subworld -> dead-mutex crash at
      // teardown).
      const bool analyze_only = parser.key_exists("es-analyze-only");
      set_response_protocol(world, header.L,
                            analyze_only ? protocol.back() : protocol.front());

      // Molecule (for natom -> nuclear expansion; harmless for pure alpha).
      Molecule molecule;
      auto archive_dir = std::filesystem::path(archive_path).parent_path();
      for (const auto &name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
        auto candidate = archive_dir / name;
        if (std::filesystem::exists(candidate)) {
          std::ifstream ifs(candidate);
          nlohmann::json j;
          ifs >> j;
          nlohmann::json mol_json;
          if (j.contains("tasks") && j["tasks"].is_array() &&
              !j["tasks"].empty())
            mol_json = j["tasks"][0]["molecule"];
          else if (j.contains("molecule"))
            mol_json = j["molecule"];
          if (!mol_json.is_null())
            molecule.from_json(mol_json);
          break;
        }
      }
      auto gs = GroundState::from_archive(world, archive_path, molecule);
      if (world.rank() == 0)
        gs.print_info();

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

      ConvergencePolicy policy;
      policy.dconv_user = dconv;
      // --fd-tensor: gate the FD theta build onto the tensor-exchange layer
      // (exch::; ClosedShell Static/Full only). Default off = per-op reference.
      if (parser.key_exists("fd-tensor")) policy.exchange_tensor = true;
      // --fd-tensor-tile=N: block the tensor-exchange Tx/Ty/g0 builds over the
      // φ-row index to bound the n² memory peak (doc 33 Inc-3b). 0 = no tiling
      // (default). Bit-identical to tile 0. Only meaningful with --fd-tensor.
      if (parser.key_exists("fd-tensor-tile"))
        policy.exchange_tile = std::stoi(parser.value("fd-tensor-tile"));
      // --explosion-guard=VAL caps the diverging-residual bail-out;
      // --no-explosion- guard disables it (test whether the guard just fires
      // early for stiff nuclear-displacement FD states).
      if (parser.key_exists("explosion-guard"))
        policy.explosion_guard = std::stod(parser.value("explosion-guard"));
      if (parser.key_exists("no-explosion-guard"))
        policy.explosion_guard = 1.0e30;
      // --step-restrict=orbital|state  (default per-state; see ConvergencePolicy).
      if (parser.key_exists("step-restrict") &&
          parser.value("step-restrict") == "state")
        policy.step_restrict_mode =
            ConvergencePolicy::StepRestrictMode::PerState;
      if (parser.key_exists("omega-factor"))
        policy.omega_residual_factor = std::stod(parser.value("omega-factor"));
      if (parser.key_exists("lock-min-pass"))
        policy.lock_min_pass = std::stoi(parser.value("lock-min-pass"));
      // Convergence-gate factors (default 5.0 each, SCF convention). The
      // effective target is factor*max(thresh,dconv). Loosen these — rather
      // than inflating --dconv — to pass a coarse protocol whose residual
      // floors above 5*dconv (e.g. the stiff nuclear-displacement Raman FD
      // channel). --conv-factor sets both; --bsh-factor / --density-factor
      // override each gate independently.
      if (parser.key_exists("conv-factor")) {
        const double f = std::stod(parser.value("conv-factor"));
        policy.bsh_residual_factor = f;
        policy.density_residual_factor = f;
      }
      if (parser.key_exists("bsh-factor"))
        policy.bsh_residual_factor = std::stod(parser.value("bsh-factor"));
      if (parser.key_exists("density-factor"))
        policy.density_residual_factor =
            std::stod(parser.value("density-factor"));

      const std::string calc_dir = parser.key_exists("calc-dir")
                                       ? parser.value_raw("calc-dir")
                                       : std::string("calc_manager_run");

      // ---- --es-analyze-only: load a converged ES bundle and print the
      // transition-property report WITHOUT solving. Restart precedence picks
      // the bundle at (or coarser than) the finest requested protocol; roots
      // are reprojected to the active k/thresh (same path the solver takes).
      // Add --es-full to analyze a Full (X,Y) bundle instead of TDA. This is
      // also an isolated harness for the analysis collectives.
      if (parser.key_exists("es-analyze-only")) {
        if (world.rank() == 0)
          print("1:");
        const double th = protocol.back();
        set_response_protocol(world, header.L, th);
        const double t0 = FunctionDefaults<3>::get_thresh();
        auto cop =
            poperatorT(CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
        gs.prepare(world, 0.001 * t0, cop, fock_json);
        if (world.rank() == 0)
          print("1:");

        auto reproject = [&](auto &st) {
          const int k = FunctionDefaults<3>::get_k();
          const double tt = FunctionDefaults<3>::get_thresh();
          for (auto &root : st.roots)
            for (auto *blk : root.blocks())
              for (auto &fn : *blk)
                fn = madness::project(fn, k, tt);
        };

        const bool full = parser.key_exists("es-full");
        bool found = false;
        if (full) {
          if (auto loaded =
                  try_load_es_bundle<Full, ClosedShell>(world, calc_dir)) {
            reproject(loaded->state);
            report_es_analysis<Full, ClosedShell>(world, gs, loaded->state,
                                                  print_level);
            found = true;
          }
        } else {
          if (auto loaded =
                  try_load_es_bundle<TDA, ClosedShell>(world, calc_dir)) {
            reproject(loaded->state);
            report_es_analysis<TDA, ClosedShell>(world, gs, loaded->state,
                                                 print_level);
            found = true;
          }
        }
        if (!found && world.rank() == 0)
          print("[ANALYZE] no", (full ? "Full" : "TDA"),
                "ES bundle found under", calc_dir);
        if (world.rank() == 0)
          print("finalize:");
        finalize();
        return 0;
      }

        if (world.rank() == 0)
          print("finalize:");
      // ---- Plan: polarizability, or (--es-roots=N) resonant gradient -----
      const int es_roots = parser.key_exists("es-roots")
                               ? std::stoi(parser.value("es-roots"))
                               : 0;
      const bool do_beta = parser.key_exists("beta");
      ResponsePropertyRequest req;
      req.axes = axes;
      req.protocol_thresholds = protocol;
      if (do_beta) {
        // Hyperpolarizability beta via 2n+1 / VBC contraction (SHG by default).
        req.kind = ResponsePropertyKind::Hyperpolarizability;
        req.beta_process = parser.key_exists("beta-static")
                               ? BetaProcess::Static
                               : BetaProcess::SHG;
        req.frequencies = freqs;
      } else if (parser.key_exists("raman")) {
        // Vibrational Raman: β(A=dipole; B=dipole@ω, C=nuclear@0). SINGLE
        // component for now (--nuc-atom/--nuc-axis; default atom 0, axis z).
        req.kind = ResponsePropertyKind::PolarizabilityGradient;
        req.gradient_mode = GradientMode::Nuclear;
        req.frequencies = freqs;
        req.raman_nuc_atom = parser.key_exists("nuc-atom")
                                 ? std::stoi(parser.value("nuc-atom"))
                                 : 0;
        req.raman_nuc_axis = parser.key_exists("nuc-axis")
                                 ? std::stoi(parser.value("nuc-axis"))
                                 : 2;
      } else if (es_roots > 0) {
        // ES(TDA) bundle + symbolic derived dipole FD; the calc manager solves
        // the bundle, then expands to FD at the converged excitation energies.
        req.kind = ResponsePropertyKind::PolarizabilityGradient;
        req.gradient_mode = GradientMode::Resonant;
        req.n_roots = es_roots;
      } else {
        req.kind = ResponsePropertyKind::Polarizability;
        req.frequencies = freqs;
      }
      ResponsePlan plan = plan_one(req);
      // --es-full: drive the Full (paired X,Y) closed-shell ES path instead of
      // TDA.
      if (es_roots > 0 && parser.key_exists("es-full"))
        for (auto &e : plan.es)
          e.tda = false;

      if (world.rank() == 0) {
        print("\nRUN CONFIG:");
        print("  archive    =", archive_path);
        print("  shell      =",
              gs.is_spin_restricted() ? "ClosedShell" : "OpenShell");
        print("  omega      =", freqs);
        print("  protocol   =", protocol);
        print("  calc_dir   =", calc_dir);
        if (parser.key_exists("dalton-dir"))
          print("  dalton_dir =", parser.value_raw("dalton-dir"),
                " (FD seed import)");
        print("  mode       =", (do_beta ? "hyperpolarizability (beta/VBC)"
                                 : es_roots > 0 ? "resonant-gradient (ES)"
                                                : "polarizability"));
        if (es_roots > 0)
          print("  es_roots   =", es_roots);
        print("  FD requests=", (int)plan.fd.size(),
              "  ES requests=", (int)plan.es.size());
        print("  maxiter    =", max_iters,
              "  conv: bsh_factor=", policy.bsh_residual_factor,
              " density_factor=", policy.density_residual_factor,
              " dconv=", policy.dconv_user);
        print("  accept_at_maxiter =",
              (parser.key_exists("accept-at-maxiter") ? "ON" : "off"));
        print("  fd_tensor  =", (policy.exchange_tensor ? "ON" : "off"),
              "  tile =", policy.exchange_tile,
              "  es_tensor =", (parser.key_exists("es-tensor") ? "ON" : "off"));
      }

      // ---- Analyze-only: load a converged ES bundle + print the report ----
      // (no solve). Prepares the GS at the top protocol, then loads + reports
      // the bundle from calc_dir. Isolates the analysis path from the solver.
      if (analyze_only) {
        // gs was already prepared once at protocol.back() above (the protocol
        // was set to the top rung up front for analyze-only) — do NOT
        // re-prepare here; a second back-to-back prepare reruns
        // exchange/MacroTaskQ and corrupts ground-state-orbital world ownership
        // at teardown.
        const bool full = parser.key_exists("es-full");
        // --es-vectors enables the heavy per-orbital AO-population pass
        // (analyze_response_orbitals); OFF by default while the teardown crash
        // is bisected — it's the current suspect.
        const bool do_vec = parser.key_exists("es-vectors");
        // --es-load-only: bisection — load+reproject, skip the report entirely.
        const bool load_only = parser.key_exists("es-load-only");
        rc = full ? analyze_es_bundle_from_disk<Full, ClosedShell>(
                        world, gs, calc_dir, print_level, do_vec, load_only)
                  : analyze_es_bundle_from_disk<TDA, ClosedShell>(
                        world, gs, calc_dir, print_level, do_vec, load_only);
        world.gop.broadcast(rc, 0);
        // NB: do NOT finalize()/return here — that would tear down the runtime
        // while gs (and its MADNESS Functions) are still in scope; their
        // destructors would then hit a dead mutex. Fall through to the single
        // finalize() after the try block, so locals are destroyed first.
      } else {

        // ---- dalton.dir import (showcase W3): seed the planned FD states
        // from a DALTON linear-response directory BEFORE the manager runs, so
        // reconcile sees the seed bundles as coarser-or-equal restart sources.
        // Fingerprint (geometry vs the active Molecule) is a HARD error; so is
        // any requested frequency absent from the RSPVEC (exact 1-to-1). The
        // GS is already prepared at protocol.front() above — the precondition
        // run_dalton_import documents.
        if (parser.key_exists("dalton-dir")) {
          auto get_opt = [&](const char *k) {
            return parser.key_exists(k) ? parser.value_raw(k) : std::string{};
          };
          run_dalton_import(world, gs, molecule, plan, calc_dir,
                            parser.value_raw("dalton-dir"),
                            get_opt("dalton-molden"), get_opt("dalton-rspvec"),
                            get_opt("dalton-out"));
        }

        // ---- Drive the calc manager ----------------------------------------
        CalcManager::Policy mgr_policy;
        mgr_policy.max_iters_per_step = max_iters;
        CalcManager mgr(plan, calc_dir, mgr_policy);
        mgr.build(molecule.natom());

        // ExecutorContext now = {world, gs, L, fock_json} + inherited
        // ExecutorSettings (R0a split). The ctx.<knob> setters below still work
        // (inherited members). FUTURE (R0b): build the Input + call run_response.
        ExecutorContext ctx(world, gs, header.L, fock_json,
                            ExecutorSettings{policy, print_level, calc_dir,
                                             max_iters});
        if (parser.key_exists("es-seed"))
          ctx.seed_derived_from_es_root = true;
        if (parser.key_exists("no-es-seed"))
          ctx.seed_derived_from_es_root = false;
        // --es-tensor: Inc-3c tensor-layer ES γ (cached g0 + batched Coulomb;
        // TDA/ClosedShell only). Separate gate from --fd-tensor (FD θ path).
        if (parser.key_exists("es-tensor"))
          ctx.es_gamma_tensor = true;
        // --accept-at-maxiter: accept a non-diverged FD solve that hits maxiter
        // without meeting the strict target (records converged + an `accepted`
        // marker) so a stiff channel climbs the protocol ladder and VBC
        // proceeds. Pair with --maxiter (iters/protocol budget) and --protocol
        // (the finest entry is the de-facto final rung; its bsh_residual is the
        // verdict).
        if (parser.key_exists("accept-at-maxiter"))
          ctx.accept_at_maxiter = true;
        // ES/KAIN experiment knobs (workstreams A + C; sweepable without
        // rebuild).
        if (parser.key_exists("kain-maxsub"))
          ctx.es_kain_maxsub = std::stoi(parser.value("kain-maxsub"));
        if (parser.key_exists("kain-delay"))
          ctx.es_main_kain_delay = std::stoi(parser.value("kain-delay"));
        if (parser.key_exists("es-warmup-iters"))
          ctx.es_tda_warmup_iters = std::stoi(parser.value("es-warmup-iters"));
        if (parser.key_exists("es-oversample"))
          ctx.es_warmup_oversample = std::stod(parser.value("es-oversample"));
        if (parser.key_exists("maxrotn"))
          ctx.es_maxrotn = std::stod(parser.value("maxrotn"));
        if (parser.key_exists("es-guess"))
          ctx.es_guess = parse_es_guess_mode(parser.value("es-guess"));
        // --tpa-residue [--tpa-prefactor=X]: corrected single-residue 2PA
        // contraction (kernels/tpa.hpp tpa_moment_residue) instead of the
        // legacy beta-reuse candidate.
        if (parser.key_exists("tpa-residue")) ctx.tpa_residue = true;
        if (parser.key_exists("tpa-prefactor"))
          ctx.tpa_prefactor = std::stod(parser.value("tpa-prefactor"));
        if (parser.key_exists("tpa-decompose")) ctx.tpa_decompose = true;
        if (parser.key_exists("tpa-diag-only")) ctx.tpa_diag_only = true;
        // --tpa-roots=0,2 : per-root filter (parallel verification; read-only)
        if (parser.key_exists("tpa-roots")) {
          std::stringstream rss(parser.value("tpa-roots"));
          std::string rtok;
          while (std::getline(rss, rtok, ','))
            if (!rtok.empty()) ctx.tpa_roots.push_back(std::stoi(rtok));
          // concurrent per-root processes share one calc dir: no persistence
          ResponseMetadata::read_only() = true;
        }
          // accepts: random | solid[_harmonics] | virtual[_ao]
        if (parser.key_exists("lock-converged"))
          ctx.es_lock_converged = true;
        if (parser.key_exists("no-lock-converged"))
          ctx.es_lock_converged = false;
        if (parser.key_exists("es-warmup-cache"))
          ctx.es_warmup_cache = true;
        if (parser.key_exists("es-warmup-cache-dir")) {
          ctx.es_warmup_cache_dir = parser.value_raw("es-warmup-cache-dir");
          ctx.es_warmup_cache = true; // a cache dir implies caching
        }
        if (world.rank() == 0 && es_roots > 0) {
          print("ES/KAIN knobs: guess=", to_string(ctx.es_guess),
                " kain_maxsub=", ctx.es_kain_maxsub,
                " kain_delay=", ctx.es_main_kain_delay,
                " warmup_iters=", ctx.es_tda_warmup_iters,
                " oversample=", ctx.es_warmup_oversample,
                " maxrotn=", ctx.es_maxrotn, " step_restrict=",
                (policy.step_restrict_mode ==
                         ConvergencePolicy::StepRestrictMode::PerState
                     ? "state"
                     : "orbital"),
                " lock_converged=", (ctx.es_lock_converged ? "on" : "off"),
                " warmup_cache=", (ctx.es_warmup_cache ? "on" : "off"));
        }
        FdResponseExecutor exec(ctx);
        mgr.run(world, exec);

        // beta: Tier-A property assembly (contraction) after the manager run.
        // beta OR raman both fill plan.vbc + share assemble_beta (it tags raman
        // vs beta per pair). alpha only for the plain-FD path.
        if (!plan.vbc.empty())
          assemble_beta(ctx, plan, protocol.back());
        else if (es_roots == 0)
          assemble_alpha(ctx, plan, protocol.back());

        // --tpa: two-photon-absorption assembly after the ES bundle + derived
        // dipole FD at omega_f/2 converge (validates kernels/tpa.hpp vs Dalton).
        if (es_roots > 0 && parser.key_exists("tpa"))
          assemble_tpa(ctx, plan, protocol.back());

        // ---- Validate at the top protocol ----------------------------------
        const std::string top_key = protocol_key_at(protocol.back());
        if (world.rank() == 0) {
          auto meta = ResponseMetadata::load_or_create(
              calc_dir + "/response_metadata.json");
          const auto &j = meta.json();
          if (do_beta) {
            int want =
                3 * static_cast<int>(plan.vbc.size()); // A axes x VBC pairs
            int got = 0;
            if (j.contains("properties") && j["properties"].contains("beta") &&
                j["properties"]["beta"].contains(top_key))
              got = static_cast<int>(j["properties"]["beta"][top_key].size());
            print("  VBC pairs=", (int)plan.vbc.size(),
                  "  beta elements recorded=", got, "/", want);
            bool ok = (want > 0) && (got == want);
            print("\n", ok ? "PASSED" : "FAILED", " (beta: ", got, "/", want,
                  ")");
            rc = ok ? 0 : 1;
          } else if (es_roots > 0) {
            // ES bundle converged with n_roots, plus the promoted derived FDs.
            bool es_ok = j["excited_states"].contains(top_key) &&
                         j["excited_states"][top_key].value("converged", false);
            int n_es = (es_ok && j["excited_states"][top_key].contains("roots"))
                           ? (int)j["excited_states"][top_key]["roots"].size()
                           : 0;
            // Expected derived frequencies = the DISTINCT freq_key(0.5 * ωₙ)
            // over the converged roots — degenerate roots collapse to one
            // frequency, so the count is unique-ωₙ/2 × axes, not n_roots ×
            // axes. (0.5 mirrors the resonant DerivedFDRequest::es_freq_factor
            // default.)
            std::set<std::string> want_fkeys;
            if (es_ok && j["excited_states"][top_key].contains("roots"))
              for (const auto &r : j["excited_states"][top_key]["roots"]) {
                const double w =
                    r.value("omega", std::numeric_limits<double>::quiet_NaN());
                if (w == w)
                  want_fkeys.insert(ResponseMetadata::freq_key(0.5 * w));
              }
            const int want_derived = static_cast<int>(want_fkeys.size()) *
                                     static_cast<int>(axes.size());
            // Count converged derived FDs at the top protocol whose freq
            // matches an expected ωₙ/2 key (ignores any pre-fix ωₙ orphans on a
            // reused dir).
            int derived = 0;
            if (j.contains("fd_states"))
              for (auto &kv : j["fd_states"].items())
                if (kv.value().contains(top_key))
                  for (auto &fe : kv.value()[top_key].items())
                    if (want_fkeys.count(fe.key()) &&
                        fe.value().value("converged", false))
                      ++derived;
            print("  ES bundle  converged=", es_ok, "  roots=", n_es, "/",
                  es_roots, "  distinct ωₙ/2 freqs=", (int)want_fkeys.size());
            print("  derived FD converged=", derived, "/", want_derived,
                  " key=", top_key);
            bool ok = es_ok && n_es == es_roots && derived == want_derived;
            print("\n", ok ? "PASSED" : "FAILED", " (ES + ", derived, "/",
                  want_derived, " derived FD)");
            rc = ok ? 0 : 1;
          } else {
            int expected = 0, converged = 0;
            for (const auto &r : plan.fd) {
              ++expected;
              const std::string pert = r.pert.description();
              const std::string fk = ResponseMetadata::freq_key(r.freq);
              bool ok =
                  j["fd_states"].contains(pert) &&
                  j["fd_states"][pert].contains(top_key) &&
                  j["fd_states"][pert][top_key].contains(fk) &&
                  j["fd_states"][pert][top_key][fk].value("converged", false);
              if (ok)
                ++converged;
              print("  ", ok ? "[PASS]" : "[FAIL]", pert, "@", r.freq,
                    " key=", top_key);
            }
            print("\n", (converged == expected) ? "PASSED" : "FAILED", " (",
                  converged, "/", expected, " FD states converged)");
            rc = (converged == expected) ? 0 : 1;
          }

          // DAG identity invariant (mode-independent): every FD node's id must
          // encode its OWN frequency. Guards the ES-expansion id/freq mismatch
          // for es_freq_factor != 1.0 — a derived node lives at ωₙ/2 and must
          // not carry an id built from the unscaled root energy ωₙ.
          int id_bad = 0;
          for (const auto &n : mgr.dag())
            if (n.kind == CalcKind::FD && n.id != fd_node_id(n.pert, n.freq)) {
              ++id_bad;
              print("  [FAIL] FD node id/freq mismatch: id=", n.id,
                    " freq=", n.freq);
            }
          if (id_bad) {
            print("FAILED (", id_bad, " FD node id/freq mismatch(es))");
            rc = 1;
          }
        }
        world.gop.broadcast(rc, 0);
      } // end else (calc-manager drive path)
    }
  } catch (const std::exception &e) {
    if (world.rank() == 0)
      print("EXCEPTION:", e.what());
    rc = 2;
  }
  finalize();
  return rc;
}
