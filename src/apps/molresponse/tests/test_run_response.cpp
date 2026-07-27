// test_run_response — R0a: exercise the public run_response seam (doc 16, L1).
//
// The clean Input-builder reference: parse a minimal CLI, build a Tier-A plan +
// ExecutorSettings, call run_response once, and validate the returned Output.
// This proves the orchestrator reproduces the existing CalcManager numbers
// through a single self-contained entry point (the same one madqc will use, R3).
// The richer (entangled) test_calc_manager_run driver migrates onto this seam in
// R0b alongside the main.cpp / header consolidation.

#include "../orchestrator/response_workflow.hpp"

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <algorithm>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

std::vector<char> parse_axes(const std::string &s) {
  std::vector<char> ax;
  for (char c : s) {
    char l = std::tolower(c);
    if (l == 'x' || l == 'y' || l == 'z') ax.push_back(l);
  }
  if (ax.empty()) ax = {'x', 'y', 'z'};
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
      print(" molresponse_v3 run_response seam test  (R0a)");
      print("=========================================================\n");
    }
    {
      commandlineparser parser(argc, argv);
      if (!parser.key_exists("archive")) {
        if (world.rank() == 0)
          print("Usage: test_run_response --archive=<path> [--omega=0.0,0.057] "
                "[--axes=xyz] [--protocol=1e-4,1e-6] [--maxiter=N] [--dconv=X] "
                "[--calc-dir=DIR] [--print-level=0..3] [--beta [--beta-static]]");
        finalize();
        return 1;
      }

      // ---- Parse ----------------------------------------------------------
      ResponseWorkflowInput in;
      in.archive_file = parser.value_raw("archive");
      in.protocols = parser.key_exists("protocol")
                         ? parse_protocol_csv(parser.value("protocol"))
                         : std::vector<double>{1e-4, 1e-6};
      const std::vector<double> freqs =
          parser.key_exists("omega") ? parse_protocol_csv(parser.value("omega"))
                                     : std::vector<double>{0.0};
      const std::vector<char> axes =
          parse_axes(parser.key_exists("axes") ? parser.value("axes") : "xyz");

      // ---- Settings (World-free) -----------------------------------------
      ConvergencePolicy policy;
      policy.dconv_user = parser.key_exists("dconv")
                              ? std::stod(parser.value("dconv"))
                              : 1e-4;
      if (parser.key_exists("conv-factor")) {
        const double f = std::stod(parser.value("conv-factor"));
        policy.bsh_residual_factor = f;
        policy.density_residual_factor = f;
      }
      in.settings.policy = policy;
      in.settings.max_iters =
          parser.key_exists("maxiter") ? std::stoi(parser.value("maxiter")) : 25;
      in.settings.calc_dir = parser.key_exists("calc-dir")
                                 ? parser.value_raw("calc-dir")
                                 : std::string("run_response_test");
      const int pl = parser.key_exists("print-level")
                         ? std::stoi(parser.value("print-level"))
                         : 1;
      in.settings.print_level =
          static_cast<PrintLevel>(std::max(0, std::min(3, pl)));
      if (parser.key_exists("accept-at-maxiter"))
        in.settings.accept_at_maxiter = true;

      // ---- Tier-A plan (alpha by default; --beta for SHG/static) ----------
      const bool do_beta = parser.key_exists("beta");
      ResponsePropertyRequest req;
      req.axes = axes;
      req.protocol_thresholds = in.protocols;
      req.frequencies = freqs;
      if (do_beta) {
        req.kind = ResponsePropertyKind::Hyperpolarizability;
        req.beta_process = parser.key_exists("beta-static") ? BetaProcess::Static
                                                            : BetaProcess::SHG;
      } else {
        req.kind = ResponsePropertyKind::Polarizability;
      }
      in.plan = plan_one(req);

      if (world.rank() == 0) {
        print("RUN CONFIG (run_response):");
        print("  archive    =", in.archive_file);
        print("  omega      =", freqs);
        print("  protocol   =", in.protocols);
        print("  calc_dir   =", in.settings.calc_dir);
        print("  mode       =", (do_beta ? "hyperpolarizability" : "polarizability"));
        print("  FD requests=", (int)in.plan.fd.size(),
              "  VBC pairs=", (int)in.plan.vbc.size());
      }

      // ---- The seam -------------------------------------------------------
      ResponseWorkflowOutput out = run_response(world, in);

      // ---- Validate (rank 0) ---------------------------------------------
      if (world.rank() == 0) {
        // R1c: confirm the scheduler-trace diagnostics block populated.
        if (out.diagnostics.contains("schedule"))
          print("diagnostics: passes=", out.diagnostics.value("passes", -1),
                "  stop_reason=",
                out.diagnostics.value("stop_reason", std::string("?")),
                "  waves=", (int)out.diagnostics["schedule"].size());
        const std::string top_key = protocol_key_at(in.protocols.back());
        const auto &j = out.metadata;
        if (do_beta) {
          int want = 3 * static_cast<int>(in.plan.vbc.size());
          int got = (j.contains("properties") &&
                     j["properties"].contains("beta") &&
                     j["properties"]["beta"].contains(top_key))
                        ? static_cast<int>(j["properties"]["beta"][top_key].size())
                        : 0;
          bool ok = want > 0 && got == want;
          print("  beta elements recorded=", got, "/", want);
          print("\n", ok ? "PASSED" : "FAILED", " (beta ", got, "/", want, ")");
          rc = ok ? 0 : 1;
        } else {
          int expected = 0, converged = 0;
          for (const auto &r : in.plan.fd) {
            ++expected;
            const std::string pert = r.pert.description();
            const std::string fk = ResponseMetadata::freq_key(r.freq);
            const bool ok =
                j.contains("fd_states") && j["fd_states"].contains(pert) &&
                j["fd_states"][pert].contains(top_key) &&
                j["fd_states"][pert][top_key].contains(fk) &&
                j["fd_states"][pert][top_key][fk].value("converged", false);
            if (ok) ++converged;
            print("  ", ok ? "[PASS]" : "[FAIL]", pert, "@", r.freq,
                  " key=", top_key);
          }
          const bool alpha_ok =
              j.contains("properties") && j["properties"].contains("alpha") &&
              j["properties"]["alpha"].contains(top_key);
          print("  alpha recorded =", alpha_ok);
          bool ok = (converged == expected) && alpha_ok;
          print("\n", ok ? "PASSED" : "FAILED", " (", converged, "/", expected,
                " FD converged, alpha ", (alpha_ok ? "present" : "missing"), ")");
          rc = ok ? 0 : 1;
        }
      }
      world.gop.broadcast(rc, 0);
    }
  } catch (const std::exception &e) {
    if (world.rank() == 0) print("EXCEPTION:", e.what());
    rc = 2;
  }
  finalize();
  return rc;
}
