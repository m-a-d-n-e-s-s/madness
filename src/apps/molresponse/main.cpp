// molresponse_v3 — frequency-dependent response solver (app entry).
//
// Thin driver over the single orchestrator seam (R0b): parse a CLI, build a
// ResponseWorkflowInput, call run_response, print what landed. All solver,
// scheduling, persistence, timing, and diagnostics logic lives behind
// orchestrator/response_workflow.hpp -> CalcManager -> the kernels/solvers.
//
// The legacy direct-FD path (top-level FDSolver/ESSolver headers + a hand-rolled
// protocol×direction loop) was retired in R0b; this binary now exercises exactly
// the same path as the tests and (later, R3) madqc.

#include "orchestrator/response_workflow.hpp"

#include <nlohmann/json.hpp>
#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <madness/world/worldprofile.h>

#include <algorithm>
#include <cstdlib>
#include <sstream>
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

void print_usage() {
  print("molresponse — frequency-dependent response solver\n");
  print("Usage: molresponse --archive=<path> [options]");
  print("  --archive=<path>     ground-state checkpoint (moldft restartdata base)");
  print("  --omega=0.0,0.057    CSV of frequencies (default 0.0)");
  print("  --axes=xyz           Cartesian axes (default xyz)");
  print("  --protocol=1e-4,1e-6 truncation-threshold ladder (default 1e-4,1e-6)");
  print("  --maxiter=N          iteration budget per protocol step (default 25)");
  print("  --fd-subworlds=P     fan FD/VBC states across P subworlds PER NODE "
        "(0=single-World default; 1=node-aligned; >1=sub-node/NUMA packing)");
  print("  --hdf5               opt-in HDF5-backed restart I/O "
        "(requires -DMADNESS_ENABLE_HDF5 build; readers auto-detect .h5)");
  print("  --dconv=X            convergence target (default 1e-4)");
  print("  --calc-dir=DIR       output dir (response_metadata.json + archives)");
  print("  --print-level=0..3   verbosity (default 1)");
  print("  --beta [--beta-static]  hyperpolarizability β (SHG, or static) via VBC");
}

} // namespace

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    commandlineparser parser(argc, argv);

    if (parser.key_exists("help") || !parser.key_exists("archive")) {
      if (world.rank() == 0) print_usage();
      finalize();
      return parser.key_exists("help") ? 0 : 1;
    }

    if (world.rank() == 0) {
      print("\n====================================");
      print("  molresponse_v3 — response solver");
      print("====================================");
      print(info::print_revision_information());
    }

    {
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

      ConvergencePolicy policy;
      policy.dconv_user =
          parser.key_exists("dconv") ? std::stod(parser.value("dconv")) : 1e-4;
      in.settings.policy = policy;
      in.settings.max_iters =
          parser.key_exists("maxiter") ? std::stoi(parser.value("maxiter")) : 25;
      // F2 (doc 32 §5): fan independent FD states across node-aligned subworlds.
      // 0 = single-World reference path. >0 = requested # of subworlds.
      in.settings.fd_subworlds =
          parser.key_exists("fd-subworlds")
              ? std::stoi(parser.value("fd-subworlds")) : 0;
      if (parser.key_exists("hdf5")) {
#ifdef MADNESS_HAS_HDF5
        set_hdf5_io_enabled(true);
#else
        if (world.rank() == 0)
          print("ERROR: --hdf5 requested but this build has no HDF5 support "
                "(configure with -DMADNESS_ENABLE_HDF5=ON)");
        finalize();
        return 1;
#endif
      }
      in.settings.calc_dir = parser.key_exists("calc-dir")
                                 ? parser.value_raw("calc-dir")
                                 : std::string("molresponse_v3_calc");
      const int pl = parser.key_exists("print-level")
                         ? std::stoi(parser.value("print-level"))
                         : 1;
      in.settings.print_level =
          static_cast<PrintLevel>(std::max(0, std::min(3, pl)));

      ResponsePropertyRequest req;
      req.axes = axes;
      req.protocol_thresholds = in.protocols;
      req.frequencies = freqs;
      if (parser.key_exists("beta")) {
        req.kind = ResponsePropertyKind::Hyperpolarizability;
        req.beta_process = parser.key_exists("beta-static") ? BetaProcess::Static
                                                            : BetaProcess::SHG;
      } else {
        req.kind = ResponsePropertyKind::Polarizability;
      }
      in.plan = plan_one(req);

      if (world.rank() == 0) {
        print("\nPARAMETERS:");
        print("  archive    =", in.archive_file);
        print("  omega      =", freqs);
        print("  protocol   =", in.protocols);
        print("  calc_dir   =", in.settings.calc_dir);
        print("  mode       =",
              (parser.key_exists("beta") ? "hyperpolarizability (beta)"
                                         : "polarizability (alpha)"));
        print("");
      }

      ResponseWorkflowOutput out = run_response(world, in);

      if (world.rank() == 0) {
        const std::string top = protocol_key_at(in.protocols.back());
        const auto &p = out.properties;
        for (const char *name : {"alpha", "beta", "raman"})
          if (p.contains(name) && p[name].contains(top))
            print("[RESULT]", name, "recorded @", top);
        if (out.timing.contains("total"))
          print("[RESULT] total_wall_s =",
                out.timing["total"].value("wall_s", 0.0));
      }
    }
  } catch (const std::exception &e) {
    if (world.rank() == 0) print("EXCEPTION:", e.what());
    rc = 2;
  }
  finalize();
  return rc;
}
