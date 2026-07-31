// seed_es_from_hdf5 — build a molresponse TDA excited-state SEED bundle from
// per-occupied-orbital response functions x_i(r) stored as pymra structured
// .mad.h5 files (schema=1, the /meta+/keys+/coeffs layout written by
// pymra.io.write_function == solvers/function_hdf5_io.hpp save_function_hdf5).
//
// It loads the N functions into a ResponseStateX<ClosedShell> (one Function per
// occupied orbital), wraps them as a 1-root ESSolver<TDA,ClosedShell>::State,
// and calls save_es_roots(..., converged=false). That writes
//   <calc-dir>/es__<protocol_key>/{root_0.*, roots.json}
// and upserts <calc-dir>/response_metadata.json (excited_states/<key>). When the
// normal executor runs its ES node at that protocol, reconcile_protocol finds a
// non-converged exact-key entry -> Resume -> try_load_es_bundle sets s0 to this
// guess -> iterate. Zero solver-code changes; this is the "warm start MADNESS ES
// from a DALTON CIS guess" seam (dalton-cis-seed thread).
//
// Constraints:
//   * NP=1 only. load_function_hdf5 is single-rank (local container == whole
//     tree), and the native seed archive is nio=1, so the ES run that Resumes
//     from it must ALSO be NP=1 (else the writer_nproc np-guard trips).
//   * cell/k must match the solver: load_function_hdf5 calls
//     FunctionDefaults::set_cell(<file cell>), so the pymra projection MUST use
//     the solver's cell (moldft L=200 -> [-200,200]^3) and the protocol's k.
//   * DALTON TDA normalizes 2(||X||^2-||Y||^2)=1 -> ||X||^2=0.5; MADNESS's
//     response metric has no spin factor (<x|x>=1), so we scale by sqrt(2) by
//     default (--scale to override). Magnitude of a guess is not critical, but
//     matching the convention makes it a faithful warm start.
//
// Usage:
//   seed_es_from_hdf5 --calc-dir=DIR --omega=<au> --x=x0.mad.h5,x1.mad.h5,...
//                     [--seed-dir=DIR] [--thresh=1e-4] [--scale=1.41421356]

#include "../ResponseProtocol.hpp"          // protocol_key, default_k_for_thresh
#include "../kernels/tags.hpp"              // TDA, ClosedShell
#include "../solvers/response_state.hpp"    // ResponseStateX<ClosedShell>
#include "../solvers/es_solver.hpp"         // ESSolver<TDA,ClosedShell>::State
#include "../solvers/es_save_load.hpp"      // save_es_roots
#include "../solvers/function_hdf5_io.hpp"  // load_function_hdf5

#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <sstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {
std::vector<std::string> split_csv(const std::string &s) {
  std::vector<std::string> out;
  std::stringstream ss(s);
  std::string t;
  while (std::getline(ss, t, ','))
    if (!t.empty()) out.push_back(t);
  return out;
}
} // namespace

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  startup(world, argc, argv, true);

#ifndef MADNESS_HAS_HDF5
  if (world.rank() == 0)
    print("ERROR: seed_es_from_hdf5 requires an HDF5-enabled build "
          "(-DMADNESS_ENABLE_HDF5).");
  finalize();
  return 2;
#else
  commandlineparser parser(argc, argv);

  if (world.size() != 1) {
    if (world.rank() == 0)
      print("ERROR: seed_es_from_hdf5 must run on a single rank (NP=1).");
    finalize();
    return 2;
  }
  if (!parser.key_exists("calc-dir") || !parser.key_exists("x") ||
      !parser.key_exists("omega")) {
    if (world.rank() == 0) {
      print("Usage: seed_es_from_hdf5 --calc-dir=DIR --omega=<au> "
            "--x=x0.mad.h5,x1.mad.h5,...");
      print("  [--seed-dir=DIR] [--thresh=1e-4] [--scale=<sqrt2>]");
    }
    finalize();
    return 2;
  }

  const std::string calc_dir = parser.value_raw("calc-dir");
  const std::string seed_dir =
      parser.key_exists("seed-dir") ? parser.value_raw("seed-dir") : "";
  const double omega = std::stod(parser.value("omega"));
  const double thresh =
      parser.key_exists("thresh") ? std::stod(parser.value("thresh")) : 1e-4;
  const double scale =
      parser.key_exists("scale") ? std::stod(parser.value("scale"))
                                 : std::sqrt(2.0);
  const auto xfiles = split_csv(parser.value_raw("x"));

  // Protocol the seed belongs to: reconcile_protocol looks the seed up by
  // protocol_key(thresh, default_k_for_thresh(thresh)); save_es_roots records
  // k/thresh into roots.json from FunctionDefaults.
  const int k = default_k_for_thresh(thresh);
  const std::string key = protocol_key(thresh, k);
  const std::string bundle = calc_dir + "/es__" + key;

  // All MADNESS Function/State objects live in this inner scope so their
  // destructors (which touch WorldContainer mutexes) run BEFORE finalize();
  // otherwise teardown hits a dead-mutex abort (RecursiveMutex::lock EINVAL,
  // functions outliving the World).
  {
    // Load each per-occupied-orbital response function. load_function_hdf5 sets
    // FunctionDefaults<3>::set_cell() from the file's stored cell, so all N
    // share one cell (must equal the solver's — see header note).
    ResponseStateX<ClosedShell> root;
    for (const auto &f : xfiles) {
      const std::string path =
          (!seed_dir.empty() && !f.empty() && f.front() != '/')
              ? (seed_dir + "/" + f)
              : f;
      Function<double, 3> fn = load_function_hdf5<double, 3>(world, path);
      if (scale != 1.0) fn.scale(scale);
      root.x_alpha.push_back(fn);
    }
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);

    ESSolver<TDA, ClosedShell>::State s;
    s.roots.push_back(std::move(root));
    s.omega = Tensor<double>(1L);
    s.omega(0L) = omega;
    s.iter = 0;
    // A never-iterated seed has no residuals; save_es_roots writes NaN (-> JSON
    // null) for empty residual vectors, and load_es_roots' value(key, default)
    // THROWS json type_error.302 on a present-but-null number (the default is
    // only used for a MISSING key, not a null one). Stamp a finite,
    // clearly-unconverged sentinel so the bundle round-trips and the ES node
    // Resumes + iterates from this guess.
    s.last_bsh_residual     = std::vector<double>(1, 1.0);
    s.last_density_residual = std::vector<double>(1, 1.0);
    s.last_omega_residual   = std::vector<double>(1, 1.0);

    save_es_roots<TDA, ClosedShell>(world, s, bundle, /*converged=*/false);

    if (world.rank() == 0) {
      print("seed_es_from_hdf5: wrote",
            static_cast<int>(s.roots[0].x_alpha.size()), "orbital functions");
      print("  omega(au) =", omega, "  protocol_key =", key, "  k =", k,
            "  thresh =", thresh, "  scale =", scale);
      print("  bundle    =", bundle);
      print("  Run (NP=1): test_calc_manager_run --archive=<gs> --calc-dir=",
            calc_dir, " --es-roots=1 --protocol=", thresh,
            " -> ES node Resumes from seed.");
    }
  }  // s / root / loaded functions destruct here, before finalize()

  finalize();
  return 0;
#endif
}
