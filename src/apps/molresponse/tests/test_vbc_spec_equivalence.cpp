// ===========================================================================
// test_vbc_spec_equivalence.cpp — THE gate for source-unification step 2:
// vbc::compute_vbc_spec (the declarative source-spec path over
// kernels/source_spec.hpp) must reproduce vbc::compute_vbc (the v2-ported,
// Dalton-validated bespoke builder) on real MRA functions to summation-order
// precision (per-channel difference norms ~1e-12, i.e. bit-identical up to
// floating-point regrouping — four+ orders below the MRA thresh, so any
// algebraic/index/truncation-point mismatch fails loudly).
//
// Three stages:
//   1. STAND-IN equality: B, C from non-Q-projected dipole+phi admixtures
//      (the pq-test recipe — projected or symmetry-pure vectors silently
//      zero terms and weaken the test). Assert per-channel diff norms.
//   2. CONVERGED equality + beta: solve the static dipole FD states for
//      axes x,z (protocol ramp), rebuild V^{xx} and V^{xz} via BOTH paths,
//      assert the same per-channel equality on converged states, and assert
//      beta_zxx / beta_xxz computed via both paths agree.
//   3. KLEINMAN pair: for static responses beta_zxx == beta_xxz analytically;
//      the residual asymmetry (solver-convergence sized) must MATCH between
//      the two paths — this exercises exactly the (B,C) index plumbing the
//      spec re-encodes.
//
//   test_vbc_spec_equivalence --archive=<moldft restartdata>
//       [--thresh=X] [--k=N] [--protocol=1e-4,1e-6] [--maxiter=N]
//       [--calc-dir=DIR] [--skip-beta]
// ===========================================================================

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../calc/calc_executor.hpp"
#include "../kernels/beta.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/vbc.hpp"
#include "../solvers/build_response_ground_state.hpp"
#include "../solvers/fd_save_load.hpp"
#include "../solvers/response_state.hpp"

#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <optional>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

using vecfuncT = std::vector<real_function_3d>;

namespace {

// Summation-order tolerance for ||V_spec - V_vbc|| per channel: rounding-level
// accumulation over ~10 orbitals x O(10) MRA ops sits at 1e-13..1e-12 for
// O(1)-normed sources; an algebraic or truncation-point mismatch shows at
// >= thresh scale (1e-6 here) — four orders of separation.
constexpr double kSpecTol = 1.0e-10;

struct DiffReport {
  double dx, dy, nx, ny;
};

DiffReport channel_diff(World &world, const ResponseStateXY<ClosedShell> &A,
                        const ResponseStateXY<ClosedShell> &B) {
  vecfuncT dx = madness::copy(world, A.x_alpha);
  gaxpy(world, 1.0, dx, -1.0, B.x_alpha);
  vecfuncT dy = madness::copy(world, A.y_alpha);
  gaxpy(world, 1.0, dy, -1.0, B.y_alpha);
  return {norm2(world, dx), norm2(world, dy),
          norm2(world, A.x_alpha), norm2(world, A.y_alpha)};
}

} // namespace

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    {
      commandlineparser parser(argc, argv);
      if (!parser.key_exists("archive")) {
        if (world.rank() == 0)
          print("Usage: test_vbc_spec_equivalence --archive=<path> "
                "[--thresh=X] [--k=N] [--protocol=1e-4,1e-6] [--maxiter=N] "
                "[--calc-dir=DIR] [--skip-beta]");
        finalize();
        return 1;
      }
      const std::string archive_path = parser.value_raw("archive");
      auto header = GroundState::read_archive_header(world, archive_path);
      const int override_k = parser.key_exists("k") ? std::stoi(parser.value("k")) : -1;
      const double thresh = parser.key_exists("thresh")
                                ? std::stod(parser.value("thresh"))
                                : default_thresh_for_k(header.k);
      set_response_protocol(world, header.L, thresh, override_k);

      Molecule molecule;
      auto dir = std::filesystem::path(archive_path).parent_path();
      for (const auto &name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
        auto cand = dir / name;
        if (std::filesystem::exists(cand)) {
          std::ifstream ifs(cand);
          nlohmann::json j;
          ifs >> j;
          nlohmann::json mj;
          if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
            mj = j["tasks"][0]["molecule"];
          else if (j.contains("molecule")) mj = j["molecule"];
          if (!mj.is_null()) molecule.from_json(mj);
          break;
        }
      }
      auto gs = GroundState::from_archive(world, archive_path, molecule);
      std::string fock_json;
      for (const auto &name : {"moldft.fock.json", "mad.fock.json"}) {
        auto cand = dir / name;
        if (std::filesystem::exists(cand)) { fock_json = cand.string(); break; }
      }
      const double t = FunctionDefaults<3>::get_thresh();
      auto coulop = poperatorT(CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t));
      gs.prepare(world, 0.001 * t, coulop, fock_json);
      auto g0 = build_response_ground_state_closed_shell(
          world, gs, gs.hf_exchange_coefficient(), gs.params().lo());
      const vecfuncT &phi = g0.amo;

      const auto mu_x = dipole_operator(world, 0);
      const auto mu_z = dipole_operator(world, 2);

      bool ok = true;

      // ================= stage 1: stand-in equality ======================
      // Non-Q-projected, phi-admixed stand-ins (see test_tpa_pq_vs_vbc header
      // note): Q-projected or symmetry-pure vectors zero <phi|..> factors and
      // let a wrong spec pass.
      {
        auto d0 = dipole_perturbation(world, gs, 0);
        auto d1 = dipole_perturbation(world, gs, 1);
        auto d2 = dipole_perturbation(world, gs, 2);
        auto comb = [&](double a, const vecfuncT &u, double b, const vecfuncT &v,
                        double cp) {
          vecfuncT r = madness::copy(world, u);
          gaxpy(world, a, r, b, v);
          gaxpy(world, 1.0, r, cp, phi);
          truncate(world, r);
          return r;
        };
        ResponseStateXY<ClosedShell> B, C;
        B.x_alpha = comb(1.0, d0,  0.23, d1,  0.31);
        B.y_alpha = comb(1.0, d1, -0.41, d2,  0.17);
        C.x_alpha = comb(1.0, d2,  0.13, d0, -0.29);
        C.y_alpha = comb(1.0, d0,  0.53, d1,  0.11);

        auto V_ref  = vbc::compute_vbc<ClosedShell>(world, g0, B, C, mu_x, mu_z);
        auto V_spec = vbc::compute_vbc_spec<ClosedShell>(world, g0, B, C, mu_x, mu_z);
        auto d = channel_diff(world, V_ref, V_spec);
        const bool pass = d.dx < kSpecTol && d.dy < kSpecTol;
        ok = ok && pass;
        if (world.rank() == 0) {
          print("\n=== stage 1: compute_vbc_spec == compute_vbc (stand-ins) ===");
          printf("  ||Vx||=%12.8f  ||Vy||=%12.8f\n", d.nx, d.ny);
          printf("  ||dVx||=%.3e  ||dVy||=%.3e   tol=%.1e   %s\n",
                 d.dx, d.dy, kSpecTol, pass ? "PASS" : "FAIL");
        }
      }

      // ================= stage 2+3: converged states, beta, Kleinman =====
      if (!parser.key_exists("skip-beta")) {
        const std::vector<double> protocols =
            parser.key_exists("protocol")
                ? parse_protocol_csv(parser.value("protocol"))
                : std::vector<double>{1e-4, 1e-6};
        const int max_iters =
            parser.key_exists("maxiter") ? std::stoi(parser.value("maxiter")) : 30;
        const std::string calc_dir = parser.key_exists("calc-dir")
                                         ? parser.value_raw("calc-dir")
                                         : std::string("vbc_spec_equiv_calc");
        if (world.rank() == 0)
          std::filesystem::create_directories(calc_dir);
        world.gop.fence();

        ConvergencePolicy policy;
        ExecutorContext ctx(world, gs, header.L, fock_json,
                            ExecutorSettings{policy, PrintLevel::Normal,
                                             calc_dir, max_iters});
        // static dipole responses for the Kleinman pair (zxx vs xxz): axes x,z
        for (int ax : {0, 2}) {
          bool first = true;
          for (double p : protocols) {
            solve_fd_protocol<Static, ClosedShell>(
                ctx, Perturbation::dipole(ax), 0.0, p,
                first ? NodeAction::Fresh : NodeAction::Restart);
            first = false;
          }
        }
        world.gop.fence();

        // reload at the top protocol; rebuild g0 at that protocol
        const double tf = protocols.back();
        set_response_protocol(world, header.L, tf, override_k);
        {
          const double t0 = FunctionDefaults<3>::get_thresh();
          auto cop = poperatorT(CoulombOperatorPtr(world, gs.params().lo(), 0.001 * t0));
          gs.prepare(world, 0.001 * t0, cop, fock_json);
        }
        auto g0f = build_response_ground_state_closed_shell(
            world, gs, gs.hf_exchange_coefficient(), gs.params().lo());

        auto load_static_xy = [&](int ax) -> std::optional<ResponseStateXY<ClosedShell>> {
          auto r = try_load_fd_state<Static, ClosedShell>(
              world, calc_dir, Perturbation::dipole(ax), 0.0);
          if (!r) return std::nullopt;
          ResponseStateXY<ClosedShell> xy;
          xy.x_alpha = madness::copy(world, r->state.responses[0].x_alpha);
          xy.y_alpha = madness::copy(world, r->state.responses[0].x_alpha);  // static: Y = X
          return xy;
        };
        auto X = load_static_xy(0);
        auto Z = load_static_xy(2);
        if (!X || !Z) {
          if (world.rank() == 0)
            print("  FAIL: could not load solved static FD states from", calc_dir);
          ok = false;
        } else {
          // V^{xx} (for beta_zxx) and V^{xz} (for beta_xxz), both paths
          auto Vxx_ref  = vbc::compute_vbc<ClosedShell>(world, g0f, *X, *X, mu_x, mu_x);
          auto Vxx_spec = vbc::compute_vbc_spec<ClosedShell>(world, g0f, *X, *X, mu_x, mu_x);
          auto Vxz_ref  = vbc::compute_vbc<ClosedShell>(world, g0f, *X, *Z, mu_x, mu_z);
          auto Vxz_spec = vbc::compute_vbc_spec<ClosedShell>(world, g0f, *X, *Z, mu_x, mu_z);
          auto dxx = channel_diff(world, Vxx_ref, Vxx_spec);
          auto dxz = channel_diff(world, Vxz_ref, Vxz_spec);
          const bool pass_conv = dxx.dx < kSpecTol && dxx.dy < kSpecTol &&
                                 dxz.dx < kSpecTol && dxz.dy < kSpecTol;

          const double b_zxx_ref  = beta::beta_abc<ClosedShell>(world, g0f, *Z, Vxx_ref,  *X, *X, mu_z);
          const double b_zxx_spec = beta::beta_abc<ClosedShell>(world, g0f, *Z, Vxx_spec, *X, *X, mu_z);
          const double b_xxz_ref  = beta::beta_abc<ClosedShell>(world, g0f, *X, Vxz_ref,  *X, *Z, mu_x);
          const double b_xxz_spec = beta::beta_abc<ClosedShell>(world, g0f, *X, Vxz_spec, *X, *Z, mu_x);
          const double dbeta = std::max(std::abs(b_zxx_spec - b_zxx_ref),
                                        std::abs(b_xxz_spec - b_xxz_ref));
          const bool pass_beta = dbeta < 1e-8;

          // Kleinman: static beta_zxx == beta_xxz analytically; the residual
          // is solver-convergence sized and must be path-independent.
          const double scale_k = std::max({std::abs(b_zxx_ref), std::abs(b_xxz_ref), 1e-6});
          const double kle_ref  = std::abs(b_zxx_ref  - b_xxz_ref)  / scale_k;
          const double kle_spec = std::abs(b_zxx_spec - b_xxz_spec) / scale_k;
          const bool pass_kle = kle_ref < 1e-2 && kle_spec < 1e-2 &&
                                std::abs(kle_ref - kle_spec) < 1e-6;

          ok = ok && pass_conv && pass_beta && pass_kle;
          if (world.rank() == 0) {
            print("\n=== stage 2: equality on CONVERGED static responses ===");
            printf("  V^{xx}: ||dVx||=%.3e  ||dVy||=%.3e   (||Vx||=%.6f, ||Vy||=%.6f)\n",
                   dxx.dx, dxx.dy, dxx.nx, dxx.ny);
            printf("  V^{xz}: ||dVx||=%.3e  ||dVy||=%.3e   (||Vx||=%.6f, ||Vy||=%.6f)\n",
                   dxz.dx, dxz.dy, dxz.nx, dxz.ny);
            printf("  per-channel tol=%.1e   %s\n", kSpecTol, pass_conv ? "PASS" : "FAIL");
            print("\n=== stage 2b: beta via both paths ===");
            printf("  beta_zxx: ref=%+.10f  spec=%+.10f  |d|=%.3e\n",
                   b_zxx_ref, b_zxx_spec, std::abs(b_zxx_spec - b_zxx_ref));
            printf("  beta_xxz: ref=%+.10f  spec=%+.10f  |d|=%.3e\n",
                   b_xxz_ref, b_xxz_spec, std::abs(b_xxz_spec - b_xxz_ref));
            printf("  max path diff = %.3e  (tol 1e-8)   %s\n", dbeta,
                   pass_beta ? "PASS" : "FAIL");
            print("\n=== stage 3: Kleinman pair (zxx vs xxz, static) ===");
            printf("  rel asym ref=%.3e  spec=%.3e  (tol 1e-2 each; paths must match)  %s\n",
                   kle_ref, kle_spec, pass_kle ? "PASS" : "FAIL");
          }
        }
      }

      if (world.rank() == 0) {
        print("\nVBC_SPEC_EQUIVALENCE:", ok ? "PASSED" : "FAILED");
        rc = ok ? 0 : 1;
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
