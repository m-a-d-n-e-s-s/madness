// ===========================================================================
// test_vbc_spec_equivalence.cpp — THE gate that the two builders of the one
// second-order source agree: vbc::compute_vbc (kernels/vbc.hpp, the spec of
// Hurtado eq 19) and tpa::quadratic_source (kernels/tpa_source_spec.hpp, the
// (P,Q) builder validated against DALTON for 2PA and SHG). They differ only in
// bookkeeping — (P,Q) leaves its apply moves un-Q-projected and builds the
// pair-density family leg by leg — so the comparison is Q(P,Q) vs V^{BC}, at
// MRA thresh scale (per-leg untruncated densities regroup at ~thresh, not at
// roundoff; an orientation or index error shows at O(0.1-1) relative).
//
// Stages:
//   1.  STAND-IN equality on NON-HERMITIAN legs (B.x != B.y): the
//       orientation-sensitive check. Mixed-axis dipole stand-ins, Q-PROJECTED:
//       response vectors live in the virtual space, and both builders assume
//       it — V^{BC}'s occupied-matrix term [M] is Sum_k x^C_k F^B_kp, i.e. in
//       the span of x^C, so a phi admixture in x^C survives in V but is
//       removed by the Q applied to (P,Q) (job 2161695: rel 0.62 / 0.15 and
//       0.77 in the Hermitian limit on un-projected stand-ins, while the
//       production recontraction of converged Q-space legs agreed to 1e-6).
//       Symmetry-mixing (x, y, z admixed) keeps every term alive; x != y
//       keeps the orientation visible.
//   1b. HERMITIAN LIMIT (y := x on both legs): must also agree — before
//       2026-09-09 the two builders agreed ONLY here (HANDOFF §2c).
//   2.  CONVERGED equality + beta: solve the static dipole FD states for
//       axes x,z (protocol ramp), rebuild V^{xx} and V^{xz} via BOTH builders,
//       assert equality on converged states and that beta_zxx / beta_xxz
//       agree (static: orientation-blind, checks the rest of the plumbing).
//   3.  KLEINMAN pair: for static responses beta_zxx == beta_xxz analytically;
//       the residual asymmetry (solver-convergence sized) must MATCH between
//       the two builders.
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
#include "../kernels/tpa_source_spec.hpp"
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

// Per-channel RELATIVE tolerance for ||Q(P) - Vx|| / ||Vx||. The two builders
// share the engine but not the truncation points of the pair-density family
// (vbc: one truncated density; (P,Q): four untruncated ones), so they regroup
// at ~thresh, not at roundoff. An orientation or index error shows at
// O(0.1-1) relative (HANDOFF §2a: 0.41 / 0.23 before the fix) — the gate
// sits at 10 x thresh, four+ orders below that.
constexpr double kSpecRelTolFactor = 10.0;

struct DiffReport {
  double dx, dy, nx, ny;
};

/// ||A - B|| per channel and ||A|| per channel.
DiffReport channel_diff(World &world, const ResponseStateXY<ClosedShell> &A,
                        const ResponseStateXY<ClosedShell> &B) {
  vecfuncT dx = madness::copy(world, A.x_alpha);
  gaxpy(world, 1.0, dx, -1.0, B.x_alpha);
  vecfuncT dy = madness::copy(world, A.y_alpha);
  gaxpy(world, 1.0, dy, -1.0, B.y_alpha);
  return {norm2(world, dx), norm2(world, dy),
          norm2(world, A.x_alpha), norm2(world, A.y_alpha)};
}

/// Q-project both channels of a source (the part every contraction sees).
ResponseStateXY<ClosedShell> project_source(World &world, const ResponseGroundState &g0,
                                            ResponseStateXY<ClosedShell> s) {
  s.x_alpha = g0.Qa(s.x_alpha);
  s.y_alpha = g0.Qa(s.y_alpha);
  truncate(world, s.x_alpha);
  truncate(world, s.y_alpha);
  return s;
}

/// Q-projected (P,Q) source ((P,Q) leaves its apply moves un-Q-projected).
ResponseStateXY<ClosedShell>
projected_pq(World &world, const ResponseGroundState &g0,
             const ResponseStateXY<ClosedShell> &B,
             const ResponseStateXY<ClosedShell> &C,
             const real_function_3d &VB_op, const real_function_3d &VC_op) {
  return project_source(world, g0,
                        tpa::quadratic_source(world, g0, B, C, VB_op, VC_op));
}

/// Q-projected V^{BC} (a no-op on Q-space inputs; keeps the gate symmetric).
ResponseStateXY<ClosedShell>
projected_vbc(World &world, const ResponseGroundState &g0,
              const ResponseStateXY<ClosedShell> &B,
              const ResponseStateXY<ClosedShell> &C,
              const real_function_3d &VB_op, const real_function_3d &VC_op) {
  return project_source(world, g0,
                        vbc::compute_vbc<ClosedShell>(world, g0, B, C, VB_op, VC_op));
}

/// Q-project a response-shaped pair in place (stand-ins into the virtual space).
void project_state(const ResponseGroundState &g0, ResponseStateXY<ClosedShell> &s) {
  s.x_alpha = g0.Qa(s.x_alpha);
  s.y_alpha = g0.Qa(s.y_alpha);
}

/// Gate: V (vbc builder) vs Q(P,Q), both channels, relative to ||V||.
bool builders_agree(World &world, const char *label,
                    const ResponseStateXY<ClosedShell> &V,
                    const ResponseStateXY<ClosedShell> &PQ, double tol) {
  auto d = channel_diff(world, V, PQ);
  const double rx = d.nx > 1e-14 ? d.dx / d.nx : d.dx;
  const double ry = d.ny > 1e-14 ? d.dy / d.ny : d.dy;
  const bool pass = rx < tol && ry < tol;
  if (world.rank() == 0) {
    printf("  %-22s ||Vx||=%12.8f  ||Vy||=%12.8f\n", label, d.nx, d.ny);
    printf("  %-22s rel ||Q(P)-Vx||=%.3e  rel ||Q(Q)-Vy||=%.3e   tol=%.1e   %s\n",
           "", rx, ry, tol, pass ? "PASS" : "FAIL");
  }
  return pass;
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
      // Symmetry-mixed dipole stand-ins, Q-projected (see header: the
      // builders are equal on the response-vector domain, i.e. Q-space).
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
        const double tol = kSpecRelTolFactor * t;
        ResponseStateXY<ClosedShell> B, C;
        B.x_alpha = comb(1.0, d0,  0.23, d1,  0.31);
        B.y_alpha = comb(1.0, d1, -0.41, d2,  0.17);
        C.x_alpha = comb(1.0, d2,  0.13, d0, -0.29);
        C.y_alpha = comb(1.0, d0,  0.53, d1,  0.11);
        project_state(g0, B);
        project_state(g0, C);

        if (world.rank() == 0)
          print("\n=== stage 1: Q(compute_vbc) == Q(tpa::quadratic_source), "
                "NON-Hermitian Q-space stand-ins (x != y) ===");
        auto V  = projected_vbc(world, g0, B, C, mu_x, mu_z);
        auto PQ = projected_pq(world, g0, B, C, mu_x, mu_z);
        const bool pass1 = builders_agree(world, "x!=y", V, PQ, tol);

        // 1b: the Hermitian limit — the only regime where the two agreed
        // before the orientation fix; must still agree.
        ResponseStateXY<ClosedShell> Bh, Ch;
        Bh.x_alpha = madness::copy(world, B.x_alpha); Bh.y_alpha = madness::copy(world, B.x_alpha);
        Ch.x_alpha = madness::copy(world, C.x_alpha); Ch.y_alpha = madness::copy(world, C.x_alpha);
        if (world.rank() == 0)
          print("\n=== stage 1b: Hermitian limit (y := x) ===");
        auto Vh  = projected_vbc(world, g0, Bh, Ch, mu_x, mu_z);
        auto PQh = projected_pq(world, g0, Bh, Ch, mu_x, mu_z);
        const bool pass1b = builders_agree(world, "y:=x", Vh, PQh, tol);
        ok = ok && pass1 && pass1b;
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
          // V^{xx} (for beta_zxx) and V^{xz} (for beta_xxz), both builders
          const double tolf = kSpecRelTolFactor * tf;
          auto Vxx_ref  = projected_vbc(world, g0f, *X, *X, mu_x, mu_x);
          auto Vxx_spec = projected_pq(world, g0f, *X, *X, mu_x, mu_x);
          auto Vxz_ref  = projected_vbc(world, g0f, *X, *Z, mu_x, mu_z);
          auto Vxz_spec = projected_pq(world, g0f, *X, *Z, mu_x, mu_z);
          if (world.rank() == 0)
            print("\n=== stage 2: equality on CONVERGED static responses ===");
          const bool pass_conv = builders_agree(world, "V^{xx}", Vxx_ref, Vxx_spec, tolf) &&
                                 builders_agree(world, "V^{xz}", Vxz_ref, Vxz_spec, tolf);
          auto dxx = channel_diff(world, Vxx_ref, Vxx_spec);
          auto dxz = channel_diff(world, Vxz_ref, Vxz_spec);

          const double b_zxx_ref  = beta::beta_abc<ClosedShell>(world, g0f, *Z, Vxx_ref,  *X, *X, mu_z);
          const double b_zxx_spec = beta::beta_abc<ClosedShell>(world, g0f, *Z, Vxx_spec, *X, *X, mu_z);
          const double b_xxz_ref  = beta::beta_abc<ClosedShell>(world, g0f, *X, Vxz_ref,  *X, *Z, mu_x);
          const double b_xxz_spec = beta::beta_abc<ClosedShell>(world, g0f, *X, Vxz_spec, *X, *Z, mu_x);
          // Contractions are Q-blind and see only the (thresh-sized) regrouping
          // of the pair-density family: gate at thresh relative to |beta|.
          const double dbeta = std::max(std::abs(b_zxx_spec - b_zxx_ref),
                                        std::abs(b_xxz_spec - b_xxz_ref));
          const double beta_scale = std::max({std::abs(b_zxx_ref), std::abs(b_xxz_ref), 1.0});
          const bool pass_beta = dbeta < tolf * beta_scale;

          // Kleinman: static beta_zxx == beta_xxz analytically; the residual
          // is solver-convergence sized and must be path-independent.
          const double scale_k = std::max({std::abs(b_zxx_ref), std::abs(b_xxz_ref), 1e-6});
          const double kle_ref  = std::abs(b_zxx_ref  - b_xxz_ref)  / scale_k;
          const double kle_spec = std::abs(b_zxx_spec - b_xxz_spec) / scale_k;
          const bool pass_kle = kle_ref < 1e-2 && kle_spec < 1e-2 &&
                                std::abs(kle_ref - kle_spec) < 1e-6;

          ok = ok && pass_conv && pass_beta && pass_kle;
          if (world.rank() == 0) {
            printf("  V^{xx}: ||dVx||=%.3e  ||dVy||=%.3e   (||Vx||=%.6f, ||Vy||=%.6f)\n",
                   dxx.dx, dxx.dy, dxx.nx, dxx.ny);
            printf("  V^{xz}: ||dVx||=%.3e  ||dVy||=%.3e   (||Vx||=%.6f, ||Vy||=%.6f)\n",
                   dxz.dx, dxz.dy, dxz.nx, dxz.ny);
            printf("  stage 2 %s\n", pass_conv ? "PASS" : "FAIL");
            print("\n=== stage 2b: beta via both builders (vbc vs (P,Q)) ===");
            printf("  beta_zxx: vbc=%+.10f  pq=%+.10f  |d|=%.3e\n",
                   b_zxx_ref, b_zxx_spec, std::abs(b_zxx_spec - b_zxx_ref));
            printf("  beta_xxz: vbc=%+.10f  pq=%+.10f  |d|=%.3e\n",
                   b_xxz_ref, b_xxz_spec, std::abs(b_xxz_spec - b_xxz_ref));
            printf("  max builder diff = %.3e  (tol %.1e)   %s\n", dbeta,
                   tolf * beta_scale, pass_beta ? "PASS" : "FAIL");
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
