// ===========================================================================
// test_tpa_pq_vs_vbc.cpp — gate for the PROMOTED state-free two-photon
// residue source (P,Q) (kernels/tpa_source_spec.hpp, source-unification
// step 3), and the measured comparison against beta's quadratic source V^{bc}.
//
// The (P,Q) assembly was derived and validated IN this test (see git history:
// the regrouping identities are pinned by test_tpa_regroup_identities, and the
// original inline build agreed with the c-grouped production composition to
// 3e-7 on stored FD responses). It now lives in kernels/tpa_source_spec.hpp
// as data over the declarative engine (kernels/source_spec.hpp); this test is
// the equality gate that promotion must keep passing:
//
//   (1) TWO-ELECTRON gate: <x^f|P> + <y^f|Q> from assemble_tpa_pq(B,C) must
//       equal -(T1+T2+T3) = tpa::tpa_e3_residue(B,C,F)/sqrt(2) to MRA
//       precision. The only use of a state vector; any stand-in works
//       because the regrouping identities are algebraic.
//   (2) ONE-ELECTRON gate: the full assemble_tpa_pq(B,C,vB,vC) adds exactly
//       tpa::tpa_moment_residue_1e's half(b,c)+half(c,b) content:
//       contraction == e3/sqrt(2) + S_1e.
//   (3) COMPARISON (informational): (P,Q) vs V^{bc}'s two-electron part,
//       family norms + difference — the measured statement that the residue
//       source is NOT the VBC source (paper "One Right-Hand Side" §obs 1-3).
//
// Stand-in vectors are used by default (CI needs no external data beyond the
// ground-state fixture). They must NOT be Q-projected and must admix phi with
// irregular coefficients: Q-projected vectors are orthogonal to phi, which
// silently zeroes every <phi|...> factor, and pure Cartesian dipoles on a
// symmetric molecule zero further terms by selection rule. Either makes the
// test pass without testing anything.
//
//   test_tpa_pq_vs_vbc --archive=<moldft restartdata> [--thresh=X] [--k=N]
//       [--calc-dir=<dir> --freq=<omega_f/2> [--baxis=N --caxis=N]]
// ===========================================================================

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/common_ops.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/tpa.hpp"
#include "../kernels/tpa_source_spec.hpp"
#include "../kernels/two_electron.hpp"
#include "../kernels/vbc.hpp"
#include "../solvers/build_response_ground_state.hpp"
#include "../solvers/fd_save_load.hpp"
#include "../solvers/response_state.hpp"

#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

using vecfuncT = std::vector<real_function_3d>;

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    {
      commandlineparser parser(argc, argv);
      if (!parser.key_exists("archive")) {
        if (world.rank() == 0)
          print("Usage: test_tpa_pq_vs_vbc --archive=<path> [--thresh=X] [--k=N]");
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
      const size_t n = phi.size();
      const double cxc = g0.c_xc;

      // ---- stand-in b, c, f vectors (see header note on why not Q-projected)
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
      vecfuncT xb = comb(1.0, d0,  0.23, d1,  0.31);
      vecfuncT yb = comb(1.0, d1, -0.41, d2,  0.17);
      vecfuncT xc = comb(1.0, d2,  0.13, d0, -0.29);
      vecfuncT yc = comb(1.0, d0,  0.53, d1,  0.11);
      const vecfuncT xf = comb(1.0, d1, -0.37, d2, -0.19);
      const vecfuncT yf = comb(1.0, d2,  0.71, d0,  0.43);

      // ---- optionally replace b,c with REAL stored FD responses ----------
      // --calc-dir=<dir> --freq=<omega_f/2> [--baxis=N --caxis=N]
      // The stand-ins above exercise the algebra; real responses give physical
      // magnitudes for the (p,q)-vs-V^{bc} comparison. Falls back to stand-ins
      // with a loud note if the metadata has no record at (pert, freq).
      std::string src = "stand-in (dipole combinations + phi admixture)";
      if (parser.key_exists("calc-dir") && parser.key_exists("freq")) {
        const std::string cdir = parser.value_raw("calc-dir");
        const double freq = std::stod(parser.value("freq"));
        const int ba = parser.key_exists("baxis") ? std::stoi(parser.value("baxis")) : 2;
        const int ca = parser.key_exists("caxis") ? std::stoi(parser.value("caxis")) : 2;
        auto rb = try_load_fd_state<Full, ClosedShell>(
            world, cdir, Perturbation::dipole(ba), freq);
        auto rcc = try_load_fd_state<Full, ClosedShell>(
            world, cdir, Perturbation::dipole(ca), freq);
        if (rb && rcc) {
          xb = madness::copy(world, rb->state.responses[0].x_alpha);
          yb = madness::copy(world, rb->state.responses[0].y_alpha);
          xc = madness::copy(world, rcc->state.responses[0].x_alpha);
          yc = madness::copy(world, rcc->state.responses[0].y_alpha);
          src = "STORED FD responses  dipole_" + std::string(1, "xyz"[ba]) +
                " / dipole_" + std::string(1, "xyz"[ca]) +
                "  @ omega=" + parser.value("freq") +
                "  [b exact=" + (rb->exact ? "yes" : "no") +
                " (" + rb->source_protocol_key + "), c exact=" +
                (rcc->exact ? "yes" : "no") + " (" + rcc->source_protocol_key + ")]";
        } else if (world.rank() == 0) {
          print("  !! no stored FD record for the requested (pert, freq) --"
                " falling back to stand-ins");
        }
      }

      ResponseStateXY<ClosedShell> B, C, F;
      B.x_alpha = madness::copy(world, xb); B.y_alpha = madness::copy(world, yb);
      C.x_alpha = madness::copy(world, xc); C.y_alpha = madness::copy(world, yc);
      F.x_alpha = madness::copy(world, xf); F.y_alpha = madness::copy(world, yf);

      // ============ (P,Q) from the PROMOTED kernel (2e part) ==============
      auto PQ = tpa::assemble_tpa_pq(world, g0, B, C);

      // ============ gate 1: does it reproduce -(T1+T2+T3)? ================
      // tpa_e3_residue returns -sqrt2*(T1+T2+T3), so +e3/sqrt2 = -(T1+T2+T3),
      // which is the quantity (P,Q) was constructed to reproduce.
      const double e3 = tpa::tpa_e3_residue(world, g0, B, C, F);
      const double cgrouped = e3 / std::sqrt(2.0);    // = -(T1+T2+T3)
      const double sgrouped = inner(xf, PQ.x_alpha) + inner(yf, PQ.y_alpha);

      // ============ gate 2: full (1e + 2e) vs residue_1e + e3 =============
      const auto mu_b = dipole_operator(world, 0);
      const auto mu_c = dipole_operator(world, 2);
      auto PQ_full = tpa::assemble_tpa_pq(world, g0, B, C, mu_b, mu_c);
      // arrays laid out so S_1e(0,1) = half(b,c) + half(c,b) with (vB, vC)
      const std::array<ResponseStateXY<ClosedShell>, 3> resp{
          B, C, C};
      const std::array<real_function_3d, 3> ops{mu_b, mu_c, mu_c};
      const auto S1e = tpa::tpa_moment_residue_1e(world, g0, F, resp, ops, 1.0);
      const double ref_full = cgrouped + S1e(0L, 1L);
      const double got_full =
          inner(xf, PQ_full.x_alpha) + inner(yf, PQ_full.y_alpha);

      // ============ gate 4: symmetrized builder == two-ordering sum =======
      // tpa_pq_spec_sym exploits the family-D collapse (D(B,C)==B(C,B), the
      // adjoint identity Gbar=G^T); this gate keeps the CURRENT per-ordering
      // builder as the validation reference for the simplified one.
      {
        auto PQ_cb = tpa::assemble_tpa_pq(world, g0, C, B);
        vecfuncT refx = madness::copy(world, PQ.x_alpha);
        gaxpy(world, 1.0, refx, 1.0, PQ_cb.x_alpha);
        vecfuncT refy = madness::copy(world, PQ.y_alpha);
        gaxpy(world, 1.0, refy, 1.0, PQ_cb.y_alpha);
        auto symr = source_spec::assemble_source(
            world, g0, tpa::tpa_pq_spec_sym(world, g0, B, C));
        truncate(world, symr[0]);
        truncate(world, symr[1]);
        auto nrm = [&](const vecfuncT &v) {
          return std::sqrt(std::abs(inner(world, v, v).sum()));
        };
        vecfuncT dx = madness::copy(world, symr[0]);
        gaxpy(world, 1.0, dx, -1.0, refx);
        vecfuncT dy = madness::copy(world, symr[1]);
        gaxpy(world, 1.0, dy, -1.0, refy);
        const double ex = nrm(dx) / std::max(1e-30, nrm(refx));
        const double ey = nrm(dy) / std::max(1e-30, nrm(refy));
        const bool ok4 = ex < 1.0e-9 && ey < 1.0e-9;
        if (world.rank() == 0)
          print("gate 4 (sym == sum of orderings):  rel dx=", ex,
                " rel dy=", ey, ok4 ? " PASS" : " FAIL");
        if (!ok4) rc = 1;
      }

      // ============ (3) V^{bc} two-electron part (informational) =========
      // Zero one-electron operators isolate the two-electron content.
      real_function_3d zop = madness::copy(phi[0]); zop.scale(0.0);
      auto V = vbc::compute_vbc<ClosedShell>(world, g0, B, C, zop, madness::copy(zop));

      auto nrm = [&](const vecfuncT &v) { return norm2(world, v); };
      vecfuncT dpx = madness::copy(world, PQ.x_alpha); gaxpy(world, 1.0, dpx, -1.0, V.x_alpha);
      vecfuncT dqy = madness::copy(world, PQ.y_alpha); gaxpy(world, 1.0, dqy, -1.0, V.y_alpha);

      if (world.rank() == 0) {
        print("\n=== promoted (P,Q) source (tpa_source_spec)  vs  references ===");
        print("  thresh =", t, "  n_occ =", n, "  c_xc =", cxc);
        print("  b,c source:", src);
        print("\n--- gate 1: two-electron regrouping vs tpa_e3_residue ---");
        printf("  c-grouped  -(T1+T2+T3)          = %+.10f\n", cgrouped);
        printf("  state-grouped <x^f|P>+<y^f|Q>   = %+.10f\n", sgrouped);
        const double rel = std::abs(cgrouped) > 1e-12
                             ? std::abs(sgrouped - cgrouped) / std::abs(cgrouped)
                             : std::abs(sgrouped - cgrouped);
        printf("  relative difference             = %.3e   %s\n", rel,
               rel < 200.0 * t ? "PASS" : "FAIL");
        print("\n--- gate 2: full (1e+2e) vs tpa_moment_residue_1e + e3 ---");
        printf("  reference  e3/sqrt2 + S_1e      = %+.10f\n", ref_full);
        printf("  contraction of full (P,Q)       = %+.10f\n", got_full);
        const double relf = std::abs(ref_full) > 1e-12
                              ? std::abs(got_full - ref_full) / std::abs(ref_full)
                              : std::abs(got_full - ref_full);
        printf("  relative difference             = %.3e   %s\n", relf,
               relf < 200.0 * t ? "PASS" : "FAIL");
        print("\n--- (P,Q) vs V^{bc} two-electron (informational) ---");
        printf("  %-28s ||P||=%10.6f  ||Q||=%10.6f\n", "(P,Q) total", nrm(PQ.x_alpha), nrm(PQ.y_alpha));
        printf("  %-28s ||x||=%10.6f  ||y||=%10.6f\n", "V^{bc} 2-electron", nrm(V.x_alpha), nrm(V.y_alpha));
        printf("  %-28s ||  ||=%10.6f  ||  ||=%10.6f\n", "(P,Q) - V^{bc}", nrm(dpx), nrm(dqy));
        const double relx = nrm(PQ.x_alpha) > 1e-12 ? nrm(dpx) / nrm(PQ.x_alpha) : 0.0;
        const double rely = nrm(PQ.y_alpha) > 1e-12 ? nrm(dqy) / nrm(PQ.y_alpha) : 0.0;
        printf("  relative:  ||P-Vx||/||P|| = %.4f     ||Q-Vy||/||Q|| = %.4f\n", relx, rely);
        rc = (rel < 200.0 * t && relf < 200.0 * t) ? 0 : 1;
        print("\n", rc == 0 ? "PQ_SOURCE_SPEC VALIDATED" : "PQ_SOURCE_SPEC FAILED");
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
