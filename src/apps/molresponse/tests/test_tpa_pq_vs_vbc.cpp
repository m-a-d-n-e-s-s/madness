// ===========================================================================
// test_tpa_pq_vs_vbc.cpp — build the two-photon residue source (p,q) as an
// EXPLICITLY STATE-FREE vector and compare it, term by term, against beta's
// quadratic source V^{bc}.
//
// Both objects are vectors built from phi and the two photon-leg responses
// (b,c) only. That is the whole point: the naive assumption was that they are
// the same object, and this test measures exactly how they differ without any
// excited-state vector entering the comparison.
//
// The c-grouped residue of kernels/tpa.hpp writes the two-electron part as
// -(T1+T2+T3) with x^f/y^f appearing inside densities and exchange pairs. The
// regrouping moves the state out using one master identity, read four ways.
// With (K[a,b]t)_i = sum_j b_j P(a_j t_i), P the Coulomb convolution:
//
//     <u|K[a,b]v> = sum_ij  int int  u_i(r) b_j(r) a_j(r') v_i(r') / |r-r'|
//
// The four vectors sit in distinguishable slots -- u,b at r; a,v at r'; u,v
// carry index i; a,b carry index j -- so ANY one of them can be isolated into
// the bra:
//     isolate v :  <u|K[a,b]v> = <v|K[b,a]u>          (M1, operator transpose)
//     isolate a :  <u|K[a,b]v> = <a|K[u,v]b>          (M3)
//     isolate b :  <u|K[a,b]v> = <b|K[v,u]a>          (M3')
// plus the Coulomb analogue  <u|J[rho]v> = int rho * J[sum_i u_i v_i], and the
// occupied-block move  <z^{uv}|w> = <u|transform(v, <w|phi>)> for
// z^{uv}_i = sum_j phi_j <u_j|v_i>.
//
// Applying these to T1, T2, T3 yields p,q containing no x^f or y^f. The test
// then does two things:
//   (1) VALIDATES the regrouping: <x^f|p> + <y^f|q> must equal -(T1+T2+T3)
//       from the existing tpa_e3_residue, to MRA precision. This is what
//       proves the hand derivation; it is the only place a state vector is
//       used, and any stand-in works because the identities are algebraic.
//   (2) COMPARES (p,q) against V^{bc}'s two-electron part family by family,
//       reporting norms of each and of their differences.
//
// Stand-in vectors are used by default. They must NOT be Q-projected and must
// admix phi with irregular coefficients: Q-projected vectors are orthogonal to
// phi, which silently zeroes every <phi|...> factor, and pure Cartesian dipoles
// on a symmetric molecule zero further terms by selection rule. Either makes
// the test pass without testing anything.
//
//   test_tpa_pq_vs_vbc --archive=<moldft restartdata> [--thresh=X] [--k=N]
// ===========================================================================

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/common_ops.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/tpa.hpp"
#include "../kernels/two_electron.hpp"
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
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;
using molresponse_v3::common_ops::apply_exchange;

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
      const double lo = gs.params().lo();
      const double cxc = g0.c_xc;
      const size_t n = phi.size();

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

      // ---- small helpers -------------------------------------------------
      auto vdot = [&](const vecfuncT &a, const vecfuncT &b) {  // sum_i a_i b_i
        vecfuncT p(a.size());
        for (size_t i = 0; i < a.size(); ++i) p[i] = a[i] * b[i];
        truncate(world, p);
        return sum(world, p);
      };
      auto Jof = [&](const real_function_3d &rho) { return apply(*coulop, rho); };
      auto K = [&](const vecfuncT &bra, const vecfuncT &ket, const vecfuncT &to) {
        return apply_exchange(world, bra, ket, to, lo);
      };
      auto scal = [&](vecfuncT v, double s) { scale(world, v, s); return v; };
      auto add = [&](vecfuncT &dst, const vecfuncT &src, double s) {
        gaxpy(world, 1.0, dst, s, src);
      };
      auto zero = [&]() { return zero_functions<double, 3>(world, (int)n); };
      // vector whose i-th component is phi_i * F
      auto phiF = [&](const real_function_3d &F) { return mul(world, F, phi, true); };
      auto zblk = [&](const vecfuncT &u, const vecfuncT &v) {  // z^{uv}
        return transform(world, phi, matrix_inner(world, u, v), true);
      };

      // ================= (p,q) from the regrouping =======================
      // --- family B : from -T2, the photon-b kernel g_B -------------------
      real_function_3d rhoB = vdot(phi, xb); rhoB += vdot(phi, yb);
      rhoB.scale(2.0); rhoB.truncate();
      auto JB = Jof(rhoB);
      auto gB  = [&](const vecfuncT &tt) {
        return two_electron::apply_gamma_raw(world, JB, tt, {{xb, phi}, {phi, yb}}, cxc, lo);
      };
      auto gBT = [&](const vecfuncT &tt) {   // transposed exchange pairs
        return two_electron::apply_gamma_raw(world, JB, tt, {{phi, xb}, {yb, phi}}, cxc, lo);
      };
      auto MB  = matrix_inner(world, gB(phi),  phi);
      auto MBT = matrix_inner(world, gBT(phi), phi);
      // -T2 = A - B - C + D
      vecfuncT pB = transform(world, xc, MB, true);  add(pB, gBT(xc), -1.0);
      vecfuncT qB = scal(gB(yc), -1.0);              add(qB, transform(world, yc, MBT, true), 1.0);

      // --- family F : from -T3, the regrouped state kernel ----------------
      // Each T3 term becomes  2*phi_j*J[rho] - cxc*K[.,.]phi , state-free.
      auto z1 = zblk(yb, xc);          // z^{y^b x^c}
      auto z2 = zblk(xb, yc);          // z^{x^b y^c}
      auto Jcb = Jof(vdot(xc, yb));
      auto Jbc = Jof(vdot(xb, yc));
      auto Jz1 = Jof(vdot(z1, phi));
      auto Jz2 = Jof(vdot(phi, z2));
      vecfuncT pF = zero(), qF = zero();
      // 3a (T3 sign +) -> subtract
      add(pF, phiF(Jcb), -2.0); add(pF, K(yb, xc, phi), +cxc);
      add(qF, phiF(Jcb), -2.0); add(qF, K(xc, yb, phi), +cxc);
      // 3b (T3 sign -) -> add
      add(pF, phiF(Jz1), +2.0); add(pF, K(phi, z1, phi), -cxc);
      add(qF, phiF(Jz1), +2.0); add(qF, K(z1, phi, phi), -cxc);
      // 3c (T3 sign +) -> subtract
      add(pF, phiF(Jbc), -2.0); add(pF, K(yc, xb, phi), +cxc);
      add(qF, phiF(Jbc), -2.0); add(qF, K(xb, yc, phi), +cxc);
      // 3d (T3 sign -) -> add
      add(pF, phiF(Jz2), +2.0); add(pF, K(z2, phi, phi), -cxc);
      add(qF, phiF(Jz2), +2.0); add(qF, K(phi, z2, phi), -cxc);

      // --- family D : from -T1, the b-with-state transition density -------
      auto W  = Jof(vdot(xc, phi) + vdot(yc, phi));
      auto Nm = matrix_inner(world, phi, mul(world, W, phi, true));
      auto wE3 = K(phi, xc, phi);      // for the zeta^D piece, x-channel
      auto wE6 = K(yc, phi, phi);      // for the zeta^D piece, y-channel
      auto Q3 = matrix_inner(world, phi, wE3);
      auto Q6 = matrix_inner(world, phi, wE6);
      vecfuncT pD = zero(), qD = zero();
      // Coulomb (T1 sign +) -> subtract
      add(pD, mul(world, W, xb, true), -2.0);
      add(pD, transform(world, xb, Nm, true), +2.0);
      add(qD, mul(world, W, yb, true), -2.0);
      add(qD, transform(world, yb, Nm, true), +2.0);
      // exchange, dyad pieces  (T1 carries -cxc; -T1 flips to +cxc)
      add(pD, K(phi, xc, xb), +cxc);
      add(pD, K(yc, phi, xb), +cxc);
      add(qD, K(xc, phi, yb), +cxc);
      add(qD, K(phi, yc, yb), +cxc);
      // exchange, zeta^D pieces (extra minus from nzetaD = -zeta^D)
      add(pD, transform(world, xb, Q3, true), -cxc);
      add(pD, transform(world, xb, Q6, true), -cxc);
      add(qD, transform(world, yb, transpose(Q3), true), -cxc);
      add(qD, transform(world, yb, transpose(Q6), true), -cxc);

      vecfuncT p = madness::copy(world, pB); add(p, pF, 1.0); add(p, pD, 1.0);
      vecfuncT q = madness::copy(world, qB); add(q, qF, 1.0); add(q, qD, 1.0);
      truncate(world, p); truncate(world, q);

      // ================= validation: does it reproduce -(T1+T2+T3)? ======
      ResponseStateXY<ClosedShell> B, C, F;
      B.x_alpha = madness::copy(world, xb); B.y_alpha = madness::copy(world, yb);
      C.x_alpha = madness::copy(world, xc); C.y_alpha = madness::copy(world, yc);
      F.x_alpha = madness::copy(world, xf); F.y_alpha = madness::copy(world, yf);
      // tpa_e3_residue returns -sqrt2*(T1+T2+T3), so +e3/sqrt2 = -(T1+T2+T3),
      // which is the quantity (p,q) was constructed to reproduce.
      const double e3 = tpa::tpa_e3_residue(world, g0, B, C, F);
      const double cgrouped = e3 / std::sqrt(2.0);    // = -(T1+T2+T3)
      const double sgrouped = inner(xf, p) + inner(yf, q);

      // ================= V^{bc} two-electron part =========================
      // Zero one-electron operators isolate the two-electron content.
      real_function_3d zop = madness::copy(phi[0]); zop.scale(0.0);
      auto V = vbc::compute_vbc<ClosedShell>(world, g0, B, C, zop, madness::copy(zop));

      auto nrm = [&](const vecfuncT &v) { return norm2(world, v); };
      vecfuncT dpx = madness::copy(world, p); add(dpx, V.x_alpha, -1.0);
      vecfuncT dqy = madness::copy(world, q); add(dqy, V.y_alpha, -1.0);

      if (world.rank() == 0) {
        print("\n=== (p,q) from the regrouping  vs  V^{bc} two-electron ===");
        print("  thresh =", t, "  n_occ =", n, "  c_xc =", cxc);
        print("  b,c source:", src);
        print("\n--- validation of the regrouping (the only use of a state vector) ---");
        printf("  c-grouped  -(T1+T2+T3)          = %+.10f\n", cgrouped);
        printf("  state-grouped <x^f|p>+<y^f|q>   = %+.10f\n", sgrouped);
        const double rel = std::abs(cgrouped) > 1e-12
                             ? std::abs(sgrouped - cgrouped) / std::abs(cgrouped)
                             : std::abs(sgrouped - cgrouped);
        printf("  relative difference             = %.3e   %s\n", rel,
               rel < 200.0 * t ? "PASS" : "FAIL");
        print("\n--- family norms of (p,q), all state-free ---");
        printf("  %-28s ||p||=%10.6f  ||q||=%10.6f\n", "B  (photon-b kernel g_B)", nrm(pB), nrm(qB));
        printf("  %-28s ||p||=%10.6f  ||q||=%10.6f\n", "F  (regrouped state kern)", nrm(pF), nrm(qF));
        printf("  %-28s ||p||=%10.6f  ||q||=%10.6f\n", "D  (b-with-state density)", nrm(pD), nrm(qD));
        printf("  %-28s ||p||=%10.6f  ||q||=%10.6f\n", "total (p,q)", nrm(p), nrm(q));
        print("\n--- V^{bc} and the difference ---");
        printf("  %-28s ||x||=%10.6f  ||y||=%10.6f\n", "V^{bc} 2-electron", nrm(V.x_alpha), nrm(V.y_alpha));
        printf("  %-28s ||  ||=%10.6f  ||  ||=%10.6f\n", "(p,q) - V^{bc}", nrm(dpx), nrm(dqy));
        const double relx = nrm(p) > 1e-12 ? nrm(dpx) / nrm(p) : 0.0;
        const double rely = nrm(q) > 1e-12 ? nrm(dqy) / nrm(q) : 0.0;
        printf("  relative:  ||p-Vx||/||p|| = %.4f     ||q-Vy||/||q|| = %.4f\n", relx, rely);
        rc = (rel < 200.0 * t) ? 0 : 1;
        print("\n", rc == 0 ? "REGROUPING VALIDATED" : "REGROUPING FAILED");
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
