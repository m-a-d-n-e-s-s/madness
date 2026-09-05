// ===========================================================================
// test_vbc_pq_terms.cpp — TERM-RESOLVED comparison of the two quadratic
// sources: V^{BC} (kernels/vbc.hpp, spec table vbc_half_spec) against the 2PA
// residue (P,Q) (kernels/tpa_source_spec.hpp), entry by entry.
//
// Motivation (2026-09-04 review): the contracted scalars agree with each
// other and with DALTON, but a contracted scalar can hide compensation
// between terms. This test evaluates EVERY SourceEntry of both tables
// separately (singleton specs over the same engine), aligns the slots that
// correspond under the compact-equation dictionary, and prints per term:
//     ||term||_2, the aligned difference, and the contraction with a
//     stand-in (or real) state vector <x^f|term> / <y^f|term>.
// So "the sums of the individual terms" are visible, not just their total.
//
// Alignment (P-channel vs V^{BC} X-channel; Q vs Y is the conjugate):
//     V.fphi  [+Σ_k x^C F^B_kp]    <-> P.B_mat [+Σ_k x^C F̄^B_kp]  (dagger)
//     V.fb    [-Q̂ F^B x^C]         <-> P.B_app [-g'[γ^{B†}]x^C]   (dagger,
//                                       and P omits Q̂ + carries no v — legal
//                                       only under a Q̂-projected contraction)
//     V.gzeta [-Q̂ g'[γ_ζ]φ]        <-> P.F_1..4 [-g'[D^{BC}]φ]     (reorient)
//     (nothing)                    <-> P.D_app + P.D_mat            (R family)
// The 1e family (v^B / v^C content) is evaluated separately so the question
// "where did v^C go relative to Σ x^B G^C_kp" gets an explicit numeric row.
//
// Also probes the beta contraction convention: b1 = -(<xA|Vx> + <yA|Vy>)
// (production, kernels/beta.hpp) against the x<->y-swapped pairing
// <yA|Vx> + <xA|Vy> that the negative-frequency exchange rule would suggest —
// printing both makes the convention knot measurable instead of argued.
//
// Symmetry caveat: on a symmetric molecule (h2o fixture) selection rules zero
// several rows even with the phi-admixed stand-ins. Run --archive against a
// C1 (no-symmetry) ground state for fully dense rows; see the reconciliation
// report ledger (reports/2026-09-04_report_pq_vbc_reconciliation).
//
//   test_vbc_pq_terms --archive=<moldft restartdata> [--thresh=X] [--k=N]
//       [--calc-dir=<dir> --freq=<omega_f/2> [--baxis=N --caxis=N]]
// ===========================================================================

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/common_ops.hpp"
#include "../kernels/source_spec.hpp"
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
#include <string>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

using vecfuncT = std::vector<real_function_3d>;

namespace {

double vnorm(World &world, const vecfuncT &v) {
  return std::sqrt(std::abs(inner(world, v, v).sum()));
}
double vinner(World &world, const vecfuncT &a, const vecfuncT &b) {
  return inner(world, a, b).sum();
}
vecfuncT vdiff(World &world, const vecfuncT &a, const vecfuncT &b) {
  vecfuncT d = madness::copy(world, a);
  gaxpy(world, 1.0, d, -1.0, b);
  return d;
}
vecfuncT vsum(World &world, const std::vector<vecfuncT> &terms,
              const std::vector<int> &idx) {
  vecfuncT s = zero_functions_compressed<double, 3>(
      world, static_cast<int>(terms.at(idx.at(0)).size()));
  for (int i : idx) gaxpy(world, 1.0, s, 1.0, terms.at(static_cast<size_t>(i)));
  return s;
}

/// Evaluate every entry of every channel separately (singleton specs).
/// out[channel][entry] = the entry's vecfuncT, sign included.
std::vector<std::vector<vecfuncT>>
eval_per_entry(World &world, const ResponseGroundState &g0,
               const std::vector<source_spec::SourceSpec> &channels) {
  std::vector<std::vector<vecfuncT>> out;
  for (const auto &ch : channels) {
    std::vector<vecfuncT> terms;
    for (const auto &e : ch.entries) {
      source_spec::SourceSpec one;
      one.entries.push_back(e);
      auto r = source_spec::assemble_source(world, g0, {one});
      terms.push_back(std::move(r[0]));
    }
    out.push_back(std::move(terms));
  }
  return out;
}

void row(World &world, const char *label, const vecfuncT *v, const vecfuncT *p,
         const vecfuncT &f) {
  if (v && p) {
    auto d = vdiff(world, *p, *v);
    if (world.rank() == 0)
      printf("  %-34s |V|=%11.4e |P|=%11.4e |P-V|=%11.4e  <f|V>=%+13.6e <f|P>=%+13.6e\n",
             label, vnorm(world, const_cast<vecfuncT &>(*v)),
             vnorm(world, const_cast<vecfuncT &>(*p)), vnorm(world, d),
             vinner(world, f, *v), vinner(world, f, *p));
  } else if (p) {
    if (world.rank() == 0)
      printf("  %-34s %-13s |P|=%11.4e %-19s %-21s <f|P>=%+13.6e\n", label, "",
             vnorm(world, const_cast<vecfuncT &>(*p)), "", "",
             vinner(world, f, *p));
  }
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
          print("Usage: test_vbc_pq_terms --archive=<path> [--thresh=X] [--k=N]");
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

      // ---- stand-in b, c, f vectors (same recipe as test_tpa_pq_vs_vbc:
      // NOT Q-projected, irregular phi admixture — see that test's header)
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
      vecfuncT xf = comb(1.0, d1, -0.37, d2, -0.19);
      vecfuncT yf = comb(1.0, d2,  0.71, d0,  0.43);

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
          src = "STORED FD responses @ omega=" + parser.value("freq");
        } else if (world.rank() == 0) {
          print("  !! no stored FD record at (pert,freq) — using stand-ins");
        }
      }
      if (world.rank() == 0) {
        print("\n[TERMS] source vectors:", src);
        print("[TERMS] archive:", archive_path, " thresh:", t,
              " k:", FunctionDefaults<3>::get_k(), " c_x:", g0.c_xc);
      }

      // ================= build BOTH tables, per-entry =====================
      // V^{BC}: both ordering halves; zero 1e op isolates two-electron parts.
      real_function_3d zop = madness::copy(phi[0]); zop.scale(0.0);
      auto zeta_bc = vbc::make_zeta(world, yb, xc, phi);
      auto zeta_cb = vbc::make_zeta(world, yc, xb, phi);
      auto Vbc = eval_per_entry(world, g0,
          vbc::vbc_half_spec(world, g0, xb, yb, xc, yc, zeta_bc, zop));
      auto Vcb = eval_per_entry(world, g0,
          vbc::vbc_half_spec(world, g0, xc, yc, xb, yb, zeta_cb, zop));
      // (P,Q): both orderings, 2e only (no 1e ops) — entry layout per
      // tpa_source_spec.hpp: P[0]=B_mat P[1]=B_app P[2..5]=F legs
      // P[6]=D_app P[7]=D_mat;   Q[0]=B_app Q[1]=B_mat Q[2..5]=F Q[6]=D_app
      // Q[7]=D_mat.
      ResponseStateXY<ClosedShell> B, C;
      B.x_alpha = madness::copy(world, xb); B.y_alpha = madness::copy(world, yb);
      C.x_alpha = madness::copy(world, xc); C.y_alpha = madness::copy(world, yc);
      auto Pbc = eval_per_entry(world, g0, tpa::tpa_pq_spec(world, g0, B, C));
      auto Pcb = eval_per_entry(world, g0, tpa::tpa_pq_spec(world, g0, C, B));

      // Ordering-summed slots. V channels: [0]=X entries {gzeta, fb, fphi},
      // [1]=Y. Sum (B,C)+(C,B) per slot.
      auto vslot = [&](int ch, int e) {
        vecfuncT s = madness::copy(world, Vbc[ch][e]);
        gaxpy(world, 1.0, s, 1.0, Vcb[ch][e]);
        return s;
      };
      auto pslot = [&](int ch, std::vector<int> idx) {
        auto a = vsum(world, Pbc[static_cast<size_t>(ch)], idx);
        auto b = vsum(world, Pcb[static_cast<size_t>(ch)], idx);
        gaxpy(world, 1.0, a, 1.0, b);
        return a;
      };

      // ===================== the aligned table ============================
      if (world.rank() == 0)
        print("\n[TERMS] ===== P-channel vs V^BC X-channel (both orderings summed, 2e only) =====");
      {
        auto v_fphi  = vslot(0, 2);              // +Σ x^C F^B_kp (+ BC image)
        auto p_bmat  = pslot(0, {0});            // +Σ x^C F̄^B_kp (+ image)
        row(world, "property-Fock MATRIX (fphi|B_mat)", &v_fphi, &p_bmat, xf);

        auto v_fb    = vslot(0, 1);              // -Q̂ F^B x^C (+ image)
        auto p_bapp  = pslot(0, {1});            // -g'[γ^B†] x^C (+ image)
        row(world, "property-Fock APPLY  (fb|B_app)", &v_fb, &p_bapp, xf);

        auto v_gz    = vslot(0, 0);              // -Q̂ g'[γ_ζ] φ (+ image)
        auto p_F     = pslot(0, {2, 3, 4, 5});   // -g'[D^{BC}] φ (+ image)
        row(world, "pair density (gzeta|D-legs)", &v_gz, &p_F, xf);

        auto p_Dapp  = pslot(0, {6});            // -g'[γ^C†] x^B (+ image)
        auto p_Dmat  = pslot(0, {7});            // +Σ x^B G^C_kp (+ image)
        row(world, "R family APPLY  -g'[γ^C†]x^B", nullptr, &p_Dapp, xf);
        row(world, "R family MATRIX +Σ x^B G^C_kp", nullptr, &p_Dmat, xf);
        auto p_R = pslot(0, {6, 7});
        row(world, "R family TOTAL", nullptr, &p_R, xf);

        // totals + engine-sum sanity
        auto v_tot = vsum(world, {v_gz, v_fb, v_fphi}, {0, 1, 2});
        auto p_tot = vsum(world, {p_bmat, p_bapp, p_F, p_R}, {0, 1, 2, 3});
        row(world, "TOTAL (V^BC.x | P)", &v_tot, &p_tot, xf);
      }

      if (world.rank() == 0)
        print("\n[TERMS] ===== Q-channel vs V^BC Y-channel =====");
      {
        auto v_fphi = vslot(1, 2);
        auto q_bmat = pslot(1, {1});             // Q entry order: [0]=B_app
        row(world, "property-Fock MATRIX", &v_fphi, &q_bmat, yf);
        auto v_fb   = vslot(1, 1);
        auto q_bapp = pslot(1, {0});
        row(world, "property-Fock APPLY", &v_fb, &q_bapp, yf);
        auto v_gz   = vslot(1, 0);
        auto q_F    = pslot(1, {2, 3, 4, 5});
        row(world, "pair density", &v_gz, &q_F, yf);
        auto q_R = pslot(1, {6, 7});
        row(world, "R family TOTAL", nullptr, &q_R, yf);
        auto v_tot = vsum(world, {v_gz, v_fb, v_fphi}, {0, 1, 2});
        auto q_tot = vsum(world, {q_bmat, q_bapp, q_F, q_R}, {0, 1, 2, 3});
        row(world, "TOTAL (V^BC.y | Q)", &v_tot, &q_tot, yf);
      }

      // ============ the v^C question: the 1e family, explicitly ==========
      // Full spec with real dipole ops; entries beyond the 8 2e ones are the
      // 1e family: P[8]=-Q̂v^B x^C, P[9]=+Σx^C <φ|v^B|φ>, P[10]=-Q̂v^C x^B,
      // P[11]=+Σx^B <φ|v^C|φ>. Rows 10/11 are the (C,B) 1e image — the v^C
      // content whose ABSENCE from the R family (G^C vs F^C) was the puzzle.
      if (world.rank() == 0)
        print("\n[TERMS] ===== 1e family (where v^C actually lives) =====");
      {
        const auto mu_b = dipole_operator(world, 0);
        const auto mu_c = dipole_operator(world, 2);
        auto Pfull = eval_per_entry(world, g0,
            tpa::tpa_pq_spec(world, g0, B, C, mu_b, mu_c));
        const auto &PT = Pfull[0];
        if (PT.size() >= 12) {
          row(world, "1e: -Q̂ v^B x^C", nullptr, &PT[8], xf);
          row(world, "1e: +Σ x^C <φ|v^B|φ>_kp", nullptr, &PT[9], xf);
          row(world, "1e: -Q̂ v^C x^B  (the missing v^C)", nullptr, &PT[10], xf);
          row(world, "1e: +Σ x^B <φ|v^C|φ>_kp", nullptr, &PT[11], xf);
          // R + 1e(C,B) = would-be full property-Fock at mixed transposition
          vecfuncT rfull = madness::copy(world, PT[10]);
          gaxpy(world, 1.0, rfull, 1.0, PT[11]);
          auto p_R = pslot(0, {6, 7});
          gaxpy(world, 1.0, rfull, 1.0, p_R);
          row(world, "R + 1e(C,B) combined", nullptr, &rfull, xf);
        } else if (world.rank() == 0) {
          print("  !! unexpected 1e entry layout (", PT.size(), "entries)");
        }
      }

      // ============ conjugation check: Q == S[P] (x<->y swap on legs) ======
      // Claim (2026-09-05, user observation): the Q channel is EXACTLY the P
      // channel with every response leg's x and y halves swapped (the
      // negative-frequency exchange rule at the source level) — one generator
      // functional, evaluated at +/- frequency legs. Verified here numerically;
      // it also holds entry-for-entry in tpa_pq_spec's tables.
      if (world.rank() == 0)
        print("\n[TERMS] ===== conjugation: Q vs P at swapped legs =====");
      {
        ResponseStateXY<ClosedShell> Bs, Cs;   // x<->y swapped legs
        Bs.x_alpha = madness::copy(world, yb); Bs.y_alpha = madness::copy(world, xb);
        Cs.x_alpha = madness::copy(world, yc); Cs.y_alpha = madness::copy(world, xc);
        auto PQswap = eval_per_entry(world, g0, tpa::tpa_pq_spec(world, g0, Bs, Cs));
        auto Pswap_tot = vsum(world, PQswap[0], {0,1,2,3,4,5,6,7});
        auto Q_tot     = vsum(world, Pbc[1],    {0,1,2,3,4,5,6,7});
        auto d = vdiff(world, Pswap_tot, Q_tot);
        if (world.rank() == 0)
          printf("  ||P[swapped legs] - Q|| = %.3e   (||Q|| = %.3e)\n",
                 vnorm(world, d), vnorm(world, Q_tot));
        // Reversed contraction pairing <y|P>+<x|Q>: the -omega_f (emission-
        // side) residue pairing. Informational — differs on stand-ins; on a
        // true eigenvector its magnitude should match the absorption pairing.
        auto P_tot = vsum(world, Pbc[0], {0,1,2,3,4,5,6,7});
        if (world.rank() == 0)
          printf("  pairing  <x|P>+<y|Q> = %+.8e    reversed <y|P>+<x|Q> = %+.8e\n",
                 vinner(world, xf, P_tot) + vinner(world, yf, Q_tot),
                 vinner(world, yf, P_tot) + vinner(world, xf, Q_tot));
      }

      // ============ beta contraction convention probe =====================
      // Production (kernels/beta.hpp): b1 = -(<xA|Vx> + <yA|Vy>) with the A
      // leg SOLVED AT +omega_sigma. The x<->y exchange-rule pairing would be
      // <yA|Vx> + <xA|Vy>. Both printed; (xf,yf) stands in for the A leg.
      if (world.rank() == 0)
        print("\n[TERMS] ===== beta contraction convention probe =====");
      {
        auto Vfull = vbc::compute_vbc_spec<ClosedShell>(world, g0, B, C, zop,
                                                        madness::copy(zop));
        const double neg  = -(vinner(world, xf, Vfull.x_alpha) +
                              vinner(world, yf, Vfull.y_alpha));
        const double swap =  (vinner(world, yf, Vfull.x_alpha) +
                              vinner(world, xf, Vfull.y_alpha));
        if (world.rank() == 0) {
          printf("  production  -(<xA|Vx>+<yA|Vy>) = %+.8e\n", neg);
          printf("  swapped      (<yA|Vx>+<xA|Vy>) = %+.8e\n", swap);
          printf("  (equal only if the source/contraction conventions match; "
                 "difference = %+.3e)\n", neg - swap);
        }
      }

      if (world.rank() == 0)
        print("\nVBC_PQ_TERMS DONE (informational — no gate; the equality "
              "gates live in test_vbc_spec_equivalence / test_tpa_pq_vs_vbc)");
    }
    world.gop.fence();
  } catch (const std::exception &e) {
    if (world.rank() == 0) print("EXCEPTION:", e.what());
    rc = 1;
  }
  finalize();
  return rc;
}
