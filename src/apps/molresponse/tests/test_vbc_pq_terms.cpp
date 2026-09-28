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
//     V.[M]  [+Σ_k x^C F^B_kp]     <-> P.B_mat  (same matrix, built as the
//                                                transposed daggered block)
//     V.[A]  [-Q̂ F^B x^C]          <-> P.B_app  (P omits Q̂ — legal under a
//                                                Q̂-projected contraction)
//     V.[L]  [-Q̂ g'[γ_L]φ]         <-> P.F_1..4 (leg by leg)
//     (the (C,B) half of V)        <-> P.D_app + P.D_mat
// HISTORY: until 2026-09-09 kernels/vbc.hpp wrote its exchange legs in
// reading order, which under madness::Exchange's convention builds every
// response density transposed; the "dagger" this test was written to hunt
// was in V, not in (P,Q) (reports/2026-09-09_orientation_derivation). With
// vbc.hpp in the equation's orientation the aligned slots now agree, and the
// VERDICT line below reads "REPRODUCES". The rows are kept as the term-level
// audit trail; the equality gates live in test_vbc_spec_equivalence /
// test_tpa_pq_vs_vbc.
// The 1e family (v^B / v^C content) is evaluated separately.
//
// Also probes the beta contraction convention: b1 = -(<xA|Vx> + <yA|Vy>)
// (production, kernels/beta.hpp) against the x<->y-swapped pairing
// <yA|Vx> + <xA|Vy> — printing both keeps the convention measurable.
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
#include "../kernels/tpa.hpp"
#include "../kernels/tpa_source_spec.hpp"
#include "../kernels/vbc.hpp"
#include "../solvers/build_response_ground_state.hpp"
#include "../solvers/es_save_load.hpp"
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

// Full contraction grid: every term against BOTH eigenvector halves.
// Columns: |V| <x|V> <y|V> | |P| <x|P> <y|P>  (V columns blank for P-only rows).
void row2(World &world, const char *label, const vecfuncT *v, const vecfuncT *p,
          const vecfuncT &fx, const vecfuncT &fy) {
  if (world.rank() != 0 && v && p) { /* inners are collective */ }
  double vn=0, vx=0, vy=0, pn=0, px=0, py=0;
  if (v) { vn=vnorm(world,const_cast<vecfuncT&>(*v)); vx=vinner(world,fx,*v); vy=vinner(world,fy,*v); }
  if (p) { pn=vnorm(world,const_cast<vecfuncT&>(*p)); px=vinner(world,fx,*p); py=vinner(world,fy,*p); }
  if (world.rank() == 0) {
    if (v && p)
      printf("  %-30s |V|=%10.3e <x|V>=%+12.5e <y|V>=%+12.5e | |P|=%10.3e <x|P>=%+12.5e <y|P>=%+12.5e\n",
             label, vn, vx, vy, pn, px, py);
    else if (p)
      printf("  %-30s %-56s | |P|=%10.3e <x|P>=%+12.5e <y|P>=%+12.5e\n",
             label, "", pn, px, py);
  }
}
void row(World &world, const char *label, const vecfuncT *v, const vecfuncT *p,
         const vecfuncT &f) {
  // kept for the leg-resolved section; single-contraction variant
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
      // --es-root=N (+ --calc-dir): use the CONVERGED eigenvector of root N as
      // (x^f,y^f), and default --freq to omega_N/2 (the physical leg freq).
      double freq_from_root = -1.0;
      if (parser.key_exists("calc-dir") && parser.key_exists("es-root")) {
        const std::string cdir = parser.value_raw("calc-dir");
        const int nroot = std::stoi(parser.value("es-root"));
        auto es = try_load_es_bundle<Full, ClosedShell>(world, cdir);
        if (es && nroot < static_cast<int>(es->state.roots.size())) {
          xf = madness::copy(world, es->state.roots[nroot].x_alpha);
          yf = madness::copy(world, es->state.roots[nroot].y_alpha);
          freq_from_root = 0.5 * es->state.omega(nroot);
          if (world.rank() == 0)
            print("[TERMS] eigenvector: CONVERGED root", nroot,
                  " omega=", es->state.omega(nroot),
                  " -> leg freq", freq_from_root,
                  " (", es->source_protocol_key, ")");
        } else if (world.rank() == 0) {
          print("  !! --es-root requested but no usable ES bundle — stand-in f");
        }
      }
      if (parser.key_exists("calc-dir") &&
          (parser.key_exists("freq") || freq_from_root > 0.0)) {
        const std::string cdir = parser.value_raw("calc-dir");
        const double freq = parser.key_exists("freq")
                                ? std::stod(parser.value("freq"))
                                : freq_from_root;
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
          char fbuf[64];
          std::snprintf(fbuf, sizeof fbuf, "%.5f", freq);
          src = std::string("STORED FD responses @ omega=") + fbuf;
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
        row2(world, "property-Fock MATRIX", &v_fphi, &p_bmat, xf, yf);

        auto v_fb    = vslot(0, 1);              // -Q̂ F^B x^C (+ image)
        auto p_bapp  = pslot(0, {1});            // -g'[γ^B†] x^C (+ image)
        row2(world, "property-Fock APPLY", &v_fb, &p_bapp, xf, yf);

        auto v_gz    = vslot(0, 0);              // -Q̂ g'[γ_ζ] φ (+ image)
        auto p_F     = pslot(0, {2, 3, 4, 5});   // -g'[D^{BC}] φ (+ image)
        row2(world, "pair density", &v_gz, &p_F, xf, yf);

        auto p_Dapp  = pslot(0, {6});            // -g'[γ^C†] x^B (+ image)
        auto p_Dmat  = pslot(0, {7});            // +Σ x^B G^C_kp (+ image)
        row2(world, "R family APPLY", nullptr, &p_Dapp, xf, yf);
        row2(world, "R family MATRIX", nullptr, &p_Dmat, xf, yf);
        auto p_R = pslot(0, {6, 7});
        row2(world, "R family TOTAL", nullptr, &p_R, xf, yf);

        // totals + engine-sum sanity
        auto v_tot = vsum(world, {v_gz, v_fb, v_fphi}, {0, 1, 2});
        auto p_tot = vsum(world, {p_bmat, p_bapp, p_F, p_R}, {0, 1, 2, 3});
        row2(world, "TOTAL (V.x | P)", &v_tot, &p_tot, xf, yf);
      }

      if (world.rank() == 0)
        print("\n[TERMS] ===== Q-channel vs V^BC Y-channel =====");
      {
        auto v_fphi = vslot(1, 2);
        auto q_bmat = pslot(1, {1});             // Q entry order: [0]=B_app
        row2(world, "property-Fock MATRIX", &v_fphi, &q_bmat, xf, yf);
        auto v_fb   = vslot(1, 1);
        auto q_bapp = pslot(1, {0});
        row2(world, "property-Fock APPLY", &v_fb, &q_bapp, xf, yf);
        auto v_gz   = vslot(1, 0);
        auto q_F    = pslot(1, {2, 3, 4, 5});
        row2(world, "pair density", &v_gz, &q_F, xf, yf);
        auto q_R = pslot(1, {6, 7});
        row2(world, "R family TOTAL", nullptr, &q_R, xf, yf);
        auto v_tot = vsum(world, {v_gz, v_fb, v_fphi}, {0, 1, 2});
        auto q_tot = vsum(world, {q_bmat, q_bapp, q_F, q_R}, {0, 1, 2, 3});
        row2(world, "TOTAL (V.y | Q)", &v_tot, &q_tot, xf, yf);
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
          row2(world, "1e: -Qv^B x^C", nullptr, &PT[8], xf, yf);
          row2(world, "1e: +Sum x^C <p|v^B|p>", nullptr, &PT[9], xf, yf);
          row2(world, "1e: -Qv^C x^B", nullptr, &PT[10], xf, yf);
          row2(world, "1e: +Sum x^B <p|v^C|p>", nullptr, &PT[11], xf, yf);
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

      // ============ V^BC pair density: LEG-RESOLVED + the relabeling identity
      // (2026-09-05, user request): expand V's gzeta into its individual legs
      // with norms, and verify numerically that the "opposite orientation"
      // relaxation block of D is the SAME kernel:
      //     sum_i phi_i(r) zbar_i(r')  ==  sum_i zeta_i(r) phi_i(r')
      // (index relabeling; zbar_i = sum_j phi_j <y_j|x_i>, zeta_i =
      //  sum_j phi_j <y_i|x_j>), hence D == gamma_L^dagger EXACTLY.
      if (world.rank() == 0)
        print("\n[TERMS] ===== V^BC pair density, leg-resolved + relabeling identity =====");
      {
        using source_spec::apply_entry;
        auto dot2 = [&](const vecfuncT &a, const vecfuncT &b) {
          auto r = common_ops::dot(world, a, b); r.scale(2.0); r.truncate();
          return r;
        };
        // zbar_i = sum_j phi_j <yb_j|xc_i>   (the Parker orientation)
        // zeta_i = sum_j phi_j <yb_i|xc_j>   (the VBC orientation)
        auto zbar = tpa::pq_detail::zblk(world, phi, yb, xc);
        auto zeta = tpa::pq_detail::zblk(world, phi, xc, yb);
        // V's gzeta legs, one entry each (J splits linearly with per-leg rho):
        source_spec::SourceSpec L1, L2a, L2b;
        L1.entries.push_back(apply_entry(dot2(xb, yc), {{xb, yc}}, phi, -1.0, true));
        L2a.entries.push_back(apply_entry(dot2(phi, zeta), {{phi, zeta}}, phi, -1.0, true));
        L2b.entries.push_back(apply_entry(dot2(zbar, phi), {{zbar, phi}}, phi, -1.0, true));
        auto l1  = source_spec::assemble_source(world, g0, {L1})[0];
        auto l2a = source_spec::assemble_source(world, g0, {L2a})[0];
        auto l2b = source_spec::assemble_source(world, g0, {L2b})[0];
        if (world.rank() == 0) {
          printf("  V gzeta leg (x^B,y^C):            |leg|=%11.4e  <f|leg>=%+13.6e\n",
                 vnorm(world, l1), vinner(world, xf, l1));
          printf("  V gzeta leg (phi,zeta)  [VBC or.] |leg|=%11.4e  <f|leg>=%+13.6e\n",
                 vnorm(world, l2a), vinner(world, xf, l2a));
          printf("  D leg      (zbar,phi) [Parker or.]|leg|=%11.4e  <f|leg>=%+13.6e\n",
                 vnorm(world, l2b), vinner(world, xf, l2b));
          auto d = vdiff(world, l2a, l2b);
          printf("  RELABELING IDENTITY  ||(phi,zeta)-(zbar,phi)|| = %.3e  (0 => D == gamma_L^T, no new math)\n",
                 vnorm(world, d));
        }
      }

      // ============ BUNDLING CHECK (2026-09-07): is the ~6x pair-density
      // norm ratio purely an accounting difference?  Compare, for ONE
      // ordering: (a) the four unbundled family-F entries summed, (b) the
      // same physics as ONE bundled entry (untruncated -> must EQUAL (a)),
      // (c) bundled + truncated like vbc's gzeta -> the apples-to-apples
      // partner for |V_gzeta|.
      if (world.rank() == 0)
        print("\n[TERMS] ===== pair-density bundling check =====");
      {
        auto F_unb = vsum(world, Pbc[0], {2, 3, 4, 5});
        source_spec::SourceSpec sb, sbt;
        sb.entries.push_back(tpa::pq_family_F_bundled(world, g0, B, C, false));
        sbt.entries.push_back(tpa::pq_family_F_bundled(world, g0, B, C, true));
        auto F_bun  = source_spec::assemble_source(world, g0, {sb})[0];
        auto F_bunt = source_spec::assemble_source(world, g0, {sbt})[0];
        auto v_gz1  = source_spec::assemble_source(
            world, g0,
            vbc::vbc_half_spec(world, g0, xb, yb, xc, yc, zeta_bc, zop))[0];
        auto d = vdiff(world, F_bun, F_unb);
        if (world.rank() == 0) {
          printf("  (a) 4 unbundled entries summed : |F|=%11.4e  <x|F>=%+13.6e\n",
                 vnorm(world, F_unb), vinner(world, xf, F_unb));
          printf("  (b) 1 bundled entry (untrunc.) : |F|=%11.4e  <x|F>=%+13.6e\n",
                 vnorm(world, F_bun), vinner(world, xf, F_bun));
          printf("      ||(b)-(a)|| = %.3e   <== 0 proves bundling is EXACT\n",
                 vnorm(world, d));
          printf("  (c) 1 bundled entry (truncated): |F|=%11.4e  <x|F>=%+13.6e\n",
                 vnorm(world, F_bunt), vinner(world, xf, F_bunt));
          printf("  NOTE |V.x whole-channel| = %11.4e (gzeta+fb+fphi, one ordering)\n",
                 vnorm(world, v_gz1));
        }
      }

      // ============ DENSITY-LEVEL probe (2026-09-07): the residual pair-
      // density ratio. Bundling was proven exact and is norm-NEUTRAL, so it
      // cannot explain the ~6x. Decompose what's left: compare the DENSITIES
      // themselves (rho_D vs rho_gammaL) and the applied results with the
      // exchange switched off (J-only), which separates Coulomb from exchange.
      if (world.rank() == 0)
        print("\n[TERMS] ===== density-level probe =====");
      {
        auto dot2t = [&](const vecfuncT &a, const vecfuncT &b) {
          auto r = common_ops::dot(world, a, b); r.scale(2.0); r.truncate(); return r;
        };
        // V's FULL gamma_L density (both halves, as the two calls sum to)
        auto rho_gL = common_ops::dot(world, xb, yc);
        rho_gL += common_ops::dot(world, xc, yb);
        rho_gL += common_ops::dot(world, phi, zeta_bc);
        rho_gL += common_ops::dot(world, phi, zeta_cb);
        rho_gL.scale(2.0); rho_gL.truncate();
        // P's FULL D density (one call already carries both leg orders)
        auto z1 = tpa::pq_detail::zblk(world, phi, yb, xc);
        auto z2 = tpa::pq_detail::zblk(world, phi, xb, yc);
        auto rho_D = common_ops::dot(world, xc, yb);
        rho_D += common_ops::dot(world, xb, yc);
        rho_D -= common_ops::dot(world, z1, phi);
        rho_D -= common_ops::dot(world, phi, z2);
        rho_D.scale(2.0); rho_D.truncate();
        auto dd = rho_D - rho_gL;
        if (world.rank() == 0) {
          printf("  ||rho_gammaL|| = %11.4e   trace = %+13.6e\n",
                 rho_gL.norm2(), rho_gL.trace());
          printf("  ||rho_D||      = %11.4e   trace = %+13.6e\n",
                 rho_D.norm2(), rho_D.trace());
          printf("  ||rho_D - rho_gammaL|| = %.4e  <== 0 if the two densities agree\n",
                 dd.norm2());
        }
        // number of exchange legs actually applied in each entry family
        if (world.rank() == 0)
          printf("  legs: V gzeta per call = 2 (x2 calls = 4);  P family F per call = 4\n"
                 "        (P's family F is ordering-COMPLETE in one call -> summing both\n"
                 "         calls in the aligned table double-counts it: factor 2 of the ratio)\n");
      }

      // ============ UNIFICATION CANDIDATE: P == V^BC built at swapped legs
      // (2026-09-05, user question): can 2PA be computed with the V^BC
      // builder itself? The conjugation + D==gamma_L^T theorems say YES with
      // one amendment: swap the x/y halves of the PHOTON LEGS fed to the
      // builder (legs at -omega), keep the natural state pairing. This row
      // measures || vbc_half_spec(swapped legs) - P^{BC}(2e) || — if ~0 (up
      // to the Q-projection difference on the fb slot, which is null against
      // a Q-projected state), tpa_pq_spec can be RETIRED in favor of the one
      // V^BC builder for all three quadratic properties.
      if (world.rank() == 0)
        print("\n[TERMS] ===== unification: V^BC(swapped legs) vs P(2e) =====");
      {
        real_function_3d zop2 = madness::copy(phi[0]); zop2.scale(0.0);
        // swap (x,y)->(y,x) on both legs; V's zeta_bc = make_zeta(y_B, x_C)
        // becomes make_zeta(x_B, y_C) under the swap.
        auto zeta_sw = vbc::make_zeta(world, xb, yc, phi);
        auto Vsw = source_spec::assemble_source(
            world, g0,
            vbc::vbc_half_spec(world, g0, yb, xb, yc, xc, zeta_sw, zop2));
        auto p2e = vsum(world, Pbc[0], {0, 1, 2, 3, 4, 5});   // famB+famF, one ordering
        auto d = vdiff(world, Vsw[0], p2e);
        if (world.rank() == 0)
          printf("  ||V^BC.x(swapped) - P(2e)|| = %.3e   (|Vsw|=%.3e |P2e|=%.3e)\n",
                 vnorm(world, d), vnorm(world, Vsw[0]), vnorm(world, p2e));
        auto q2e = vsum(world, Pbc[1], {0, 1, 2, 3, 4, 5});
        auto dy = vdiff(world, Vsw[1], q2e);
        if (world.rank() == 0)
          printf("  ||V^BC.y(swapped) - Q(2e)|| = %.3e   (|Vsw|=%.3e |Q2e|=%.3e)\n",
                 vnorm(world, dy), vnorm(world, Vsw[1]), vnorm(world, q2e));
      }

      // ===== FLAG-TOGGLE: is (P,Q) the ADJOINT of V^BC? ====================
      // 2e-only, ONE ordering (B,C). The test builds V with a ZERO one-electron
      // operator, so v^B is already neutralized and the ONLY differences
      // between V's {fb,fphi} and P's family B are three flags:
      //   APPLY : exchange bra/ket order, and project_Q
      //   MATRIX: transpose_matrix
      // Toggle them one at a time and see which combination maps V onto P
      // EXACTLY (norm of the difference, not just a contraction).
      if (world.rank() == 0)
        print("\n[TERMS] ===== FLAG TOGGLE: does V^BC daggered == P? =====");
      {
        using source_spec::apply_entry;
        using source_spec::occupied_matrix_entry;
        // rho_B exactly as BOTH builders make it (they agree by construction)
        real_function_3d rho_B = common_ops::dot(world, xb, phi);
        rho_B += common_ops::dot(world, phi, yb);
        rho_B.scale(2.0); rho_B.truncate();

        auto one = [&](source_spec::SourceEntry e) {
          source_spec::SourceSpec S; S.entries.push_back(std::move(e));
          source_spec::SourceSpec E;
          return source_spec::assemble_source(world, g0, {S, E})[0];
        };
        auto rpt = [&](const char *tag, const vecfuncT &a, const vecfuncT &ref_) {
          vecfuncT d = sub(world, a, ref_);
          if (world.rank() == 0)
            printf("   %-34s |.|=%10.4e  <x|.>=%+12.5e  <y|.>=%+12.5e  ||.-P||=%9.3e\n",
                   tag, vnorm(world, const_cast<vecfuncT&>(a)),
                   vinner(world, xf, a), vinner(world, yf, a), vnorm(world, d));
        };
        // ---- APPLY: P's family-B reference
        vecfuncT P_app = one(apply_entry(rho_B, {{phi, xb}, {yb, phi}},
                                         xc, -1.0, /*Q=*/false));
        if (world.rank() == 0) print("  -- APPLY family (P ref = fam-B apply) --");
        rpt("P famB apply  [swap,Q=0] (REF)", P_app, P_app);
        rpt("V fb as-built [orig,Q=1]",
            one(apply_entry(rho_B, {{xb, phi}, {phi, yb}}, xc, -1.0, true)), P_app);
        rpt("V fb  [orig,Q=0]",
            one(apply_entry(rho_B, {{xb, phi}, {phi, yb}}, xc, -1.0, false)), P_app);
        rpt("V fb  [SWAP,Q=1]",
            one(apply_entry(rho_B, {{phi, xb}, {yb, phi}}, xc, -1.0, true)), P_app);
        rpt("V fb  [SWAP,Q=0]",
            one(apply_entry(rho_B, {{phi, xb}, {yb, phi}}, xc, -1.0, false)), P_app);
        // ---- MATRIX: P's family-B reference
        vecfuncT P_mat = one(occupied_matrix_entry(rho_B, {{xb, phi}, {phi, yb}},
                                                   phi, xc, +1.0, /*T=*/true));
        if (world.rank() == 0) print("  -- MATRIX family (P ref = fam-B matrix) --");
        rpt("P famB matrix [T=1] (REF)", P_mat, P_mat);
        rpt("V fphi as-built [T=0]",
            one(occupied_matrix_entry(rho_B, {{xb, phi}, {phi, yb}},
                                      phi, xc, +1.0, false)), P_mat);
        rpt("V fphi [T=1]",
            one(occupied_matrix_entry(rho_B, {{xb, phi}, {phi, yb}},
                                      phi, xc, +1.0, true)), P_mat);
      }

      // ===== HERMITIAN LIMIT: the falsifiable prediction ===================
      // F^B = v^B + g'[gamma^B], gamma^B = x^B phi^T + phi y^B^T, is NOT a
      // symmetric operator when x^B != y^B, so the exchange bra/ket
      // orientation is physical, not conventional. V^BC uses F^B (correct for
      // the SOLVE: F^(1) acting right on phi^(1)); (P,Q) uses g'[gamma^Bdag].
      // Those differ ONLY through the non-Hermiticity.
      // PREDICTION: set y := x on both legs (gamma Hermitian, F symmetric,
      // adjoint = no-op) => the two sources must become IDENTICAL.
      if (world.rank() == 0)
        print("\n[TERMS] ===== HERMITIAN LIMIT (y:=x): must make V == P =====");
      {
        ResponseStateXY<ClosedShell> Bh, Ch;
        Bh.x_alpha = madness::copy(world, xb); Bh.y_alpha = madness::copy(world, xb);
        Ch.x_alpha = madness::copy(world, xc); Ch.y_alpha = madness::copy(world, xc);
        auto zeta_h = vbc::make_zeta(world, Bh.y_alpha, Ch.x_alpha, phi);
        auto Vh = source_spec::assemble_source(
            world, g0, vbc::vbc_half_spec(world, g0, Bh.x_alpha, Bh.y_alpha,
                                          Ch.x_alpha, Ch.y_alpha, zeta_h, zop));
        auto zeta_h_cb = vbc::make_zeta(world, Ch.y_alpha, Bh.x_alpha, phi);
        auto Vh2 = source_spec::assemble_source(
            world, g0, vbc::vbc_half_spec(world, g0, Ch.x_alpha, Ch.y_alpha,
                                          Bh.x_alpha, Bh.y_alpha, zeta_h_cb, zop));
        auto Ph = tpa::assemble_tpa_pq(world, g0, Bh, Ch);
        vecfuncT Vfull_x = add(world, Vh[0], Vh2[0]);
        vecfuncT Vfull_y = add(world, Vh[1], Vh2[1]);
        if (world.rank() == 0) {
          printf("   ||V_half.x||       = %10.4e   <x|.>=%+12.5e\n",
                 vnorm(world, Vh[0]), vinner(world, xf, Vh[0]));
          printf("   ||V_full.x||       = %10.4e   <x|.>=%+12.5e\n",
                 vnorm(world, Vfull_x), vinner(world, xf, Vfull_x));
          printf("   ||P.x||            = %10.4e   <x|.>=%+12.5e\n",
                 vnorm(world, Ph.x_alpha), vinner(world, xf, Ph.x_alpha));
          printf("   ||V_full.x - P.x|| = %10.4e   ||V_full.y - Q.y|| = %10.4e\n",
                 vnorm(world, sub(world, Vfull_x, Ph.x_alpha)),
                 vnorm(world, sub(world, Vfull_y, Ph.y_alpha)));
          printf("   ||V_half.x - P.x|| = %10.4e\n",
                 vnorm(world, sub(world, Vh[0], Ph.x_alpha)));
        }
      }

      // ============ THE CLEAN TEST (2026-09-08): does the DRIVEN source
      // V^BC, contracted with a TRUE eigenvector, equal the residue answer?
      // The first-principles 2n+1 elimination, read literally, says the
      // residue hands the eigenvector the SAME V^BC that beta uses. Our code
      // instead uses the daggered (P,Q). Gate 1 already pins
      // <x^f|P>+<y^f|Q> == e3/sqrt2 (one ordering), so the NEW information is
      // whether the V^BC contraction also equals e3/sqrt2. Requires a REAL
      // eigenvector: with stand-ins the residue identity need not hold,
      // because it leans on the homogeneous equation (F0-eps-w_f)x^f =
      // -Q g'[gamma^f]phi.
      if (world.rank() == 0)
        print("\n[TERMS] ===== CLEAN TEST: V^BC vs (P,Q) at the eigenvector =====");
      {
        // reference: validated c-grouped 2e answer, ONE ordering (B,C)
        ResponseStateXY<ClosedShell> F;   // the eigenvector as a state
        F.x_alpha = madness::copy(world, xf); F.y_alpha = madness::copy(world, yf);
        const double e3   = tpa::tpa_e3_residue(world, g0, B, C, F);
        const double ref  = e3 / std::sqrt(2.0);
        // (P,Q), one ordering
        auto PQ1 = tpa::assemble_tpa_pq(world, g0, B, C);
        const double s_pq = vinner(world, xf, PQ1.x_alpha) + vinner(world, yf, PQ1.y_alpha);
        // V^BC (2e only), ONE ordering — vbc_half_spec(B,C) is the (B,C) half
        auto vhalf = source_spec::assemble_source(
            world, g0, vbc::vbc_half_spec(world, g0, xb, yb, xc, yc, zeta_bc, zop));
        const double s_v  =  vinner(world, xf, vhalf[0]) + vinner(world, yf, vhalf[1]);
        const double s_vn = -s_v;
        if (world.rank() == 0) {
          printf("  reference  e3/sqrt2 (c-grouped, validated) = %+.8e\n", ref);
          printf("  (P,Q)      <x|P>+<y|Q>                     = %+.8e   dev = %+.3e\n",
                 s_pq, s_pq - ref);
          printf("  V^BC       <x|Vx>+<y|Vy>                   = %+.8e   dev = %+.3e\n",
                 s_v, s_v - ref);
          printf("  V^BC neg  -(<x|Vx>+<y|Vy>)                 = %+.8e   dev = %+.3e\n",
                 s_vn, s_vn - ref);
          // ADJOINT HYPOTHESIS: the 2n+1 step moves the resolvent left as its
          // ADJOINT. The response pencil is non-self-adjoint (paired Delta
          // metric), so its LEFT eigenvector is the x<->y swap of the right
          // one. Prediction: contracting the SWAPPED eigenvector against the
          // driven V^BC must reproduce the reference — equivalently, (P,Q) is
          // the adjoint image of V^BC.
          const double s_vsw = vinner(world, yf, vhalf[0]) + vinner(world, xf, vhalf[1]);
          printf("  V^BC swap  <y|Vx>+<x|Vy>  (LEFT eigvec)     = %+.8e   dev = %+.3e\n",
                 s_vsw, s_vsw - ref);
          printf("  V^BC swap neg                              = %+.8e   dev = %+.3e\n",
                 -s_vsw, -s_vsw - ref);
          printf("  ratio ref/V^BC = %.6f   ratio ref/V^BCswap = %.6f\n",
                 ref/s_v, ref/s_vsw);
          printf("  VERDICT: %s\n",
                 (std::abs(s_v-ref) < 1e-6*std::max(1.0,std::abs(ref)) ||
                  std::abs(s_vn-ref) < 1e-6*std::max(1.0,std::abs(ref)))
                   ? "V^BC REPRODUCES the residue answer (one source; expected since 2026-09-09)"
                   : "V^BC does NOT reproduce it -> orientation regression in vbc.hpp");
        }
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

        // ===== THE SCALAR QUESTION (user, 2026-09-05): with plain V^BC
        // (normal legs, no daggers), does ANY pairing of the state halves
        // against the two V channels reproduce <x|P>+<y|Q>?  All four inner
        // products printed, all candidate pairings compared, both sides
        // symmetrized over orderings (Vfull is (B,C)+(C,B); P/Q summed the
        // same way). Density-level theory predicts NO pairing matches at
        // c_x=1 (the dagger + transposition are real); this measures it.
        auto P_tot = vsum(world, Pbc[0], {0,1,2,3,4,5,6,7});
        {
          auto pcb = vsum(world, Pcb[0], {0,1,2,3,4,5,6,7});
          gaxpy(world, 1.0, P_tot, 1.0, pcb);
        }
        auto Q_tot = vsum(world, Pbc[1], {0,1,2,3,4,5,6,7});
        {
          auto qcb = vsum(world, Pcb[1], {0,1,2,3,4,5,6,7});
          gaxpy(world, 1.0, Q_tot, 1.0, qcb);
        }
        const double xVx = vinner(world, xf, Vfull.x_alpha);
        const double yVy = vinner(world, yf, Vfull.y_alpha);
        const double yVx = vinner(world, yf, Vfull.x_alpha);
        const double xVy = vinner(world, xf, Vfull.y_alpha);
        const double S_ours = vinner(world, xf, P_tot) + vinner(world, yf, Q_tot);
        if (world.rank() == 0) {
          print("\n[TERMS] ===== scalar pairings: plain V^BC vs <x|P>+<y|Q> =====");
          printf("  <x|Vx> = %+.8e   <y|Vy> = %+.8e\n", xVx, yVy);
          printf("  <y|Vx> = %+.8e   <x|Vy> = %+.8e\n", yVx, xVy);
          printf("  S_ours  = <x|P>+<y|Q>          = %+.8e   (both orderings)\n", S_ours);
          printf("  S_a     =  <y|Vx>+<x|Vy>       = %+.8e   dev = %+.3e\n",  yVx + xVy,  (yVx + xVy) - S_ours);
          printf("  S_b     = -(<x|Vx>+<y|Vy>)     = %+.8e   dev = %+.3e\n", -(xVx + yVy), -(xVx + yVy) - S_ours);
          printf("  S_c     = -(<y|Vx>+<x|Vy>)     = %+.8e   dev = %+.3e\n", -(yVx + xVy), -(yVx + xVy) - S_ours);
          printf("  S_d     =  <x|Vx>+<y|Vy>       = %+.8e   dev = %+.3e\n",  (xVx + yVy),  (xVx + yVy) - S_ours);
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
