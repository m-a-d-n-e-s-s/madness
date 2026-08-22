// ===========================================================================
// test_tpa_regroup_identities.cpp — numerically pins the three algebraic moves
// that underlie the state-grouped reformulation of the 2PA two-electron
// residue (paper Sec. "A State-Grouped Reformulation").
//
// The c-grouped implementation in kernels/tpa.hpp contracts E[3] on the c leg;
// Parker's form contracts on the state leg, <X_f|(p,q)>. Regrouping between
// them needs exactly three identities. This test asserts each one directly on
// real MRA functions, so the derivation rests on measurement rather than on
// symbol pushing:
//
//   M1  <w|K[a,b]f>            == <f|K[b,a]w>            (operator transpose)
//   M2  <z^{uv}|w>             == <u|T[v;M]>             (state out of a z block)
//   M3a <x^c|K[x^f,phi]y^b>    == <x^f|K[x^c,y^b]phi>    (state out of the operator)
//   M3b <x^c|K[phi,y^f]y^b>    == <y^f|K[y^b,x^c]phi>
//   M3c <x^c|J[rho^f]y^b>      == 2 sum_j <x^f_j+y^f_j| phi_j J[rho^cb]>
//
// M3 is the load-bearing move: it is what lets T1/T3 be regrouped at all, since
// there the state sits inside the density and the exchange pairs rather than in
// a contracted vector. Its right-hand side is a SINGLE exchange application
// built from the photon legs acting on phi -- the reason the reformulation is
// cheaper, not merely tidier.
//
// The six vectors are stand-ins, not physical responses: the identities are
// algebraic and must hold for arbitrary functions, so distinct combinations of
// dipole perturbations exercise them. ALLOCATION test (real MRA).
//
//   test_tpa_regroup_identities --archive=<moldft restartdata> [--thresh=X] [--k=N]
// ===========================================================================

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/common_ops.hpp"
#include "../kernels/tags.hpp"
#include "../solvers/build_response_ground_state.hpp"
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
          print("Usage: test_tpa_regroup_identities --archive=<path> "
                "[--thresh=X] [--k=N]");
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
          else if (j.contains("molecule"))
            mj = j["molecule"];
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

      // Six distinct stand-in vectors. The identities are algebraic, so any
      // functions serve. Two requirements make the test informative rather than
      // merely passing:
      //   (a) Do NOT Q-project. Q-projected vectors are orthogonal to phi, so
      //       every <phi|...> factor vanishes identically and the M2 check --
      //       and any z-block term -- collapses to 0 == 0.
      //   (b) Admix phi with irregular coefficients, which breaks the molecular
      //       point group. With pure Cartesian dipoles on a symmetric molecule
      //       several contractions vanish by selection rule, again giving
      //       0 == 0 rather than a real comparison.
      auto d0 = dipole_perturbation(world, gs, 0);
      auto d1 = dipole_perturbation(world, gs, 1);
      auto d2 = dipole_perturbation(world, gs, 2);
      auto comb = [&](double a, const vecfuncT &u, double b, const vecfuncT &v,
                      double c) {
        vecfuncT r = madness::copy(world, u);
        gaxpy(world, a, r, b, v);
        gaxpy(world, 1.0, r, c, phi);  // phi admixture: breaks (a) and (b)
        truncate(world, r);
        return r;
      };
      const vecfuncT xb = comb(1.0, d0, 0.23, d1, 0.31);
      const vecfuncT yb = comb(1.0, d1, -0.41, d2, 0.17);
      const vecfuncT xc = comb(1.0, d2, 0.13, d0, -0.29);
      const vecfuncT yc = comb(1.0, d0, 0.53, d1, 0.11);
      const vecfuncT xf = comb(1.0, d1, -0.37, d2, -0.19);
      const vecfuncT yf = comb(1.0, d2, 0.71, d0, 0.43);

      auto vecdot = [&](const vecfuncT &a, const vecfuncT &b) {
        vecfuncT p(a.size());
        for (size_t i = 0; i < a.size(); ++i) p[i] = a[i] * b[i];
        truncate(world, p);
        return sum(world, p);
      };

      struct Check { std::string name; double lhs; double rhs; };
      std::vector<Check> checks;
      auto add = [&](const std::string &n, double l, double r) {
        Check c;
        c.name = n;
        c.lhs = l;
        c.rhs = r;
        checks.push_back(c);
      };

      // --- M1: <w|K[a,b]f> == <f|K[b,a]w> -------------------------------
      add("M1  <y^c|K[x^b,phi]x^f> == <x^f|K[phi,x^b]y^c>",
          inner(yc, apply_exchange(world, xb, phi, xf, lo)),
          inner(xf, apply_exchange(world, phi, xb, yc, lo)));

      // --- M2: <z^{uv}|w> == <u|T[v;M]>, M_ji = <phi_j|w_i> --------------
      // z^{uv}_i = sum_j phi_j <u_j|v_i>; transform(w,M)_k = sum_j w_j M(j,k),
      // and for real functions <w_i|phi_j> == <phi_j|w_i>, so the matrix that
      // realizes T[v;M]_j = sum_i v_i <phi_j|w_i> is matrix_inner(w, phi).
      {
        const vecfuncT &u = xf, &v = xc, &w = yb;
        auto z = transform(world, phi, matrix_inner(world, u, v), true);
        auto T = transform(world, v, matrix_inner(world, w, phi), true);
        add("M2  <z^{x^f x^c}|y^b>    == <x^f|T[x^c;<phi|y^b>]>", inner(z, w),
            inner(u, T));
      }

      // --- M3a/M3b: pull the state out of the exchange operator ----------
      add("M3a <x^c|K[x^f,phi]y^b>  == <x^f|K[x^c,y^b]phi>",
          inner(xc, apply_exchange(world, xf, phi, yb, lo)),
          inner(xf, apply_exchange(world, xc, yb, phi, lo)));
      add("M3b <x^c|K[phi,y^f]y^b>  == <y^f|K[y^b,x^c]phi>",
          inner(xc, apply_exchange(world, phi, yf, yb, lo)),
          inner(yf, apply_exchange(world, yb, xc, phi, lo)));

      // --- M3c: the Coulomb analogue -------------------------------------
      // rho^f = 2 sum_j phi_j (x^f_j + y^f_j),  rho^cb = sum_i x^c_i y^b_i
      {
        auto xfyf = comb(1.0, xf, 1.0, yf, 0.0);  // exactly x^f + y^f
        auto rho_f = vecdot(phi, xfyf);
        rho_f.scale(2.0);
        auto rho_cb = vecdot(xc, yb);
        auto J_f = apply(*coulop, rho_f);
        auto J_cb = apply(*coulop, rho_cb);
        const double lhs = inner(xc, mul(world, J_f, yb, true));
        const double rhs = 2.0 * inner(xfyf, mul(world, J_cb, phi, true));
        add("M3c <x^c|J[rho^f]y^b>    == 2<x^f+y^f|phi J[rho^cb]>", lhs, rhs);
      }

      // Relative tolerance: these are exact algebraic identities, so the only
      // error is the MRA representation of the intermediates. Scale by the
      // magnitude of the quantity so a near-zero contraction is not judged on
      // an absolute floor it cannot meet.
      const double tol = 200.0 * t;
      bool ok = true;
      if (world.rank() == 0) {
        print("\n=== 2PA state-grouped regrouping identities ===");
        print("  thresh =", t, "  rel tol =", tol, "  n_occ =", phi.size());
        printf("\n  %-52s %14s %14s %11s\n", "identity", "lhs", "rhs", "rel.diff");
      }
      // A contraction that vanishes identically agrees with anything, so a
      // near-zero pair is not evidence. Treat it as a failure of the test
      // rather than a pass of the identity.
      const double informative = 1e-4;
      for (const auto &c : checks) {
        const double scale = std::max(std::abs(c.lhs), std::abs(c.rhs));
        const double rel = scale > 1e-12 ? std::abs(c.lhs - c.rhs) / scale
                                         : std::abs(c.lhs - c.rhs);
        const bool weak = scale < informative;
        const bool pass = rel < tol && !weak;
        ok = ok && pass;
        if (world.rank() == 0)
          printf("  %-52s %14.8f %14.8f %11.2e %s\n", c.name.c_str(), c.lhs,
                 c.rhs, rel,
                 weak ? "UNINFORMATIVE (|value| too small)"
                      : (pass ? "ok" : "FAIL"));
      }
      if (world.rank() == 0) {
        print("\n", ok ? "PASSED" : "FAILED",
              " (M1/M2/M3 hold => the state-grouped regrouping is sound)");
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
