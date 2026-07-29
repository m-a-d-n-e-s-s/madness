// ===========================================================================
// test_kernel_equivalence.cpp — pins the unified two-electron kernel against
// the Dalton-validated linear kernel. Asserts, for a closed-shell Full state χ:
//   Kernels<Full,ClosedShell>::apply_g(χ, Φ, Φ)  ==  Kernels<Full,ClosedShell>
//                                                        ::compute_gamma(χ, ρ)
// where Φ is the ground state dressed as {φ, φ}. Both apply Qa + truncate, so a
// match to ~thresh proves the unification (exchange pairing + projection +
// density factor) before any solver is switched over. ALLOCATION test (real MRA).
//
//   test_kernel_equivalence --archive=<moldft restartdata> [--thresh=1e-6] [--k=N]
// ===========================================================================

#include "../GroundState.hpp"
#include "../Perturbations.hpp"
#include "../ResponseProtocol.hpp"
#include "../kernels/full.hpp"
#include "../kernels/tags.hpp"
#include "../kernels/two_electron.hpp"
#include "../solvers/build_response_ground_state.hpp"
#include "../solvers/response_state.hpp"

#include <madness/misc/info.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

#include <cmath>
#include <filesystem>
#include <fstream>
#include <string>

using namespace madness;
using namespace molresponse_v3;

int main(int argc, char **argv) {
  World &world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);
    {
      commandlineparser parser(argc, argv);
      if (!parser.key_exists("archive")) {
        if (world.rank() == 0)
          print("Usage: test_kernel_equivalence --archive=<path> [--thresh=X] [--k=N]");
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
          std::ifstream ifs(cand); nlohmann::json j; ifs >> j;
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

      // χ: a closed-shell Full state with x != y (exercises cross-channel
      // exchange). Use distinct dipole perturbations for x and y.
      ResponseStateXY<ClosedShell> chi;
      chi.x_alpha = dipole_perturbation(world, gs, 0);  // d/dx
      chi.y_alpha = dipole_perturbation(world, gs, 1);  // d/dy
      chi.x_alpha = g0.Qa(chi.x_alpha);
      chi.y_alpha = g0.Qa(chi.y_alpha);

      // Φ: the ground state dressed as {φ, φ}.
      ResponseStateXY<ClosedShell> Phi;
      Phi.x_alpha = madness::copy(world, g0.amo);
      Phi.y_alpha = madness::copy(world, g0.amo);

      using KF = Kernels<Full, ClosedShell>;
      auto rho_ref   = KF::compute_density(world, g0, chi);
      auto gamma_ref = KF::compute_gamma(world, g0, chi, rho_ref);

      // Unified path: Kernels<Full,ClosedShell>::apply_g rebuilds rho from
      // (chi, Phi) via the two-state compute_density and returns it.
      auto [gamma_new, rho_new] = KF::apply_g(world, g0, chi, Phi, Phi);

      const double dr  = (rho_ref - rho_new).norm2();
      double dx = 0.0, dy = 0.0;
      for (size_t i = 0; i < gamma_ref.x_alpha.size(); ++i)
        dx += (gamma_ref.x_alpha[i] - gamma_new.x_alpha[i]).norm2();
      for (size_t i = 0; i < gamma_ref.y_alpha.size(); ++i)
        dy += (gamma_ref.y_alpha[i] - gamma_new.y_alpha[i]).norm2();

      const double tol = 10.0 * t;
      if (world.rank() == 0) {
        print("\n=== kernel equivalence (Full, ClosedShell) ===");
        print("  ||rho_ref - rho_new||      =", dr);
        print("  Sum ||gamma_x diff||       =", dx);
        print("  Sum ||gamma_y diff||       =", dy);
        print("  tol                        =", tol);
        bool ok = dr < tol && dx < tol && dy < tol;
        print("\n", ok ? "PASSED" : "FAILED",
              " (unified compute_g == compute_gamma<Full,ClosedShell>)");
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
