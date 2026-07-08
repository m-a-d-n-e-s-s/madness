// =========================================================================
// test_state_archive_hdf5.cpp — io-hdf5 P2a: round-trip the real restart units
// ResponseStateX<ClosedShell> (X only; static / TDA-X) AND
// ResponseStateXY<ClosedShell> (X+Y blocks; dynamic-alpha / Full-ES) — each a
// *vector* of orbitals — through BOTH the legacy archive and the opt-in HDF5
// path, and confirm the loader AUTO-DETECTS the .h5. Exercises the wiring in
// response_state.hpp:
//   - env MADRESPONSE_IO_HDF5 set  => save() writes <file>.h5
//   - load() finds <file>.h5       => reads HDF5; else legacy
// Built only when MADNESS_ENABLE_HDF5=ON. NP=1 sanity (path is multi-rank-ready).
// =========================================================================

#include "../solvers/response_state.hpp"

#include <madness/mra/mra.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

static const std::size_t D = 3;

static double g0(const Vector<double, D>& r) {
  return std::exp(-1.0 * (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]));
}
static double g1(const Vector<double, D>& r) {
  const double x = r[0] - 0.5;
  return std::exp(-1.5 * (x * x + r[1] * r[1] + r[2] * r[2]));
}
static double g2(const Vector<double, D>& r) {
  const double y = r[1] + 0.4, z = r[2] - 0.3;
  return std::exp(-0.8 * (r[0] * r[0] + y * y + z * z));
}

static double state_max_err(const ResponseStateX<ClosedShell>& a,
                            const ResponseStateX<ClosedShell>& b) {
  if (a.x_alpha.size() != b.x_alpha.size()) return 1e9;
  double e = 0;
  for (std::size_t i = 0; i < a.x_alpha.size(); ++i)
    e = std::max(e, (a.x_alpha[i] - b.x_alpha[i]).norm2());  // collective
  return e;
}

// XY state carries both x_alpha and y_alpha (dynamic-alpha / Full-ES).
static double statexy_max_err(const ResponseStateXY<ClosedShell>& a,
                              const ResponseStateXY<ClosedShell>& b) {
  if (a.x_alpha.size() != b.x_alpha.size() ||
      a.y_alpha.size() != b.y_alpha.size())
    return 1e9;
  double e = 0;
  for (std::size_t i = 0; i < a.x_alpha.size(); ++i)
    e = std::max(e, (a.x_alpha[i] - b.x_alpha[i]).norm2());  // collective
  for (std::size_t i = 0; i < a.y_alpha.size(); ++i)
    e = std::max(e, (a.y_alpha[i] - b.y_alpha[i]).norm2());  // collective
  return e;
}

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  startup(world, argc, argv);
  std::cout.precision(8);

  // both variants share the same projection defaults
  FunctionDefaults<D>::set_k(8);
  FunctionDefaults<D>::set_thresh(1e-6);
  FunctionDefaults<D>::set_cubic_cell(-20.0, 20.0);

  bool ok_x = false, ok_xy = false;

  // ---- ResponseStateX<ClosedShell> — X only (static / TDA-X) ----
  {
    ResponseStateX<ClosedShell> st;
    st.x_alpha.push_back(FunctionFactory<double, D>(world).f(g0));
    st.x_alpha.push_back(FunctionFactory<double, D>(world).f(g1));
    st.x_alpha.push_back(FunctionFactory<double, D>(world).f(g2));
    for (auto& f : st.x_alpha) f.truncate();

    // legacy path: env OFF -> save writes <file>.00000; load (no .h5) -> legacy
    unsetenv("MADRESPONSE_IO_HDF5");
    st.save(world, "state_leg");
    auto s_leg = ResponseStateX<ClosedShell>::load(world, "state_leg");
    const double e_leg = state_max_err(st, s_leg);

    // HDF5 path: env ON -> save writes state_h5.h5; load auto-detects it
    setenv("MADRESPONSE_IO_HDF5", "1", 1);
    st.save(world, "state_h5");
    auto s_h5 = ResponseStateX<ClosedShell>::load(world, "state_h5");
    const double e_h5 = state_max_err(st, s_h5);
    const bool h5_present = std::filesystem::exists("state_h5.h5");

    ok_x = (e_leg < 1e-10) && (e_h5 < 1e-10) && h5_present;
    if (world.rank() == 0) {
      std::printf("ResponseStateX<ClosedShell> round-trip (%zu orbitals, NP=%d):\n",
                  st.x_alpha.size(), world.size());
      std::printf("  legacy   max_err = %.1e\n", e_leg);
      std::printf("  hdf5     max_err = %.1e   (state_h5.h5 present: %s)\n", e_h5,
                  h5_present ? "yes" : "no");
    }
  }

  // ---- ResponseStateXY<ClosedShell> — X and Y blocks (dynamic-alpha / Full-ES) ----
  {
    ResponseStateXY<ClosedShell> st;
    st.x_alpha.push_back(FunctionFactory<double, D>(world).f(g0));
    st.x_alpha.push_back(FunctionFactory<double, D>(world).f(g1));
    st.x_alpha.push_back(FunctionFactory<double, D>(world).f(g2));
    st.y_alpha.push_back(FunctionFactory<double, D>(world).f(g2));  // distinct Y
    st.y_alpha.push_back(FunctionFactory<double, D>(world).f(g0));
    st.y_alpha.push_back(FunctionFactory<double, D>(world).f(g1));
    for (auto& f : st.x_alpha) f.truncate();
    for (auto& f : st.y_alpha) f.truncate();

    unsetenv("MADRESPONSE_IO_HDF5");
    st.save(world, "statexy_leg");
    auto s_leg = ResponseStateXY<ClosedShell>::load(world, "statexy_leg");
    const double e_leg = statexy_max_err(st, s_leg);

    setenv("MADRESPONSE_IO_HDF5", "1", 1);
    st.save(world, "statexy_h5");
    auto s_h5 = ResponseStateXY<ClosedShell>::load(world, "statexy_h5");
    const double e_h5 = statexy_max_err(st, s_h5);
    const bool h5_present = std::filesystem::exists("statexy_h5.h5");

    ok_xy = (e_leg < 1e-10) && (e_h5 < 1e-10) && h5_present;
    if (world.rank() == 0) {
      std::printf("ResponseStateXY<ClosedShell> round-trip (%zu+%zu orbitals, NP=%d):\n",
                  st.x_alpha.size(), st.y_alpha.size(), world.size());
      std::printf("  legacy   max_err = %.1e\n", e_leg);
      std::printf("  hdf5     max_err = %.1e   (statexy_h5.h5 present: %s)\n", e_h5,
                  h5_present ? "yes" : "no");
    }
  }

  const bool ok = ok_x && ok_xy;
  if (world.rank() == 0)
    std::printf(ok ? " VERDICT: PASS — X and XY round-trip exact on both paths; "
                     "HDF5 auto-detected.\n"
                   : " VERDICT: FAIL — see errors above.\n");
  const int rc = ok ? 0 : 1;

  world.gop.fence();
  finalize();
  return rc;
}
