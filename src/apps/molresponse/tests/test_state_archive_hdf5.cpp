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

#include "../solvers/fd_save_load.hpp"  // detail_fd_save_load::check_writer_nproc
#include "../solvers/response_state.hpp"

#include <madness/mra/mra.h>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <stdexcept>
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

  // ---- Backend toggle: a stale twin must never shadow the newer save ----
  // (REVIEW finding 10: save writes the opt-in backend to the same basename
  // and previously left the other backend's file behind; auto-detect then
  // preferred the STALE .h5 forever.)
  bool ok_toggle = false;
  {
    ResponseStateX<ClosedShell> a, b;
    a.x_alpha.push_back(FunctionFactory<double, D>(world).f(g0));
    b.x_alpha.push_back(FunctionFactory<double, D>(world).f(g1));  // distinct
    a.x_alpha[0].truncate();
    b.x_alpha[0].truncate();

    // Run 1: HDF5 on — checkpoint A as .h5.
    setenv("MADRESPONSE_IO_HDF5", "1", 1);
    const IoBackend b1 = a.save(world, "state_toggle");
    const bool h5_first = std::filesystem::exists("state_toggle.h5");

    // Run 2: toggle OFF — save the NEWER state B natively at the same
    // basename. The stale .h5 twin must be removed, and an auto-detect load
    // must return B, not A.
    unsetenv("MADRESPONSE_IO_HDF5");
    const IoBackend b2 = b.save(world, "state_toggle");
    const bool h5_gone = !std::filesystem::exists("state_toggle.h5");
    auto got_auto = ResponseStateX<ClosedShell>::load(world, "state_toggle");
    const double e_newer = state_max_err(b, got_auto);
    // Explicit backend hint (what the metadata layer now records/threads).
    auto got_native = ResponseStateX<ClosedShell>::load(world, "state_toggle",
                                                        IoBackend::Native);
    const double e_native = state_max_err(b, got_native);

    // Run 3: toggle back ON — save A as .h5 again; the native parts must go.
    setenv("MADRESPONSE_IO_HDF5", "1", 1);
    const IoBackend b3 = a.save(world, "state_toggle");
    const bool native_gone = !std::filesystem::exists("state_toggle.00000");
    auto got_h5 = ResponseStateX<ClosedShell>::load(world, "state_toggle",
                                                    IoBackend::Hdf5);
    const double e_h5b = state_max_err(a, got_h5);

    ok_toggle = h5_first && (b1 == IoBackend::Hdf5) &&
                h5_gone && (b2 == IoBackend::Native) &&
                (e_newer < 1e-10) && (e_native < 1e-10) &&
                (b3 == IoBackend::Hdf5) && native_gone && (e_h5b < 1e-10);
    if (world.rank() == 0) {
      std::printf("backend toggle (stale-twin removal):\n");
      std::printf("  h5 written then removed on native re-save: %s/%s\n",
                  h5_first ? "yes" : "NO", h5_gone ? "yes" : "NO");
      std::printf("  auto-detect after toggle-off loads the NEWER state: "
                  "err=%.1e\n", e_newer);
      std::printf("  native parts removed on h5 re-save: %s  (h5 err=%.1e)\n",
                  native_gone ? "yes" : "NO", e_h5b);
    }
  }

  // ---- np-guard: native np-mismatch refused, hdf5 + legacy pass ----
  bool ok_npguard = false;
  {
    using detail_fd_save_load::check_writer_nproc;
    bool threw_native = false;
    try {
      check_writer_nproc(world, IoBackend::Native,
                         static_cast<int>(world.size()) + 1, "test", "arch");
    } catch (const std::runtime_error&) { threw_native = true; }
    bool ok_hdf5 = true, ok_legacy = true, ok_match = true;
    try {  // hdf5 blob is np-portable — mismatch must NOT throw
      check_writer_nproc(world, IoBackend::Hdf5,
                         static_cast<int>(world.size()) + 1, "test", "arch");
    } catch (...) { ok_hdf5 = false; }
    try {  // writer_nproc==0 = legacy entry — proceed
      check_writer_nproc(world, IoBackend::Native, 0, "test", "arch");
    } catch (...) { ok_legacy = false; }
    try {  // matching count — proceed
      check_writer_nproc(world, IoBackend::Native,
                         static_cast<int>(world.size()), "test", "arch");
    } catch (...) { ok_match = false; }
    ok_npguard = threw_native && ok_hdf5 && ok_legacy && ok_match;
    if (world.rank() == 0)
      std::printf("np-guard: native mismatch refused=%s  hdf5 pass=%s  "
                  "legacy pass=%s  match pass=%s\n",
                  threw_native ? "yes" : "NO", ok_hdf5 ? "yes" : "NO",
                  ok_legacy ? "yes" : "NO", ok_match ? "yes" : "NO");
  }

  // ---- atomic .h5 save: a stranded truncated tmp must not break anything ----
  // (REVIEW finding: .h5 used to be written in place — a kill mid-save
  // bricked restart. Now the writer goes tmp+rename; simulate the kill by
  // planting a garbage tmp and confirming save/load ignore+replace it.)
  bool ok_atomic = false;
  {
    ResponseStateX<ClosedShell> a;
    a.x_alpha.push_back(FunctionFactory<double, D>(world).f(g2));
    a.x_alpha[0].truncate();

    setenv("MADRESPONSE_IO_HDF5", "1", 1);
    a.save(world, "state_atomic");  // good checkpoint on disk
    if (world.rank() == 0) {        // simulate a kill mid-(re)save
      std::ofstream junk("state_atomic.h5.tmp");
      junk << "truncated garbage, not an HDF5 file";
    }
    world.gop.fence();
    // The good .h5 is untouched by the stranded tmp.
    auto got = ResponseStateX<ClosedShell>::load(world, "state_atomic",
                                                 IoBackend::Hdf5);
    const double e_before = state_max_err(a, got);
    // A re-save replaces the garbage tmp and installs atomically.
    a.save(world, "state_atomic");
    const bool tmp_gone = !std::filesystem::exists("state_atomic.h5.tmp");
    auto got2 = ResponseStateX<ClosedShell>::load(world, "state_atomic",
                                                  IoBackend::Hdf5);
    const double e_after = state_max_err(a, got2);
    ok_atomic = (e_before < 1e-10) && tmp_gone && (e_after < 1e-10);
    if (world.rank() == 0)
      std::printf("atomic h5 save: load-with-stranded-tmp err=%.1e  "
                  "re-save clears tmp=%s  err=%.1e\n",
                  e_before, tmp_gone ? "yes" : "NO", e_after);
  }

  const bool ok = ok_x && ok_xy && ok_toggle && ok_npguard && ok_atomic;
  if (world.rank() == 0)
    std::printf(ok ? " VERDICT: PASS — X and XY round-trip exact on both paths; "
                     "HDF5 auto-detected; backend toggle removes the stale "
                     "twin; np-guard + atomic tmp OK.\n"
                   : " VERDICT: FAIL — see errors above.\n");
  const int rc = ok ? 0 : 1;

  world.gop.fence();
  finalize();
  return rc;
}
