// =========================================================================
// test_function_hdf5.cpp — P2 Layer A: three-way Function round-trip A/B.
//
// For each (k, thresh) in a large box, round-trips a known real_function_3d
// through THREE paths and reports size + write/read wall-time + round-trip error:
//   * legacy      — madness::save / ParallelInputArchive<BinaryFstream> (native)
//   * structured  — save/load_function_hdf5        (/meta+/keys+/coeffs datasets)
//   * archive     — save/load_function_archive_hdf5 (MADNESS-serialized HDF5 blob)
// Aggregate PASS only if every path round-trips exact (≤1e-10) for every config.
//
// Edit SWEEP / BOXL below. NP=1, double (docs/30 §4 P2 Layer A). Built only when
// MADNESS_ENABLE_HDF5=ON.
// =========================================================================

#include "../GroundState.hpp"
#include "../solvers/function_hdf5_io.hpp"

#include <madness/chem/molecule.h>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/misc/info.h>  // commandlineparser
#include <madness/mra/mra.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <string>
#include <utility>
#include <vector>

using namespace madness;
using namespace molresponse_v3;  // GroundState, etc.
namespace fs = std::filesystem;

static const std::size_t D = 3;

// --- knobs (edit here) ---------------------------------------------------
static const double BOXL = 100.0;  // cubic cell [-BOXL/2, BOXL/2] (bohr)
static const int GZ_LEVEL = 1;     // gzip level for the archive+gz path (1..9)
static const std::vector<std::pair<long, double>> SWEEP = {
    {6, 1e-4}, {8, 1e-6}, {10, 1e-8}};

static double testfn(const Vector<double, D>& r) {
  const double x = r[0] - 0.7, y = r[1] + 0.3, z = r[2] - 0.2;
  const double r2 = x * x + y * y + z * z;
  return std::exp(-1.3 * r2);
}

static long long sum_sizes(const std::string& stem) {
  long long total = 0;
  std::error_code ec;
  for (const auto& e : fs::directory_iterator(fs::current_path(), ec)) {
    const std::string name = e.path().filename().string();
    if (name.rfind(stem, 0) == 0 && fs::is_regular_file(e, ec))
      total += static_cast<long long>(fs::file_size(e, ec));
  }
  return total;
}

struct PathResult {
  long long bytes;
  double tw, tr, err;
};
struct RunResult {
  long k;
  double thresh, norm_ref;
  PathResult legacy, structured, archive, archive_gz;
  bool ok;
};

// One (k, thresh) three-way round-trip. Functions are local to this helper, so
// they destruct on return — before finalize() (avoids the shutdown mutex abort).
static RunResult run_one(World& world, long k, double thresh, double L) {
  FunctionDefaults<D>::set_k(int(k));
  FunctionDefaults<D>::set_thresh(thresh);
  FunctionDefaults<D>::set_refine(true);
  FunctionDefaults<D>::set_initial_level(2);
  FunctionDefaults<D>::set_truncate_mode(0);
  FunctionDefaults<D>::set_cubic_cell(-L / 2, L / 2);

  Function<double, D> f = FunctionFactory<double, D>(world).f(testfn);
  f.truncate();
  const double norm_ref = f.norm2();

  const std::string legstem = "func_legacy";
  const std::string strpath = "func_struct.mad.h5";
  const std::string arcpath = "func_arc.mad.h5";

  // --- legacy (native archive) ---
  double t = wall_time();
  madness::save(f, legstem);
  const double tw_leg = wall_time() - t;
  t = wall_time();
  Function<double, D> f_leg;
  {
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
        world, legstem.c_str(), 1);
    ar& f_leg;
  }
  const double tr_leg = wall_time() - t;

  // --- structured HDF5 (NP=1 only: interop path; no multi-rank gather) ---
  double tw_str = 0, tr_str = 0, err_str = -1.0;  // err < 0 => skipped
  long long sz_str = 0;
  if (world.size() == 1) {
    t = wall_time();
    molresponse_v3::save_function_hdf5(f, strpath);
    tw_str = wall_time() - t;
    t = wall_time();
    Function<double, D> f_str =
        molresponse_v3::load_function_hdf5<double, D>(world, strpath);
    tr_str = wall_time() - t;
    err_str = (f - f_str).norm2();
    sz_str = sum_sizes("func_struct");
  }

  // --- archive-backend HDF5 ---
  t = wall_time();
  molresponse_v3::save_function_archive_hdf5(f, arcpath);
  const double tw_arc = wall_time() - t;
  t = wall_time();
  Function<double, D> f_arc =
      molresponse_v3::load_function_archive_hdf5<double, D>(world, arcpath);
  const double tr_arc = wall_time() - t;

  // --- archive-backend HDF5 + gzip ---
  const std::string gzpath = "func_arcgz.mad.h5";
  t = wall_time();
  molresponse_v3::save_function_archive_hdf5(f, gzpath, GZ_LEVEL);
  const double tw_gz = wall_time() - t;
  t = wall_time();
  Function<double, D> f_gz =
      molresponse_v3::load_function_archive_hdf5<double, D>(world, gzpath);
  const double tr_gz = wall_time() - t;

  const double err_leg = (f - f_leg).norm2();
  const double err_arc = (f - f_arc).norm2();
  const double err_gz = (f - f_gz).norm2();

  const long long sz_leg = sum_sizes(legstem);
  const long long sz_arc = sum_sizes("func_arc.");   // trailing '.' excludes func_arcgz
  const long long sz_gz = sum_sizes("func_arcgz");

  const double tol = 1e-10;
  // err_str < 0 means the structured path was skipped (NP>1) -> counts as ok.
  const bool ok = (err_leg < tol) && (err_str < tol) && (err_arc < tol) &&
                  (err_gz < tol);
  return {k,
          thresh,
          norm_ref,
          {sz_leg, tw_leg, tr_leg, err_leg},
          {sz_str, tw_str, tr_str, err_str},
          {sz_arc, tw_arc, tr_arc, err_arc},
          {sz_gz, tw_gz, tr_gz, err_gz},
          ok};
}

static void print_block(const RunResult& r) {
  std::printf("\n[k=%ld  thresh=%.1e]   norm2=%.6f\n", r.k, r.thresh, r.norm_ref);
  std::printf("  %-12s %12s %10s %10s %10s\n", "path", "bytes", "write(s)",
              "read(s)", "err");
  auto row = [](const char* name, const PathResult& p) {
    std::printf("  %-12s %12lld %10.4f %10.4f %10.1e\n", name, p.bytes, p.tw,
                p.tr, p.err);
  };
  row("legacy", r.legacy);
  if (r.structured.err >= 0)
    row("structured", r.structured);
  else
    std::printf("  %-12s %12s %10s %10s %10s\n", "structured", "-", "-", "-",
                "skip(NP>1)");
  row("archive", r.archive);
  row("archive+gz", r.archive_gz);
  if (r.archive.bytes > 0)
    std::printf("  (gzip level %d: %.3f x the uncompressed archive size)\n",
                GZ_LEVEL, double(r.archive_gz.bytes) / double(r.archive.bytes));
}

// Real-orbital mode (--archive=<moldft restart>): load the actual ground state
// and round-trip EACH orbital through legacy + the archive-backend HDF5 path
// (true, irregular orbital trees — the realistic restart payload). Multi-rank
// OK (the archive path gathers/distributes). Scaffolding mirrors
// test_v3_fd_skeleton's ground-state loader.
static bool run_real_mo(World& world, const std::string& archive_path) {
  auto header = GroundState::read_archive_header(world, archive_path);
  FunctionDefaults<D>::set_cubic_cell(-header.L / 2, header.L / 2);
  FunctionDefaults<D>::set_k(header.k);
  FunctionDefaults<D>::set_thresh(1e-6);  // round-trip is exact regardless of thresh

  Molecule molecule;
  auto archive_dir = fs::path(archive_path).parent_path();
  for (const auto& name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
    auto cand = archive_dir / name;
    if (fs::exists(cand)) {
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
  const auto& mos = gs.orbitals();
  const int nmo = static_cast<int>(mos.size());

  if (world.rank() == 0)
    std::printf("\n=== real ground-state MOs: %d orbitals, k=%d, L=%.2f, NP=%d ===\n"
                "  %4s %12s %12s %12s\n",
                nmo, header.k, header.L, world.size(), "mo", "norm2",
                "err(legacy)", "err(archive)");

  double tw_leg = 0, tr_leg = 0, tw_arc = 0, tr_arc = 0, max_err = 0;
  long long sz_leg = 0, sz_arc = 0;
  for (int i = 0; i < nmo; ++i) {
    const auto& mo = mos[i];
    double t = wall_time();
    madness::save(mo, "mo_legacy");
    tw_leg += wall_time() - t;
    t = wall_time();
    Function<double, D> f_leg;
    {
      archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
          world, "mo_legacy", 1);
      ar& f_leg;
    }
    tr_leg += wall_time() - t;

    t = wall_time();
    molresponse_v3::save_function_archive_hdf5(mo, "mo_arc.mad.h5");
    tw_arc += wall_time() - t;
    t = wall_time();
    Function<double, D> f_arc =
        molresponse_v3::load_function_archive_hdf5<double, D>(world, "mo_arc.mad.h5");
    tr_arc += wall_time() - t;

    const double n0 = mo.norm2();
    const double e_leg = (mo - f_leg).norm2();
    const double e_arc = (mo - f_arc).norm2();
    max_err = std::max({max_err, e_leg, e_arc});
    sz_leg += sum_sizes("mo_legacy");
    sz_arc += sum_sizes("mo_arc.");
    if (world.rank() == 0)
      std::printf("  %4d %12.6f %12.1e %12.1e\n", i, n0, e_leg, e_arc);
  }

  const bool ok = max_err < 1e-10;
  if (world.rank() == 0) {
    std::printf("  --------------------------------------------------------\n");
    std::printf("  totals over %d MOs:  legacy %lld B (w %.3f r %.3f) | "
                "archive %lld B (w %.3f r %.3f)\n",
                nmo, sz_leg, tw_leg, tr_leg, sz_arc, tw_arc, tr_arc);
    std::printf(ok ? " VERDICT: PASS — every real MO round-trips exact "
                     "(max_err %.1e).\n"
                   : " VERDICT: FAIL — a real MO drifted (max_err %.1e).\n",
                max_err);
  }
  return ok;
}

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  startup(world, argc, argv);
  std::cout.precision(8);

  // --archive=<moldft restart> => real-orbital mode; otherwise the analytic sweep.
  commandlineparser parser(argc, argv);
  if (parser.key_exists("archive")) {
    const bool ok = run_real_mo(world, parser.value_raw("archive"));
    world.gop.fence();
    finalize();
    return ok ? 0 : 1;
  }

  if (world.rank() == 0)
    std::printf("NP=%d  (legacy + archive[+gz] run at any rank count; "
                "structured is NP=1-only)\n",
                world.size());

  std::vector<RunResult> results;
  for (const auto& cfg : SWEEP)
    results.push_back(run_one(world, cfg.first, cfg.second, BOXL));

  bool all_ok = true;
  for (const auto& r : results) all_ok = all_ok && r.ok;

  if (world.rank() == 0) {
    std::printf("\n==========================================================\n");
    std::printf(" P2 Layer A — Function round-trip A/B   (box [%.0f,%.0f]^3, NP=%d)\n",
                -BOXL / 2, BOXL / 2, world.size());
    std::printf("==========================================================");
    for (const auto& r : results) print_block(r);
    std::printf("\n==========================================================\n");
    if (all_ok)
      std::printf(" VERDICT: PASS — every path round-trips exact for every (k,thresh).\n");
    else
      std::printf(" VERDICT: FAIL — a path drifted (>1e-10).\n");
    std::printf("==========================================================\n");
  }

  world.gop.fence();
  finalize();
  return all_ok ? 0 : 1;
}
