#ifndef MOLRESPONSE_V3_SOLVERS_STATE_METRICS_HPP
#define MOLRESPONSE_V3_SOLVERS_STATE_METRICS_HPP

// ===========================================================================
// Per-state memory + convergence instrumentation (calc-manager foundation).
//
// Records, per response state at each protocol boundary:
//   coeffs  — total MRA coefficients across all blocks (Function::size,
//             a collective global sum). The molecule/protocol-specific
//             footprint the memory-scaling study needs (predict
//             high-protocol size from a measured low-protocol solve).
//   bytes   — coeffs × sizeof(double): the per-state memory footprint.
//   rss_gb  — worst-task resident set size (GiB) via gop.max. The
//             OOM-relevant number — exit 137 is the kernel killing on RSS
//             (see CLAUDE.md memory bottleneck). /proc/self/statm.
//   iters   — convergence iterations spent reaching this protocol.
//
// Wall-time per protocol is intentionally NOT here — it belongs to the
// calc-manager's scheduling loop, which owns the timer. These three
// (coeffs, rss, iters) are the low-coupling, high-value metrics that drop
// straight into the existing FD/ES save paths.
//
// measure_state() is COLLECTIVE on `world` (size_local sum + gop reduces);
// every rank must call it. One sum + one max reduce per call — cheap at
// protocol boundaries, not for the hot iteration loop.
// ===========================================================================

#include <madness/external/nlohmann_json/json.hpp>
#include <madness/world/MADworld.h>

#include <unistd.h>
#ifdef __APPLE__
#include <mach/mach.h>
#endif

#include <cstddef>
#include <fstream>

namespace molresponse_v3 {

/// Resident set size of THIS process in GiB. Linux: /proc/self/statm
/// (field 1 = resident pages). macOS: mach task_info (no procfs there —
/// returning 0.0 failed test_state_metrics' "RSS > 0" on macos CI).
/// Returns 0.0 only if the platform query itself fails.
inline double process_rss_gb() {
#ifdef __APPLE__
  mach_task_basic_info info;
  mach_msg_type_number_t count = MACH_TASK_BASIC_INFO_COUNT;
  if (task_info(mach_task_self(), MACH_TASK_BASIC_INFO,
                reinterpret_cast<task_info_t>(&info), &count) != KERN_SUCCESS)
    return 0.0;
  return static_cast<double>(info.resident_size) /
         (1024.0 * 1024.0 * 1024.0);
#else
  std::ifstream statm("/proc/self/statm");
  if (!statm) return 0.0;
  long total_pages = 0, resident_pages = 0;
  statm >> total_pages >> resident_pages;
  if (!statm) return 0.0;
  const long page = ::sysconf(_SC_PAGESIZE);
  const double bytes =
      static_cast<double>(resident_pages) * static_cast<double>(page);
  return bytes / (1024.0 * 1024.0 * 1024.0);
#endif
}

struct StateMetrics {
  std::size_t coeffs = 0;   // total MRA coefficients (all blocks, global sum)
  std::size_t bytes  = 0;   // coeffs × sizeof(double)
  double      rss_gb = 0.0; // worst-task RSS (gop.max), GiB
  int         iters  = 0;   // convergence iterations for this protocol
  double      wall_s = 0.0; // wall time for this (state,protocol) solve (R1b);
                            // 0.0 = not timed. Set by the caller, NOT by
                            // measure_state (which has no timer); uniform across
                            // FD / ES / VBC saves.

  nlohmann::json to_json() const {
    return nlohmann::json{
        {"coeffs", coeffs},
        {"bytes",  bytes},
        {"rss_gb", rss_gb},
        {"iters",  iters},
        {"wall_s", wall_s}};
  }
};

/// Collective. Sum the local coefficient counts over every block of the
/// state, reduce once; take the worst-task RSS. `iters` is carried through
/// from the caller (the solver's convergence count).
template <typename Storage>
StateMetrics measure_state(madness::World &world, const Storage &state,
                           int iters) {
  StateMetrics m;
  m.iters = iters;

  std::size_t coeffs_local = 0;
  for (const auto &f : state.flatten()) coeffs_local += f.size_local();
  world.gop.sum(coeffs_local);
  m.coeffs = coeffs_local;
  m.bytes  = coeffs_local * sizeof(double);

  double rss = process_rss_gb();
  world.gop.max(rss);
  m.rss_gb = rss;
  return m;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_STATE_METRICS_HPP
