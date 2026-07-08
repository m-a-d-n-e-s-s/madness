// ===========================================================================
// Unit tests for the non-collective parts of state_metrics.hpp:
//   - process_rss_gb() returns a plausible positive RSS on Linux
//   - StateMetrics::to_json carries the expected keys/values
// The collective measure_state() path needs a World and is exercised on the
// allocation via the FD/ES save paths.
// ===========================================================================

#include "../solvers/state_metrics.hpp"

#include <cstdio>

using namespace molresponse_v3;

namespace {
int failed = 0;
#define EXPECT(cond, label)                                                \
  do {                                                                     \
    if (cond) { std::printf("  [PASS]  %s\n", label); }                    \
    else      { std::printf("  [FAIL]  %s\n", label); ++failed; }          \
  } while (0)
} // namespace

int main() {
  std::printf("=== process_rss_gb ===\n");
  const double rss = process_rss_gb();
  // On Linux /proc/self/statm exists; this process has nonzero resident
  // memory and certainly less than a terabyte.
  EXPECT(rss > 0.0,      "RSS > 0 (this process is resident)");
  EXPECT(rss < 1024.0,   "RSS < 1 TiB (sanity upper bound)");
  std::printf("    measured rss_gb = %.4f\n", rss);

  std::printf("=== StateMetrics::to_json ===\n");
  StateMetrics m;
  m.coeffs = 123456;
  m.bytes  = 123456 * sizeof(double);
  m.rss_gb = 2.5;
  m.iters  = 7;
  auto j = m.to_json();
  EXPECT(j["coeffs"] == 123456,                 "coeffs round-trips");
  EXPECT(j["bytes"]  == 123456 * sizeof(double), "bytes = coeffs × 8");
  EXPECT(j["rss_gb"] == 2.5,                     "rss_gb round-trips");
  EXPECT(j["iters"]  == 7,                       "iters round-trips");

  std::printf("\n%s: %d failure(s)\n",
              failed == 0 ? "ALL PASS" : "FAILED", failed);
  return failed == 0 ? 0 : 1;
}
