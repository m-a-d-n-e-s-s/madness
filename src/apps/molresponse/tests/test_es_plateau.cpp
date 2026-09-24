// =========================================================================
// Per-root ES plateau detection (ConvergencePolicy::root_plateau,
// es_roots_stalled).
//
// The solve-level detector (ConvergencePolicy::plateau) watches the max over
// the active roots of drho/target and |dw|/target, and withholds a stall as
// soon as either track has met its target. A stuck ES root typically has a
// converged omega (second order) and a density change jittering above target
// (first order). |dw| is then met for every root, so plateau() never fires
// and the solve spends every remaining iteration going nowhere. The per-root
// test judges each root on its own tracks, and a met track does not block it.
//
// Pure C++; no World needed.
// =========================================================================

#include "../solvers/convergence_policy.hpp"

#include <cstdio>
#include <limits>
#include <vector>

using namespace molresponse_v3;

namespace {

int failed = 0;

void expect(bool cond, const char *label) {
  std::printf("  [%s]  %s\n", cond ? "PASS" : "FAIL", label);
  if (!cond) ++failed;
}

using Track = std::vector<double>;

Track flat(double v, int n) { return Track(static_cast<size_t>(n), v); }

// v0, v0*f, v0*f^2, ...
Track geometric(double v0, double f, int n) {
  Track t;
  double v = v0;
  for (int i = 0; i < n; ++i) { t.push_back(v); v *= f; }
  return t;
}

} // namespace

int main() {
  ConvergencePolicy p;  // stall_window 6, stall_ratio 0.10
  const int n = p.stall_window + 1;
  const double inf = std::numeric_limits<double>::infinity();

  std::printf("=== the gap: density flat above target, |dw| met ===\n");
  {
    const std::vector<Track> stuck = {flat(3.0, n), flat(0.5, n)};
    expect(!p.plateau(stuck),
           "solve-level plateau() withholds: the met |dw| track blocks it");
    expect(p.root_plateau(stuck),
           "root_plateau() fires: every unmet track is flat");
  }

  std::printf("=== root_plateau: what must not fire ===\n");
  expect(!p.root_plateau({geometric(3.0, 0.9, n), flat(0.5, n)}),
         "density still falling (10% per iter) -> not a plateau");
  expect(!p.root_plateau({flat(0.8, n), flat(0.5, n)}),
         "both tracks met -> converged, not a plateau");
  expect(!p.root_plateau({flat(3.0, n - 1), flat(0.5, n - 1)}),
         "history shorter than window+1 -> no verdict");
  {
    Track t = flat(3.0, n);
    t.back() = inf;
    expect(!p.root_plateau({t, flat(0.5, n)}),
           "newest entry not measured -> inconclusive");
    Track d = flat(3.0, n);
    d.front() = inf;
    expect(!p.root_plateau({d, flat(0.5, n)}),
           "unmet track with no measurement window iters ago -> inconclusive");
    Track m = flat(0.5, n);
    m.front() = inf;
    expect(p.root_plateau({flat(3.0, n), m}),
           "a met track needs no history: its old sentinel does not block");
  }
  {
    ConvergencePolicy off = p;
    off.stall_window = 0;
    expect(!off.root_plateau({flat(3.0, n), flat(0.5, n)}),
           "stall_window 0 disables the per-root test too");
  }
  expect(p.root_plateau({flat(3.0, n), flat(2.0, n)}),
         "both tracks flat above target -> plateau (same as plateau())");

  std::printf("=== es_roots_stalled: the solve-level verdict ===\n");
  // active, converged, plateaued per slot
  expect(es_roots_stalled({1, 1}, {1, 0}, {0, 1}),
         "one root converged, the other plateaued -> stall");
  expect(!es_roots_stalled({1, 1}, {0, 0}, {0, 1}),
         "the other root still iterating -> no stall");
  expect(!es_roots_stalled({1, 1}, {1, 1}, {0, 0}),
         "all converged -> no stall (that is convergence)");
  expect(es_roots_stalled({0, 1}, {0, 0}, {0, 1}),
         "a locked (inactive) root is ignored");
  expect(!es_roots_stalled({0, 0}, {0, 0}, {0, 0}),
         "no active root -> no stall");
  expect(es_roots_stalled({1, 1}, {0, 0}, {1, 1}),
         "every root plateaued -> stall");
  expect(!es_roots_stalled({1, 1}, {1, 0}, {1, 0}),
         "a converged root flagged plateaued does not count as stuck");

  std::printf("\n%s: %d failure(s)\n", failed == 0 ? "ALL PASS" : "FAILED",
              failed);
  return failed == 0 ? 0 : 1;
}
