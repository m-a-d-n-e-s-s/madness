// =========================================================================
// Pin the protocol_key() format contract (doc 13).
//
// protocol_key is the join key FD and ES persistence agree on so properties
// can assert "same accuracy". The exact string format matters: it appears in
// archive filenames AND as a JSON key, so a future change to the formatter
// would silently break restart/property matching. Lock it down here.
//
// Pure C++; no MPI / no MADNESS World needed. Tests only the two-arg
// overload; the no-arg overload is a trivial pass-through to FunctionDefaults.
// =========================================================================

#include "../ResponseProtocol.hpp"

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

int failed = 0;

void expect_eq(const std::string &got, const std::string &want,
               const char *label) {
  if (got == want) {
    std::printf("  [PASS]  %s -> \"%s\"\n", label, got.c_str());
  } else {
    std::printf("  [FAIL]  %s -> got \"%s\"  want \"%s\"\n",
                label, got.c_str(), want.c_str());
    ++failed;
  }
}

void expect_true(bool ok, const char *label) {
  std::printf("  [%s]  %s\n", ok ? "PASS" : "FAIL", label);
  if (!ok) ++failed;
}

} // namespace

int main() {
  using molresponse_v3::protocol_key;

  std::printf("=== protocol_key(thresh, k) — canonical ramp ===\n");
  expect_eq(protocol_key(1e-2,  4),  "1e-02_k4",  "1e-2 / k=4 ");
  expect_eq(protocol_key(1e-4,  6),  "1e-04_k6",  "1e-4 / k=6 ");
  expect_eq(protocol_key(1e-6,  8),  "1e-06_k8",  "1e-6 / k=8 ");
  expect_eq(protocol_key(1e-8,  10), "1e-08_k10", "1e-8 / k=10");
  expect_eq(protocol_key(1e-10, 12), "1e-10_k12", "1e-10 / k=12");

  std::printf("=== protocol_key(thresh, k) — off-table ===\n");
  expect_eq(protocol_key(1e-5, 6),  "1e-05_k6",  "1e-5 / k=6 (user-specified thresh)");
  expect_eq(protocol_key(5e-7, 8),  "5e-07_k8",  "5e-7 / k=8 (non-power-of-10)");

  std::printf("=== apply_seed_start_rung (seed.start_rung contract) ===\n");
  {
    using molresponse_v3::apply_seed_start_rung;
    // fine + seed: ladder collapses to the finest rung.
    std::vector<double> p{1e-4, 1e-6};
    expect_true(apply_seed_start_rung(p, "fine", true) &&
                    p.size() == 1 && p[0] == 1e-6,
                "fine+seed [1e-4,1e-6] -> [1e-6]");
    // 3-rung ladder: finest wins even if unsorted input.
    p = {1e-6, 1e-4, 1e-8};
    expect_true(apply_seed_start_rung(p, "fine", true) &&
                    p.size() == 1 && p[0] == 1e-8,
                "fine+seed unsorted [1e-6,1e-4,1e-8] -> [1e-8]");
    // coarse (default): no-op, ladder untouched.
    p = {1e-4, 1e-6};
    expect_true(!apply_seed_start_rung(p, "coarse", true) && p.size() == 2,
                "coarse+seed -> full ladder unchanged");
    // fine WITHOUT a seed: no-op (caller warns; cold fine start is just a
    // single-rung protocol request).
    p = {1e-4, 1e-6};
    expect_true(!apply_seed_start_rung(p, "fine", false) && p.size() == 2,
                "fine without seed -> ignored, full ladder");
    // single-rung ladder: nothing to skip.
    p = {1e-6};
    expect_true(!apply_seed_start_rung(p, "fine", true) &&
                    p.size() == 1 && p[0] == 1e-6,
                "fine+seed single rung -> no-op");
    // typo'd mode: loud error, not a silent default.
    p = {1e-4, 1e-6};
    bool threw = false;
    try {
      apply_seed_start_rung(p, "finest", true);
    } catch (const std::invalid_argument &) {
      threw = true;
    }
    expect_true(threw, "unknown mode 'finest' throws invalid_argument");
  }

  std::printf("\n%s: %d failure(s)\n",
              failed == 0 ? "ALL PASS" : "FAILED", failed);
  return failed == 0 ? 0 : 1;
}
