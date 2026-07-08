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
#include <string>

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

  std::printf("\n%s: %d failure(s)\n",
              failed == 0 ? "ALL PASS" : "FAILED", failed);
  return failed == 0 ? 0 : 1;
}
