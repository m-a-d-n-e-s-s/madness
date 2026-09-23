// =========================================================================
// common_ops::bsh_shift — the level shift that keeps every BSH exponent real.
//
// make_bsh_operators builds mu_p = sqrt(-2 (eps_p + omega + shift)); the shift
// must make eps_p + omega + shift < 0 for EVERY occupied orbital. With
// localized orbitals the energies are Fock diagonals in LMO order (SCF.cc,
// do_localize), not sorted, so the highest one need not be last (review
// finding C11: the shift read eps(size-1)).
//
// Pure C++ on a Tensor; no World needed.
// =========================================================================

#include "../kernels/common_ops.hpp"

#include <cstdio>

namespace {

int failed = 0;

void expect(bool cond, const char *label) {
  std::printf("  [%s]  %s\n", cond ? "PASS" : "FAIL", label);
  if (!cond) ++failed;
}

madness::Tensor<double> energies(std::initializer_list<double> v) {
  madness::Tensor<double> t(static_cast<long>(v.size()));
  long i = 0;
  for (double e : v) t(i++) = e;
  return t;
}

bool all_exponents_real(const madness::Tensor<double> &eps, double omega) {
  const double shift = molresponse_v3::common_ops::bsh_shift(eps, omega);
  for (long p = 0; p < eps.size(); ++p)
    if (!(eps(p) + omega + shift < 0.0)) return false;
  return true;
}

} // namespace

int main() {
  using molresponse_v3::common_ops::bsh_shift;

  std::printf("=== bsh_shift: canonical (ascending) energies ===\n");
  {
    const auto eps = energies({-0.5, -0.3, -0.2});
    expect(bsh_shift(eps, 0.1) == 0.0, "no shift needed when HOMO + omega < 0");
    expect(all_exponents_real(eps, 0.25), "HOMO + omega >= 0: every exponent real");
  }

  std::printf("=== bsh_shift: localized (unsorted) energies ===\n");
  {
    const auto eps = energies({-0.5, -0.2, -0.3});  // highest is NOT last
    expect(std::abs(bsh_shift(eps, 0.25) - (-0.1)) < 1e-15,
           "shift taken from the highest energy, not the last");
    expect(all_exponents_real(eps, 0.25), "every exponent real");
  }

  std::printf("\n%s  (%d failures)\n", failed ? "FAILED" : "PASSED", failed);
  return failed ? 1 : 0;
}
