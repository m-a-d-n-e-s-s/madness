#ifndef MOLRESPONSE_V3_KERNELS_COMMON_OPS_HPP
#define MOLRESPONSE_V3_KERNELS_COMMON_OPS_HPP

// =========================================================================
// Kernel-agnostic per-state operator wrappers used by EVERY (Type, Shell)
// specialization of `Kernels<…>`:
//
//   apply_kinetic       — −½∇² applied to a vector of orbitals
//   bsh_shift           — level-shift so the BSH μ² stays positive
//   make_bsh_operators  — vector of per-orbital BSH operators
//   apply_exchange      — compact wrapper around madness::Exchange
//                         + multiworld-efficient-row algorithm
//
// Previously these lived inside `namespace detail_tda` at the top of
// tda.hpp — that was a misnomer; static.hpp and full.hpp both reach
// across into them. They're shared kernel infrastructure, not
// TDA-specific. Lifted here so each kernel header includes one neutral
// location.
//
// Also hosts the `poperatorT` typedef (shared pointer to a 3D
// SeparatedConvolution) that the BSH and Coulomb operators use.
// =========================================================================

#include <madness/chem/SCFOperators.h>   // Exchange
#include <madness/mra/mra.h>              // real_function_3d, real_derivative_3d
#include <madness/mra/operator.h>         // BSHOperatorPtr3D
#include <madness/tensor/tensor.h>        // Tensor

#include <cmath>
#include <memory>
#include <vector>

namespace molresponse_v3 {

using poperatorT = std::shared_ptr<madness::real_convolution_3d>;

namespace common_ops {

/// Pointwise-product sum of two equal-length vecfuncs:
///
///     result = Σ_p left[p] · right[p]
///
/// This is the operation every `compute_density` performs: pair the
/// kernel's φ-side flat vecfunc with `state.flatten()` and reduce to a
/// single response-density function.
///
/// Currently a thin pass-through to `madness::dot(world, left, right)`
/// so the per-kernel call sites become uniform without changing the
/// numerical behavior of the original code. (A hand-rolled "batched
/// mul + single fence" variant was tried first; it produced subtly
/// different results — internal truncate ordering vs the per-element
/// `a[i]*b[i]` MADNESS dot does — and was rolled back. If we ever
/// need fence-fused multi-pair density assembly we can revisit, but
/// for correctness this passthrough is the trusted form.)
///
/// Callers should still `scale(...)` and `truncate()` after the call,
/// matching the original idiom and the post-scale truncation
/// threshold that the existing Dalton-validated runs were built on.
inline madness::real_function_3d
dot(madness::World &world,
    const std::vector<madness::real_function_3d> &left,
    const std::vector<madness::real_function_3d> &right) {
  return madness::dot(world, left, right);
}

inline std::vector<madness::real_function_3d>
apply_kinetic(madness::World &world,
              const std::vector<madness::real_function_3d> &v) {
  if (v.empty()) return {};
  std::vector<madness::real_function_3d> result;
  for (int d = 0; d < 3; ++d) {
    madness::real_derivative_3d D(world, d);
    auto dv = apply(world, D, v);
    auto dv2 = apply(world, D, dv);
    if (result.empty()) result = std::move(dv2);
    else gaxpy(world, 1.0, result, 1.0, dv2);
  }
  scale(world, result, -0.5);
  truncate(world, result);
  return result;
}

/// Pick a level-shift so the per-orbital BSH μ² stays positive.
/// `eps` carries OCCUPIED orbital energies only (HOMO last). When
/// `HOMO + omega` would push μ² = -2(ε_HOMO + omega) ≤ 0, return
/// a shift that pushes the effective energy below zero by `guard`.
inline double bsh_shift(const madness::Tensor<double> &eps, double omega) {
  constexpr double guard = 0.05;
  const double homo_shifted = eps(eps.size() - 1) + omega;
  return (homo_shifted >= 0.0) ? -guard - homo_shifted : 0.0;
}

inline std::vector<poperatorT>
make_bsh_operators(madness::World &world,
                   const madness::Tensor<double> &eps, double omega,
                   double lo) {
  const double tol = madness::FunctionDefaults<3>::get_thresh();
  const double shift = bsh_shift(eps, omega);
  std::vector<poperatorT> ops(eps.size());
  for (long p = 0; p < eps.size(); ++p) {
    const double mu = std::sqrt(-2.0 * (eps(p) + omega + shift));
    ops[p] = poperatorT(madness::BSHOperatorPtr3D(world, mu, lo, tol));
  }
  return ops;
}

/// Compact wrapper around `madness::Exchange::set_bra_and_ket(bra, ket)`
/// + apply. Each call constructs an Exchange operator with the
/// multiworld-efficient-row algorithm and applies it to `apply_to`.
/// Used 4-9× per compute_gamma / compute_V0x; one-line replacement
/// for the 5-line Exchange<double, 3> ... set_bra_and_ket ... set_algorithm
/// ... apply pattern.
inline std::vector<madness::real_function_3d>
apply_exchange(madness::World &world,
               const std::vector<madness::real_function_3d> &bra,
               const std::vector<madness::real_function_3d> &ket,
               const std::vector<madness::real_function_3d> &apply_to,
               double lo) {
  madness::Exchange<double, 3> K(world, lo);
  K.set_bra_and_ket(bra, ket);
  K.set_algorithm(madness::Exchange<double, 3>::
                      ExchangeAlgorithm::multiworld_efficient_row);
  return K(apply_to);
}

/// Build a ground-state exchange operator K = K[mos, mos] once, so callers can
/// reuse it across iterations instead of reconstructing it (which deep-copies
/// `mos` in Exchange::set_bra_and_ket on every call). Returns nullptr for an
/// empty mos (e.g. the beta block of a closed-shell run). Built at the active
/// k/thresh — rebuild when the ground state is re-prepared to a new protocol.
inline std::shared_ptr<madness::Exchange<double, 3>>
make_ground_exchange(madness::World &world,
                     const std::vector<madness::real_function_3d> &mos,
                     double lo) {
  if (mos.empty()) return nullptr;
  auto K = std::make_shared<madness::Exchange<double, 3>>(world, lo);
  K->set_bra_and_ket(mos, mos);
  K->set_algorithm(madness::Exchange<double, 3>::
                       ExchangeAlgorithm::multiworld_efficient_row);
  return K;
}

/// Apply a prebuilt ground exchange operator (from make_ground_exchange) to
/// `apply_to`. Falls back to building a fresh K[mos, mos] if the cached operator
/// is null — identical numerics to the old per-call path, so a ground state
/// built without the cache still works (defensive; the build helpers set it).
inline std::vector<madness::real_function_3d>
apply_ground_exchange(madness::World &world,
                      const std::shared_ptr<madness::Exchange<double, 3>> &K,
                      const std::vector<madness::real_function_3d> &mos,
                      const std::vector<madness::real_function_3d> &apply_to,
                      double lo) {
  if (K) return (*K)(apply_to);
  return apply_exchange(world, mos, mos, apply_to, lo);
}

} // namespace common_ops
} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_COMMON_OPS_HPP
