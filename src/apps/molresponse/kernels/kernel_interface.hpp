#ifndef MOLRESPONSE_V3_KERNELS_KERNEL_INTERFACE_HPP
#define MOLRESPONSE_V3_KERNELS_KERNEL_INTERFACE_HPP

// =========================================================================
// Enforce the "Kernel" interface contract for every Kernels<Type, Shell>
// specialization. Two refinement levels:
//
//   FDKernel<K>   — basic response-density kernel. Required by FDSolver.
//                   Members:
//                     using State                       — Storage type
//                     compute_density(w, g0, s)       → real_function_3d
//                     compute_gamma  (w, g0, s, rho)  → State
//                     apply_g(w, g0, S1, S2, S3)      → {State, rho}
//                     compute_V0x    (w, g0, s)       → State
//                     compute_E0x    (w, g0, s)       → State   (no-diag)
//                     bsh_apply      (w, g0, s, theta, omega) → State
//                     compute_residual_norm(w, s_old, s_new) → double
//
//   ESKernel<K>   — FDKernel + the extra T0x / Efull pieces the ES
//                   subspace eigenproblem needs. Required by ESSolver.
//                   Adds:
//                     compute_T0x        (w, g0, s)        → State
//                     compute_E0x_full   (w, g0, s)        → State  (with diag)
//
// (θ / Λ assembly is shell-agnostic and lives in kernels/assembly.hpp as
//  the State-generic `assemble_theta` / `assemble_lambda` free functions,
//  so it isn't part of the kernel interface.)
//
// Current specializations and their required level:
//   Kernels<Static, *>   — FDKernel  (FD-only)
//   Kernels<Full,   *>   — FDKernel  (today FD; promote to ESKernel when
//                                     ESSolver<Full, *> is added)
//   Kernels<TDA,    *>   — ESKernel  (ES-only)
//
// Implementation:
//   * C++17 path (codebase default): SFINAE `has_<method>` traits + a
//     compound `is_fd_kernel_v` / `is_es_kernel_v` bool, with per-method
//     static_asserts so a missing method gets a single surgical error.
//   * C++20 path (when CMAKE_CXX_STANDARD=20): a `FDKernel`/`ESKernel`
//     concept refinement, also wired into the per-method asserts and
//     usable directly as a template constraint on solvers.
//
// Both paths trip at compile time at the bottom of each kernel header,
// so adding a new (Type, Shell) without implementing the full interface
// is caught immediately — not as a fifty-line template error at the
// solver call site.
// =========================================================================

#include "tda.hpp"                         // ResponseGroundState

#include <madness/mra/mra.h>

#include <type_traits>
#include <utility>

#if defined(__cpp_concepts) && __cpp_concepts >= 201907L
#include <concepts>
#endif

namespace molresponse_v3::detail_kernel {

// --- detection traits: one per required method ------------------------
// Each trait probes a specific call expression via std::void_t. If the
// expression is well-formed for K, has_X<K>::value is true; otherwise
// SFINAE selects the false primary template.

#define MV3_DEFINE_KERNEL_TRAIT(name, expr)                              \
  template <typename, typename = std::void_t<>>                          \
  struct has_##name : std::false_type {};                                \
  template <typename K>                                                  \
  struct has_##name<K, std::void_t<decltype(expr)>> : std::true_type {}; \
  template <typename K>                                                  \
  inline constexpr bool has_##name##_v = has_##name<K>::value

// Probe-expression aliases keep the trait bodies readable.
//   W  = madness::World&
//   G0 = const ResponseGroundState&
//   S  = const typename K::State&
//   R  = const madness::real_function_3d&
//   D  = double
#define MV3_PROBE_W   std::declval<madness::World&>()
#define MV3_PROBE_G0  std::declval<const ResponseGroundState&>()
#define MV3_PROBE_S   std::declval<const typename K::State&>()
#define MV3_PROBE_R   std::declval<const madness::real_function_3d&>()
#define MV3_PROBE_D   std::declval<double>()

// has_state: K::State must name a type.
template <typename, typename = std::void_t<>>
struct has_state : std::false_type {};
template <typename K>
struct has_state<K, std::void_t<typename K::State>> : std::true_type {};
template <typename K>
inline constexpr bool has_state_v = has_state<K>::value;

// FDKernel methods
MV3_DEFINE_KERNEL_TRAIT(compute_density,
    K::compute_density(MV3_PROBE_W, MV3_PROBE_G0, MV3_PROBE_S));
MV3_DEFINE_KERNEL_TRAIT(compute_gamma,
    K::compute_gamma(MV3_PROBE_W, MV3_PROBE_G0, MV3_PROBE_S, MV3_PROBE_R));
MV3_DEFINE_KERNEL_TRAIT(apply_g,
    K::apply_g(MV3_PROBE_W, MV3_PROBE_G0, MV3_PROBE_S, MV3_PROBE_S,
               MV3_PROBE_S));
MV3_DEFINE_KERNEL_TRAIT(compute_V0x,
    K::compute_V0x(MV3_PROBE_W, MV3_PROBE_G0, MV3_PROBE_S));
MV3_DEFINE_KERNEL_TRAIT(compute_E0x,
    K::compute_E0x(MV3_PROBE_W, MV3_PROBE_G0, MV3_PROBE_S));
MV3_DEFINE_KERNEL_TRAIT(bsh_apply,
    K::bsh_apply(MV3_PROBE_W, MV3_PROBE_G0, MV3_PROBE_S, MV3_PROBE_S,
                 MV3_PROBE_D));
MV3_DEFINE_KERNEL_TRAIT(compute_residual_norm,
    K::compute_residual_norm(MV3_PROBE_W, MV3_PROBE_S, MV3_PROBE_S));

// ESKernel-only methods (extra subspace pieces; Λ assembly itself is
// shell-agnostic — see kernels/assembly.hpp).
MV3_DEFINE_KERNEL_TRAIT(compute_T0x,
    K::compute_T0x(MV3_PROBE_W, MV3_PROBE_G0, MV3_PROBE_S));
MV3_DEFINE_KERNEL_TRAIT(compute_E0x_full,
    K::compute_E0x_full(MV3_PROBE_W, MV3_PROBE_G0, MV3_PROBE_S));

#undef MV3_PROBE_W
#undef MV3_PROBE_G0
#undef MV3_PROBE_S
#undef MV3_PROBE_R
#undef MV3_PROBE_D
#undef MV3_DEFINE_KERNEL_TRAIT

template <typename K>
inline constexpr bool is_fd_kernel_v =
    has_state_v<K>            &&
    has_compute_density_v<K>  &&
    has_compute_gamma_v<K>    &&
    has_apply_g_v<K>          &&
    has_compute_V0x_v<K>      &&
    has_compute_E0x_v<K>      &&
    has_bsh_apply_v<K>        &&
    has_compute_residual_norm_v<K>;

template <typename K>
inline constexpr bool is_es_kernel_v =
    is_fd_kernel_v<K>           &&
    has_compute_T0x_v<K>        &&
    has_compute_E0x_full_v<K>;

} // namespace molresponse_v3::detail_kernel

// ----------------------------------------------------------------------
// C++20 concepts (when available) — parallel formulation. Useful for
// constraining solver template parameters directly:
//     template <typename T, typename S>
//       requires ESKernel<Kernels<T, S>>
//     class ESSolver { ... };
//
// Gated by __cpp_concepts; falls away silently under C++17. The
// always-on SFINAE traits above still cover correctness under C++17.
// ----------------------------------------------------------------------
#if defined(__cpp_concepts) && __cpp_concepts >= 201907L
namespace molresponse_v3 {

template <typename K>
concept FDKernel = requires(
    madness::World& w,
    const ResponseGroundState& g0,
    const typename K::State& s,
    const madness::real_function_3d& rho,
    double omega)
{
  typename K::State;
  { K::compute_density (w, g0, s)            } -> std::convertible_to<madness::real_function_3d>;
  { K::compute_gamma   (w, g0, s, rho)       } -> std::convertible_to<typename K::State>;
  { K::apply_g         (w, g0, s, s, s)      };   // -> std::pair<State, real_function_3d>
  { K::compute_V0x     (w, g0, s)            } -> std::convertible_to<typename K::State>;
  { K::compute_E0x     (w, g0, s)            } -> std::convertible_to<typename K::State>;
  { K::bsh_apply       (w, g0, s, s, omega)  } -> std::convertible_to<typename K::State>;
  { K::compute_residual_norm(w, s, s)        } -> std::convertible_to<double>;
};

template <typename K>
concept ESKernel = FDKernel<K> && requires(
    madness::World& w,
    const ResponseGroundState& g0,
    const typename K::State& s)
{
  { K::compute_T0x     (w, g0, s)            } -> std::convertible_to<typename K::State>;
  { K::compute_E0x_full(w, g0, s)            } -> std::convertible_to<typename K::State>;
};

} // namespace molresponse_v3
#endif // __cpp_concepts

// ----------------------------------------------------------------------
// Public assertion macros — emit one static_assert per required method,
// each with a focused error message so a missing method gets a single
// surgical diagnostic, not a fifty-line template error.
//
// Variadic: the kernel type expression contains a comma
// (`Kernels<Static, ClosedShell>`) which the preprocessor would
// otherwise read as an extra macro argument. `__VA_ARGS__` lets the
// caller pass the whole type unmolested.
// ----------------------------------------------------------------------
#define MV3_ASSERT_FD_KERNEL(...)                                                                    \
  static_assert(::molresponse_v3::detail_kernel::has_state_v<__VA_ARGS__>,                           \
                #__VA_ARGS__ " is missing required `using State`");                                  \
  static_assert(::molresponse_v3::detail_kernel::has_compute_density_v<__VA_ARGS__>,                 \
                #__VA_ARGS__ " is missing `compute_density(w, g0, s)`");                             \
  static_assert(::molresponse_v3::detail_kernel::has_compute_gamma_v<__VA_ARGS__>,                   \
                #__VA_ARGS__ " is missing `compute_gamma(w, g0, s, rho)`");                          \
  static_assert(::molresponse_v3::detail_kernel::has_apply_g_v<__VA_ARGS__>,                         \
                #__VA_ARGS__ " is missing `apply_g(w, g0, S1, S2, S3)`");                            \
  static_assert(::molresponse_v3::detail_kernel::has_compute_V0x_v<__VA_ARGS__>,                     \
                #__VA_ARGS__ " is missing `compute_V0x(w, g0, s)`");                                 \
  static_assert(::molresponse_v3::detail_kernel::has_compute_E0x_v<__VA_ARGS__>,                     \
                #__VA_ARGS__ " is missing `compute_E0x(w, g0, s)`");                                 \
  static_assert(::molresponse_v3::detail_kernel::has_bsh_apply_v<__VA_ARGS__>,                       \
                #__VA_ARGS__ " is missing `bsh_apply(w, g0, s, theta, omega)`");                     \
  static_assert(::molresponse_v3::detail_kernel::has_compute_residual_norm_v<__VA_ARGS__>,           \
                #__VA_ARGS__ " is missing `compute_residual_norm(w, s_old, s_new)`");                \
  static_assert(::molresponse_v3::detail_kernel::is_fd_kernel_v<__VA_ARGS__>,                        \
                #__VA_ARGS__ " does not satisfy FDKernel")

#define MV3_ASSERT_ES_KERNEL(...)                                                                    \
  MV3_ASSERT_FD_KERNEL(__VA_ARGS__);                                                                 \
  static_assert(::molresponse_v3::detail_kernel::has_compute_T0x_v<__VA_ARGS__>,                     \
                #__VA_ARGS__ " is missing `compute_T0x(w, g0, s)`");                                 \
  static_assert(::molresponse_v3::detail_kernel::has_compute_E0x_full_v<__VA_ARGS__>,                \
                #__VA_ARGS__ " is missing `compute_E0x_full(w, g0, s)`");                            \
  static_assert(::molresponse_v3::detail_kernel::is_es_kernel_v<__VA_ARGS__>,                        \
                #__VA_ARGS__ " does not satisfy ESKernel")

#endif // MOLRESPONSE_V3_KERNELS_KERNEL_INTERFACE_HPP
