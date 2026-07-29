#ifndef MOLRESPONSE_V3_KERNELS_TAGS_HPP
#define MOLRESPONSE_V3_KERNELS_TAGS_HPP

// =========================================================================
// Tag types + trait mappings for two axes:
//
//   ResponseType:  Static, Full, TDA
//                  — the math of the response (density formula, gamma
//                    formula, BSH shift). Static and TDA both use x-only
//                    storage but differ on density (y=x vs y=0).
//
//   Shell:         ClosedShell, OpenShell
//                  — whether beta channel exists. Closed-shell types
//                    do not even declare a x_beta / y_beta field; an
//                    OpenShell kernel cannot accidentally reach into a
//                    ClosedShell state's beta block because the field
//                    does not exist on that type.
//
// Six (ResponseType, Shell) pairs total; each gets its own Kernels<T,S>
// specialization with the math appropriate for that case. Storage is
// factored into two logical shapes (X-only and X-and-Y), each
// templated on Shell so the field set is shell-strict. Save/load lives
// on the storage shape — there are two "save/load logic" bodies
// (X-only / X-and-Y), each instantiated for both shells.
// =========================================================================

namespace molresponse_v3 {

/// Solver output verbosity. Shared by ESSolver<T,S>, FDSolver<T,S>,
/// and the iterate<> driver. Lives here so kernel/solver headers
/// don't need to include FDSolver.hpp just for the enum.
///
///   Silent   nothing printed.
///   Normal   per-iter banner + header + final summary.
///   Verbose  + per-root residual breakdown each iter.
///   Debug    + per-state component norms, A/S/U matrix dumps,
///            <X|V0|X> / <X|Lambda|X> / <X|Theta|X> inner-product
///            matrices each iter.
enum class PrintLevel { Silent = 0, Normal = 1, Verbose = 2, Debug = 3 };

// ---- Response-type tags ----
struct Static {};
struct Full   {};
struct TDA    {};

// ---- Shell tags ----
struct ClosedShell {};
struct OpenShell   {};

// ---- Forward declarations of the trait classes ----
// Concrete specializations live in kernels/<type>.hpp.
template <typename Type, typename Shell> struct Kernels;

// Forward-declare the two storage templates so kernel headers can use
// StorageOf_t before pulling in response_state.hpp.
template <typename Shell> struct ResponseStateX;
template <typename Shell> struct ResponseStateXY;

// Storage selection by (ResponseType, Shell):
//   Static, TDA  → ResponseStateX<Shell>     (x only)
//   Full         → ResponseStateXY<Shell>    (x and y)
template <typename Type, typename Shell> struct StorageOf;
template <typename S> struct StorageOf<Static, S> { using type = ResponseStateX<S>;  };
template <typename S> struct StorageOf<TDA,    S> { using type = ResponseStateX<S>;  };
template <typename S> struct StorageOf<Full,   S> { using type = ResponseStateXY<S>; };

template <typename Type, typename Shell>
using StorageOf_t = typename StorageOf<Type, Shell>::type;

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_TAGS_HPP
