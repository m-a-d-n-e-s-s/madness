#ifndef MOLRESPONSE_V3_KERNELS_RESPONSE_SPACE_OPS_HPP
#define MOLRESPONSE_V3_KERNELS_RESPONSE_SPACE_OPS_HPP

// =========================================================================
// Type/Shell-agnostic operations on a response_space.
//
// A response_space here is the same object as in the original
// molresponse codebase and Hurtado thesis: a collection of vectors of
// orbital functions, one per excited state (or per perturbation),
// arranged as an (n_states × n_orbitals) "matrix". Each row is one
// state's response orbitals.
//
//   using response_space = std::vector<vector_real_function_3d>;
//
// rs::inner mirrors the efficient form of molresponse's
// response_space_inner (src/apps/molresponse/response_functions.h):
// transpose to (n_orbitals × n_states) once, then accumulate
// matrix_inner per orbital index. For m states and n orbitals that's
// n batched matrix_inners of size m×m rather than m² scalar inners.
// =========================================================================

#include <madness/mra/mra.h>
#include <madness/tensor/tensor.h>
#include <madness/tensor/tensor_lapack.h>

#include <utility>
#include <algorithm>
#include <type_traits>
#include <vector>

namespace molresponse_v3 {

/// A response_space: n_states × n_orbitals matrix of 3D functions.
/// Each row is one state's response orbitals. Thesis / molresponse
/// vocabulary.
using response_space = std::vector<madness::vector_real_function_3d>;

namespace rs {

/// Inner-product matrix between two response_spaces.
///
///   result(i, j) = sum_p <a[i][p] | b[j][p]>
///
/// Efficient form (ported from molresponse::response_space_inner):
/// transpose each operand to (n_orbitals × n_states), then one
/// matrix_inner per orbital index. n_orbitals batched calls instead
/// of m² scalar inners.
inline madness::Tensor<double>
inner(const response_space &a, const response_space &b) {
  MADNESS_CHECK(!a.empty());
  MADNESS_CHECK(!b.empty());
  MADNESS_CHECK(a[0].size() == b[0].size());

  const long m_a   = static_cast<long>(a.size());
  const long m_b   = static_cast<long>(b.size());
  const long n_orb = static_cast<long>(a[0].size());

  madness::World &world = a[0][0].world();

  // matrix_inner requires compressed form.
  for (const auto &vi : a) compress(world, vi, false);
  for (const auto &vi : b) compress(world, vi, false);
  world.gop.fence();

  // Transposed views: aT[p][i] = a[i][p].
  // Function<double,3> has shared_ptr semantics; this aliases rather
  // than copies the underlying impls.
  std::vector<madness::vector_real_function_3d> aT(n_orb);
  std::vector<madness::vector_real_function_3d> bT(n_orb);
  for (long p = 0; p < n_orb; ++p) {
    aT[p].resize(m_a);
    bT[p].resize(m_b);
    for (long i = 0; i < m_a; ++i) aT[p][i] = a[i][p];
    for (long j = 0; j < m_b; ++j) bT[p][j] = b[j][p];
  }

  madness::Tensor<double> result(m_a, m_b);
  // matrix_inner is itself collective + self-synchronizing (it reduces internally),
  // and the aT/bT views are already compressed (fenced above) and not mutated in
  // the loop — so the per-orbital fence was redundant. Accumulate all n_orb
  // batched inners, fence once. Same arithmetic; drops n_orb-1 global barriers
  // per rs::inner call (called for A and S every ES iteration).
  for (long p = 0; p < n_orb; ++p)
    result += matrix_inner(world, aT[p], bT[p]);
  world.gop.fence();
  return result;
}

/// Output of `diagonalize` — includes evals + U plus diagnostic
/// data the caller can log to track slot identity across iters:
///
///   omega_ascending — eigenvalues in ascending order, indexed by the
///                     ω-rank (NOT by output slot).
///   dominance_perm  — `dominance_perm[i] = j` means output slot i was
///                     filled by the ω-rank-j eigenvector during the
///                     diagonal-dominance swap loop. Identity perm
///                     means the sygvp ascending order matched the
///                     input slot ordering — slots are "well-aligned".
///   U_diag          — diag(U) AFTER all fixups, indexed by output
///                     slot. Values near ±1 mean the rotation kept slot
///                     identity; values near 0 mean the slot picked up
///                     character from another slot via the rotation.
struct DiagonalizeResult {
  madness::Tensor<double> omega;          ///< eigenvalues, output-slot indexed
  madness::Tensor<double> U;              ///< rotation, output-slot columns
  madness::Tensor<double> omega_ascending; ///< pre-dominance ω order
  std::vector<long>       dominance_perm; ///< i ↦ source ω-rank of slot i
  madness::Tensor<double> U_diag;         ///< diag(U) post-fixups, per slot
};

/// Solve A U = S U diag(omega). Returns evals + U + diagnostic data
/// (DiagonalizeResult).
///
/// Beyond `sygvp`, this performs the same identity-preserving fixups
/// that the legacy molresponse::get_fock_transformation does — without
/// them, per-slot response densities oscillate near degeneracies (H2's
/// ω₂ ≈ ω₃ case):
///
///   0. Symmetrize: A := 0.5·(A + Aᵀ), S := 0.5·(S + Sᵀ). Numerical
///      noise (and the asymmetry in <X|Λ> when Λ isn't computed
///      self-adjointly) can push A off-symmetric by O(thresh); sygvp
///      assumes symmetric. Mirrors ExcitedResponse.cpp:835,846.
///   1. SVD of S — drop columns whose singular value < 10·thresh_degen
///      (subspace reduction; rare in practice, included for robustness).
///   2. Generalized eigenproblem on the (possibly reduced) A, S.
///   3. Diagonal-dominance sort: repeatedly swap pairs of columns so the
///      diagonal of U holds the largest entries. This keeps slot s
///      tracking the eigenvalue closest to slot s's input axis, so
///      iter-to-iter slot identity is continuous.
///   4. Phase fix: flip the sign of any column with U(i,i) < 0 so the
///      diagonal is non-negative — fixes the canonical phase choice.
///   5. Cluster unmixing: within each near-degenerate eigenvalue
///      cluster (|ω_i − ω_j| < cluster_factor · thresh_degenerate ·
///      max(|ω_i|, 1)), replace the cluster's U block with the polar-
///      decomposition orthogonal factor of that block — picks the
///      rotation closest to identity within the degenerate subspace,
///      so eigenvectors don't rotate freely between iterations.
///   6. Transform U back to the original size if step 1 reduced it.
///
/// `thresh_degenerate < 0` (default) → use FunctionDefaults<3>::thresh.
/// `cluster_factor`: TDA convention is 100 (loose), Full/RPA is 10
///   (tighter — matches legacy molresponse Full path, see
///   ExcitedResponse.cpp:1073 / TDDFT.cc:3097). The caller decides.
inline DiagonalizeResult
diagonalize(const madness::Tensor<double> &A_in,
            const madness::Tensor<double> &S_in,
            double thresh_degenerate = -1.0,
            double cluster_factor    = 100.0,
            madness::World *world_ptr = nullptr) {
  // world_ptr: the world whose ranks are calling this (sygvp broadcasts on
  // it). nullptr keeps the legacy World::get_default() — WRONG inside a
  // subworld (kernels review: cross-communicator broadcast is a hang/cross-
  // talk landmine). Subworld-capable callers MUST pass their world; new call
  // sites should always pass it explicitly.
  using madness::Slice;
  using madness::_;

  if (thresh_degenerate < 0.0) {
    thresh_degenerate = madness::FunctionDefaults<3>::get_thresh();
  }

  madness::Tensor<double> A = madness::copy(A_in);
  madness::Tensor<double> S = madness::copy(S_in);

  // ---- 0. Symmetrize (cheap, harmless if already symmetric) ------------
  A = 0.5 * (A + madness::transpose(A));
  S = 0.5 * (S + madness::transpose(S));

  // ---- 1. SVD of overlap; identify singular columns --------------------
  madness::Tensor<double> l_vecs, s_vals, r_vecs;
  {
    auto S_copy = madness::copy(S);
    madness::svd(S_copy, l_vecs, s_vals, r_vecs);
  }

  std::size_t num_sv = 0;
  for (long i = 0; i < s_vals.dim(0); ++i) {
    if (s_vals(i) < 10.0 * thresh_degenerate) ++num_sv;
  }
  const std::size_t size_l = static_cast<std::size_t>(s_vals.dim(0));
  const std::size_t size_s = size_l - num_sv;

  // ---- 2a. Reduce subspace if needed ------------------------------------
  if (num_sv > 0) {
    S = madness::Tensor<double>(size_s, size_s);
    for (std::size_t i = 0; i < size_s; ++i) S(i, i) = s_vals(i);

    auto l_vecs_s = madness::copy(l_vecs(_, Slice(0, static_cast<long>(size_s) - 1)));
    madness::Tensor<double> work(size_l, size_s);
    madness::mxm(size_l, size_s, size_l, work.ptr(), A.ptr(), l_vecs_s.ptr());
    A = madness::Tensor<double>(size_s, size_s);
    auto l_vecs_t = madness::transpose(l_vecs);
    madness::mxm(size_s, size_s, size_l, A.ptr(), l_vecs_t.ptr(), work.ptr());
  }

  // ---- 2b. Generalized eigenproblem -------------------------------------
  madness::Tensor<double> U_small, evals;
  madness::sygvp(world_ptr ? *world_ptr : madness::World::get_default(),
                 A, S, 1, U_small, evals);

  const long nmo = A.dim(0);

  // Snapshot ascending-order evals BEFORE the dominance swap so the
  // caller can compare "what sygvp returned" vs "what slots ended up
  // holding". omega_ascending(k) = k-th smallest ω.
  madness::Tensor<double> omega_ascending = madness::copy(evals);

  // Track the dominance permutation: perm[i] = which sygvp column
  // ended up in slot i after the swap loop. Identity perm means the
  // sygvp ascending order matched the input slot ordering.
  std::vector<long> perm(nmo);
  for (long i = 0; i < nmo; ++i) perm[i] = i;

  // ---- 3. Diagonal-dominance sort (preserve slot identity) -------------
  bool switched = true;
  while (switched) {
    switched = false;
    for (long i = 0; i < nmo; ++i) {
      for (long j = i + 1; j < nmo; ++j) {
        const double sold = U_small(i, i) * U_small(i, i)
                          + U_small(j, j) * U_small(j, j);
        const double snew = U_small(i, j) * U_small(i, j)
                          + U_small(j, i) * U_small(j, i);
        if (snew > sold) {
          auto tmp = madness::copy(U_small(_, i));
          U_small(_, i) = U_small(_, j);
          U_small(_, j) = tmp;
          std::swap(evals[i], evals[j]);
          std::swap(perm[i], perm[j]);
          switched = true;
        }
      }
    }
  }

  // ---- 4. Phase fix (column-wise sign so U(i,i) ≥ 0) -------------------
  for (long i = 0; i < nmo; ++i) {
    if (U_small(i, i) < 0.0) U_small(_, i).scale(-1.0);
  }

  // ---- 5. Polar-decomposition unmixing within degenerate clusters -----
  // cluster_factor=100 matches legacy TDA path; 10 matches legacy Full.
  long ilo = 0;
  while (ilo < nmo - 1) {
    long ihi = ilo;
    while (std::fabs(evals[ilo] - evals[ihi + 1]) <
           thresh_degenerate * cluster_factor *
               std::max(std::fabs(evals[ilo]), 1.0)) {
      ++ihi;
      if (ihi == nmo - 1) break;
    }
    const long nclus = ihi - ilo + 1;
    if (nclus > 1) {
      auto q = madness::copy(U_small(Slice(ilo, ihi), Slice(ilo, ihi)));
      madness::Tensor<double> W(nclus, nclus), VH(nclus, nclus), sigma(nclus);
      madness::svd(q, W, sigma, VH);
      q = madness::transpose(madness::inner(W, VH));
      U_small(_, Slice(ilo, ihi)) =
          madness::inner(U_small(_, Slice(ilo, ihi)), q);
    }
    ilo = ihi + 1;
  }

  // NOTE (no step 5.5): an earlier revision re-imposed ascending-omega
  // order here (port of legacy sort_eigenvalues). That re-sort exactly
  // INVERTS the dominance permutation of step 3 (perm collapses back to
  // identity, so the [ROT-SLOTS] REORDERED diagnostic could never fire),
  // voids the step-4 phase fix (a swapped pair can end with U(i,i) < 0
  // un-refixed), and rotates step 5's polar-unmix blocks off their slots.
  // Near-degenerate roots — whose cluster windows scale with thresh and
  // are widest at the coarse rung — then swap slots and flip signs every
  // iteration, poisoning the per-slot KAIN history. Slot order here is
  // DOMINANCE order (continuity with the caller's input slots); callers
  // wanting ascending-omega output canonicalize once at the end via
  // ESSolver::sort_state_by_omega (iterate_protocol's post-ramp hook).

  // ---- 6. Transform U back to original size if reduced -----------------
  madness::Tensor<double> U;
  if (num_sv > 0) {
    madness::Tensor<double> temp_U(size_l, size_l);
    temp_U(Slice(0, static_cast<long>(size_s) - 1),
           Slice(0, static_cast<long>(size_s) - 1)) = madness::copy(U_small);
    for (std::size_t i = size_s; i < size_l; ++i) temp_U(i, i) = 1.0;
    U = madness::Tensor<double>(size_l, size_l);
    madness::mxm(size_l, size_l, size_l,
                 U.ptr(), l_vecs.ptr(), temp_U.ptr());

    // U was expanded back to size_l but the eigenproblem ran at size_s, so
    // omega / omega_ascending / perm are SHORT. Pad them to size_l so the
    // caller can index omega(s) for every slot (previously out-of-bounds →
    // garbage BSH shift). The padded slots are the near-null-space
    // directions of S dropped in step 2a; they carry the largest retained
    // eigenvalue as a bounded placeholder (the next orthonormalize +
    // rotation re-seeds those slots), and perm maps them to themselves.
    {
      const double pad = evals(static_cast<long>(size_s) - 1);
      madness::Tensor<double> evals_full(static_cast<long>(size_l));
      madness::Tensor<double> asc_full(static_cast<long>(size_l));
      for (std::size_t i = 0; i < size_s; ++i) {
        evals_full(static_cast<long>(i)) = evals(static_cast<long>(i));
        asc_full(static_cast<long>(i))   = omega_ascending(static_cast<long>(i));
      }
      for (std::size_t i = size_s; i < size_l; ++i) {
        evals_full(static_cast<long>(i)) = pad;
        asc_full(static_cast<long>(i))   = pad;
        perm.push_back(static_cast<long>(i));
      }
      evals           = std::move(evals_full);
      omega_ascending = std::move(asc_full);
      if ((world_ptr ? *world_ptr : madness::World::get_default()).rank() == 0) {
        madness::print("[DIAG] WARNING: rank-deficient subspace overlap —",
                       num_sv, "of", size_l,
                       "slots dropped from the eigenproblem; their omega is a",
                       "placeholder copy of the top retained eigenvalue.");
      }
    }
  } else {
    U = std::move(U_small);
  }

  // Slot-wise diag(U) post-fixups, for diagnostic. After subspace
  // expansion (step 6) any slots beyond size_s are pure identity.
  madness::Tensor<double> U_diag(U.dim(0));
  for (long i = 0; i < U.dim(0); ++i) U_diag(i) = U(i, i);

  DiagonalizeResult out;
  out.omega           = std::move(evals);
  out.U               = std::move(U);
  out.omega_ascending = std::move(omega_ascending);
  out.dominance_perm  = std::move(perm);
  out.U_diag          = std::move(U_diag);
  return out;
}

/// Rotate a response_space in place: X_new[j] = sum_i U(i, j) * X_old[i].
///
/// Implementation: transpose to (n_orb × n_states), dispatch to the
/// native `madness::transform` (vmra.h) once per orbital, untranspose.
/// The native routine uses tensor-train sparsity via FunctionImpl's
/// vtransform pathway — much faster than n×m manual gaxpys on the
/// nested response_space layout, especially as n_orb grows.
inline void transform(madness::World &world, response_space &X,
                      const madness::Tensor<double> &U) {
  const long n = static_cast<long>(X.size());
  if (n == 0) return;
  const long n_orb = static_cast<long>(X[0].size());
  MADNESS_CHECK(U.dim(0) == n && U.dim(1) == n);

  // Transposed view: XT[p][i] = X[i][p].
  // Function<double,3> has shared_ptr semantics — this aliases impls
  // rather than copying tree data.
  std::vector<madness::vector_real_function_3d> XT(n_orb);
  for (long p = 0; p < n_orb; ++p) {
    XT[p].resize(n);
    for (long i = 0; i < n; ++i) XT[p][i] = X[i][p];
  }

  // Native transform per orbital: YT[p] = XT[p] · U (sum over states).
  std::vector<madness::vector_real_function_3d> YT(n_orb);
  for (long p = 0; p < n_orb; ++p) {
    YT[p] = madness::transform(world, XT[p], U, /*fence=*/false);
  }
  world.gop.fence();

  // Untranspose into Y[i][p] = YT[p][i].
  response_space Y(n);
  for (long i = 0; i < n; ++i) Y[i].resize(n_orb);
  for (long p = 0; p < n_orb; ++p)
    for (long i = 0; i < n; ++i)
      Y[i][p] = YT[p][i];

  const double thresh = madness::FunctionDefaults<3>::get_thresh();
  for (auto &v : Y) truncate(world, v, thresh);
  X = std::move(Y);
}

// ----------------------------------------------------------------------
// State-templated overloads. They flatten each State into the canonical
// component ordering (`flatten()` already does α (+β) concat per shell)
// then dispatch to the response_space versions above. The OpenShell
// bundle metric ⟨x_α|x_α⟩ + ⟨x_β|x_β⟩ and the slot-preserving rotation
// fall out of `flatten()` for free, so solvers no longer need their own
// pack/unpack scaffolding.
//
// Overload resolution: the non-template response_space variants win when
// the argument IS a response_space (`std::vector<vector_real_function_3d>`);
// the templates kick in only when the argument is a
// `std::vector<ResponseStateX<…>>` etc.
//
// Header-only — caller is responsible for having included the relevant
// `ResponseStateX/XY<…>` definitions before instantiation. No include
// of response_state.hpp here, to keep kernels/ free of a hard dependency
// on solvers/.
// ----------------------------------------------------------------------

template <typename State>
inline madness::Tensor<double>
inner(const std::vector<State> &a, const std::vector<State> &b) {
  response_space ra(a.size()), rb(b.size());
  for (std::size_t i = 0; i < a.size(); ++i) ra[i] = a[i].flatten();
  for (std::size_t j = 0; j < b.size(); ++j) rb[j] = b[j].flatten();
  return inner(ra, rb);
}

// Member-presence detection for the State storage blocks. Used by
// rs::metric to assemble the RPA overlap from whichever blocks a given
// (Type, Shell) State actually carries. void_t idiom — same pattern as
// kernel_interface.hpp's MV3_DEFINE_KERNEL_TRAIT.
namespace detail {
template <typename, typename = void> struct has_x_beta  : std::false_type {};
template <typename T> struct has_x_beta <T,
    std::void_t<decltype(std::declval<T>().x_beta )>> : std::true_type {};
template <typename, typename = void> struct has_y_alpha : std::false_type {};
template <typename T> struct has_y_alpha<T,
    std::void_t<decltype(std::declval<T>().y_alpha)>> : std::true_type {};
template <typename, typename = void> struct has_y_beta  : std::false_type {};
template <typename T> struct has_y_beta <T,
    std::void_t<decltype(std::declval<T>().y_beta )>> : std::true_type {};
}  // namespace detail

/// RPA subspace overlap S = ⟨a | M | b⟩ with the indefinite metric
/// M = diag(I_X, −I_Y):
///
///   S(i,j) = ⟨X_a_i | X_b_j⟩ − ⟨Y_a_i | Y_b_j⟩   (closed shell)
///
/// Assembled block-by-block from whichever members the State carries
/// (detected at compile time). The de-excitation (y_*) blocks enter
/// with a minus sign — that's the RPA metric (legacy
/// ExcitedResponse.cpp:871 S = ⟨X|X⟩ − ⟨Y|Y⟩). The excitation (x_*)
/// blocks enter with plus. So:
///   X-only (TDA): S = ⟨X|X⟩                       (= inner(a,b))
///   XY  (Full)  : S = ⟨X|X⟩ − ⟨Y|Y⟩
///   + β blocks for OpenShell, same signs per spin.
///
/// `if constexpr` discards the branch for any block a given State lacks
/// (the &State::y_beta etc. is dependent, so it's never instantiated for
/// types without that member). The sign lives on the scalar inner-product
/// matrix — no function copies, unlike a negate-the-Y-block approach.
///
/// The companion subspace matrix A = inner(roots, lambda) stays the
/// Euclidean (+) inner; only S carries the sign. As long as the states
/// are excitation-dominated (‖X‖ > ‖Y‖ per state — true for physical
/// singlets away from a triplet instability) S stays positive-definite
/// and the standard sygvp path in `diagonalize` applies unmodified.
template <typename State>
inline madness::Tensor<double>
metric(const std::vector<State> &a, const std::vector<State> &b) {
  // Pull one named block from every state into a response_space.
  auto block = [](const std::vector<State> &v,
                  madness::vector_real_function_3d State::*mem) {
    response_space rs(v.size());
    for (std::size_t i = 0; i < v.size(); ++i) rs[i] = v[i].*mem;
    return rs;
  };

  auto S = inner(block(a, &State::x_alpha), block(b, &State::x_alpha));
  if constexpr (detail::has_x_beta<State>::value)
    S += inner(block(a, &State::x_beta),  block(b, &State::x_beta));
  if constexpr (detail::has_y_alpha<State>::value)
    S -= inner(block(a, &State::y_alpha), block(b, &State::y_alpha));
  if constexpr (detail::has_y_beta<State>::value)
    S -= inner(block(a, &State::y_beta),  block(b, &State::y_beta));
  return S;
}

/// Scalar response-metric inner product of two single states — the
/// per-state form of `metric` above, and the convention used to
/// normalize response vectors (legacy ResponseBase.cpp::normalize):
///
///   x-only  state : ⟨a|b⟩ = ⟨Xa|Xb⟩                  (standard norm)
///   x-and-y state : ⟨a|b⟩ = ⟨Xa|Xb⟩ − ⟨Ya|Yb⟩        (indefinite RPA metric)
///
/// (+ β blocks for OpenShell, same signs per spin.) The de-excitation
/// (Y) block carries the minus sign; for an X-only state there is no Y
/// block so this is just the standard ⟨X|X⟩. Same indefinite metric as
/// the subspace overlap `metric`, so a bundle orthonormalized with this
/// goes into `diagonalize` with S ≈ I.
template <typename State>
inline double metric_inner(const State &a, const State &b) {
  double s = madness::inner(a.x_alpha, b.x_alpha);
  if constexpr (detail::has_x_beta<State>::value)
    s += madness::inner(a.x_beta,  b.x_beta);
  if constexpr (detail::has_y_alpha<State>::value)
    s -= madness::inner(a.y_alpha, b.y_alpha);
  if constexpr (detail::has_y_beta<State>::value)
    s -= madness::inner(a.y_beta,  b.y_beta);
  return s;
}

template <typename State>
inline void
transform(madness::World &world, std::vector<State> &X,
          const madness::Tensor<double> &U) {
  response_space rX(X.size());
  for (std::size_t i = 0; i < X.size(); ++i) rX[i] = X[i].flatten();
  transform(world, rX, U);
  for (std::size_t i = 0; i < X.size(); ++i) X[i].from_flat(rX[i]);
}

/// Project every state's response orbitals onto the ground-state virtual
/// space — i.e. strip occupied-orbital contamination, per spin. Alpha
/// blocks (x_alpha, y_alpha) go through `Qa`; beta blocks (x_beta,
/// y_beta) through `Qb`. Which blocks a State carries is detected at
/// compile time, so this is correct for X-only (TDA) and (X,Y) (Full)
/// states and for both shells; the beta branches are discarded for
/// closed shell, so `Qb` is never invoked there.
///
/// `Proj` is any callable mapping a vecfunc to its projected vecfunc
/// (e.g. madness::QProjector<double,3>). Templated so this header keeps
/// no chem dependency.
///
/// Used at the top of each ES iteration: BSH + KAIN leak a little
/// occupied-orbital character into the bundle, and without re-projection
/// that residual gets amplified by the next BSH apply into the
/// ω = −ε_core ghost eigenvector. (Q is idempotent, so re-projecting an
/// already-clean block — e.g. a Y block that bsh_apply just projected —
/// is a harmless no-op.) Mirrors legacy ExcitedResponse.cpp:2483-2517.
template <typename State, typename Proj>
inline void project(std::vector<State> &roots,
                    const Proj &Qa, const Proj &Qb) {
  for (auto &s : roots) {
    s.x_alpha = Qa(s.x_alpha);
    if constexpr (detail::has_y_alpha<State>::value) s.y_alpha = Qa(s.y_alpha);
    if constexpr (detail::has_x_beta <State>::value) s.x_beta  = Qb(s.x_beta);
    if constexpr (detail::has_y_beta <State>::value) s.y_beta  = Qb(s.y_beta);
  }
}

/// In-place modified Gram-Schmidt orthonormalization of a state bundle
/// in the RESPONSE metric (`metric_inner`): ⟨X|X⟩ for x-only states,
/// ⟨X|X⟩ − ⟨Y|Y⟩ for x-and-y states. This is the same indefinite metric
/// the subspace step uses, so the bundle enters `diagonalize` with
/// S ≈ I. Slot 0 is anchored and each higher slot is projected against
/// all lower slots, so the input ordering is preserved — important for
/// keeping per-slot KAIN history aligned.
///
/// Two passes by default: a single MGS pass leaves O(machine·κ) residual
/// overlap for an ill-conditioned bundle; the second cleans it up. A
/// slot whose metric norm collapses below 1e-12 (linearly dependent on
/// the lower slots) is left un-normalized rather than blown up.
template <typename State>
inline void orthonormalize(madness::World &world, std::vector<State> &roots,
                           int n_passes = 2) {
  const std::size_t M = roots.size();

  for (int pass = 0; pass < n_passes; ++pass) {
    for (std::size_t i = 0; i < M; ++i) {
      // Project out the lower (already metric-normalized) slots:
      //   |i⟩ -= ⟨j|i⟩ |j⟩,   inner products in the response metric.
      for (std::size_t j = 0; j < i; ++j) {
        const double c = metric_inner(roots[j], roots[i]);
        roots[i].axpy(world, -c, roots[j]);
      }
      // Normalize slot i in the response metric: ⟨X|X⟩ (x-only) or
      // ⟨X|X⟩ − ⟨Y|Y⟩ (x-and-y). abs() guards the (unphysical) negative-
      // norm case so a stray slot can't NaN the whole bundle.
      const double nrm = std::sqrt(std::abs(metric_inner(roots[i], roots[i])));
      if (nrm > 1.0e-12) roots[i].scale(world, 1.0 / nrm);
    }
  }
}

} // namespace rs
} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_KERNELS_RESPONSE_SPACE_OPS_HPP
