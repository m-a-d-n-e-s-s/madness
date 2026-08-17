#ifndef MOLRESPONSE_V3_KERNELS_SOURCE_SPEC_HPP
#define MOLRESPONSE_V3_KERNELS_SOURCE_SPEC_HPP

// ===========================================================================
// source_spec — declarative assembly engine for the quadratic response
// sources (V^{BC} for beta/Raman, the residue (P,Q) for 2PA, and future
// members of the same family: double residues, resonant Raman).
//
// The unifying observation (paper §"One Right-Hand Side", eqs. tpa_fbar /
// tpa_zetabar / tpa_dbc / tpa_compact_P,Q): every term of every quadratic
// source is one of exactly TWO moves built from ONE object — the perturbed
// Fock operator
//
//     F = v + g'[pair density],        g'[Σ|bra><ket|] t
//           = J[rho]*t - c_xc * Σ_pairs K(bra, ket) t          (eq:tpa_fbar)
//
// The two moves are
//
//   (1) APPLY the operator to a target vector, optionally Q̂-projected:
//           ± Q̂ ( g'[..] t + v t )            e.g.  -Q̂ F^B x^C  in
//       eq:tpa_compact_P, or vbc.hpp's gzeta / fb terms;
//
//   (2) take its OCCUPIED MATRIX and transform a second target by it:
//           ± transform( s, <phi | g'[..] t + v t > )
//       i.e.  Σ_k s_k F_kp  — the property-Fock matrix term of
//       eq:tpa_compact_P/Q, or vbc.hpp's fphi term.
//
// What distinguishes V^{BC} from (P,Q) is pure DATA: which (bra,ket) legs
// build the pair density (γ^B vs γ^{B†} vs γ_L^{BC} vs D^{BC}), which target,
// which sign, where Q̂ sits. This header holds that data model (SourceEntry /
// SourceSpec) and the evaluator assemble_source. It contains NO physics: the
// numerics are the EXISTING primitives two_electron::apply_gamma_raw,
// matrix_inner and transform, called in the same order as the bespoke
// builders they replace (so a spec evaluates bit-identically, up to
// floating-point summation order, to the hand-written source it encodes).
//
// Field-to-equation dictionary (see SourceEntry):
//   pairs + coulomb_density   the pair density inside g'[..]  — γ^B legs for
//                             F^B, transposed legs for the daggered F̄^B
//                             (eq:tpa_fbar), the ζ-carrying γ_L^{BC} of the
//                             vbc gzeta term, or the D^{BC} legs (eq:tpa_dbc).
//                             The closed-shell factor 2 lives IN the density
//                             (v3 convention, kernels/vbc.hpp header note).
//   one_electron              the raw perturbation operator v inside
//                             F = v + g'[..] (eq:tpa_fbar). Left
//                             uninitialized for pure two-electron entries.
//   target                    the vector F acts on (x^C, y^C, or phi).
//   mode = Apply              move (1):  ± Q̂ (F target)      [project_Q]
//   mode = OccupiedMatrix     move (2):  ± transform(transform_target,
//                                              <phi|F target>)
//   transpose_matrix          occupied-block dagger: use <phi|F target>^T,
//                             i.e. the F̄_kp orientation of eq:tpa_compact_P.
//   sign                      the ± in front of the whole move.
// ===========================================================================

#include "tags.hpp"
#include "tda.hpp"           // ResponseGroundState
#include "two_electron.hpp"  // two_electron::{ExchangePair, apply_gamma_raw}

#include <madness/mra/mra.h>

#include <map>
#include <stdexcept>
#include <vector>

namespace molresponse_v3::source_spec {

using vecfuncT = std::vector<madness::real_function_3d>;

/// One (bra, ket) exchange leg of a pair density, by VALUE (a vecfuncT copy is
/// a shallow vector of shared-impl Function handles, so this is cheap and
/// lifetime-safe — unlike two_electron::ExchangePair, whose references only
/// suit immediate braced-init-list calls).
struct LegPair {
  vecfuncT bra;
  vecfuncT ket;
};

/// One term of a quadratic source, in the vocabulary of the compact form
/// (eq:tpa_compact_P/Q) — see the field dictionary in the header comment.
struct SourceEntry {
  enum class Mode {
    Apply,           ///< move (1): ± Q̂( g'[pairs] target + v*target )
    OccupiedMatrix,  ///< move (2): ± transform(transform_target,
                     ///<                      <phi| g'[pairs] target + v*target>)
  };

  /// Pair density: Coulomb side. J = coulop(coulomb_density) enters as
  /// J*target. Must carry the closed-shell factor 2 itself (v3 convention).
  /// Leave uninitialized (default Function) for a pure one-electron entry.
  madness::real_function_3d coulomb_density;

  /// Pair density: exchange side. Each leg contributes -c_xc*K(bra,ket)target.
  /// Transposing every leg (bra<->ket) daggers the kernel: g'[γ] -> g'[γ†]
  /// (the F vs F̄ distinction of eq:tpa_fbar).
  std::vector<LegPair> pairs;

  /// The vector the perturbed Fock operator acts on.
  vecfuncT target;

  /// Optional one-electron operator v of F = v + g'[..]; adds v*target
  /// before projection / matrix extraction. Uninitialized => none.
  madness::real_function_3d one_electron;

  /// Overall sign (the ± of the move); applied last.
  double sign = 1.0;

  Mode mode = Mode::Apply;

  /// Mode::OccupiedMatrix only: the vector transformed by the occupied
  /// matrix — the "Σ_k s_k F_kp" second target s of move (2).
  vecfuncT transform_target;

  /// Mode::OccupiedMatrix only: use the transposed occupied block
  /// <phi_k|F|phi_p> -> <phi_p|F|phi_k> (the occupied-space dagger).
  bool transpose_matrix = false;

  /// Mode::Apply only: project the result onto the virtual space (Q̂).
  bool project_Q = false;
};

/// One output channel (e.g. the X or Y block of V^{BC}, or P or Q of the 2PA
/// residue) = an ordered list of entries, summed with their signs.
struct SourceSpec {
  std::vector<SourceEntry> entries;
};

/// Factory for a move-(1) entry:  sign * Q̂?( g'[legs] target + v*target ).
/// Pass a default-constructed Function for `rho`/`v` to omit that piece.
inline SourceEntry
apply_entry(madness::real_function_3d rho, std::vector<LegPair> legs,
            vecfuncT target, double sign, bool project_Q,
            madness::real_function_3d v = {}) {
  SourceEntry e;
  e.coulomb_density = std::move(rho);
  e.pairs           = std::move(legs);
  e.target          = std::move(target);
  e.one_electron    = std::move(v);
  e.sign            = sign;
  e.mode            = SourceEntry::Mode::Apply;
  e.project_Q       = project_Q;
  return e;
}

/// Factory for a move-(2) entry:
///   sign * transform(s, <phi| g'[legs] target + v*target >^(T?)) —
/// the property-Fock matrix move Σ_k s_k F_kp (F̄_kp when transposed).
inline SourceEntry
occupied_matrix_entry(madness::real_function_3d rho, std::vector<LegPair> legs,
                      vecfuncT target, vecfuncT s, double sign,
                      bool transpose_matrix,
                      madness::real_function_3d v = {}) {
  SourceEntry e;
  e.coulomb_density  = std::move(rho);
  e.pairs            = std::move(legs);
  e.target           = std::move(target);
  e.one_electron     = std::move(v);
  e.sign             = sign;
  e.mode             = SourceEntry::Mode::OccupiedMatrix;
  e.transform_target = std::move(s);
  e.transpose_matrix = transpose_matrix;
  return e;
}

namespace detail {

/// F(target) = g'[pairs](target) + v*target for one entry, projection-free.
/// `jcache` maps a density's FunctionImpl pointer to its applied Coulomb
/// potential so a density shared between entries/channels (e.g. γ^B feeding
/// both the fb and fphi terms of vbc, in both X and Y) is convolved once —
/// exactly like the bespoke builders, which hoist J out of their term loops.
inline vecfuncT
eval_fock_action(madness::World &world, const ResponseGroundState &g0,
                 const SourceEntry &e,
                 std::map<const void *, madness::real_function_3d> &jcache) {
  using namespace madness;
  vecfuncT val;
  if (e.coulomb_density.is_initialized()) {
    auto key = static_cast<const void *>(e.coulomb_density.get_impl().get());
    auto it = jcache.find(key);
    if (it == jcache.end())
      it = jcache.emplace(key, apply(*g0.coulop, e.coulomb_density)).first;
    // runtime-sized view onto the entry's legs for apply_gamma_raw
    std::vector<two_electron::ExchangePair> legs;
    legs.reserve(e.pairs.size());
    for (const auto &p : e.pairs) legs.push_back({p.bra, p.ket});
    val = two_electron::apply_gamma_raw(world, it->second, e.target, legs,
                                        g0.c_xc, g0.lo);
  } else {
    if (!e.pairs.empty())
      throw std::runtime_error(
          "source_spec: entry has exchange pairs but no coulomb_density — "
          "a pair density must supply both sides (put the closed-shell "
          "factor 2 in the density).");
    val = zero_functions_compressed<double, 3>(
        world, static_cast<int>(e.target.size()));
  }
  if (e.one_electron.is_initialized()) {
    auto vt = mul(world, e.one_electron, e.target, true);
    gaxpy(world, 1.0, val, 1.0, vt);
  }
  return val;
}

} // namespace detail

/// Evaluate a set of channel specs into their source vectors (one vecfuncT
/// per channel, NOT truncated — the caller owns truncation points so a spec
/// path can reproduce its bespoke builder's truncation sequence exactly).
/// The Coulomb cache is shared across channels: a density Function reused
/// between entries (same underlying impl) is convolved once.
inline std::vector<vecfuncT>
assemble_source(madness::World &world, const ResponseGroundState &g0,
                const std::vector<SourceSpec> &channels) {
  using namespace madness;
  std::map<const void *, real_function_3d> jcache;
  std::vector<vecfuncT> out;
  out.reserve(channels.size());
  for (const auto &spec : channels) {
    vecfuncT acc;
    for (const auto &e : spec.entries) {
      vecfuncT val = detail::eval_fock_action(world, g0, e, jcache);
      vecfuncT term;
      if (e.mode == SourceEntry::Mode::Apply) {
        // move (1): ± Q̂ (F target)
        term = e.project_Q ? g0.Qa(val) : std::move(val);
      } else {
        // move (2): ± transform(s, <phi|F target>) — the property-Fock
        // matrix term Σ_k s_k F_kp; transpose_matrix selects the daggered
        // occupied block F̄_kp (eq:tpa_compact_P vs Q).
        auto m = matrix_inner(world, g0.amo, val);
        if (e.transpose_matrix) m = transpose(m);
        term = transform(world, e.transform_target, m, true);
      }
      if (acc.empty())
        acc = zero_functions_compressed<double, 3>(
            world, static_cast<int>(term.size()));
      gaxpy(world, 1.0, acc, e.sign, term);
    }
    out.push_back(std::move(acc));
  }
  return out;
}

} // namespace molresponse_v3::source_spec

#endif // MOLRESPONSE_V3_KERNELS_SOURCE_SPEC_HPP
