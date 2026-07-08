#ifndef MOLRESPONSE_V3_SOLVERS_ES_ROOT_IDENTITY_HPP
#define MOLRESPONSE_V3_SOLVERS_ES_ROOT_IDENTITY_HPP

// =========================================================================
// Stable identity for excited-state roots (Increment 9, doc 03).
//
// The problem: for an excited state the "frequency" omega is an eigenvalue
// — an OUTPUT that shifts every iteration and can cross neighbours between
// protocols. So neither omega nor the array slot is a stable identity. We
// assign each root a permanent `stable_index` at first appearance and carry
// it along through the (a) discrete slot permutations done by
// `sort_state_by_omega` and (b) save/load of the bundle. The bundle — all
// roots together — is the unit of persistence; per-root identity lives
// inside it as the `slot_permutation` (slot -> stable_index).
//
//   root_id        "es_root_0001"  — derived from stable_index, never changes
//   stable_index   permanent integer identity, assigned once
//   slot_index     current array position, CAN change on reorder
//   energy         current omega, metadata only
//   display_name   "S1", "S2_a", ... assigned at protocol boundaries by
//                  energy grouping (tol = 10 * thresh), suffix a/b/c for
//                  near-degenerate groups
//
// Identity assignment is a protocol-boundary / orchestration concern, NOT a
// per-iteration one: within a protocol step() mixes roots by continuous
// subspace rotation, so a clean permutation only exists at the discrete
// sort points. display_name grouping therefore happens at save / end-of-ramp.
// =========================================================================

#include <madness/tensor/tensor.h>

#include <cstdio>
#include <numeric>
#include <string>
#include <vector>

namespace molresponse_v3 {

/// "es_root_0001" from a stable_index.
inline std::string make_root_id(int stable_index) {
  char buf[32];
  std::snprintf(buf, sizeof buf, "es_root_%04d", stable_index);
  return std::string(buf);
}

/// Stable identity descriptor for one root within a bundle. This is the
/// human/metadata view; the solver itself only needs the `stable_index`
/// vector (slot -> stable_index) threaded through State.
struct ExcitedRootDescriptor {
  int         stable_index = -1;   // permanent identity
  int         slot_index   = -1;   // current array position (mutable)
  double      energy       = 0.0;  // current omega (metadata only)
  std::string display_name;        // assigned at protocol boundaries

  std::string root_id() const { return make_root_id(stable_index); }
};

/// Assign stable_index 0..M-1 to slots 0..M-1 (first appearance). No-op if
/// `stable_index` already has the right size — restart paths load it from
/// disk, and we must not clobber a loaded identity.
inline void assign_initial_stable_index(std::vector<int> &stable_index,
                                        int M) {
  if (static_cast<int>(stable_index.size()) == M) return;
  stable_index.resize(M);
  std::iota(stable_index.begin(), stable_index.end(), 0);
}

/// Apply a slot permutation `perm` (new slot i <- old slot perm[i]) to the
/// stable_index array, keeping each root's identity attached to its data as
/// `sort_state_by_omega` shuffles slots. No-op if `stable_index` is empty.
inline void permute_stable_index(std::vector<int> &stable_index,
                                 const std::vector<long> &perm) {
  if (stable_index.empty()) return;
  const long M = static_cast<long>(perm.size());
  std::vector<int> next(stable_index.size());
  for (long i = 0; i < M; ++i) {
    const long src = perm[i];
    next[i] = (src >= 0 && src < static_cast<long>(stable_index.size()))
                  ? stable_index[src]
                  : -1;
  }
  stable_index = std::move(next);
}

/// Group roots by energy (assumed ascending) with tolerance `tol` and assign
/// display names: S1, S2, ... where a group of size > 1 (near-degenerate
/// within tol) gets suffixes _a, _b, _c. Returns one name per slot.
inline std::vector<std::string>
assign_display_names(const madness::Tensor<double> &omega, double tol) {
  const long M = omega.dim(0);
  std::vector<std::string> names(static_cast<size_t>(M));
  if (M == 0) return names;

  // Partition consecutive slots into groups where the energy gap <= tol.
  std::vector<std::pair<long, long>> groups;  // [begin, end) per group
  long begin = 0;
  for (long i = 1; i < M; ++i) {
    if (omega(i) - omega(i - 1) > tol) {
      groups.emplace_back(begin, i);
      begin = i;
    }
  }
  groups.emplace_back(begin, M);

  int state_number = 1;
  for (const auto &g : groups) {
    const long n = g.second - g.first;
    for (long k = 0; k < n; ++k) {
      std::string name = "S" + std::to_string(state_number);
      if (n > 1) name += std::string("_") + static_cast<char>('a' + k);
      names[static_cast<size_t>(g.first + k)] = name;
    }
    ++state_number;
  }
  return names;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_ES_ROOT_IDENTITY_HPP
