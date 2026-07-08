#ifndef MOLRESPONSE_V3_SOLVERS_RESPONSE_STATE_HPP
#define MOLRESPONSE_V3_SOLVERS_RESPONSE_STATE_HPP

// =========================================================================
// Storage shapes — two templates, each specialized on Shell.
//
//   ResponseStateX<ClosedShell>   : { x_alpha }
//   ResponseStateX<OpenShell>     : { x_alpha, x_beta }
//   ResponseStateXY<ClosedShell>  : { x_alpha, y_alpha }
//   ResponseStateXY<OpenShell>    : { x_alpha, y_alpha, x_beta, y_beta }
//
// Four distinct types via two structures × two shells. The TYPE
// determines exactly which fields are accessible — no empty-vector
// guard, no "is this populated" check; a closed-shell kernel cannot
// even spell `state.x_beta` because that member does not exist on
// `ResponseStateX<ClosedShell>`.
//
// The mechanical per-block operations (axpy, copy, truncate, flatten,
// from_flat, flat_size) are IDENTICAL modulo "which blocks exist", so
// each struct exposes a single `blocks()` accessor (pointers to its
// member vecfuncs, in flatten order) and delegates the operations to
// the shared `state_ops::` helpers below. Adding or changing one of
// these operations is then a one-place edit, and a new block layout
// only has to declare its `blocks()` list. The structs stay plain
// aggregates (no base class) so the `State{...}` brace-init used all
// over the kernels keeps working.
//
// Save / load stay hand-written per shape: the on-disk format is
// shape-specific (and predates this header), so it's kept explicit
// rather than routed through blocks().
// =========================================================================

#include "../kernels/tags.hpp"
#include "function_hdf5_io.hpp"  // opt-in HDF5 restart (no-op unless MADNESS_HAS_HDF5)

#include <madness/mra/mra.h>
#include <madness/world/parallel_archive.h>

#include <array>
#include <cstddef>
#include <filesystem>
#include <string>
#include <vector>

namespace molresponse_v3 {

using namespace madness;

// ---- Shared per-block operations ----------------------------------------
// Generic over any State that exposes:
//   blocks()        -> std::array<std::vector<real_function_3d>*, N>
//   blocks() const  -> std::array<const std::vector<real_function_3d>*, N>
// with the blocks listed in the canonical flatten order. The flatten/
// from_flat round-trip preserves all storage; that order is stable per
// State type so KAIN's stored iterates index the same components every
// iteration.
namespace state_ops {

template <class S>
inline void axpy(S &s, World &world, double alpha, const S &other) {
  auto mine = s.blocks();
  auto theirs = other.blocks();
  for (std::size_t b = 0; b < mine.size(); ++b)
    gaxpy(world, 1.0, *mine[b], alpha, *theirs[b]);
}

/// Deep copy of every block. The default copy shares
/// shared_ptr<FunctionImpl> handles; callers that mutate in place
/// (assemble_*, axpy) need an independent buffer.
template <class S>
inline S copy(const S &s, World &world) {
  S out = s;  // shallow — shares handles, fixed up block-by-block below
  for (auto *blk : out.blocks()) *blk = madness::copy(world, *blk);
  return out;
}

template <class S>
inline void truncate_all(S &s, World &world, double thresh) {
  for (auto *blk : s.blocks()) madness::truncate(world, *blk, thresh);
}

template <class S>
inline void scale(S &s, World &world, double factor) {
  for (auto *blk : s.blocks()) madness::scale(world, *blk, factor);
}

template <class S>
inline std::size_t flat_size(const S &s) {
  std::size_t n = 0;
  for (const auto *blk : s.blocks()) n += blk->size();
  return n;
}

template <class S>
inline std::vector<real_function_3d> flatten(const S &s) {
  std::vector<real_function_3d> v;
  v.reserve(flat_size(s));
  for (const auto *blk : s.blocks())
    v.insert(v.end(), blk->begin(), blk->end());
  return v;
}

/// Unpack a flat vector back into the blocks. Each block must already
/// be sized (this reads blk->size() to partition) — true on every call
/// site, where from_flat follows a flatten of an already-shaped state.
template <class S>
inline void from_flat(S &s, const std::vector<real_function_3d> &v) {
  std::size_t off = 0;
  for (auto *blk : s.blocks()) {
    const std::size_t n = blk->size();
    *blk = std::vector<real_function_3d>(v.begin() + off, v.begin() + off + n);
    off += n;
  }
}

}  // namespace state_ops

// ---- ResponseStateX<ClosedShell> ----------------------------------------
template <>
struct ResponseStateX<ClosedShell> {
  std::vector<real_function_3d> x_alpha;

  std::size_t num_alpha() const { return x_alpha.size(); }

  // Blocks in flatten order. Sole per-type customization point — every
  // mechanical op below is generic over this list (see state_ops::).
  std::array<std::vector<real_function_3d>*, 1> blocks() {
    return {&x_alpha};
  }
  std::array<const std::vector<real_function_3d>*, 1> blocks() const {
    return {&x_alpha};
  }

  void axpy(World &world, double alpha,
            const ResponseStateX<ClosedShell> &other) {
    state_ops::axpy(*this, world, alpha, other);
  }
  ResponseStateX<ClosedShell> copy(World &world) const {
    return state_ops::copy(*this, world);
  }
  void truncate_all(World &world, double thresh) {
    state_ops::truncate_all(*this, world, thresh);
  }
  void scale(World &world, double factor) {
    state_ops::scale(*this, world, factor);
  }
  std::size_t flat_size() const { return state_ops::flat_size(*this); }
  std::vector<real_function_3d> flatten() const {
    return state_ops::flatten(*this);
  }
  void from_flat(const std::vector<real_function_3d>& v) {
    state_ops::from_flat(*this, v);
  }

  void save(World &world, const std::string &filename) const {
#ifdef MADNESS_HAS_HDF5
    if (hdf5_io_enabled()) {  // opt-in (env MADRESPONSE_IO_HDF5); writes <file>.h5
      save_parallel_archive_hdf5(world, filename + ".h5", 0, [&](auto &ar) {
        const std::size_t na = x_alpha.size();
        ar & na;
        for (const auto &f : x_alpha) ar & f;
      });
      return;
    }
#endif
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
        world, filename.c_str(), 1);
    const std::size_t na = x_alpha.size();
    ar & na;
    for (const auto &f : x_alpha) ar & f;
  }

  static ResponseStateX load(World &world, const std::string &filename) {
    ResponseStateX s;
#ifdef MADNESS_HAS_HDF5
    if (std::filesystem::exists(filename + ".h5")) {  // auto-detect HDF5 checkpoint
      load_parallel_archive_hdf5(world, filename + ".h5", [&](auto &ar) {
        std::size_t na;
        ar & na;
        s.x_alpha.resize(na);
        for (auto &f : s.x_alpha) ar & f;
      });
      return s;
    }
#endif
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
        world, filename.c_str(), 1);
    std::size_t na;
    ar & na;
    s.x_alpha.resize(na);
    for (auto &f : s.x_alpha) ar & f;
    return s;
  }
};

// ---- ResponseStateX<OpenShell> ------------------------------------------
template <>
struct ResponseStateX<OpenShell> {
  std::vector<real_function_3d> x_alpha;
  std::vector<real_function_3d> x_beta;

  std::size_t num_alpha() const { return x_alpha.size(); }
  std::size_t num_beta()  const { return x_beta.size(); }

  std::array<std::vector<real_function_3d>*, 2> blocks() {
    return {&x_alpha, &x_beta};
  }
  std::array<const std::vector<real_function_3d>*, 2> blocks() const {
    return {&x_alpha, &x_beta};
  }

  void axpy(World &world, double alpha,
            const ResponseStateX<OpenShell> &other) {
    state_ops::axpy(*this, world, alpha, other);
  }
  ResponseStateX<OpenShell> copy(World &world) const {
    return state_ops::copy(*this, world);
  }
  void truncate_all(World &world, double thresh) {
    state_ops::truncate_all(*this, world, thresh);
  }
  void scale(World &world, double factor) {
    state_ops::scale(*this, world, factor);
  }
  std::size_t flat_size() const { return state_ops::flat_size(*this); }
  std::vector<real_function_3d> flatten() const {
    return state_ops::flatten(*this);
  }
  void from_flat(const std::vector<real_function_3d>& v) {
    state_ops::from_flat(*this, v);
  }

  void save(World &world, const std::string &filename) const {
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
        world, filename.c_str(), 1);
    const std::size_t na = x_alpha.size();
    const std::size_t nb = x_beta.size();
    ar & na & nb;
    for (const auto &f : x_alpha) ar & f;
    for (const auto &f : x_beta)  ar & f;
  }

  static ResponseStateX load(World &world, const std::string &filename) {
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
        world, filename.c_str(), 1);
    ResponseStateX s;
    std::size_t na, nb;
    ar & na & nb;
    s.x_alpha.resize(na);
    s.x_beta.resize(nb);
    for (auto &f : s.x_alpha) ar & f;
    for (auto &f : s.x_beta)  ar & f;
    return s;
  }
};

// ---- ResponseStateXY<ClosedShell> ---------------------------------------
template <>
struct ResponseStateXY<ClosedShell> {
  std::vector<real_function_3d> x_alpha;
  std::vector<real_function_3d> y_alpha;

  std::size_t num_alpha() const { return x_alpha.size(); }

  std::array<std::vector<real_function_3d>*, 2> blocks() {
    return {&x_alpha, &y_alpha};
  }
  std::array<const std::vector<real_function_3d>*, 2> blocks() const {
    return {&x_alpha, &y_alpha};
  }

  void axpy(World &world, double alpha,
            const ResponseStateXY<ClosedShell> &other) {
    state_ops::axpy(*this, world, alpha, other);
  }
  ResponseStateXY<ClosedShell> copy(World &world) const {
    return state_ops::copy(*this, world);
  }
  void truncate_all(World &world, double thresh) {
    state_ops::truncate_all(*this, world, thresh);
  }
  void scale(World &world, double factor) {
    state_ops::scale(*this, world, factor);
  }
  std::size_t flat_size() const { return state_ops::flat_size(*this); }
  std::vector<real_function_3d> flatten() const {
    return state_ops::flatten(*this);
  }
  void from_flat(const std::vector<real_function_3d>& v) {
    state_ops::from_flat(*this, v);
  }

  void save(World &world, const std::string &filename) const {
#ifdef MADNESS_HAS_HDF5
    if (hdf5_io_enabled()) {  // opt-in (env MADRESPONSE_IO_HDF5); writes <file>.h5
      save_parallel_archive_hdf5(world, filename + ".h5", 0, [&](auto &ar) {
        const std::size_t na = x_alpha.size();
        ar & na;
        for (const auto &f : x_alpha) ar & f;
        for (const auto &f : y_alpha) ar & f;
      });
      return;
    }
#endif
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
        world, filename.c_str(), 1);
    const std::size_t na = x_alpha.size();
    ar & na;
    for (const auto &f : x_alpha) ar & f;
    for (const auto &f : y_alpha) ar & f;
  }

  static ResponseStateXY load(World &world, const std::string &filename) {
    ResponseStateXY s;
#ifdef MADNESS_HAS_HDF5
    if (std::filesystem::exists(filename + ".h5")) {  // auto-detect HDF5 checkpoint
      load_parallel_archive_hdf5(world, filename + ".h5", [&](auto &ar) {
        std::size_t na;
        ar & na;
        s.x_alpha.resize(na);
        s.y_alpha.resize(na);
        for (auto &f : s.x_alpha) ar & f;
        for (auto &f : s.y_alpha) ar & f;
      });
      return s;
    }
#endif
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
        world, filename.c_str(), 1);
    std::size_t na;
    ar & na;
    s.x_alpha.resize(na);
    s.y_alpha.resize(na);
    for (auto &f : s.x_alpha) ar & f;
    for (auto &f : s.y_alpha) ar & f;
    return s;
  }
};

// ---- ResponseStateXY<OpenShell> -----------------------------------------
template <>
struct ResponseStateXY<OpenShell> {
  std::vector<real_function_3d> x_alpha;
  std::vector<real_function_3d> y_alpha;
  std::vector<real_function_3d> x_beta;
  std::vector<real_function_3d> y_beta;

  std::size_t num_alpha() const { return x_alpha.size(); }
  std::size_t num_beta()  const { return x_beta.size(); }

  std::array<std::vector<real_function_3d>*, 4> blocks() {
    return {&x_alpha, &y_alpha, &x_beta, &y_beta};
  }
  std::array<const std::vector<real_function_3d>*, 4> blocks() const {
    return {&x_alpha, &y_alpha, &x_beta, &y_beta};
  }

  void axpy(World &world, double alpha,
            const ResponseStateXY<OpenShell> &other) {
    state_ops::axpy(*this, world, alpha, other);
  }
  ResponseStateXY<OpenShell> copy(World &world) const {
    return state_ops::copy(*this, world);
  }
  void truncate_all(World &world, double thresh) {
    state_ops::truncate_all(*this, world, thresh);
  }
  void scale(World &world, double factor) {
    state_ops::scale(*this, world, factor);
  }
  std::size_t flat_size() const { return state_ops::flat_size(*this); }
  std::vector<real_function_3d> flatten() const {
    return state_ops::flatten(*this);
  }
  void from_flat(const std::vector<real_function_3d>& v) {
    state_ops::from_flat(*this, v);
  }

  void save(World &world, const std::string &filename) const {
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
        world, filename.c_str(), 1);
    const std::size_t na = x_alpha.size();
    const std::size_t nb = x_beta.size();
    ar & na & nb;
    for (const auto &f : x_alpha) ar & f;
    for (const auto &f : y_alpha) ar & f;
    for (const auto &f : x_beta)  ar & f;
    for (const auto &f : y_beta)  ar & f;
  }

  static ResponseStateXY load(World &world, const std::string &filename) {
    archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
        world, filename.c_str(), 1);
    ResponseStateXY s;
    std::size_t na, nb;
    ar & na & nb;
    s.x_alpha.resize(na);
    s.y_alpha.resize(na);
    s.x_beta.resize(nb);
    s.y_beta.resize(nb);
    for (auto &f : s.x_alpha) ar & f;
    for (auto &f : s.y_alpha) ar & f;
    for (auto &f : s.x_beta)  ar & f;
    for (auto &f : s.y_beta)  ar & f;
    return s;
  }
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_RESPONSE_STATE_HPP
