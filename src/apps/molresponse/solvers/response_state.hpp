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
#include <cstdio>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

namespace molresponse_v3 {

using namespace madness;

// ---- Restart-archive backend selection -----------------------------------
// Which physical file family holds a state archive at a basename:
//   Native = MADNESS parallel archive parts (<basename>.00000, ...)
//   Hdf5   = single <basename>.h5 blob (opt-in, MADNESS_HAS_HDF5 builds)
// The metadata layer records the backend PER ENTRY at save time (fd_states /
// vbc_states `backend` field, roots.json `backend`) and threads it back into
// load — so a stale twin left by a backend toggle can never shadow the state
// the metadata actually describes. AutoDetect (= legacy entry with no backend
// recorded) falls back to existence-based detection, decided on rank 0 and
// broadcast so every rank takes the same branch.
enum class IoBackend { AutoDetect, Native, Hdf5 };

inline const char *io_backend_tag(IoBackend b) {
  switch (b) {
    case IoBackend::Native: return "native";
    case IoBackend::Hdf5:   return "hdf5";
    default:                return "auto";
  }
}
inline IoBackend io_backend_from_tag(const std::string &tag) {
  if (tag == "hdf5")   return IoBackend::Hdf5;
  if (tag == "native") return IoBackend::Native;
  return IoBackend::AutoDetect;  // "", "auto", or anything legacy
}

/// True when archives written by `backend` can be reloaded by a different
/// process count. Native nio=1 parallel archives are np-locked (worlddc
/// assumes #writers == #readers — a mismatch silently corrupts the container);
/// the HDF5 blob is gathered to a single client at write time, so it is
/// np-portable by construction. AutoDetect (legacy, backend unknown) is
/// treated as np-locked — conservative.
inline bool io_backend_np_portable(IoBackend b) { return b == IoBackend::Hdf5; }

namespace detail_state_io {

/// Rank-0 decides whether the .h5 twin exists, then broadcasts — per-rank
/// exists() on a weak-coherence filesystem can take different branches on
/// different ranks (HDF5 collective vs native ParallelInputArchive), which
/// deadlocks on mismatched collectives. Legacy fallback only.
inline bool h5_twin_exists_collective(World &world, const std::string &filename) {
  int have = 0;
  if (world.rank() == 0)
    have = std::filesystem::exists(filename + ".h5") ? 1 : 0;
  world.gop.broadcast(have, 0);
  return have != 0;
}

/// Rank-0-only: drop native parallel-archive parts (<basename>.00000, ...) at
/// this basename. Parts are contiguous from 0, so stop at the first missing.
inline void remove_native_twin_rank0(const std::string &filename) {
  for (int p = 0;; ++p) {
    char part[16];
    std::snprintf(part, sizeof part, ".%5.5d", p);
    std::error_code ec;
    if (!std::filesystem::remove(filename + part, ec)) break;
  }
}

/// Rank-0-only: drop the .h5 twin (and a stranded atomic-write tmp).
inline void remove_h5_twin_rank0(const std::string &filename) {
  std::error_code ec;
  std::filesystem::remove(filename + ".h5", ec);
  std::filesystem::remove(filename + ".h5.tmp", ec);
}

/// Shared save plumbing for every Storage shape. `cb(ar)` runs the shape's
/// `ar & ...` ops identically on either backend. Writes the opt-in backend,
/// then removes the OTHER backend's twin at the same basename so a later
/// backend toggle can never resurrect stale state via auto-detect. Collective
/// (archive write + fence); returns the backend written so the metadata layer
/// can record it.
template <class StoreCb>
inline IoBackend save_state(World &world, const std::string &filename,
                            StoreCb &&cb) {
#ifdef MADNESS_HAS_HDF5
  if (hdf5_io_enabled()) {  // opt-in (deck response.hdf5 / env MADRESPONSE_IO_HDF5)
    save_parallel_archive_hdf5(world, filename + ".h5", 0, cb);  // atomic (tmp+rename)
    if (world.rank() == 0) remove_native_twin_rank0(filename);
    world.gop.fence();
    return IoBackend::Hdf5;
  }
#endif
  // Atomic native write (review io HIGH): the native default path used to
  // stream directly into <filename>.00000 in truncate mode, so a kill mid-save
  // left a truncated archive that BinaryFstreamInputArchive reads back SILENTLY
  // SHORT (its ifstream sets no exceptions() and nothing checks fail()). Write
  // to a .tmp basename, then rank-0 renames the part file(s) into place — same
  // tmp+rename contract the HDF5 branch already had. (nio is 1 here, so a
  // single <base>.tmp.00000; the glob-rename stays correct if nio ever grows.)
  {
    archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
        world, (filename + ".tmp").c_str(), 1);
    cb(ar);
  }
  world.gop.fence();
  if (world.rank() == 0) {
    namespace fs = std::filesystem;
    const fs::path fp(filename);
    const fs::path dir = fp.has_parent_path() ? fp.parent_path() : fs::path(".");
    const std::string base = fp.filename().string();
    const std::string tmp_prefix = base + ".tmp";   // parts are "<base>.tmp.NNNNN"
    std::error_code ec;
    for (const auto &e : fs::directory_iterator(dir, ec)) {
      if (ec) break;
      const std::string n = e.path().filename().string();
      if (n.rfind(tmp_prefix, 0) != 0) continue;    // not one of our tmp parts
      const std::string suffix = n.substr(tmp_prefix.size());  // ".NNNNN" (or "")
      fs::rename(e.path(), dir / (base + suffix));
    }
    remove_h5_twin_rank0(filename);
  }
  world.gop.fence();
  return IoBackend::Native;
}

/// Shared load plumbing. `backend` comes from the state's metadata entry
/// (io_backend_from_tag); AutoDetect = legacy entry → rank-0-decided
/// existence fallback. The requested-backend branch is rank-uniform by
/// construction (callers broadcast the tag), so a throw here is collective.
template <class LoadCb>
inline void load_state(World &world, const std::string &filename,
                       IoBackend backend, LoadCb &&cb) {
  [[maybe_unused]] bool use_h5 = false;
  if (backend == IoBackend::Hdf5) {
#ifdef MADNESS_HAS_HDF5
    use_h5 = true;
#else
    throw std::runtime_error(
        "response_state load: metadata records backend=hdf5 for " + filename +
        " but this build has no HDF5 support — configure with "
        "-DMADNESS_ENABLE_HDF5=ON or recompute the state");
#endif
  } else if (backend == IoBackend::AutoDetect) {
#ifdef MADNESS_HAS_HDF5
    use_h5 = h5_twin_exists_collective(world, filename);
#endif
  }
#ifdef MADNESS_HAS_HDF5
  if (use_h5) {
    load_parallel_archive_hdf5(world, filename + ".h5", cb);
    return;
  }
#endif
  archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(
      world, filename.c_str(), 1);
  cb(ar);
}

}  // namespace detail_state_io

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

  IoBackend save(World &world, const std::string &filename) const {
    return detail_state_io::save_state(world, filename, [&](auto &ar) {
      const std::size_t na = x_alpha.size();
      ar & na;
      for (const auto &f : x_alpha) ar & f;
    });
  }

  static ResponseStateX load(World &world, const std::string &filename,
                             IoBackend backend = IoBackend::AutoDetect) {
    ResponseStateX s;
    detail_state_io::load_state(world, filename, backend, [&](auto &ar) {
      std::size_t na;
      ar & na;
      s.x_alpha.resize(na);
      for (auto &f : s.x_alpha) ar & f;
    });
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

  IoBackend save(World &world, const std::string &filename) const {
    return detail_state_io::save_state(world, filename, [&](auto &ar) {
      const std::size_t na = x_alpha.size();
      const std::size_t nb = x_beta.size();
      ar & na & nb;
      for (const auto &f : x_alpha) ar & f;
      for (const auto &f : x_beta)  ar & f;
    });
  }

  static ResponseStateX load(World &world, const std::string &filename,
                             IoBackend backend = IoBackend::AutoDetect) {
    ResponseStateX s;
    detail_state_io::load_state(world, filename, backend, [&](auto &ar) {
      std::size_t na, nb;
      ar & na & nb;
      s.x_alpha.resize(na);
      s.x_beta.resize(nb);
      for (auto &f : s.x_alpha) ar & f;
      for (auto &f : s.x_beta)  ar & f;
    });
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

  IoBackend save(World &world, const std::string &filename) const {
    return detail_state_io::save_state(world, filename, [&](auto &ar) {
      const std::size_t na = x_alpha.size();
      ar & na;
      for (const auto &f : x_alpha) ar & f;
      for (const auto &f : y_alpha) ar & f;
    });
  }

  static ResponseStateXY load(World &world, const std::string &filename,
                              IoBackend backend = IoBackend::AutoDetect) {
    ResponseStateXY s;
    detail_state_io::load_state(world, filename, backend, [&](auto &ar) {
      std::size_t na;
      ar & na;
      s.x_alpha.resize(na);
      s.y_alpha.resize(na);
      for (auto &f : s.x_alpha) ar & f;
      for (auto &f : s.y_alpha) ar & f;
    });
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

  IoBackend save(World &world, const std::string &filename) const {
    return detail_state_io::save_state(world, filename, [&](auto &ar) {
      const std::size_t na = x_alpha.size();
      const std::size_t nb = x_beta.size();
      ar & na & nb;
      for (const auto &f : x_alpha) ar & f;
      for (const auto &f : y_alpha) ar & f;
      for (const auto &f : x_beta)  ar & f;
      for (const auto &f : y_beta)  ar & f;
    });
  }

  static ResponseStateXY load(World &world, const std::string &filename,
                              IoBackend backend = IoBackend::AutoDetect) {
    ResponseStateXY s;
    detail_state_io::load_state(world, filename, backend, [&](auto &ar) {
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
    });
    return s;
  }
};

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_RESPONSE_STATE_HPP
