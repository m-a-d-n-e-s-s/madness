#ifndef MOLRESPONSE_V3_SOLVERS_FUNCTION_HDF5_IO_HPP
#define MOLRESPONSE_V3_SOLVERS_FUNCTION_HDF5_IO_HPP

// =========================================================================
// function_hdf5_io.hpp — HDF5-backed Function I/O, Layer A (restart).
//
// TWO paths (compared in tests/test_function_hdf5.cpp):
//   * STRUCTURED  save/load_function_hdf5        — /meta + /keys + /coeffs
//                 datasets; interop-capable (MRChem/ParaView); double, NP=1.
//   * ARCHIVE     save/load_function_archive_hdf5 — opaque MADNESS-serialized
//                 blob in one dataset; universal + multi-rank; restart-only.
//
// --- STRUCTURED ---
// A structured (contiguous) HDF5 alternative to the native MADNESS archive
// (`madness::save/load` → ParallelOutputArchive<BinaryFstream...>). Stores the
// RAW coefficient tensors + per-node tree structure, mirroring what
// `FunctionImpl::store` serializes (`ar & k & thresh & ... & tree_state` then
// `ar & coeffs`) — so it round-trips reconstructed AND compressed trees.
//
// OPT-IN: compiled only when MADNESS_HAS_HDF5 is defined (the io-hdf5 thread's
// -DMADNESS_ENABLE_HDF5=ON path). The native archive remains the default; this
// header is a no-op in a normal build. See docs/30_io_hdf5_survey_and_plan.md.
//
// On-disk layout (<stem>.mad.h5):
//   /meta   attrs: schema, k, thresh, ndim, tree_state, initial_level,
//                  truncate_mode, autorefine, n_nodes, n_coeff_nodes,
//                  cell[NDIM][2]
//   /keys   int64  [n_nodes][3+NDIM] = (level, transl_0..transl_{NDIM-1},
//                                       has_children, has_coeff)
//   /coeffs double [n_coeff_nodes][k^NDIM]  contiguous; one row per coeff-bearing
//           node, same order as its /keys row (chunk only if adding partial
//           reads / compression — restart reads the whole function)
//
// SCOPE (P2 Layer A, first increment): T == double, single rank (NP=1). At NP=1
// the local container IS the whole tree, so no gather is needed. Multi-rank
// rank-0 gather and complex<T> are explicit follow-ups (guarded below).
// =========================================================================

#ifdef MADNESS_HAS_HDF5

#include <madness/mra/mra.h>
#include <madness/world/parallel_archive.h>
#include <madness/world/vector_archive.h>

#include <hdf5.h>

#include <algorithm>
#include <cstdlib>
#include <string>
#include <type_traits>
#include <vector>

namespace molresponse_v3 {
using namespace madness;

inline constexpr int kFunctionHdf5Schema = 1;

namespace detail_function_hdf5 {

// HDF5's default handler auto-prints the error stack to stderr on any failed
// call. On this build (1.14.3, threadsafety=no) that printer can recurse
// (H5E_printf_stack -> H5I_inc_ref -> H5E__push_stack -> ...) and overflow the
// stack, turning a benign/recoverable error into a SIGSEGV. Disable it around
// our I/O so failures surface as checked return codes, not a crash. RAII:
// restores the prior handler on scope exit.
struct H5ErrorScope {
  H5E_auto2_t old_func = nullptr;
  void* old_data = nullptr;
  H5ErrorScope() {
    H5Eget_auto2(H5E_DEFAULT, &old_func, &old_data);
    H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
  }
  ~H5ErrorScope() { H5Eset_auto2(H5E_DEFAULT, old_func, old_data); }
  H5ErrorScope(const H5ErrorScope&) = delete;
  H5ErrorScope& operator=(const H5ErrorScope&) = delete;
};

// --- thin core-HDF5-C-API helpers (no HL dependency) ---------------------
inline void write_attr(hid_t loc, const char* name, hid_t htype, const void* v) {
  hid_t sp = H5Screate(H5S_SCALAR);
  hid_t at = H5Acreate2(loc, name, htype, sp, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(at, htype, v);
  H5Aclose(at);
  H5Sclose(sp);
}
inline void write_attr_i(hid_t loc, const char* name, long long v) {
  write_attr(loc, name, H5T_NATIVE_LLONG, &v);
}
inline void write_attr_d(hid_t loc, const char* name, double v) {
  write_attr(loc, name, H5T_NATIVE_DOUBLE, &v);
}
inline void write_attr_darray(hid_t loc, const char* name, const double* v, hsize_t n) {
  hid_t sp = H5Screate_simple(1, &n, nullptr);
  hid_t at = H5Acreate2(loc, name, H5T_NATIVE_DOUBLE, sp, H5P_DEFAULT, H5P_DEFAULT);
  H5Awrite(at, H5T_NATIVE_DOUBLE, v);
  H5Aclose(at);
  H5Sclose(sp);
}
inline long long read_attr_i(hid_t loc, const char* name) {
  long long v = 0;
  hid_t at = H5Aopen(loc, name, H5P_DEFAULT);
  H5Aread(at, H5T_NATIVE_LLONG, &v);
  H5Aclose(at);
  return v;
}
inline double read_attr_d(hid_t loc, const char* name) {
  double v = 0.0;
  hid_t at = H5Aopen(loc, name, H5P_DEFAULT);
  H5Aread(at, H5T_NATIVE_DOUBLE, &v);
  H5Aclose(at);
  return v;
}
inline void read_attr_darray(hid_t loc, const char* name, double* out) {
  hid_t at = H5Aopen(loc, name, H5P_DEFAULT);
  H5Aread(at, H5T_NATIVE_DOUBLE, out);
  H5Aclose(at);
}
// 2-D dataset; chunk_by_row => chunk = [1, ncol] (one node's tensor per chunk).
inline void write_dataset_2d(hid_t file, const char* name, hid_t htype,
                             const void* data, hsize_t nrow, hsize_t ncol,
                             bool chunk_by_row) {
  hsize_t dims[2] = {nrow, ncol};
  hid_t sp = H5Screate_simple(2, dims, nullptr);
  hid_t dcpl = H5P_DEFAULT;
  hid_t plist = -1;
  if (chunk_by_row && nrow > 0) {
    plist = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk[2] = {1, ncol};
    H5Pset_chunk(plist, 2, chunk);
    dcpl = plist;
  }
  hid_t ds = H5Dcreate2(file, name, htype, sp, H5P_DEFAULT, dcpl, H5P_DEFAULT);
  if (nrow > 0) H5Dwrite(ds, htype, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
  H5Dclose(ds);
  if (plist >= 0) H5Pclose(plist);
  H5Sclose(sp);
}
inline void read_dataset(hid_t file, const char* name, hid_t htype, void* out) {
  hid_t ds = H5Dopen2(file, name, H5P_DEFAULT);
  H5Dread(ds, htype, H5S_ALL, H5S_ALL, H5P_DEFAULT, out);
  H5Dclose(ds);
}

}  // namespace detail_function_hdf5

// Save `f` to a single structured HDF5 file (raw coeffs + tree metadata).
template <typename T, std::size_t NDIM>
void save_function_hdf5(const Function<T, NDIM>& f, const std::string& path) {
  static_assert(std::is_same<T, double>::value,
                "function_hdf5_io Layer A is double-only for now (complex = follow-up)");
  namespace h = detail_function_hdf5;

  const auto impl = f.get_impl();
  World& world = impl->world;
  // First increment: single-rank only — local container == whole tree.
  MADNESS_CHECK(world.size() == 1);

  const auto& coeffs = impl->get_coeffs();
  const long k = impl->get_k();
  std::size_t npts = 1;
  for (std::size_t d = 0; d < NDIM; ++d) npts *= std::size_t(k);

  // --- collect every container node (faithful to `ar & coeffs`) ---
  std::vector<int64_t> keyrows;  // n_nodes * (3 + NDIM)
  std::vector<double> coeffbuf;  // n_coeff_nodes * npts (ragged: coeff nodes only)
  int64_t n_nodes = 0, n_coeff_nodes = 0;
  for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
    const Key<NDIM>& key = it->first;
    const auto& node = it->second;
    const bool hc = node.has_children();
    const bool hcoeff = node.has_coeff();
    keyrows.push_back(int64_t(key.level()));
    for (std::size_t d = 0; d < NDIM; ++d)
      keyrows.push_back(int64_t(key.translation()[d]));
    keyrows.push_back(hc ? 1 : 0);
    keyrows.push_back(hcoeff ? 1 : 0);
    ++n_nodes;
    if (hcoeff) {
      Tensor<T> t = node.coeff().full_tensor_copy();
      MADNESS_CHECK(std::size_t(t.size()) == npts);
      const T* p = t.ptr();
      coeffbuf.insert(coeffbuf.end(), p, p + npts);
      ++n_coeff_nodes;
    }
  }

  // --- write the file ---
  hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  MADNESS_CHECK(file >= 0);

  hid_t meta = H5Gcreate2(file, "meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  h::write_attr_i(meta, "schema", kFunctionHdf5Schema);
  h::write_attr_i(meta, "k", k);
  h::write_attr_d(meta, "thresh", impl->get_thresh());
  h::write_attr_i(meta, "ndim", int(NDIM));
  h::write_attr_i(meta, "tree_state", int(impl->get_tree_state()));
  h::write_attr_i(meta, "initial_level", impl->get_initial_level());
  h::write_attr_i(meta, "truncate_mode", impl->get_truncate_mode());
  h::write_attr_i(meta, "autorefine", impl->get_autorefine() ? 1 : 0);
  h::write_attr_i(meta, "n_nodes", n_nodes);
  h::write_attr_i(meta, "n_coeff_nodes", n_coeff_nodes);
  {
    const Tensor<double>& cell = FunctionDefaults<NDIM>::get_cell();  // [NDIM][2]
    h::write_attr_darray(meta, "cell", cell.ptr(), hsize_t(NDIM * 2));
  }
  H5Gclose(meta);

  h::write_dataset_2d(file, "keys", H5T_NATIVE_INT64, keyrows.data(),
                      hsize_t(n_nodes), hsize_t(3 + NDIM), /*chunk_by_row=*/false);
  // Contiguous: restart reads the whole function, so chunk-per-leaf is pure
  // per-call overhead (~2x slower I/O). Switch to chunked only if/when we add
  // per-leaf partial reads or a compression filter (gzip needs chunking).
  h::write_dataset_2d(file, "coeffs", H5T_NATIVE_DOUBLE, coeffbuf.data(),
                      hsize_t(n_coeff_nodes), hsize_t(npts), /*chunk_by_row=*/false);

  H5Fclose(file);
}

// Load a Function previously written by save_function_hdf5.
template <typename T, std::size_t NDIM>
Function<T, NDIM> load_function_hdf5(World& world, const std::string& path) {
  static_assert(std::is_same<T, double>::value,
                "function_hdf5_io Layer A is double-only for now (complex = follow-up)");
  namespace h = detail_function_hdf5;
  MADNESS_CHECK(world.size() == 1);

  hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  MADNESS_CHECK(file >= 0);

  hid_t meta = H5Gopen2(file, "meta", H5P_DEFAULT);
  const int schema = int(h::read_attr_i(meta, "schema"));
  MADNESS_CHECK(schema == kFunctionHdf5Schema);
  const long k = long(h::read_attr_i(meta, "k"));
  const double thresh = h::read_attr_d(meta, "thresh");
  const int ndim = int(h::read_attr_i(meta, "ndim"));
  MADNESS_CHECK(ndim == int(NDIM));
  const int tree_state = int(h::read_attr_i(meta, "tree_state"));
  const int64_t n_nodes = h::read_attr_i(meta, "n_nodes");
  const int64_t n_coeff_nodes = h::read_attr_i(meta, "n_coeff_nodes");
  double cellflat[NDIM * 2];
  h::read_attr_darray(meta, "cell", cellflat);
  H5Gclose(meta);

  std::size_t npts = 1;
  for (std::size_t d = 0; d < NDIM; ++d) npts *= std::size_t(k);

  std::vector<int64_t> keyrows(std::size_t(n_nodes) * (3 + NDIM));
  if (n_nodes > 0) h::read_dataset(file, "keys", H5T_NATIVE_INT64, keyrows.data());
  std::vector<double> coeffbuf(std::size_t(n_coeff_nodes) * npts);
  if (n_coeff_nodes > 0) h::read_dataset(file, "coeffs", H5T_NATIVE_DOUBLE, coeffbuf.data());
  H5Fclose(file);

  // --- rebuild the function ---
  Tensor<double> cell(NDIM, 2);
  std::copy(cellflat, cellflat + NDIM * 2, cell.ptr());
  FunctionDefaults<NDIM>::set_cell(cell);

  FunctionFactory<T, NDIM> factory(world);
  Function<T, NDIM> f(factory.k(int(k)).thresh(thresh).empty());
  world.gop.fence();

  auto& C = f.get_impl()->get_coeffs();
  long dims[NDIM];
  for (std::size_t d = 0; d < NDIM; ++d) dims[d] = k;

  std::size_t ci = 0;
  for (int64_t i = 0; i < n_nodes; ++i) {
    const std::size_t base = std::size_t(i) * (3 + NDIM);
    Level n = Level(keyrows[base]);
    Vector<Translation, NDIM> l;
    for (std::size_t d = 0; d < NDIM; ++d) l[d] = Translation(keyrows[base + 1 + d]);
    const bool hc = keyrows[base + 1 + NDIM] != 0;
    const bool hcoeff = keyrows[base + 2 + NDIM] != 0;
    Key<NDIM> key(n, l);
    if (hcoeff) {
      Tensor<T> t(static_cast<long>(NDIM), dims);
      std::copy(coeffbuf.begin() + ci * npts, coeffbuf.begin() + (ci + 1) * npts, t.ptr());
      ++ci;
      C.replace(key, FunctionNode<T, NDIM>(GenTensor<T>(t), hc));
    } else {
      C.replace(key, FunctionNode<T, NDIM>(GenTensor<T>(), hc));
    }
  }
  world.gop.fence();
  f.get_impl()->set_tree_state(TreeState(tree_state));
  f.verify_tree();
  return f;
}

// =========================================================================
// Archive-backend HDF5 (the "vector-archive" path).
//
// An opaque whole-Function blob produced by MADNESS's OWN parallel
// serialization, persisted to one HDF5 uint8 dataset. Reuses the optimized
// `ParallelOutputArchive<VectorOutputArchive>` gather (worlddc.h: thread-
// parallel serialize + MPI_Gatherv to rank 0) — so it is:
//   * universal — any tree state, GenTensor, ANY serializable T (incl. complex);
//   * multi-rank — works at nio=1 for NP>1 (rank 0 gathers/distributes), unlike
//     the structured path above (NP=1 only);
//   * a near drop-in for madness::save/load with the byte sink swapped to HDF5.
// Trade-off: OPAQUE (MADNESS bytes in an HDF5 container) — for restart/
// checkpoint, NOT interchange. Use the structured path for MRChem/plotting.
// =========================================================================

namespace detail_function_hdf5 {
// rank-0 byte-blob write: one uint8 dataset "madness_archive" (+ tiny /meta).
// deflate_level 0 = contiguous; 1..9 = chunked + gzip (transparent on read).
inline void write_byte_dataset(const std::string& path,
                               const std::vector<unsigned char>& buf,
                               int deflate_level) {
  H5ErrorScope h5errs;  // no recursive auto-print; codes checked below
  hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  MADNESS_CHECK_THROW(file >= 0, "write_byte_dataset: H5Fcreate failed (path/permissions/disk?)");
  hid_t meta = H5Gcreate2(file, "meta", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  MADNESS_CHECK_THROW(meta >= 0, "write_byte_dataset: H5Gcreate2(/meta) failed");
  write_attr_i(meta, "schema", kFunctionHdf5Schema);
  write_attr_i(meta, "archive_format", 1);
  write_attr_i(meta, "deflate", deflate_level);
  H5Gclose(meta);
  hsize_t n = buf.size();
  hid_t sp = H5Screate_simple(1, &n, nullptr);
  MADNESS_CHECK_THROW(sp >= 0, "write_byte_dataset: H5Screate_simple failed");
  hid_t dcpl = H5P_DEFAULT, plist = -1;
  if (deflate_level > 0 && n > 0) {  // gzip needs chunked storage
    plist = H5Pcreate(H5P_DATASET_CREATE);
    hsize_t chunk = std::min<hsize_t>(n, hsize_t(1) << 20);  // ~1 MiB chunks
    H5Pset_chunk(plist, 1, &chunk);
    H5Pset_deflate(plist, unsigned(deflate_level));
    dcpl = plist;
  }
  hid_t ds = H5Dcreate2(file, "madness_archive", H5T_NATIVE_UCHAR, sp, H5P_DEFAULT,
                        dcpl, H5P_DEFAULT);
  MADNESS_CHECK_THROW(ds >= 0, "write_byte_dataset: H5Dcreate2(madness_archive) failed");
  if (n) {
    herr_t st = H5Dwrite(ds, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    MADNESS_CHECK_THROW(st >= 0, "write_byte_dataset: H5Dwrite failed");
  }
  H5Dclose(ds);
  if (plist >= 0) H5Pclose(plist);
  H5Sclose(sp);
  H5Fclose(file);
}
inline void read_byte_dataset(const std::string& path,
                              std::vector<unsigned char>& buf) {
  H5ErrorScope h5errs;
  hid_t file = H5Fopen(path.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  MADNESS_CHECK_THROW(file >= 0, "read_byte_dataset: H5Fopen failed");
  hid_t ds = H5Dopen2(file, "madness_archive", H5P_DEFAULT);
  MADNESS_CHECK_THROW(ds >= 0, "read_byte_dataset: H5Dopen2(madness_archive) failed");
  hid_t sp = H5Dget_space(ds);
  hssize_t n = H5Sget_simple_extent_npoints(sp);
  buf.resize(n > 0 ? std::size_t(n) : 0);
  if (n > 0) {
    herr_t st = H5Dread(ds, H5T_NATIVE_UCHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, buf.data());
    MADNESS_CHECK_THROW(st >= 0, "read_byte_dataset: H5Dread failed");
  }
  H5Sclose(sp);
  H5Dclose(ds);
  H5Fclose(file);
}
}  // namespace detail_function_hdf5

/// Programmatic override for the HDF5 opt-in (deck parameter `response.hdf5`
/// or the standalone --hdf5 flag). -1 = unset (fall back to the env var).
inline int& hdf5_io_override() {
  static int v = -1;
  return v;
}
inline void set_hdf5_io_enabled(bool on) { hdf5_io_override() = on ? 1 : 0; }

/// Opt-in switch for HDF5-backed restart I/O: deck/flag override first, else
/// env MADRESPONSE_IO_HDF5 set and not "0". Controls only what is WRITTEN;
/// readers auto-detect a `.h5` file.
inline bool hdf5_io_enabled() {
  if (hdf5_io_override() >= 0) return hdf5_io_override() == 1;
  const char* e = std::getenv("MADRESPONSE_IO_HDF5");
  return e && e[0] && e[0] != '0';
}

/// Generic: serialize arbitrary parallel-archive content to one HDF5 dataset via
/// the optimized VectorOutputArchive gather. The callback runs the `ar & ...`
/// ops exactly as the BinaryFstream path would (e.g. `ar & na; for(f) ar & f;`).
template <class StoreCb>
void save_parallel_archive_hdf5(World& world, const std::string& path,
                                int deflate_level, StoreCb&& cb) {
  std::vector<unsigned char> buf;
  {
    archive::VectorOutputArchive var(buf);
    archive::ParallelOutputArchive<archive::VectorOutputArchive> par(world, var, 1);
    cb(par);  // 2067: thread-parallel serialize + MPI_Gatherv to rank 0's buf
    par.flush();
  }
  if (world.rank() == 0)
    detail_function_hdf5::write_byte_dataset(path, buf, deflate_level);
  world.gop.fence();
}

/// Generic inverse of save_parallel_archive_hdf5 (rank 0 reads; the callback's
/// `ar & ...` distributes — header broadcast, containers to owners).
template <class LoadCb>
void load_parallel_archive_hdf5(World& world, const std::string& path, LoadCb&& cb) {
  std::vector<unsigned char> buf;  // only the io-node (rank 0) needs the bytes
  if (world.rank() == 0) detail_function_hdf5::read_byte_dataset(path, buf);
  archive::VectorInputArchive var(buf);
  archive::ParallelInputArchive<archive::VectorInputArchive> par(world, var, 1);
  cb(par);  // 2328
  world.gop.fence();
}

// Save f as a single HDF5 dataset of MADNESS parallel-archive bytes.
// deflate_level 0 = uncompressed; 1..9 = gzip (opt-in size lever).
template <typename T, std::size_t NDIM>
void save_function_archive_hdf5(const Function<T, NDIM>& f, const std::string& path,
                                int deflate_level = 0) {
  save_parallel_archive_hdf5(f.world(), path, deflate_level,
                             [&](auto& ar) { ar & f; });
}

// Load a Function previously written by save_function_archive_hdf5.
template <typename T, std::size_t NDIM>
Function<T, NDIM> load_function_archive_hdf5(World& world, const std::string& path) {
  Function<T, NDIM> f;  // null impl OK: parallel Function load takes world from ar
  load_parallel_archive_hdf5(world, path, [&](auto& ar) { ar & f; });
  return f;
}

}  // namespace molresponse_v3

#endif  // MADNESS_HAS_HDF5
#endif  // MOLRESPONSE_V3_SOLVERS_FUNCTION_HDF5_IO_HPP
