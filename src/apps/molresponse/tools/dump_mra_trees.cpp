// Offline MRA tree-skeleton dumper  (molresponse_v3 — R2/L4 export, structural target)
//
// Loads already-converged ground-state orbitals from a moldft restart archive
// and dumps each orbital's adaptive multiresolution octree as JSON, in TWO
// representations:
//
//   mo_<i>_recon.json  reconstructed: leaves hold scaling coefficients, so the
//                      per-node coeff norm = the function's L2 content in that
//                      box  ->  magnitude / spatial-support map.
//   mo_<i>_comp.json   compressed: internal nodes hold difference (wavelet)
//                      coefficients, so the per-node coeff norm = the local
//                      detail/error added at that scale  ->  error-at-resolution
//                      map (the essence of MRA).
//
// Output is consumed by the gecko `mra_tree` Python module, which turns the JSON
// into a flat box table + HDF5/HyperTreeGrid for ParaView and the analysis/viz
// layer (depth histogram, wavelet energy-per-scale spectrum, 2D slices, and the
// N x N product-overlap matrix used to study Coulomb/exchange partitioning).
//
// This tool LOADS ONLY; it never triggers an electronic solve (Tier-B sibling,
// doc 18). It reuses MADNESS's own tree dumper (print_tree_jsonfile), which
// already emits per-node {has_coeff, has_children, norm, norm_tree, snorm,
// dnorm, rank, dim} + owner, so no new tree-walking code is needed.
//
// All MADNESS Function ops below (reconstruct/compress, print_tree_jsonfile,
// norm2) are COLLECTIVE -- they run on every rank; only printing is gated on
// rank 0.
//
// Usage:
//   dump_mra_trees --archive=<prefix>.restartdata [--out=DIR]
//                  [--maxlevel=N] [--k=N] [--thresh=X] [--max-orbitals=N]

#include "../GroundState.hpp"
#include "../Perturbations.hpp"        // Perturbation (FD identity)
#include "../ResponseProtocol.hpp"
#include "../kernels/tags.hpp"         // Static / Full / TDA, ClosedShell
#include "../kernels/static.hpp"       // Kernels<Static,ClosedShell>::compute_density
#include "../kernels/full.hpp"         // Kernels<Full,ClosedShell>::compute_density
#include "../kernels/tda.hpp"          // Kernels<TDA,ClosedShell>::compute_density (ES tdens)
#include "../solvers/fd_save_load.hpp" // try_load_fd_state<Type,Shell>
#include "../solvers/es_save_load.hpp" // load_es_roots<TDA,ClosedShell> (ES bundle)

#include <madness/chem/atomutil.h>
#include <madness/external/nlohmann_json/json.hpp>
#include <madness/misc/info.h>
#include <madness/mra/funcplot.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>

// Native HyperTreeGrid (.htg) export is OPTIONAL, behind MADNESS_ENABLE_VTK
// (CMake -> defines MADNESS_HAS_VTK). A no-VTK build compiles unchanged and
// --htg becomes a no-op with a one-line notice; the JSON/.vtk path is untouched.
#ifdef MADNESS_HAS_VTK
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkHyperTreeGrid.h>
#include <vtkHyperTreeGridNonOrientedGeometryCursor.h>
#include <vtkNew.h>
#include <vtkXMLHyperTreeGridWriter.h>
// Overlapping-AMR (function values on the adaptive octree -> .vthb)
#include <vtkAMRBox.h>
#include <vtkFloatArray.h>
#include <vtkOverlappingAMR.h>
#include <vtkUniformGrid.h>
#include <vtkXMLUniformGridAMRWriter.h>
#include <array>
#endif

// Native MRA coefficient export (.mad.h5) is OPTIONAL, behind MADNESS_ENABLE_HDF5
// (CMake -> defines MADNESS_HAS_HDF5). It reuses the io-hdf5 structured Function
// writer (solvers/function_hdf5_io.hpp); a no-HDF5 build compiles unchanged and
// --coeffs becomes a one-line notice (the header is itself a no-op without the macro).
#ifdef MADNESS_HAS_HDF5
#include "../solvers/function_hdf5_io.hpp"
#endif

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <functional>
#include <map>
#include <sstream>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

using namespace madness;
using namespace molresponse_v3;

namespace {

// Mirror the molecule sourcing the FD/ES skeletons use: read the geometry from
// the moldft calc_info JSON sitting next to the archive. from_archive needs a
// Molecule to derive electron counts; an empty one starts the closed-shell load
// path but matching the canonical path avoids surprises.
Molecule load_molecule_near(const std::string& archive_path) {
  Molecule molecule;
  auto archive_dir = std::filesystem::path(archive_path).parent_path();
  for (const auto& name : {"moldft.calc_info.json", "mad.calc_info.json"}) {
    auto candidate = archive_dir / name;
    if (std::filesystem::exists(candidate)) {
      std::ifstream ifs(candidate);
      nlohmann::json j;
      ifs >> j;
      nlohmann::json mol_json;
      if (j.contains("tasks") && j["tasks"].is_array() && !j["tasks"].empty())
        mol_json = j["tasks"][0]["molecule"];
      else if (j.contains("molecule"))
        mol_json = j["molecule"];
      if (!mol_json.is_null()) molecule.from_json(mol_json);
      break;
    }
  }
  return molecule;
}

// Write one Function's octree as JSON, MPI-safe. `Function::print_tree_json` is
// COLLECTIVE -- every rank must call it (rank 0 walks the tree and pulls remote
// nodes behind internal fences) -- but it emits content only on rank 0. So we
// buffer into a stringstream on every rank and let ONLY rank 0 open and write
// the file. (MADNESS's own print_tree_jsonfile opens the ofstream on every rank
// and writes the wrapper on every rank, so under -np>1 all ranks race on the
// same path and produce corrupt JSON.)
void write_tree_json(World& world, const Function<double, 3>& f,
                     const std::string& path) {
  std::ostringstream oss;
  const Tensor<double> cell = FunctionDefaults<3>::get_cell();
  oss << "{\"cell\":[";
  for (int d = 0; d < 3; ++d) {
    oss << "[" << cell(d, 0) << "," << cell(d, 1) << "]";
    if (d != 2) oss << ",";
  }
  oss << "],\"tree\":{";
  f.print_tree_json(oss);  // collective: call on every rank; fills oss on rank 0
  oss << "}}";
  if (world.rank() == 0) {
    std::ofstream os(path);
    os << oss.str();
  }
}

// Write the nuclear geometry as VTK POLYDATA points in BOHR -- the same frame as
// the box dumps -- with per-atom atomic number Z and covalent radius, so ParaView
// can glyph the atoms as element-sized spheres overlaid on the octree. The
// Molecule is a small replicated object (not a distributed Function), so this is
// pure local data: write on rank 0 only, no collective ops, no fence.
void write_geometry_vtk(World& world, const Molecule& mol,
                        const std::string& path) {
  if (world.rank() != 0) return;
  const std::size_t n = mol.natom();
  std::ostringstream pts, zs, rs;
  for (std::size_t i = 0; i < n; ++i) {
    const Atom& a = mol.get_atom(i);
    pts << a.x << " " << a.y << " " << a.z << "\n";
    zs << a.atomic_number << "\n";
    rs << get_atomic_data(a.atomic_number).covalent_radius << "\n";
  }
  std::ofstream os(path);
  os << "# vtk DataFile Version 3.0\nMRA molecule geometry (bohr)\nASCII\n"
     << "DATASET POLYDATA\nPOINTS " << n << " float\n"
     << pts.str() << "VERTICES " << n << " " << 2 * n << "\n";
  for (std::size_t i = 0; i < n; ++i) os << "1 " << i << "\n";
  os << "POINT_DATA " << n << "\n"
     << "SCALARS Z float 1\nLOOKUP_TABLE default\n" << zs.str()
     << "SCALARS radius float 1\nLOOKUP_TABLE default\n" << rs.str();
}

// Sample one orbital on a regular grid and write a Gaussian .cube (geometry
// baked in) so ParaView/VMD can render it as signed isosurfaces -- the actual
// orbital lobes/nodes, which the box dumps (one norm per box) cannot show.
//
// The grid is a cubic box snug around the molecule (bbox +/- `pad` bohr), NOT the
// full ~200-bohr cell, so the points land where the orbital has support. Sampling
// uses Function::eval_cube -- the same batched, distributed grid evaluator behind
// plotvtk/plotdx (each box fills its grid cells), which is COLLECTIVE (all ranks
// call it) and orders of magnitude faster than per-point evaluation. Only rank 0
// then writes the file.
void write_cube_file(World& world, Function<double, 3> f, const Molecule& mol,
                     const std::string& path, int npoints, double pad) {
  double lo[3] = {1e9, 1e9, 1e9}, hi[3] = {-1e9, -1e9, -1e9};
  for (std::size_t a = 0; a < mol.natom(); ++a) {
    const Atom& at = mol.get_atom(a);
    const double c[3] = {at.x, at.y, at.z};
    for (int d = 0; d < 3; ++d) {
      lo[d] = std::min(lo[d], c[d]);
      hi[d] = std::max(hi[d], c[d]);
    }
  }
  Tensor<double> cell(3, 2);
  for (int d = 0; d < 3; ++d) {
    lo[d] -= pad;
    hi[d] += pad;
    cell(d, 0) = lo[d];
    cell(d, 1) = hi[d];
  }

  // Collective, batched grid evaluation -> full grid (rank 0 holds it).
  const std::vector<long> npt(3, npoints);
  const Tensor<double> grid = f.eval_cube(cell, npt);

  if (world.rank() != 0) return;
  double delta[3];
  for (int d = 0; d < 3; ++d) delta[d] = (hi[d] - lo[d]) / (npoints - 1);
  FILE* file = std::fopen(path.c_str(), "w");
  if (!file) return;
  std::fprintf(file, "MADNESS MRA orbital (bohr)\n%s\n", path.c_str());
  std::fprintf(file, "%d %12.6f %12.6f %12.6f\n", int(mol.natom()),
               lo[0], lo[1], lo[2]);
  std::fprintf(file, "%d %12.6f %12.6f %12.6f\n", npoints, delta[0], 0.0, 0.0);
  std::fprintf(file, "%d %12.6f %12.6f %12.6f\n", npoints, 0.0, delta[1], 0.0);
  std::fprintf(file, "%d %12.6f %12.6f %12.6f\n", npoints, 0.0, 0.0, delta[2]);
  for (std::size_t a = 0; a < mol.natom(); ++a) {
    const Atom& at = mol.get_atom(a);
    std::fprintf(file, "%d %12.6f %12.6f %12.6f %12.6f\n", at.atomic_number,
                 double(at.atomic_number), at.x, at.y, at.z);
  }
  long per_line = 0;
  for (int i = 0; i < npoints; ++i)
    for (int j = 0; j < npoints; ++j)
      for (int k = 0; k < npoints; ++k) {
        std::fprintf(file, "%13.5e ", grid(i, j, k));
        if (++per_line == 6) {
          std::fprintf(file, "\n");
          per_line = 0;
        }
      }
  std::fprintf(file, "\n");
  std::fclose(file);
}

// ---- native HyperTreeGrid (.htg) export -----------------------------------
//
// Build a vtkHyperTreeGrid from the LIVE octree (a single root tree over the
// whole simulation cell, branch factor 2) instead of the legacy one-hexahedron-
// per-leaf VTK. Geometry is implicit in the tree path, so the file is ~20-30x
// smaller and ParaView gets a native Maximum-Level depth control + per-level
// opacity. Mirrors the validated reference write_htg() in
// madness_studies/scripts/view_boxes.py (child order, rms formula, the cursor
// SetGlobalIndexStart gotcha, the bohr frame, the array names) -- but builds
// from the live tree rather than the dumped _boxes.vtk.
#ifdef MADNESS_HAS_VTK
struct HtgNode {
  double norm = 0.0;
  double dnorm = -1.0;
  int owner = 0;
  bool has_children = false;
};
// Keyed by (level, lx, ly, lz) so std::map ordering needs no Vector operator<.
using HtgKey = std::tuple<int, Translation, Translation, Translation>;
using HtgMap = std::map<HtgKey, HtgNode>;

// rank-0-only recursion mirroring FunctionImpl::do_print_tree_json (mraimpl.h):
// coeffs.find(key).get() pulls remote nodes while the OTHER ranks service the
// request inside the trailing gop.fence() in collect_live_tree(). The norm
// matches FunctionNode::print_json (funcimpl.h); owner = coeffs.owner(key), the
// live-container pmap owner -- the same source do_print_tree_json writes and the
// legacy _boxes.vtk uses (the recon-JSON owner can be unreliable; this is not).
template <typename ImplT>
void collect_live_tree_r(const ImplT& impl, const Key<3>& key, HtgMap& out) {
  const auto& coeffs = impl.get_coeffs();
  auto it = coeffs.find(key).get();
  if (it == coeffs.end()) return;  // absent child of a has_children node: skip
  const auto& node = it->second;
  HtgNode rec;
  rec.norm = node.has_coeff() ? node.coeff().normf() : 0.0;
  if (rec.norm < 1e-12) rec.norm = 0.0;
  rec.dnorm = node.get_dnorm();
  rec.owner = coeffs.owner(key);
  rec.has_children = node.has_children();
  const auto& l = key.translation();
  out[{int(key.level()), l[0], l[1], l[2]}] = rec;
  if (node.has_children())
    for (KeyChildIterator<3> kit(key); kit; ++kit)
      collect_live_tree_r(impl, kit.key(), out);
}

// Collective: rank 0 walks the tree; the fence lets its remote finds complete.
void collect_live_tree(World& world, const Function<double, 3>& f, HtgMap& out) {
  const auto impl = f.get_impl();
  if (world.rank() == 0)
    collect_live_tree_r(*impl, impl->get_cdata().key0, out);
  world.gop.fence();
}

// Build the single-root HTG from the node map and write it (rank 0 only). DFS
// tracks (level, translation) explicitly; child order c -> (c&1,(c>>1)&1,(c>>2)&1)
// with l_child = 2*l + (i,j,k); rms = norm / sqrt(box_volume),
// box_volume = cell_volume / 2^(3*level). Cell arrays match the legacy .vtk
// names so existing ParaView coloring (level/norm/rms/owner) works unchanged.
void write_htg(World& world, const HtgMap& nodes, const std::string& path) {
  if (world.rank() != 0 || nodes.empty()) return;
  const Tensor<double> cell = FunctionDefaults<3>::get_cell();
  double cell_volume = 1.0;
  for (int d = 0; d < 3; ++d) cell_volume *= (cell(d, 1) - cell(d, 0));

  vtkNew<vtkHyperTreeGrid> htg;
  htg->Initialize();
  htg->SetDimensions(2, 2, 2);  // 1x1x1 cells -> a single root tree
  htg->SetBranchFactor(2);
  for (int d = 0; d < 3; ++d) {
    vtkNew<vtkDoubleArray> coord;
    coord->SetNumberOfValues(2);
    coord->SetValue(0, cell(d, 0));
    coord->SetValue(1, cell(d, 1));
    if (d == 0) htg->SetXCoordinates(coord.Get());
    else if (d == 1) htg->SetYCoordinates(coord.Get());
    else htg->SetZCoordinates(coord.Get());
  }

  vtkNew<vtkDoubleArray> a_level, a_norm, a_rms, a_dnorm, a_owner;
  a_level->SetName("level");
  a_norm->SetName("norm");
  a_rms->SetName("rms");
  a_dnorm->SetName("dnorm");
  a_owner->SetName("owner");

  vtkNew<vtkHyperTreeGridNonOrientedGeometryCursor> cur;
  htg->InitializeNonOrientedGeometryCursor(cur.Get(), 0, /*create=*/true);
  cur->SetGlobalIndexStart(0);  // without this GetGlobalNodeIndex() == -1

  std::function<void(int, Translation, Translation, Translation)> rec =
      [&](int level, Translation lx, Translation ly, Translation lz) {
        const vtkIdType gid = cur->GetGlobalNodeIndex();
        auto it = nodes.find({level, lx, ly, lz});
        const HtgNode hn = (it != nodes.end()) ? it->second : HtgNode{};
        const double box_volume = cell_volume / std::pow(2.0, 3.0 * level);
        a_level->InsertValue(gid, double(level));
        a_norm->InsertValue(gid, hn.norm);
        a_rms->InsertValue(
            gid, box_volume > 0.0 ? hn.norm / std::sqrt(box_volume) : 0.0);
        a_dnorm->InsertValue(gid, hn.dnorm);
        a_owner->InsertValue(gid, double(hn.owner));
        if (it != nodes.end() && it->second.has_children) {
          cur->SubdivideLeaf();
          for (int c = 0; c < 8; ++c) {
            cur->ToChild(c);
            rec(level + 1, 2 * lx + (c & 1), 2 * ly + ((c >> 1) & 1),
                2 * lz + ((c >> 2) & 1));
            cur->ToParent();
          }
        }
      };
  rec(0, 0, 0, 0);  // root = (level 0, {0,0,0}) = FunctionImpl cdata.key0

  htg->GetCellData()->AddArray(a_level.Get());
  htg->GetCellData()->AddArray(a_norm.Get());
  htg->GetCellData()->AddArray(a_rms.Get());
  htg->GetCellData()->AddArray(a_dnorm.Get());
  htg->GetCellData()->AddArray(a_owner.Get());

  vtkNew<vtkXMLHyperTreeGridWriter> writer;
  writer->SetFileName(path.c_str());
  // Inline ASCII, NOT the default Appended mode: in VTK 9.1 the appended-data
  // offset-fixup pass mangles the header (writes a RangeMax/offset fragment over
  // byte 0, clobbering <VTKFile>), so ParaView/VTK 9.5 rejects it as not
  // well-formed. ASCII writes data inline (no seek-back). Switch to
  // SetDataModeToBinary() later for the size win once validated -- it is also
  // inline (base64) and avoids the buggy appended path.
  writer->SetDataModeToAscii();
  writer->SetInputData(htg.Get());
  writer->Write();
}

// Collect (collective) + write (rank 0). Safe to call on every rank.
void write_htg_from_function(World& world, const Function<double, 3>& f,
                             const std::string& path) {
  HtgMap nodes;
  collect_live_tree(world, f, nodes);
  write_htg(world, nodes, path);
}

// ---- adaptive function-value export: vtkOverlappingAMR (.vthb) -------------
//
// PROTOTYPE (see AMR_EXPORT_NOTES.md). HTG stores ONE scalar per box (skeleton);
// OverlappingAMR stores the FUNCTION sampled on a small m^3 uniform grid PER LEAF
// BOX across refinement levels, so ParaView's AMR Contour gives smooth isosurfaces
// at sub-box resolution while keeping adaptivity (fine blocks only where the tree
// refined) -- the adaptive replacement for the dense .cube. Octree => refinement
// ratio 2; origin/spacing in bohr (aligns with geometry.vtk / the .htg).
//
// The MADNESS-side mapping below (leaf walk, box geometry, level grouping,
// index-space AMR box, per-level spacing) is the load-bearing logic. The VTK
// OverlappingAMR API is version-sensitive -- the assembly calls are a best-effort
// first cut to verify at the first -DMADNESS_ENABLE_VTK build (NOTES §VTK-API).

// Leaf (level,lx,ly,lz) list replicated on ALL ranks: rank 0 walks the tree
// (collect_live_tree_r), keeps the leaves (no children), and broadcasts the flat
// list so the eval_cube sampling loop stays collective + in lockstep everywhere.
std::vector<std::array<long, 4>>
collect_leaf_keys(World& world, const Function<double, 3>& f) {
  std::vector<long> flat;
  if (world.rank() == 0) {
    HtgMap nodes;
    const auto impl = f.get_impl();
    collect_live_tree_r(*impl, impl->get_cdata().key0, nodes);
    for (const auto& [key, nd] : nodes)
      if (!nd.has_children) {
        flat.push_back(std::get<0>(key));
        flat.push_back(std::get<1>(key));
        flat.push_back(std::get<2>(key));
        flat.push_back(std::get<3>(key));
      }
  }
  world.gop.fence();  // pair rank 0's remote coeffs.find()s
  std::size_t nflat = flat.size();
  world.gop.broadcast(&nflat, 1, 0);
  flat.resize(nflat);
  if (nflat) world.gop.broadcast(flat.data(), nflat, 0);
  std::vector<std::array<long, 4>> leaves(nflat / 4);
  for (std::size_t i = 0; i < leaves.size(); ++i)
    leaves[i] = {flat[4 * i], flat[4 * i + 1], flat[4 * i + 2], flat[4 * i + 3]};
  return leaves;
}

// Sample f on every leaf box (m^3 uniform grid via the collective eval_cube; rank
// 0 receives each grid) and assemble a vtkOverlappingAMR; write <path> on rank 0.
void write_amr(World& world, Function<double, 3> f, const std::string& path,
               int m) {
  f.reconstruct();
  const auto leaves = collect_leaf_keys(world, f);
  const Tensor<double> cell = FunctionDefaults<3>::get_cell();
  double lo[3], width[3];
  for (int d = 0; d < 3; ++d) {
    lo[d] = cell(d, 0);
    width[d] = cell(d, 1) - cell(d, 0);
  }
  long Lmax = 0, Lmin = leaves.empty() ? 0 : leaves.front()[0];
  for (const auto& k : leaves) {
    Lmax = std::max(Lmax, k[0]);
    Lmin = std::min(Lmin, k[0]);
  }

  // collective sampling: ALL ranks call eval_cube per leaf; rank 0 keeps the grids
  struct Block { long n, lx, ly, lz; Tensor<double> g; };
  std::vector<Block> blocks;
  const std::vector<long> npt(3, m);
  for (const auto& k : leaves) {
    Tensor<double> bc(3, 2);
    for (int d = 0; d < 3; ++d) {
      const double h = width[d] / std::pow(2.0, double(k[0]));  // leaf box width
      const double dx = h / m;                                  // m-grid sub-cell width
      // Sample at the m CELL CENTERS, not box corners: the outer corner of a
      // boundary leaf sits exactly on the sim-cell edge, where eval_cube ->
      // plot_cube_kernel maps it to x=1.0(+fp eps) and trips its x in [0,1]
      // assertion. Centers are strictly interior (x in (0,1)) and also align the
      // sampled values with the m-cell block's cell-data positions (fixes the
      // half-cell offset noted in AMR_EXPORT_NOTES.md simplification #2).
      bc(d, 0) = lo[d] + double(k[d + 1]) * h + 0.5 * dx;       // first cell center
      bc(d, 1) = lo[d] + double(k[d + 1] + 1) * h - 0.5 * dx;   // last cell center
    }
    Tensor<double> g = f.eval_cube(bc, npt);  // collective; rank 0 holds the grid
    if (world.rank() == 0) blocks.push_back({k[0], k[1], k[2], k[3], copy(g)});
  }
  if (world.rank() != 0) return;

  // --- assemble on rank 0 (VTK OverlappingAMR) ---
  // RE-BASE the hierarchy to the populated level range [Lmin, Lmax]: MADNESS
  // leaves never sit at the top tree levels (the root box is always subdivided),
  // so AMR level 0 must map to Lmin, not the MADNESS root. An empty AMR level 0
  // segfaults vtkXMLUniformGridAMRReader::UpdateInformation (it cannot build the
  // metadata). AMR level index = n - Lmin; spacing/origin stay physical (bohr).
  const int nlev = int(Lmax - Lmin + 1);
  std::vector<int> per_level(nlev, 0);
  for (const auto& b : blocks) per_level[b.n - Lmin]++;
  // Review MED: the rebase guarantees a populated level 0 (Lmin), but an
  // INTERIOR level can still be empty if the leaf depths have a gap — that can
  // trip the same vtkXMLUniformGridAMRReader metadata path as an empty level 0.
  // Warn here; the full fix (compact to only populated levels with the correct
  // per-level refinement ratio) needs ParaView validation and belongs to the
  // viz/pymra thread.
  for (int lev = 1; lev + 1 < nlev; ++lev)
    if (per_level[lev] == 0)
      print("[--amr] WARNING: interior AMR level", lev, "(MADNESS n =",
            Lmin + lev, ") has NO leaf boxes; the .vthb may fail to load in "
            "some readers. Report to the viz thread if ParaView rejects it.");

  vtkNew<vtkOverlappingAMR> amr;
  amr->Initialize(nlev, per_level.data());
  double origin[3] = {lo[0], lo[1], lo[2]};
  amr->SetOrigin(origin);
  amr->SetGridDescription(VTK_XYZ_GRID);
  for (int lev = 0; lev < nlev; ++lev) {
    const long n = Lmin + lev;
    double sp[3];
    for (int d = 0; d < 3; ++d) sp[d] = (width[d] / std::pow(2.0, double(n))) / m;
    amr->SetSpacing(lev, sp);
    if (lev > 0) amr->SetRefinementRatio(lev - 1, 2);  // octree => 2x per level
  }
  std::vector<int> next(nlev, 0);
  for (const auto& b : blocks) {
    const long n = b.n;
    const int lev = int(n - Lmin), j = next[lev]++;
    const long l[3] = {b.lx, b.ly, b.lz};
    int boxlo[3], boxhi[3];
    for (int d = 0; d < 3; ++d) {
      boxlo[d] = int(l[d] * m);
      boxhi[d] = boxlo[d] + m - 1;
    }
    vtkAMRBox box(boxlo, boxhi);
    amr->SetAMRBox(lev, j, box);

    vtkNew<vtkUniformGrid> grid;
    double gorigin[3], sp[3];
    for (int d = 0; d < 3; ++d) {
      const double h = width[d] / std::pow(2.0, double(n));
      gorigin[d] = lo[d] + double(l[d]) * h;
      sp[d] = h / m;
    }
    grid->SetOrigin(gorigin);
    grid->SetSpacing(sp);
    grid->SetDimensions(m + 1, m + 1, m + 1);  // m^3 cells
    vtkNew<vtkFloatArray> vals;
    vals->SetName("psi");
    vals->SetNumberOfValues(vtkIdType(m) * m * m);
    for (int iz = 0; iz < m; ++iz)       // eval_cube g(ix,iy,iz); VTK cell x-fastest
      for (int iy = 0; iy < m; ++iy)
        for (int ix = 0; ix < m; ++ix)
          vals->SetValue((vtkIdType(iz) * m + iy) * m + ix, float(b.g(ix, iy, iz)));
    grid->GetCellData()->AddArray(vals.Get());
    amr->SetDataSet(lev, j, grid.Get());
  }
  vtkNew<vtkXMLUniformGridAMRWriter> w;
  w->SetFileName(path.c_str());
  w->SetDataModeToAscii();  // see write_htg: default Appended mode mangles VTK 9.1 output
  w->SetInputData(amr.Get());
  w->Write();
}

void write_amr_from_function(World& world, const Function<double, 3>& f,
                             const std::string& path, int m) {
  write_amr(world, f, path, m);
}
#else
// No-VTK build: collective-safe no-ops (nothing happens on any rank).
inline void write_htg_from_function(World&, const Function<double, 3>&,
                                    const std::string&) {}
inline void write_amr_from_function(World&, const Function<double, 3>&,
                                    const std::string&, int) {}
#endif  // MADNESS_HAS_VTK

// ---- native MRA coefficient export: .mad.h5 (io-hdf5 structured writer) -----
//
// Reuse the io-hdf5 structured Function writer (solvers/function_hdf5_io.hpp) to emit
// each function's RECONSTRUCTED octree as one <stem>.mad.h5: /meta (k, thresh, cell,
// tree_state, ...), /keys [n_nodes x (3+NDIM)] = (level, lx,ly,lz, has_children,
// has_coeff) per node, /coeffs [n_coeff_nodes x k^3] = the scaling-coefficient tensor
// per leaf. Unlike --amr (which point-SAMPLES each leaf onto an m^3 grid), this is the
// function's native multiresolution basis, so a reader can reconstruct the MADNESS
// expansion exactly at any point/zoom. "Don't roll a new serializer"
// we write through the same writer io-hdf5 built.
//
// The structured writer is double-only and MADNESS_CHECKs world.size()==1 (the local
// container == the whole tree only at NP=1) and t.size()==k^3 per coeff node (so the
// tree MUST be reconstructed -- a compressed (2k)^3 node would trip the check). We
// therefore reconstruct first and skip cleanly under NP>1. The skip is uniform across
// ranks (every rank evaluates world.size()>1 identically) -> no collective divergence.
#ifdef MADNESS_HAS_HDF5
void write_coeffs_hdf5(World& world, Function<double, 3> f,
                       const std::string& path) {
  if (world.size() > 1) {
    if (world.rank() == 0)
      print("  [HDF5] --coeffs needs NP=1 (structured writer is single-rank); "
            "skipping", path);
    return;
  }
  f.reconstruct();  // leaves -> k^3 scaling coefficients (the shape /coeffs expects)
  molresponse_v3::save_function_hdf5(f, path);
}
#else
// No-HDF5 build: collective-safe no-op (nothing happens on any rank).
inline void write_coeffs_hdf5(World&, Function<double, 3>, const std::string&) {}
#endif  // MADNESS_HAS_HDF5

// Dump one orbital's octree in both representations. `stem` is the output path
// prefix (no extension). With `htg`, also write the native HyperTreeGrid next to
// each JSON (<stem>_boxes.htg from the reconstructed tree, <stem>_error.htg from
// the compressed tree). Collective: must be called on every rank.
void dump_function_trees(World& world, Function<double, 3> f,
                         const std::string& stem, bool htg = false) {
  f.reconstruct();
  write_tree_json(world, f, stem + "_recon.json");
  if (htg) write_htg_from_function(world, f, stem + "_boxes.htg");
  f.compress();
  write_tree_json(world, f, stem + "_comp.json");
  if (htg) write_htg_from_function(world, f, stem + "_error.htg");
}

// Inverse of Perturbation::description() (Perturbations.hpp): turn a
// response_metadata.json <pert> key back into a Perturbation so we can drive
// try_load_fd_state. Mirrors the only formats description() emits:
//   dipole_<a>   magnetic_<a>   nuc_<atom>_<a>   (atom "*" = all-atoms, -1).
Perturbation parse_perturbation(const std::string& desc) {
  auto axis_of = [](char c) -> int {
    return c == 'x' ? 0 : c == 'y' ? 1 : c == 'z' ? 2 : -1;
  };
  if (desc.rfind("dipole_", 0) == 0)
    return Perturbation::dipole(axis_of(desc.back()));
  if (desc.rfind("magnetic_", 0) == 0)
    return Perturbation::magnetic(axis_of(desc.back()));
  if (desc.rfind("nuc_", 0) == 0) {
    const auto last = desc.rfind('_');                  // before the axis char
    const std::string mid = desc.substr(4, last - 4);   // atom tag between _ _
    return Perturbation::nuclear(mid == "*" ? -1 : std::stoi(mid),
                                 axis_of(desc.back()));
  }
  throw std::runtime_error("parse_perturbation: unrecognized '" + desc + "'");
}

// Load one converged FD response point at the ACTIVE protocol (caller sets
// FunctionDefaults via set_response_protocol) using the SAME collective
// try_load_fd_state the solver's restart path uses, then dump each response
// orbital's octree (recon + compressed) exactly like a ground MO -- and, with
// `cube`, the signed .cube isosurfaces too, so the response orbitals and rho^(1)
// render as lobes just like the ground MOs. Static (static α) carries x only;
// Full (dynamic α) carries x and y. Collective.
template <typename Type, typename Shell>
int dump_fd_point(World& world, const std::string& calc_dir,
                  const Perturbation& pert, double freq,
                  const std::string& out_dir, const std::string& label,
                  int max_orbitals, const vecfuncT& gs_amo,
                  const Molecule& mol, bool cube, int cube_npoints,
                  double cube_pad, bool htg, bool amr, int amr_m, bool coeffs) {
  auto loaded = try_load_fd_state<Type, Shell>(world, calc_dir, pert, freq);
  if (!loaded) {
    if (world.rank() == 0) print("  [FD] skip (no loadable bundle):", label);
    return 0;
  }
  const auto& store = loaded->state.responses.at(0);

  auto dump_block = [&](const std::vector<real_function_3d>& blk,
                        const char* tag) {
    const std::size_t n =
        (max_orbitals >= 0)
            ? std::min<std::size_t>(blk.size(),
                                    static_cast<std::size_t>(max_orbitals))
            : blk.size();
    for (std::size_t i = 0; i < n; ++i) {
      const std::string stem =
          out_dir + "/" + label + "_" + tag + std::to_string(i);
      dump_function_trees(world, blk[i], stem, htg);
      // Same in-memory Function -> .cube as the ground MOs (write_cube_file),
      // same stem as the _boxes.vtk so the local tools pick it up automatically.
      if (cube)
        write_cube_file(world, blk[i], mol, stem + ".cube", cube_npoints,
                        cube_pad);
      // vtkOverlappingAMR (.vthb) value export, same in-memory Function -> one
      // adaptive-isosurface dataset per response orbital (parity with --cube).
      if (amr) write_amr_from_function(world, blk[i], stem + ".vthb", amr_m);
      // Native MRA coefficients (.mad.h5): the function's basis, not a sample.
      if (coeffs) write_coeffs_hdf5(world, blk[i], stem + ".mad.h5");
      const double nrm = blk[i].norm2();  // collective -- all ranks
      if (world.rank() == 0)
        print("  dumped", label, tag, i, " norm2=", nrm);
    }
  };

  dump_block(store.x_alpha, "x");
  if constexpr (std::is_same_v<Type, Full>) dump_block(store.y_alpha, "y");

  // Response density rho^(1). Kernels::compute_density is the single source of
  // truth for the closed-shell convention (Static folds Y=X -> factor 4; Full
  // uses x+y -> factor 2) and reads only g0.amo, so a 1-field ResponseGroundState
  // suffices (no prepare()/Coulomb/Q). The GS orbitals were loaded at the GS
  // native k; reproject a LOCAL copy to the ACTIVE (FD) k so the phi*x products
  // conform (the dumped GS trees stay native -- this copy is for the density only).
  const int    kfd = FunctionDefaults<3>::get_k();
  const double tfd = FunctionDefaults<3>::get_thresh();
  vecfuncT amo_k;
  amo_k.reserve(gs_amo.size());
  for (const auto& o : gs_amo)
    amo_k.push_back(o.k() == kfd ? o : project(o, kfd, tfd, false));
  world.gop.fence();

  ResponseGroundState rgs;
  rgs.amo = amo_k;
  real_function_3d rho1 = Kernels<Type, Shell>::compute_density(world, rgs, store);
  const std::string rho1_stem = out_dir + "/" + label + "_rho1";
  dump_function_trees(world, rho1, rho1_stem, htg);
  if (cube)
    write_cube_file(world, rho1, mol, rho1_stem + ".cube", cube_npoints,
                    cube_pad);
  if (amr) write_amr_from_function(world, rho1, rho1_stem + ".vthb", amr_m);
  if (coeffs) write_coeffs_hdf5(world, rho1, rho1_stem + ".mad.h5");
  const double rnorm = rho1.norm2();  // collective -- all ranks
  if (world.rank() == 0) print("  dumped", label, "rho1  norm2=", rnorm);
  return 1;
}

// Walk response_metadata.json and dump every saved CLOSED-SHELL FD response
// point's trees into `out_dir`, alongside the ground orbitals, so a single
// mra_tree.py run yields a combined ground + response product-overlap matrix.
//
// The metadata file is small and local: every rank reads it to build the
// IDENTICAL job list (nlohmann objects are key-ordered), so the downstream
// collective loads + dumps stay in lockstep. `L` (box length, from the
// ground-state archive header) sets each point's protocol before loading.
void dump_fd_states(World& world, double L, const std::string& calc_dir,
                    const std::string& out_dir, int max_orbitals,
                    const vecfuncT& gs_amo, const Molecule& mol, bool cube,
                    int cube_npoints, double cube_pad, bool htg, bool amr,
                    int amr_m, bool coeffs) {
  const std::string meta_path = calc_dir + "/response_metadata.json";
  if (!std::filesystem::exists(meta_path)) {
    if (world.rank() == 0)
      print("  [FD] no response_metadata.json in", calc_dir, "-- nothing to dump");
    return;
  }
  nlohmann::json j;
  {
    std::ifstream ifs(meta_path);
    ifs >> j;
  }
  if (!j.contains("fd_states")) {
    if (world.rank() == 0) print("  [FD] no fd_states in", meta_path);
    return;
  }
  const auto protocols = j.value("protocols", nlohmann::json::object());

  int npoints = 0;
  for (auto& [pert, by_key] : j["fd_states"].items()) {
    Perturbation p;
    try {
      p = parse_perturbation(pert);
    } catch (const std::exception& e) {
      if (world.rank() == 0) print("  [FD] skip:", e.what());
      continue;
    }
    for (auto& [key, by_fkey] : by_key.items()) {
      // (thresh, k) recorded by save_fd_state -> set the protocol; no key parse.
      if (!protocols.contains(key)) continue;
      const double thresh = protocols[key].value("thresh", 0.0);
      const int    k      = protocols[key].value("k", 0);
      if (thresh <= 0.0 || k <= 0) continue;
      set_response_protocol(world, L, thresh, k);

      for (auto& [fkey, entry] : by_fkey.items()) {
        if (entry.value("shell", "") != "closed_shell") {
          if (world.rank() == 0)
            print("  [FD] skip non-closed-shell:", pert, key, fkey);
          continue;
        }
        const double freq      = entry.value("freq", 0.0);
        const std::string type = entry.value("type", "");
        const std::string label = "fd_" + pert + "__" + key + "__" + fkey;
        // FD points are static (α(0)) or full (dynamic α(ω)); TDA is an ES-bundle
        // type and never appears in fd_states (FDPerturbationOf<TDA> is undefined).
        if (type == "static")
          npoints += dump_fd_point<Static, ClosedShell>(
              world, calc_dir, p, freq, out_dir, label, max_orbitals, gs_amo,
              mol, cube, cube_npoints, cube_pad, htg, amr, amr_m, coeffs);
        else if (type == "full")
          npoints += dump_fd_point<Full, ClosedShell>(
              world, calc_dir, p, freq, out_dir, label, max_orbitals, gs_amo,
              mol, cube, cube_npoints, cube_pad, htg, amr, amr_m, coeffs);
        else if (world.rank() == 0)
          print("  [FD] skip unknown type:", type, "for", label);
      }
    }
  }
  if (world.rank() == 0)
    print("DUMP_MRA_FD  calc_dir=", calc_dir, "  points_dumped=", npoints);
}

// Walk response_metadata.json excited_states and dump every converged CLOSED-SHELL
// TDA excited-state's response orbitals (x_i) + transition density into `out_dir`,
// alongside the ground + FD dumps. Mirrors dump_fd_states: the ES loader
// (load_es_roots) and all export writers already exist. Collective.
void dump_es_roots(World& world, double L, const std::string& calc_dir,
                   const std::string& out_dir, int max_orbitals,
                   const vecfuncT& gs_amo, const Molecule& mol, bool cube,
                   int cube_npoints, double cube_pad, bool htg, bool amr,
                   int amr_m, bool coeffs) {
  const std::string meta_path = calc_dir + "/response_metadata.json";
  if (!std::filesystem::exists(meta_path)) {
    if (world.rank() == 0)
      print("  [ES] no response_metadata.json in", calc_dir, "-- nothing to dump");
    return;
  }
  nlohmann::json j;
  { std::ifstream ifs(meta_path); ifs >> j; }
  if (!j.contains("excited_states")) {
    if (world.rank() == 0) print("  [ES] no excited_states in", meta_path);
    return;
  }
  const auto protocols = j.value("protocols", nlohmann::json::object());

  int nroots = 0;
  for (auto& [key, blk] : j["excited_states"].items()) {
    if (!protocols.contains(key)) continue;
    const double thresh = protocols[key].value("thresh", 0.0);
    const int    k      = protocols[key].value("k", 0);
    if (thresh <= 0.0 || k <= 0) continue;
    set_response_protocol(world, L, thresh, k);

    ESSolver<TDA, ClosedShell>::State state;
    try {
      state = load_es_roots<TDA, ClosedShell>(world, calc_dir + "/es__" + key);
    } catch (const std::exception& e) {
      if (world.rank() == 0) print("  [ES] skip", key, "--", e.what());
      continue;
    }

    // GS orbitals reprojected to the ES k for the transition density (same
    // discipline as dump_fd_point's rho1: reproject a LOCAL copy).
    const int    kfd = FunctionDefaults<3>::get_k();
    const double tfd = FunctionDefaults<3>::get_thresh();
    vecfuncT amo_k;
    amo_k.reserve(gs_amo.size());
    for (const auto& o : gs_amo)
      amo_k.push_back(o.k() == kfd ? o : project(o, kfd, tfd, false));
    world.gop.fence();
    ResponseGroundState rgs;
    rgs.amo = amo_k;

    const int nr = static_cast<int>(state.roots.size());
    for (int f = 0; f < nr; ++f) {
      const auto& root = state.roots[f];
      // Name to pymra.web's field contract: es_<state>__<protocol>_<comp>
      // (web.py _RE_ES). Comp suffixes _x<i>/_tdens are appended below, matching
      // the FD path's fd_<pert>_<dir>__<protocol>__f<omega>_<comp> convention.
      const std::string label = "es_" + std::to_string(f) + "__" + key;
      const std::size_t no =
          (max_orbitals >= 0)
              ? std::min<std::size_t>(root.x_alpha.size(),
                                      static_cast<std::size_t>(max_orbitals))
              : root.x_alpha.size();
      for (std::size_t i = 0; i < no; ++i) {
        const std::string stem = out_dir + "/" + label + "_x" + std::to_string(i);
        dump_function_trees(world, root.x_alpha[i], stem, htg);
        if (cube)
          write_cube_file(world, root.x_alpha[i], mol, stem + ".cube",
                          cube_npoints, cube_pad);
        if (amr) write_amr_from_function(world, root.x_alpha[i], stem + ".vthb", amr_m);
        if (coeffs) write_coeffs_hdf5(world, root.x_alpha[i], stem + ".mad.h5");
        const double nrm = root.x_alpha[i].norm2();  // collective
        if (world.rank() == 0) print("  dumped", label, "x", i, " norm2=", nrm);
      }
      // Transition density (TDA, closed-shell): the single-source-of-truth
      // Kernels density of the root against the ground state.
      real_function_3d tdens =
          Kernels<TDA, ClosedShell>::compute_density(world, rgs, root);
      const std::string tstem = out_dir + "/" + label + "_tdens";
      dump_function_trees(world, tdens, tstem, htg);
      if (cube)
        write_cube_file(world, tdens, mol, tstem + ".cube", cube_npoints, cube_pad);
      if (amr) write_amr_from_function(world, tdens, tstem + ".vthb", amr_m);
      if (coeffs) write_coeffs_hdf5(world, tdens, tstem + ".mad.h5");
      const double tn = tdens.norm2();  // collective
      if (world.rank() == 0) print("  dumped", label, "tdens  norm2=", tn);
      ++nroots;
    }
  }
  if (world.rank() == 0)
    print("DUMP_MRA_ES  calc_dir=", calc_dir, "  roots_dumped=", nroots);
}

}  // namespace

int main(int argc, char** argv) {
  World& world = initialize(argc, argv);
  int rc = 0;
  try {
    startup(world, argc, argv, true);

    commandlineparser parser(argc, argv);
    if (!parser.key_exists("archive")) {
      if (world.rank() == 0) {
        print("Usage: dump_mra_trees --archive=<prefix>.restartdata "
              "[--out=DIR] [--maxlevel=N] [--k=N] [--thresh=X] "
              "[--max-orbitals=N] [--cube] [--cube-npoints=N] [--cube-pad=X] "
              "[--htg] [--amr] [--amr-m=N] [--coeffs] [--fd] [--fd-calc-dir=DIR]");
        print("  Loads ground-state orbitals at their native (k, thresh) and");
        print("  dumps per-orbital MRA octree JSON (reconstructed + compressed),");
        print("  geometry.vtk, and (with --cube) per-orbital .cube isosurfaces.");
        print("  With --htg (needs a -DMADNESS_ENABLE_VTK=ON build), also writes");
        print("  each tree as a native VTK HyperTreeGrid (<stem>_boxes.htg /");
        print("  <stem>_error.htg) next to the legacy JSON/.vtk. With --amr,");
        print("  writes vtkOverlappingAMR (<stem>.vthb; function values on an");
        print("  --amr-m^3 grid/leaf box -> adaptive isosurfaces) for the ground");
        print("  MOs and (with --fd) the response orbitals + rho^(1).");
        print("  With --coeffs (needs -DMADNESS_ENABLE_HDF5=ON, NP=1), writes each");
        print("  function's NATIVE MRA coefficients as <stem>.mad.h5 (/meta+/keys+");
        print("  /coeffs: per-leaf n,l,s[k^3] + cell,k,thresh) -- the exact basis");
        print("  for offline reconstruction, not a point sample like .cube/.vthb.");
        print("  With --fd, also dumps the converged closed-shell FD response");
        print("  orbitals from --fd-calc-dir (default: the archive's directory)");
        print("  into the same out dir -> one combined ground+response analysis.");
      }
      finalize();
      return 1;
    }

    // Inner scope so MADNESS objects destruct before finalize().
    {
      const std::string archive_path = parser.value_raw("archive");
      const std::string out_dir =
          parser.key_exists("out") ? parser.value_raw("out")
                                   : std::string("mra_trees");
      const int max_orbitals = parser.key_exists("max-orbitals")
                                   ? std::stoi(parser.value("max-orbitals"))
                                   : -1;
      // Optional orbital isosurface export (.cube). Off by default -- it is the
      // expensive step (per-point grid evaluation).
      const bool cube = parser.key_exists("cube");
      const int cube_npoints = parser.key_exists("cube-npoints")
                                   ? std::stoi(parser.value("cube-npoints")) : 80;
      const double cube_pad = parser.key_exists("cube-pad")
                                  ? std::stod(parser.value("cube-pad")) : 6.0;
      // FD response orbitals: dump the converged closed-shell FD response states
      // from a calc dir (response_metadata.json) into the same out dir, so one
      // mra_tree.py run yields a combined ground + response overlap matrix.
      const bool dump_fd = parser.key_exists("fd");
      // ES excited-state response orbitals + transition densities from the same
      // calc dir (response_metadata.json excited_states + es__<key>/ bundles).
      const bool dump_es = parser.key_exists("es");
      const std::string fd_calc_dir =
          parser.key_exists("fd-calc-dir")
              ? parser.value_raw("fd-calc-dir")
              : std::filesystem::path(archive_path).parent_path().string();
      // Native HyperTreeGrid (.htg) export of every tree, in addition to the
      // JSON/.vtk. Needs a -DMADNESS_ENABLE_VTK=ON build; otherwise a no-op.
      const bool htg = parser.key_exists("htg");
      // Adaptive function-value export (.vthb, vtkOverlappingAMR) -- samples each
      // leaf box on an m^3 grid (--amr-m, default 8). Covers the ground MOs and
      // (with --fd) the response orbitals + rho^(1), matching --cube/--htg.
      const bool amr = parser.key_exists("amr");
      const int amr_m = parser.key_exists("amr-m")
                            ? std::stoi(parser.value("amr-m")) : 8;
      // Native MRA coefficient export (.mad.h5) via the io-hdf5 structured writer:
      // per-leaf (n, l, s[k^3]) + header (cell, k, thresh). Needs a
      // -DMADNESS_ENABLE_HDF5=ON build and NP=1; otherwise skipped with a notice.
      // Covers the ground MOs and (with --fd) the response orbitals + rho^(1),
      // matching --cube/--amr.
      const bool coeffs = parser.key_exists("coeffs");
#ifndef MADNESS_HAS_VTK
      if ((htg || amr) && world.rank() == 0)
        print("  [VTK] --htg/--amr requested but this binary was built without "
              "VTK; reconfigure with -DMADNESS_ENABLE_VTK=ON. Skipping .htg/.vthb "
              "(JSON/.vtk/.cube output is unaffected).");
#endif
#ifndef MADNESS_HAS_HDF5
      if (coeffs && world.rank() == 0)
        print("  [HDF5] --coeffs requested but this binary was built without "
              "HDF5; reconfigure with -DMADNESS_ENABLE_HDF5=ON. Skipping .mad.h5 "
              "(JSON/.vtk/.cube/.vthb output is unaffected).");
#endif

      auto header = GroundState::read_archive_header(world, archive_path);

      // Load at NATIVE resolution unless explicitly overridden. Matching the
      // archive's (k, thresh) means SCF::load_mos performs no reprojection, so
      // the tree we walk is exactly the one the solver produced. Reprojecting
      // would rebuild the tree and we'd be visualizing a different object.
      const int k =
          parser.key_exists("k") ? std::stoi(parser.value("k")) : header.k;
      const double thresh = parser.key_exists("thresh")
                                ? std::stod(parser.value("thresh"))
                                : header.converged_for_thresh;
      set_response_protocol(world, header.L, thresh, k);

      auto molecule = load_molecule_near(archive_path);
      auto gs = GroundState::from_archive(world, archive_path, molecule);
      if (world.rank() == 0) gs.print_info();

      if (world.rank() == 0) std::filesystem::create_directories(out_dir);
      world.gop.fence();

      write_geometry_vtk(world, gs.molecule(), out_dir + "/geometry.vtk");
      if (world.rank() == 0) {
        print("  wrote geometry.vtk  (", gs.molecule().natom(), "atoms )");
      }

      const vecfuncT& mos = gs.orbitals_alpha();
      const std::size_t n =
          (max_orbitals >= 0)
              ? std::min<std::size_t>(mos.size(),
                                      static_cast<std::size_t>(max_orbitals))
              : mos.size();

      if (world.rank() == 0) {
        print("DUMP_MRA  archive=", archive_path, "  out=", out_dir,
              "  k=", k, "  thresh=", thresh, "  n_orbitals=", n);
      }

      for (std::size_t i = 0; i < n; ++i) {
        const std::string stem = out_dir + "/mo_" + std::to_string(i);
        dump_function_trees(world, mos[i], stem, htg);
        if (cube) {
          write_cube_file(world, mos[i], gs.molecule(), stem + ".cube",
                          cube_npoints, cube_pad);
        }
        if (amr) write_amr_from_function(world, mos[i], stem + ".vthb", amr_m);
        if (coeffs) write_coeffs_hdf5(world, mos[i], stem + ".mad.h5");
        const double nrm = mos[i].norm2();  // collective -- all ranks
        if (world.rank() == 0) print("  dumped mo_", i, "  norm2=", nrm);
      }
      world.gop.fence();

      // Combined ground + response analysis: dump the FD response orbitals into
      // the same out dir. Runs after the GS dump (it mutates FunctionDefaults to
      // each FD point's protocol). header.L is the box length from the GS header.
      if (dump_fd) {
        dump_fd_states(world, header.L, fd_calc_dir, out_dir, max_orbitals, mos,
                       gs.molecule(), cube, cube_npoints, cube_pad, htg, amr,
                       amr_m, coeffs);
        world.gop.fence();
      }

      // ES excited-state orbitals + transition densities (same calc dir as --fd).
      if (dump_es) {
        dump_es_roots(world, header.L, fd_calc_dir, out_dir, max_orbitals, mos,
                      gs.molecule(), cube, cube_npoints, cube_pad, htg, amr,
                      amr_m, coeffs);
        world.gop.fence();
      }
    }

    finalize();
    return rc;
  } catch (const std::exception& e) {
    if (world.rank() == 0) print("Error:", e.what());
    finalize();
    return 1;
  }
}
