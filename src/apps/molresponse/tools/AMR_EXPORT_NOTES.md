# `dump_mra_trees --amr` — Overlapping-AMR function-value export (PROTOTYPE)

**Status:** **BUILT + VALIDATED (2026-06-22, xm047, NP=1, c2h4).** Both `--htg`
(HyperTreeGrid) and `--amr` (Overlapping-AMR `.vthb`) build with
`-DMADNESS_ENABLE_VTK=ON` against the cluster VTK 9.1 and **parse in ParaView 6.0.1
(VTK 9.5)** via pvpython: HTG = 1833 nodes + arrays [level,norm,rms,dnorm,owner];
AMR = `vtkOverlappingAMR` 17 levels / 1604 blocks, `psi` present, reader OK.
Getting there fixed three prototype bugs (all in `dump_mra_trees.cpp`):
1. **Writer corruption** — the default `Appended` data mode mangles the file header
   in VTK 9.1 (offset fragment overwrites `<VTKFile>` → ParaView "not well-formed");
   fix = `writer->SetDataModeToAscii()` on both XML writers.
2. **`--amr` crash** — per-leaf `eval_cube` sampled box corners; a boundary leaf's
   outer corner maps to `x=1.0+ε` and trips `plot_cube_kernel`'s `x∈[0,1]` assert;
   fix = sample at the m cell **centers** (strictly interior; also fixes offset #2).
3. **`.vthb` reader segfault** — MADNESS leaves never sit at the top tree levels, so
   AMR levels 0/1 were empty and `vtkXMLUniformGridAMRReader::UpdateInformation`
   segfaulted; fix = **re-base** the AMR hierarchy to `[Lmin, Lmax]` (AMR level 0 =
   coarsest populated MADNESS level).
Remaining: (a) multi-rank (NP>1) launch is blocked in a nested `srun` step on this
OpenMPI 5.0.2/SLURM (`mpirun --mca plm isolated` and `srun --mpi=pmix` both fail) —
NP=1 fully exercises the rank-0 export, so the multi-rank `owner` view is deferred;
(b) ParaView visual pass (AMR Contour smoothness, blanking) still to eyeball in the
GUI; (c) `SetDataModeToBinary()` is the size-optimized follow-up to ASCII.

**Build prerequisite RESOLVED (2026-06-19, env check, no alloc):** the buildable
VTK on the cluster is **VTK 9.1.0** at
`VTK_DIR=/gpfs/software/fsl/6.0.6.5/lib/cmake/vtk-9.1` (FSL package). A standalone
`find_package(VTK 9.0 REQUIRED COMPONENTS CommonCore CommonDataModel IOXML)`
configure resolves against it cleanly, so the current `CMakeLists.txt` component
list is sufficient — **`FiltersAMR` is NOT needed** (resolves the §VTK-API hedge
below). Every HTG + AMR VTK call in this prototype was checked against the 9.1
headers and the signatures match (see §VTK-API). The `paraview/6.0.1` module is a
**view-only** runtime binary (no dev headers / `*Config.cmake`) — usable to open
`.htg`/`.vthb`, not to build. Only-unverified step left = link/runtime ABI of the
conda-built VTK 9.1 against the GCC 13.2.0 static MADNESS build (settles at the
alloc build; put `/gpfs/software/fsl/6.0.6.5/lib` on `LD_LIBRARY_PATH` to run).

## Why a third export format

Three complementary views of one MADNESS `Function`:

| format | what it stores | isosurfaces | scales? |
|---|---|---|---|
| `.cube` (have) | function on a **dense uniform grid** | smooth | no — uniform fine grid everywhere |
| `.htg` (have) | **one scalar per leaf box** (skeleton: level/norm/rms/owner) | blocky (cell-resolution) | yes (compact) |
| **`.vthb` (this)** | function sampled on an **m³ grid per leaf box, across levels** | **smooth + adaptive** | yes (fine blocks only where refined) |

`.vthb` (`vtkOverlappingAMR`) is the adaptive replacement for the dense cube: ParaView's
AMR Contour gives smooth isosurfaces at sub-box resolution, but only the refined regions
carry fine blocks — so it stays small for large molecules where a uniform cube blows up.
This is the natural home for the "k³ quadrature values per box" idea (HTG, being one-value-
per-cell, cannot hold sub-box grids).

## MADNESS → OverlappingAMR mapping (the load-bearing logic — verified)

- **Leaf walk:** rank 0 walks the reconstructed tree (`collect_live_tree_r`), keeps leaves
  (`!has_children`), broadcasts the flat `(level,lx,ly,lz)` list to all ranks so the
  sampling loop is collective + in lockstep.
- **Per-leaf sampling:** every rank calls `f.eval_cube(box_cell, {m,m,m})` for each leaf
  (collective; rank 0 receives the m³ grid). `box_cell` (bohr) = `[lo+l·h, lo+(l+1)·h]`,
  `h = cell_width / 2^level` per axis.
- **AMR assembly (rank 0):** refinement ratio **2** (octree). Global origin = `cell_lo`.
  Level-n cell spacing = `(cell_width/2^n)/m`. Each leaf box → one `vtkUniformGrid` block
  (origin = box lower corner bohr, m³ cells) at AMR level n, with index-space
  `vtkAMRBox` lo = `l·m`, hi = `l·m + m − 1`. Values attached as cell data `"psi"`.
- Output `<stem>.vthb` via `vtkXMLUniformGridAMRWriter`. Wired for the ground MOs
  **and** (with `--fd`) the FD response orbitals + ρ⁽¹⁾ — same call-chain threading
  as `--cube` (`dump_fd_states → dump_fd_point → dump_block`), parity with
  `--cube`/`--htg`.

## Known prototype simplifications (intentional — fix during iteration)

1. **Per-box `eval_cube` is O(#leaves) collective calls** — fine for h2o/c2h4, slow at
   scale (each call fences). Production path: transform each leaf's scaling coefficients
   to point values **locally** (no per-box collective) — `node.coeff()` → values via the
   quadrature transform. Prototype favors the simple `eval_cube` reuse.
2. **Cell-vs-point half-cell offset:** the m³ `eval_cube` samples (grid points, inclusive
   endpoints) are attached as **cell** data on an m-cell block — a half-cell offset.
   Acceptable for a first look; tighten by sampling at cell centers (offset the cell by
   h/2m) or by using point data with matching `vtkAMRBox` dims.
3. **Blanking** (coarse cells covered by finer blocks) is left to ParaView's AMR overlap
   resolution from the box geometry; we do not pre-blank. Verify it renders correctly;
   if not, call `vtkOverlappingAMR::GenerateParentChildInformation()` / set blanking.
4. ~~**Ground-state only** — FD response/ρ⁽¹⁾ not yet wired (trivial follow-up).~~
   **DONE:** `--amr` now threads through `dump_fd_states → dump_fd_point →
   dump_block` (alongside `--cube`), so `--fd --amr` emits `fd_*_x*.vthb`,
   `fd_*_y*.vthb` (Full), and `fd_*_rho1.vthb`. Pending alloc rebuild + ParaView
   eyeball of the FD `.vthb` outputs.

## VTK API — CONFIRMED against cluster VTK 9.1.0 (2026-06-19, header check)

All calls below were checked against `/gpfs/software/fsl/6.0.6.5/include/vtk-9.1`
and the signatures match — no version drift from the 9.5-era prototype:
- `vtkOverlappingAMR::Initialize(int, const int*)` — `vtkUniformGridAMR.h:55`. ✓
- `SetGridDescription(int)` — on the parent `vtkUniformGridAMR.h:60`; `VTK_XYZ_GRID`
  (= 8) in `vtkStructuredData.h:46`. (`<vtkStructuredData.h>` not yet `#include`d —
  it transitively resolves today, but add it if the constant goes unfound.) ✓
- `SetSpacing(unsigned int, const double[3])` — `vtkOverlappingAMR.h:70`. ✓
- `SetRefinementRatio(unsigned int, int)` — `vtkOverlappingAMR.h:114`. ✓
- `SetAMRBox(unsigned int, unsigned int, const vtkAMRBox&)` — `vtkOverlappingAMR.h:78`;
  `vtkAMRBox(const int lo[3], const int hi[3])` — `vtkAMRBox.h:61`. ✓
- `SetDataSet(unsigned int, unsigned int, vtkUniformGrid*)` — `vtkUniformGridAMR.h:89`. ✓
- `vtkXMLUniformGridAMRWriter` / `vtkXMLHyperTreeGridWriter` — headers present (IOXML). ✓
- HTG: `SetDimensions(int,int,int)`/`SetBranchFactor`/`SetXCoordinates` and the cursor's
  `SetGlobalIndexStart`/`GetGlobalNodeIndex`/`SubdivideLeaf`/`ToChild`/`ToParent`. ✓
- **CMake components:** the existing `find_package(VTK 9.0 REQUIRED COMPONENTS CommonCore
  CommonDataModel IOXML)` is **sufficient** — it resolves against 9.1 (the AMR data-model
  classes live in CommonDataModel, the XML AMR writer in IOXML). **`FiltersAMR`/`IOAMR`
  are NOT needed.** Re-check only if the *link/ABI* step (conda VTK vs GCC 13.2 static
  MADNESS) fails — that, not the component list, is the remaining unknown.

## Build + validate (on the alloc)

```bash
# via the cm.sh harness (note: export VTK_DIR BEFORE cm_use — it bakes the flag in):
source /gpfs/projects/rjh/adrian/madness_studies/es_bench/cm.sh
cm_arch <arch>                  # e.g. the xeonmax/amd96 you build on
export VTK_DIR=/gpfs/software/fsl/6.0.6.5/lib/cmake/vtk-9.1
cm_use feat/amr-export          # -> CM_CMAKE_EXTRA="-DMADNESS_ENABLE_VTK=ON -DVTK_DIR=$VTK_DIR"
cm_build                        # configures fresh + builds dump_mra_trees (+ v3 tests)
# run (FSL VTK is conda-built -> put its lib on the runtime path):
export LD_LIBRARY_PATH=/gpfs/software/fsl/6.0.6.5/lib:$LD_LIBRARY_PATH
mpirun -np 4 $BUILD/.../dump_mra_trees --archive=.../h2o/mad.restartdata --out=amr --htg --amr --amr-m=8
```
Then in ParaView: open `amr/mo_0.vthb` → **Contour** on `psi` → expect a **smooth**
isosurface (vs the blocky `.htg`), with **finer blocks near the nuclei**. Cross-check the
lobe shape against `mo_0.cube`'s isosurface (should match) and compare file size +
ParaView load time vs the `.cube` (the adaptivity payoff grows with molecule size).
```
