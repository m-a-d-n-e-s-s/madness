#ifndef MOLRESPONSE_V3_SOLVERS_DALTON_GS_SEED_HPP
#define MOLRESPONSE_V3_SOLVERS_DALTON_GS_SEED_HPP

// ===========================================================================
// GROUND-STATE SEED FROM A DALTON MOLDEN FILE (2026-09-09) — the library form
// of tools/seed_moldft_from_dalton, so the madqc response workflow can seed
// moldft itself when the deck carries `dalton.dir` (the seeding showcase:
// every stage starts from the basis-set calculation).
//
// Projects the n_occ occupied DALTON/molden MOs onto the MRA basis at the
// run's box L and (thresh, k), Loewdin-orthonormalizes them (independent MRA
// projections leave O(1e-4) mutual overlaps that SCF::load_mos does not clean,
// and an unorthonormal first density collapsed water to -106.9 Ha), and writes
// <prefix>.restartdata through RestartMetadata (chem/Restart.h) — the same
// header SCF::save_mos writes, so SCF::load_mos / `restart auto` resume from it
// as from any other archive. converged_for_dconv is left "unknown", so the
// plan is `iterate`, never `read_only`: the DALTON orbitals are a guess.
//
// Conventions (pinned by the tool's A/B, water 10 -> 5 SCF iterations):
//   * aocc = 1.0 per occupied spatial orbital (moldft doubles for closed
//     shell; the molden's chemist occ=2.0 would double the density).
//   * spin_restricted, alpha block only (closed-shell).
//   * L MUST equal the deck's box (load_mos throws otherwise); k/thresh
//     mismatches are re-projected by load_mos, so seeding at the first
//     protocol rung is fine.
// ===========================================================================

#include "../tools/dalton_gto.hpp"
#include "../ResponseProtocol.hpp"
#include "function_hdf5_io.hpp"   // HDF5 twin of the seed archive (no-op without MADNESS_HAS_HDF5)

#include <madness/chem/Restart.h>
#include <madness/chem/molecule.h>
#include <madness/mra/mra.h>
#include <madness/world/MADworld.h>
#include <madness/world/parallel_archive.h>

#include <algorithm>
#include <cctype>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace molresponse_v3 {
namespace dalton_gs_seed_detail {

/// One DALTON MO (its AO-coefficient column) as a real MRA functor.
class DaltonMOFunctor : public madness::FunctionFunctorInterface<double, 3> {
  const DaltonMoldenBasis &basis;
  std::vector<double> c;
  std::vector<madness::coord_3d> centers;
public:
  DaltonMOFunctor(const DaltonMoldenBasis &b, std::vector<double> w)
      : basis(b), c(std::move(w)) {
    for (const auto &sh : basis.shells) centers.push_back({sh.cx, sh.cy, sh.cz});
  }
  double operator()(const madness::coord_3d &r) const override {
    double val = 0.0, bf[9];
    for (size_t s = 0; s < basis.shells.size(); s++) {
      const auto &sh = basis.shells[s];
      sh.evaluate(r[0], r[1], r[2], bf);
      const int off = basis.ao_offsets[s];
      for (int k = 0; k < sh.n_ao; k++) val += c[static_cast<size_t>(off + k)] * bf[k];
    }
    return val;
  }
  std::vector<madness::coord_3d> special_points() const override { return centers; }
};

/// Element symbol (molden labels carry index suffixes: H1, H2, O) -> Z.
inline int symbol_to_Z(const std::string &s) {
  std::string elem;
  for (char c : s) { if (std::isalpha(static_cast<unsigned char>(c))) elem += c; else break; }
  static const std::vector<std::pair<std::string, int>> tab = {
      {"H", 1}, {"He", 2}, {"Li", 3}, {"Be", 4}, {"B", 5}, {"C", 6},
      {"N", 7}, {"O", 8}, {"F", 9}, {"Ne", 10}, {"Na", 11}, {"Mg", 12},
      {"Al", 13}, {"Si", 14}, {"P", 15}, {"S", 16}, {"Cl", 17}, {"Ar", 18}};
  for (const auto &[sym, z] : tab) if (sym == elem) return z;
  throw std::runtime_error("dalton gs seed: unknown element symbol '" + s + "'");
}

} // namespace dalton_gs_seed_detail

struct GsSeedOptions {
  double      L        = 200.0;   ///< box half-width; MUST match the run's `l`
  double      thresh   = 1e-4;    ///< projection thresh (k from the protocol table)
  double      energy   = 0.0;     ///< provenance only
  std::string xc       = "hf";    ///< provenance only (load_mos discards)
  std::string localize = "canon"; ///< provenance only
  int         nio      = 1;
  /// Extra archive prefixes to write the SAME seed to (e.g. "<prefix>.gs_seed"):
  /// moldft's save_mos overwrites <prefix>.restartdata on its first save, so a
  /// preserved copy is the only record of what the SCF started from.
  std::vector<std::string> extra_prefixes;
  /// Molecule to stamp into the RestartMetadata. RestartPlan matches the
  /// archive geometry against the deck at 1e-8 bohr and its eprec exactly;
  /// the molden coordinates differ from the deck at ~1e-8 (DALTON prints
  /// 10 digits) and carry no eprec, so a seed stamped with them was rejected
  /// as "geometry moved" and moldft fell back to the initial guess (closeout
  /// job 2162152, h2o). The caller has already fingerprinted molden vs the
  /// active molecule at 1e-4; stamp the active molecule.
  const madness::Molecule *active_molecule = nullptr;
};

struct GsSeedReport {
  int    n_occ = 0, n_ao = 0, n_mo = 0;
  double max_offdiag_pre = 0.0;   ///< max |S_ij|, i!=j, before Loewdin
  double max_norm_dev    = 0.0;   ///< max | ||phi_i|| - 1 | after
  std::string archive;            ///< <out_prefix>.restartdata
};

/// Write <out_prefix>.restartdata seeded from the molden's lowest n_occ MOs.
/// Collective. Sets FunctionDefaults<3> cell/k/thresh for the projection
/// (callers re-set their own protocol afterwards if it differs).
inline GsSeedReport
write_gs_seed_from_molden(madness::World &world, const std::string &molden_path,
                          int n_occ, const std::string &out_prefix,
                          const GsSeedOptions &opt = {}) {
  using namespace madness;
  using dalton_gs_seed_detail::DaltonMOFunctor;
  using dalton_gs_seed_detail::symbol_to_Z;
  GsSeedReport rep;

  DaltonMoldenResult molden = read_molden(molden_path);
  const int n_ao = molden.n_ao;
  if (n_occ <= 0 || n_occ > molden.n_mo)
    throw std::runtime_error("dalton gs seed: n_occ=" + std::to_string(n_occ) +
                             " out of range (molden n_mo=" + std::to_string(molden.n_mo) + ")");
  rep.n_occ = n_occ; rep.n_ao = n_ao; rep.n_mo = molden.n_mo;

  Molecule molecule;
  for (size_t a = 0; a < molden.atom_symbols.size(); ++a) {
    const int Z = symbol_to_Z(molden.atom_symbols[a]);
    molecule.add_atom(molden.coords[a][0], molden.coords[a][1], molden.coords[a][2],
                      static_cast<double>(Z), Z);
  }

  const int k = default_k_for_thresh(opt.thresh);
  Tensor<double> cell(3L, 2L);
  for (int i = 0; i < 3; i++) { cell(i, 0) = -opt.L; cell(i, 1) = opt.L; }
  FunctionDefaults<3>::set_cell(cell);
  FunctionDefaults<3>::set_k(k);
  FunctionDefaults<3>::set_thresh(opt.thresh);

  std::vector<real_function_3d> amo;
  for (int i = 0; i < n_occ; i++) {
    std::vector<double> col(molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i) * n_ao,
                            molden.mo_coeffs.begin() + static_cast<ptrdiff_t>(i + 1) * n_ao);
    std::shared_ptr<FunctionFunctorInterface<double, 3>> ff =
        std::make_shared<DaltonMOFunctor>(molden.basis, std::move(col));
    amo.push_back(FunctionFactory<double, 3>(world).functor(ff).thresh(opt.thresh)
                      .truncate_on_project());
  }
  {
    Tensor<double> S = matrix_inner(world, amo, amo, true);
    for (int i = 0; i < n_occ; i++)
      for (int j = 0; j < n_occ; j++)
        if (i != j) rep.max_offdiag_pre = std::max(rep.max_offdiag_pre, std::abs(S(i, j)));
  }
  amo = orthonormalize_symmetric(amo);

  Tensor<double> aeps(static_cast<long>(n_occ)), aocc(static_cast<long>(n_occ));
  for (int i = 0; i < n_occ; i++) {
    aeps(i) = (i < static_cast<int>(molden.mo_energies.size())) ? molden.mo_energies[i] : 0.0;
    aocc(i) = 1.0;
    rep.max_norm_dev = std::max(rep.max_norm_dev, std::abs(amo[i].norm2() - 1.0));
  }
  std::vector<int> aset(static_cast<size_t>(n_occ), 0);

  RestartMetadata meta;
  meta.current_energy       = opt.energy;
  meta.spin_restricted      = true;
  meta.L                    = opt.L;
  meta.k                    = k;
  const Molecule &stamp     = opt.active_molecule ? *opt.active_molecule : molecule;
  meta.molecule             = stamp;
  meta.xc                   = opt.xc;
  meta.localize             = opt.localize;
  meta.converged_for_thresh = opt.thresh;
  meta.representation       = Representation::mo;
  meta.eprec                = stamp.parameters.eprec();
  meta.madness_version      = MADNESS_PACKAGE_VERSION;
  // one id for the seed and its preserved copy: they are the same orbitals
  meta.archive_id           = new_archive_id(world);
  meta.origin               = "dalton-seed";
  meta.nalpha               = n_occ;
  meta.nbeta                = n_occ;
  meta.nmo_beta             = n_occ;
  auto emit = [&](auto &ar) {
    meta.write(ar);
    ar & static_cast<unsigned int>(amo.size());
    ar & aeps & aocc & aset;
    for (unsigned int i = 0; i < amo.size(); ++i) ar & amo[i];
  };
  auto write_to = [&](const std::string &prefix) {
    const std::string name = prefix + ".restartdata";
    {
      archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(
          world, name.c_str(), opt.nio);
      emit(ar);
    }
    world.gop.fence();
#ifdef MADNESS_HAS_HDF5
    // HDF5 twin of the seed (opt-in, same switch as the response states), so the
    // DALTON starting point can be visualized next to the MADNESS result.
    if (hdf5_io_enabled())
      save_parallel_archive_hdf5(world, name + ".h5", /*deflate=*/0,
                                 [&](auto &ar) { emit(ar); });
#endif
    return name;
  };
  rep.archive = write_to(out_prefix);
  for (const auto &px : opt.extra_prefixes) write_to(px);

  if (world.rank() == 0)
    print("[DALTON-SEED] GS: wrote", n_occ, "occupied MOs ->", rep.archive,
          " (molden", molden_path, " n_ao =", n_ao, " n_mo =", molden.n_mo,
          " L =", opt.L, " k =", k, " thresh =", opt.thresh,
          " max|S_ij| pre-Loewdin =", rep.max_offdiag_pre, ")");
  return rep;
}

} // namespace molresponse_v3

#endif // MOLRESPONSE_V3_SOLVERS_DALTON_GS_SEED_HPP
