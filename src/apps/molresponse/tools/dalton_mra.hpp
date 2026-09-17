#pragma once

// dalton_mra.hpp — shared DALTON(molden) -> MADNESS MRA projection machinery.
//
// Extracted from tpa_from_dalton.cpp / seed_from_dalton.cpp (which carried
// private copies) so the FD seed path (solvers/dalton_import.hpp) can reuse it
// instead of duplicating a third copy. Pure projection: a Gaussian-LC functor
// over a DaltonMoldenBasis, a single-function projector, and the occ-vir
// block projector
//   x_i(r) = sum_a blk[i,a] * phi_{n_occ+a}(r)
//          = sum_mu (C_vir · blkᵀ)[mu,i] chi_mu(r)
// used for both response vectors (RSPVEC blocks) and plain MO columns.

#include "dalton_gto.hpp"

#include <madness/mra/mra.h>
#include <madness/tensor/tensor_lapack.h>   // syev (Loewdin gauge polish)

#include <algorithm>
#include <cmath>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

namespace molresponse_v3 {

/// Gaussian linear-combination functor: evaluates sum_mu w[mu] chi_mu(r) over
/// a parsed molden basis. special_points() = the shell centers, so the initial
/// refinement finds the nuclei.
class DaltonResponseFunctor : public madness::FunctionFunctorInterface<double, 3> {
  const DaltonMoldenBasis &basis;
  std::vector<double> weights;  // AO weights for this function
  std::vector<madness::coord_3d> centers;

public:
  DaltonResponseFunctor(const DaltonMoldenBasis &b, std::vector<double> w)
      : basis(b), weights(std::move(w)) {
    for (const auto &sh : basis.shells) centers.push_back({sh.cx, sh.cy, sh.cz});
  }
  double operator()(const madness::coord_3d &r) const override {
    double val = 0.0;
    double bf[9];  // max components for a g shell
    for (size_t s = 0; s < basis.shells.size(); s++) {
      const auto &sh = basis.shells[s];
      sh.evaluate(r[0], r[1], r[2], bf);
      const int off = basis.ao_offsets[s];
      for (int k = 0; k < sh.n_ao; k++)
        val += weights[static_cast<size_t>(off + k)] * bf[k];
    }
    return val;
  }
  std::vector<madness::coord_3d> special_points() const override { return centers; }
};

/// Project one AO-weight vector to an MRA function at `thresh` (active k/cell).
inline madness::real_function_3d
project_dalton_weights(madness::World &world, const DaltonMoldenBasis &basis,
                       std::vector<double> w, double thresh) {
  std::shared_ptr<madness::FunctionFunctorInterface<double, 3>> f =
      std::make_shared<DaltonResponseFunctor>(basis, std::move(w));
  return madness::FunctionFactory<double, 3>(world)
      .functor(f).thresh(thresh).truncate_on_project();
}

/// Project a flat (occ-outer, vir-inner) coefficient block into n_occ MRA
/// functions: out[i] = scale * sum_a blk[i,a] phi_{n_occ+a}. `C` is the
/// column-major molden MO matrix C[mu + n_ao*mo].
inline std::vector<madness::real_function_3d>
project_dalton_ov_block(madness::World &world, const DaltonMoldenBasis &basis,
                        const std::vector<double> &C, int n_ao, int n_mo,
                        int n_occ, int n_vir, const std::vector<double> &blk,
                        double thresh, double scale) {
  (void)n_mo;
  std::vector<madness::real_function_3d> out;
  for (int i = 0; i < n_occ; ++i) {
    std::vector<double> w(static_cast<size_t>(n_ao), 0.0);
    for (int a = 0; a < n_vir; ++a) {
      const double coef = blk[static_cast<size_t>(i) * n_vir + a] * scale;
      if (coef == 0.0) continue;
      const double *col = &C[static_cast<size_t>(n_ao) * (n_occ + a)];
      for (int mu = 0; mu < n_ao; ++mu) w[static_cast<size_t>(mu)] += col[mu] * coef;
    }
    out.push_back(project_dalton_weights(world, basis, std::move(w), thresh));
  }
  madness::truncate(world, out);
  return out;
}

/// Occupied-orbital GAUGE ROTATION (shared by the FD dalton.dir import,
/// solvers/dalton_import.hpp, and the ES seed tool, tools/seed_from_dalton).
///
/// DALTON's response/excitation vectors are indexed by its CANONICAL occupied
/// MOs; the MADNESS ground-state orbitals are typically LOCALIZED (moldft
/// `localize new`). Response functions transform covariantly with the occupied
/// orbitals, so without this the per-orbital pairing is scrambled and the seed
/// is worthless (observed on the h2o FD A/B: seed-implied alpha_zz 3.56 vs
/// DALTON's 8.39, zero iteration savings). Build the occupied-occupied overlap
/// M(i,j) = <phi^DAL_i | phi^MAD_j> in MRA space, Loewdin-polish it to an
/// exact unitary U = M (M^T M)^{-1/2}, and rotate every seed block:
///     x^MAD_j = sum_i x^DAL_i U(i,j)     (same for y),
/// i.e. x = transform(world, x, U), applied BEFORE any Q-projection.
///
/// The eigenvalues of M^T M are ALSO the quantitative ground-state import
/// fidelity: ~1 when the two occupied spaces coincide. An eigenvalue below
/// 0.5 means a MADNESS occupied orbital has little support in the DALTON
/// occupied space — a different state entirely — hard error.
inline madness::Tensor<double>
occupied_gauge_rotation(madness::World &world,
                        const std::vector<madness::real_function_3d> &phi_dal_mra,
                        const std::vector<madness::real_function_3d> &phi_mad,
                        bool verbose) {
  using madness::Tensor;
  Tensor<double> M = matrix_inner(world, phi_dal_mra, phi_mad);
  Tensor<double> S = inner(transpose(M), M);   // M^T M, n_occ x n_occ
  Tensor<double> V, ev;
  syev(S, V, ev);
  double ev_min = ev(0L);
  for (long k = 0; k < ev.dim(0); ++k) ev_min = std::min(ev_min, ev(k));
  if (verbose && world.rank() == 0)
    madness::print("[DALTON-SEED] occupied-gauge overlap M^T M eigenvalues:", ev,
                   "  (import fidelity; 1 = perfect span match)");
  if (ev_min < 0.5) {
    std::ostringstream os;
    os << "dalton import: occupied-space mismatch — min eigenvalue of the "
          "DALTON/MADNESS occupied overlap M^T M is " << ev_min
       << " (want ~1). The DALTON ground state does not span the MADNESS "
          "occupied orbitals (different state/charge/geometry?)";
    throw std::runtime_error(os.str());
  }
  const long n = ev.dim(0);
  Tensor<double> Sinvhalf(n, n);
  for (long a = 0; a < n; ++a)
    for (long b = 0; b < n; ++b) {
      double s = 0.0;
      for (long k = 0; k < n; ++k)
        s += V(a, k) * V(b, k) / std::sqrt(ev(k));
      Sinvhalf(a, b) = s;
    }
  return inner(M, Sinvhalf);   // Loewdin: exactly unitary
}

} // namespace molresponse_v3
