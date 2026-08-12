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

#include <memory>
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

} // namespace molresponse_v3
