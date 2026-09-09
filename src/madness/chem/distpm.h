#ifndef MADNESS_DISTPM_H
#define MADNESS_DISTPM_H

#include<madness/chem/MolecularOrbitals.h>
#include<madness/chem/molecularbasis.h>
#include<madness/chem/molecule.h>
#include <madness/tensor/distributed_matrix.h>


namespace madness {

extern DistributedMatrix<double> distributed_localize_PM(World& world,
                                                         const std::vector<Function<double, 3>>& mo,
                                                         const std::vector<Function<double, 3>>& ao,
                                                         const std::vector<int>& set,
                                                         const std::vector<int>& at_to_bf,
                                                         const std::vector<int>& at_nbf,
                                                         const double thresh = 1e-9,
                                                         const double thetamax = 0.5,
                                                         const bool randomize = true,
                                                         const bool doprint = false);

/// Distributed optimization for the "new" localization method.

/// The "new" objective is Pipek-Mezey in an orthonormal basis, so it is optimized with
/// the same systolic Jacobi sweeps as distributed_localize_PM. C is the replicated
/// (nmo x nao) coefficient matrix already transformed to the orthonormal
/// atomic-eigenfunction basis; at_to_bf/at_nbf hold the per-atom 1s / 2s2p / rest
/// blocks built by Localizer::prepare_new_basis.
extern DistributedMatrix<double> distributed_localize_new(World& world,
                                                          const Tensor<double>& C,
                                                          const std::vector<int>& set,
                                                          const std::vector<int>& at_to_bf,
                                                          const std::vector<int>& at_nbf,
                                                          const double thresh = 1e-9,
                                                          const double thetamax = 0.5);
}

#endif // MADNESS_DISTPM_H
