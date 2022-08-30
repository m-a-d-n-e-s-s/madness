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
}

#endif // MADNESS_DISTPM_H
