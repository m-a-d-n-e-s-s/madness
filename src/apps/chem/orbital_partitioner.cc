//
// Created by Florian Bischoff on 12/10/20.
//

#include <chem/orbital_partitioner.h>
#include<math.h>

/// partition the orbitals for the exchange operator evenly

/// \param min_ntask_per_world  minimum number of tasks per world for work balancing
/// \param nsubworld            number of subworlds working on the tasks
/// \param nocc                 number of occupied orbitals
/// \return                     a vector for pairs<long,long> [start,end) for each batch
std::vector<std::pair<long,long> > OrbitalPartitioner::partition_for_exchange(
        long min_ntask_per_world, long nsubworld, long nocc) {

    // split up the exchange matrix in chunks
    long ntile_target=min_ntask_per_world*nsubworld;
    long nbatch=long(ceil(sqrt(double(ntile_target))));
    long batchsize=ceil(double(nocc)/nbatch);
    std::vector<std::pair<long,long> > ranges(nbatch);
    for (long i=0; i<nbatch; ++i) {
        ranges[i].first=i*batchsize;
        ranges[i].second=std::min(nocc,(i+1)*batchsize);
    }

    // make one range large and one range smaller for better work balancing after sorting the tasks
    if (batchsize>8) {
        ranges[nbatch-2].second+=batchsize/2;
        ranges[nbatch-1].first+=batchsize/2;
    }

    return ranges;
}
