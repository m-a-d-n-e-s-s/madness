//
// Created by Florian Bischoff on 12/10/20.
//

#ifndef MADNESS_ORBITAL_PARTITIONER_H
#define MADNESS_ORBITAL_PARTITIONER_H

#include<vector>

/// partition orbitals into sets for the multiworld/macrotask algorithms
class OrbitalPartitioner {

public:
    static std::vector<std::pair<long,long> > partition_for_exchange(
            long min_ntask_per_world, long nsubworld, long nocc);



};


#endif //MADNESS_ORBITAL_PARTITIONER_H
