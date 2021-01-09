//
// Created by Florian Bischoff on 12/10/20.
//

#ifndef MADNESS_ORBITAL_PARTITIONER_H
#define MADNESS_ORBITAL_PARTITIONER_H

#include<vector>
#include<string>

/// partition orbitals into sets for the multiworld/macrotask algorithms
class OrbitalPartitioner {

public:
    static std::vector<std::pair<long,long> > partition_for_exchange(
            long min_batch_size, long nsubworld, long nocc, std::string policy="guided");



};


#endif //MADNESS_ORBITAL_PARTITIONER_H
