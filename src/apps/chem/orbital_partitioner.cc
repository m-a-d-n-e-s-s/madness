//
// Created by Florian Bischoff on 12/10/20.
//

#include<chem/orbital_partitioner.h>
#include<cmath>
#include<madness/world/madness_exception.h>

/// partition the orbitals for the exchange operator evenly

/// \param min_ntask_per_world  minimum number of tasks per world for work balancing
/// \param nsubworld            number of subworlds working on the tasks
/// \param nocc                 number of occupied orbitals
/// \param policy               static, guided, akin to OMP scheduling
/// \return                     a vector for pairs<long,long> [start,end) for each batch
std::vector<std::pair<long, long> > OrbitalPartitioner::partition_for_exchange(
        long min_ntask_per_world, long nsubworld, long nocc, std::string policy) {

    std::vector<std::pair<long, long>> ranges;

    /* static
     * Divide the loop into equal-sized chunks or as equal as possible
     *
     * guided
     * Similar to dynamic scheduling, but the chunk size starts off large and decreases to better handle
     * load imbalance between iterations. The optional chunk parameter specifies them minimum size chunk to use.
     * By default the chunk size is approximately loop_count/number_of_threads.
     */

    // split up the exchange matrix in chunks
    if (policy == "static") {
        long ntile_target = min_ntask_per_world * nsubworld;
        long nbatch = long(ceil(sqrt(double(ntile_target))));
        long batchsize = ceil(double(nocc) / double(nbatch));
        ranges.resize(nbatch);
        for (long i = 0; i < nbatch; ++i) {
            ranges[i].first = i * batchsize;
            ranges[i].second = std::min(nocc, (i + 1) * batchsize);
        }
        // make one range large and one range smaller for better work balancing after sorting the tasks
        if (batchsize > 8) {
            ranges[nbatch - 2].second += batchsize / 2;
            ranges[nbatch - 1].first += batchsize / 2;
        }

    } else if (policy == "guided") {
        long begin=0;
        long end=0;
        while (end<nocc) {
            end+=std::max(min_ntask_per_world,((nocc-end)/nsubworld));
            end=std::min(end,nocc);
            ranges.push_back(std::pair<long,long>(begin,end));
            begin=end;
        }


    } else {
        std::string msg="unknown partitioning policy: "+policy;
        MADNESS_EXCEPTION(msg.c_str(),1);
    }

    return ranges;
}
