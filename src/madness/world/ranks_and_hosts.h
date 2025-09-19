//
// Created by Florian Bischoff on 18.09.25.
//

#ifndef MADNESS_RANKS_AND_HOSTS_H
#define MADNESS_RANKS_AND_HOSTS_H

#include <madness/world/world.h>
#include <madness/world/worldgop.h>
#include <madness/misc/misc.h>


/// declare all functions in namespace madness, from ranks_and_hosts.cpp
namespace madness {
    double get_rss_usage_in_GB();
    std::map<long,std::pair<std::string,double>> rank_to_host_and_rss_map(World& universe);
    std::map<std::string,std::vector<long>> ranks_per_host(World& universe);
    long lowest_rank_on_host_of_rank(const std::map<std::string,std::vector<long>> ranks_per_host1, int rank);
    std::vector<int> primary_ranks_per_host(World& world, const std::map<std::string,std::vector<long>> ranks_per_host1);
    std::string get_hostname();

}

#endif //MADNESS_RANKS_AND_HOSTS_H