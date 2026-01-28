//
// Created by Florian Bischoff on 19.09.25.
//

#include<ranks_and_hosts.h>

namespace madness {


    double get_rss_usage_in_GB() {
        double kb_to_GB=1.0/(1024*1024);
#ifdef __APPLE__
        double b_to_GB=kb_to_GB/1024;
        struct rusage usage;
        if (getrusage(RUSAGE_SELF, &usage) == -1) {
            std::cerr << "Unable to get RSS usage" << std::endl;
            return -1;
        }
        return usage.ru_maxrss*b_to_GB;
#else
        std::ifstream statm_file("/proc/self/statm");
        unsigned long size, resident;
        if (statm_file.is_open()) {
            statm_file >> size >> resident;
            statm_file.close();
        } else {
            std::cerr << "Unable to open /proc/self/statm" << std::endl;
            return -1;
        }
        long page_size_kb = sysconf(_SC_PAGE_SIZE) / 1024;
        return resident * page_size_kb*kb_to_GB;
#endif
    }

    /// return a mapping rank to hostname, return value on rank 0 only
    std::map<long,std::pair<std::string,double>> rank_to_host_and_rss_map(World& universe) {
        std::vector<std::pair<long,std::pair<std::string,double>>> rank_to_host;
        // rank-local
        auto hostname_and_rss=std::pair<std::string,double>(get_hostname(),get_rss_usage_in_GB());
        rank_to_host.push_back(std::pair<long,std::pair<std::string,double>>(universe.rank(),hostname_and_rss));
        // gather on rank 0
        rank_to_host=universe.gop.concat0(rank_to_host);
        // turn into map
        std::map<long,std::pair<std::string,double>> map;
        for (const auto& [rank,hostname] : rank_to_host) {
            map[rank]=hostname;
        }
        return map;
    }

    /// for each host, return a list of its ranks
    std::map<std::string,std::vector<long>> ranks_per_host(World& universe) {

        std::map<std::string,std::vector<long>> result;
        // Get mapping from rank to (hostname, rss), returns on rank 0 only
        auto rank_to_host = rank_to_host_and_rss_map(universe);
        universe.gop.broadcast_serializable(rank_to_host,0);

        // loop over all ranks and reorder to: hostname -> list of ranks
        for (const auto& [rank, host_rss] : rank_to_host) {
            const std::string& hostname = host_rss.first;
            result[hostname].push_back(rank);
        }
        // Sort the ranks for each host
        for (auto& [hostname, ranks] : result) {
            std::sort(ranks.begin(), ranks.end());
        }
        return result;
    }

    /// for a given rank return the lowest rank on its host
    /// @param[in]  ranks_per_host1 from ranks_per_host()
    long lowest_rank_on_host_of_rank(const std::map<std::string,std::vector<long>> ranks_per_host1, int rank) {
        for (const auto& [hostname,ranks] : ranks_per_host1) {
            for (const auto& r : ranks) {
                if (r==rank) return ranks.front();
            }
        }
        return -1l;
    };

    /**
     * Returns a vector of ranks such that exactly one rank is chosen per host (the lowest rank on each host).
     * Only valid on rank 0. If called on other ranks, returns an empty vector.
     */
    std::vector<int> primary_ranks_per_host(World& world, const std::map<std::string,std::vector<long>>& ranks_per_host1) {
        std::vector<int> result;
        // loop over all ranks, but only keep the lowest rank on each host
        for (int rank = 0; rank < world.size(); ++rank) {
            int lowest_rank = lowest_rank_on_host_of_rank(ranks_per_host1, rank);
            if (lowest_rank == rank) result.push_back(rank);
        }
        std::sort(result.begin(), result.end());
        return result;
    }


    std::string get_hostname() {
        char hostname[256];
#if defined(HAVE_UNISTD_H)
        gethostname(hostname, 256);
#endif
        return std::string(hostname);
    }
}
