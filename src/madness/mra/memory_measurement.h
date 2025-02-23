//
// Created by Florian Bischoff on 2/16/25.
//

#ifndef MEMORY_MEASUREMENT_H
#define MEMORY_MEASUREMENT_H


#include<madness/mra/memory_measurement.h>
#include<madness/mra/mra.h>

namespace madness {
    /// measure the memory usage of all FunctionImpl objects of all worlds

    /// Assuming FunctionImpl are the largest objects in a calculation
    /// data is kept as key-value pairs in a map: key=world_id,rank,hostname,NDIM, value=#functions,memory_GB
    class MemoryMeasurer {
    public:

        /// measure the memory usage of all objects of all worlds
        static void measure_and_print(World& world) {
            MemoryMeasurer mm;
            world.gop.fence();
            mm.search_all_worlds();
            world.gop.fence();
            mm.print_memory_map(world);
            world.gop.fence();
            mm.clear_map();
        }

    private:
        struct MemKey {
            unsigned long world_id=1;
            unsigned long rank=0;
            std::string hostname="localhost";
            std::size_t DIM=0;
            MemKey() = default;
            MemKey(World& world) : world_id(world.id()), rank(world.rank()) {
                hostname=MemoryMeasurer::get_hostname();
            }

            template<typename T, std::size_t NDIM>
            MemKey(const FunctionImpl<T,NDIM>& fimpl) : MemKey(fimpl.world) {
                DIM=NDIM;
            }
            MemKey(const MemKey& other) = default;

            template<typename Archive>
            void serialize(Archive& ar) const {
                ar & world_id & rank & hostname & DIM;
            }
        };
        friend bool operator<(const MemKey& lhs, const MemKey& other) {
            if (lhs.hostname!=other.hostname) return lhs.hostname<other.hostname;
            if (lhs.world_id!=other.world_id) return lhs.world_id<other.world_id;
            if (lhs.rank!=other.rank) return lhs.rank<other.rank;
            return lhs.DIM<other.DIM;
        }

        struct MemInfo {
            MemInfo() =default;
            MemInfo(const MemInfo& other) =default;
            long num_functions=0;
            double memory_GB=0.0;
            template <typename Archive>
            void serialize(Archive& ar) const {
                ar & num_functions & memory_GB;
            }
        };

        typedef std::map<MemKey,MemInfo> MemInfoMapT;

        template<typename T, std::size_t NDIM>
        const FunctionImpl<T,NDIM>* cast_to_funcimpl_ptr(const uniqueidT obj_id) {
            World& world=*World::world_from_id(obj_id.get_world_id());
            auto ptr_opt = world.ptr_from_id< WorldObject< FunctionImpl<T,NDIM> > >(obj_id);
            if (!ptr_opt)
                MADNESS_EXCEPTION("FunctionImpl: remote operation attempting to use a locally uninitialized object",0);
            return (dynamic_cast< const FunctionImpl<T,NDIM>*>(*ptr_opt));
        }

        /// keeps track of the memory usage of all objects of one or many worlds **on this rank**
        MemInfoMapT world_memory_map;
        bool debug=false;

        template<typename T, std::size_t NDIM>
        void add_memory_to_map(const FunctionImpl<T,NDIM>& f) {
            const double toGB=double(sizeof(T))/(1024*1024*1024); // convert to GB
            auto sz=f.size_local();
            if (debug) print("funcimpl<T,",NDIM,"> id",f.id(), "rank",f.world.rank(),"size in GB",sz*toGB);

            // accumulate the sizes into the world_memory_map
            world_memory_map[MemKey(f)].num_functions++;
            world_memory_map[MemKey(f)].memory_GB+=sz*toGB;
        }

    public:

        /// add all FunctionImpl<T,NDIM> objects of the given world to the memory map
        /// the memory map is a rank-local object
        void search_world(World& world) {

            auto all_objects=world.get_object_ids();
            if (debug and (world.rank()==0)) print("objects in this world ",all_objects);

            for (const auto& obj : all_objects) {
                if (auto funcimpl=cast_to_funcimpl_ptr<double,1>(obj)) add_memory_to_map(*funcimpl);
                if (auto funcimpl=cast_to_funcimpl_ptr<double,2>(obj)) add_memory_to_map(*funcimpl);
                if (auto funcimpl=cast_to_funcimpl_ptr<double,3>(obj)) add_memory_to_map(*funcimpl);
                if (auto funcimpl=cast_to_funcimpl_ptr<double,4>(obj)) add_memory_to_map(*funcimpl);
                if (auto funcimpl=cast_to_funcimpl_ptr<double,5>(obj)) add_memory_to_map(*funcimpl);
                if (auto funcimpl=cast_to_funcimpl_ptr<double,6>(obj)) add_memory_to_map(*funcimpl);
            }
        }

        /// add all FunctionImpl<T,NDIM> objects **of all worlds** to the memory map
        /// the memory map is a rank-local object
        void search_all_worlds() {
            auto all_worlds=World::get_world_ids(); // all worlds but the default world
            all_worlds.push_back(World::get_default().id());  // add the default world
            if (debug) print("searching worlds",all_worlds);
            for (auto world_id : all_worlds) {
                if (debug) print("searching world",world_id);
                World* thisworld=World::world_from_id(world_id);
                if (World::exists(thisworld)) search_world(*thisworld);
            }
        }

        /// reset the memory map
        void clear_map() {
            world_memory_map.clear();
        }

        /// gather all information of the map on rank 0 of the universe
        void reduce_map(World& universe) {
            // turn map into vector
            std::vector<std::pair<MemKey,MemInfo>> memory_vec(world_memory_map.begin(),world_memory_map.end());
            // gather all data on rank 0
            memory_vec=universe.gop.concat0(memory_vec);
            // turn back into map
            clear_map();
            for (const auto& [memkey,memval] : memory_vec) {
                world_memory_map[memkey]=memval;
            }
        }

        /// get the hostname of this machine, rank-local
        static std::string get_hostname() {
            char buffer[256];
            gethostname(buffer, 256);
            return std::string(buffer);
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

        /// given the hostname, return number of ranks and total rss on that node
        std::map<std::string,std::pair<int,double>> host_to_nrank_and_rss_map(World& universe) {
            auto accumulate_left =[](std::pair<int,double>& a, const std::pair<int,double>& b) {
                a.first++;
                a.second+=b.second;
            };
            auto rank_to_host=rank_to_host_and_rss_map(universe);
            std::map<std::string,std::pair<int,double>> host_to_rank;
            for (const auto& [rank,hostname_and_rss] : rank_to_host) {
                accumulate_left(host_to_rank[hostname_and_rss.first],std::pair<int,double>(rank,hostname_and_rss.second));
            }
            return host_to_rank;
        }


        /// accumulate the memory usage of all objects of all worlds for this rank per host

        /// integrate out world and dim from MemKey, result lives on rank 0 only
        std::vector<std::pair<std::pair<std::string,long>,double>> memory_per_host_and_rank(World& world) const {

            std::map<std::pair<std::string,long>,double> memory_per_host;
            for (const auto& [memkey,memval] : world_memory_map) {
                memory_per_host[{memkey.hostname,memkey.rank}]+=memval.memory_GB;
            }

            // turn map into vector and sort
            std::vector<std::pair<std::pair<std::string,long>,double>> memory_per_host_vec(memory_per_host.begin(),memory_per_host.end());
            std::sort(memory_per_host_vec.begin(),memory_per_host_vec.end(),[](const auto& a, const auto& b){return a.first<b.first;});

            return memory_per_host_vec;
        }

        /// accumulate the memory usage of all objects of all worlds over all ranks per host

        /// integrate out world, dim and rank, only hostname is left
        std::vector<std::pair<std::string,double>> memory_per_host_all_ranks(
            const std::vector<std::pair<std::pair<std::string,long>,double>>& mem_per_host_and_rank) const {
            std::map<std::string,double> mem_per_host;
            for (auto& [hostname_and_rank,memory] : mem_per_host_and_rank) {
                auto hostname=hostname_and_rank.first;
                mem_per_host[hostname]+=memory;
            }
            // turn map into vector
            std::vector<std::pair<std::string,double>> mem_per_host_vec(mem_per_host.begin(),mem_per_host.end());
            return mem_per_host_vec;
        }

        /// return the total memory usage over all hosts
        double total_memory(World& world) const {
            double total_memory=0.0;
            for (const auto& [memkey,memval] : world_memory_map) {
                total_memory+=memval.memory_GB;
            }
            return total_memory;
        }

        /// @param[in] msg a message to print before the memory map
        /// @param[in] world used only for clean printing
        void print_memory_map(World& world, std::string msg="") {
            reduce_map(world);
            world.gop.fence();
            if (world.rank()==0) {
                print("final memory map:",msg);
                print("hostname                      world rank  DIM  #funcs  memory_GB");
            }
            std::size_t bufsize=256;
            char line[bufsize];

            // print all information
            world.gop.fence();
            // turn into vector
            std::vector<std::pair<MemKey,MemInfo>> memory_vec(world_memory_map.begin(),world_memory_map.end());
            std::sort(memory_vec.begin(),memory_vec.end(),[](const std::pair<MemKey,MemInfo>& a, const std::pair<MemKey,MemInfo>& b){return a.first<b.first;});
            for (const auto& [memkey,memval] : memory_vec) {
                snprintf(line, bufsize, "%20s %12lu %5lu %5lu %5lu    %e", memkey.hostname.c_str(), memkey.world_id, memkey.rank, memkey.DIM, memval.num_functions, memval.memory_GB);
                print(std::string(line));
            }
            world.gop.fence();

            // print memory on each host
            auto mem_per_host_and_rank=memory_per_host_and_rank(world);
            auto host_to_nrank_and_rss=host_to_nrank_and_rss_map(world);
            if (world.rank()==0) {
                print("memory per host");
                auto info=memory_per_host_all_ranks(mem_per_host_and_rank);
                print("hostname                      memory_GB     nrank(universe)  rss_GB/host");
                // print("hostname                      memory_GB");
                for (const auto& [hostname,memory] : info) {
                    snprintf(line, bufsize, "%20s          %e       %d           %e", hostname.c_str(), memory,
                    host_to_nrank_and_rss[hostname].first, host_to_nrank_and_rss[hostname].second);
                    print(std::string(line));
                }
            }
            if (world.rank()==0) {
                auto info=memory_per_host_all_ranks(mem_per_host_and_rank);
                double total_mem=total_memory(world);
                double total_rss=0.0;
                for (auto& [hostname,memory] : info) {
                    total_rss+=host_to_nrank_and_rss[hostname].second;
                }
                std::string word="all hosts";
                snprintf(line, bufsize, "%20s          %e       %d           %e",
                    word.c_str(), total_mem, world.size(), total_rss);
                print(std::string(line));
            }

        }

        static double get_rss_usage_in_GB() {
            double kb_to_GB=1.0/(1024*1024);
            double b_to_GB=kb_to_GB/1024;
#ifdef __APPLE__
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

    };
}

#endif //MEMORY_MEASUREMENT_H
