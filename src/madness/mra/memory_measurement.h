//
// Created by Florian Bischoff on 2/16/25.
//

#ifndef MEMORY_MEASUREMENT_H
#define MEMORY_MEASUREMENT_H


#include<madness/world/ranks_and_hosts.h>
#include<madness/mra/mra.h>
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
#include <madness/mra/stacktrace_util.h>
#endif

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
        /// get the hostname of this machine, rank-local
        static std::string get_hostname() {
            char buffer[256];
            gethostname(buffer, 256);
            return std::string(buffer);
        }

        struct MemKey {
            unsigned long world_id=1;
            unsigned long rank=0;
            std::string hostname="localhost";
            std::size_t DIM=0;
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
            uniqueidT func_id;  ///< Unique ID of the FunctionImpl (avoids hash collisions)
#endif
            MemKey() = default;
            MemKey(World& world) : world_id(world.id()), rank(world.rank()) {
                hostname=MemoryMeasurer::get_hostname();
            }

            template<typename T, std::size_t NDIM>
            MemKey(const FunctionImpl<T,NDIM>& fimpl) : MemKey(fimpl.world) {
                DIM=NDIM;
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
                func_id = fimpl.id();
#endif
            }
            MemKey(const MemKey& other) = default;

            template<typename Archive>
            void serialize(Archive& ar) const {
                ar & world_id & rank & hostname & DIM;
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
                ar & func_id;
#endif
            }
        };
        friend bool operator<(const MemKey& lhs, const MemKey& other) {
            if (lhs.hostname!=other.hostname) return lhs.hostname<other.hostname;
            if (lhs.world_id!=other.world_id) return lhs.world_id<other.world_id;
            if (lhs.rank!=other.rank) return lhs.rank<other.rank;
            if (lhs.DIM!=other.DIM) return lhs.DIM<other.DIM;
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
            return lhs.func_id<other.func_id;
#else
            return false;
#endif
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
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
        /// Map from unique FunctionImpl ID to stacktrace string (no hash collisions)
        std::map<uniqueidT, std::string> stacktrace_strings_;
#endif

        template<typename T, std::size_t NDIM>
        void add_memory_to_map(const FunctionImpl<T,NDIM>& f) {
            const double toGB=double(sizeof(T))/(1024*1024*1024); // convert to GB
            auto sz=f.nCoeff_local();
            if (debug) print("funcimpl<T,",NDIM,"> id",f.id(), "rank",f.world.rank(),"size in GB",sz*toGB);

            // accumulate the sizes into the world_memory_map
            MemKey key(f);
            world_memory_map[key].num_functions++;
            world_memory_map[key].memory_GB+=sz*toGB;
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
            // Store stacktrace string indexed by unique ID (no collisions)
            if (key.func_id) {
                stacktrace_strings_[key.func_id] = f.get_creation_stacktrace();
            }
#endif
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
            all_worlds.push_back(0);  // the default world's ID
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
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
            stacktrace_strings_.clear();
#endif
        }

        /// gather all information of the map on rank 0 of the universe
        void reduce_map(World& universe) {
            // turn map into vector
            std::vector<std::pair<MemKey,MemInfo>> memory_vec(world_memory_map.begin(),world_memory_map.end());
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
            // Gather stacktrace strings across ranks (BEFORE clear_map!)
            std::vector<std::pair<uniqueidT, std::string>> stacktrace_vec(stacktrace_strings_.begin(), stacktrace_strings_.end());
#endif
            // gather all data on rank 0
            memory_vec=universe.gop.concat0(memory_vec);
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
            stacktrace_vec = universe.gop.concat0(stacktrace_vec);
#endif
            // turn back into map
            clear_map();
            for (const auto& [memkey,memval] : memory_vec) {
                world_memory_map[memkey]=memval;
            }
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
            for (const auto& [func_id, trace] : stacktrace_vec) {
                stacktrace_strings_[func_id] = trace;
            }
#endif
        }


        /// given the hostname, return number of ranks and total rss on that node
        static std::map<std::string,std::pair<int,double>> host_to_nrank_and_rss_map(World& universe) {
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
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
                print("hostname                      world rank  DIM  #funcs  memory_GB  func_id");
#else
                print("hostname                      world rank  DIM  #funcs  memory_GB");
#endif
            }
            constexpr std::size_t bufsize=512;
            char line[bufsize];

            // print all information
            world.gop.fence();
            // turn into vector
            std::vector<std::pair<MemKey,MemInfo>> memory_vec(world_memory_map.begin(),world_memory_map.end());
            std::sort(memory_vec.begin(),memory_vec.end(),[](const std::pair<MemKey,MemInfo>& a, const std::pair<MemKey,MemInfo>& b){return a.first<b.first;});
            for (const auto& [memkey,memval] : memory_vec) {
#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
                snprintf(line, bufsize, "%20s %12lu %5lu %5lu %5lu    %e  {%lu,%lu}", memkey.hostname.c_str(), memkey.world_id, memkey.rank, memkey.DIM, memval.num_functions, memval.memory_GB, memkey.func_id.get_world_id(), memkey.func_id.get_obj_id());
#else
                snprintf(line, bufsize, "%20s %12lu %5lu %5lu %5lu    %e", memkey.hostname.c_str(), memkey.world_id, memkey.rank, memkey.DIM, memval.num_functions, memval.memory_GB);
#endif
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

#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
            // Print detailed stacktrace breakdown
            print_stacktrace_details(world);
#endif
        }

#ifdef MADNESS_HAS_MEM_PROFILE_STACKTRACE
        /// Print top allocation sites grouped by stacktrace (aggregated by string, no hash collisions)
        void print_stacktrace_details(World& world) {
            if (world.rank() != 0) return;

            // Group memory by stacktrace string (eliminates hash collision issue)
            std::map<std::string, double> stacktrace_memory;
            std::map<std::string, long> stacktrace_count;

            for (const auto& [memkey, memval] : world_memory_map) {
                if (memkey.func_id) {
                    // Look up the stacktrace for this function's unique ID
                    auto it = stacktrace_strings_.find(memkey.func_id);
                    if (it != stacktrace_strings_.end()) {
                        const std::string& trace = it->second;
                        stacktrace_memory[trace] += memval.memory_GB;
                        stacktrace_count[trace] += memval.num_functions;
                    }
                }
            }

            if (stacktrace_memory.empty()) {
                print("\nNo stacktrace information available");
                return;
            }

            // Convert to vector and sort by memory (descending)
            std::vector<std::pair<std::string, double>> sorted_traces;
            for (const auto& [trace, memory] : stacktrace_memory) {
                sorted_traces.push_back({trace, memory});
            }
            std::sort(sorted_traces.begin(), sorted_traces.end(),
                     [](const auto& a, const auto& b) { return a.second > b.second; });

            // Print top 20 allocation sites
            print("\n========================================");
            print("Top allocation sites by stacktrace:");
            print("========================================");

            const int max_sites = 20;
            int count = 0;
            for (const auto& [trace, memory] : sorted_traces) {
                if (count >= max_sites) break;

                print("\n----------------------------------------");
                print("Rank:", count + 1);
                print("Total memory:", memory, "GB");
                print("Number of functions:", stacktrace_count[trace]);
                print("\nStacktrace:");
                print(trace);

                count++;
            }
            print("\n========================================");
        }
#endif

    };
}

#endif //MEMORY_MEASUREMENT_H
