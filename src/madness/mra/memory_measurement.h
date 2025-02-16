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
        struct MemKey {
            unsigned long world_id=1;
            unsigned long rank=0;
            std::string hostname="localhost";
            std::size_t DIM=0;
            MemKey(World& world) : world_id(world.id()), rank(world.rank()) {
                char buffer[256];
                gethostname(buffer, 256);
                hostname=std::string(buffer);
            }

            template<typename T, std::size_t NDIM>
            MemKey(const FunctionImpl<T,NDIM>& fimpl) : MemKey(fimpl.world) {
                DIM=NDIM;
            }
            bool operator<(const MemKey& other) const {
                return std::tie(world_id,rank,hostname,DIM) < std::tie(other.world_id,other.rank,other.hostname,other.DIM);
            }
        };

        struct MemInfo {
            long num_functions=0;
            double memory_GB=0.0;
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
            auto all_worlds=World::get_world_ids();
            all_worlds.push_back(World::get_default().id());
            print("searching worlds",all_worlds);
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

        /// measure the memory usage of all objects of all worlds
        void measure_and_print(World& world) {
            world.gop.fence();
            search_all_worlds();
            world.gop.fence();
            print_memory_map(world);
            world.gop.fence();
        }

        /// accumulate the memory usage of all objects of all worlds for this rank per host
        std::map<std::string,double> memory_per_host_local() const {
            std::map<std::string,double> memory_per_host;
            for (const auto& [memkey,memval] : world_memory_map) {
                memory_per_host[memkey.hostname]+=memval.memory_GB;
            }
            return memory_per_host; // different on each rank
        }

        /// accumulate the memory usage of all objects of all worlds over all ranks per host
        std::map<std::string,double> memory_per_host(World& world) const {
            auto memory_per_host=memory_per_host_local();
            for (auto& [host, mem] : memory_per_host)  world.gop.sum(mem);
            return memory_per_host; // same on each rank
        }

        /// return the total memory usage over all hosts
        double total_memory(World& world) const {
            double total_memory=0.0;
            for (const auto& [memkey,memval] : world_memory_map) {
                total_memory+=memval.memory_GB;
            }
            world.gop.sum(total_memory);
            return total_memory;
        }

        /// @param[in] msg a message to print before the memory map
        /// @param[in] world used only for clean printing
        void print_memory_map(World& world, std::string msg="") {
            world.gop.fence();
            if (world.rank()==0) {
                print("final memory map:",msg);
                print("hostname                      world rank  DIM  #funcs  memory_GB");
            }
            std::size_t bufsize=256;
            char line[bufsize];

            // print all information
            world.gop.fence();
            for (const auto& [memkey,memval] : world_memory_map) {
                snprintf(line, bufsize, "%20s %12lu %5lu %5lu %5lu    %e", memkey.hostname.c_str(), memkey.world_id, memkey.rank, memkey.DIM, memval.num_functions, memval.memory_GB);
                print(std::string(line));
            }
            world.gop.fence();

            // print memory on each host and rank
            if (world.rank()==0) print("memory per host and rank");
            std::map<std::string,double> memory_per_host=memory_per_host_local();
            world.gop.fence();
            // turn map into vector
            std::vector<std::pair<std::string,double>> memory_per_host_vec(memory_per_host.begin(),memory_per_host.end());
            // concatenate all vectors, lives on rank 0 of the world
            world.gop.concat0(memory_per_host_vec);
            world.gop.fence();

            if (world.rank()==0) {
                print("hostname               rank   memory_GB");
                for (const auto& [hostname,memory] : memory_per_host_vec) {
                    snprintf(line, bufsize, "%20s %5d    %e", hostname.c_str(), world.rank(), memory);
                    print(std::string(line));
                }
            }
            world.gop.fence();

            // reduce memory_per_host
            // if (world.rank()==0) print("memory per host");
            for (auto& [host, mem] : memory_per_host)  world.gop.sum(mem);
            world.gop.fence();
            if (world.rank()==0) {
                // print("hostname                      memory_GB");
                for (const auto& [hostname,memory] : memory_per_host) {
                    snprintf(line, bufsize, "%20s   all    %e", hostname.c_str(), memory);
                    print(std::string(line));
                }
            }
            world.gop.fence();
            print("RSS usage on rank ",world.rank(),"      ",get_rss_usage_in_GB());
            world.gop.fence();
            if (world.rank()==0) print("");


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
