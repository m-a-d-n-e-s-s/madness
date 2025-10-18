
/**
 \file cloud.h
 \brief Declares the \c Cloud class for storing data and transfering them between worlds
 \ingroup world

*/

/**
 * TODO:  - delete container record upon caching if container is replicated
 */

#ifndef SRC_MADNESS_WORLD_CLOUD_H_
#define SRC_MADNESS_WORLD_CLOUD_H_


#include <madness/world/parallel_dc_archive.h>
#include<any>
#include<iomanip>


/*!
  \file cloud.h
  \brief Defines and implements most of madness cloud storage

  TODO: check use of preprocessor directives
  TODO: clear cache in destructor won't work because no subworld is present -> must be explicitly called, error prone/
*/

namespace madness {

    /// \brief A utility to get the name of a type as a string from chatGPT
    template<typename T>
    struct type_name {
        static const char* value() { return typeid(T).name();}
    };

    template<>
    struct type_name<Function<double,1>> { static const char* value() { return "Function<double,1>"; } };
    template<>
    struct type_name<Function<double,2>> { static const char* value() { return "Function<double,2>"; } };
    template<>
    struct type_name<Function<double,3>> { static const char* value() { return "Function<double,3>"; } };
    template<>
    struct type_name<Function<double,4>> { static const char* value() { return "Function<double,4>"; } };
    template<>
    struct type_name<Function<double,5>> { static const char* value() { return "Function<double,5>"; } };
    template<>
    struct type_name<Function<double,6>> { static const char* value() { return "Function<double,6>"; } };

    template<>
    struct type_name<std::vector<Function<double,1>>> { static const char* value() { return "std::vector<Function<double,1>>"; } };
    template<>
    struct type_name<std::vector<Function<double,2>>> { static const char* value() { return "std::vector<Function<double,2>>"; } };
    template<>
    struct type_name<std::vector<Function<double,3>>> { static const char* value() { return "std::vector<Function<double,3>>"; } };
    template<>
    struct type_name<std::vector<Function<double,4>>> { static const char* value() { return "std::vector<Function<double,4>>"; } };
    template<>
    struct type_name<std::vector<Function<double,5>>> { static const char* value() { return "std::vector<Function<double,5>>"; } };
    template<>
    struct type_name<std::vector<Function<double,6>>> { static const char* value() { return "std::vector<Function<double,6>>"; } };

template<typename keyT>
struct Recordlist {
    std::list<keyT> list;

    Recordlist() : list() {};

    explicit Recordlist(const keyT &key) : list{key} {};

    Recordlist(const Recordlist &other) : list(other.list) {};

    Recordlist &operator+=(const Recordlist &list2) {
        for (auto &l2 : list2.list) list.push_back(l2);
        return *this;
    }

    Recordlist &operator+=(const keyT &key) {
        list.push_back(key);
        return *this;
    }

    keyT pop_front_and_return() {
        keyT key = list.front();
        list.pop_front();
        return key;
    }

    std::size_t size() const {
        return list.size();
    }

    // if type provides id() member function (i.e. WorldObject) use that for hashing, otherwise use hash_value() for
    // fundamental types (see worldhash.h)
    template <typename T>
    using member_id_t = decltype(std::declval<T>().id());

    template <typename T>
    using has_member_id = madness::meta::is_detected<member_id_t, T>;

    // if type provides a hashing function use that, intrusive hashing, see worldhash.h
    template <typename T>
    using member_hash_t = decltype(std::declval<T>().hash());

    template <typename T>
    using has_member_hash = madness::meta::is_detected<member_hash_t, T>;

    template<typename T, std::size_t NDIM>
    static keyT compute_record(const Function<T,NDIM>& arg) {return hash_value(arg.get_impl()->id());}

    template<typename T, std::size_t NDIM>
    static keyT compute_record(const FunctionImpl<T,NDIM>* arg) {return hash_value(arg->id());}

    template<typename keyQ, typename valueT>
    static keyT compute_record(const WorldContainer<keyQ,valueT>& arg) {return hash_value(arg.id());}

    template<typename keyQ, typename valueT>
    static keyT compute_record(const std::shared_ptr<WorldContainer<keyQ,valueT>>& arg) {return hash_value(arg->id());}

    template<typename T, std::size_t NDIM>
    static keyT compute_record(const std::shared_ptr<madness::FunctionImpl<T, NDIM>>& arg) {return hash_value(arg->id());}

    template<typename T>
    static keyT compute_record(const std::vector<T>& arg) {return hash_range(arg.begin(), arg.end());}

    template<typename T>
    static keyT compute_record(const Tensor<T>& arg) {return hash_value(arg.normf());}

    template<typename T>
    static keyT compute_record(const std::shared_ptr<T>& arg) {return compute_record(*arg);}

    template<typename T>
    static keyT compute_record(const T& arg) {
        if constexpr (has_member_id<T>::value) {
            return hash_value(arg.id());
        } else if constexpr (std::is_pointer_v<T> && has_member_id<std::remove_pointer_t<T>>::value) {
            return hash_value(arg->id());
        } else {
            // compute hash_code for fundamental types
            std::size_t hashtype = typeid(T).hash_code();
            hash_combine(hashtype,hash_value(arg));
            return hashtype;
        }
    }


    friend std::ostream &operator<<(std::ostream &os, const Recordlist &arg) {
        using namespace madness::operators;
        os << arg.list;
        return os;
    }

};

/// cloud class

/// store and load data to/from the cloud into arbitrary worlds
///
/// Distributed data is always bound to a certain world. If it needs to be
/// present in another world it can be serialized to the cloud and deserialized
/// from there again. For an example see test_cloud.cc
///
/// Data is stored into a distributed container living in the universe.
/// During storing a (replicated) list of records is returned that can be used to find the data
/// in the container. If a combined object (a vector, tuple, etc) is stored a list of records
/// will be generated. When loading the data from the world the record list will be used to
/// deserialize all stored objects.
///
/// Note that there must be a fence after the destruction of subworld containers, as in:
///
///  create subworlds
///  {
///      dcT(subworld)
///      do work
///  }
///  subworld.gop.fence();
class Cloud {

    bool debug = false;       ///< prints debug output
    bool is_replicated=false;   ///< if contents of the container are replicated
    bool dofence = true;      ///< fences after load/store
    bool force_load_from_cache = false;       ///< forces load from cache (mainly for debugging)
    bool use_cache=true;

public:

    typedef std::any cached_objT;
    using keyT = madness::archive::ContainerRecordOutputArchive::keyT;
    using valueT = std::vector<unsigned char>;
    typedef std::map<keyT, cached_objT> cacheT;
    typedef Recordlist<keyT> recordlistT;

    enum StoragePolicy {
        StoreFunction,          ///< store a madness function in the cloud  -- can have a large memory impact
                                ///< equivalent to a deep copy
        StoreFunctionPointer,   ///< store the pointer to the function in the cloud.
                                ///< Return type still is a Function<T,NDIM> with a pointer to the universe function impl.
                                ///< equivalent to a shallow copy
    };


    friend std::ostream& operator<<(std::ostream& os, const StoragePolicy& sp) {
        switch(sp) {
            case StoreFunction: os << "StoreFunction"; break;
            case StoreFunctionPointer: os << "StoreFunctionPointer"; break;
            default: os << "UnknownStoragePolicy"; break;
        }
        return os;
    }

    friend std::string to_string(const StoragePolicy sp) {
        std::ostringstream os;
        os << sp;
        return os.str();
    }

private:
    /// are the functions (WorldObjects) stored in the cloud or only pointers to them
    StoragePolicy storage_policy = StoreFunctionPointer;

    /// cloud is a container: replication policy for the cloud container: distributed, node-replicated, rank-replicated
    DistributionType cloud_replication_policy = Distributed;

    mutable madness::WorldContainer<keyT, valueT> container;
    cacheT cached_objects;
    recordlistT local_list_of_container_keys;   // a world-local list of keys occupied in container

public:
    std::list<WorldObjectBase*> world_object_base_list; // list of world objects stored in the cloud

    template <typename T>
    using member_cloud_serialize_t = decltype(std::declval<T>().cloud_store(std::declval<World&>(), std::declval<Cloud&>()));

    template <typename T>
    using has_cloud_serialize = madness::meta::is_detected<member_cloud_serialize_t, T>;

public:

    /// @param[in]	universe	the universe world
    Cloud(madness::World &universe) : container(universe), reading_time(0l), copy_time(0l), writing_time(0l),
        cache_reads(0l), cache_stores(0l) {
    }

    ~Cloud() {
        if ((not cached_objects.empty()) or (not local_list_of_container_keys.list.empty())) {
            print("\nCloud::~Cloud(): cached_objects not empty, size=", cached_objects.size());
            print("You need to call clear_cache(subworld) before destroying the cloud");
            print("\n------------------------------\n");
            std::string msg="deferred destruction of cloud with non-empty cache";
            std::cerr << msg << std::endl;
        }
    }

    void set_debug(bool value) {
        debug = value;
    }

    void set_fence(bool value) {
        dofence = value;
    }

    void set_force_load_from_cache(bool value) {
        force_load_from_cache = value;
    }

    /// is the cloud container replicated: per rank, per node, or distributed
    void set_replication_policy(const DistributionType value) {
        cloud_replication_policy = value;
        use_cache=false;
        if (value == RankReplicated) use_cache=true;
    }

    /// is the cloud container replicated: per rank, per node, or distributed
    DistributionType get_replication_policy() const {
        return cloud_replication_policy;
    }

    bool validate_replication_policy() const {
        auto disttype=validate_distribution_type(container);
        if (disttype!=cloud_replication_policy) {
            std::cout << "Cloud::validate_distribution(): distribution type mismatch, container is " << disttype
                      << " but cloud_replication_policy is " << cloud_replication_policy << std::endl;
            return false;
        }
        return true;
    }


    /// storing policy refers to storing functions or pointers to functions
    void set_storing_policy(const StoragePolicy value) {
        storage_policy = value;
    }

    /// storing policy refers to storing functions or pointers to functions
    StoragePolicy get_storing_policy() const {
        return storage_policy;
    }

    void print_size(World& universe) {
        nlohmann::json stats=gather_memory_statistics(universe);
        double byte2gbyte=1.0/(1024*1024*1024);
        double global_memsize=stats["memory_global"].template get<double>();
        double max_record_size=stats["max_record_size"].template get<double>();
        double min_memsize=stats["memory_min"].template get<double>();
        double max_memsize=stats["memory_max"].template get<double>();
        double global_size=stats["container_size_global"].template get<double>();

        if (universe.rank()==0) {
            print("Cloud memory:");
            print("  replicated:",is_replicated);
            print("size of cloud (total)");
            print("  number of records:        ",global_size);
            print("  memory in GBytes:         ",global_memsize*byte2gbyte);
            print("size of cloud (average per node)");
            print("  number of records:        ",double(global_size)/universe.size());
            print("  memory in GBytes:         ",global_memsize*byte2gbyte/universe.size());
            print("min/max of node");
            print("  memory in GBytes:         ",min_memsize*byte2gbyte,max_memsize*byte2gbyte);
            print("  max record size in GBytes:",max_record_size*byte2gbyte);

        }
    }

    /// return a json object with the cloud settings and statistics
    nlohmann::json get_statistics(World& world) const {
        nlohmann::json j;
        {   // settings
            j["storage_policy"]=to_string(storage_policy);
            j["cloud_replication_policy"]=to_string(cloud_replication_policy);
            j["is_replicated"]=is_replicated;
            j["local_cached_objects_size"]=cached_objects.size();
        }
        // timings
        j.update(gather_timings(world));
        j.update(gather_memory_statistics(world));
        return j;

    }

    /// get size of the cloud container
    nlohmann::json gather_memory_statistics(World &universe) const {

        std::size_t memsize=0;
        std::size_t max_record_size=0;
        for (auto& item : container) {
            memsize+=item.second.size();
            max_record_size=std::max(max_record_size,item.second.size());
        }
        std::size_t global_memsize=memsize;
        std::size_t max_memsize=memsize;
        std::size_t min_memsize=memsize;
        double rss=madness::get_rss_usage_in_GB();
        universe.gop.sum(global_memsize);
        universe.gop.max(max_memsize);
        universe.gop.max(max_record_size);
        universe.gop.min(min_memsize);
        universe.gop.max(rss);

        auto local_size=container.size();
        auto global_size=local_size;
        universe.gop.sum(global_size);
        nlohmann::json j;
        j["container_size_global"] = global_size;
        j["memory_global"] = global_memsize;
        j["memory_min"] = min_memsize;
        j["memory_max"] = max_memsize;
        j["memory_rss_GB_max"] = rss;
        j["max_record_size"] = max_record_size;
        return j;
    }

    nlohmann::json gather_timings(World &universe) const {
        double rtime_max = double(reading_time)*1.e-6;
        double rtime_acc = double(reading_time)*1.e-6;
        double rtime_av = double(reading_time)*1.e-6;
        double ctime_max = double(copy_time)*1.e-6;
        double ctime_acc = double(copy_time)*1.e-6;
        double ctime_av = double(copy_time)*1.e-6;
        double wtime = double(writing_time)*1.e-6;
        double ptime = double(replication_time)*1.e-6;
        double tptime = double(target_replication_time)*1.e-6;
        universe.gop.max(rtime_max);
        universe.gop.sum(rtime_acc);
        rtime_av = rtime_acc/universe.size();
        universe.gop.max(ctime_max);
        universe.gop.sum(ctime_acc);
        ctime_av = ctime_acc/universe.size();
        universe.gop.max(wtime);
        universe.gop.max(ptime);
        universe.gop.max(tptime);
        long creads = long(cache_reads);
        long cstores = long(cache_stores);
        universe.gop.sum(creads);
        universe.gop.sum(cstores);
        nlohmann::json j;
        j["reading_time_max_s"] = rtime_max;
        j["reading_time_acc_s"] = rtime_acc;
        j["reading_time_av_s"] = rtime_av;
        j["copy_time_max_s"] = ctime_max;
        j["copy_time_acc_s"] = ctime_acc;
        j["copy_time_av_s"] = ctime_av;
        j["writing_time_s"] = wtime;
        j["replication_time_s"] = ptime;
        j["target_replication_time_s"] = tptime;
        j["cache_reads"] = creads;
        j["cache_stores"] = cstores;
        return j;
    }

    /// backwards compatibility
    void print_timings(World& universe) const {
        print_timings(gather_timings(universe));
    }

    static void print_timings(const nlohmann::json timings) {
        double rtime_max=timings["reading_time_max_s"].template get<double>();
        double rtime_av=timings["reading_time_av_s"].template get<double>();
        double rtime_acc=timings["reading_time_acc_s"].template get<double>();
        // double ctime_max=timings["copy_time_max_s"].template get<double>();
        // double ctime_av=timings["copy_time_av_s"].template get<double>();
        // double ctime_acc=timings["copy_time_acc_s"].template get<double>();
        double wtime=timings["writing_time_s"].template get<double>();
        double ptime=timings["replication_time_s"].template get<double>();
        double tptime=timings["target_replication_time_s"].template get<double>();
        long creads=timings["cache_reads"].template get<long>();
        long cstores=timings["cache_stores"].template get<long>();

        auto precision = std::cout.precision();
        std::cout << std::fixed << std::setprecision(1);
        print("cloud storing wall time                        ", wtime);
        print("cloud replication wall time                    ", ptime);
        print("target replication wall time                   ", tptime);
        print("cloud max reading time (all procs)             ", rtime_max, std::defaultfloat);
        print("cloud average reading cpu time (all procs)     ", rtime_av, std::defaultfloat);
        print("cloud accumulated reading cpu time (all procs) ", rtime_acc, std::defaultfloat);
        std::cout << std::setprecision(precision) << std::scientific;
        print("cloud cache stores                             ", long(cstores));
        print("cloud cache loads                              ", long(creads));
    }

    static void print_memory_statistics(const nlohmann::json stats) {
        double byte2gbyte=1.0/(1024*1024*1024);
        double global_memsize=stats["memory_global"].template get<double>();
        double max_record_size=stats["max_record_size"].template get<double>();
        double min_memsize=stats["memory_min"].template get<double>();
        double max_memsize=stats["memory_max"].template get<double>();
        double global_size=stats["container_size_global"].template get<double>();

        print("Cloud memory:");
        print("  size of cloud (total)");
        print("    number of records:        ",global_size);
        print("    memory in GBytes:         ",global_memsize*byte2gbyte);
        // print("  size of cloud (average per node)");
        // print("    number of records:        ",double(global_size)/madness::world().size());
        // print("    memory in GBytes:         ",global_memsize*byte2gbyte/madness::world().size());
        print("  min/max of node");
        print("    memory in GBytes:         ",min_memsize*byte2gbyte,max_memsize*byte2gbyte);
        print("    max record size in GBytes:",max_record_size*byte2gbyte);
    }

    void clear_cache(World &subworld) {
        cached_objects.clear();
        local_list_of_container_keys.list.clear();
        subworld.gop.fence();
    }

    void clear() {
        container.clear();
    }

    void clear_timings() {
        reading_time=0l;
        copy_time=0l;
        writing_time=0l;
        writing_time1=0l;
        replication_time=0l;
        target_replication_time=0l;
        cache_stores=0l;
        cache_reads=0l;
    }

    /// functor to distribute/rank/node-replicate a function, passed in as a pointer to WorldObjectBase
    template<typename T, std::size_t NDIM>
    struct DistributeFunctor {
        DistributionType dt= Distributed;
        explicit DistributeFunctor(const DistributionType dt) : dt(dt) {}
        int operator()(WorldObjectBase* wo) const {
            // figure out if wo is a FunctionImpl and do the distribution
            if (auto fimpl=dynamic_cast<FunctionImpl<T, NDIM>*>(wo)) {
                // fimpl->get_pmap()->print_data_sizes(world,"before distribution of function in cloud");
                if (dt==RankReplicated) {
                    fimpl->replicate(false);
                } else if (dt==NodeReplicated) {
                    print("replicating function per node",fimpl);;
                    // fimpl->print_size("fimpl");
                    fimpl->replicate_on_hosts(true);
                } else if (dt==Distributed) {
                    fimpl->undo_replicate(false);
                } else {
                    MADNESS_EXCEPTION("unknown distribution type",1);
                }
                // fimpl->get_pmap()->print_data_sizes(world,"after distribution of function in cloud");
            }
            return 0;
        }
    };

    /// distribute/node/rank replicate the targets of all world objects stored in the cloud
    void distribute_targets(const DistributionType dt= Distributed) {
        if (world_object_base_list.empty()) return;

        World& world=world_object_base_list.front()->get_world();
        for (auto wo : world_object_base_list) {
            loop_types<DistributeFunctor, double, float, double_complex, float_complex>(std::tuple<DistributionType>(dt),wo);
            world.gop.fence();
        }
        world.gop.fence();

    }

    /// @param[in]  world the subworld the objects are loaded to
    /// @param[in]  recordlist the list of records where the objects are stored

    /// load a single object from the cloud, recordlist is kept unchanged
    template<typename T>
    T load(madness::World &world, const recordlistT recordlist) const {
        recordlistT rlist = recordlist;
        cloudtimer t(world, reading_time);

        // forward_load will consume the recordlist while loading elements
        return forward_load<T>(world, rlist);
    }

    /// similar to load, but will consume the recordlist

    /// @param[in]  world the subworld the objects are loaded to
    /// @param[in]  recordlist the list of records where the objects are stored
    template<typename T>
    T consuming_load(madness::World &world, recordlistT& recordlist) const {
        cloudtimer t(world, reading_time);

        // forward_load will consume the recordlist while loading elements
        return forward_load<T>(world, recordlist);
    }

    /// load a single object from the cloud, recordlist is consumed while loading elements
    template<typename T>
    T forward_load(madness::World &world, recordlistT& recordlist) const {
        // different objects are stored in different ways
        // - tuples are split up into their components
        // - classes with their own cloud serialization are stored using that
        // - everything else is stored using their usual serialization
        if constexpr (is_tuple<T>::value) {
            return load_tuple<T>(world, recordlist);
        } else if constexpr (has_cloud_serialize<T>::value) {
            T target = allocator<T>(world);
            target.cloud_load(world, *this, recordlist);
            return target;
        } else {
            return do_load<T>(world, recordlist);
        }
    }

    /// @param[in]  world presumably the universe
    template<typename T>
    recordlistT store(madness::World &world, const T &source) {
        if (is_replicated) {
            print("Cloud contents are replicated and read-only!");
            MADNESS_EXCEPTION("cloud error",1);
        }
        cloudtimer t(world,writing_time);

        // different objects are stored in different ways
        // - tuples are split up into their components
        // - classes with their own cloud serialization are stored using that
        // - everything else is stored using their usual serialization
        recordlistT recordlist;
        if constexpr (is_tuple<T>::value) {
            recordlist+=store_tuple(world,source);
        } else if constexpr (has_cloud_serialize<T>::value) {
            recordlist+=source.cloud_store(world,*this);
        } else {
            recordlist+=store_other(world,source);
        }
        if (dofence) world.gop.fence();
        return recordlist;
    }

    void replicate_according_to_policy(const std::size_t chunk_size=INT_MAX) {
        if (cloud_replication_policy == Distributed) {
            // if (debug and (container.size() > 0)) print("no replication of container");
            return;
        }
        else if (cloud_replication_policy == RankReplicated) {
            replicate(chunk_size);
        }
        else if (cloud_replication_policy == NodeReplicated) {
            replicate_per_node(chunk_size);
        }
        else {
            MADNESS_EXCEPTION("unknown replication policy",1);
        }
        container.get_world().gop.fence();
    }

    void replicate_per_node(const std::size_t chunk_size=INT_MAX) {
        // this will fail if the container values are larger that 2GB
        // need to reimplement that at some point
        try {
            double cpu0=cpu_time();
            World& world=container.get_world();
            world.gop.fence();
            cloudtimer t(world,replication_time);
            MADNESS_CHECK_THROW(not is_replicated,"cloud::replicate_per_node: container is already replicated");
            container.replicate_on_hosts(true);
            is_replicated=true;
            world.gop.fence();
            double cpu1=cpu_time();
            if (debug and (world.rank()==0)) print("replication_per_node ended after ",cpu1-cpu0," seconds");
        } catch (...) {
            MADNESS_EXCEPTION("cloud replication_per_node failed, presumably because some data is larger than 2GB",1);
        }
    }

    // replicates the contents of the container
    void replicate(const std::size_t chunk_size=INT_MAX) {
        MADNESS_CHECK_THROW(not is_replicated,"cloud::replicate_per_node: container is already replicated");

        double cpu0=cpu_time();
        World& world=container.get_world();
        world.gop.fence();
        cloudtimer t(world,replication_time);
        container.reset_pmap_to_local();
        is_replicated=true;

        std::list<keyT> keylist;
        for (auto it=container.begin(); it!=container.end(); ++it) {
            keylist.push_back(it->first);
        }

        for (ProcessID rank=0; rank<world.size(); rank++) {
            if (rank == world.rank()) {
                std::size_t keylistsize = keylist.size();
                world.mpi.Bcast(&keylistsize,sizeof(keylistsize),MPI_BYTE,rank);

                for (auto key : keylist) {
                    madness::WorldContainer<keyT, std::vector<unsigned char> >::const_accessor acc;
                    bool found=container.find(acc,key);
                    MADNESS_CHECK(found);
                    auto data = acc->second;
                    std::size_t sz=data.size();

                    world.mpi.Bcast(&key,sizeof(key),MPI_BYTE,rank);
                    world.mpi.Bcast(&sz,sizeof(sz),MPI_BYTE,rank);

                    // if data is too large for MPI_INT break it into pieces to avoid integer overflow
                    for (std::size_t start=0; start<sz; start+=chunk_size) {
                        std::size_t remainder = std::min(sz - start, chunk_size);
                        world.mpi.Bcast(&data[start], remainder, MPI_BYTE, rank);
                    }

                }
            }
            else {
                std::size_t keylistsize;
                world.mpi.Bcast(&keylistsize,sizeof(keylistsize),MPI_BYTE,rank);
                for (size_t i=0; i<keylistsize; i++) {
                    keyT key;
                    world.mpi.Bcast(&key,sizeof(key),MPI_BYTE,rank);
                    std::size_t sz;
                    world.mpi.Bcast(&sz,sizeof(sz),MPI_BYTE,rank);
                    valueT data(sz);
//                    world.mpi.Bcast(&data[0],sz,MPI_BYTE,rank);
                    for (std::size_t start=0; start<sz; start+=chunk_size) {
                        std::size_t remainder=std::min(sz-start,chunk_size);
                        world.mpi.Bcast(&data[start],remainder,MPI_BYTE,rank);
                    }

                    container.replace(key,data);
                }
            }
        }
        world.gop.fence();
        double cpu1=cpu_time();
        if (debug and (world.rank()==0)) print("replication ended after ",cpu1-cpu0," seconds");
    }

private:

    mutable std::atomic<long> reading_time=0l;     // in microseconds
public:
    mutable std::atomic<long> copy_time=0l;        // if pointers are stored in cloud, time to copy from universe to subworld
    mutable std::atomic<long> target_replication_time=0l;     // if pointers are stored in cloud, time to replicate targets
private:
    mutable std::atomic<long> writing_time=0l;     // in microseconds
    mutable std::atomic<long> writing_time1=0l;    // in microseconds
    mutable std::atomic<long> replication_time=0l;    // in microseconds
    mutable std::atomic<long> cache_reads=0l;
    mutable std::atomic<long> cache_stores=0l;

    template<typename> struct is_tuple : std::false_type { };
    template<typename ...T> struct is_tuple<std::tuple<T...>> : std::true_type { };

    template<typename Q> struct is_vector : std::false_type { };
    template<typename Q> struct is_vector<std::vector<Q>> : std::true_type { };

    template<typename T> using is_parallel_serializable_object = std::is_base_of<archive::ParallelSerializableObject,T>;

    template<typename T> using is_world_constructible = std::is_constructible<T, World &>;

public:
    struct cloudtimer {
        World& world;
        double wall0;
        std::atomic<long> &rtime;

        cloudtimer(World& world, std::atomic<long> &readtime) : world(world), wall0(wall_time()), rtime(readtime) {}

        ~cloudtimer() {
            long deltatime=long((wall_time() - wall0) * 1000000l);
            rtime += deltatime;
        }
    };
private:

    template<typename T>
    void cache(madness::World &world, const T &obj, const keyT &record) const {
        const_cast<cacheT &>(cached_objects).insert({record,std::make_any<T>(obj)});
    }

    /// load an object from the cache, record is unchanged
    template<typename T>
    T load_from_cache(madness::World &world, const keyT &record) const {
        if (world.rank()==0) cache_reads++;
        if (debug) print("loading", type_name<T>::value(), "from cache record", record, "to world", world.id());
        if (auto obj = std::any_cast<T>(&cached_objects.find(record)->second)) return *obj;
        MADNESS_EXCEPTION("failed to load from cloud-cache", 1);
        T target = allocator<T>(world);
        return target;
    }

    bool is_cached(const keyT &key) const {
        return (cached_objects.count(key) == 1);
    }

    /// checks if a (universe) container record is used

    /// currently implemented with a local copy of the recordlist, might be
    /// reimplemented with container.find(), which would include blocking communication.
    bool is_in_container(const keyT &key) const {
        auto it = std::find(local_list_of_container_keys.list.begin(),
                            local_list_of_container_keys.list.end(), key);
        return it!=local_list_of_container_keys.list.end();
    }

    template<typename T>
    T allocator(World &world) const {
        if constexpr (is_world_constructible<T>::value) {
            return T(world);
        } else {
            return T();
        }
    }

    template<typename T>
    recordlistT store_other(madness::World &world, const T &source) {
        auto record = Recordlist<keyT>::compute_record(source);
        bool is_already_present= is_in_container(record);
        if (debug and world.rank()==0) {
            if (is_already_present) std::cout << "skipping ";
            if constexpr (Recordlist<keyT>::has_member_id<T>::value) {
                std::cout << "storing world object of " << type_name<T>::value() << "id " << source.id()
                << " to record " << record << std::endl;
            }
            std::cout << "storing object of " << type_name<T>::value() << " to record " << record << std::endl;
        }
        if constexpr (is_madness_function<T>::value) {
            if (source.is_compressed() and T::dimT>3) print("WARNING: storing compressed hi-dim `function");
        }

        // scope is important because of destruction ordering of world objects and fence
        if (is_already_present) {
            if (world.rank()==0) cache_stores++;
        } else {
            cloudtimer t(world,writing_time1);
            madness::archive::ContainerRecordOutputArchive ar(world, container, record);
            madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
            if (storage_policy==StoreFunctionPointer) {
                if constexpr (is_madness_function<T>::value) {
                    // store the pointer to the function, not the function itself
                    par & source.get_impl();
                    // store the pointer to the WorldObject in a list for later reference (replication/redistribution)
                    WorldObjectBase* wobj=source.get_impl().get();
                    world_object_base_list.push_back(wobj);
                } else {
                    // store everything else
                    par & source;
                }
            } else {
                // store everything else
                par & source;
            }
            local_list_of_container_keys+=record;
        }
        if (dofence) world.gop.fence();
        return recordlistT{record};
    }

public:
    /// load a vector from the cloud, pop records from recordlist
    ///
    /// @param[inout]    world	destination world
    /// @param[inout]    recordlist	list of records to load from (reduced by the first few elements)
    template<typename T>
    typename std::enable_if<is_vector<T>::value, T>::type
    do_load(World &world, recordlistT &recordlist) const {
        std::size_t sz = do_load<std::size_t>(world, recordlist);
        T target(sz);
        for (std::size_t i = 0; i < sz; ++i) {
            target[i] = do_load<typename T::value_type>(world, recordlist);
        }
        return target;
    }

    /// load a single object from the cloud, pop record from recordlist
    ///
    /// @param[inout]    world	destination world
    /// @param[inout]    recordlist	list of records to load from (reduced by the first element)
    template<typename T>
    typename std::enable_if<!is_vector<T>::value, T>::type
    do_load(World &world, recordlistT &recordlist) const {
        keyT record = recordlist.pop_front_and_return();
        if (force_load_from_cache) MADNESS_CHECK(is_cached(record));

        if (is_cached(record)) return load_from_cache<T>(world, record);
        if (debug) print("loading", type_name<T>::value(), "from container record", record, "to world", world.id());
        T target = allocator<T>(world);
        madness::archive::ContainerRecordInputArchive ar(world, container, record);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        if constexpr (is_madness_function<T>::value) {
            if (storage_policy==StoreFunctionPointer) {
                // load the pointer to the function, not the function itself
                // this is important for large functions, as they are not replicated
                // and only copied to subworlds when needed
                try {
                    typedef madness::FunctionImpl<typename T::typeT, T::dimT> implT;
                    std::shared_ptr<implT> impl;
                    par & impl;
                    target.set_impl(impl); // target now points to a universe function impl
                } catch (...) {
                    {
                        io_redirect_cout redirect;
                        print("failed to load function pointer from cloud, maybe the target is out of scope?");
                        print("record:", record, "world:", world.id());
                        print("function type:", type_name<T>::value());
                        print("\n");
                    }
                    MADNESS_EXCEPTION("load/store error of pointers in cloud", 1);
                }
            } else {
                // load everything else
                par & target;
            }
        } else {
            // load everything else
            par & target;
        }

        if (use_cache) {
            cache(world, target, record);
            if (is_replicated) container.erase(record);
        }

        return target;
    }

public:

    // overloaded
    template<typename T>
    recordlistT store_other(madness::World& world, const std::vector<T>& source) {
        if (debug and world.rank()==0)
            std::cout << "storing vector of " << type_name<T>::value() << " of size " << source.size() << std::endl;
        recordlistT l = store_other(world, source.size());
        for (const auto& s : source) l += store_other(world, s);
        if (dofence) world.gop.fence();
        if (debug and world.rank()==0) std::cout << "done with vector storing; container size "
                << container.size() << std::endl;
        return l;
    }

    /// store a tuple in multiple records
    template<typename... Ts>
    recordlistT store_tuple(World &world, const std::tuple<Ts...> &input) {
        recordlistT v;
        auto storeaway = [&](const auto &arg) {
            v += this->store(world, arg);
        };
        auto l = [&](Ts const &... arg) {
            ((storeaway(arg)), ...);
        };
        std::apply(l, input);
        return v;
    }

    /// load a tuple from the cloud, pop records from recordlist
    ///
    /// @param[inout]    world	destination world
    /// @param[inout]    recordlist	list of records to load from (reduced by the first few elements)
    template<typename T>
    T load_tuple(madness::World &world, recordlistT &recordlist) const {
        if (debug) std::cout << "loading tuple of type " << typeid(T).name() << " to world " << world.id() << std::endl;
        T target;
        std::apply([&](auto &&... args) {
            ((args = forward_load<typename std::remove_reference<decltype(args)>::type>(world, recordlist)), ...);
        }, target);
        return target;
    }
};

} /* namespace madness */

#endif /* SRC_MADNESS_WORLD_CLOUD_H_ */
