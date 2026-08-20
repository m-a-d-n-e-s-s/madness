
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
#include<algorithm>
#include<any>
#include<atomic>
#include<iomanip>
#include<limits>
#include<list>
#include<mutex>


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

/// Process map for the cloud's batch container

/// Routes a record to an explicitly assigned owner, falling back to a hash for any
/// record that was never registered. `set_owner` is collective, so all ranks agree on
/// routing without further communication. Unlike the default cloud container this map
/// is never reset to local, which is what keeps owner pinning stable.
template <typename keyT, typename hashfunT = Hash<keyT>>
class CloudOwnerPmap : public WorldDCPmapInterface<keyT> {
private:
    const int nproc;
    hashfunT hashfun;
    std::shared_ptr<std::map<keyT, ProcessID>> table;

public:
    CloudOwnerPmap(World& world, const hashfunT& hf = hashfunT())
        : nproc(world.mpi.nproc()), hashfun(hf), table(new std::map<keyT, ProcessID>()) {}

    /// collective: every rank must register the same (key, owner) pair
    void set_owner(const keyT& key, const ProcessID owner) { (*table)[key] = owner; }

    ProcessID owner(const keyT& key) const override {
        auto it = table->find(key);
        if (it != table->end()) return it->second;
        if (nproc == 1) return 0;
        return hashfun(key) % nproc;
    }

    void print_table(const std::string& tag = "") const {
        print("CloudOwnerPmap::table", tag, "size=", table->size(), "(nproc=", nproc, ")");
        for (const auto& kv : *table) {
            std::ostringstream os;
            os << "  key=0x" << std::hex << kv.first << std::dec << "  owner=" << kv.second;
            print(os.str());
        }
    }
};

/// BatchTransport holds a back-reference to its Cloud; its bodies are defined below the class
class Cloud;

using batch_keyT = madness::archive::ContainerRecordOutputArchive::keyT;
/// the serialized batch travels by value through a Future, so it must be
/// archive-serializable: a plain vector is, a shared_ptr is not
using batch_bytesT = std::vector<unsigned char>;

/// Point-to-point transfer of serialized function batches between universe ranks

/// The bytes stream straight from the owner's local batch store to the requester by
/// MPI point-to-point; they never ride inside an active-message payload, so there is
/// no eager-buffer limit and no extra copy on the wire.
///
/// Both endpoints are posted from **comm-thread AM handlers** (`WorldObject::send` runs
/// the member inline on the RMI receiver thread), never from worker tasks. That is what
/// buys overlap under worker saturation: at a tight protocol every worker sits in the
/// exchange kernel for a long time, so an endpoint posted as a *task* would queue behind
/// it and the MPI op would not be posted until compute ended. On the comm thread the RMI
/// loop's Testsome drives the rendezvous to completion during compute instead, leaving
/// only the final await for the worker.
///
/// Wire protocol:
///   1. requester (worker): record a pending slot keyed by `tag`, send `on_trigger`
///   2. owner `on_trigger` (comm thread): Isend the local bytes, reply `on_reply` with the size
///   3. requester `on_reply` (comm thread): size the buffer, post the Irecv, enqueue `finish_recv`
///   4. requester `finish_recv` (worker): await the Irecv and set the future
///
/// The size travels in the reply of step 2 rather than a separate Isend so that step 3 can
/// post the payload Irecv *during* compute; posting it at consume time would move the data
/// transfer to post-compute and lose the overlap.
class BatchTransport : public WorldObject<BatchTransport> {
public:
    /// tags live in the range MADNESS does not manage (safempi.h); 32767 is the
    /// conservative MPI_TAG_UB floor
    static constexpr int BATCH_TAG_BASE = 8192;
    static constexpr int BATCH_TAG_CAP = 32767;

    /// reply size meaning "the owner does not hold this record"; see on_trigger
    static constexpr std::size_t BATCH_NOT_FOUND = ~std::size_t(0);

    /// bytes per MPI message in a batch transfer; a larger payload is split into several

    /// MPI byte counts are `int`, and batches do exceed 2 GiB: a 161-orbital run at k=10
    /// measured 3.03 GiB, whose count narrowed to a negative int. MPI rejects that and
    /// SafeMPI throws it on the comm thread, where nothing catches it.
    static std::size_t batch_chunk_bytes() { return batch_chunk_bytes_; }

    /// Set the chunk size, for tests only.

    /// Collective and not thread-safe: call it on every rank before any transfer. It exists so
    /// a unit test can reach the multi-chunk path with an affordable payload.
    static void set_batch_chunk_bytes(const std::size_t n) {
        MADNESS_CHECK_THROW(n > 0 and n <= std::size_t(std::numeric_limits<int>::max()),
                            "batch chunk size must be positive and fit an MPI count");
        batch_chunk_bytes_ = n;
    }

    /// @param[in] universe  the world the cloud lives in (collective construction)
    /// @param[in] cloud     back-reference used to read owner-local batch bytes
    BatchTransport(World& universe, Cloud* cloud)
        : WorldObject<BatchTransport>(universe), cloud_(cloud), next_tag_(0) {
        this->process_pending();
    }

    /// Future to the serialized bytes of `record`, fetched from its owner

    /// Resolves locally without MPI when this rank owns the record. The trigger is in
    /// flight on return, so the round trip overlaps work until the future is consumed.
    Future<batch_bytesT> request(batch_keyT record);

private:
    Cloud* cloud_;                 ///< back-reference (not owned)
    std::atomic<int> next_tag_;

    /// requester-side state for one outstanding receive. Held by shared_ptr so the buffer
    /// address is stable for the Irecv and the entry can leave pending_ while finish_recv
    /// still owns its reference.
    struct PendingRecv {
        ProcessID owner = -1;
        int tag = -1;
        bool not_found = false;             ///< owner reported it does not hold the record
        batch_bytesT buf;                   ///< sized in on_reply
        std::vector<SafeMPI::Request> reqs; ///< one per chunk, posted in on_reply
        Future<batch_bytesT> fut;           ///< set by finish_recv
    };
    std::mutex pending_mtx_;
    std::map<int, std::shared_ptr<PendingRecv>> pending_;

    /// owner-side in-flight Isends, reaped lazily. Their buffers live in the cloud's
    /// batch container and stay valid for the duration, so an un-reaped Isend is harmless.
    std::mutex sends_mtx_;
    std::list<SafeMPI::Request> sends_;

    int alloc_tag() {
        const int span = BATCH_TAG_CAP - BATCH_TAG_BASE + 1;
        const int t = next_tag_.fetch_add(1);
        return BATCH_TAG_BASE + ((t % span) + span) % span;
    }

    static inline std::size_t batch_chunk_bytes_ = std::size_t(1) << 30;   ///< 1 GiB

    void reap_sends();

    /// owner side, comm thread: Isend the record's bytes, then reply with the count.
    /// Must not throw: an exception here escapes every task-level handler.
    void on_trigger(batch_keyT record, ProcessID requester, int tag);

    /// requester side, comm thread: size the buffer, post the Irecvs, enqueue finish_recv

    /// \param chunk  the chunking the *owner* used, echoed back so the receiver never has to
    ///               infer it. See on_trigger.
    void on_reply(int tag, std::size_t size, std::size_t chunk);

    /// requester side, worker task: await the background-progressed Irecv and set the future
    void finish_recv(int tag);
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
            case StoreFunction: os << "Function"; break;
            case StoreFunctionPointer: os << "FunctionPointer"; break;
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
    /// dedicated container for owner-pinned function batches; see store_batch / fetch_batch_p2p.
    /// Uses CloudOwnerPmap so each batch record lives on an explicitly chosen rank.
    std::shared_ptr<CloudOwnerPmap<keyT>> batch_pmap;
    mutable madness::WorldContainer<keyT, valueT> batch_container;
    /// constructed after batch_container so it is destroyed first, as WorldObject lifetimes require
    std::unique_ptr<BatchTransport> batch_transport_;
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
    Cloud(madness::World &universe) : container(universe),
        batch_pmap(new CloudOwnerPmap<keyT>(universe)),
        batch_container(universe, batch_pmap),
        batch_transport_(new BatchTransport(universe, this)),
        reading_time(0l), copy_time(0l), writing_time(0l),
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

    /// dump the record->owner table on rank 0; pair with the caller's own task-to-rank
    /// print to check that task assignment and batch routing agree
    void print_batch_owner_map(World& universe, const std::string& tag = "") const {
        if (universe.rank() != 0) return;
        batch_pmap->print_table(tag);
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
        double global_memsize=stats["memory_global_GB"].template get<double>();
        double max_record_size=stats["max_record_size"].template get<double>();
        double min_memsize=stats["memory_min_GB"].template get<double>();
        double max_memsize=stats["memory_max_GB"].template get<double>();
        double global_size=stats["container_size_global"].template get<double>();

        if (universe.rank()==0) {
            print("Cloud memory:");
            print("  replicated:",is_replicated);
            print("size of cloud (total)");
            print("  number of records:        ",global_size);
            print("  memory in GBytes:         ",global_memsize);
            print("size of cloud (average per node)");
            print("  number of records:        ",double(global_size)/universe.size());
            print("  memory in GBytes:         ",global_memsize/universe.size());
            print("min/max of node");
            print("  memory in GBytes:         ",min_memsize,max_memsize);
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
        // batch records are held separately, and for a caller that stores batches they are
        // the larger half; leaving them out would understate the cloud's footprint
        std::size_t batch_memsize=0;
        for (auto& item : batch_container) {
            batch_memsize+=item.second.size();
            max_record_size=std::max(max_record_size,item.second.size());
        }
        memsize+=batch_memsize;
        std::size_t global_memsize=memsize;
        std::size_t max_memsize=memsize;
        std::size_t min_memsize=memsize;
        double rss=madness::get_rss_usage_in_GB();
        double rss_av=rss;
        universe.gop.sum(global_memsize);
        universe.gop.max(max_memsize);
        universe.gop.max(max_record_size);
        universe.gop.min(min_memsize);
        universe.gop.max(rss);
        universe.gop.sum(rss_av);
        double byte2gbyte=1.0/(1024*1024*1024);

        // convert type(container item).second to GB, i.e. number of bytes in the container to GB
        double uchar2gbyte=byte2gbyte*sizeof(unsigned char);


        auto local_size=container.size();
        auto global_size=local_size;
        universe.gop.sum(global_size);
        auto batch_global_size=batch_container.size();
        universe.gop.sum(batch_global_size);
        std::size_t batch_global_memsize=batch_memsize;
        universe.gop.sum(batch_global_memsize);
        nlohmann::json j;
        j["container_size_global"] = global_size;
        j["batch_container_size_global"] = batch_global_size;
        j["batch_memory_global_GB"] = batch_global_memsize*uchar2gbyte;
        j["memory_global_GB"] = global_memsize*uchar2gbyte;
        j["memory_min_GB"] = min_memsize*uchar2gbyte;
        j["memory_max_GB"] = max_memsize*uchar2gbyte;
        j["memory_rss_GB_max"] = rss;
        j["memory_rss_GB_av"] = rss_av/universe.size();
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
        double global_memsize=stats["memory_global_GB"].template get<double>();
        double max_record_size=stats["max_record_size"].template get<double>();
        double min_memsize=stats["memory_min_GB"].template get<double>();
        double max_memsize=stats["memory_max_GB"].template get<double>();
        double global_size=stats["container_size_global"].template get<double>();

        print("Cloud memory:");
        print("  size of cloud (total)");
        print("    number of records:        ",global_size);
        print("    memory in GBytes:         ",global_memsize);
        // print("  size of cloud (average per node)");
        // print("    number of records:        ",double(global_size)/madness::world().size());
        // print("    memory in GBytes:         ",global_memsize*byte2gbyte/madness::world().size());
        print("  min/max of node");
        print("    memory in GBytes:         ",min_memsize,max_memsize);
        print("    max record size in GBytes:",max_record_size*byte2gbyte);
        // the owner-pinned batches, reported separately because they are a different lifetime and
        // usually the bulk of it. Zero unless some task stored batches.
        const double b_size = stats.value("batch_container_size_global", 0.0);
        if (b_size > 0.0) {
            print("  owner-pinned batches");
            print("    number of records:        ", b_size);
            print("    memory in GBytes:         ", stats.value("batch_memory_global_GB", 0.0));
        }
    }

    void clear_cache(World &subworld) {
        cached_objects.clear();
        local_list_of_container_keys.list.clear();
        subworld.gop.fence();
    }

    void clear() {
        container.clear();
        // The owner-pinned batches too. They are not leaked without this -- an SCF run derives the
        // same record keys each time, so the next application overwrites them -- but they are the
        // largest thing the cloud holds, and without this they stay resident through everything
        // that follows the exchange in an iteration, which is where the memory ceiling actually is.
        batch_container.clear();
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
                    // print("replicating function per node",fimpl);;
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
            // world.gop.fence();
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

    /// Register the owner of a batch record; local map insert, no communication

    /// Collective in the same sense as store_batch: every rank must call it with an
    /// identical (record, owner) pair or fetches will route inconsistently. Separating
    /// registration from the payload lets all ranks replicate the routing while each
    /// owner stores only its own bytes, over a size-1 subworld.
    void register_batch_owner(const keyT record, const ProcessID owner) {
        batch_pmap->set_owner(record, owner);
    }

    /// Store a batch of functions as one owner-pinned record

    /// The whole vector -- its size and each function -- is serialized into a single
    /// record in the batch container and routed to `owner`, so one batch is one record
    /// with one owner. Must be called with identical (owner, record) on every rank of
    /// `world`.
    ///
    /// @param[in] fence  false lets a caller storing many batches emit one fence for all
    ///                   of them; the collective gather inside the archive self-synchronizes
    template<typename T, std::size_t NDIM>
    keyT store_batch(madness::World& world, const std::vector<Function<T, NDIM>>& batch,
                     const ProcessID owner, const keyT record, const bool fence = true) {
        if (is_replicated) {
            print("Cloud contents are replicated and read-only!");
            MADNESS_EXCEPTION("cloud error", 1);
        }
        batch_pmap->set_owner(record, owner);
        cloudtimer t_batch(world, batch_store_time);
        {
            madness::archive::ContainerRecordOutputArchive ar(world, batch_container, record);
            madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
            par.set_dofence(false);
            std::size_t fsize = batch.size();
            par & fsize;
            for (std::size_t i = 0; i < fsize; ++i) par & batch[i];
        }
        if (dofence and fence) world.gop.fence();
        return record;
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
                std::size_t keylistsize = 0;
                world.mpi.Bcast(&keylistsize,sizeof(keylistsize),MPI_BYTE,rank);
                for (size_t i=0; i<keylistsize; i++) {
                    keyT key;
                    world.mpi.Bcast(&key,sizeof(key),MPI_BYTE,rank);
                    std::size_t sz = 0;
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
    mutable std::atomic<long> batch_store_time=0l;       ///< store_batch wall time, microseconds
    mutable std::atomic<long> batch_find_time=0l;        ///< waiting on the p2p transfer, microseconds
    mutable std::atomic<long> batch_deserialize_time=0l; ///< deserializing the bytes, microseconds
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

public:

    /// the owner of a batch record; a pmap lookup, no communication
    ProcessID batch_owner(const keyT record) const {
        return batch_pmap->owner(record);
    }

private:

    /// only the transport reads owner-local bytes; callers go through fetch_batch_p2p
    friend class BatchTransport;

    /// bytes of a batch record held by this rank, empty if it holds none

    /// Never throws: the callers are BatchTransport's comm-thread handlers, where an
    /// escaping exception bypasses every task-level handler and surfaces as an
    /// unattributable abort. A miss is reported to the requester instead, and raised
    /// there in task context.
    ///
    /// Emptiness is a sound "not found" marker because a stored batch always begins with
    /// its serialized element count, so a present record is never zero bytes.
    valueT try_get_local_batch_bytes(const keyT record) const {
        typename madness::WorldContainer<keyT, valueT>::const_accessor acc;
        if (batch_container.find(acc, record)) return acc->second;
        return valueT();
    }

    /// stable pointer and size of a local batch record, {nullptr,0} if this rank holds none

    /// Lets the comm thread Isend without copying the payload. The accessor lock is
    /// released on return, but the address stays valid because batch records are neither
    /// erased nor mutated between the store and the end of the consuming operation.
    /// Never throws, for the reason given on try_get_local_batch_bytes.
    std::pair<const unsigned char*, std::size_t> try_get_local_batch_ptr(const keyT record) const {
        typename madness::WorldContainer<keyT, valueT>::const_accessor acc;
        if (batch_container.find(acc, record))
            return {acc->second.data(), acc->second.size()};
        return {nullptr, 0};
    }

public:

    /// start fetching `record` from its owner; the trigger is in flight on return
    Future<batch_bytesT> request_batch_bytes_async(const keyT record) const {
        return batch_transport_->request(record);
    }

    /// turn the bytes of a p2p transfer into the batch of functions

    /// Blocks on `fut` only if the transfer has not landed yet. Runs in a task, which is
    /// where a missing record is reported so the failure is attributable.
    /// @param[in] cache_result  default false: the cloud-side cache is not safe to keep
    ///                          across changes of the calling world, so opting in is the
    ///                          caller's decision
    template<typename T, std::size_t NDIM>
    std::vector<Function<T, NDIM>> deserialize_batch_p2p(madness::World& subworld,
            Future<batch_bytesT> fut, const keyT record,
            const bool cache_result = false) const {
        typedef std::vector<Function<T, NDIM>> vecfuncT;
        if (is_cached(record)) return load_from_cache<vecfuncT>(subworld, record);
        cloudtimer t(subworld, reading_time);
        const double w0 = wall_time();
        batch_bytesT bytes = fut.get();
        const double w1 = wall_time();
        batch_find_time += long((w1 - w0) * 1.e6);
        MADNESS_CHECK_THROW(not bytes.empty(),
                            "deserialize_batch_p2p: the owner does not hold this batch record");
        vecfuncT batch;
        {
            madness::archive::VectorInputArchive var(bytes);
            madness::archive::ParallelInputArchive<madness::archive::VectorInputArchive> par(subworld, var);
            std::size_t fsize = 0;
            par & fsize;
            batch.resize(fsize);
            for (std::size_t i = 0; i < fsize; ++i) par & batch[i];
        }
        batch_deserialize_time += long((wall_time() - w1) * 1.e6);
        if (use_cache and cache_result) cache(subworld, batch, record);
        return batch;
    }

    /// fetch a batch stored by store_batch; resolves without MPI when this rank owns it
    template<typename T, std::size_t NDIM>
    std::vector<Function<T, NDIM>> fetch_batch_p2p(madness::World& subworld,
            const keyT record, const bool cache_result = false) const {
        typedef std::vector<Function<T, NDIM>> vecfuncT;
        if (is_cached(record)) return load_from_cache<vecfuncT>(subworld, record);
        return deserialize_batch_p2p<T, NDIM>(subworld, request_batch_bytes_async(record),
                                              record, cache_result);
    }

private:

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
        auto print_debug = [&]() {
             if (is_already_present) std::cout << "skipping ";
             if constexpr (Recordlist<keyT>::has_member_id<T>::value) {
                 std::cout << "storing world object of " << type_name<T>::value() << "id " << source.id()
                 << " to record " << record << std::endl;
             }
             std::cout << "storing object of " << type_name<T>::value() << " to record " << record << std::endl;
        };
        if (debug and world.rank()==0) print_debug();

        try {
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
        } catch (std::exception& e) {
            {
                print("exception during storing object in cloud: ", e.what());
                print("failed to store object in cloud, maybe the object is not serializable or out of scope?");
                print("record:", record, "world:", world.id());
                print("object type:", type_name<T>::value());
                print("\n");
                io_redirect_cout redirect;
                print("exception during storing object in cloud: ", e.what());
                print("failed to store object in cloud, maybe the object is not serializable or out of scope?");
                print("record:", record, "world:", world.id());
                print("object type:", type_name<T>::value());
                print("\n");
                std::rethrow_exception(std::current_exception());
            }
            MADNESS_EXCEPTION("load/store error in cloud", 1);
        }
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

// ---- BatchTransport bodies; they need the complete Cloud ----

inline Future<batch_bytesT> BatchTransport::request(batch_keyT record) {
    World& u = this->get_world();
    const ProcessID owner = cloud_->batch_owner(record);
    if (owner == u.rank()) {
        // local: no MPI. An absent record yields empty bytes, which the consumer
        // reports in task context, exactly as for the remote path.
        return Future<batch_bytesT>(cloud_->try_get_local_batch_bytes(record));
    }
    const int tag = alloc_tag();
    auto p = std::make_shared<PendingRecv>();
    p->owner = owner;
    p->tag = tag;
    {
        std::lock_guard<std::mutex> g(pending_mtx_);
        pending_[tag] = p;
    }
    // send, not add: the trigger runs inline on the owner's comm thread, so the owner
    // posts its Isend without queueing behind its own saturated workers
    this->send(owner, &BatchTransport::on_trigger, record, u.rank(), tag);
    return p->fut;
}

inline void BatchTransport::reap_sends() {
    std::lock_guard<std::mutex> g(sends_mtx_);
    for (auto it = sends_.begin(); it != sends_.end(); ) {
        if (it->Test()) it = sends_.erase(it);
        else ++it;
    }
}

inline void BatchTransport::on_trigger(batch_keyT record, ProcessID requester, int tag) {
    World& u = this->get_world();
    reap_sends();
    auto ptr_size = cloud_->try_get_local_batch_ptr(record);
    if (ptr_size.first == nullptr) {
        // Not held here -- reply with the sentinel and post nothing. Throwing instead
        // would escape this comm-thread handler past every task-level handler and abort
        // the run with no attribution; the requester raises it in task context.
        this->send(requester, &BatchTransport::on_reply, tag, BATCH_NOT_FOUND, std::size_t(0));
        return;
    }
    const std::size_t n = ptr_size.second;
    // Sent along rather than read again by the requester: the size is a process-local static, so
    // the two ends can disagree, and a receiver that guessed would post Irecvs that do not match
    // the messages -- MPI_ERR_TRUNCATE on the comm thread, where nothing catches it.
    const std::size_t chunk = batch_chunk_bytes_;
    // MPI does not overtake between messages sharing source, tag and communicator, so posting
    // the chunks in ascending offset order matches the requester's Irecvs pairwise without a
    // per-chunk tag. That does rely on one transfer per tag, which alloc_tag's span makes true.
    {
        std::lock_guard<std::mutex> g(sends_mtx_);
        for (std::size_t off = 0; off < n; off += chunk) {
            const std::size_t len = std::min(chunk, n - off);
            sends_.push_back(u.mpi.Isend(ptr_size.first + off, int(len), MPI_BYTE, requester, tag));
        }
    }
    // the size rides in the reply so the requester can post its Irecvs now, during its
    // own compute, and let the rendezvous finish in the background
    this->send(requester, &BatchTransport::on_reply, tag, n, chunk);
}

inline void BatchTransport::on_reply(int tag, std::size_t size, std::size_t chunk) {
    World& u = this->get_world();
    std::shared_ptr<PendingRecv> p;
    {
        std::lock_guard<std::mutex> g(pending_mtx_);
        auto it = pending_.find(tag);
        MADNESS_CHECK(it != pending_.end());
        p = it->second;
    }
    if (size == BATCH_NOT_FOUND) {
        // no Isend was posted, so post no Irecv; empty bytes mark the miss
        {
            std::lock_guard<std::mutex> g(pending_mtx_);
            pending_.erase(tag);
        }
        p->not_found = true;
        p->fut.set(batch_bytesT());
        return;
    }
    p->buf.resize(size);
    // the owner's framing, not ours; see on_trigger. Positive whenever the record was found, and
    // unguarded because a throw on this thread is what the framing exists to avoid
    for (std::size_t off = 0; off < size; off += chunk) {
        const std::size_t len = std::min(chunk, size - off);
        p->reqs.push_back(u.mpi.Irecv(p->buf.data() + off, int(len), MPI_BYTE, p->owner, tag));
    }
    u.taskq.add(this, &BatchTransport::finish_recv, tag);
}

inline void BatchTransport::finish_recv(int tag) {
    std::shared_ptr<PendingRecv> p;
    {
        std::lock_guard<std::mutex> g(pending_mtx_);
        auto it = pending_.find(tag);
        MADNESS_CHECK(it != pending_.end());
        p = it->second;
        pending_.erase(it);
    }
    // the comm thread has been progressing these Irecvs all along, so the awaits are short
    for (auto& r : p->reqs) World::await(r, true);
    p->fut.set(std::move(p->buf));
}

} /* namespace madness */

#endif /* SRC_MADNESS_WORLD_CLOUD_H_ */
