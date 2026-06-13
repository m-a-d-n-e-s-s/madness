
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
#include<atomic>
#include<iomanip>
#include<memory>
#include<map>
#include<mutex>
#include<list>
#include<utility>


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

/// Process map for the cloud's dedicated batch container.

/// Routes a record to an explicitly assigned owner rank, falling back to a
/// hash for any record that was never registered. The owner table is
/// populated collectively by Cloud::store_batch (every rank registers the
/// same record->owner pair), so all ranks agree on routing without extra
/// communication. Unlike the default cloud container, this map is never
/// reset to local (no replication), which keeps owner pinning stable.
template <typename keyT, typename hashfunT = Hash<keyT>>
class CloudOwnerPmap : public WorldDCPmapInterface<keyT> {
private:
    const int nproc;
    hashfunT hashfun;
    std::shared_ptr<std::map<keyT, ProcessID>> table;

public:
    CloudOwnerPmap(World& world, const hashfunT& hf = hashfunT())
        : nproc(world.mpi.nproc()), hashfun(hf), table(new std::map<keyT, ProcessID>()) {}

    /// Register the owner of a record. Collective: must be called with
    /// identical (key, owner) on every rank so all ranks route consistently.
    void set_owner(const keyT& key, const ProcessID owner) { (*table)[key] = owner; }

    ProcessID owner(const keyT& key) const override {
        auto it = table->find(key);
        if (it != table->end()) return it->second;
        if (nproc == 1) return 0;
        return hashfun(key) % nproc;
    }

    /// dump the (key -> owner) registrations to stdout in a stable order.
    /// Intended for diagnostic printing on universe rank 0 at high print levels.
    void print_table(const std::string& tag = "") const {
        print("CloudOwnerPmap::table", tag, "size=", table->size(), "(nproc=", nproc, ")");
        for (const auto& kv : *table) {
            std::ostringstream os;
            os << "  key=0x" << std::hex << kv.first << std::dec << "  owner=" << kv.second;
            print(os.str());
        }
    }
};

// forward declaration; BatchTransport holds a back-reference to its Cloud and
// its method bodies (which touch Cloud members) are defined after the Cloud class.
class Cloud;

/// types used by the owner-pinned batch point-to-point transfer (see BatchTransport)
using batch_keyT     = madness::archive::ContainerRecordOutputArchive::keyT;   // == long
/// the serialized batch travels by value through a Future; a plain
/// std::vector<unsigned char> is archive-serializable (shared_ptr is not, and
/// Future<T> instantiates the serialization path unconditionally).
using batch_bytesT   = std::vector<unsigned char>;

/// Point-to-point transfer of serialized function batches between universe ranks.
///
/// Replaces the WorldContainer::find round-trip used by the row-owner exchange
/// fetch path. The raw serialized bytes stream straight from the owner's local
/// batch store to the requester via MPI point-to-point; the blob never rides
/// inside an active-message payload, so there is no eager-buffer limit and no
/// extra serialize/deserialize copy on the wire (pure memcpy over wire).
///
/// Both transfer endpoints are posted from **comm-thread AM handlers** (via
/// WorldObject::send, whose handler runs the member inline on the RMI receiver
/// thread — see world_object.h), never from worker tasks. This is the key to
/// overlap under worker saturation: at high protocol every worker runs the
/// exchange kernel for ~10^3 s, so a transfer endpoint posted as a *task* (the
/// previous do_send/receive_task design) would queue behind the kernel and the
/// MPI op would not even be posted until compute ends — zero overlap. Posting on
/// the comm thread sidesteps the worker queue entirely, and the RMI loop's
/// continuous Testsome (worldrmi.cc) then drives the rendezvous to completion in
/// the background *during* compute. The only worker-side step is the final await
/// at consume, which returns immediately when the transfer already landed.
///
/// Wire protocol (no MPI op ever touches a worker task):
///   1. requester (worker, pre-compute): record pending slot keyed by `tag`,
///      send(owner, &on_trigger, record, requester, tag).  Returns a Future.
///   2. owner on_trigger (comm thread): look up the local bytes (stable pointer,
///      no copy), Isend(bytes) on `tag`, send(requester, &on_reply, tag, size).
///   3. requester on_reply (comm thread): allocate the receive buffer, post
///      Irecv(bytes) on `tag`, enqueue finish_recv.
///   4. requester finish_recv (worker, post-compute): await the Irecv (already
///      progressed in the background) and set the Future.
///
/// The size travels in the reply AM (step 2) rather than a separate Isend so the
/// requester can post the bytes Irecv *during* compute (step 3) — posting it only
/// at consume would defer the rendezvous data movement to post-compute and lose
/// the overlap on the large payload.
class BatchTransport : public WorldObject<BatchTransport> {
public:
    /// tags live in [8192, MPI_TAG_UB], the range MADNESS does not manage
    /// (see safempi.h); 32767 is the conservative MPI_TAG_UB floor.
    static constexpr int BATCH_TAG_BASE = 8192;
    static constexpr int BATCH_TAG_CAP  = 32767;

    /// @param[in] universe  the world the cloud lives in (collective construction)
    /// @param[in] cloud     back-reference used to read owner-local batch bytes
    BatchTransport(World& universe, Cloud* cloud)
        : WorldObject<BatchTransport>(universe), cloud_(cloud), next_tag_(0) {
        this->process_pending();
    }

    /// Requester entry point: returns a future to the serialized bytes of
    /// `record`, fetched from its owner. Resolves locally (no MPI) when this
    /// rank owns the record. The trigger message is in flight on return, so the
    /// round-trip overlaps subsequent work until the future is consumed.
    Future<batch_bytesT> request(batch_keyT record);

private:
    Cloud* cloud_;                 ///< back-reference (not owned)
    std::atomic<int> next_tag_;    ///< per-rank round-robin tag counter

    /// requester-side state for one outstanding receive, created in request()
    /// and completed in finish_recv. Held by shared_ptr so the buffer address is
    /// stable for the Irecv and the entry can be erased from pending_ while
    /// finish_recv still owns its reference.
    struct PendingRecv {
        ProcessID owner = -1;
        int tag = -1;
        batch_bytesT buf;              ///< receive buffer (sized in on_reply)
        SafeMPI::Request req;          ///< the Irecv request (posted in on_reply)
        Future<batch_bytesT> fut;      ///< set by finish_recv at completion
    };
    std::mutex pending_mtx_;
    std::map<int, std::shared_ptr<PendingRecv>> pending_;   ///< keyed by tag

    /// owner-side in-flight Isend requests, reaped lazily on the comm thread.
    /// The send buffers live in the cloud's batch_container (stable for the
    /// duration of the exchange), so an un-reaped Isend is harmless.
    std::mutex sends_mtx_;
    std::list<SafeMPI::Request> sends_;

    /// allocate the next point-to-point tag in [BATCH_TAG_BASE, BATCH_TAG_CAP]
    int alloc_tag() {
        const int span = BATCH_TAG_CAP - BATCH_TAG_BASE + 1;
        const int t = next_tag_.fetch_add(1);
        return BATCH_TAG_BASE + ((t % span) + span) % span;
    }

    /// drop completed Isend requests from sends_ (comm thread, amortized)
    void reap_sends();

    /// owner side, comm-thread handler: Isend the record's local bytes to
    /// `requester` on `tag`, then reply with the byte count. Non-blocking.
    void on_trigger(batch_keyT record, ProcessID requester, int tag);

    /// requester side, comm-thread handler: allocate the buffer for `size` bytes,
    /// post the Irecv on `tag`, and enqueue finish_recv. Non-blocking.
    void on_reply(int tag, std::size_t size);

    /// requester side, worker task: await the (background-progressed) Irecv and
    /// set the future. Runs post-compute, so its await returns promptly.
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
    bool batch_debug = false; ///< per-call diagnostics on store_batch / fetch_batch / deserialize_batch (narrower than debug)
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
    /// dedicated container for owner-pinned function batches; see store_batch / fetch_batch.
    /// uses CloudOwnerPmap so each batch record lives on an explicitly chosen rank.
    std::shared_ptr<CloudOwnerPmap<keyT>> batch_pmap;
    mutable madness::WorldContainer<keyT, valueT> batch_container;
    /// point-to-point fetch of owner-pinned batches (see BatchTransport).
    /// Constructed after batch_container so it is destroyed before it (reverse
    /// of construction, as WorldObject lifetimes require).
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

    /// per-call diagnostics on store_batch / fetch_batch / deserialize_batch /
    /// fetch_batch_record_async (BATCH_REQ / BATCH_DESER lines). Narrower than
    /// set_debug, so the existing chatty store-prints don't fire too.
    void set_batch_debug(bool value) {
        batch_debug = value;
    }

    /// counters exposed so external code (the exchange task) can tally
    /// algorithm-specific prefetch outcomes into the cloud diagnostics bucket.
    void tally_batch_prefetch_hit()  const { ++batch_prefetch_hits;  }
    void tally_batch_prefetch_miss() const { ++batch_prefetch_misses; }

    /// dump the CloudOwnerPmap's record->owner table (rank 0 only).
    /// Diagnostic — pair with the task's owner_map_ print to verify the
    /// algorithm's task-to-rank assignment aligns with the cloud's batch routing.
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

        // batch_container (StoreFunctionBatched payloads). Zero unless any
        // task called store_batch this invocation.
        std::size_t b_memsize = 0;
        std::size_t b_max_record_size = 0;
        for (auto& item : batch_container) {
            b_memsize += item.second.size();
            b_max_record_size = std::max(b_max_record_size, item.second.size());
        }
        std::size_t b_global_memsize = b_memsize;
        std::size_t b_max_memsize = b_memsize;
        std::size_t b_min_memsize = b_memsize;
        universe.gop.sum(b_global_memsize);
        universe.gop.max(b_max_memsize);
        universe.gop.max(b_max_record_size);
        universe.gop.min(b_min_memsize);
        auto b_local_size = batch_container.size();
        auto b_global_size = b_local_size;
        universe.gop.sum(b_global_size);

        nlohmann::json j;
        j["container_size_global"] = global_size;
        j["memory_global_GB"] = global_memsize*uchar2gbyte;
        j["memory_min_GB"] = min_memsize*uchar2gbyte;
        j["memory_max_GB"] = max_memsize*uchar2gbyte;
        j["memory_rss_GB_max"] = rss;
        j["memory_rss_GB_av"] = rss_av/universe.size();
        j["max_record_size"] = max_record_size;
        j["batch_container_size_global"] = b_global_size;
        j["batch_memory_global_GB"] = b_global_memsize*uchar2gbyte;
        j["batch_memory_min_GB"]    = b_min_memsize*uchar2gbyte;
        j["batch_memory_max_GB"]    = b_max_memsize*uchar2gbyte;
        j["batch_max_record_size"]  = b_max_record_size;
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

        // writing_time1 is the per-record inner write (now also covers store_batch);
        // surface it so the cost of serializing each record is visible.
        double wtime1_max = double(writing_time1)*1.e-6;
        universe.gop.max(wtime1_max);

        // StoreFunctionBatched-specific timers (zero for non-batched runs).
        double bstime_max = double(batch_store_time)*1.e-6;
        double bftime_max = double(batch_find_time)*1.e-6;
        double bftime_acc = double(batch_find_time)*1.e-6;
        double bdtime_max = double(batch_deserialize_time)*1.e-6;
        double bdtime_acc = double(batch_deserialize_time)*1.e-6;
        universe.gop.max(bstime_max);
        universe.gop.max(bftime_max);
        universe.gop.sum(bftime_acc);
        universe.gop.max(bdtime_max);
        universe.gop.sum(bdtime_acc);
        long bphits  = long(batch_prefetch_hits);
        long bpmiss  = long(batch_prefetch_misses);
        universe.gop.sum(bphits);
        universe.gop.sum(bpmiss);

        nlohmann::json j;
        j["reading_time_max_s"] = rtime_max;
        j["reading_time_acc_s"] = rtime_acc;
        j["reading_time_av_s"] = rtime_av;
        j["copy_time_max_s"] = ctime_max;
        j["copy_time_acc_s"] = ctime_acc;
        j["copy_time_av_s"] = ctime_av;
        j["writing_time_s"] = wtime;
        j["writing_time1_max_s"] = wtime1_max;
        j["replication_time_s"] = ptime;
        j["target_replication_time_s"] = tptime;
        j["cache_reads"] = creads;
        j["cache_stores"] = cstores;
        j["batch_store_time_max_s"]       = bstime_max;
        j["batch_find_time_max_s"]        = bftime_max;
        j["batch_find_time_acc_s"]        = bftime_acc;
        j["batch_deserialize_time_max_s"] = bdtime_max;
        j["batch_deserialize_time_acc_s"] = bdtime_acc;
        j["batch_prefetch_hits"]   = bphits;
        j["batch_prefetch_misses"] = bpmiss;
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
        print("cloud per-record storing wall time max         ", timings.value("writing_time1_max_s", 0.0));
        print("cloud replication wall time                    ", ptime);
        print("target replication wall time                   ", tptime);
        print("cloud max reading time (all procs)             ", rtime_max, std::defaultfloat);
        print("cloud average reading cpu time (all procs)     ", rtime_av, std::defaultfloat);
        print("cloud accumulated reading cpu time (all procs) ", rtime_acc, std::defaultfloat);
        std::cout << std::setprecision(precision) << std::scientific;
        print("cloud cache stores                             ", long(cstores));
        print("cloud cache loads                              ", long(creads));

        // StoreFunctionBatched-specific block. Zero for non-batched runs;
        // prints to keep both backends comparable.
        std::cout << std::fixed << std::setprecision(3);
        print("batch store wall time max (all procs)          ", timings.value("batch_store_time_max_s", 0.0));
        print("batch find wall time max (all procs)           ", timings.value("batch_find_time_max_s", 0.0));
        print("batch find wall time acc (all procs)           ", timings.value("batch_find_time_acc_s", 0.0));
        print("batch deserialize wall time max (all procs)    ", timings.value("batch_deserialize_time_max_s", 0.0));
        print("batch deserialize wall time acc (all procs)    ", timings.value("batch_deserialize_time_acc_s", 0.0));
        std::cout << std::setprecision(precision) << std::scientific;
        print("batch prefetch hits                            ", long(timings.value("batch_prefetch_hits", 0l)));
        print("batch prefetch misses                          ", long(timings.value("batch_prefetch_misses", 0l)));
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
        // batch_container (StoreFunctionBatched payloads). Zero unless any
        // task stored owner-pinned batches this invocation.
        const double b_global_size   = stats.value("batch_container_size_global", 0.0);
        if (b_global_size > 0.0) {
            const double b_global_memsize    = stats.value("batch_memory_global_GB", 0.0);
            const double b_min_memsize       = stats.value("batch_memory_min_GB", 0.0);
            const double b_max_memsize       = stats.value("batch_memory_max_GB", 0.0);
            const double b_max_record_size   = stats.value("batch_max_record_size", 0.0);
            print("Cloud batch_container memory (StoreFunctionBatched):");
            print("  size of batch container (total)");
            print("    number of records:        ", b_global_size);
            print("    memory in GBytes:         ", b_global_memsize);
            print("  min/max of node");
            print("    memory in GBytes:         ", b_min_memsize, b_max_memsize);
            print("    max record size in GBytes:", b_max_record_size*byte2gbyte);
        }
    }

    void clear_cache(World &subworld) {
        cached_objects.clear();
        local_list_of_container_keys.list.clear();
        subworld.gop.fence();
    }

    void clear() {
        container.clear();
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
        batch_store_time=0l;
        batch_find_time=0l;
        batch_deserialize_time=0l;
        batch_prefetch_hits=0l;
        batch_prefetch_misses=0l;
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

    /// Register the owner of a batch record WITHOUT storing any payload.

    /// Local map insert, no communication. Like store_batch's internal routing,
    /// it must be called with identical (record, owner) on EVERY rank so all
    /// ranks route fetches consistently. Used by the distributed per-owner store
    /// path: all ranks register routing for all records, then each owner calls
    /// store_batch only for its own records (over its size-1 subworld). Separating
    /// registration from payload keeps routing replicated while the byte-store is
    /// owner-local.
    void register_batch_owner(const keyT record, const ProcessID owner) {
        batch_pmap->set_owner(record, owner);
    }

    /// Store a batch of functions as a single owner-pinned record.

    /// The whole vector (size + each function) is serialized into one record
    /// in the dedicated batch container, routed to `owner` via CloudOwnerPmap.
    /// Mirrors the save_function / load_function idiom (vmra.h) so a batch is
    /// exactly one cloud record: one record -> one owner.
    ///
    /// Collective and must be called with identical (owner, record) on every
    /// rank of `world` (typically the universe), like the existing store():
    /// every rank registers the same record->owner so all ranks route finds
    /// consistently, while only rank 0 writes the bytes.
    ///
    /// @param[in] world   the world the batch currently lives in (universe)
    /// @param[in] batch   the functions to store
    /// @param[in] owner   the rank that should physically hold the record
    /// @param[in] record  the (caller-chosen, globally consistent) record key
    /// @param[in] fence   if true (default), fence on `world` before returning.
    ///                    Pass false when emitting many batches back-to-back so
    ///                    one outer fence at the end of the batching loop
    ///                    suffices; cuts setup from O(nsubworld) fences to one.
    /// @return the record key (for manifest bookkeeping)
    template<typename T, std::size_t NDIM>
    keyT store_batch(madness::World& world, const std::vector<Function<T, NDIM>>& batch,
                     const ProcessID owner, const keyT record, const bool fence = true) {
        if (is_replicated) {
            print("Cloud contents are replicated and read-only!");
            MADNESS_EXCEPTION("cloud error", 1);
        }
        // collective: every rank registers the routing so owner(record) agrees everywhere
        batch_pmap->set_owner(record, owner);
        const double wall0 = wall_time();
        cloudtimer t_batch(world, batch_store_time);   // dedicated bucket for batch payload writes (separate from writing_time1)
        {
            madness::archive::ContainerRecordOutputArchive ar(world, batch_container, record);
            madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
            // Phase 4a-i: skip the per-function internal gop.fences inside the gather.
            // The inputs are complete (read-only) and the collective MPI_Gatherv
            // self-synchronizes; store_batches() emits a single fence after all
            // batches. Removes the per-function fence0/fence1 (~45% of store cost).
            par.set_dofence(false);
            std::size_t fsize = batch.size();
            par & fsize;
            for (std::size_t i = 0; i < fsize; ++i) par & batch[i];
        }
        if (dofence and fence) world.gop.fence();
        if (batch_debug and world.rank() == 0) {
            std::ostringstream os;
            os << "BATCH_STORE rank=" << world.rank()
               << " t=" << std::fixed << std::setprecision(3) << wall_time()
               << " key=0x" << std::hex << record << std::dec
               << " owner=" << owner
               << " nfuncs=" << batch.size()
               << " wall=" << std::fixed << std::setprecision(3) << (wall_time() - wall0);
            print(os.str());
        }
        return record;
    }

    /// Load a batch of functions stored by store_batch (synchronous).

    /// Blocks on the find round-trip to the record's owner; resolves locally
    /// when the calling rank owns the record (the held-vf case). Caches the
    /// deserialized batch by record (per-world), like do_load. The batch is
    /// reconstructed into `subworld`.
    ///
    /// NOTE (Phase A): synchronous. The overlapped (async / Future) variant
    /// that hides the find latency behind the previous task's compute is
    /// added together with its caller (the exchange prefetch hooks).
    ///
    /// @param[in] subworld      destination world (size 1 for the owner algorithm)
    /// @param[in] record        the record key returned by store_batch
    /// @param[in] cache_result  if true (default), cache the deserialized batch
    ///                          for subsequent same-record lookups (dedup). Pass
    ///                          false for transient batches that are used once
    ///                          and released, to avoid unbounded cache growth.
    template<typename T, std::size_t NDIM>
    std::vector<Function<T, NDIM>> fetch_batch(madness::World& subworld, const keyT record,
                                               const bool cache_result = true) const {
        typedef std::vector<Function<T, NDIM>> vecfuncT;
        if (is_cached(record)) return load_from_cache<vecfuncT>(subworld, record);
        const double wall0 = wall_time();
        cloudtimer t(subworld, reading_time);
        const ProcessID owner = batch_pmap->owner(record);
        const bool local = (owner == subworld.rank());  // strictly: subworld is size 1 so this is "owner == this rank"
        vecfuncT batch;
        double find_dt = 0.0, deser_dt = 0.0;
        {
            if (batch_debug) {
                fprintf(stderr, "FIXB_DBG cloud.fetch_batch BEFORE_FIND key=0x%lx owner=%d local=%d sw_rank=%d sw_id=%lu t=%.6f\n",
                        static_cast<unsigned long>(record), owner, int(local),
                        subworld.rank(), static_cast<unsigned long>(subworld.id()), wall_time());
                fflush(stderr);
            }
            // The find round-trip lives inside ContainerRecordInputArchive's ctor
            // (it does dc.find(record).get()). Time it separately from the
            // byte-fed deserialize loop so we can see the cost of the AM
            // independently from the local CPU work.
            const double w0 = wall_time();
            madness::archive::ContainerRecordInputArchive ar(subworld, batch_container, record);
            const double w1 = wall_time();
            if (batch_debug) {
                fprintf(stderr, "FIXB_DBG cloud.fetch_batch AFTER_FIND  key=0x%lx find_ms=%.3f t=%.6f\n",
                        static_cast<unsigned long>(record), (w1-w0)*1e3, wall_time());
                fflush(stderr);
            }
            find_dt = w1 - w0;
            batch_find_time += long(find_dt * 1e6);
            madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(subworld, ar);
            std::size_t fsize = 0;
            par & fsize;
            batch.resize(fsize);
            for (std::size_t i = 0; i < fsize; ++i) par & batch[i];
            deser_dt = wall_time() - w1;
            batch_deserialize_time += long(deser_dt * 1e6);
        }
        if (use_cache and cache_result) cache(subworld, batch, record);
        if (batch_debug) {
            std::ostringstream os;
            os << "BATCH_REQ rank=" << subworld.rank()
               << " sw=" << subworld.id()
               << " t=" << std::fixed << std::setprecision(3) << wall_time()
               << " key=0x" << std::hex << record << std::dec
               << " owner=" << owner
               << " kind=" << (local ? "sync_local" : "sync_remote")
               << " nfuncs=" << batch.size()
               << " find_ms=" << std::fixed << std::setprecision(3) << (find_dt*1e3)
               << " deser_ms=" << std::fixed << std::setprecision(3) << (deser_dt*1e3)
               << " wall_ms=" << std::fixed << std::setprecision(3) << ((wall_time()-wall0)*1e3);
            print(os.str());
        }
        return batch;
    }

    /// iterator type of the batch container (for prefetch futures)
    typedef typename madness::WorldContainer<keyT, valueT>::iterator batch_iterator;

    /// serialized-bytes payload type for the p2p batch fetch (see BatchTransport).
    /// Re-exported so callers can spell prefetch futures as Future<Cloud::batch_bytesT>.
    using batch_bytesT = madness::batch_bytesT;

    /// Issue the non-blocking find for a batch record on its owner (Step 2b).

    /// Returns a future to the container iterator; the find AM is in flight on
    /// return, so this can be issued during the previous task's compute to hide
    /// the round-trip. Resolves immediately when this rank owns the record.
    /// Pair with deserialize_batch to obtain the vecfunc once the future is set.
    Future<batch_iterator> fetch_batch_record_async(const keyT record) const {
        if (batch_debug) {
            const ProcessID owner = batch_pmap->owner(record);
            fprintf(stderr, "FIXB_DBG cloud.fetch_batch_record_async ENTER key=0x%lx owner=%d inflight_pre=%ld t=%.6f\n",
                    static_cast<unsigned long>(record), owner, inflight_finds_.load(), wall_time());
            fflush(stderr);
        }
        ++inflight_finds_;
        Future<batch_iterator> f = batch_container.find(record);
        if (batch_debug) {
            fprintf(stderr, "FIXB_DBG cloud.fetch_batch_record_async EXIT key=0x%lx inflight_now=%ld t=%.6f\n",
                    static_cast<unsigned long>(record), inflight_finds_.load(), wall_time());
            fflush(stderr);
        }
        return f;
    }

    /// Whether per-call BATCH_* diagnostic prints are enabled. Exposed so the
    /// caller (e.g., exchange's prefetch hook) can emit its own log lines that
    /// include subworld/rank context not visible inside the cloud.
    bool is_batch_debug() const { return batch_debug; }

    /// Look up the owner of a batch record (no AM). Used by the exchange's
    /// prefetch hook so its BATCH_PREFETCH log line can include the target rank.
    ProcessID batch_owner(const keyT record) const {
        return batch_pmap->owner(record);
    }

    /// Snapshot of in-flight async finds (issued by fetch_batch_record_async,
    /// not yet consumed by deserialize_batch). Diagnostic only.
    long inflight_finds() const { return inflight_finds_.load(); }

    /// Deserialize a batch from a prefetched find future (Step 2b).

    /// Blocks on `fut` only if the find round-trip has not completed yet, then
    /// deserializes the fetched bytes locally. The stored bytes are a plain
    /// vector-archive stream (see store_batch), so a ParallelInputArchive over a
    /// VectorInputArchive reproduces the functions exactly as fetch_batch's
    /// ContainerRecordInputArchive would. Caches by record like fetch_batch.
    template<typename T, std::size_t NDIM>
    std::vector<Function<T, NDIM>> deserialize_batch(madness::World& subworld,
            Future<batch_iterator> fut, const keyT record,
            const bool cache_result = true) const {
        typedef std::vector<Function<T, NDIM>> vecfuncT;
        if (is_cached(record)) return load_from_cache<vecfuncT>(subworld, record);
        const double wall0 = wall_time();
        cloudtimer t(subworld, reading_time);
        const ProcessID owner = batch_pmap->owner(record);
        if (batch_debug) {
            fprintf(stderr, "FIXB_DBG cloud.deserialize_batch BEFORE_GET key=0x%lx owner=%d sw_rank=%d sw_id=%lu inflight=%ld t=%.6f\n",
                    static_cast<unsigned long>(record), owner, subworld.rank(),
                    static_cast<unsigned long>(subworld.id()), inflight_finds_.load(), wall_time());
            fflush(stderr);
        }
        // fut.get() blocks only if the find round-trip is still outstanding;
        // with proper prefetch overlap this should be ~0 for remote records.
        const double w0 = wall_time();
        batch_iterator it = fut.get();
        const double w1 = wall_time();
        --inflight_finds_;
        if (batch_debug) {
            fprintf(stderr, "FIXB_DBG cloud.deserialize_batch AFTER_GET key=0x%lx owner=%d find_ms=%.3f inflight_after=%ld t=%.6f\n",
                    static_cast<unsigned long>(record), owner, (w1-w0)*1e3,
                    inflight_finds_.load(), wall_time());
            fflush(stderr);
        }
        const double find_dt = w1 - w0;
        batch_find_time += long(find_dt * 1e6);
        MADNESS_CHECK_THROW(it != batch_container.end(), "deserialize_batch: record not found");
        valueT bytes = it->second;       // copy out the serialized stream
        vecfuncT batch;
        {
            madness::archive::VectorInputArchive var(bytes);
            madness::archive::ParallelInputArchive<madness::archive::VectorInputArchive> par(subworld, var);
            std::size_t fsize = 0;
            par & fsize;
            batch.resize(fsize);
            for (std::size_t i = 0; i < fsize; ++i) par & batch[i];
        }
        const double deser_dt = wall_time() - w1;
        batch_deserialize_time += long(deser_dt * 1e6);
        if (use_cache and cache_result) cache(subworld, batch, record);
        if (batch_debug) {
            std::ostringstream os;
            os << "BATCH_DESER rank=" << subworld.rank()
               << " sw=" << subworld.id()
               << " t=" << std::fixed << std::setprecision(3) << wall_time()
               << " key=0x" << std::hex << record << std::dec
               << " owner=" << owner
               << " kind=prefetched"
               << " nfuncs=" << batch.size()
               << " find_ms=" << std::fixed << std::setprecision(3) << (find_dt*1e3)
               << " deser_ms=" << std::fixed << std::setprecision(3) << (deser_dt*1e3)
               << " wall_ms=" << std::fixed << std::setprecision(3) << ((wall_time()-wall0)*1e3);
            print(os.str());
        }
        return batch;
    }

    // ---- owner-pinned batch point-to-point fetch (BatchTransport) ----
    // New fetch path that replaces the WorldContainer::find round-trip
    // (fetch_batch_record_async + deserialize_batch). Not yet wired into the
    // exchange; the old path is retained until the swap.

    /// Read a batch record's serialized bytes from this rank's local batch store.
    /// Precondition: this rank owns `record` (batch_pmap->owner(record)==rank).
    /// Used by BatchTransport::do_send on the owner side.
    valueT get_local_batch_bytes(const keyT record) const {
        typename madness::WorldContainer<keyT, valueT>::const_accessor acc;
        if (batch_container.find(acc, record)) return acc->second;   // copy out of the local map
        MADNESS_EXCEPTION("get_local_batch_bytes: record not held locally", int(record));
    }

    /// Stable pointer + size of a local batch record's bytes, for a zero-copy
    /// Isend from BatchTransport::on_trigger on the comm thread (copying 0.6 GB
    /// on the comm thread would stall RMI). The accessor lock is released on
    /// return, but the pointer stays valid because batch_container entries are
    /// not erased or mutated between store_batches and the end of the exchange.
    /// Precondition: this rank owns `record`.
    /// (Confirmation-prototype shortcut for the not-yet-built owner_store_, which
    /// will formalize this stable-address guarantee.)
    std::pair<const unsigned char*, std::size_t> get_local_batch_ptr(const keyT record) const {
        typename madness::WorldContainer<keyT, valueT>::const_accessor acc;
        if (batch_container.find(acc, record))
            return {acc->second.data(), acc->second.size()};
        MADNESS_EXCEPTION("get_local_batch_ptr: record not held locally", int(record));
    }

    /// Issue the point-to-point fetch of `record` from its owner (Step 2b, p2p).
    /// Returns a future to the serialized bytes; the trigger is in flight on
    /// return so the round-trip overlaps the caller's compute. Pair with
    /// deserialize_batch_p2p to obtain the vecfunc. Mirror of
    /// fetch_batch_record_async, but transports raw bytes instead of an iterator.
    Future<batch_bytesT> request_batch_bytes_async(const keyT record) const {
        return batch_transport_->request(record);
    }

    /// Deserialize a batch from a p2p byte future (Step 2b, p2p).
    /// Blocks on `fut` only if the transfer has not landed yet, then deserializes
    /// the received bytes locally. Mirror of deserialize_batch, fed by raw bytes.
    template<typename T, std::size_t NDIM>
    std::vector<Function<T, NDIM>> deserialize_batch_p2p(madness::World& subworld,
            Future<batch_bytesT> fut, const keyT record,
            const bool cache_result = true) const {
        typedef std::vector<Function<T, NDIM>> vecfuncT;
        if (is_cached(record)) return load_from_cache<vecfuncT>(subworld, record);
        cloudtimer t(subworld, reading_time);
        const double w0 = wall_time();
        batch_bytesT bytes = fut.get();                     // await the transfer
        const double w1 = wall_time();
        batch_find_time += long((w1 - w0) * 1e6);
        vecfuncT batch;
        {
            madness::archive::VectorInputArchive var(bytes);
            madness::archive::ParallelInputArchive<madness::archive::VectorInputArchive> par(subworld, var);
            std::size_t fsize = 0;
            par & fsize;
            batch.resize(fsize);
            for (std::size_t i = 0; i < fsize; ++i) par & batch[i];
        }
        batch_deserialize_time += long((wall_time() - w1) * 1e6);
        if (use_cache and cache_result) cache(subworld, batch, record);
        return batch;
    }

    /// Synchronous p2p batch fetch (Step 2b, p2p). Drop-in for fetch_batch: cache
    /// hit short-circuits; otherwise request the bytes from the owner and
    /// deserialize (the request is skipped entirely on a cache hit). Resolves
    /// locally without MPI when this rank owns the record.
    template<typename T, std::size_t NDIM>
    std::vector<Function<T, NDIM>> fetch_batch_p2p(madness::World& subworld,
            const keyT record, const bool cache_result = true) const {
        typedef std::vector<Function<T, NDIM>> vecfuncT;
        if (is_cached(record)) return load_from_cache<vecfuncT>(subworld, record);
        return deserialize_batch_p2p<T, NDIM>(subworld, request_batch_bytes_async(record),
                                              record, cache_result);
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
public:
    mutable std::atomic<long> copy_time=0l;        // if pointers are stored in cloud, time to copy from universe to subworld
    mutable std::atomic<long> target_replication_time=0l;     // if pointers are stored in cloud, time to replicate targets
private:
    mutable std::atomic<long> writing_time=0l;     // in microseconds
    mutable std::atomic<long> writing_time1=0l;    // in microseconds
    mutable std::atomic<long> replication_time=0l;    // in microseconds
    mutable std::atomic<long> cache_reads=0l;
    mutable std::atomic<long> cache_stores=0l;
    // ---- StoreFunctionBatched diagnostics (the cloud-batch fetch path) ----
    mutable std::atomic<long> batch_store_time=0l;       ///< store_batch wall time (microseconds; subset of writing_time1)
    mutable std::atomic<long> batch_find_time=0l;        ///< wall time inside fetch_batch / deserialize_batch waiting on the find round-trip (microseconds)
    mutable std::atomic<long> batch_deserialize_time=0l; ///< wall time deserializing the bytes (microseconds; the CPU half of reading)
    mutable std::atomic<long> batch_prefetch_hits=0l;    ///< exchange's row-owner k-batch prefetch hits (set by tally_batch_prefetch_hit)
    mutable std::atomic<long> batch_prefetch_misses=0l;  ///< exchange's row-owner k-batch prefetch misses (set by tally_batch_prefetch_miss)
    mutable std::atomic<long> inflight_finds_=0l;        ///< async finds issued by fetch_batch_record_async but not yet consumed

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

// ---- BatchTransport out-of-line definitions (need the complete Cloud) ----

inline Future<batch_bytesT> BatchTransport::request(batch_keyT record) {
    World& u = this->get_world();
    const ProcessID owner = cloud_->batch_owner(record);
    if (owner == u.rank()) {
        // local: no MPI, just hand back a copy of the local bytes
        return Future<batch_bytesT>(cloud_->get_local_batch_bytes(record));
    }
    const int tag = alloc_tag();
    auto p = std::make_shared<PendingRecv>();
    p->owner = owner;
    p->tag = tag;
    {
        std::lock_guard<std::mutex> g(pending_mtx_);
        pending_[tag] = p;
    }
    // trigger runs inline on the owner's comm thread (send, not task), so the
    // owner posts its Isend without queueing behind its saturated workers.
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
    // comm-thread handler: post the bytes Isend straight from the stable local
    // buffer (no copy), then reply with the size. Must not block.
    World& u = this->get_world();
    reap_sends();
    auto ptr_size = cloud_->get_local_batch_ptr(record);
    const std::size_t n = ptr_size.second;
    SafeMPI::Request rd = u.mpi.Isend(ptr_size.first, int(n), MPI_BYTE, requester, tag);
    {
        std::lock_guard<std::mutex> g(sends_mtx_);
        sends_.push_back(rd);
    }
    // reply carries the size so the requester can post its Irecv now (during the
    // requester's compute), enabling the rendezvous to complete in the background.
    this->send(requester, &BatchTransport::on_reply, tag, n);
}

inline void BatchTransport::on_reply(int tag, std::size_t size) {
    // comm-thread handler: allocate the buffer and post the Irecv, then enqueue
    // the worker-side completion. Must not block.
    World& u = this->get_world();
    std::shared_ptr<PendingRecv> p;
    {
        std::lock_guard<std::mutex> g(pending_mtx_);
        auto it = pending_.find(tag);
        MADNESS_CHECK(it != pending_.end());
        p = it->second;
    }
    p->buf.resize(size);
    p->req = u.mpi.Irecv(p->buf.data(), int(size), MPI_BYTE, p->owner, tag);
    u.taskq.add(this, &BatchTransport::finish_recv, tag);
}

inline void BatchTransport::finish_recv(int tag) {
    // worker task (post-compute): the Irecv has progressed in the background on
    // the comm thread, so this await returns promptly; then set the future.
    std::shared_ptr<PendingRecv> p;
    {
        std::lock_guard<std::mutex> g(pending_mtx_);
        auto it = pending_.find(tag);
        MADNESS_CHECK(it != pending_.end());
        p = it->second;
        pending_.erase(it);
    }
    World::await(p->req, true);
    p->fut.set(std::move(p->buf));
}

} /* namespace madness */

#endif /* SRC_MADNESS_WORLD_CLOUD_H_ */
