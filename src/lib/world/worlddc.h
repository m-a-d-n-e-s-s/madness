#ifndef WORLD_DC
#define WORLD_DC

/// \file worlddc.h
/// \brief Implements WorldContainer



/*
  
Current class provides hash_map-like functinonality with caching.
Other containers will follow when this design stabilizes.

To do (in rough order of importance):

c) Straighten out semantics for caching and test it

e) Clarification of exception handling

f) Do we need a global iterator?  Not yet for MRA

g) Global clear operation?  Can just assign from a new container so no
rush.

*/


namespace madness {
        
    template <typename keyT, 
              typename valueT, 
              typename procmapT,
              typename attrT>
    class WorldContainer;
    
    template <typename keyT, 
              typename valueT, 
              typename procmapT,
              typename attrT>
    class WorldContainerImpl;
    
    
    /// Defines attributes of the container, esp. caching
    struct DCDefaultAttr {
        // There are only three sequentially consistent policies supported.
        // a) No read or write caching.  Every access is remote => false,false.
        // b) Cache reads. Writes invalidate locally, write remotely => true,false.
        // c) Cache reads. Writes write locally and remotely => true,true.
        
        
        // !! CURRENTLY, only false,false is supported
        
        static const bool CacheReadPolicy = false;
        static const bool CacheWritePolicy = false;
        static const size_t CacheMaxEntries = 1024;
    };
    
    
    /// Default process hash is "random" using madness::hash(key)
    template <typename keyT>
    class DCDefaultProcmap {
    private:
        int nproc;
        unsigned int mask;
        
    public:
        DCDefaultProcmap(World& world) {
            nproc = world.mpi.nproc();
            // Find bit mask to eliminate use of modulo operator which relies
            // upon having a VERY GOOD hash, which we do.
            mask = 1;
            while (mask < (unsigned) nproc) mask <<= 1;
            mask--;
        };
        
        // Assumes key is a contiguous byte stream
        ProcessID operator()(const keyT& key) const {
            if (nproc == 1) return 0;
            hashT hh = hash(key);
            int h = hh&mask;
            if (h >= nproc) h -= nproc;
            return h;
        };
    };
    
    
    // Used to store info about stuff in the cache
    struct CacheInfo {
        mutable int nref;  // Reference count
        mutable std::vector<CallbackInterface*> callbacks;
        
        CacheInfo(int nref = 0) : nref(nref), callbacks() {};
        
        // Iterators register their interest here
        void register_callback(CallbackInterface* callback) {
            nref++;
            callbacks.push_back(callback);
        };
        
        // Iterators unregister their interest here
        void unregister_callback(CallbackInterface* callback) {
            nref--;
            for (std::size_t i=0; i<callbacks.size(); i++) {
                if (callbacks[i] == callback) {
                    callbacks[i] = 0;
                    break;
                }
            }
        };
        
        // Destructor invalidates referencing iterators
        ~CacheInfo() {
            for (std::size_t i=0; i<callbacks.size(); i++) {
                if (callbacks[i]) callbacks[i]->notify();
            }
        };
    };
    
    /// Iterator for distributed container wraps the local iterator
    
    /// Reference counting is performed for data in the cache so that
    /// if you copy the iterator its reference count is also
    /// incremented.  This is so that we can autoclean the cache
    /// without invalidating your references.  Thus, do NOT take a
    /// direct pointer/reference to a cached value since if you then
    /// free the iterator your pointer may well be junk.
    ///
    /// However, the reference counting still leaves the usual problem
    /// (present in even the STL containers) that if you erase an
    /// entry (either via the iterator or via the key) that your
    /// iterators are now invalid.  Just the way it is.
    template <class implT, class internal_iteratorT, class pairT>
    class WorldContainerIterator : private CallbackInterface {
    private:
        typedef typename implT::cacheinfo_iteratorT cacheinfo_iteratorT;
        typedef typename implT::attributes attrT;
        
        implT* impl;                   //< Pointer to container implementation
        internal_iteratorT  it;         //< Iterator from local/cache container
        cacheinfo_iteratorT cacheit;   //< Iterator for reference counting cache
        bool fromcache;                //< True if is a cached remote value
        
        // Called by the cache to notify iterators to invalidate themselves
        void notify() {
            print("Invalidating iterator!",fromcache);
            impl = 0;
            it = internal_iteratorT();
        };
        
        inline void register_cache_callback() {
            if (impl && fromcache) cacheit->second.register_callback(static_cast<CallbackInterface*>(this));
        };
        
        inline void unregister_cache_callback() {
            if (impl && fromcache) {
                CacheInfo& info = cacheit->second;
                info.unregister_callback(static_cast<CallbackInterface*>(this));
                MADNESS_ASSERT(info.nref>=0);
                if (!attrT::CacheReadPolicy && info.nref==0) {
                    //print("Iterator deleting cached data");
                    impl->cache_erase(it,cacheit);
                }
            }
        };
        
        
    public:
        /// Default constructor (same as end())
        WorldContainerIterator() 
            : impl(0), it(), cacheit(), fromcache(false)
        {};
        
        /// Private: Constructor used internally
        WorldContainerIterator(implT* impl, 
                                     const internal_iteratorT& it, 
                                     const cacheinfo_iteratorT& cacheit,
                                     bool fromcache) 
            : impl(impl), it(it), cacheit(cacheit), fromcache(fromcache) 
        {
            register_cache_callback();
        };
        
        /// Copy constructor (increments reference count)
        WorldContainerIterator(const WorldContainerIterator& other) 
            : impl(other.impl), it(other.it), cacheit(other.cacheit), fromcache(other.fromcache)
        {
            register_cache_callback();
        };
        
        /// Assignment (increments reference count)
        WorldContainerIterator& operator=(const WorldContainerIterator& other) {
            if (this != &other) {
                unregister_cache_callback();  // <---
                impl = other.impl;
                it = other.it;
                cacheit = other.cacheit;
                fromcache = other.fromcache;
                register_cache_callback();   // <---
            }
            return *this;
        };
        
        /// Determines if two iterators are identical
        bool operator==(const WorldContainerIterator& other) const {
            return it==other.it;
        };
        
        /// Determines if two iterators are different
        bool operator!=(const WorldContainerIterator& other) const {
            return it != other.it;
        };
        
        /// Pre-increment of an iterator (i.e., ++it) --- \em local iterators only
        
        /// Trying to increment a remote iterator will throw
        WorldContainerIterator& operator++() {
            if (impl) {
                MADNESS_ASSERT(!fromcache);
                if (++it == impl->local.end()) *this = impl->end();
            }
            return *this;
        };
        
        /// Iterators dereference to std::pair<keyT,valueT>
        const pairT* operator->() const {
            return it.operator->();
        };
        
        /// Iterators dereference to std::pair<keyT,valueT>
        pairT* operator->() {
            return it.operator->();
        };
        
        /// Iterators dereference to std::pair<keyT,valueT>
        const pairT& operator*() const {
            return *it;
        };
        
        /// Iterators dereference to std::pair<keyT,valueT>
        pairT& operator*() {
            return *it;
        };
        
        /// Private: (or should be) Returns iterator of internal container
        const internal_iteratorT& get_internal_iterator() const {
            return it;
        };
        
        /// Returns true if this is a cached remote value
        bool is_cached() const {
            return fromcache;
        };
        
        ~WorldContainerIterator() {
            unregister_cache_callback();
        };
        
        template <typename Archive>
        void serialize(const Archive& ar) {
            throw "Serializing DC iterator ... why?";
        }
    };
    
    
    
    
    /// Implementation of distributed container to enable PIMPL
    template <typename keyT, 
              typename valueT, 
              typename procmapT,
              typename attrT>
    class WorldContainerImpl 
        : public WorldObject< WorldContainerImpl<keyT, valueT, procmapT, attrT> >
          , private NO_DEFAULTS 
    {
    public:
        typedef std::pair<const keyT,valueT> pairT;
        typedef const std::pair<const keyT,valueT> const_pairT;
        typedef WorldContainerImpl<keyT,valueT,procmapT,attrT> implT;
        
#ifdef WORLDDC_USES_GNU_HASH_MAP
        template <typename T> 
        struct DCLocalHash {
            std::size_t operator()(const T& t) const {
                return hash(t);
            };
        };
        typedef HASH_MAP_NAMESPACE::hash_map< keyT,CacheInfo,DCLocalHash<keyT> > cacheinfoT;
        typedef HASH_MAP_NAMESPACE::hash_map< keyT,valueT,DCLocalHash<keyT> > internal_containerT;
#else
        typedef std::map<keyT,int> cacheinfoT;
        typedef std::map<keyT,valueT> internal_containerT;
#endif
        
        typedef typename cacheinfoT::iterator cacheinfo_iteratorT;
        typedef typename internal_containerT::iterator internal_iteratorT;
        typedef typename internal_containerT::const_iterator internal_const_iteratorT;
        typedef WorldContainerIterator<implT,internal_iteratorT,pairT> iteratorT;
        typedef WorldContainerIterator<implT,internal_iteratorT,pairT> iterator;
        typedef WorldContainerIterator<const implT, internal_const_iteratorT, const_pairT> const_iteratorT;
        typedef WorldContainerIterator<const implT, internal_const_iteratorT, const_pairT> const_iterator;
        typedef attrT attributes;
        
        friend class WorldContainer<keyT,valueT,procmapT,attrT>;
        friend class WorldContainerIterator<implT,internal_iteratorT,pairT>;
        friend class WorldContainerIterator<const implT,internal_const_iteratorT,const_pairT>;
        
    private:
        
        WorldContainerImpl();   // Inhibit default constructor
        
        World& world;
        const uniqueidT theid;    //< Universe-wide unique ID for this instance
        const procmapT procmap;       //< Function/class to map from keys to owning process
        const ProcessID me;           //< My MPI rank
        internal_containerT local;    //< Locally owned data
        mutable internal_containerT cache;    //< Local cache for remote data
        mutable cacheinfoT cacheinfo;         //< Maps cached keys to info
        const iterator end_iterator;          //< For fast return of end
        const const_iterator end_const_iterator; //< For fast return of end
        
        /// Removes (remote) item from local cache
        void cache_erase (const internal_const_iteratorT& it, const cacheinfo_iteratorT& cacheit) const {
            cache.erase(it->first);
            cacheinfo.erase(cacheit);
        };
        
        
        /// Handles find request
        void find_handler(ProcessID requestor, const keyT& key, const RemoteReference< FutureImpl<iterator> >& ref) {
            if (owner(key) != me) {
                send(owner(key), &implT::find_handler, requestor, key, ref);
            }
            else {
                internal_iteratorT r = local.find(key);
                if (r == local.end()) send(requestor, &implT::find_failure_handler, ref);
                else send(requestor, &implT::find_success_handler, ref, *r);
            }
        };
        
        /// Handles successful find response
        void find_success_handler(const RemoteReference< FutureImpl<iterator> >& ref, const pairT& datum) {
            FutureImpl<iterator>* f = ref.get();
            // Check if the value is already in the cache.  If so, overwrite it with 
            // the new value ... but don't invalidate any old iterators in the process.
            cacheinfo_iteratorT cacheit;
            internal_iteratorT it = cache.find(datum.first);
            if (it == cache.end()) {
                it = cache.insert(datum).first;
                cacheit = cacheinfo.insert(std::pair<keyT,CacheInfo>(datum.first,CacheInfo(0))).first;
            }
            else {
                it->second = datum.second;
                cacheit = cacheinfo.find(datum.first);
            }
            f->set(iterator(this, it, cacheit, true));
            ref.dec(); // Matching inc() in find() where ref was made
        };
        
        /// Handles unsuccessful find response
        void find_failure_handler(const RemoteReference< FutureImpl<iterator> >& ref) {
            FutureImpl<iterator>* f = ref.get();
            f->set(end());
            ref.dec(); // Matching inc() in find() where ref was made
        };
        
    public:
        
        WorldContainerImpl(World& world, const procmapT& procmap)
            : WorldObject< WorldContainerImpl<keyT, valueT, procmapT, attrT> >(world)
            , world(world)
            , theid(world.register_ptr(this))
            , procmap(procmap)
            , me(world.mpi.rank())
            , local()
            , cache()
            , cacheinfo()
        {
            this->process_pending();
        };
        
        
        ~WorldContainerImpl() {
            print("In DCImpl destructor");
        };
        
        bool is_local(const keyT& key) const {
            return owner(key) == me;
        };
        
        ProcessID owner(const keyT& key) const {
            return procmap(key);
        };
        
        bool probe(const keyT& key) const {
            ProcessID dest = owner(key);
            if (dest == me) 
                return local.find(key) != local.end();
            else if (attrT::CacheReadPolicy) 
                return cache.find(key) != cache.end();
            else 
                return false;
        };
        
        std::size_t size() const {
            return local.size();
        };
        
        void handle_cache_overflow() {};
        
        void insert(const pairT& datum) {
            ProcessID dest = owner(datum.first);
            if (dest == me) {
                local.insert(datum);
            }
            else {
                if (attrT::CacheWritePolicy) {  // Remote writes also hit cache
                    if (cache.size() >= attrT::CacheMaxEntries) 
                        handle_cache_overflow();
                    cache.insert(datum);
                    cacheinfo.insert(std::pair<keyT,CacheInfo>(datum.first,0));
                }
                else if (attrT::CacheReadPolicy) { // Remote writes only invalidate cache
                    cache.erase(datum.first);
                    cacheinfo.erase(datum.first);
                }
                send(dest, &implT::insert, datum);
            }
        };
        
        
        void clear() {
            local.clear();
            cache.clear();
            cacheinfo.clear();
        };
        
        
        void erase(const keyT& key) {
            ProcessID dest = owner(key);
            if (dest == me) {
                local.erase(key);
            }
            else {
                cache.erase(key);
                cacheinfo.erase(key);
                void (implT::*eraser)(const keyT&);
                send(dest, eraser, key);
            }                
        };
        
        
        void erase(const iterator& it) {
            if (it == end()) {
                return;
            }
            else if (it.is_cached()) {
                erase(it->first);
            }
            else { 
                local.erase(it.get_internal_iterator());
            }
        };
        
        void erase(const iterator& start, const iterator& finish) {
            iterator it = start;
            while (it!=finish && it!=end()) {
                iterator last=it;
                ++it;
                erase(last);
            }
        };
        
        iterator begin() {
            internal_iteratorT it = local.begin();
            if (it == local.end()) return end();
            return iterator(this,it,cacheinfo_iteratorT(),false);
        };
        
        const_iterator begin() const {
            internal_const_iteratorT it = local.begin();
            if (it == local.end()) return end();
            return const_iterator(const_cast<const implT*>(this),
                                  it, cacheinfo_iteratorT(),false);
        };
        
        const iterator& end() {
            return end_iterator;
        };
        
        const const_iterator& end() const {
            return end_const_iterator;
        };
        
        
        Future<iterator> find(const keyT& key) {
            ProcessID dest = owner(key);
            if (dest == me) {  // Local read
                internal_iteratorT it = local.find(key);
                if (it == local.end()) return Future<iterator>(end());
                else return Future<iterator>(iterator(this,it,cacheinfo_iteratorT(),false));
            }
            else {             // Remote read
                if (attrT::CacheReadPolicy) {  // Try cache
                    internal_iteratorT it = cache.find(key);
                    if (it != cache.end()) {
                        return Future<iterator>(iterator(this,it,cacheinfo.find(key),true));
                    }
                }
                Future<iterator> result;
                send(dest, &implT::find_handler, me, key, result.remote_ref(world));
                return result;
            }
        };
        
        // Used to forward call to item member function
        template <typename memfunT>
        RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT))
            itemfun(const keyT& key, memfunT memfun) {
            return detail::MemfunRun<MEMFUN_RETURNT(memfunT)>::run(local[key], memfun);
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T>
        RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT))        
            itemfun(const keyT& key, memfunT memfun, const arg1T& arg1) {
            return detail::MemfunRun<MEMFUN_RETURNT(memfunT)>::run(local[key], memfun, arg1);
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T>
        RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT))        
            itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            return detail::MemfunRun<MEMFUN_RETURNT(memfunT)>::run(local[key], memfun, arg1, arg2);
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT))        
            itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            return detail::MemfunRun<MEMFUN_RETURNT(memfunT)>::run(local[key], memfun, arg1, arg2, arg3);
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T>
        RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT))        
            itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg3T& arg4) {
            return detail::MemfunRun<MEMFUN_RETURNT(memfunT)>::run(local[key], memfun, arg1, arg2, arg3, arg4);
        }
        
        // Used to forward call to item member function
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T, typename arg4T, typename arg5T>
        RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT))        
            itemfun(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3, const arg3T& arg4, const arg5T& arg5) {
            return detail::MemfunRun<MEMFUN_RETURNT(memfunT)>::run(local[key], memfun, arg1, arg2, arg3, arg4, arg5);
        }
        
    };
    
    
    /// Makes a distributed container with specified attributes
    
    /// There is no communication or syncronization associated with
    /// making a new container, but every process must invoke the
    /// constructor for each container in the same order.  This is so
    /// that we can assign each container a unique ID without any
    /// communication.  Since remotely invoked operations may start
    /// happening before local construction, messages on not yet 
    /// constructed containers are buffered pending construction.
    ///
    /// Similarly, when a container is destroyed, the actual
    /// destruction is deferred until a synchronization point
    /// (world.gop.fence()) in order to eliminate the need to fence
    /// before and after destroying every container.
    ///
    /// The caching behavior is controlled by the attributes class.
    /// Currently, this is untested and no caching is enabled.
    ///
    /// The distribution of data between processes is controlled by
    /// the process map (ProcMap) class.  The default is uniform
    /// hashing based upon a strong (Bob Jenkins, lookup3) bytewise
    /// hash of the key.
    ///
    /// All operations, including constructors and destructors, are
    /// non-blocking and return immediately.  If communication occurs
    /// it is asynchronous, otherwise operations are local.
    ///
    /// !! The key is presently assumed to be small, probably less
    /// than 100 bytes.  This can be relaxed as the need arises.
    template <typename keyT, 
              typename valueT, 
              typename procmapT = DCDefaultProcmap<keyT>,
              typename attrT = DCDefaultAttr>
    class WorldContainer {
    public:
        typedef WorldContainer<keyT,valueT,procmapT,attrT> containerT;
        typedef WorldContainerImpl<keyT,valueT,procmapT,attrT> implT;
        typedef typename implT::pairT pairT;
        typedef typename implT::iterator iterator;
        typedef typename implT::const_iterator const_iterator;
        typedef Future<iterator> futureT;
        typedef Future<iterator> future;
        typedef Future<const_iterator> const_future;
        
    private:
        SharedPtr<implT> p;
        
        inline void check_initialized() const {
            MADNESS_ASSERT(p);
        };
    public:
        
        /// Makes an uninitialized container (no communication)
        
        /// The container is useless until assigned to from a fully
        /// constructed container.  There is no need to worry about
        /// default constructors being executed in order.
        WorldContainer() 
            : p(0)
        {};
        
        
        /// Makes a container with default data distribution (no communication)
        
        /// A unique ID is associated with every distributed container
        /// within a world.  In order to avoid synchronization when
        /// making a container, we have to assume that all processes
        /// execute this constructor in the same order (does not apply
        /// to the non-initializing, default constructor).
        WorldContainer(World& world) 
            : p(new implT(world,procmapT(world)))
        {
            world.deferred_cleanup(p);
        };
        
        /// Makes an initialized, empty container (no communication)
        
        /// A unique ID is associated with every distributed container
        /// within a world.  In order to avoid synchronization when
        /// making a container, we have to assume that all processes
        /// execute this constructor in the same order (does not apply
        /// to the non-initializing, default constructor).
        WorldContainer(World& world, const procmapT& procmap) 
            : p(new implT(world,procmap))
        {
            world.deferred_cleanup(p);
        };
        
        
        /// Copy constructor is shallow (no communication)
        
        /// The copy refers to exactly the same container as other
        /// which must be initialized.
        WorldContainer(const WorldContainer& other) 
            : p(other.p)
        {
            check_initialized();
        };
        
        /// Assignment is shallow (no communication)
        
        /// The copy refers to exactly the same container as other
        /// which must be initialized.
        containerT& operator=(const containerT& other) {
            if (this != &other) {
                other.check_initialized();
                p = other.p;
            }
            return *this;
        };
        
        /// Returns the world associated with this container
        World& world() {
            check_initialized();
            return p->world;
        };
        
        
        /// Inserts key+value pair (non-blocking communication if key not local)
        void insert(const std::pair<keyT,valueT>& datum) {
            check_initialized();
            p->insert(datum);
        };
        
        
        /// Inserts key+value pair (non-blocking communication if key not local)
        void insert(const keyT& key, const valueT& value) {
            insert(pairT(key,value));
        };
        
        
        /// Inserts pairs (non-blocking communication if key(s) not local)
        
        /// See insert(pair) for cache semantics.
        template <typename input_iterator>
        void insert(input_iterator& start, input_iterator& end) {
            check_initialized();
            for_each(start,end,bind1st(mem_fun(&containerT::insert),this));
        }
        
        
        
        /// Returns true if local data is immediately available (no communication)
        void probe(const keyT& key) const {
            check_initialized();
            return p->probe(key);
        };
        
        
        /// Returns processor that logically owns key (no communication)
        
        /// Local remapping may have changed its physical location, but all
        /// operations should forward correctly. 
        inline ProcessID owner(const keyT& key) const {
            check_initialized();
            return p->owner(key);
        };
        
        
        /// Returns true if the key maps to the local processor (no communication)
        bool is_local(const keyT& key) const {
            check_initialized();
            return p->is_local();
        };
        
        
        /// Returns a future iterator (non-blocking communication if key not local)
        
        /// Like an std::map an iterator "points" to an std::pair<const keyT,valueT>.
        ///
        /// Refer to Future for info on how to avoid blocking.
        Future<iterator> find(const keyT& key) const {
            check_initialized();
            return p->find(key);
        };
        
        
        /// Returns an iterator to the beginning of the \em local data (no communication)
        iterator begin() {
            check_initialized();
            return p->begin();
        };
        
        /// Returns an iterator to the beginning of the \em local data (no communication)
        const_iterator begin() const {
            check_initialized();
            return const_cast<const implT*>(p.get())->begin();
        };
        
        /// Returns an iterator past the end of the \em local data (no communication)
        const iterator& end() {
            check_initialized();
            return p->end();
        };
        
        /// Returns an iterator past the end of the \em local data (no communication)
        const const_iterator& end() const {
            check_initialized();
            return const_cast<const implT*>(p.get())->end();
        };
        
        /// Erases entry from container (non-blocking comm if remote)
        
        /// Missing keys are quietly ignored.
        ///
        /// Note that erasing an entry may invalidate iterators on the
        /// remote end.  This is just the same as what happens when
        /// using STL iterators on an STL container in a sequential
        /// algorithm.  To protect against such an intrusion you
        /// should consider creating around your use of local
        /// iterators a critical section either by locally suspending
        /// processing of AM and/or tasks, or by bracketting with a
        /// global_fence.
        void erase(const keyT& key) {
            check_initialized();
            p->erase(key);
        };
        
        /// Erases entry corresponding to \em local iterator (no communication)
        void erase(const iterator& it) {
            check_initialized();
            p->erase(it);
        };            
        
        /// Erases range defined by \em local iterators (no communication)
        void erase(const iterator& start, const iterator& finish) {
            check_initialized();
            p->erase(start,finish);
        };
        
        
        /// Clears all \em local data (no communication)
        
        /// Invalidates all iterators
        void clear() {
            check_initialized();
            p->clear();
        };
        
        /// Returns the number of \em local entries (no communication)
        std::size_t size() const {
            check_initialized();
            return p->size();
        };
        

        /// Sends message "resultT memfun()" to item (non-blocking comm if remote)
        
        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for remote messages.
        /// 
        /// Returns a future result (Future<void> may be ignored).
        template <typename memfunT>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(const keyT& key, memfunT memfun) {
            check_initialized();
            RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT)) (implT::*itemfun)(const keyT&, memfunT) = &implT:: template itemfun<memfunT>;
            return p->send(owner(key), itemfun, key, memfun);
        }
        
        
        /// Sends message "resultT memfun(arg1T)" to item (non-blocking comm if remote)
        
        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for remote messages.
        /// 
        /// Returns a future result (Future<void> may be ignored).
        template <typename memfunT, typename arg1T>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(const keyT& key, memfunT memfun, const arg1T& arg1) {
            check_initialized();
            RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT)) (implT::*itemfun)(const keyT&, memfunT, const arg1T&) = &implT:: template itemfun<memfunT,arg1T>;
            return p->send(owner(key), itemfun, key, memfun, arg1);
        }
        
        
        /// Sends message "resultT memfun(arg1T,arg2T)" to item (non-blocking comm if remote)
        
        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for both local and remote messages.
        /// 
        /// Returns a future result (Future<void> may be ignored).
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            check_initialized();
            RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT)) (implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&) = &implT:: template itemfun<memfunT,arg1T,arg2T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2);
        }
        
        
        /// Sends message "resultT memfun(arg1T,arg2T)" to item (non-blocking comm if remote)
        
        /// If item does not exist it is made with the default constructor.
        ///
        /// Future arguments must be ready for both local and remote messages.
        /// 
        /// Returns a future result (Future<void> may be ignored).
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< MEMFUN_RETURNT(memfunT) > 
        send(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            check_initialized();
            RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT)) (implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T>;
            return p->send(owner(key), itemfun, key, memfun, arg1, arg2, arg3);
        }
        
        
        /// Adds task "resultT memfun()" in process owning item (non-blocking comm if remote)
        
        /// If item does not exist it is made with the default constructor.  
        /// 
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        /// 
        /// Returns a future result (Future<void> may be ignored).
        template <typename memfunT>
        Future< MEMFUN_RETURNT(memfunT) >
        task(const keyT& key, memfunT memfun) {
            check_initialized();
            RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT)) (implT::*itemfun)(const keyT&, memfunT) = &implT:: template itemfun<memfunT>;
            return p->task(owner(key), itemfun, key, memfun);
        }
        
        /// Adds task "resultT memfun(arg1T)" in process owning item (non-blocking comm if remote)
        
        /// If item does not exist it is made with the default constructor.  
        /// 
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        /// 
        /// Returns a future result (Future<void> may be ignored).
        template <typename memfunT, typename arg1T>
        Future< MEMFUN_RETURNT(memfunT) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1) {
            check_initialized();
            RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT)) (implT::*itemfun)(const keyT&, memfunT, const arg1T&) = &implT:: template itemfun<memfunT,arg1T>;
            return p->task(owner(key), itemfun, key, memfun, arg1);
        }
        
        /// Adds task "resultT memfun(arg1T,arg2T)" in process owning item (non-blocking comm if remote)
        
        /// If item does not exist it is made with the default constructor.  
        /// 
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        /// 
        /// Returns a future result (Future<void> may be ignored).
        template <typename memfunT, typename arg1T, typename arg2T>
        Future< MEMFUN_RETURNT(memfunT) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2) {
            check_initialized();
            RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT)) (implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&) = &implT:: template itemfun<memfunT,arg1T,arg2T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2);
        }
        
        /// Adds task "resultT memfun(arg1T,arg2T,arg3T)" in process owning item (non-blocking comm if remote)
        
        /// If item does not exist it is made with the default constructor.  
        /// 
        /// Future arguments for local tasks can generate dependencies, but for remote
        /// tasks all futures must be ready.
        /// 
        /// Returns a future result (Future<void> may be ignored).
        template <typename memfunT, typename arg1T, typename arg2T, typename arg3T>
        Future< MEMFUN_RETURNT(memfunT) >
        task(const keyT& key, memfunT memfun, const arg1T& arg1, const arg2T& arg2, const arg3T& arg3) {
            check_initialized();
            RETURN_WRAPPERT(MEMFUN_RETURNT(memfunT)) (implT::*itemfun)(const keyT&, memfunT, const arg1T&, const arg2T&, const arg3T&) = &implT:: template itemfun<memfunT,arg1T,arg2T,arg3T>;
            return p->task(owner(key), itemfun, key, memfun, arg1, arg2, arg3);
        }
        
        
        /// (de)Serialize --- !! ONLY for purpose of interprocess communication
        
        /// This just writes/reads the unique id to/from the archive.  If you want
        /// to serialize the actual contents you'll have to write your own
        /// specialization of archive::ArchiveLoadImpl and ditto for store
        /// for that specific type of archive, or recast this in that form.
        template <typename Archive>
        void serialize(const Archive& ar) {
            if (Archive::is_output_archive) {
                check_initialized();
                ar & static_cast<WorldObject<implT>*>(p.get());
            }
            else {
                WorldObject<implT>* ptr;
                ar & ptr;
                MADNESS_ASSERT(ptr);
                p = SharedPtr<implT>(static_cast<implT*>(ptr),false,false);  // use_count will be 0, which is good
            }
        }
        
        /// Destructor passes ownership of implementation to world for deferred cleanup
        virtual ~WorldContainer() {};
        
    };
}

#endif
