#ifndef WORLD_FUT_H
#define WORLD_FUT_H

/// \file worldfut.h
/// \brief Implements Future


namespace madness {


    template <typename T> class Future;


    /// Boost-type-trait-like testing of if a type is a future
    template <typename T>
    struct is_future {
        static const bool value = false;
    };

    /// Boost-type-trait-like testing of if a type is a future
    template <typename T>
    struct is_future< Future<T> > {
        static const bool value = true;
    };

    /// Boost-type-trait-like mapping of Future<T> to T
    template <typename T>
    struct remove_future {
        typedef T type;
    };

    /// Boost-type-trait-like mapping of Future<T> to T
    template <typename T>
    struct remove_future< Future<T> > {
        typedef T type;
    };

#define REMFUTURE(T) typename remove_future< T >::type    


    /// Implements the functionality of Futures
    template <typename T>
    class FutureImpl : NO_DEFAULTS {
        friend class Future<T>;
    private:
        std::vector<CallbackInterface*> callbacks; 
        std::vector< SharedPtr< FutureImpl<T> > > assignment_list;
        volatile bool assigned;
        World * const world;
        RemoteReference< FutureImpl<T> > remote_ref;
        T t;
        
        /// Private: AM handler for remote set operations
        static void set_handler(World& world, ProcessID src, void *buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg*) buf;
            BufferInputArchive ar(arg->buf,nbyte);
            RemoteReference< FutureImpl<T> > ptr;
            ar & ptr;
            FutureImpl<T>* f = ptr.get();
            ar & f->t;                    // Assigns value
            f->set_assigned();
            ptr.dec();                    // Releases reference
        };
        

        /// Private:  invoked locally by set routine after assignment
        inline void set_assigned() {
            MADNESS_ASSERT(!assigned);
            assigned = true;
            for (int i=0; i<(int)callbacks.size(); i++)  
                callbacks[i]->notify();
            callbacks.clear();
            for (int i=0; i<(int)assignment_list.size(); i++) 
                assignment_list[i]->set(t);
            assignment_list.clear();
        };

        inline void add_to_assignment_list(const SharedPtr< FutureImpl<T> >& f) {
            assignment_list.push_back(f);
        };
            
        
    public:
        
        /// Local unassigned value
        FutureImpl() 
            : callbacks()
            , assignment_list()
            , assigned(false)
            , world(0)
            , remote_ref()
            , t()
        {};
        
        
        /// Local assigned value
        FutureImpl(const T& t) 
            : callbacks()
            , assignment_list()
            , assigned(false)
            , world(0)
            , remote_ref()
            , t(t)
        {
            set_assigned();
        };
        

        /// Wrapper for a remote future
        FutureImpl(const RemoteReference< FutureImpl<T> >& remote_ref)
            : callbacks()
            , assignment_list()
            , assigned(false)
            , world(World::world_from_id(remote_ref.world_id()))
            , remote_ref(remote_ref)
            , t()
        {};


        /// Returns true if the value has been assigned
        inline bool probe() const {
            return assigned;
        };
        

        /// Registers a function to be invoked when future is assigned

        /// Callbacks are invoked in the order registered.  If the
        /// future is already assigned the callback is immediately
        /// invoked.
        inline void register_callback(CallbackInterface* callback) {
            if (assigned) callback->notify();
            else callbacks.push_back(callback);
        };

         
        /// Sets the value of the future (assignment)
        void set(const T& value) {
            if (world) { 
                LongAmArg* arg = new LongAmArg;
                size_t nbyte = arg->stuff(remote_ref, value);
                world->am.send_long_managed(remote_ref.owner(), 
                                            FutureImpl<T>::set_handler, 
                                            arg, nbyte);
            }
            else {
                t = value;
            }
            set_assigned();
        };
        

        /// Gets/forces the value, waiting if necessary (error if not local)
        T& get() {
            MADNESS_ASSERT(!world);  // Only for local futures
            World::await(bind_nullary_mem_fun(this,&FutureImpl<T>::probe));
            return t;
        };

        /// Gets/forces the value, waiting if necessary (error if not local)
        const T& get() const {
            MADNESS_ASSERT(!world);  // Only for local futures
            World::await(bind_nullary_mem_fun(this,&FutureImpl<T>::probe));
            return t;
        };

        T& operator->() {
            return get();
        };

        const T& operator->() const {
            return get();
        };

        virtual ~FutureImpl() {
            if (!assigned && world) {
                print("Future: unassigned remote future being destroyed?");
                remote_ref.dec();
            }
            if (callbacks.size()) {
                print("Future: uninvoked callbacks being destroyed?");
            }
        };
    };
                        
    /// A future is a possibly yet unevaluated value

    /// Uses delegation to FutureImpl to provide desired
    /// copy/assignment semantics as well as safe reference counting
    /// for remote futures.
    ///
    /// Since we are using Futures a lot to store local values coming
    /// from containers and inside task wrappers for messages, we
    /// included in this class a value.  If a future is assigned
    /// before a copy/remote-reference is taken, the shared ptr is
    /// never made.  The point of this is to eliminate the two mallocs
    /// that must be peformed for every new shared_ptr.
    template <typename T> 
    class Future {
    private:
        mutable SharedPtr< FutureImpl<T> > f;
        mutable T value;
        mutable bool value_set;
        bool is_the_default_initializer;

        /// Private:  ensures internal implementation is initialized before use
        inline void initialize(const T* value = 0) const {
            if (!f) {
                if (value) {
                    f = SharedPtr< FutureImpl<T> >(new FutureImpl<T>(*value));
                }
                else {
                    f = SharedPtr< FutureImpl<T> >(new FutureImpl<T>());
                }
            }
        };

    public:

        /// Makes an unassigned local future
        Future() 
            : f()
            , value()
            , value_set(false)
            , is_the_default_initializer(false)
        {};

        /// Makes an assigned local future
        /// get rid of these?
        explicit Future(const T& t) 
            : f()
            , value(t)
            , value_set(true)
            , is_the_default_initializer(false)
        {};


        /// Makes a future wrapping a remote reference
        explicit Future(const RemoteReference< FutureImpl<T> >& remote_ref)
            : f(new FutureImpl<T>(remote_ref))
            , value()
            , value_set(false)
            , is_the_default_initializer(false)
        {};


        /// Copy constructor is shallow

        /// Trying to copy an unconstructed future forces its
        /// construction unless it is already locally assigned.  The
        /// exception is if other is default initializer which makes
        /// the copy constructor behave like the default constructor
        /// (see Gotchas on the main documentation page).
        Future(const Future<T>& other) 
            : f()
            , value()
            , value_set(other.value_set)
            , is_the_default_initializer(false)
        {
            if (!other.is_the_default_initializer) {
                if (other.value_set) {
                    value = other.value;
                }
                else {
                    other.initialize();
                    f = other.f;
                }
            }
        };


        /// See Gotchas on the documentation mainpage about why this exists and how to use it.
        static const Future<T>& default_initializer() {
            static Future<T> result;
            result.is_the_default_initializer = true;
            return result;
        };


        /// Assignment future = value is the same as future.set(value)
        inline Future<T>& operator=(const T& value) {
            this->set(value);
            return *this;
        };


        /// Assignment future = future makes a shallow copy

        /// The future A is replaced by a shallow copy of B.  If B is not
        /// yet assigned, both A and B will refer to the same underlying
        /// implementation.  Note that this is NOT the same as set.  
        /// A.set(B) assigns the (future) value of B to A whereas A=B
        /// destroys the existing future A and replaces it either with
        /// the assigned value of B or a reference to its future value.
        Future<T>& operator=(const Future<T>& other) {
            if (this != &other) {
                value_set = other.value_set;
                if (value_set) {
                    value = other.value;
                    f = SharedPtr< FutureImpl<T> >(0);
                }
                else {
                    other.initialize();
                    f = other.f;
                }
            }
            return *this;
        };


        /// A.set(B) where A & B are futures ensures A has/will-have the same value as B. 

        /// An exception is thrown if A is already assigned since a
        /// Future is a single assignment variable.
        ///
        /// If B is already assigned, this is the same as A.set(B.get()).
        ///
        /// If B has not yet been assigned, the behavior is to ensure
        /// that when B is assigned that both A and B will be assigned
        /// and have the same value (though they may/may not refer to 
        /// the same underlying copy of the value and indeed may even
        /// be in different processes).
        void set(const Future<T>& other) {
            if (this != &other) { 
                MADNESS_ASSERT(!probe());
                if (other.probe()) {
                    set(other.get());     // The easy case
                }
                else if (f) {
                    if (other.f) {
                        // Both are initialized but not assigned.
                        // Store a reference to this inside other
                        if (f != other.f) other.f->add_to_assignment_list(f);
                    }
                    else {
                        // This is initialized but not assigned, and
                        // other is not initialized.  Make other refer
                        // to the same future as this.
                        const_cast<Future<T>&>(other).f = f;
                    }
                }
                else { 
                    // This is not initialized so just make it refer to other
                    other.initialize();
                    f = other.f;
                }
            }
        };


        /// Assigns the value ... it can only be set ONCE.
        inline void set(const T& value) {
            if (f) {
                f->set(value);
            }
            else {
                MADNESS_ASSERT(!value_set);
                this->value = value;
                value_set = true;
            }
        };


        /// Gets the value, waiting if necessary (error if not a local future)
        inline T& get() {
            if (value_set) {
                return value;
            }
            else {
                initialize();
                return f->get();
            }
        };


        /// Gets the value, waiting if necessary (error if not a local future)
        inline const T& get() const {
            if (value_set) {
                return value;
            }
            else {
                initialize();
                return f->get();
            }
        };


        /// Returns true if the future has been assigned
        inline bool probe () const {
            if (f) return f->probe();
            else return value_set;
        };


        /// Same as get()
        inline operator T () {
            return get();
        };

        /// Same as get() const
        inline operator const T () const {
            return get();
        };


        /// Private: Returns a structure used to pass references to another process.

        /// This is used for passing pointers/references to another
        /// process.  To make remote references completely safe, the
        /// RemoteReference increments the internal reference count of
        /// the Future.  The counter is decremented by either
        /// assigning to the remote Future or its destructor if it is
        /// never assigned.  The remote Future is ONLY useful for
        /// setting the future.  It will NOT be notified if the value
        /// is set elsewhere.
        inline RemoteReference< FutureImpl<T> > remote_ref(World& world) const {
            MADNESS_ASSERT(!probe());
            initialize();
            if (f->world) 
                return f->remote_ref;
            else 
                return RemoteReference< FutureImpl<T> >(world, f);
        };


        /// Registers an object to be called when future is assigned

        /// Callbacks are invoked in the order registered.  If the
        /// future is already assigned the callback is immediately
        /// invoked. 
        inline void register_callback(CallbackInterface* callback) {
            if (probe()) callback->notify();
            else {
                initialize();
                f->register_callback(callback);
            }
        };
    };


    /// Specialization of FutureImpl<void> for internal convenience ... does nothing useful!
    template <> class FutureImpl<void> {};


    /// Specialization of Future<void> for internal convenience ... does nothing useful!
    template <> class Future<void> {
        
    public:
        RemoteReference< FutureImpl<void> > remote_ref(World& world) const {
            return RemoteReference< FutureImpl<void> >();
        };

        Future(){};

        Future(const RemoteReference< FutureImpl<void> >& ref) {};

        
        inline void set(const Future<void>& f) {};

        inline Future<void>& operator=(const Future<void>& f) {
            return *this;
        };

        inline void set() {};

        template <class Archive>
        void serialize(const Archive& ar) {}

        virtual ~Future(){};
    };

    /// Specialization of Future for vector of Futures

    /// Enables passing a vector of futures into a task and having
    /// the dependencies correctly tracked.  Does not directly
    /// support most operations that other futures do ... that is
    /// the responsiblility of the individual futures in the vector.
    template <typename T>
    class Future< std::vector< Future<T> > > : public DependencyInterface {
    private:
        typedef std::vector< Future<T> > vectorT;
        vectorT v;

    public:
        Future() : v() {};
        Future(const vectorT& v) : v(v) { // Register dependencies as created
            for (int i=0; i<(int)v.size(); i++) {
                if (!this->v[i].probe()) {
                    inc();
                    this->v[i].register_callback(this);
                }
            }
        };
        vectorT& get() {return v;};
        const vectorT& get() const {return v;};
        operator vectorT& () {return get();};
        operator const vectorT& () const {return get();};
    };
        


    /// Probes a future for readiness, other types are always ready
    template <typename T>
    ENABLE_IF(is_future<T>,bool) future_probe(const T& t) {return t.probe();}

    /// Probes a future for readiness, other types are always ready
    template <typename T>
    DISABLE_IF(is_future<T>,bool) future_probe(const T& t) {return true;}
        

    /// Friendly I/O to streams for futures
    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f);

    /// Factory for vectors of futures (see section Gotchas on the mainpage)
    template <typename T>
    std::vector< Future<T> > future_vector_factory(std::size_t n) {
        return std::vector< Future<T> >(n, Future<T>::default_initializer());
    }



    namespace archive {
        /// Serialize an assigned future
        template <class Archive, typename T>
        struct ArchiveStoreImpl< Archive, Future<T> > {
            static inline void store(const Archive& ar, const Future<T>& f) {
                MAD_ARCHIVE_DEBUG(std::cout << "serializing future" << std::endl);
                MADNESS_ASSERT(f.probe());
                ar & f.get();
            }
        };
        
        
        /// Deserialize a future into an unassigned future
        template <class Archive, typename T>
        struct ArchiveLoadImpl< Archive, Future<T> > {
            static inline void load(const Archive& ar, Future<T>& f) {
                MAD_ARCHIVE_DEBUG(std::cout << "deserializing future" << std::endl);
                MADNESS_ASSERT(!f.probe());
                T value;
                ar & value;
                f.set(value);
            }
        };
    }

#ifdef WORLD_INSTANTIATE_STATIC_TEMPLATES
    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f) {
        if (f.probe()) out << f.get();
        else out << "<unassigned>";
        return out;
    }
#endif

}


#endif
