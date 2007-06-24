/*
  This file is part of MADNESS.
  
  Copyright (C) <2007> <Oak Ridge National Laboratory>
  
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
  
  For more information please contact:

  Robert J. Harrison
  Oak Ridge National Laboratory
  One Bethel Valley Road
  P.O. Box 2008, MS-6367

  email: harrisonrj@ornl.gov 
  tel:   865-241-3937
  fax:   865-572-0680

  $LastChangedDate$
  $Rev$
*/

  
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
        static const int MAXCALLBACKS = 4;
        //std::vector<CallbackInterface*> callbacks; 
        //std::vector< SharedPtr< FutureImpl<T> > > assignments;
        Stack<CallbackInterface*,MAXCALLBACKS> callbacks;
        Stack<SharedPtr< FutureImpl<T> >,MAXCALLBACKS> assignments;
        volatile bool assigned;
        World * world;
        RemoteReference< FutureImpl<T> > remote_ref;
        T t;
        
        /// Private: AM handler for remote set operations
        static void set_handler(World& world, ProcessID src, void *buf, size_t nbyte) {
            LongAmArg* arg = (LongAmArg*) buf;
            BufferInputArchive ar(arg->buf,nbyte);
            RemoteReference< FutureImpl<T> > ref;
            ar & ref;
            FutureImpl<T>* f = ref.get();
            ar & f->t;                    // Assigns value
            f->set_assigned();
            ref.dec();                    // Releases reference
        };
        

        /// Private:  invoked locally by set routine after assignment
        inline void set_assigned() {
            MADNESS_ASSERT(!assigned);
            assigned = true;
//             for (int i=0; i<(int)callbacks.size(); i++)  
//                 callbacks[i]->notify();
//             callbacks.clear();
//             for (int i=0; i<(int)assignments.size(); i++) 
//                 assignments[i]->set(t);
//             assignments.clear();
            while (callbacks.size()) callbacks.pop()->notify();
            while (assignments.size()) {
                SharedPtr< FutureImpl<T> >& p = assignments.pop();
                p->set(t);
                p = SharedPtr< FutureImpl<T> >();
            };
        };

        inline void add_to_assignments(const SharedPtr< FutureImpl<T> >& f) {
            //assignments.push_back(f);
            assignments.push(f);
        };
            
        
    public:
        
        /// Local unassigned value
        FutureImpl() 
            : callbacks()
            , assignments()
            , assigned(false)
            , world(0)
            , remote_ref()
            , t()
        {};
        
        
        /// Local assigned value
        FutureImpl(const T& t) 
            : callbacks()
            , assignments()
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
            , assignments()
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
            else callbacks.push(callback);
            //else callbacks.push_back(callback);
        };

         
        /// Sets the value of the future (assignment)
        void set(const T& value) {
            if (world) { 
                if (remote_ref.owner() == world->rank()) {
                    remote_ref.get()->set(value); 
                    remote_ref.dec(); // Releases reference
                }
                else {
                    LongAmArg* arg = new LongAmArg;
                    size_t nbyte = arg->stuff(remote_ref, value);
                    world->am.send_long_managed(remote_ref.owner(), 
                                                FutureImpl<T>::set_handler, 
                                                arg, nbyte);
                }
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

        bool is_local() const {
            return world == 0;
        };

        bool replace_with(FutureImpl<T>* f) {
            if (world) return false;
            MADNESS_ASSERT(!assigned || f->assigned);
            if (f->world) {
                world = f->world;
                remote_ref = f->remote_ref;
                f->world = 0;
            }
            while(f->callbacks.size()) callbacks.push(f->callbacks.pop());
            while(f->assignments.size()) assignments.push(f->assignments.pop());
            return true;
        };

        virtual ~FutureImpl() {
            if (!assigned && world) {
                print("Future: unassigned remote future being destroyed?");
                remote_ref.dec();
            }
            if (callbacks.size()) {
                print("Future: uninvoked callbacks being destroyed?");
            }
            if (assignments.size()) {
                print("Future: uninvoked assignment being destroyed?");
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
        typedef RemoteReference< FutureImpl<T> > remote_refT;

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
         explicit Future(const remote_refT& remote_ref)
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


//         /// Assignment future = value is the same as future.set(value)
//         inline Future<T>& operator=(const T& value) {
//             this->set(value);
//             return *this;
//         };


        /// Assignment future = future makes a shallow copy

        /// The future A is replaced by a shallow copy of B.  If B is
        /// assigned, A=B causes A to have the same value as B, though
        /// it does not necessarily refer to the same actual data.  If
        /// B is not assigned, A=B makes both A and B refer to the
        /// same underlying future so that both variables will "see"
        /// the eventual assigned value.
        ///
        /// Note that this is \em not the same as A.set(B) which will throw
        /// if A is already assigned and otherwise ensures that A will
        /// take on the (future) value of B.
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
        /// Future is a single assignment variable.  We don't yet
        /// track multiple assignments from unassigned futures.
        ///
        /// If B is already assigned, this is the same as A.set(B.get())
        /// which sets A to the value of B.
        ///
        /// If B has not yet been assigned, the behavior is to ensure
        /// that when B is assigned that both A and B will be assigned
        /// and have the same value (though they may/may not refer to
        /// the same underlying copy of the data and indeed may even
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
                        // Used to just put this on the assignment
                        // list of other, but now first check to see
                        // if the internal state of either can be
                        // harmlessly replaced.  This requires that
                        // one of this or other (but not both) have a
                        // reference of one.  The main optimization
                        // enabled by this is shortcircuiting
                        // communication when forwarding futures that
                        // are actually remote references.  However,
                        // it also avoids most uses of the assignment
                        // list. 
                        if (f == other.f) {
                            print("future.set(future): you are setting this with this?");
                            return;
                        }
                        if (f.use_count()==1 && other.f->is_local()) {
                            other.f->replace_with(f);
                            f = other.f;
                            //print("future.set(future): replaced other with this");
                            return;
                        }
                        if (other.f.use_count()==1 && f->is_local()) {
                            f->replace_with(other.f);
                            const_cast<Future<T>&>(other).f = f;
                            print("future.set(future): replaced this with other");
                            return;
                        }
                        print("future.set(future): adding to assignment list");
                        other.f->add_to_assignments(f);
                        return;
                    }
                    else {
                        // This is initialized but not assigned, and
                        // other is not initialized.  Make other refer
                        // to the same future as this.
                        //print("future.set(future): this initialized and other not initialized");
                        const_cast<Future<T>&>(other).f = f;
                        return;
                    }
                }
                else { 
                    // Neither initialized ... make this refer to other
                    other.initialize();
                    f = other.f;
                    return;
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


        /// Returns a structure used to pass references to another process.

        /// This is used for passing pointers/references to another
        /// process.  To make remote references completely safe, the
        /// RemoteReference increments the internal reference count of
        /// the Future.  The counter is decremented by either
        /// assigning to the remote Future or its destructor if it is
        /// never assigned.  The remote Future is ONLY useful for
        /// setting the future.  It will NOT be notified if the value
        /// is set elsewhere.
        ///
        /// If this is already a reference to a remote future, the
        /// actual remote reference is returned ... i.e., \em not a
        /// a reference to the local future.
        inline remote_refT remote_ref(World& world) const {
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


    /// A future of a future is forbidden (by private constructor)
    template <typename T> class Future< Future<T> >{
        Future(){};
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

        Future(const Future<Void>& f) {};
        
        inline void set(const Future<void>& f) {};

        inline Future<void>& operator=(const Future<void>& f) {
            return *this;
        };

        inline void set() {};

        template <class Archive>
        void serialize(const Archive& ar) {}

        virtual ~Future(){};
    };

    /// Specialization of FutureImpl<Void> for internal convenience ... does nothing useful!
    template <> class FutureImpl<Void> {};

    /// Specialization of Future<Void> for internal convenience ... does nothing useful!
    template <> class Future<Void> {
    public:
        RemoteReference< FutureImpl<Void> > remote_ref(World& world) const {
            return RemoteReference< FutureImpl<Void> >();
        };

        Future(){};

        Future(const RemoteReference< FutureImpl<Void> >& ref) {};

        Future(const Future<void>& f) {};
        
        inline void set(const Future<Void>& f) {};

        inline Future<Void>& operator=(const Future<Void>& f) {
            return *this;
        };

        inline void set(const Void& f) {};

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
        typedef typename std::vector< Future<T> > vectorT;
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

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<void>& f);

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<Void>& f);

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


        /// Serialize an assigned future
        template <class Archive>
        struct ArchiveStoreImpl< Archive, Future<void> > {
            static inline void store(const Archive& ar, const Future<void>& f) {
            }
        };
        
        
        /// Deserialize a future into an unassigned future
        template <class Archive>
        struct ArchiveLoadImpl< Archive, Future<void> > {
            static inline void load(const Archive& ar, Future<void>& f) {
            }
        };

        /// Serialize an assigned future
        template <class Archive>
        struct ArchiveStoreImpl< Archive, Future<Void> > {
            static inline void store(const Archive& ar, const Future<Void>& f) {
            }
        };
        
        
        /// Deserialize a future into an unassigned future
        template <class Archive>
        struct ArchiveLoadImpl< Archive, Future<Void> > {
            static inline void load(const Archive& ar, Future<Void>& f) {
            }
        };
    }

    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f) ;

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<void>& f) ;

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<Void>& f) ;

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
