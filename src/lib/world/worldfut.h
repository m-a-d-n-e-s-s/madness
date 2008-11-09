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

  
  $Id$
*/

  
#ifndef WORLD_FUT_H
#define WORLD_FUT_H

/// \file worldfut.h
/// \brief Implements Future

#include <vector>
#include <world/nodefaults.h>
#include <world/worlddep.h>
#include <world/array.h>
#include <world/sharedptr.h>
#include <world/worldref.h>
#include <world/typestuff.h>

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

    template <typename T>
    std::ostream& operator<<(std::ostream& out, const Future<T>& f);


    /// Implements the functionality of Futures
    template <typename T>
    class FutureImpl : private Spinlock {
        friend class Future<T>;
        friend std::ostream& operator<< <T>(std::ostream& out, const Future<T>& f);

    private:
        static const int MAXCALLBACKS = 4;
        typedef Stack<CallbackInterface*,MAXCALLBACKS> callbackT;
        typedef Stack<SharedPtr< FutureImpl<T> >,MAXCALLBACKS> assignmentT;
        volatile callbackT callbacks;
        volatile assignmentT assignments;
        volatile bool assigned;
        World * world;
        RemoteReference< FutureImpl<T> > remote_ref;
        volatile T t;
        
        /// Private: AM handler for remote set operations
        static void set_handler(const AmArg& arg) {
            RemoteReference< FutureImpl<T> > ref;
            T t;
            arg & ref & t;
            FutureImpl<T>* f = ref.get();
            f->set(t);
            ref.dec();                    // Releases reference
        };
        

        /// Private:  invoked locally by set routine after assignment
        inline void set_assigned() {
            // ASSUME THAT WE HAVE THE MUTEX WHILE IN HERE (so that
            // destructor is not invoked by another thread as a result
            // of our actions)
            MADNESS_ASSERT(!assigned);
            assigned = true;

            callbackT& cb = const_cast<callbackT&>(callbacks);
            assignmentT& as = const_cast<assignmentT&>(assignments);
            while (as.size()) {
                SharedPtr< FutureImpl<T> >& p = as.pop();
                MADNESS_ASSERT(p);
                p->set(const_cast<T&>(t));
                p = SharedPtr< FutureImpl<T> >();
            };
            while (cb.size()) {
                CallbackInterface* p = cb.pop();
                MADNESS_ASSERT(p);
                p->notify();
            }
        };

        inline void add_to_assignments(const SharedPtr< FutureImpl<T> >& f) {
            // ASSUME lock is already acquired by Future<T>::set()
            assignmentT& nvas = const_cast<assignmentT&>(assignments);
            nvas.push(f);
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
        {
            //print("FUTCON(a)",(void*) this);
        };
        
        
        /// Local assigned value
        FutureImpl(const T& t) 
            : callbacks()
            , assignments()
            , assigned(false)
            , world(0)
            , remote_ref()
            , t(t)
        {
            //print("FUTCON(b)",(void*) this);
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
        {
            //print("FUTCON(c)",(void*) this);
        };


        /// Returns true if the value has been assigned
        inline bool probe() const {
            return assigned;
        };
        

        /// Registers a function to be invoked when future is assigned

        /// Callbacks are invoked in the order registered.  If the
        /// future is already assigned the callback is immediately
        /// invoked.
        inline void register_callback(CallbackInterface* callback) {
            ScopedMutex<Spinlock> fred(this);
            if (assigned) callback->notify();
            else const_cast<callbackT&>(callbacks).push(callback);
        };

         
        /// Sets the value of the future (assignment)
        void set(const T& value) {
            ScopedMutex<Spinlock> fred(this);
            if (world) { 
                if (remote_ref.owner() == world->rank()) {
                    remote_ref.get()->set(value); 
                    remote_ref.dec(); // Releases reference
                }
                else {
                    world->am.send(remote_ref.owner(), 
                                   FutureImpl<T>::set_handler, 
                                   new_am_arg(remote_ref, value));
                }
            }
            else {
                const_cast<T&>(t) = value;
            }
            set_assigned();
        };
        

        /// Gets/forces the value, waiting if necessary (error if not local)
        T& get() {
            MADNESS_ASSERT(!world);  // Only for local futures
            World::await(bind_nullary_mem_fun(this,&FutureImpl<T>::probe));
            return *const_cast<T*>(&t);
        };

        /// Gets/forces the value, waiting if necessary (error if not local)
        const T& get() const {
            MADNESS_ASSERT(!world);  // Only for local futures
            World::await(bind_nullary_mem_fun(this,&FutureImpl<T>::probe));
            return *const_cast<const T*>(&t);
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
            throw "IS THIS WORKING? maybe now we have the mutex";
//            ScopedMutex<Spinlock> fred(this);
//             MADNESS_ASSERT(!world); // was return false;
//             MADNESS_ASSERT(!assigned || f->assigned);
//             if (f->world) {
//                 world = f->world;
//                 remote_ref = f->remote_ref;
//                 f->world = 0;
//             }
//             while(f->callbacks.size()) callbacks.push(f->callbacks.pop());
//             while(f->assignments.size()) assignments.push(f->assignments.pop());
            return true;
        };

        virtual ~FutureImpl() {
            ScopedMutex<Spinlock> fred(this);
            //print("FUTDEL",(void*) this);
            if (!assigned && world) {
                print("Future: unassigned remote future being destroyed?");
                //remote_ref.dec();
                //abort();
            }
            if (const_cast<callbackT&>(callbacks).size()) {
                print("Future: uninvoked callbacks being destroyed?");
                //abort();
            }
            if (const_cast<assignmentT&>(assignments).size()) {
                print("Future: uninvoked assignment being destroyed?");
                //abort();
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
        friend std::ostream& operator<< <T>(std::ostream& out, const Future<T>& f);

    private:
        mutable SharedPtr< FutureImpl<T> > f;
        T value;
        const bool value_set;
        const bool is_the_default_initializer;

        class dddd{};
        explicit Future(const dddd&) 
            : f()
            , value()
            , value_set(false)
            , is_the_default_initializer(true)
        {}

    public:
        typedef RemoteReference< FutureImpl<T> > remote_refT;

        /// Makes an unassigned future
        Future() 
            : f(new FutureImpl<T>())
            , value()
            , value_set(false)
            , is_the_default_initializer(false)
        {}

        /// Makes an assigned future
        explicit Future(const T& t) 
            : f()
            , value(t)
            , value_set(true)
            , is_the_default_initializer(false)
        {}


        /// Makes a future wrapping a remote reference
        explicit Future(const remote_refT& remote_ref)
            : f(new FutureImpl<T>(remote_ref))
            , value()
            , value_set(false)
            , is_the_default_initializer(false)
        {}


        /// Copy constructor is shallow
        Future(const Future<T>& other) 
            : f(0)
            , value(other.value)
            , value_set(other.value_set)
            , is_the_default_initializer(false)
        {
            if (other.is_the_default_initializer) {
                f = SharedPtr< FutureImpl<T> >(new FutureImpl<T>());
            }
            else if (!other.value_set) {
                f = other.f;
            }
        }


        /// See Gotchas on the documentation mainpage about why this exists and how to use it.
        static const Future<T> default_initializer() {
            return Future<T>(dddd());
        }


        /// Assignment future = future makes a shallow copy just like copy constructor
        Future<T>& operator=(const Future<T>& other) {
            if (this != &other) {
                MADNESS_ASSERT(!value_set);
                f = other.f;
            }
            return *this;
        }

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
            //
            // THIS NEEDS THOROUGH CHECKING FOR RACE NASTIES
            //

            // In a correct program a future is only assigned once so
            // it is safe-ish to assume that 'this' is single threaded
            // ... however other may be having stuff done to it
            // concurrently by other threads.

            if (this != &other) { 
                MADNESS_ASSERT(!probe());
                if (other.probe()) {
                    set(other.get());     // The easy case
                }
                else {
                    other.f->lock();
                    // Both are initialized but not assigned.  Used to
                    // just put this on the assignment list of other,
                    // but now first check to see if the internal
                    // state of either can be harmlessly replaced.
                    // This requires that one of this or other (but
                    // not both) have a reference of one.  The main
                    // optimization enabled by this is shortcircuiting
                    // communication when forwarding futures that are
                    // actually remote references.  However, it also
                    // avoids most uses of the assignment list.
                    //
                    if (f == other.f) {
                        print("future.set(future): you are setting this with this?");
                    }
                    else if (f.use_count()==1 && other.f.use_count()==1 && other.f->is_local()) {
                        other.f->replace_with(f);
                        f = other.f;
                        //print("future.set(future): replaced other with this");
                    }
                    else if (f.use_count()==1 && other.f.use_count()==1 && f->is_local()) {
                        f->replace_with(other.f);
                        const_cast<Future<T>&>(other).f = f;
                        //print("future.set(future): replaced this with other");
                    }
                    else {
                        //print("future.set(future): adding to assignment list");
                        other.f->add_to_assignments(f);
                    }
                    other.f->unlock();
                    return;
                }
            }
        }


        /// Assigns the value ... it can only be set ONCE.
        inline void set(const T& value) {
            MADNESS_ASSERT(!value_set);
            f->set(value);
        }


        /// Gets the value, waiting if necessary (error if not a local future)
        inline T& get() {
            if (value_set) {
                return value;
            }
            else {
                return f->get();
            }
        }


        /// Gets the value, waiting if necessary (error if not a local future)
        inline const T& get() const {
            if (value_set) {
                return value;
            }
            else {
                return f->get();
            }
        }


        /// Returns true if the future has been assigned
        inline bool probe () const {
            if (f) return f->probe();
            else return value_set;
        }


        /// Same as get()
        inline operator T () {
            return get();
        }

        /// Same as get() const
        inline operator const T () const {
            return get();
        }


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
        /// a reference to the local future.  Therefore, the local
        /// future will not be notified when the result is set
        /// (i.e., the communication is short circuited).
        inline remote_refT remote_ref(World& world) const {
            MADNESS_ASSERT(!probe());
            if (f->world) 
                return f->remote_ref;
            else 
                return RemoteReference< FutureImpl<T> >(world, f);
        }


        inline bool is_local() const {
            return ((!f) || f->is_local() || probe());
        }

        inline bool is_remote() const {
            return !is_local();
        }


        /// Registers an object to be called when future is assigned

        /// Callbacks are invoked in the order registered.  If the
        /// future is already assigned the callback is immediately
        /// invoked. 
        inline void register_callback(CallbackInterface* callback) {
            if (probe()) callback->notify();
            else {
                f->register_callback(callback);
            }
        }
    };


    /// A future of a future is forbidden (by private constructor)
    template <typename T> class Future< Future<T> >{
        Future(){}
    };


    /// Specialization of FutureImpl<void> for internal convenience ... does nothing useful!
    template <> class FutureImpl<void> {};

    /// Specialization of Future<void> for internal convenience ... does nothing useful!
    template <> class Future<void> {
    public:
        RemoteReference< FutureImpl<void> > remote_ref(World& world) const {
            return RemoteReference< FutureImpl<void> >();
        }

        Future(){}

        Future(const RemoteReference< FutureImpl<void> >& ref) {}

        Future(const Future<Void>& f) {}
        
        inline void set(const Future<void>& f) {}

        inline Future<void>& operator=(const Future<void>& f) {
            return *this;
        }

        inline void set() {}

        template <class Archive>
        void serialize(const Archive& ar) {}

        virtual ~Future(){}
    };

    /// Specialization of FutureImpl<Void> for internal convenience ... does nothing useful!
    template <> class FutureImpl<Void> {};

    /// Specialization of Future<Void> for internal convenience ... does nothing useful!
    template <> class Future<Void> {
    public:
        RemoteReference< FutureImpl<Void> > remote_ref(World& world) const {
            return RemoteReference< FutureImpl<Void> >();
        }

        Future(){}

        Future(const RemoteReference< FutureImpl<Void> >& ref) {}

        Future(const Future<void>& f) {}
        
        inline void set(const Future<Void>& f) {}

        inline Future<Void>& operator=(const Future<Void>& f) {
            return *this;
        }

        inline void set(const Void& f) {}

        template <class Archive>
        void serialize(const Archive& ar) {}

        virtual ~Future(){}
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
        Future() : v() {}

        Future(const vectorT& v) : v(v) { // Register dependencies as created
            for (int i=0; i<(int)v.size(); i++) {
                if (!this->v[i].probe()) {
                    inc();
                    this->v[i].register_callback(this);
                }
            }
        }
        vectorT& get() {return v;}
        const vectorT& get() const {return v;}
        operator vectorT& () {return get();}
        operator const vectorT& () const {return get();}
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
        else if (f.is_remote()) out << f.f->remote_ref;
        else if (f.f) out << "<unassigned refcnt=" << f.f.use_count() << ">";
        else out << "<unassigned>";
        return out;
    }

#endif

}


#endif
