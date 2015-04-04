/*
  This file is part of MADNESS.

  Copyright (C) 2007,2010 Oak Ridge National Laboratory

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

/**
 \file worldfwd.h
 \brief Implements World.
 \addtogroup world



 \addtogroup world
@{

*/
  
#ifndef MADNESS_WORLD_WORLDFWD_H__INCLUDED
#define MADNESS_WORLD_WORLDFWD_H__INCLUDED

#include <madness/madness_config.h>

// #ifdef SEEK_SET
// #undef SEEK_SET
// #endif
// #ifdef SEEK_CUR
// #undef SEEK_CUR
// #endif
// #ifdef SEEK_END
// #undef SEEK_END
// #endif

// Standerd C++ header files needed by world.h
#include <iostream>
#include <list>
#include <utility>
#include <cstddef>

#ifdef HAVE_RANDOM
#include <stdlib.h>
#endif

// Madness world header files needed by world
#include <madness/world/worldmpi.h>
#include <madness/world/worldhashmap.h>
#include <madness/world/worldprofile.h>
#include <madness/world/worldthread.h>
#include <madness/world/uniqueid.h>
#include <madness/world/nodefaults.h>

namespace madness {

    class World;
    class WorldTaskQueue;
    class WorldAmInterface;
    class WorldGopInterface;

    void redirectio(World& world);

    /// Initialize the MADNESS runtime

    /// Call this once at the very top of your main program to initialize the
    /// MADNESS runtime. Call this function instead of \c MPI_Init() or
    /// \c MPI_Init_thread() .
    /// \param argc Application argument count
    /// \param argv Application argument values
    /// \return A reference to the default world which is constructed with
    /// \c MPI_COMM_WORLD .
    World& initialize(int& argc, char**& argv);

    /// Initialize the MADNESS runtime

    /// Call this once at the very top of your main program to initialize the
    /// MADNESS runtime. Call this function instead of \c MPI_Init() or
    /// \c MPI_Init_thread() .
    /// \param argc Application argument count
    /// \param argv Application argument values
    /// \param comm The communicator that should be used to construct the
    /// default \c World object.
    /// \return A reference to the default world which is constructed with
    /// \c comm .
    World& initialize(int& argc, char**& argv, const SafeMPI::Intracomm& comm);

    /// Initialize the MADNESS runtime

    /// Call this once at the very top of your main program to initialize the
    /// MADNESS runtime. Call this function instead of \c MPI_Init() or
    /// \c MPI_Init_thread() .
    /// \param argc Application argument count
    /// \param argv Application argument values
    /// \param comm The MPI communicator that should be used to construct the
    /// default \c World object.
    /// \return A reference to the default world which is constructed with
    /// \c comm .
    World& initialize(int& argc, char**& argv, const MPI_Comm& comm);

    /// Call this once at the very end of your main program instead of calling MPI_Finalize
    void finalize();

    /// @return true if madness::initialize() had been called more recently than finalize(), false otherwise
    bool initialized();

    /// Call this to print misc. stats ... collective
    void print_stats(World& world);

    extern void xterm_debug(const char* path, const char* display);

    void error(const char *msg);

    template <typename T>
    static void error(const char *msg, const T& data) {
        std::cerr << "MADNESS: fatal error: " << msg << " " << data << std::endl;
        SafeMPI::COMM_WORLD.Abort();
    }


    /// A parallel world with full functionality wrapping an MPI communicator

    /// Multiple worlds with different communicators can co-exist.
    class World : private NO_DEFAULTS {
    private:
        friend class WorldAmInterface;
        friend class WorldGopInterface;
        friend World& initialize(int&, char**&, const SafeMPI::Intracomm&);
        friend void finalize();

        // Static member variables
        static unsigned long idbase;        ///< Base for unique world ID range for this process
        static World* default_world;        ///< Default world
        static std::list<World*> worlds;    ///< Maintains list of active worlds

        struct hashuniqueT {
            inline std::size_t operator()(const uniqueidT& id) const {
                return id.objid;    // The object id's are guaranteed to be unique
            }
        };

        struct hashvoidp {
            inline std::size_t operator()(const void* p) const {
                return std::size_t(p);    // The ptr's are guaranteed to be unique
            }
        };

//        Mutex globalmutex;  ///< Worldwide mutex
        typedef madness::ConcurrentHashMap<uniqueidT, void *, hashuniqueT> map_id_to_ptrT;
        typedef madness::ConcurrentHashMap<void *, uniqueidT, hashvoidp> map_ptr_to_idT;
        map_id_to_ptrT map_id_to_ptr;
        map_ptr_to_idT map_ptr_to_id;


        unsigned long _id;                  ///< Universe wide unique ID of this world
        unsigned long obj_id;               ///< Counter to generate unique IDs within this world
        void* user_state;                   ///< Holds user defined & managed local state

        // Default copy constructor and assignment won't compile
        // (which is good) due to reference members.

    public:
        // Here we use Pimpl to both hide implementation details and also
        // to partition the namespace for users as world.mpi, world.am, etc.
        // We also embed a reference to this instance in the am and task
        // instances so that they have access to everything.
        //
        // The downside is we cannot do much of anything here without
        // using wrapper functions to foward the calls to the hidden
        // class methods.

        // !!! Order of declaration is important for correct order of initialization !!!
        WorldMpiInterface& mpi;  ///< MPI interface
        WorldAmInterface& am;    ///< AM interface
        WorldTaskQueue& taskq;   ///< Task queue
        WorldGopInterface& gop;  ///< Global operations

    private:
        unsigned int myrand_next;///< State of crude internal random number generator

    public:
        /// Give me a communicator and I will give you the world.
        /// Does not check if another world using the same comm already exists (use instance() to check that)
        World(const SafeMPI::Intracomm& comm);

        /// Find the World corresponding to the given communicator

        /// \param comm the communicator
        /// \return nonzero pointer to the World that was constructed from
        /// \c comm ; if it does not exist, return 0.
        static World* find_instance(const SafeMPI::Intracomm& comm) {
            typedef std::list<World*>::const_iterator citer;
            for(citer it = worlds.begin(); it != worlds.end(); ++it) {
                if ((*it)->mpi.comm() == comm)
                    return *it;
            }
            return 0;
        }

        /// Default \c World object accessor

        /// This function returns a reference to the default world object; this
        /// is the same \c World object that was returned by
        /// \c madness::initialize().
        /// \return A reference to the default world.
        static World& get_default() {
            MADNESS_ASSERT(default_world);
            return *default_world;
        }

        /// Sets a pointer to user-managed local state

        /// Rather than having all remotely invoked actions carry all
        /// of their data with them, they can access local state thru
        /// their world instance.  The user is responsible for
        /// consistently managing and freeing this data.
        ///
        /// A more PC C++ style would be for the app to put state in
        /// a singleton.
        void set_user_state(void* state) { user_state = state; }

        /// Returns pointer to user-managed state set by set_user_state()

        /// Will be NULL if set_user_state() has not been invoked.
        void* get_user_state() { return user_state; }

        /// Clears user-defined state ... same as set_user_state(0)
        void clear_user_state() { user_state = NULL; }

        /// Processes command line arguments

        /// Mostly for world test codes but most usefully provides -dx option
        /// to start x debugger.
        void args(int argc, char**argv);

        /// Returns the system-wide unique integer ID of this world
        unsigned long id() const { return _id; }

        /// Returns the process rank in this world (same as MPI_Comm_rank()))
        ProcessID rank() const { return mpi.rank(); }


        /// Returns the number of processes in this world (same as MPI_Comm_size())
        ProcessID nproc() const { return mpi.nproc(); }

        /// Returns the number of processes in this world (same as MPI_Comm_size())
        ProcessID size() const { return mpi.size(); }

        /// Returns new universe-wide unique ID for objects created in this world.  No comms.

        /// You should consider using register_ptr(), unregister_ptr(),
        /// id_from_ptr() and ptr_from_id() rather than using this directly.
        ///
        /// Currently relies on this being called in the same order on
        /// every process within the current world in order to avoid
        /// synchronization.
        ///
        /// The value objid=0 is guaranteed to be invalid.
        uniqueidT unique_obj_id() { return uniqueidT(_id,obj_id++); }



        /// Associate a local pointer with a universe-wide unique id

        /// Use the routines register_ptr(), unregister_ptr(),
        /// id_from_ptr() and ptr_from_id() to map distributed data
        /// structures identified by the unique id to/from
        /// process-local data.
        ///
        /// !! The pointer will be internally cast to a (void *)
        /// so don't attempt to shove member pointers in here.
        ///
        /// !! ALL unique objects of any type within a world must
        /// presently be created in the same order on all processes so
        /// as to provide the uniquess property without global
        /// communication.
        template <typename T>
        uniqueidT register_ptr(T* ptr) {
            MADNESS_ASSERT(sizeof(T*) == sizeof(void *));
            uniqueidT id = unique_obj_id();
            map_id_to_ptr.insert(std::pair<uniqueidT,void*>(id,static_cast<void*>(ptr)));
            map_ptr_to_id.insert(std::pair<void*,uniqueidT>(static_cast<void*>(ptr),id));
            return id;
        }


        /// Unregister a unique id for a local pointer
        template <typename T>
        void unregister_ptr(T* ptr) {
            uniqueidT id = id_from_ptr(ptr);  // Will be zero if invalid
            map_id_to_ptr.erase(id);
            map_ptr_to_id.erase((void *) ptr);
        }


        /// Unregister a unique id for a local pointer based on id

        /// Same as world.unregister_ptr(world.ptr_from_id<T>(id));
        template <typename T>
        void unregister_ptr(const uniqueidT id) {
            T* const ptr = ptr_from_id<T>(id);
            map_id_to_ptr.erase(id);
            map_ptr_to_id.erase((void *) ptr);
        }


        /// Look up local pointer from world-wide unique id.

        /// Returns NULL if the id was not found.
        template <typename T>
        T* ptr_from_id(uniqueidT id) const {
            map_id_to_ptrT::const_iterator it = map_id_to_ptr.find(id);
            if (it == map_id_to_ptr.end())
                return NULL;
            else
                return (T*)(it->second);
        }


        /// Look up id from local pointer

        /// Returns invalid id if the ptr was not found
        template <typename T>
        const uniqueidT& id_from_ptr(T* ptr) const {
            static uniqueidT invalidid(0,0);
            map_ptr_to_idT::const_iterator it = map_ptr_to_id.find(ptr);
            if (it == map_ptr_to_id.end())
                return invalidid;
            else
                return it->second;
        }

#ifndef MADNESS_DISABLE_SHARED_FROM_THIS

        /// Look up local pointer from world-wide unique id.

        /// \return A default constructed std::shared_ptr if the id was not found.
        template <typename T>
        std::shared_ptr<T> shared_ptr_from_id(uniqueidT id) const {
            T* ptr = ptr_from_id<T>(id);
            return (ptr ? ptr->shared_from_this() : std::shared_ptr<T>());
        }

        /// Look up id from local pointer

        /// \return invalid id if the ptr was not found
        template <typename T>
        const uniqueidT& id_from_ptr(std::shared_ptr<T>& ptr) const {
            return id_from_ptr(ptr.get());
        }

#endif // MADNESS_DISABLE_SHARED_FROM_THIS


        /// Convert world id to world pointer

        /// The id will only be valid if the process calling this routine
        /// is a member of that world.  Thus a null return value does not
        /// necessarily mean the world does not exist --- it could just
        /// not include the calling process.
        /// \param id The world id
        /// \return A pointer to the world associated with \c id
        static World* world_from_id(unsigned long id) {
            // This is why C++ iterators are stupid, stupid, stupid, ..., gack!
            for (std::list<World *>::iterator it=worlds.begin(); it != worlds.end(); ++it) {
                if ((*it) && (*it)->_id == id) return *it;
            }
            return NULL;
        }

    private:

        // Cannot use bind_nullary here since SafeMPI::Request::Test is non-const
        struct MpiRequestTester {
            mutable SafeMPI::Request* r;
            MpiRequestTester(SafeMPI::Request& r) : r(&r) {}
            bool operator()() const {
                return r->Test();
            }
        };

    public:

        /// Wait for MPI request to complete
        static void inline await(SafeMPI::Request& request, bool dowork = true) {
            ThreadPool::await(MpiRequestTester(request), dowork);
        }

        /// Gracefully wait for a condition to become true ... executes tasks if any in queue

        /// Probe should be an object that when called returns the status.
        template <typename Probe>
        static void inline await(const Probe& probe, bool dowork = true) {
            ThreadPool::await(probe, dowork);
        }

        void srand(unsigned long seed = 0ul) {
            if (seed == 0) seed = rank();
#ifdef HAVE_RANDOM
            srandom(seed);
#else
            myrand_next = seed;
            for (int i=0; i<1000; ++i) rand(); // Warmup
#endif // HAVE_RANDOM
        }


        /// Returns a CRUDE, LOW-QUALITY, random number uniformly distributed in [0,2**24).

        /// Each process has a distinct seed for the generator.
        int rand() {
#ifdef HAVE_RANDOM
            return int(random() & 0xfffffful);
#else
            myrand_next = myrand_next * 1103515245UL + 12345UL;
            return int((myrand_next>>8) & 0xfffffful);
#endif // HAVE_RANDOM
        }


        /// Returns a CRUDE, LOW-QUALITY, random number uniformly distributed in [0,1).
        double drand() { return rand()/16777216.0; }

        /// Returns a random process number [0,world.size())
        ProcessID random_proc() { return rand()%size(); }

        /// Returns a random process number [0,world.size()) != current process

        /// Makes no sense to call this with just one process, but just in case you
        /// do it returns -1 in the hope that you won't actually use the result.
        ProcessID random_proc_not_me() {
            if (size() == 1) return -1;
            ProcessID p;
            do {
                p = rand()%size();
            } while (p == rank());
            return p;
        }

        ~World();
    }; // class World

    namespace archive {

        template <typename, typename>
        struct ArchiveLoadImpl;
        template <typename, typename>
        struct ArchiveStoreImpl;

        template <class Archive>
        struct ArchiveLoadImpl<Archive,World*> {
            static inline void load(const Archive& ar, World*& wptr) {
                unsigned long id = 0ul;
                ar & id;
                wptr = World::world_from_id(id);
                MADNESS_ASSERT(wptr);
            }
        }; // struct ArchiveLoadImpl<Archive,World*>

        template <class Archive>
        struct ArchiveStoreImpl<Archive,World*> {
            static inline void store(const Archive& ar, World* const & wptr) {
                ar & wptr->id();
            }
        }; // struct ArchiveStoreImpl<Archive,World*>
    } // namespace archive
} // namespace madness

/*@}*/

#endif // MADNESS_WORLD_WORLDFWD_H__INCLUDED
