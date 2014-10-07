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

/*!
  \file worldfwd.h
  \brief Implements World
  \addtogroup parallel_runtime

THIS IS VERY OUT OF DATE AND SOME CONTENT IS NOT APPROPRIATE HERE.
NEEDS FIXING BEFORE THE RELEASE.

The MADNESS parallel programming environment combines several
successful elements from other models and aims to provide a
rich and scalable framework for massively parallel computing while
seamlessly integrating with legacy applications and libraries.
It includes
 - Distributed sparse containers with one-sided access to items,
   transparent remote method invocation, an owner-computes task model,
   and optional user control over placement/distribution.
 - Distributed objects that can be globally addressed.
 - Futures (results of unevaluated expressions) for composition of latency tolerant
   algorithms and expression of dependencies between tasks.
 - Globally accessible task queues in each process which
   can be used individually or collectively to provide a single global
   task queue.
 - Work stealing for dynamic load balancing (prototype being tested)
 - Facile management of computations on processor sub-groups.
 - Integration with MPI
 - Optional integration with Global Arrays (J. Nieplocha,
   http://www.emsl.pnl.gov/docs/global).
 - Active messages to items in a container, distributed objects,
   and processes.
 - Efficient use of multicore processors using pthreads.

\section motivations Motivations and attributions for the parallel runtime

There were several motivations for developing this environment.
 -  The rapid evolution of machines from hundreds (pre-2000), to
    millions (post-2008) of processors demonstrates the need to abandon
    process-centric models of computation and move to paradigms that
    virtualize or even hide the concept of a process.
    The success of applications using the
    Charm++ environment to scale rapidly to 30+K processes and the enormous effort
    required to scale most process-centric applications are the central examples.
 -  The arrival of multi-core processes and the associated needs of
    expressing much more concurrency and adopting techniques for
    latency hiding motivate the use of light weight work queues to
    capture much more concurrency and the use of futures for
    latency hiding.
 -  The complexity of composing irregular applications in partitioned, global-address space
    (PGAS) models using only MPI and/or one-sided memory access (GA, UPC, SHMEM, co-Array)
    motivates the use of an object-centric active-message or remote method invocation (RMI) model
    so that computation may be moved to the data with the same ease as
    which data can be moved.  This greatly simplifies the task of maintaining
    and using distributed data structures.
 -  Interoperability with existing programming models to leverage existing
    functionality and to provide an evolutionary path forward.

The two main early influences for this work were Cilk (Kuszmaul,
http://supertech.csail.mit.edu/cilk) and Charm++ (Kale,
http://charm.cs.uiuc.edu).  Subsequently, ACE (Schmidt,
http://www.cs.wustl.edu/~schmidt/ACE.html), STAPL (Rauchwerger and
Amato, http://parasol.tamu.edu/groups/rwergergroup/research/stapl), and
the HPCS language projects (X10,
http://domino.research.ibm.com/comm/research_projects.nsf/pages/x10.index.html
; Chapel, http://chapel.cs.washington.edu ; Fortress, http://fortress.sunsource.net )
and the amazingly talented teams and individuals developing these.

\par Introduction to the parallel runtime

The entire parallel environment is encapsulated in an instance of the
class \c World which is instantiated by wrapping an MPI communicator.
Multiple worlds may exist, overlap, or be dynamically created and
destroyed.  Distributed containers (currently associative arrays or
hash tables) and distributed objects may be constructed from a world
instance.

The recommended approaches to develop scalable and latency tolerant
parallel algorithms are either object- or task-centric decompositions
rather than the process-centric approach usually forced upon MPI
applications.  The object-centric approach uses distributed containers
(or distributed objects) to store application data.  Computation is
expressed by sending tasks or messages to objects, using the task
queue to automatically manage dependencies expressed via futures.
Placement of data and scheduling/placement of computation can be
delgated to the container and task queue, unless there are spefic
performance concerns in which case the application can have full
knowledge and control of these.

Items in a container may be accessed largely as if in a standard STL
container, but instead of returning an iterator, accessors instead
return a Future<iterator>. A future is a container for the result of a
possibly unevaluated expression.  In the case of an accessor, if the
requested item is local then the result is immediately
available. However, if the item is remote, it may take some time
before the data is made available locally.  You could immediately try
to use the future, which would work but with the downside of
internally waiting for all of the communication to occur.  Much better
is to keep on working and only use the future when it is ready.

By far the best way to compute with futures is to pass them as
arguments to a new task.  Once the futures are ready, the task will be
automatically scheduled for execution.  Tasks that produce a
result also return it as a future, so this same mechanism may be used
to express dependencies between tasks.

Thus, a very natural expression of a parallel algorithm is as a
sequence of dependent tasks.  For example, in MADNESS many of the
algorithms working on distributed, multidimension trees start with
just a single task working on the root of the tree, with all other
processes waiting for something to do.  That one task starts
recursively (depth or breadth first) traversing the tree and
generating new tasks for each node.  These in turn generate more tasks
on their sub-trees.

The \c World.am member provides inter-process active message functionality, which is
the foundation on which everything else is built.  We do not recommend
that applications make routine or direct use of inter-process active messages.
Instead, try to compose applications using messaging
to/between items in distributed containers and the local
task queue(s).

The \c World.mpi member is the preferred way to use MPI since it has a growing
amount of instrumentation and debugging capability, though MPI
routines may be called directly if necessary.  However, MPI is again a
low-level model and we do not encourage its direct use.  It is there
since it is the portable standard for communication and to facilitate
integration with legacy applications.

The \c World.gop member provides global operations that are internally
non-blocking, enabling the invoking thread to continue working.

The execution model is sequentially consistent.  That is,
from the perspective of a single thread of execution, operations
on the same local/remote object behave as if executed sequentially
in the same order as programmed.   This means that performing
a read after a write/modify returns the modified value, as expected.
Such behavior applies only to the view of a single thread ---
the execution of multiple threads and active messages from different
threads may be interleaved arbitrarily.

Creating, executing, and reaping a local, null task with
no arguments or results presently takes about 350ns (Centos 4, 3GHz
Core2, Pathscale 3.0 compiler, -Ofast).  The time
is dominated by \c new and and \c delete of the
task structure, and as such is unlikely to get any faster
except by the application caching and reusing the task structures.
Creating and then executing a chain of
dependent tasks with the result of one task fed as the argument
of the next task (i.e., the input argument is an unevaluated future
which is assigned by the next task) requires about 2000ns per
task, which we believe can be redcued to about 1us (3 GHz Core2).

Creating a remote task adds the
overhead of interprocess communication which is on the scale of 1-3us
(Cray XT).  Note that this is not the actual wall-time latency since
everything is presently performed using asynchronous messaging and
polling via MPI.  The wall-time latency, which is largely irrelevant
to the application if it has expressed enough parallelism, is mostly
determined by the polling interval which is dynamically adjusted
depending upon the amount of local work available to reduce the
overhead from polling.  We can improve the runtime software through better
agregation of messages and use of deeper message queues to reduce the
overhead of remote task creation to essentially that of a local task.

Thus, circa 1us defines the ganularity above which it is worth
considering encapsulating work (c.f., Hockney's n1/2).  However, this
is just considering the balance between overhead incurred v.s. useful
work performed.  The automatic scheduling of tasks dependent upon
future arguments confers many benefits, including
 - hiding the wall-time latency of remote data access,
 - removing from the programmer the burden of correct scheduling
   of dependent tasks,
 - expressing all parallelism at all scales of the algorithm
   for facile scaling to heavily multi-core architectures and
   massively parallel computers, and
 - virtualizing the system resources for maximum
   future portability and scalability.

Available memory limits the number of tasks that can be generated
before any are consumed.  In addition to application specific data,
each task consumes circa 64 bytes on a 64-bit computer.  Thus, a few
hundred thousand outstanding tasks per processor are eminently
feasible even on the IBM BG/L.  Rather than making the application
entirely responsible for throttling it's own task production (which it
can), if the system exceeds more than a user-settable number of
outstanding tasks, it starts to run ready tasks before accepting new
tasks.  The success of this strategy presupposes that there are ready
tasks and that these tasks on average produce less than one new task
with unsatisfied dependencies per task run.  Ultimately, similar to
the Cilk execution model, safe algorithms (in the same sense as safe
MPI programs) must express tasks so that dependencies can be satisfied
without unreasonable expectation of buffering.

In a multiscale approach to parallelism, coarse gain tasks are
first enqueued, and these generate finer-grain tasks, which
in turn generate finer and finer grain work.   [Expand this discussion
and include examples along with work stealing discussion]

\par Distributed Containers (WorldContainer)

The only currently provided containers are associative arrays or maps
that are almost directly equivalent to the STL map or the GNU
hash_map.  Indeed, the implementation can use either of these for the
local storage, though the GNU hash_map is to be preferred for
performance reasons and is the only one discussed here.

A map generalizes the concept of an array (which maps an integer index
in a dense range to a value) by mapping an arbitrary key to a value.
This is a very natural, general and efficient mechanism for storing
sparse data structures.  The distribution of items in the container
between processes is based upon a function which maps the key
to a process.  There is a default mapping which is essentially
a pseudo-random uniform mapping, but the user can provide their own
(possibly data-dependent) operator to control the distribution.

The keys and values associated with containers must be serializble
by the MADNESS archive mechanism.
Please refer to world/archive/archive.h and documentation therein for
information about this.  In addition, the keys must support
 - testing for equality, either by overloading \c == or by
   specializing \c std::equal_to<key_type>, and
 - computing a hash value. See worldhash.h for details.

Here is an example of a key that might be used in an octtree.
\code
struct Key {
   typedef unsigned long ulong;
   ulong n, i, j, k;
   hashT hashval;

   Key() {}

   // Precompute the hash function for speed
   Key(ulong n, ulong i, ulong j, ulong k)
       : n(n), i(i), j(j), k(k), hashval(0)
   {
       madness::hash_combine(hashval, n);
       madness::hash_combine(hashval, i);
       madness::hash_combine(hashval, j);
       madness::hash_combine(hashval, k);
   }

   hashT hash() const {
       return hashval;
   }

   template <typename Archive>
   void serialize(const Archive& ar) {
       ar & n & i & j & k & hashval;
   }

   bool operator==(const Key& b) const {
       // Different keys will probably have a different hash
       return hashval==b.hashval && n==b.n && i==b.i && j==b.j && k==b.k;
   }
};
\endcode

\par Distributed Objects (WorldObject)

Distributed objects (WorldObject) provide all of the communication
and other resources necessary to build new distributed capabilities.
The distributed container class (WorldContainer) actually inherits
most of its functionality from the WorldObject.


\par Static data, etc., for templated classes

Several of the templated classes (currently just the
DistributedContainer, Future and RemoteReference classes) have static
data or helper functions associated with them.  These must be defined
in one and only one file.  To facilitate this definition, the
necessary templates have been wrapped in C-preprocessor conditional
block so that they are only enabled if \c
WORLD_INSTANTIATE_STATIC_TEMPLATES is defined.  In one of your source
(not header) files,
define this macro \em before including \c world.h, and then
instantiate the templates that you are using.

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

        struct hashvoidp {
            inline std::size_t operator()(const void* p) const {
                return std::size_t(p);    // The ptr's are guaranteed to be unique
            }
        };

//        Mutex globalmutex;  ///< Worldwide mutex
        typedef madness::ConcurrentHashMap<uniqueidT, void *, uniqueidT> map_id_to_ptrT;
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
