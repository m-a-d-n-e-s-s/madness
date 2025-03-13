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
*/

/**
 \file world.h
 \brief Declares the \c World class for the parallel runtime environment.
 \ingroup world

 \todo More detailed description of this file.

 \todo Are some of the forward declarations in this file necessary? A quick inspection suggests most of the functions before the World class don't need to be declared first...
*/

#ifndef MADNESS_WORLD_WORLD_H__INCLUDED
#define MADNESS_WORLD_WORLD_H__INCLUDED

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

// Standard C++ header files needed by MADworld.h
#include <cstddef>
#include <cstdint>  // std::uint64_t
#include <iostream>
#include <list>     // std::list
#include <optional> // std::optional
#include <utility>  // std::pair

#ifdef HAVE_RANDOM
#include <cstdlib>
#endif

// Madness world header files needed by world
#include <madness/world/worldinit.h>
#include <madness/world/worldmpi.h>
#include <madness/world/worldhashmap.h>
#include <madness/world/worldprofile.h>
#include <madness/world/thread.h>
#include <madness/world/uniqueid.h>
#include <madness/world/nodefaults.h>

/// \addtogroup world
/// @{

namespace madness {

    class World;
    class WorldTaskQueue;
    class WorldAmInterface;
    class WorldGopInterface;

    /// Print miscellaneous stats on a World.

    /// This requires collective operations within the World.
    /// \param[in] world The World to analyze.
    void print_stats(World& world);

    /// \todo Brief description needed.

    /// \todo Detailed description needed.
    /// \param path Description.
    /// \param display Description.
    extern void xterm_debug(const char* path, const char* display);

    /// \todo Brief description needed.

    /// \todo Detailed description needed. I imagine this prints/logs the supplied
    /// error message and subsequently aborts the program.
    /// \param[in] msg Message associated with the error.
    void error(const char *msg);

    /// Prints an error message, along with some data, and aborts the program.

    /// \todo Does this function need to be static? (Esp. since it's in a header file...)
    ///
    /// \tparam T The type of data.
    /// \param[in] msg The error message.
    /// \param[in] data The data to be printed.
    template <typename T>
    static void error(const char *msg, const T& data);


    /// A parallel world class.

    /// World wraps a MPI communicator. Multiple worlds with different
    /// communicators can co-exist.
    ///
    /// Here we use Pimpl to both hide implementation details and also
    /// to partition the namespace for users as \c world.mpi, \c world.am, etc.
    /// We also embed a reference to this instance in the \c am and \c task
    /// instances so that they have access to everything.
    ///
    /// The downside is we cannot do much of anything here without
    /// using wrapper functions to forward the calls to the hidden
    /// class methods.
    class World : private NO_DEFAULTS {
    private:
        friend class WorldAmInterface;
        friend class WorldGopInterface;
        friend World& initialize(int&, char**&, const SafeMPI::Intracomm&, int, bool);
        friend void finalize();

        // Static member variables
        static std::pair<std::uint64_t, std::uint64_t> world_id__next_last; ///< Unique {next, last} world IDs to be used by this rank.
        static World* default_world; ///< Default world.
        static std::list<World*> worlds; ///< Maintains list of active worlds, EXCLUDES default_world

        /// \todo Brief description needed.
        struct hashuniqueT {
            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] id The ID.
            /// \return Return description needed.
            inline std::size_t operator()(const uniqueidT& id) const {
                return id.objid;    // The object id's are guaranteed to be unique
            }
        };

        /// \todo Brief description needed.
        struct hashvoidp {
            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param[in] p Missing description.
            /// \return Return description needed.
            inline std::size_t operator()(const void* p) const {
                return std::size_t(p);    // The ptr's are guaranteed to be unique
            }
        };

        // Mutex globalmutex;  ///< Worldwide mutex
        /// \todo Brief description of typedef needed.
        typedef madness::ConcurrentHashMap<uniqueidT, void *, hashuniqueT> map_id_to_ptrT;
        /// \todo Brief description of typedef needed.
        typedef madness::ConcurrentHashMap<void *, uniqueidT, hashvoidp> map_ptr_to_idT;

        /// Maps unique ID of a global object to its local pointer
        ///
        /// \note Unique IDs are created by increasing `obj_id`. To be able to assign unique IDs to objects created in
        /// different worlds, they must be created in exactly the same sequence in every process. Hence in general
        /// object construction happens in main thread. Object construction first updates `obj_id`, then inserts ID
        /// into this. The rest of code only knows about this object once it's been added to this map, hence
        /// this map (NOT `obj_id`) defines the known objects. The only writers to map_id_to_ptr are
        /// `{register,unregister}_ptr()`.
        map_id_to_ptrT map_id_to_ptr;
        /// inverse of map_id_to_ptr
        map_ptr_to_idT map_ptr_to_id; ///< \todo Verify: Map from a pointer to its unique hash ID.


        std::uint64_t _id; ///< Universe wide unique ID of this world.
        unsigned long obj_id; ///< Counter for generating unique IDs within this world.
        void* user_state; ///< Holds a user-defined and managed local state.

        // Default copy constructor and assignment won't compile
        // (which is good) due to reference members.

        /// initializes the value of world_id__next_last to {global_rank*2^32, (global_rank+1)*2^32-1}
        /// \param global_rank rank of this process in COMM_WORLD
        static void initialize_world_id_range(int global_rank);
        /// \return next unique World id
        static std::uint64_t next_world_id();



    public:
        // !!! Order of declaration is important for correct order of initialization !!!
        WorldMpiInterface& mpi; ///< MPI interface.
        WorldAmInterface& am; ///< AM interface.
        WorldTaskQueue& taskq; ///< Task queue.
        WorldGopInterface& gop; ///< Global operations.

    private:
        unsigned int myrand_next; ///< State of crude internal random number generator.

    public:
        /// Constructs a \c World from a communicator.

        /// This function does not check if another \c World exists that uses
        /// the same communicator. Use instance() to check this.
        ///
        /// \param[in] comm The communicator.
        /// \param[in] fence if true, will synchronize ranks before exiting;
        ///                  setting to false removes the extra synchronization
        ///                  but may cause premature arrival of RMI messages
        ///                  that refer to this world while it's being
        ///                  registered
        World(const SafeMPI::Intracomm& comm,
              bool fence = true);

        /// Find the World (if it exists) corresponding to the given communicator.

        /// \param[in] comm The communicator.
        /// \return Pointer to the World that was constructed from \c comm;
        ///     if such a World does not exist, return 0.
        static World* find_instance(const SafeMPI::Intracomm& comm) {
          if (default_world->mpi.comm() == comm)
            return default_world;
          else {
            typedef std::list<World *>::const_iterator citer;
            for (citer it = worlds.begin(); it != worlds.end(); ++it) {
              if ((*it)->mpi.comm() == comm)
                return *it;
            }
            return nullptr;
          }
        }

        /// Check if the World exists in the registry.

        /// \param[in] world pointer to a World object
        /// \return true if \c world exists
        static bool exists(World* world) {
          return world == default_world || std::find(worlds.begin(), worlds.end(), world) != worlds.end();
        }

        /// Default World object accessor.

        /// This function returns a reference to the default world object; this
        /// is the same World object that is returned by
        /// madness::initialize().
        /// \return A reference to the default World.
        /// \throw madness::Exception if the MADNESS_DISABLE_WORLD_GET_DEFAULT preprocessor macro is defined
        static World& get_default() {
#ifdef MADNESS_DISABLE_WORLD_GET_DEFAULT
          MADNESS_EXCEPTION("World::get_default() was called while disabled", 0);
#endif
          MADNESS_ASSERT(default_world);
          return *default_world;
        }

        /// Checks if the default World object corresponds to the given Intracomm

        /// \param[in] comm The communicator.
        /// \return true if \c comm is the default World object's communicator
        /// \sa World::get_default()
        static bool is_default(const SafeMPI::Intracomm& comm) {
            auto* comm_world = find_instance(comm);
            if (comm_world)
              return default_world == comm_world;
            else
              return false;
        }

        /// Sets the user-managed local state.

        /// Rather than having all remotely invoked actions carry all
        /// of their data with them, they can access local state through
        /// their \c World instance.  The user is responsible for
        /// consistently managing and freeing this data.
        ///
        /// A more PC C++ style would be for the app to put state in
        /// a singleton.
        /// \param[in] state The user-managed local state.
        void set_user_state(void* state) { user_state = state; }

        /// Returns a pointer to the user-managed local state set by set_user_state().

        /// \return A pointer to the user-managed local state; NULL if
        ///     set_user_state() has not been invoked.
        void* get_user_state() { return user_state; }

        /// Clears the user-defined state.

        /// This has the same effect as `set_user_state(0)`.
        void clear_user_state() { user_state = nullptr; }

        /// Processes command line arguments.

        /// Mostly intended for \c World test codes, but also provides the
        /// `-dx option` to start `x` debugger.
        /// \param[in] argc The number of command-line arguments.
        /// \param[in,out] argv The command-line arguments.
        void args(int argc, char**argv);

        /// Returns the system-wide unique integer ID of this \c World.
        ///
        /// \return The system-wide unique integer ID of this \c World.
        unsigned long id() const { return _id; }

        /// Returns the process rank in this \c World (same as MPI_Comm_rank()).

        /// \return The process rank in this \c World.
        ProcessID rank() const { return mpi.rank(); }

        /// Returns the number of processes in this \c World (same as MPI_Comm_size()).

        /// \return The number of processes in this \c World.
        ProcessID nproc() const { return mpi.nproc(); }

        /// Returns the number of processes in this \c World (same as MPI_Comm_size()).

        /// \return The number of processes in this \c World.
        ProcessID size() const { return mpi.size(); }

        /// Creates a new universe-wide unique ID for objects created in this \c World. No comms.

        /// You should consider using \c register_ptr(), \c unregister_ptr(),
        /// \c id_from_ptr(), or \c ptr_from_id() rather than using this directly.
        ///
        /// Currently relies on this being called in the same order on
        /// every process within the current \c World in order to avoid
        /// synchronization.
        ///
        /// Unique IDs are created by increasing static counter.
        /// The value obj_id=0 is invalid.
        /// \return A new universe-wide unique ID for objects created in this \c World.
        /// \warning This is not re-entrant, should be called from a single (typically, main) thread
        uniqueidT make_unique_obj_id() {
          return uniqueidT(_id,obj_id++);
        }

        /// Reports the next universe-wide unique ID generated by make_unique_obj_id() objects created in this \c World. No comms.

        /// \return A the next universe-wide unique ID for objects created in this \c World.
        /// \warning This is not re-entrant, should be called from a single (typically, main) thread
        uniqueidT next_unique_obj_id() const {
          return uniqueidT(_id,obj_id);
        }

        /// Associate a local pointer with a universe-wide unique ID.

        /// Use the routines \c register_ptr(), \c unregister_ptr(),
        /// \c id_from_ptr(), \c ptr_from_id() to map distributed data
        /// structures identified by the unique ID to/from
        /// process-local data.
        ///
        /// \note The pointer will be internally cast to a (void *),
        /// so don't use member pointers here.
        ///
        /// \note All unique objects of any type within a \c World must
        /// presently be created in the same order on all processes so
        /// as to guarantee the uniqueness without global
        /// communication.
        /// \tparam T The type of data to be associated.
        /// \param[in] ptr Pointer to the data that will be associated
        ///     with a unique ID.
        /// \return The unique ID associated with the supplied data.
        template <typename T>
        [[nodiscard]] uniqueidT register_ptr(T* ptr) {
            MADNESS_ASSERT(sizeof(T*) == sizeof(void *));
            uniqueidT id = make_unique_obj_id();
            [[maybe_unused]] auto&& [it1, inserted1] = map_id_to_ptr.insert(std::pair<uniqueidT,void*>(id,static_cast<void*>(ptr)));
            MADNESS_ASSERT(inserted1);
            [[maybe_unused]] auto&& [it2, inserted2] = map_ptr_to_id.insert(std::pair<void*,uniqueidT>(static_cast<void*>(ptr),id));
            MADNESS_ASSERT(inserted2);
            return id;
        }

        /// Unregister a global object via its the unique ID.

        /// \param[in] id The unique ID of the data to unregister.
        void unregister_ptr(const uniqueidT id) {
          // START critical section
          map_id_to_ptrT::accessor acc;
          [[maybe_unused]] auto found = map_id_to_ptr.find(acc, id);
          MADNESS_ASSERT(found); // make sure id is valid

          void* ptr = static_cast<void*>(acc->second);

          // if NDEBUG is not defined, keep unregistered ptr IDs in map_id_to_ptr else remove them
#ifdef NDEBUG
          map_id_to_ptr.erase(acc);  // END critical section
#else
          MADNESS_ASSERT(ptr);  // ensure ID is not expired
          acc->second = nullptr;
          acc.release();  // END critical section
#endif
          [[maybe_unused]] auto erased2 = map_ptr_to_id.try_erase(ptr);
          MADNESS_ASSERT(erased2);
        }

        /// Look up a local pointer from a world-wide unique ID.

        /// \tparam T The type of the data to look up.
        /// \param[in] id The unique ID of the data.
        /// \return The local pointer (can be \c nullptr if \p id has been deregistered and
        /// `NDEBUG` is not `#define`d); if \p id has not been registered then returns null optional
        template <typename T>
        [[nodiscard]] std::optional<T*> ptr_from_id(uniqueidT id) const {
          map_id_to_ptrT::const_accessor acc;
          [[maybe_unused]] auto found = map_id_to_ptr.find(acc, id);
          if (!found)
            return {};
          else
            return {static_cast<T *>(acc->second)};
        }

        /// Look up an ID from a local pointer.

        /// \tparam T The type of the data to look up.
        /// \param[in] ptr The local pointer.
        /// \return Optional containing the unique ID if the pointer is registered, or a null optional if the pointer is not registered.
        template <typename T>
        [[nodiscard]] std::optional<uniqueidT> id_from_ptr(T* ptr) const {
          map_ptr_to_idT::const_accessor acc;
          [[maybe_unused]] auto found = map_ptr_to_id.find(acc, ptr);
          if (!found)
            return {};
          else
            return {acc->second};
        }

#ifndef MADNESS_DISABLE_SHARED_FROM_THIS

        /// Look up a local pointer from a world-wide unique ID.

        /// \tparam T The type of the data to look up.
        /// \param[in] id The unique ID.
        /// \return The pointer or a default constructed `std::shared_ptr` if
        ///     the ID is not found.
        template <typename T>
        [[nodiscard]] std::shared_ptr<T> shared_ptr_from_id(uniqueidT id) const {
            std::optional<T*> ptr_opt = ptr_from_id<T>(id);
            if (ptr_opt) MADNESS_ASSERT(*ptr_opt);
            return (ptr_opt ? *ptr_opt->shared_from_this() : std::shared_ptr<T>());
        }

        /// Look up a unique ID from a local pointer.

        /// \tparam T The type of the data to look up.
        /// \param[in] ptr The local pointer.
        /// \return Optional containing the unique ID, if the pointer is registered, or a null optional, if the pointer is not registered.
        template <typename T>
        [[nodiscard]] std::optional<uniqueidT> id_from_ptr(std::shared_ptr<T>& ptr) const {
            return id_from_ptr(ptr.get());
        }

#endif // MADNESS_DISABLE_SHARED_FROM_THIS

        /// Returns a vector of all unique IDs in this \c World.
        std::vector<uniqueidT> get_object_ids() const {
            std::vector<uniqueidT> ids;
            for (auto it = map_id_to_ptr.begin(); it != map_id_to_ptr.end(); ++it)
                ids.push_back(it->first);
            return ids;
        }

        /// return a vector containing all world ids
        static std::vector<unsigned long> get_world_ids() {
            std::vector<unsigned long> ids;
            for (auto it = worlds.begin(); it != worlds.end(); ++it)
                ids.push_back((*it)->id());
            return ids;

        }

        /// Convert a \c World ID to a \c World pointer.

        /// The ID will only be valid if the process calling this routine
        /// is a member of that \c World. Thus a \c NULL return value does not
        /// necessarily mean that the \c World does not exist --- it could just
        /// not include the calling process.
        /// \param[in] id The ID of the \c World.
        /// \return A pointer to the world associated with \c id, or \c NULL.
        static World* world_from_id(std::uint64_t id) {
          if (id != 0) {
            for (std::list<World *>::iterator it = worlds.begin();
                 it != worlds.end(); ++it) {
              if ((*it) && (*it)->_id == id)
                return *it;
            }
            return nullptr;
          }
          else
            return default_world;
        }

    private:

        // Cannot use bind_nullary here since SafeMPI::Request::Test is non-const
        /// \todo Brief description needed.
        struct MpiRequestTester {
            mutable SafeMPI::Request* r; ///< \todo Brief description needed.

            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \param r Description needed.
            MpiRequestTester(SafeMPI::Request& r) : r(&r) {}

            /// \todo Brief description needed.

            /// \todo Descriptions needed.
            /// \return Description needed.
            bool operator()() const {
                return r->Test();
            }
        };

    public:

        /// Wait for a MPI request to complete.

        /// \todo Descriptions needed.
        /// \param[in,out] request The MPI request on which to wait.
        /// \param dowork Work while waiting - default is true
        static void inline await(SafeMPI::Request& request, bool dowork = true) {
	  ThreadPool::await(MpiRequestTester(request), dowork, true); // Last arg is sleep=true --- don't hard spin on MPI requests
        }

        /// Gracefully wait for a condition to become true.

        /// In the mean time, execute any tasks in the queue.
        /// \todo Descriptions needed.
        /// \tparam Probe An object that, when called, returns the status.
        /// \param[in] probe The conditional's test.
        /// \param dowork Work while waiting - default is true
	/// \param sleep Sleep instead of spin while waiting - default is false
        template <typename Probe>
	  static void inline await(const Probe& probe, bool dowork = true, bool sleep=false) {
            ThreadPool::await(probe, dowork);
        }

        /// Crude seed function for random number generation.

        /// \param[in] seed The seed.
        /// \todo Since we're switching to C++11, would it be worth using the new C++11 random number generation capabilities?
        void srand(unsigned long seed = 0ul) {
            if (seed == 0) seed = rank();
#ifdef HAVE_RANDOM
            srandom(seed);
#else
            myrand_next = seed;
            for (int i=0; i<1000; ++i) rand(); // Warmup
#endif // HAVE_RANDOM
        }


        /// Returns a CRUDE, LOW-QUALITY, random number (integer) uniformly distributed in [0,2**24).

        /// Each process has a distinct seed for the generator.
        /// \return The random number.
        /// \todo Since we're switching to C++11, would it be worth using the new C++11 random number generation capabilities?
        int rand() {
#ifdef HAVE_RANDOM
            return int(random() & 0xfffffful);
#else
            myrand_next = myrand_next * 1103515245UL + 12345UL;
            return int((myrand_next>>8) & 0xfffffful);
#endif // HAVE_RANDOM
        }


        /// Returns a CRUDE, LOW-QUALITY, random number (real) uniformly distributed in [0,1).

        /// Each process has a distinct seed for the generator.
        /// \return The random number.
        /// \todo Since we're switching to C++11, would it be worth using the new C++11 random number generation capabilities?
        double drand() { return rand()/16777216.0; }

        /// Returns a random process number; that is, an integer in [0,`world.size()`).

        /// \return The random process number.
        ProcessID random_proc() { return rand()%size(); }

        /// Returns a random process number [0,`world.size()`) that is not the calling process.

        /// It doesn't make any sense to call this with just one process, but just in
        /// case, this returns -1 in the hope that you won't actually use the result.
        /// \return The random process number, or -1.
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

        /// Specialization of \c ArchiveLoadImpl for \c World pointers.

        /// Helps in archiving (reading) \c World objects.
        /// \tparam Archive The archive type.
        template <class Archive>
        struct ArchiveLoadImpl<Archive,World*> {
            /// Loads a \c World from the specified archive.

            /// \note Aborts the program if a \c World cannot be read from the archive.
            /// \param[in,out] ar The archive.
            /// \param[out] wptr Pointer to the \c World.
            static inline void load(const Archive& ar, World*& wptr) {
                unsigned long id = 0ul;
                ar & id;
                wptr = World::world_from_id(id);
                MADNESS_ASSERT(wptr);
            }
        }; // struct ArchiveLoadImpl<Archive,World*>

        /// Specialization of \c ArchiveStoreImpl for \c World pointers.

        /// Helps in archiving (writing) \c World objects.
        /// \tparam Archive The archive type.
        template <class Archive>
        struct ArchiveStoreImpl<Archive,World*> {
            /// Writes a \c World to the specified archive.

            /// \param[in,out] ar The archive.
            /// \param[in] wptr Pointer to the \c World.
            static inline void store(const Archive& ar, World* const & wptr) {
                ar & wptr->id();
            }
        }; // struct ArchiveStoreImpl<Archive,World*>
    } // namespace archive



    // implementation of templated functions
    template <typename T>
    static void error(const char *msg, const T& data) {
        std::cerr << "MADNESS: fatal error: " << msg << " " << data << std::endl;
        SafeMPI::COMM_WORLD.Abort();
    }
} // namespace madness

/// @}

#endif // MADNESS_WORLD_WORLD_H__INCLUDED
