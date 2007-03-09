#ifndef WORLD_H
#define WORLD_H

/// \file world.h
/// \brief Implements World and includes pretty much every header you'll need

// must include mpi before io stuff
#include <mpi.h>
#include <iostream>
#include <cstdio>
#include <unistd.h>
#include <assert.h>
#include <algorithm>
#include <functional>
#include <list>
#include <map>
#include <vector>

#ifdef UINT64_T
typedef UINT64_T uint64_t;
#endif

#include <world/typestuff.h>
#include <world/worldhash.h>
#include <world/array.h>
#include <world/print.h>
#include <world/worldexc.h>
#include <world/sharedptr.h>
#include <world/typestuff.h>
#include <world/nodefaults.h>
#include <world/worldmpi.h>
#include <world/worldser.h>
#include <world/worldtime.h>

namespace madness {
    
    extern void xterm_debug(const char* path, const char* display);

    class WorldTaskQueue;
    class WorldAmInterface;
    class WorldGopInterface;
    class World; 
   
    static void world_do_poll(World* world);
    static void world_do_run_task(World* world, bool* status);
    static WorldAmInterface* world_am_interface_factory(World* world);
    static void world_am_interface_unfactory(WorldAmInterface* am);
    static WorldGopInterface* world_gop_interface_factory(World* world);
    static void world_gop_interface_unfactory(WorldGopInterface* gop);
    static WorldTaskQueue* world_taskq_factory(World* world);
    static void world_taskq_unfactory(WorldTaskQueue* taskq);
    static void world_assign_id(World* world);
    

    /// For purpose of deferring cleanup to synchronization points
    struct DeferredCleanupInterface {
        virtual ~DeferredCleanupInterface(){};
    };

    static void error(const char *msg) {
        fprintf(stderr,"fatal error: %s\n",msg);
        MPI_Abort(MPI_COMM_WORLD,1);
    }
    
    // This needs cleaning up and perhaps reference counting 
    // via sharedptr.

    /// A parallel world with full functionality wrapping an MPI communicator

    /// Multiple worlds with different communicators can co-exist.
    class World {
    private:
        friend class WorldAmInterface;
        friend class WorldGopInterface;
        friend void world_assign_id(World* world);

        static unsigned long idbase;        //< Base for unique world ID range for this process
        static std::list<World*> worlds;    //< Maintains list of active worlds for polling, etc.
        static std::list<void (*)()> polls; //< List of routines to invoke while polling
        static uint64_t poll_delay;//< Min. no. of instructions between calls to poll if working
        static uint64_t last_poll;//< Instruction count at last poll
        
        unsigned long _id;                  //< Universe wide unique ID of this world
        unsigned long obj_id;               //< Counter to generate unique IDs within this world
        void* user_state;                   //< Holds user defined & managed local state
        std::list< SharedPtr<DeferredCleanupInterface> > deferred; //< List of stuff to delete at next sync point

        // Default copy constructor and assignment won't compile
        // (which is good) due to reference members.

        /// Does any deferred cleanup and returns true if cleaning was necessary
        bool do_deferred_cleanup() {
            if (deferred.empty()) {
                return false;
            }
            else {
                print("do_deferred_cleanup: cleaning",deferred.size(),"items");
                deferred.clear();
                return true;
            }
        };

        // Private: tries to run a task in each world
        static bool run_tasks() {
            bool status = false;
            for_each(worlds.begin(), worlds.end(), std::bind2nd(std::ptr_fun(world_do_run_task),&status));
            return status;
        };


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
        WorldMpiInterface& mpi;  //< MPI interface
        WorldAmInterface& am;    //< AM interface
        WorldTaskQueue& taskq;   //< Task queue
        WorldGopInterface& gop;  //< Global operations

    private:
        const ProcessID me;      //< My rank ... needs to be declared after MPI
        int nprocess;            //< No. of processes ... ditto

    public:
        /// Give me a communicator and I will give you the world
        World(MPI::Intracomm& comm) 
            : obj_id(0)
            , user_state(0)
            , deferred()
            , mpi(*(new WorldMpiInterface(comm)))
            , am(*world_am_interface_factory(this)) 
            , taskq(*world_taskq_factory(this))
            , gop(*world_gop_interface_factory(this))
            , me(mpi.rank())
            , nprocess(mpi.nproc())
        {
            worlds.push_back(this); 
            
            // Assign a globally (within COMM_WORLD) unique ID to this
            // world by assigning to each processor a unique range of indices
            // and broadcasting from node 0 of the current communicator.
            world_assign_id(this);  // Also acts as barrier 
            
            // Determine cost of polling and from this limit the
            // frequency with which poll_all will be run while there
            // is work in the task queue.
            uint64_t ins = cycle_count();
            for (int i=0; i<32; i++) World::poll_all();
            poll_delay = (cycle_count()-ins)>>5; // Actual cost per poll
            poll_delay = poll_delay<<3;  // *8 --> no more than 12.5% of time in polling
            ///world.mpi.Bcast(poll_delay,0); // For paranoia make sure all have same value?
            if (rank()==0) print("poll_delay",poll_delay);
        };


        /// Sets a pointer to user-managed local state

        /// Rather than having all remotely invoked actions carry all
        /// of their data with them, they can access local state thru
        /// their world instance.  The user is responsible for
        /// consistently managing and freeing this data.
        void set_user_state(void* state) {
            user_state = state;
        };


        /// Returns pointer to user-managed state set by set_user_state()

        /// Will be NULL if set_user_state() has not been invoked.
        void* get_user_state() {
            return user_state;
        };


        /// Clears user-defined state ... same as set_user_state(0)
        void clear_user_state() {
            set_user_state(0);
        };


        /// Invokes any necessary polling for all existing worlds
        static void poll_all(bool working = false) {
            if (working  &&  (cycle_count() < last_poll+poll_delay)) return;
            for_each(worlds.begin(), worlds.end(), world_do_poll);
            last_poll = cycle_count();
        };


        /// Returns the system-wide unique integer ID of this world
        unsigned long id() const {
            return _id;
        };


        /// Returns the process rank in this world (same as MPI::Get_rank()))
        ProcessID rank() const {return me;};


        /// Returns the number of processes in this world (same as MPI::Get_size())
        ProcessID nproc() const {return nprocess;};


        /// Returns unique ID for objects created in this world

        /// Currently relies on this being called in the same order
        /// on every process in order to avoid synchronization.
        /// Ideally, would like to relax this to only in order
        /// for the same type.
        unsigned long unique_obj_id() {
            return obj_id++;
        };


        /// Returns a pointer to the world with given ID or null if not found

        /// The id will only be valid if the process calling this routine
        /// is a member of that world.  Thus a null return value does not
        /// necessarily mean the world does not exist --- it could just
        /// not include the calling process.
        static World* world_from_id(unsigned long id) {
            // This is why C++ iterators are stupid, stupid, stupid, ..., gack!
            for (std::list<World *>::iterator it=worlds.begin(); it != worlds.end(); ++it) {
                if ((*it) && (*it)->_id == id) return *it;
            }
            return 0;
        };


        // Cannot use bind_nullary here since MPI::Request::Test is non-const
        struct MpiRequestTester {
            mutable MPI::Request& r;
            MpiRequestTester(MPI::Request& r) : r(r) {};
            bool operator()() const {return r.Test();};
        };


        /// Wait for MPI request to complete while polling and processing tasks
        static void inline await(MPI::Request& request) {
            await(MpiRequestTester(request));
        };


        /// Wait for a condition to become true while polling and processing tasks

        /// Probe should be an object that when called returns the status.
        ///
        /// Ensures progress is made in all worlds.
        template <typename Probe>
        static void inline await(const Probe& probe) {
            // Critical here is that poll() is NOT called after a
            // successful test of the request since polling may
            // trigger an activity that invalidates the condition.
            bool working = false;
            while (!probe()) {
                poll_all(working);  // If working poll_all will increase polling interval
                working = run_tasks();
            }
        }


        /// Adds item to list of stuff to be deleted at next global_fence()
        void deferred_cleanup(const SharedPtr<DeferredCleanupInterface>& item) {
            deferred.push_back(item);
        };


        ~World() {
            worlds.remove(this);
            do_deferred_cleanup();
            world_taskq_unfactory(&taskq);
            world_gop_interface_unfactory(&gop);
            world_am_interface_unfactory(&am);
            delete &mpi;
        };
    };
}

// Order of these is important
#include <world/worldam.h>
#include <world/worldref.h>
#include <world/worlddep.h>
#include <world/worldfut.h>
#include <world/worldtask.h>
#include <world/worldgop.h>
#include <world/worlddc.h>

namespace madness {

    // This nonsense needs cleaning up and probably eliminating
    // now that the class interfaces have stabilized.

    void redirectio(World& world);

    static inline void world_do_poll(World* world) {
        if (world) world->am.poll();
    }
    static inline void world_do_run_task(World* world, bool* status) {
        if (world) *status = *status || world->taskq.run_next_ready_task();
    }
    static WorldAmInterface* world_am_interface_factory(World* world) {
        return new WorldAmInterface(*world);
    }
    static void world_am_interface_unfactory(WorldAmInterface* am) {
        delete am;
    }
    static WorldGopInterface* world_gop_interface_factory(World* world) {
        return new WorldGopInterface(*world);
    }
    static void world_gop_interface_unfactory(WorldGopInterface* gop) {
        delete gop;
    }
    static WorldTaskQueue* world_taskq_factory(World* world) {
        return new WorldTaskQueue(*world);
    }
    static void world_taskq_unfactory(WorldTaskQueue* taskq) {
        delete taskq;
    }
    static void world_assign_id(World* world) {
        if (world->idbase == 0) world->idbase = MPI::COMM_WORLD.Get_rank()*10000;
        if (world->mpi.rank() == 0) world->_id = World::idbase++;
        world->gop.broadcast(world->_id);
        world->gop.barrier();
    }

    namespace archive {
        template <class Archive>
        struct ArchiveLoadImpl<Archive,World*> {
            static inline void load(const Archive& ar, World*& wptr) {
                unsigned long id;
                ar & id;
                wptr = World::world_from_id(id);
                MADNESS_ASSERT(wptr);
            };
        };
        
        template <class Archive>
        struct ArchiveStoreImpl<Archive,World*> {
            static inline void store(const Archive& ar, World* const & wptr) {
                ar & wptr->id();
            };
        };
    }
}




#endif
