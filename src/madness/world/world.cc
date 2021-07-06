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
 \file world.cc
 \brief Implementation of the \c World class and associated functions.
 \ingroup world
*/

#include <madness/world/world.h>
#include <madness/world/worldmem.h>
#include <madness/world/timers.h>
#include <madness/world/worldam.h>
#include <madness/world/world_task_queue.h>
#include <madness/world/worldgop.h>
#include <cstdlib>
#include <sstream>

#ifdef MADNESS_HAS_ELEMENTAL
#if defined(HAVE_EL_H)
# include <El.hpp>
namespace elem = El;
#elif defined(HAVE_ELEMENTAL_H)
# include <elemental.hpp>
using namespace elem;
#else
# error "MADNESS_HAS_ELEMENTAL set but neither HAVE_EL_H nor HAVE_ELEMENTAL_H set: file an issue at " MADNESS_PACKAGE_URL
#endif
#endif

namespace madness {

    // File scope variables
    namespace {
        double start_cpu_time; ///< \todo Documentation needed.
        double start_wall_time; ///< \todo Documentation needed.
        bool madness_initialized_ = false;  ///< Tracks if MADNESS has been initialized.
        bool madness_quiet_ = false;  ///< Tracks if madness::initialize() requested quiet operation
    } // namespace

    // World static member variables
    std::list<World*> World::worlds; ///< List of \c World pointers in the parallel runtime.
    World* World::default_world = nullptr; ///< The default \c World.
    unsigned long World::idbase = 0; ///< \todo Verify: Base unique ID for objects in the runtime.

    bool initialized() {
      return madness_initialized_;
    }
    bool quiet() {
      return madness_quiet_;
    }

    World::World(const SafeMPI::Intracomm& comm)
            : obj_id(1)          ///< start from 1 so that 0 is an invalid id
            , user_state(0)
            , mpi(*(new WorldMpiInterface(comm)))
            , am(* (new WorldAmInterface(*this)))
            , taskq(*(new WorldTaskQueue(*this)))
            , gop(* (new WorldGopInterface(*this)))
            , myrand_next(0)
    {
        worlds.push_back(this);
        srand();  // Initialize random number generator
        cpu_frequency();

        // Assign a globally (within COMM_WORLD) unique ID to this
        // world by assigning to each processor a unique range of indices
        // and broadcasting from node 0 of the current communicator.
        // Each process in COMM_WORLD is given unique ids for 10K new worlds
        if(idbase == 0 && rank()) {
            idbase = rank()*10000;
        }
        // The id of a new world is taken from the unique range of ids
        // assigned to the process with rank=0 in the sub-communicator
        if(rank() == 0) {
            _id = idbase++;
        }
        // Use MPI for broadcast as incoming messages may try to access an
        // uninitialized world.
        mpi.Bcast(_id, 0);
//        gop.broadcast(_id);
//        gop.barrier();
        am.worldid = _id;

        //std::cout << "JUST MADE WORLD " << id() << " with "
        //          << comm.Get_size() << " members (I am # "
        //          << comm.Get_rank() << ")"<< std::endl;
    }


    void World::args(int argc, char** argv) {
        for (int arg=1; arg<argc; ++arg) {
            if (strcmp(argv[arg],"-dx")==0) xterm_debug("objtest", 0);
//             if (strcmp(argv[arg],"-dam")==0) am.set_debug(true);
//            if (strcmp(argv[arg],"-dmpi")==0) mpi.set_debug(true);
//             if (strcmp(argv[arg],"-dref")==0) mpi.set_debug(true);
        }
    }

    World::~World() {
//        stray WorldObjects are allowed as long as they outlive madness::finalize() :(
//        MADNESS_ASSERT_NOEXCEPT(map_ptr_to_id.size() == 0);
//        MADNESS_ASSERT_NOEXCEPT(map_id_to_ptr.size() == 0);
        worlds.remove(this);
        delete &taskq;
        delete &gop;
        delete &am;
        delete &mpi;
    }

    void error(const char *msg) {
        std::cerr << "MADNESS: fatal error: " << msg << std::endl;
        SafeMPI::COMM_WORLD.Abort(1);
    }


    World& initialize(int& argc, char**& argv, bool quiet) {
        return initialize(argc, argv, SafeMPI::COMM_WORLD, quiet);
    }

    World& initialize(int& argc, char**& argv, const MPI_Comm& comm, bool quiet) {
        return initialize(argc, argv, SafeMPI::Intracomm(comm), quiet);
    }

    World& initialize(int& argc, char**& argv, const SafeMPI::Intracomm& comm, bool quiet) {
        madness_quiet_ = quiet;

#ifdef HAVE_PAPI
        initialize_papi();
#endif

        bool bind[3];
        int cpulo[3];

        const char* sbind = getenv("MAD_BIND");
        if (!sbind) sbind = MAD_BIND_DEFAULT;
        std::istringstream s(sbind);
        for (int i=0; i<3; ++i) {
            int t;
            s >> t;
            if (t < 0) {
                bind[i] = false;
                cpulo[i] = 0;
            }
            else {
                bind[i] = true;
                cpulo[i] = t;
            }
        }

        ThreadBase::set_affinity_pattern(bind, cpulo); // Decide how to locate threads before doing anything
        ThreadBase::set_affinity(0);         // The main thread is logical thread 0

#if defined(HAVE_IBMBGQ) and defined(HPM)
        // HPM Profiler
        // Convention for thread IDs is a bit odd, but reflects their
        // internal labeling in the code.
        // HPM_THREAD_ID = -10, all threads and aggregates
        // HPM_THREAD_ID =  -2, main thread
        // HPM_THREAD_ID =  -1, threads not in pool, i.e. communication
        // HPM_THREAD_ID = 0..MAD_NUM_THREADS - 2, threads in the pool
        int hpm_thread_id;
        char *chpm_thread_id = getenv("HPM_THREAD_ID");
        if (chpm_thread_id) {
            int result = sscanf(chpm_thread_id, "%d", &hpm_thread_id);
            if (result != 1)
                MADNESS_EXCEPTION("HPM_THREAD_ID is not an integer", result);
        }
        ThreadBase::set_hpm_thread_env(hpm_thread_id);
#endif
        detail::WorldMpi::initialize(argc, argv, MADNESS_MPI_THREAD_LEVEL);

        // Construct the default world before starting RMI so that incoming active messages can find this world
        World::default_world = new World(comm);

        start_cpu_time = cpu_time();
        start_wall_time = wall_time();
        ThreadPool::begin();        // Must have thread pool before any AM arrives
        if(comm.Get_size() > 1) {
            RMI::begin(comm);           // Must have RMI while still running single threaded
            // N.B. sync everyone up before messages start flying
            // this is needed to avoid hangs with some MPIs, e.g. Intel MPI on commodity hardware
            comm.Barrier();
        }

#ifdef HAVE_PAPI
        begin_papi_measurement();
#endif // HAVE_PAPI

#ifdef MADNESS_HAS_ELEMENTAL
        elem::Initialize(argc,argv);
#endif // HAVE_ELEMENTAL

        // mark this thread as part of MADNESS pool
        set_thread_tag(ThreadTag_MADNESS | ThreadTag_Main);
        madness_initialized_ = true;
        if(!quiet && comm.Get_rank() == 0)
            std::cout << "MADNESS runtime initialized with " << ThreadPool::size()
                << " threads in the pool and affinity " << sbind << "\n";

        return * World::default_world;
    }

    void finalize() {
        World::default_world->gop.fence();
        const auto rank = World::default_world->rank();
        const auto world_size = World::default_world->size();

        // Destroy the default world
        delete World::default_world;
        World::default_world = nullptr;
        // warn if there are surviving worlds
        if (!World::worlds.empty() && !quiet() && rank == 0) {
            const auto nworlds = World::worlds.size();
            std::cerr << "MADNESS runtime finalized but " << nworlds << " world"
                      << (nworlds > 1 ? "s" : "") << " still exist"
                      << (nworlds > 1 ? "" : "s") << std::endl;
        }

#ifdef MADNESS_HAS_ELEMENTAL
        elem::Finalize();
#endif

        if(world_size > 1)
            RMI::end();
        ThreadPool::end();
        detail::WorldMpi::finalize();
        madness_initialized_ = false;
        madness_quiet_ = false;
    }

    void print_stats(World& world) {
        world.gop.fence();

        double total_wall_time = wall_time()-start_wall_time;
        double total_cpu_time = cpu_time()-start_cpu_time;
        RMIStats rmi = RMI::get_stats();
        DQStats q = ThreadPool::get_stats();
#ifdef HAVE_PAPI
        // For papi ... this only make sense if done once after all
        // other worker threads have exited
        end_papi_measurement();
        const long long* values = get_papi_measurement();
#endif

        double nmsg_sent = rmi.nmsg_sent;
        double nmsg_recv = rmi.nmsg_recv;
        double nbyte_sent = rmi.nbyte_sent;
        double nbyte_recv = rmi.nbyte_recv;
        double server_q = rmi.max_serv_send_q;
        world.gop.sum(nmsg_sent);
        world.gop.sum(nmsg_recv);
        world.gop.sum(nbyte_sent);
        world.gop.sum(nbyte_recv);
        world.gop.sum(server_q);

        double max_nmsg_sent = rmi.nmsg_sent;
        double max_nmsg_recv = rmi.nmsg_recv;
        double max_nbyte_sent = rmi.nbyte_sent;
        double max_nbyte_recv = rmi.nbyte_recv;
        double max_server_q = rmi.max_serv_send_q;
        world.gop.max(max_nmsg_sent);
        world.gop.max(max_nmsg_recv);
        world.gop.max(max_nbyte_sent);
        world.gop.max(max_nbyte_recv);
        world.gop.max(max_server_q);

        double min_nmsg_sent = rmi.nmsg_sent;
        double min_nmsg_recv = rmi.nmsg_recv;
        double min_nbyte_sent = rmi.nbyte_sent;
        double min_nbyte_recv = rmi.nbyte_recv;
        double min_server_q = rmi.max_serv_send_q;
        world.gop.min(min_nmsg_sent);
        world.gop.min(min_nmsg_recv);
        world.gop.min(min_nbyte_sent);
        world.gop.min(min_nbyte_recv);
        world.gop.min(min_server_q);

        double npush_back = q.npush_back;
        double npush_front = q.npush_front;
        double npop_front = q.npop_front;
        double ntask = q.npush_back + q.npush_front;
        double nmax = q.nmax;
        world.gop.sum(npush_back);
        world.gop.sum(npush_front);
        world.gop.sum(npop_front);
        world.gop.sum(ntask);
        world.gop.sum(nmax);

        double max_npush_back = q.npush_back;
        double max_npush_front = q.npush_front;
        double max_npop_front = q.npop_front;
        double max_ntask = q.npush_back + q.npush_front;
        double max_nmax = q.nmax;
        world.gop.max(max_npush_back);
        world.gop.max(max_npush_front);
        world.gop.max(max_npop_front);
        world.gop.max(max_ntask);
        world.gop.max(max_nmax);

        double min_npush_back = q.npush_back;
        double min_npush_front = q.npush_front;
        double min_npop_front = q.npop_front;
        double min_ntask = q.npush_back + q.npush_front;
        double min_nmax = q.nmax;
        world.gop.min(min_npush_back);
        world.gop.min(min_npush_front);
        world.gop.min(min_npop_front);
        world.gop.min(min_ntask);
        world.gop.min(min_nmax);

#ifdef HAVE_PAPI
        double val[NUMEVENTS], max_val[NUMEVENTS], min_val[NUMEVENTS];
        for (int i=0; i<NUMEVENTS; ++i) {
            val[i] = max_val[i] = min_val[i] = values[i];
        }
        world.gop.sum(val, NUMEVENTS);
        world.gop.max(max_val, NUMEVENTS);
        world.gop.min(min_val, NUMEVENTS);
#endif

        if (world.rank() == 0) {
            printf("\n");
            printf("    Parallel environment\n");
            printf("    --------------------\n");
            printf("                  #nodes    %d\n", world.size());
            if (world.size() == 1) {
                printf("       #threads per node    %d+main = %d\n", int(ThreadPool::size()), int(ThreadPool::size()+1));
                printf("          #total threads    %d\n", int(ThreadPool::size()+1));
            }
            else {
                printf("       #threads per node    %d+main+server = %d\n", int(ThreadPool::size()), int(ThreadPool::size()+2));
                printf("          #total threads    %d\n", int(ThreadPool::size()+2)*world.size());
            }
            printf("\n");

            printf("  RMI message statistics (min / avg / max)\n");
            printf("  ----------------------\n");
            printf("   #messages in server q    %.2e / %.2e / %.2e\n",
                   min_server_q, server_q/world.size(), max_server_q);
            printf(" #messages sent per node    %.2e / %.2e / %.2e\n",
                   min_nmsg_sent, nmsg_sent/world.size(), max_nmsg_sent);
            printf("    #bytes sent per node    %.2e / %.2e / %.2e\n",
                   min_nbyte_sent, nbyte_sent/world.size(), max_nbyte_sent);
            printf(" #messages recv per node    %.2e / %.2e / %.2e\n",
                   min_nmsg_recv, nmsg_recv/world.size(), max_nmsg_recv);
            printf("    #bytes recv per node    %.2e / %.2e / %.2e\n",
                   min_nbyte_recv, nbyte_recv/world.size(), max_nbyte_recv);
            printf("        #msgs systemwide    %.2e\n", nmsg_sent);
            printf("       #bytes systemwide    %.2e\n", nbyte_sent);
            printf("\n");
            printf("  Thread pool statistics (min / avg / max)\n");
            printf("  ----------------------\n");
            printf("         #tasks per node    %.2e / %.2e / %.2e\n",
                   min_ntask, ntask/world.size(), max_ntask);
            printf("     #max q len per node    %.2e / %.2e / %.2e\n",
                   min_nmax, nmax/world.size(), max_nmax);
            printf("  #hi-pri tasks per node    %.2e / %.2e / %.2e\n",
                   min_npush_front, npush_front/world.size(), max_npush_front);
            printf("\n");
#ifdef HAVE_PAPI
            printf("         PAPI statistics (min / avg / max)\n");
            printf("         ---------------\n");
            for (int i=0; i<NUMEVENTS; ++i) {
                printf("  %3d   #events per node    %.2e / %.2e / %.2e\n",
                       i, min_val[i], val[i]/world.size(), max_val[i]);
            }
            for (int i=0; i<NUMEVENTS; ++i) {
                printf("  %3d #events systemwide    %.2e\n", i, val[i]);
            }
            if (total_wall_time > 0) {
                for (int i=0; i<NUMEVENTS; ++i) {
                    printf("  %3d   #op/s systemwide    %.2e\n", i, val[i]/total_wall_time);
                }
            }
            printf("\n");
#endif
#ifdef WORLD_GATHER_MEM_STATS
            world_mem_info()->print();
#endif

            printf("         Total wall time    %.1fs\n", total_wall_time);
            printf("         Total  cpu time    %.1fs\n", total_cpu_time);
            printf("\n");
        }
        world.gop.fence();
#ifdef WORLD_PROFILE_ENABLE
        WorldProfile::print(world);
#endif
        world.gop.fence();
    }

} // namespace madness
