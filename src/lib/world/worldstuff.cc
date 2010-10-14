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


#include <madness_config.h>
#include <world/world.h>
#include <world/worldmem.h>
#include <cstdlib>
#include <sstream>

/// \file worldstuff.cc
/// \brief Static variables/functions that must be linked in

#ifdef __CYGWIN__
#include <windows.h>
#endif

namespace madness {


    //    SharedCounter future_count;

    static double start_cpu_time;
    static double start_wall_time;
    const int WorldAmInterface::NSEND;

    void error(const char *msg) {
        std::cerr << "MADNESS: fatal error: " << msg << std::endl;
        MPI_Abort(MPI_COMM_WORLD,1);
    }

    void initialize(int argc, char** argv) {
        start_cpu_time = cpu_time();
        start_wall_time = wall_time();
#ifdef HAVE_PAPI
        initialize_papi();
#endif

        bool bind[3];
        int cpulo[3];

        const char* sbind = getenv("MAD_BIND");
        if (!sbind) sbind = MAD_BIND_DEFAULT;
        std::istringstream s(sbind);
        for (int i=0; i<3; i++) {
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

#ifdef SERIALIZE_MPI
        int required = MPI::THREAD_SERIALIZED;
#else
        int required = MPI::THREAD_MULTIPLE;
#endif
        int provided = MPI::Init_thread(argc, argv, required);
        int me = MPI::COMM_WORLD.Get_rank();
        if (provided < required && me == 0) {
            std::cout << "!! Warning: MPI::Init_thread did not provide requested functionality " << required << " " << provided << std::endl;
        }

        ThreadPool::begin();        // Must have thread pool before any AM arrives
        RMI::begin();               // Must have RMI while still running single threaded
        if (me == 0) std::cout << "Runtime initialized with " << ThreadPool::size() << " threads in the pool and affinity " << sbind << "\n";

#ifdef HAVE_PAPI
        begin_papi_measurement();
#endif
    }

    void finalize() {
        RMI::end();
        ThreadPool::end(); // 8/Dec/08 : II added this line as trial
        MPI::Finalize();
    }

    std::list<World*> World::worlds;
    unsigned long World::idbase = 0;
    bool TaskInterface::debug = false;

    // Enables easy printing of MadnessExceptions
    std::ostream& operator<<(std::ostream& out, const MadnessException& e) {
        out << "MadnessException : ";
        if (e.msg) out << "msg=" << e.msg << " : ";
        if (e.assertion) out << "assertion=" << e.assertion << " : ";
        out << "value=" << e.value << " : ";
        if (e.line) out << "line=" << e.line << " : ";
        if (e.function) out << "function=" << e.function << " : ";
        if (e.filename) out << "filename='" << e.filename << "'";
        out << std::endl;
        return out;
    }


    std::ostream& operator<<(std::ostream& s, const uniqueidT& id) {
        s << "{" << id.get_world_id() << "," << id.get_obj_id() << "}";
        return s;
    }



    double wall_time() {
#ifdef __CYGWIN__
        static bool initialized = false;
        static double rfreq;
        if (!initialized) {
            _LARGE_INTEGER freq;
            if (QueryPerformanceFrequency(&freq))
                rfreq = 1.0/double(freq.QuadPart);
            else
                rfreq = 0.0;
            initialized = true;
        }
        _LARGE_INTEGER ins;
        QueryPerformanceCounter(&ins);
        return double(ins.QuadPart)*rfreq;
#else
        static bool first_call = true;
        static double start_time;

        struct timeval tv;
        gettimeofday(&tv,0);
        double now = tv.tv_sec + 1e-6*tv.tv_usec;

        if (first_call) {
            first_call = false;
            start_time = now;
        }
        return now - start_time;
#endif
    }

    double cpu_frequency() {
        static double freq = -1.0;
        if (freq == -1.0) {
            double used = wall_time();
            uint64_t ins = cycle_count();
            if (ins == 0) return 0;
            while ((cycle_count()-ins) < 100000000);  // 100M cycles at 1GHz = 0.1s
            ins = cycle_count() - ins;
            used = wall_time() - used;
            freq = ins/used;
        }
        return freq;
    }

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<void>& f) {
        out << "<void>";
        return out;
    }

    template <>
    std::ostream& operator<<(std::ostream& out, const Future<Void>& f) {
        out << "<Void>";
        return out;
    }

    void World::args(int argc, char** argv) {
        for (int arg=1; arg<argc; arg++) {
            if (strcmp(argv[arg],"-dx")==0) xterm_debug("world", 0);
//             if (strcmp(argv[arg],"-dam")==0) am.set_debug(true);
//            if (strcmp(argv[arg],"-dmpi")==0) mpi.set_debug(true);
//             if (strcmp(argv[arg],"-dref")==0) mpi.set_debug(true);
        }
    }


    __thread WorldProfileObj* WorldProfileObj::call_stack = 0;

    Spinlock WorldProfile::mutex;
    volatile std::vector<WorldProfileEntry> WorldProfile::items;
    double WorldProfile::cpu_start = madness::cpu_time();
    double WorldProfile::wall_start = madness::wall_time();


#ifdef WORLD_PROFILE_ENABLE
    static void profile_do_print(World& world, const std::vector<WorldProfileEntry>& v) {
        double cpu_total = 0.0;
        for (unsigned int i=0; i<v.size(); i++)
            cpu_total += v[i].xcpu.sum;

        double cpu_sum = 0.0;
        std::printf(" cum%% cpu%%   cpu/s   cpu-min  cpu-avg  cpu-max  cpu-eff   inc/s   inc-min  inc-avg  inc-max  inc-eff   calls  call-min call-avg call-max call-eff name\n");
        std::printf(" ---- ---- -------- -------- -------- -------- -------- -------- -------- -------- -------- --------  ------- -------- -------- -------- -------- --------------------\n");

        //
        for (unsigned int i=0; i<v.size(); i++) {
            double cpu = v[i].xcpu.sum;
            double inc = v[i].icpu.sum;
            double count = v[i].count.sum;

            cpu_sum += cpu;

            double cum_cpu_percent = cpu_total ? 100.0*cpu_sum/cpu_total : 0.0;
            double cpu_percent = cpu_total ? 100.0*cpu/cpu_total : 0.0;

            double cpu_mean = cpu/world.size();
            double count_mean = count/world.size();
            double count_eff = v[i].count.max ? count_mean/v[i].count.max : 1.0;
            double cpu_eff = v[i].xcpu.max ? cpu_mean/v[i].xcpu.max : 1.0;

            double inc_mean = inc/world.size();
            double inc_eff = v[i].icpu.max ? inc_mean/v[i].icpu.max : 1.0;

            printf("%5.1f%5.1f%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e%9.2e %s\n",
                   cum_cpu_percent,
                   cpu_percent,
                   cpu, v[i].xcpu.min, cpu_mean, v[i].xcpu.max, cpu_eff,
                   inc, v[i].icpu.min, inc_mean, v[i].icpu.max, inc_eff,
                   double(count), double(v[i].count.min), count_mean, double(v[i].count.max), count_eff,
                   v[i].name.c_str());
            printf("                %9d         %9d                  %9d         %9d                  %9d         %9d\n",
                   v[i].xcpu.pmin, v[i].xcpu.pmax,
                   v[i].icpu.pmin, v[i].icpu.pmax,
                   v[i].count.pmin, v[i].count.pmax);
        }
    }
#endif

    static void est_profile_overhead() {
        PROFILE_MEMBER_FUNC(WorldProfile);
    }

    void print_stats(World& world) {
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
        world.gop.sum(nmsg_sent);
        world.gop.sum(nmsg_recv);
        world.gop.sum(nbyte_sent);
        world.gop.sum(nbyte_recv);

        double max_nmsg_sent = rmi.nmsg_sent;
        double max_nmsg_recv = rmi.nmsg_recv;
        double max_nbyte_sent = rmi.nbyte_sent;
        double max_nbyte_recv = rmi.nbyte_recv;
        world.gop.max(max_nmsg_sent);
        world.gop.max(max_nmsg_recv);
        world.gop.max(max_nbyte_sent);
        world.gop.max(max_nbyte_recv);

        double min_nmsg_sent = rmi.nmsg_sent;
        double min_nmsg_recv = rmi.nmsg_recv;
        double min_nbyte_sent = rmi.nbyte_sent;
        double min_nbyte_recv = rmi.nbyte_recv;
        world.gop.min(min_nmsg_sent);
        world.gop.min(min_nmsg_recv);
        world.gop.min(min_nbyte_sent);
        world.gop.min(min_nbyte_recv);

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
        for (int i=0; i<NUMEVENTS; i++) {
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
            for (int i=0; i<NUMEVENTS; i++) {
                printf("  %3d   #events per node    %.2e / %.2e / %.2e\n",
                       i, min_val[i], val[i]/world.size(), max_val[i]);
            }
            for (int i=0; i<NUMEVENTS; i++) {
                printf("  %3d #events systemwide    %.2e\n", i, val[i]);
            }
            if (total_wall_time > 0) {
                for (int i=0; i<NUMEVENTS; i++) {
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
}
