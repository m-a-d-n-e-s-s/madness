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

  
#include <world/world.h>

/// \file worldstuff.cc
/// \brief Static variables/functions that must be linked in

#ifdef __CYGWIN__
#include <windows.h>
#endif

namespace madness {
    void initialize(int argc, char** argv) {
        bool bind[3] = {true, true, true};
        int cpulo[3] = {0, 1, 2};
        ThreadBase::set_affinity_pattern(bind, cpulo); // Decide how to locate threads before doing anything
        ThreadBase::set_affinity(0);         // The main thread is logical thread 0
        
#ifdef SERIALIZE_MPI    
        int required = MPI::THREAD_SERIALIZED;
#else
        int required = MPI::THREAD_MULTIPLE;
#endif
        int provided = MPI::Init_thread(argc, argv, required);
        if (provided < required && MPI::COMM_WORLD.Get_rank() == 0) {
            std::cout << "!! Warning: MPI::Init_thread did not provide requested functionality" << std::endl;
        }
        
        ThreadPool::begin();        // Must have thread pool before any AM arrives
        RMI::begin();               // Must have RMI while still running single threaded
    }

    void finalize() {
        RMI::end();
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
            while ((cycle_count()-ins) < 10000000);  // 10M cycles at 1GHz = 0.01s
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

#define WORLD_PROFILE_ENABLE
    WorldProfileObj* WorldProfileObj::call_stack = 0;
    std::vector<WorldProfileEntry> WorldProfile::items;
    double WorldProfile::cpu_start = madness::cpu_time();
    double WorldProfile::wall_start = madness::wall_time();


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

    static void est_profile_overhead() {
        PROFILE_MEMBER_FUNC(WorldProfile);
    }

    void WorldProfile::print(World& world) {
#ifdef WORLD_PROFILE_ENABLE
        for (int i=0; i<100; i++) est_profile_overhead();

        ProcessID me = world.rank();
        for (unsigned int i=0; i<items.size(); i++) {
            items[i].init_par_stats(me);
        }
        
        recv_stats(world, 2*me+1);
        recv_stats(world, 2*me+2);

        if (me) {
            MPIOutputArchive ar(world, (me-1)/2);
            ar & items;
        }
        else {
            double overhead = 0.0;
            int overid = find("WorldProfile::est_profile_overhead");
            if (overid != -1) {
                overhead = get_entry(overid).xcpu.sum/get_entry(overid).count.sum;
            }

            std::printf("\n    MADNESS global parallel profile\n");
            std::printf("    -------------------------------\n\n");
            std::printf("    o  estimated profiling overhead %.1e seconds per call\n", overhead);
            std::printf("    o  total  cpu time on process zero %.1e seconds\n", madness::cpu_time()-WorldProfile::cpu_start);
            std::printf("    o  total wall time on process zero %.1e seconds\n", madness::wall_time()-WorldProfile::wall_start);
            std::printf("    o  exclusive cpu time excludes called profiled routines\n");
            std::printf("    o  inclusive cpu time includes called profiled routines and\n");
            std::printf("       does not double count recursive calls\n");
            std::printf("    o  process with max/min value is printed under the entry\n");
            std::printf("    o  first printed with items sorted in descending order by total exclusive\n");
            std::printf("       cpu time and then sorted by total inclusive cpu time\n");
            std::printf("    o  in emacs use toggle-truncate-lines to toggle wrapping long lines\n");
            std::printf("\n");
            std::printf("      cum%% - percent cumulative exclusive cpu time (summed over all nodes)\n");
            std::printf("      cpu%% - percent exclusive cpu time (summed over all nodes)\n");
            std::printf("       cpu - total exclusive cpu time (summed over all nodes)\n");
            std::printf("   cpu-min - minimum exclusive cpu time on any processor\n");
            std::printf("   cpu-avg - mean exclusive cpu time per processor\n");
            std::printf("   cpu-max - maximum exclusive cpu time on any processor\n");
            std::printf("   cpu-eff - cpu efficiency = avg/max\n");
            std::printf("       inc - total inclusive cpu time (summed over all nodes)\n");
            std::printf("   inc-min - minimum inclusive cpu time on any processor\n");
            std::printf("   inc-avg - mean inclusive cpu time per processor\n");
            std::printf("   inc-max - maximum inclusive cpu time on any processor\n");
            std::printf("   inc-eff - inclusive cpu efficiency = avg/max\n");
            std::printf("     calls - total number calls time\n");
            std::printf(" calls-min - minimum number calls on any processor\n");
            std::printf(" calls-avg - mean number calls per processor\n");
            std::printf(" calls-max - maximum number calls on any processor\n");
            std::printf(" calls-eff - calls efficiency = avg/max\n");
            std::printf("\n");

            std::vector<WorldProfileEntry> v(items);
            std::sort(v.begin(), v.end(), &WorldProfileEntry::exclusivecmp);
            std::printf("  ** sorted by exclusive cpu time **\n");
            profile_do_print(world, v);

            std::sort(v.begin(), v.end(), &WorldProfileEntry::inclusivecmp);
            std::printf("  ** sorted by inclusive cpu time **\n");
            profile_do_print(world, v);
            
        }
        world.gop.fence();
        
#endif
    }
    
    void WorldProfile::recv_stats(World& world, ProcessID p) {
        if (p >= world.size()) return;
        MPIInputArchive ar(world, p);
        const std::vector<WorldProfileEntry> v;
        ar & v;
        for (unsigned int i=0; i<v.size(); i++) {
            int id = find(v[i].name);
            if (id != -1) {
                WorldProfileEntry& d = get_entry(id);
                d.par_reduce(v[i]);
            }
            else {
                id = register_id(v[i].name.c_str());
                WorldProfileEntry& d = get_entry(id);
                d = v[i];
            }
        }
    }      

}
