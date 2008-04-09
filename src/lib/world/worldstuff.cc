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

#ifdef _CRAY
#include <catamount/dclock.h>
#endif

#ifdef __CYGWIN__
#include <windows.h>
#endif

namespace madness {
    
    std::list<World*> World::worlds;
    unsigned long World::idbase = 0;
    uint64_t World::poll_delay = 0;
    uint64_t World::last_poll = 0;
    
    Tag WorldMpiInterface::dynamic_tag_general  = DYNAMIC_TAG_BASE;
    Tag WorldMpiInterface::dynamic_tag_reserved = DYNAMIC_TAG_BASE+1;
    
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
    

    const TaskAttributes& task_attr_generator() {
        static const TaskAttributes attr(true,false);
        return attr;
    }


    // Internal: To handle order of class definitions
    void task_ready_callback_function(WorldTaskQueue* taskq, TaskInterface* task) {
        taskq->add_ready_task(task);
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
    
    double cpu_time() {
#ifdef _CRAY
        return dclock();
#else
        return double(clock())/CLOCKS_PER_SEC;
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
            if (strcmp(argv[arg],"-dam")==0) am.set_debug(true);
            if (strcmp(argv[arg],"-dmpi")==0) mpi.set_debug(true);
            if (strcmp(argv[arg],"-dref")==0) mpi.set_debug(true);
        }
    }

    void watchdog_bark(World& world, double time) {
        std::cerr << world.rank() << ": World: watchdog: I've been idle for " << time << "s\n";
        world.am.print_stats();
        world.taskq.print_stats();
    }

}
