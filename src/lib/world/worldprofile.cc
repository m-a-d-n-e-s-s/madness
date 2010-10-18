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

  $Id: $
*/

#include <world/worldprofile.h>
#include <world/mpiar.h>
//#include <world/worldtypes.h>
//#include <world/worldtime.h>
//#include <world/worldthread.h>
//#include <world/worldhashmap.h>
#include <world/world.h>

namespace madness {

    WorldProfileEntry::WorldProfileEntry(const char* name)
            : name(name), depth(0) {};

    WorldProfileEntry::WorldProfileEntry(const WorldProfileEntry& other)
            : Spinlock() {
        *this = other;
    }

    WorldProfileEntry& WorldProfileEntry::operator=(const WorldProfileEntry& other) {
        name = other.name;
        depth = other.depth;
        count = other.count;
        xcpu = other.xcpu;
        icpu = other.icpu;
        return *this;
    }

    bool WorldProfileEntry::exclusivecmp(const WorldProfileEntry&a, const WorldProfileEntry& b) {
        return a.xcpu.sum > b.xcpu.sum;
    }

    bool WorldProfileEntry::inclusivecmp(const WorldProfileEntry&a, const WorldProfileEntry& b) {
        return a.icpu.sum > b.icpu.sum;
    }

    void WorldProfileEntry::init_par_stats(ProcessID me) {
        count.init_par_stats(me);
        xcpu.init_par_stats(me);
        icpu.init_par_stats(me);
    }

    void WorldProfileEntry::par_reduce(const WorldProfileEntry& other) {
        count.par_reduce(other.count);
        xcpu.par_reduce(other.xcpu);
        icpu.par_reduce(other.icpu);
    }

    void WorldProfileEntry::clear() {
        count.clear();
        xcpu.clear();
        icpu.clear();
    }


    std::vector<WorldProfileEntry>& WorldProfile::nvitems() {
        return const_cast<std::vector<WorldProfileEntry>&>(items);
    }


    /// Returns id of the entry associated with the name.  Returns -1 if not found;
    int WorldProfile::find(const std::string& name) {
        // ASSUME WE HAVE THE MUTEX ALREADY
        std::vector<WorldProfileEntry>& nv = nvitems();
        size_t sz = nv.size();
        if (sz == 0) nv.reserve(1000); // Avoid resizing during execution ... stupid code somewhere below not thread safe?
        if (sz >=1000) MADNESS_EXCEPTION("WorldProfile: did not reserve enough space!", sz);
        for (unsigned int i=0; i<nv.size(); i++) {
            if (name == nv[i].name) return i;
        }
        return -1;
    }

    /// Returns id for the name, registering if necessary.
    int WorldProfile::register_id(const char* name) {
        ScopedMutex<Spinlock> fred(mutex);
        int id = find(name);
        if (id < 0) {
            std::vector<WorldProfileEntry>& nv = nvitems();
            id = nv.size();
            nv.push_back(name);
        }
        return id;
    }

    /// Returns id for the name, registering if necessary.
    int WorldProfile::register_id(const char* classname, const char* function) {
        ScopedMutex<Spinlock> fred(mutex);
        std::string name = std::string(classname) + std::string("::") + std::string(function);
        int id = find(name.c_str());
        if (id < 0) {
            std::vector<WorldProfileEntry>& nv = nvitems();
            id = nv.size();
            nv.push_back(name.c_str());
        }
        return id;
    }

    /// Clears all profiling information
    void WorldProfile::clear() {
        ScopedMutex<Spinlock> fred(mutex);
        cpu_start = madness::cpu_time();
        wall_start = madness::wall_time();
        std::vector<WorldProfileEntry>& nv = nvitems();
        for (unsigned int i=0; i<nv.size(); i++) {
            nv[i].clear();
        }
    }

    /// Returns a reference to the specified entry.  Throws if id is invalid.
    WorldProfileEntry& WorldProfile::get_entry(int id) {
        std::vector<WorldProfileEntry>& nv = nvitems();
        if (id<0 || id >= int(nv.size())) MADNESS_EXCEPTION("WorldProfileEntry: get_entry: invalid id", id);
        return nv[id];
    }


    void WorldProfile::print(World& world) {
#ifdef WORLD_PROFILE_ENABLE
        for (int i=0; i<100; i++) est_profile_overhead();

        std::vector<WorldProfileEntry>& nv = const_cast<std::vector<WorldProfileEntry>&>(items);

        ProcessID me = world.rank();
        for (unsigned int i=0; i<nv.size(); i++) {
            nv[i].init_par_stats(me);
        }

        recv_stats(world, 2*me+1);
        recv_stats(world, 2*me+2);

        if (me) {
            MPIOutputArchive ar(world, (me-1)/2);
            ar & nv;
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
            std::printf("    o  total  cpu time on process zero %.1f seconds\n", madness::cpu_time()-WorldProfile::cpu_start);
            std::printf("    o  total wall time on process zero %.1f seconds\n", madness::wall_time()-WorldProfile::wall_start);
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

            std::vector<WorldProfileEntry> v(nv);
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
        archive::MPIInputArchive ar(world, p);
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

    WorldProfileObj::WorldProfileObj(int id) : prev(call_stack), id(id), cpu_base(madness::cpu_time()) {
        cpu_start = cpu_base;
        call_stack = this;
        WorldProfile::get_entry(id).depth++; // Keep track of recursive calls to avoid double counting time in self
        if (prev) prev->pause(cpu_start);
    }

    /// Pause profiling while we are not executing ... accumulate time in self
    void WorldProfileObj::pause(double now) {
        ScopedMutex<Spinlock> martha(WorldProfile::get_entry(id));
        WorldProfile::get_entry(id).xcpu.value += (now - cpu_start);
    }

    /// Resume profiling
    void WorldProfileObj::resume(double now) {
        cpu_start = now;
    }

    WorldProfileObj::~WorldProfileObj() {
        // if (call_stack != this) throw "WorldProfileObject: call stack confused\n"; // destructors should not throw
        double now = madness::cpu_time();
        WorldProfileEntry& d = WorldProfile::get_entry(id);
        {
            ScopedMutex<Spinlock> martha(d);
            d.count.value++;
            d.xcpu.value += (now - cpu_start);
            d.depth--;
            if (d.depth == 0) d.icpu.value += (now - cpu_base); // Don't double count recursive calls
        }
        call_stack = prev;
        if (call_stack) call_stack->resume(now);
    }

} // namespace madness
