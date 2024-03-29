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
#ifndef MADNESS_WORLD_WORLDPROFILE_H__INCLUDED
#define MADNESS_WORLD_WORLDPROFILE_H__INCLUDED

#include <madness/madness_config.h>
#include <madness/world/worldrmi.h>
#include <madness/world/worldtypes.h>
#include <madness/world/worldmutex.h>
#include <string>
#include <vector>

// NEED TO ADD ATTRIBUTION TO SHINY ON SOURCE FORGE

namespace madness {

    class World;

    /// Simple container for parallel profile statistic
    template <typename T>
    struct ProfileStat {
        T value, max, min, sum;  // local value, parallel max, min, sum
        ProcessID pmax, pmin;    // processor with max, min values

        ProfileStat() = default;

        /// Copies local stats into parallel stats in prep for global reduction
        void init_par_stats(ProcessID me) {
            max = min = sum = value;
            pmax = pmin = me;
        }

        /// Reduction of parallel data (max, min, sum)
        void par_reduce(const ProfileStat<T>& other) {
            if (other.max > max) {
                max = other.max;
                pmax = other.pmax;
            }
            if (other.min < min) {
                min = other.min;
                pmin = other.pmin;
            }
            sum += other.sum;
        }

        /// Zeros all data
        void clear() {
            value = max = min = sum = 0;
            pmax = pmin = 0;
        }
    }; // struct ProfileStat

    /// Used to store profiler info
    struct WorldProfileEntry : public Spinlock {
        std::string name;          ///< name of the entry
        static const int MAX_NTHREAD=64;
        int depth[MAX_NTHREAD];             ///< depth of recursive calls by thread (0 if no active calls)

        ProfileStat<unsigned long> count;   ///< count of times called
        ProfileStat<double> xcpu; ///< exclusive cpu time (i.e., excluding calls)
        ProfileStat<double> icpu; ///< inclusive cpu call (i.e., including calls)
        ProfileStat<unsigned long> xnmsg_sent; ///< No. of active messages sent ... exclusive
        ProfileStat<unsigned long> inmsg_sent; ///< No. of active messages sent ... inclusive
        ProfileStat<unsigned long> xnmsg_recv; ///< No. of active messages recv ... exclusive
        ProfileStat<unsigned long> inmsg_recv; ///< No. of active messages recv ... inclusive
        ProfileStat<unsigned long> xnbyt_sent; ///< No. of bytes sent ... exclusive
        ProfileStat<unsigned long> inbyt_sent; ///< No. of bytes sent ... inclusive
        ProfileStat<unsigned long> xnbyt_recv; ///< No. of bytes recv ... exclusive
        ProfileStat<unsigned long> inbyt_recv; ///< No. of bytes recv ... inclusive

        WorldProfileEntry(const char* name = "");

        WorldProfileEntry(const WorldProfileEntry& other);

        WorldProfileEntry& operator=(const WorldProfileEntry& other);

        static bool exclusivecmp(const WorldProfileEntry&a, const WorldProfileEntry& b);

        static bool inclusivecmp(const WorldProfileEntry&a, const WorldProfileEntry& b);

        static bool exclusivebytcmp(const WorldProfileEntry&a, const WorldProfileEntry& b);

        static bool inclusivebytcmp(const WorldProfileEntry&a, const WorldProfileEntry& b);

        void init_par_stats(ProcessID me);

        void par_reduce(const WorldProfileEntry& other);

        void clear();

        template <class Archive>
        void serialize(const Archive& ar) {
            ar & name & depth & count & xcpu & icpu & xnmsg_sent & inmsg_sent & xnmsg_recv & inmsg_recv & xnbyt_sent & inbyt_sent & xnbyt_recv & inbyt_recv;
        }
    }; // struct WorldProfileEntry


    /// Singleton-like class for holding profiling data and functionality

    /// Use the macros PROFILE_FUNC, PROFILE_BLOCK, PROFILE_MEMBER_FUNC
    class WorldProfile {
        //static ConcurrentHashMap<std::string,WorldProfileEntry> items;
        static std::vector<WorldProfileEntry> items;
        static Spinlock mutex;
        static double cpu_start;
        static double wall_start;

        static std::vector<WorldProfileEntry>& nvitems();


        /// Returns id of the entry associated with the name.  Returns -1 if not found;
        static int find(const std::string& name);


    public:
        /// Returns id for the name, registering if necessary.
        static int register_id(const char* name);

        /// Returns id for the name, registering if necessary.
        static int register_id(const char* classname, const char* function);

        /// Clears all profiling information
        static void clear();

        /// Returns a reference to the specified entry.  Throws if id is invalid.
        static WorldProfileEntry& get_entry(int id);

        /// Prints global profiling information.  Global fence involved.  Implemented in worldstuff.cc
        static void print(World& world);

    private:
        /// Private.  Accumlates data from process into parallel statistics.  Implemented in worldstuff.cc
        static void recv_stats(World& world, ProcessID p);
    };


    class WorldProfileObj {
        static thread_local WorldProfileObj* call_stack;  ///< Current top of this thread's call stack
        static thread_local int mythreadid; ///< My unique thread id
        WorldProfileObj* const prev; ///< Pointer to the entry that called me
        const int id;                ///< My entry in the world profiler
        const double cpu_base;       ///< Time that I started executing
        RMIStats stats_base;         ///< Msg stats when I start executing
        double cpu_start;            ///< Time that I was at top of stack
        RMIStats stats_start;        ///< Msg stats when I was at top of stack;
    public:

        WorldProfileObj(int id);

        /// Pause profiling while we are not executing ... accumulate time in self
        void pause(double now, const RMIStats& stats);

        /// Resume profiling
        void resume(double now, const RMIStats& stats);

        ~WorldProfileObj();
    };
}

#ifdef WORLD_PROFILE_ENABLE
#  define PROFILE_STRINGIFY(s) #s

#  define PROFILE_BLOCK(name)                                             \
    static const int __name##_id=madness::WorldProfile::register_id(PROFILE_STRINGIFY(name)); \
    madness::WorldProfileObj name(__name##_id)

#  define PROFILE_FUNC                                                    \
    static const int __profile_id=madness::WorldProfile::register_id(__FUNCTION__); \
    madness::WorldProfileObj __profile_obj(__profile_id)

#  define PROFILE_MEMBER_FUNC(classname)                                       \
    static const int __profile_id=madness::WorldProfile::register_id(PROFILE_STRINGIFY(classname),  __FUNCTION__); \
    madness::WorldProfileObj __profile_obj(__profile_id)


#else

#  define PROFILE_BLOCK(name)
#  define PROFILE_FUNC
#  define PROFILE_MEMBER_FUNC(classname)

#endif

#endif // MADNESS_WORLD_WORLDPROFILE_H__INCLUDED
