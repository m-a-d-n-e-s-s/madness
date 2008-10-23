#ifndef WORLD_PROFILE_H
#define WORLD_PROFILE_H

#include <world/worldtypes.h>
#include <world/worldtime.h>
#include <world/worldthread.h>
#include <world/worldhashmap.h>

namespace madness {

    /// Simple container for parallel profile statistic
    template <typename T>
    struct ProfileStat {
        T value, max, min, sum;  // local value, parallel max, min, sum
        ProcessID pmax, pmin;    // processor with max, min values

        /// Constructor initializes all members to zero
        ProfileStat() : value(0), max(0), min(0), sum(0), pmax(0), pmin(0) {}

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

        template <class Archive>
        void serialize(const Archive& ar) {
            ar & value & max & min & sum & pmax & pmin;
        }
    };

    /// Used to store profiler info
    struct WorldProfileEntry {
        std::string name;          ///< name of the entry
        int depth;                 ///< depth of recursive calls (0 if no active calls)
        
        ProfileStat<unsigned long> count;   ///< count of times called
        ProfileStat<double> xcpu; ///< exclusive cpu time (i.e., excluding calls)
        ProfileStat<double> icpu; ///< inclusive cpu call (i.e., including calls)
        
        WorldProfileEntry(const char* name = "") 
            : name(name), depth(0)
        {};

        static bool exclusivecmp(const WorldProfileEntry&a, const WorldProfileEntry& b) {
            return a.xcpu.sum > b.xcpu.sum;
        }

        static bool inclusivecmp(const WorldProfileEntry&a, const WorldProfileEntry& b) {
            return a.icpu.sum > b.icpu.sum;
        }

        void init_par_stats(ProcessID me) {
            count.init_par_stats(me);
            xcpu.init_par_stats(me);
            icpu.init_par_stats(me);
        }

        void par_reduce(const WorldProfileEntry& other) {
            count.par_reduce(other.count);
            xcpu.par_reduce(other.xcpu);
            icpu.par_reduce(other.icpu);
        }
        
        void clear() {
            count.clear();
            xcpu.clear();
            icpu.clear();
        }
        
        template <class Archive>
        void serialize(const Archive& ar) {
            ar & name & depth & count & xcpu & icpu;
        }
    };
    
    /// Singleton-like class for holding profiling data and functionality

    /// Use the macros PROFILE_FUNC, PROFILE_BLOCK, PROFILE_MEMBER_FUNC
    class WorldProfile : private Mutex {
        //static ConcurrentHashMap<std::string,WorldProfileEntry> items;
        static std::vector<WorldProfileEntry> items;
        static double cpu_start;
        static double wall_start;
    public:
        /// Registers the name and returns index of the entry
        static int register_id(const char* name) {
            int id = find(name);
            if (id < 0) {
                id = items.size();
                items.push_back(name);
            }
            return id;
        }

        static int register_id(const char* classname, const char* function) {
            std::string name = std::string(classname) + std::string("::") + std::string(function);
            int id = find(name.c_str());
            if (id < 0) {
                id = items.size();
                items.push_back(name.c_str());
            }
            return id;
        }

        /// Returns id of the entry associated with the name.  Returns -1 if not found;
        static int find(const std::string& name) {
            for (unsigned int i=0; i<items.size(); i++) {
                if (name == items[i].name) return i;
            }
            return -1;
        }
        
        /// Clears all profiling information
        static void clear() {
            cpu_start = madness::cpu_time();
            wall_start = madness::wall_time();
            for (unsigned int i=0; i<items.size(); i++) {
                items[i].clear();
            }
        }
        
        /// Returns a reference to the specified entry.  Throws if id is invalid.
        static WorldProfileEntry& get_entry(int id) {
            // MADNESS_ASSERT not defined where this is needed?????
            if(id<0 || id >= int(items.size())) throw "WorldProfileEntry: get_entry: invalid id";
            return items[id];
        }
        
        /// Prints global profiling information.  Global fence involved.  Implemented in worldstuff.cc
        static void print(World& world);

    private:
        /// Private.  Accumlates data from process into parallel statistics.  Implemented in worldstuff.cc
        static void recv_stats(World& world, ProcessID p);
    };
    
    
    class WorldProfileObj {
        static WorldProfileObj* call_stack;  ///< Current top of the call stack


        //????????????????????????????
        //... need a thread-specific stack ... 


        WorldProfileObj* const prev; ///< Pointer to the entry that called me 
        const int id;                ///< My entry in the world profiler
        const double cpu_base;       ///< Time that I started executing
        double cpu_start;            ///< Time that I was at top of stack
    public:

        WorldProfileObj(int id) : prev(call_stack), id(id), cpu_base(madness::cpu_time()) 
        {
            cpu_start = cpu_base;
            call_stack = this;
            WorldProfile::get_entry(id).depth++; // Keep track of recursive calls to avoid double counting time in self
            if (prev) prev->pause(cpu_start);
        }

        /// Pause profiling while we are not executing ... accumulate time in self
        void pause(double now) {
            WorldProfile::get_entry(id).xcpu.value += (now - cpu_start);
        }

        /// Resume profiling 
        void resume(double now) {
            cpu_start = now;
        }
        
        ~WorldProfileObj() {
            if (call_stack != this) throw "WorldProfileObject: call stack confused\n";
            double now = madness::cpu_time();
            WorldProfileEntry& d = WorldProfile::get_entry(id);
            d.count.value++;
            d.xcpu.value += (now - cpu_start);
            d.depth--;
            if (d.depth == 0) d.icpu.value += (now - cpu_base); // Don't double count recursive calls
            call_stack = prev;
            if (call_stack) call_stack->resume(now);
        }
    };
}

//#define WORLD_PROFILE_ENABLE

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

#endif
