#include <madness_config.h>
#ifdef HAVE_PAPI
#include <papi.h>
#include <world/worldthread.h>
#include <world/worldpapi.h>
namespace madness {
    
    static int events[NUMEVENTS] = {PAPI_FP_OPS};
    static Mutex papi_mutex;
    static volatile long long total_values[NUMEVENTS] = {0};
    static __thread long long values[NUMEVENTS];

    void begin_papi_measurement() {
        // Ignore PAPI errors since don't want to kill the calculation
        PAPI_start_counters(events, NUMEVENTS);
    }

    void end_papi_measurement() {
        PAPI_stop_counters(values, NUMEVENTS);
        papi_mutex.lock();
        for (int i=0; i<NUMEVENTS; i++) total_values[i] += values[i];
        papi_mutex.unlock();
    }

    void reset_papi_measurement() {
        for (int i=0; i<NUMEVENTS; i++) values[i] = 0;
    }

    const long long* get_papi_measurement() {
        return const_cast<long long*>(values);
    }
}

#endif
