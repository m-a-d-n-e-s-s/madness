#include <madness_config.h>
#ifdef HAVE_PAPI
#include <papi.h>
#include <world/worldthread.h>
#include <world/worldpapi.h>
namespace madness {
    

    static int events[NUMEVENTS] = {PAPI_FP_OPS};
    static Mutex papi_mutex;
    static volatile long long total_values[NUMEVENTS] = {0};
    static volatile int th=0;

    void initialize_papi() {
        if (PAPI_library_init(PAPI_VER_CURRENT) != PAPI_VER_CURRENT)
            throw "Could not init PAPI";
        if (PAPI_thread_init(pthread_self) != PAPI_OK)
            throw "Could not init PAPI thread API";

        reset_papi_measurement();
    }

    void begin_papi_measurement() {
        PAPI_start_counters(events, NUMEVENTS);
    }

    void end_papi_measurement() {
        // Ignore PAPI errors since don't want to kill the calculation
        // at its very end
        long long values[NUMEVENTS];
        PAPI_stop_counters(values, NUMEVENTS);
        papi_mutex.lock();
        th++;
        std::cout << "THREAD PAPI " << values[0] << std::endl;
        for (int i=0; i<NUMEVENTS; i++) total_values[i] += values[i];
        papi_mutex.unlock();
    }

    void reset_papi_measurement() {
        for (int i=0; i<NUMEVENTS; i++) total_values[i] = 0;
    }

    const long long* get_papi_measurement() {
        std::cout << "PAPI THREAD COUNT = " << th <<std::endl;
        return const_cast<long long*>(total_values);
    }
}

#endif
