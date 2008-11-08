#ifndef WORLD_PAPI_H
#define WORLD_PAPI_H

#include <madness_config.h>

#ifdef HAVE_PAPI
namespace madness {
    const int NUMEVENTS = 1;

    void initialize_papi();
    void begin_papi_measurement();
    void end_papi_measurement();
    void reset_papi_measurement();
    const long long* get_papi_measurement();
}
#else
namespace madness {
    const int NUMEVENTS = 0;

    inline void initialize_papi(){};
    inline void begin_papi_measurement(){};
    inline void end_papi_measurement(){};
    inline void reset_papi_measurement(){};
    inline const long long* get_papi_measurement(){return 0;};
}
#endif

#endif
