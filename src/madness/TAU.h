#ifndef __TAU_MADNESS_H__
#define __TAU_MADNESS_H__

#ifndef PROFILING_ON

#define TAU_START(a) 
#define TAU_STOP(a)

#else 

#define TAU_ENABLED
#include <Profile/Profiler.h>

#endif /* PROFILING_ON */

#endif /* __TAU_MADNESS_H__ */
