
#ifndef SRC_APPS_ADRIAN_TIMER_H
#define SRC_APPS_ADRIAN_TIMER_H

#include <vector>

#include "TDDFT.h"
#include "madness/mra/mra.h"
// Needed for timers
static double pop(std::vector<double>& v);
// Pulled from SCF.cc, starts a timer
static std::vector<double> ttt, sss;
static void start_timer(World& world);  // Stops a timer
static void end_timer(World& world, const char* msg);
#endif
