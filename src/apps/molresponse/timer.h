
#ifndef SRC_APPS_molresponse_TIMER_H
#define SRC_APPS_molresponse_TIMER_H

#include <vector>

#include "TDDFT.h"
#include "madness/mra/mra.h"
// Needed for timers
double pop(std::vector<double>& v);
// Pulled from SCF.cc, starts a timer
std::vector<double> ttt, sss;
void start_timer(World& world);

// Stops a timer
void end_timer(World& world, const char* msg);


#endif
