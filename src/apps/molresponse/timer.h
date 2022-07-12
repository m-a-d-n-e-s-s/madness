// Copyright 2021 Adrian Hurtado
#ifndef SRC_APPS_MOLRESPONSE_TIMER_H_
#define SRC_APPS_MOLRESPONSE_TIMER_H_

#include <vector>

#include "TDDFT.h"
#include "madness/mra/mra.h"
// Needed for timers
namespace molresponse {
double pop(std::vector<double>& v);
// Pulled from SCF.cc, starts a timer
static std::vector<double> ttt, sss;
void start_timer(World& world);

// Stops a timer
void end_timer(World& world, const char* msg);
}  // namespace molresponse

#endif  // SRC_APPS_MOLRESPONSE_TIMER_H_
