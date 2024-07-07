// Copyright 2021 Adrian Hurtado
#ifndef SRC_APPS_MOLRESPONSE_TIMER_H_
#define SRC_APPS_MOLRESPONSE_TIMER_H_

#include <vector>

#include "madness/mra/mra.h"
// Needed for timers
namespace molresponse {
double pop(std::vector<double>& v);
// Pulled from SCF.cc, starts a timer
static std::vector<double> ttt, sss;
void start_timer(madness::World& world);

// Stops a timer
void end_timer(madness::World& world, const char* msg);

void end_timer(madness::World& world, const char* msg, const std::string& key,
               std::map<std::string, std::pair<double, double>>& time);

}  // namespace molresponse

#endif  // SRC_APPS_MOLRESPONSE_TIMER_H_
