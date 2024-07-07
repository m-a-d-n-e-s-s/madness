
// Copyright 2021 Adrian Hurtado

#include "molresponse/timer.h"

#include <vector>

#include "madness/mra/mra.h"
// Needed for timers
namespace molresponse {
    double pop(std::vector<double> &v) {
        double x = v.back();
        v.pop_back();
        return x;
    }
    // Pulled from SCF.cc, starts a timer
    void start_timer(madness::World &world) {
        world.gop.fence();
        ttt.push_back(madness::wall_time());
        sss.push_back(madness::cpu_time());
    }
    // Stops a timer
    void end_timer(madness::World &world, const char *msg) {
        world.gop.fence();
        MADNESS_CHECK(!ttt.empty());

        double wall = madness::wall_time() - pop(ttt);
        double cpu = madness::cpu_time() - pop(sss);


        if (world.rank() == 0) printf("   timer: %20.20s %8.2fs %8.2fs\n", msg, cpu, wall);
    }
    void end_timer(madness::World &world, const char *msg, const std::string &key,
                   std::map<std::string, std::pair<double, double>> &time) {
        world.gop.fence();
        double wall = madness::wall_time() - pop(ttt);
        double cpu = madness::cpu_time() - pop(sss);

        std::pair<double, double> timings;
        // first is wall
        // second is cpu
        timings.first = wall;
        timings.second = cpu;

        time[key] = timings;


        if (world.rank() == 0) printf("   timer: %20.20s %8.2fs %8.2fs\n", msg, cpu, wall);
    }
}// namespace molresponse
