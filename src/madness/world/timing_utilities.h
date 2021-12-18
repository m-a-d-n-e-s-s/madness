//
// Created by Florian Bischoff on 5/15/21.
//

#ifndef MADNESS_TIMING_UTILITIES_H
#define MADNESS_TIMING_UTILITIES_H

namespace madness {
struct timer {
    World &world;
    double ttt, sss;
    bool do_print = true;

    timer(World &world, bool do_print = true) : world(world), do_print(do_print) {
        world.gop.fence();
        ttt = wall_time();
        sss = cpu_time();
    }

    void tag(const std::string msg) {
        world.gop.fence();
        double tt1 = wall_time() - ttt;
        double ss1 = cpu_time() - sss;
        if (world.rank() == 0 and do_print) {
            std::stringstream ss;
            ss << "timer:" << std::setw(30) << msg << std::setw(8) << std::setprecision(2)
                      << std::fixed << ss1 << "s " << tt1 <<"s";
            std::cout << ss.str() << std::endl;
        }
        ttt = wall_time();
        sss = cpu_time();
    }

    void end(const std::string msg) {
        tag(msg);
//        world.gop.fence();
//        double tt1 = wall_time() - ttt;
//        double ss1 = cpu_time() - sss;
//        if (world.rank() == 0 and do_print) printf("timer: %20.20s %8.2fs %8.2fs\n", msg.c_str(), ss1, tt1);
    }
};
}

#endif //MADNESS_TIMING_UTILITIES_H
