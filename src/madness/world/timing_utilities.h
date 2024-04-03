//
// Created by Florian Bischoff on 5/15/21.
//

#ifndef MADNESS_TIMING_UTILITIES_H
#define MADNESS_TIMING_UTILITIES_H

namespace madness {
struct timer {
    World &world;
    double ttt=0.0, sss=0.0;       // duration
    bool do_print = true;
    bool is_running=false;

    timer(World &world, bool do_print = true) : world(world), do_print(do_print) {
        world.gop.fence();
        resume();
    }

    void resume() {
        if (is_running) print("timer was already running!");
        world.gop.fence();
        ttt-=wall_time();
        sss-=cpu_time();
        is_running=true;
    }

    double interrupt() {
        world.gop.fence();
        ttt+=wall_time();
        sss+=cpu_time();
        is_running=false;
        return sss;
    }

    void print(const std::string msg) const {
        if (world.rank() == 0 and do_print) {
            std::stringstream ss;
            ss << "timer:" << std::setw(30) << msg << std::setw(8) << std::setprecision(2)
                      << std::fixed << sss << "s " << ttt <<"s";
            std::cout << ss.str() << std::endl;
        }
    }

    double tag(const std::string msg) {
        world.gop.fence();
        interrupt();
        print(msg);
        double cpu=sss;
        ttt=0.0;
        sss=0.0;
        resume();
        return cpu;
    }

    double end(const std::string msg) {
        return tag(msg);
//        world.gop.fence();
//        double tt1 = wall_time() - ttt;
//        double ss1 = cpu_time() - sss;
//        if (world.rank() == 0 and do_print) printf("timer: %20.20s %8.2fs %8.2fs\n", msg.c_str(), ss1, tt1);
    }
};
}

#endif //MADNESS_TIMING_UTILITIES_H
