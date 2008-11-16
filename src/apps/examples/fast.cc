#define WORLD_INSTANTIATE_STATIC_TEMPLATES  
#include <mra/mra.h>
#include <mra/mra.h>
#include <mra/operator.h>
#include <constants.h>
#include <mra/lbdeux.h>

using namespace madness;
using namespace std;

Void doit(int niter) {
    Tensor<double> c(20,20,20), d(20,20,20), e(20,20,20);
    Tensor<double> a(20,20);
    for (int i=0; i<niter; i++) {
        fast_transform(c,a,d,e);
    }
    return None;
}


int main(int argc, char** argv) {
    initialize(argc, argv);
    try {
        World world(MPI::COMM_WORLD);
        cout.precision(8);

        int niter = 100000;

        world.gop.fence();
        double start = wall_time();
        if (world.rank() == 0) print("starting at", start);
        world.gop.fence();
        for (unsigned int i=0; i<=ThreadPool::size(); i++) world.taskq.add(doit, niter); // Note +1
        world.taskq.fence();
        world.gop.fence();
        double end = wall_time();
        if (world.rank() == 0) print("starting at", end);
        world.gop.fence();

        ThreadPool::end();
        print_stats(world);

        if (world.rank() == 0) {
            double nflop = (ThreadPool::size()+1.0)*20.0*20.0*20.0*20.0*2.0*3.0*niter;
            print("");
            print("NFLOPS ", nflop);
            print("  USED ", end-start);
            print("  RATE ", nflop/(end-start));
            print("");
        }
    } catch (const MPI::Exception& e) {
        //        print(e);
        error("caught an MPI exception");
    } catch (const madness::MadnessException& e) {
        print(e);
        error("caught a MADNESS exception");
    } catch (const madness::TensorException& e) {
        print(e);
        error("caught a Tensor exception");
    } catch (const char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (char* s) {
        print(s);
        error("caught a c-string exception");
    } catch (const std::string& s) {
        print(s);
        error("caught a string (class) exception");
    } catch (const std::exception& e) {
        print(e.what());
        error("caught an STL exception");
    } catch (...) {
        error("caught unhandled exception");
    }
    
    finalize();
    return 0;
}
