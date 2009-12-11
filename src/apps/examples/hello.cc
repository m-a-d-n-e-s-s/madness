
/// \file examples/hello.cc
/// \brief Simplest example program for MADNESS

#include <world/world.h>

using namespace madness;

int main(int argc, char**argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    
    print("Hello from processor",world.rank());

    finalize();
    return 0;
}
