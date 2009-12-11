
/// \file ii/hello.cc
/// \brief Simplest example program for MADNESS

#include <mra/mra.h>

using namespace madness;

int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    World world(MPI::COMM_WORLD);
    
    print("Hello from processor",world.rank());

    MPI::Finalize();

    return 0;
}
