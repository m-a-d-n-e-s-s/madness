
/// \file hello.cc
/// \brief Simplest example program for MADNESS

#include <mra/mra.h>

using namespace madness;

int main(int argc, char**argv) {
    MPI::Init(argc, argv);
    ThreadPool::begin();
    RMI::begin();
    MPI::COMM_WORLD.Barrier();
    World world(MPI::COMM_WORLD);
    
    print("Hello from processor",world.rank());

    RMI::end();
    MPI::Finalize();

    return 0;
}
