

/// \file examples/hello.cc
/// \brief Simplest example program for MADNESS
/// \defgroup hellowworldmad Hello world MADNESS style
/// \ingroup examples
///
/// Simplest program that initializes the MADNESS parallel runtime
/// using initialize(), makes a madness::World object, prints 
/// a greeting, and then cleans up.
///
/// To initialize the MADNESS numerical environment you also need
/// \c startup(world,argc,argv) and should include mra/mra.h rather
/// than world/world.h .

#include <world/world.h>

using namespace madness;

int main(int argc, char**argv) {
    initialize(argc, argv);
    World world(MPI::COMM_WORLD);
    
    print("Hello from processor",world.rank());

    finalize();
    return 0;
}
