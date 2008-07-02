#include <world/world.h>

using namespace madness;


void dave(int i) {
    PROFILE_FUNC;
    if (i&1) {PROFILE_BLOCK(mary);}
    if (i&3) {PROFILE_BLOCK(calvin);}
}

void fred(int i) {
    PROFILE_FUNC;
    dave(i);
}

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    
    World world(MPI::COMM_WORLD);
    for (int i=0; i<1000; i++)
        fred(i);
    
    WorldProfile::print(world);
    
    MPI::Finalize();
    return 0;
}
