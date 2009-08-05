#include <world/world.h>

int main(int argc, char** argv) {
    madness::initialize(argc,argv);
    madness::World world(MPI::COMM_WORLD);

    std::cout << "Hello from " << world.rank() << std::endl;

    madness::finalize();
    return 0;
}
