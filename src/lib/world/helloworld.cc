#include <world/world.h>

int main(int argc, char** argv) {
    MPI::Init(argc, argv);
    madness::World world(MPI::COMM_WORLD);

    std::cout << "Hello from " << world.rank() << std::endl;

    MPI::Finalize();
}
