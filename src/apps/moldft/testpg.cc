#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <world/world.h>
#include <moldft/pointgroup.h>

using namespace madness;

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(MPI::COMM_WORLD);
    world.gop.fence();

    PointGroup::test();

    madness::finalize();
    return 0;
}
