#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/MADworld.h>
#include <madness/world/parallel_dc_archive.h>

int main(int argc, char** argv) {
    madness::World& world = madness::initialize(argc, argv);

    madness::archive::xxxtest(world);

    world.gop.fence();
    madness::finalize();

    return 0;
}
