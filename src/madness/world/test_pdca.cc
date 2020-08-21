#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/world/MADworld.h>
#include <madness/world/worlddc.h>
#include <madness/world/parallel_dc_archive.h>

void test(madness::World& world)
{
    madness::WorldContainer<long,std::vector<unsigned char>> container(world);
    {
        madness::archive::ContainerRecordOutputArchive ar(world,container,7);
        madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(world, ar);
        par & 1;
    }
    
    {
        madness::archive::ContainerRecordInputArchive ar(world,container,7);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(world, ar);
        int value;
        par & value;
        std::cout << value << std::endl;
    }
    world.gop.fence();
}
void test_multi_work(madness::World& universe, madness::World& subworld, madness::WorldContainer<long,std::vector<unsigned char>>& container) {
    long key = 1+(universe.rank()%2);
    int value;
    {
        madness::archive::ContainerRecordInputArchive ar(subworld,container,key);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(subworld, ar);
        par & value;
        std::cout << subworld.id() << " : subworld read process " << subworld.rank() << " read " << value << std::endl;
    }
    {
        long newkey = key*1000;
        madness::archive::ContainerRecordOutputArchive ar(subworld,container,newkey);
        madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(subworld, ar);
        par & value*1000;
    }
}

void test_multi(madness::World& universe) {

    // From universe put data into container
    madness::WorldContainer<long,std::vector<unsigned char>> container(universe);
    if (universe.rank() == 0) {
        {
            madness::archive::ContainerRecordOutputArchive ar(universe,container,1);
            madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(universe, ar);
            par & 1;
        }
        {
            madness::archive::ContainerRecordOutputArchive ar(universe,container,2);
            madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(universe, ar);
            par & 2;
        }
    }
    universe.gop.fence(); // Rank 0 must finish writing before workers start

    // Make subworlds
    int color = universe.rank() % 2;
    SafeMPI::Intracomm comm = universe.mpi.comm().Split(color, universe.rank() / 2);
    madness::World subworld(comm);
    
    // In subworlds read data from universe container and put results back
    test_multi_work(universe, subworld, container);
    universe.gop.fence(); // Workers must finish before reading the results

    // In universe read back results from subworld computations
    {
        madness::archive::ContainerRecordInputArchive ar(universe,container,1000);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(universe, ar);
        int value;
        par & value;
        madness::print("first result", value);
    }
    {
        madness::archive::ContainerRecordInputArchive ar(universe,container,2000);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(universe, ar);
        int value;
        par & value;
        madness::print("second result", value);
    }
    universe.gop.fence();

        
}

int main(int argc, char** argv) {
    madness::World& world = madness::initialize(argc, argv);

    //madness::archive::xxxtest(world);
    //test(world);
    test_multi(world);
    
    madness::finalize();

    return 0;
}
