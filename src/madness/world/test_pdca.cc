
//#define MAD_ARCHIVE_DEBUG_ENABLE
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


void test_multi_work(madness::World& universe, madness::World& subworld,
		madness::WorldContainer<long,std::vector<unsigned char>>& container) {

	long key = 1+(universe.rank()%2);
    int value;
    {
        madness::WorldContainer<int,double> data(subworld);
        madness::archive::ContainerRecordInputArchive ar(subworld,container,key);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(subworld, ar);
        par & value & data;
        std::cout << subworld.id() << " : subworld read process " << subworld.rank() << " read " << value << std::endl;
        std::cout << "data.size() " << subworld.id() << " " << subworld.rank() << " " << data.size() << std::endl;
    }
    {
        madness::WorldContainer<int,double> subdata(subworld);
        subdata.replace(subworld.rank(), subworld.rank()*100.0);
        subworld.gop.fence();
        
        long newkey = key*1000;
        madness::archive::ContainerRecordOutputArchive ar(subworld,container,newkey);
        madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(subworld, ar);
        par & value*1000 & subdata;
    }
}

void test_multi(madness::World& universe) {

    // For testing passing of containers from universe to subworld
    madness::WorldContainer<int,double> data(universe);
    data.replace(universe.rank(), universe.rank()*100.0);
    universe.gop.fence();
    int sz=data.size();
    std::cout << "universe data size in rank " << universe.rank() << " "  << sz << std::endl;
    universe.gop.sum(sz);
    std::cout << "universe data size total   "   << sz << std::endl;
    for (auto& a : data) {
    	std::cout << " universe content " << a.first << " " << a.second << std::endl;
    }

    // From universe put data into container
    madness::WorldContainer<long,std::vector<unsigned char>> container(universe);
    {
        madness::archive::ContainerRecordOutputArchive ar(universe,container,1);
        madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(universe, ar);
        par & 1 & data;
    }
    {
        madness::archive::ContainerRecordOutputArchive ar(universe,container,2);
        madness::archive::ParallelOutputArchive<madness::archive::ContainerRecordOutputArchive> par(universe, ar);
        par & 2 & data;
    }
    universe.gop.fence(); // Rank 0 must finish writing before workers start

    // Make subworlds
    int color = universe.rank() % 2;
    SafeMPI::Intracomm comm = universe.mpi.comm().Split(color, universe.rank() / 2);
    madness::World subworld(comm);
    
    // In subworlds read data from universe container and put results back
    test_multi_work(universe, subworld, container);
    subworld.gop.fence();
    universe.gop.fence(); // Workers must finish before reading the results

    // In universe read back results from subworld computations
    {
        madness::WorldContainer<int,double> subdata(universe);
        madness::archive::ContainerRecordInputArchive ar(universe,container,1000);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(universe, ar);
        int value;
        par & value & subdata;
        int sz=subdata.size();
        universe.gop.fence(); // Workers must finish before reading the results
        madness::print("first result", value, sz);
        for (auto& a : subdata) {
        	std::cout << " universe content after work " << a.first << " " << a.second << std::endl;
        }

    }
    {
        madness::WorldContainer<int,double> subdata(universe);
        madness::archive::ContainerRecordInputArchive ar(universe,container,2000);
        madness::archive::ParallelInputArchive<madness::archive::ContainerRecordInputArchive> par(universe, ar);
        int value;
        par & value & subdata;
        int sz=subdata.size();
        madness::print("second result", value,sz);
    }
    universe.gop.fence();
}

int main(int argc, char** argv) {
    madness::World& world = madness::initialize(argc, argv);

    //madness::archive::xxxtest(world);
//    test(world);
//    test_multi(world);
    test_cloud(world);
    
    madness::finalize();

    return 0;
}
