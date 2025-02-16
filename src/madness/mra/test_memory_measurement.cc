//
// Created by Florian Bischoff on 1/9/25.
//

#include<madness/mra/memory_measurement.h>
#include<madness/mra/mra.h>

using namespace madness;

/// for each process create a world using a communicator shared with other processes by round-robin
/// copy-paste from test_world.cc
static std::shared_ptr<World> create_worlds(World& universe, const std::size_t nsubworld) {

    int color = universe.rank() % nsubworld;
    SafeMPI::Intracomm comm = universe.mpi.comm().Split(color, universe.rank() / nsubworld);

    std::shared_ptr<World> all_worlds;
    all_worlds.reset(new World(comm));

    universe.gop.fence();
    return all_worlds;
}

template<std::size_t NDIM>
int test_size(World& world) {

    // create a slater function
    auto slater=[](const Vector<double,2*NDIM>& r){return exp(-r.normf());};
    Function<double,2*NDIM> f2=FunctionFactory<double,2*NDIM>(world).functor(slater);

    if (world.rank()==0) print_header2("1 function in the universe");
    MemoryMeasurer mm;
    mm.measure_and_print(world);
    double total_memory=mm.total_memory(world);
    if (world.rank()==0) print("total memory in universe",total_memory);


    // create functions in all worlds
    {
        if (world.rank()==0) print_header2("1 function per subworld");
        std::shared_ptr<World> subworld=create_worlds(world,world.size());

        {
            // Function<double,2*NDIM> g2_universe=FunctionFactory<double,2*NDIM>(world).functor(slater);
            FunctionDefaults<2*NDIM>::set_default_pmap(*subworld);
            Function<double,2*NDIM> g2=FunctionFactory<double,2*NDIM>(*subworld).functor(slater);

            print("\n---\n");
            MemoryMeasurer mm1;
            mm1.search_world(*subworld);
            mm1.print_memory_map(*subworld,"subworld"+std::to_string(subworld->id()));

            if (world.rank()==0) print("\n---\n");
            MemoryMeasurer mm2;
            mm2.search_all_worlds();
            mm2.print_memory_map(world,"all worlds");
            FunctionDefaults<2*NDIM>::set_default_pmap(world);

            // print success
            double total_memory1=mm2.total_memory(world);
            double total_memory_ref=world.size()*total_memory+total_memory;
            if (world.rank()==0) {
                print("total memory in universe",total_memory1);
                print("should be (nsubworld+1)*total_mem",total_memory_ref);
                print("difference",total_memory1-total_memory_ref);
            }
        }
        subworld->gop.fence();
    }


    return 0;
}

int main(int argc, char** argv) {
    madness::initialize(argc, argv);

    madness::World world(SafeMPI::COMM_WORLD);
    world.gop.fence();
    startup(world,argc,argv);
    const int k=7;
    const double thresh=1.e-5;
    const double L=24.0;
    FunctionDefaults<1>::set_cubic_cell(-L,L);
    FunctionDefaults<2>::set_cubic_cell(-L,L);
    FunctionDefaults<3>::set_cubic_cell(-L,L);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<2>::set_thresh(thresh);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<2>::set_k(k);
    FunctionDefaults<3>::set_k(k);
    int result=0;

    result+=test_size<2>(world);

    print("result",result);
    madness::finalize();
    return result;

}
