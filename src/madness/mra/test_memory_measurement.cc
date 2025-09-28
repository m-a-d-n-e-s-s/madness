//
// Created by Florian Bischoff on 1/9/25.
//

#include<madness/mra/memory_measurement.h>
#include<madness/world/test_utilities.h>
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

    test_output t1("test_size");
    t1.set_cout_to_terminal();
    // create a slater function, slightly offset to create an uneven distribution
    auto slater=[](const Vector<double,2*NDIM>& r){return exp(-(r-0.1).normf());};
    Function<double,2*NDIM> f2=FunctionFactory<double,2*NDIM>(world).functor(slater);

    if (world.rank()==0) print_header2("1 function in the universe");
    MemoryMeasurer::measure_and_print(world);


    // create functions in all worlds
    {
        if (world.rank()==0) print_header2("1 function per subworld");
        std::shared_ptr<World> subworld=create_worlds(world,world.size());

        {
            // Function<double,2*NDIM> g2_universe=FunctionFactory<double,2*NDIM>(world).functor(slater);
            FunctionDefaults<2*NDIM>::set_default_pmap(*subworld);
            Function<double,2*NDIM> g2=FunctionFactory<double,2*NDIM>(*subworld).functor(slater);
            MemoryMeasurer::measure_and_print(world);
            FunctionDefaults<2*NDIM>::set_default_pmap(world);
        }
        subworld->gop.fence();
    }


    return t1.end();
}

/// test functions to get a map of ranks to hosts
/// and to get the lowest rank on each host
int test_host_rank_map(World& world) {
    test_output t1("testing host-rank map");
    t1.set_cout_to_terminal();
    // print a list: hostname, rank, lowest rank
    if (world.rank()==0) print_header2("ranks per host");
    auto ranks_per_host1=ranks_per_host(world);
    auto rank=world.rank();
    auto rank_host_map=rank_to_host_and_rss_map(world);
    world.gop.broadcast_serializable(rank_host_map,0);
    // ordered printing
    for (int r=0; r<world.size(); r++, world.gop.fence()) {
        if (world.rank()==r) print("rank",rank,"hostname",rank_host_map[rank].first);
    }

    world.gop.fence();
    // print unique ranks per host
    if (world.rank()==0) {
        print("unique ranks per host (the lowest rank on each host):");
        auto primary_ranks=primary_ranks_per_host(world,ranks_per_host1);
        for (const auto& r : primary_ranks) print("rank",r);
    }
    // ordered printing
    for (int r=0; r<world.size(); r++, world.gop.fence()) {
        if (world.rank()==r) print("lowest rank on host of rank",rank,
            lowest_rank_on_host_of_rank(ranks_per_host1,rank));
    }
    world.gop.fence();
    return t1.end();
}


/// replicate a function on all nodes (not ranks)
/// check that all keys are on the expected rank
/// check that the norm of a specific key is the same as in the distributed function
int test_node_replicated_function(World& world) {
    test_output t1("testing node-replicated function");
    t1.set_cout_to_terminal();

    // some information for the user
    if (world.rank()==0) {
        print("world.size()",world.size());
        print("ranks and hosts");
    }

    // which ranks are on which host
    auto rank_host_map=rank_to_host_and_rss_map(world);
    world.gop.broadcast_serializable(rank_host_map,0);
    // ordered printing
    for (int r=0; r<world.size(); r++, world.gop.fence()) {
        if (world.rank()==r) print("rank",r,"hostname",rank_host_map[r].first);
    }

    // a 2d Gaussian -- reasonably small
    real_function_2d f=real_factory_2d(world).f([](const coord_2d& r){return exp(-r.normf());});

    // print out the tree size of the distributed function, global sum over all ranks
    std::size_t sz_distributed=f.tree_size();
    print("tree size of the distributed function",sz_distributed);

    // compute the norm of a specific key
    Key<2> key(4, {8,9});
    double n=f.get_impl()->get_coeffs().find(key).get()->second.coeff().normf();

    // check if all keys of functimpl are on a given rank
    auto check_keys_are_on_rank=[&](const FunctionImpl<double,2>* fimpl, const ProcessID rank) {
        for (const auto& [key,node] : fimpl->get_coeffs()) {
            int owner=fimpl->get_coeffs().owner(key);
            if (owner!=rank) {
                print("key",key,"owner",owner,"expected owner",rank);
                return false;
            }
        }
        return true;
    };

    // replicate f on all ranks
    {
        auto f1=copy(f);
        f1.replicate();
        std::size_t sz_replicated=f1.tree_size();
        t1.checkpoint(sz_replicated==sz_distributed*world.size(),"rank "+std::to_string(world.rank())+
            ": tree size replicated = size distributed * nproc = "+std::to_string(sz_replicated));

        // check that all keys are on the expected rank, ie. the local rank
        for (int r=0; r<world.size(); r++, world.gop.fence()) {
            if (world.rank()==r) {
                print("rank ",r,"hostname",rank_host_map[r].first);
                // f1.get_impl()->do_print_tree(f1.get_impl()->get_cdata().key0,std::cout,4);
                t1.checkpoint(check_keys_are_on_rank(f1.get_impl().get(),r),"all keys on rank "+std::to_string(r));
            }
            world.gop.fence();
        }
    }

    // replicate f on all hosts, i.e. the lowest rank on each host only
    {
        auto f1=copy(f);
        f1.replicate_on_hosts();
        // f1.replicate_on_nodes();
        std::size_t sz_replicated=f1.tree_size();
        auto ranks_per_host1=ranks_per_host(world);
        auto primary_ranks=primary_ranks_per_host(world,ranks_per_host1);

        print("unique ranks per host (the lowest rank on each host):");
        for (const auto& r : primary_ranks) print("rank",r);
        print("tree size",sz_replicated);
        world.gop.fence();

        t1.set_do_print(true);  // print on all ranks for debugging

        // check that all keys are on the expected rank, ie. the corresponding primary rank
        long myowner = lowest_rank_on_host_of_rank(ranks_per_host1, world.rank());
        for (int r=0; r<world.size(); r++, world.gop.fence()) {
            if (world.rank()==r) {
                print("rank ",r,"hostname",rank_host_map[r].first);
                // f1.get_impl()->do_print_tree(f1.get_impl()->get_cdata().key0,std::cout,4);
                t1.checkpoint(check_keys_are_on_rank(f1.get_impl().get(),myowner),"all keys on rank "+std::to_string(myowner));
            }
            world.gop.fence();
        }

        // compute the norm of a specific key, should be the same as in the distributed function
        double n1=f1.get_impl()->get_coeffs().find(key).get()->second.coeff().normf();
        t1.checkpoint(n,n1,1.e-12,"norm of key is the same as in the distributed function");
        t1.checkpoint(sz_replicated==sz_distributed*primary_ranks.size(),"tree size replicated on nodes = size distributed * n_nodes = "
            +std::to_string(sz_replicated));

    }

    t1.set_do_print(world.rank()==0);
    return t1.end();
}

/// test creating 2 distinct MPI groups and communicators
/// then create a world on each communicator and do some work there
int test_mpi_group(World& world) {
    test_output t1("testing MPI groups",world.rank()==0);
    t1.set_cout_to_terminal();

    /// print in rank-order
    auto oprint = [&](World& world, auto &&... args) {
        world.gop.fence();
        for (int r=0; r<world.size(); ++r) {
            if (r==world.rank()) {
                std::cout << "rank " << world.rank() << ": ";
                print(std::forward<decltype(args)>(args)...);
            }
            world.gop.fence();
        }
    };


    // find unique ranks per host
    auto ranks_per_host1=ranks_per_host(world);
    std::vector<int> primary_ranks=primary_ranks_per_host(world,ranks_per_host1);
    world.gop.broadcast_serializable(primary_ranks,0);
    world.gop.fence();
    oprint(world,"primary_ranks: ", primary_ranks);
    // get a list of all other ranks that are not unique
    std::vector<int> secondary_ranks;
    for (int r=0; r<world.size(); ++r) {
        if (std::find(primary_ranks.begin(),primary_ranks.end(),r)==primary_ranks.end())
            secondary_ranks.push_back(r);
    }
    oprint(world,"secondary_ranks: ", secondary_ranks);
    t1.checkpoint(secondary_ranks.size()+primary_ranks.size() == world.size(),"total number of ranks");


    long myowner = lowest_rank_on_host_of_rank(ranks_per_host1, world.rank());
    // check if this rank is in the unique list
    bool i_am_in_unique_list=world.rank()==myowner;
    oprint(world," myowner ",myowner," i_am_in_unique_list ",i_am_in_unique_list);
    // check that myowner is in the list of unique ranks
    t1.checkpoint(std::find(primary_ranks.begin(),primary_ranks.end(),myowner)!=primary_ranks.end(),
        "myowner is in the list of unique ranks");

    // step 2-1: create a world with only the unique ranks and replicate there (see test_world.cc)
    const int color= i_am_in_unique_list ? 1 : 0;
    if (color) {
        SafeMPI::Group g_unique = world.mpi.comm().Get_group().Incl(int(primary_ranks.size()), &primary_ranks[0]);
        SafeMPI::Intracomm comm_unique = world.mpi.comm().Create(g_unique);
        print(world.rank(),"in unique list");
        {
            World world_unique(comm_unique);
            world_unique.gop.fence();
            // get the hostname of all ranks in this unique world
            std::string hostname=madness::get_hostname();
            auto hostnames=world_unique.gop.concat0(std::vector<std::string>(1,hostname));
            if (world_unique.rank()==0) {
                print("hostnames of unique ranks:");
                for (const auto& h : hostnames) print("hostname",h);
            }
            world_unique.gop.fence();

            // check if there are duplicates in hostnames
            std::set<std::string> hostnames_set(hostnames.begin(),hostnames.end());
            t1.checkpoint(hostnames_set.size()==hostnames.size(),"all hostnames in unique world are different");


        }
    } else {
        // need this to avoid deadlock in MPI_Comm_create (why??)
        print(world.rank(),"not in unique list");
        SafeMPI::Group g_nonunique = world.mpi.comm().Get_group().Incl(secondary_ranks.size(), &secondary_ranks[0]);
        SafeMPI::Intracomm comm_nonunique = world.mpi.comm().Create(g_nonunique);
    }
    world.gop.fence();
    return t1.end();
}

int main(int argc, char** argv) {
    madness::World& world=madness::initialize(argc, argv);

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
    result+=test_host_rank_map(world);
    result+=test_mpi_group(world);
    result+=test_node_replicated_function(world);


    print("result",result);
    madness::finalize();
    return result;

}
