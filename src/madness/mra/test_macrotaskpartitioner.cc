//
// Created by Florian Bischoff on 5/17/21.
//

#include<madness/mra/mra.h>
#include<madness/mra/macrotaskpartitioner.h>
#include<madness/world/test_utilities.h>

using namespace madness;

int test_batch_1D(World& world) {
    Batch_1D b;
    Batch_1D b1(0,-1);
    Batch_1D b2(2,8);
    Batch_1D b3(2,5);

    MADNESS_CHECK(b.is_full_size());
    MADNESS_CHECK(b1.is_full_size());
    MADNESS_CHECK(not b2.is_full_size());
    MADNESS_CHECK(b2.size()==6);

    vector_real_function_1d vf= zero_functions<double,1>(world,10);
    vector_real_function_1d vf3= zero_functions<double,1>(world,3);
    auto tuple=std::make_tuple(vf,3.0,vf);

    auto tuple1=b2.copy_batch<decltype(tuple),0>(tuple);
    MADNESS_CHECK(std::get<0>(tuple1).size()==6);
    auto tuple3=b3.copy_batch<decltype(tuple),1>(tuple);
    MADNESS_CHECK(std::get<2>(tuple3).size()==3);
    MADNESS_CHECK(std::get<2>(tuple3)[0].get_impl() == std::get<2>(tuple)[2].get_impl());

    vf=b3.insert_batch(vf,vf3);
    MADNESS_CHECK(vf[2].get_impl()==vf3[0].get_impl());
    MADNESS_CHECK(vf[3].get_impl()==vf3[1].get_impl());

    print("b ",b);
    print("b3",b3);
    return 0;
}

int test_batch(World& world) {

    Batch batch0;

    Batch_1D input1(0,1), input2(3,8), result(1,4 );
    Batch batch1(input1, input2, result);
    Batch batch2(input1, result);

    print("batch0",batch0);
    print("batch1",batch1);
    print("batch2",batch2);


    vector_real_function_1d vf= zero_functions<double,1>(world,10);
    vector_real_function_1d vf3= zero_functions<double,1>(world,3);

    {
        auto tuple2=std::make_tuple(vf,3.0,vf);
        auto tuple2_parts = batch1.copy_input_batch(tuple2);
        auto result2_parts = batch1.insert_result_batch(vf, vf3);
        MADNESS_CHECK(std::get<0>(tuple2_parts).size() == 1);
        MADNESS_CHECK(result2_parts.size() == 10);
        MADNESS_CHECK(result2_parts[0].get_impl() == vf[0].get_impl());
        MADNESS_CHECK(result2_parts[1].get_impl() == vf3[0].get_impl());
        MADNESS_CHECK(result2_parts[2].get_impl() == vf3[1].get_impl());
        MADNESS_CHECK(result2_parts[3].get_impl() == vf3[2].get_impl());
        MADNESS_CHECK(result2_parts[4].get_impl() == vf[4].get_impl());
    }

    {
        auto tuple2=std::make_tuple(vf,3.0);
        auto tuple2_parts = batch2.copy_input_batch(tuple2);
        auto result2_parts = batch2.insert_result_batch(vf, vf3);
        MADNESS_CHECK(std::get<0>(tuple2_parts).size() == 1);
        MADNESS_CHECK(result2_parts.size() == 10);
        MADNESS_CHECK(result2_parts[0].get_impl() == vf[0].get_impl());
        MADNESS_CHECK(result2_parts[1].get_impl() == vf3[0].get_impl());
        MADNESS_CHECK(result2_parts[2].get_impl() == vf3[1].get_impl());
        MADNESS_CHECK(result2_parts[3].get_impl() == vf3[2].get_impl());
        MADNESS_CHECK(result2_parts[4].get_impl() == vf[4].get_impl());
    }

    return 0;
}

int test_partitioner(World& world) {

    using partitionT  = MacroTaskPartitioner::partitionT;

    vector_real_function_1d vf= zero_functions<double,1>(world,40);
    auto tuple2=std::make_tuple(vf,3.0);
    auto tuple3=std::make_tuple(vf,3.0,vf);

    MacroTaskPartitioner mtp;

    {
        print("\n1D partitioning");
        partitionT partition = mtp.partition_tasks(tuple3);
        for (auto p : partition) print(p);
    }

    {
        mtp.set_nsubworld(2);
        print("\n1D partitioning, nsubworld=2");
        partitionT partition = mtp.partition_tasks(tuple3);
        for (auto p : partition) print(p);
    }

    {
        mtp.set_nsubworld(2);
        mtp.set_dimension(2);
        print("\n2D partitioning, nsubworld=2");
        partitionT partition = mtp.partition_tasks(tuple3);
        for (auto p : partition) print(p);
    }

    {
        mtp.set_nsubworld(2);
        mtp.set_dimension(1);
        print("\n2D partitioning, nsubworld=2");
        partitionT partition = mtp.partition_tasks(tuple2);
        for (auto p : partition) print(p);
    }

    return 0;
}

int main(int argc, char **argv) {

    madness::World &universe = madness::initialize(argc, argv);
    startup(universe, argc, argv);
    int success=1;

    success+=test_batch_1D(universe);
    success+=test_batch(universe);
    success+=test_partitioner(universe);
//    success+=test_partitioner(universe);
//    success+=test_partitioner(universe);


    universe.gop.fence();
    madness::finalize();
    return success;
}



