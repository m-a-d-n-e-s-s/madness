/*
 * test_cloud.cc
 *
 *  Created on: Aug 28, 2020
 *      Author: fbischoff
 */

#include<madness/mra/mra.h>
#include<madness/world/cloud.h>
#include<madness/mra/macrotaskq.h>
#include<madness/world/test_utilities.h>

using namespace madness;


struct gaussian {
    double a;

    gaussian() : a() {};

    gaussian(double aa) : a(aa) {}

    double operator()(const coord_4d &r) const {
        double x = r[0], y = r[1], z = r[2], aa = r[3];
        return exp(-a * (x * x + y * y + z * z * aa * aa));//*abs(sin(abs(2.0*x))) *cos(y);
    }

    double operator()(const coord_3d &r) const {
        double x = r[0], y = r[1], z = r[2];
        return exp(-a * (x * x + y * y + z * z));//*abs(sin(abs(2.0*x))) *cos(y);
    }
};


/// this class stores different member variables in different records of the cloud
class custom_serialize_tester {
public:
    int i;
    double d;

    custom_serialize_tester() : i(0), d(0.0) {}
    bool operator==(const custom_serialize_tester& other) const {
        return i == other.i && d == other.d;
    }

    /// customized function to store this to the cloud

    /// functions and constant_part can be very large and we want to split them and store them in differenc records
    Recordlist<Cloud::keyT> cloud_store(World& world, Cloud& cloud) const {
        // save bookkeeping stuff in a vector
        std::vector<unsigned char> v;
        archive::VectorOutputArchive arout(v);
        arout &  i;

        Recordlist<Cloud::keyT> records;
        records+=cloud.store(world,v);
        records+=cloud.store(world,d);
        return records;
    }

    void cloud_load(World& world, const Cloud& cloud, Recordlist<Cloud::keyT>& recordlist) {
        std::vector<unsigned char> v=cloud.forward_load<std::vector<unsigned char>>(world,recordlist);
        archive::VectorInputArchive arin(v);
        arin & i;
        d=cloud.forward_load<double>(world,recordlist);
    }


};

template<typename T>
double norm(const T i1) { return fabs(i1); }

template<typename T, std::size_t NDIM>
double norm(const Function<T, NDIM> &f1) { return (f1).norm2(); }

template<typename T>
double norm(const Tensor<T> &t) { return t.normf(); }

template<typename T, std::size_t NDIM>
double norm(const std::vector<Function<T, NDIM>> &f1) { return norm2(f1.front().world(), f1); }

template<typename T>
double norm(const std::vector<T> &v) {
    double result=0.0;
    for (const auto &e : v) result += double(e);
    return result;
}

/// a simple example for how to use the cloud for inter-world communication

/// Notes:
/// - during the subworld section no universe-wide fence must be invoked, including the
///   creation of universe WorldObjects -- they must be constructed before
/// - Certain operation can be performed between objects living in different world (e.g. Function::operator+=)
///    -- store/lood pointers to these universe-wide world objects
/// - when the subworld is destroyed all subworld objects must have been destroyed
/// - subworld objects will be destroyed only at subworld fences
void simple_example(World &universe) {

    // this function lives in the universe
    real_function_3d f_universe = real_factory_3d(universe).functor(gaussian(1.0));

    // create the cloud
    {
        Cloud cloud(universe);

        // store f_universal into the cloud, the return value holds the record to find the function again.
        auto recordlist = cloud.store(universe, f_universe);

        // begin subworld section
        auto subworld_ptr = MacroTaskQ::create_worlds(universe, universe.size());
        World &subworld = *subworld_ptr;

        // from now on there must be no universe-wide fences!
        //
        // scopes are important because of deferred destruction:
        // when the subworld is destroyed all objects living in it must
        // have been destroyed before
        {
            // reset process map to subworlds
            MacroTaskQ::set_pmap(subworld);

            // load f into the worlds by passing in the recordlist
            real_function_3d f_subworld = real_factory_3d(subworld);
            if (universe.rank() == 0) {
                f_subworld = cloud.load<real_function_3d>(subworld, recordlist); // has a subworld fence, that's ok
            }
            double norm = f_subworld.norm2();
            // this will print 0 often and the actual norm once
            print("norm of f in subworld", subworld.id(), ":", norm);

            // end subworld section
            MacroTaskQ::set_pmap(universe);
            cloud.clear_cache(subworld);        // includes subworld fence
        } // f_subworld goes out of scope here
        subworld.gop.fence();   // f_subworld is destroyed here
        universe.gop.fence();
    }   // subworld is destroyed here
}

/// test the cloud with message larger than chunk size set in cloud.replicate()
int chunk_example(World &universe) {
    int test_size = 100;
    std::vector<int> testvec;
    for (int i=0; i<test_size; i++) {
        testvec.push_back(i);
    }

    {
        test_output bla("testing replication");
        Cloud cloud(universe);
        auto recordlist = cloud.store(universe, testvec);
        cloud.replicate(50);

        std::vector<int> cloud_vector = cloud.load<std::vector<int>>(universe, recordlist);
        int sum = 0;
        for (int i = 0; i < testvec.size(); i++){
            sum += std::abs(testvec[i] - cloud_vector[i]);
        }
        bla.end((sum==0));
        return sum;
    }
}


template<typename T> using is_world_constructible = std::is_constructible<T, World &>;


/// test storing and loading a custom WorldObject, used e.g. for the scalar output of a macrotask
int test_custom_worldobject(World& universe, World& subworld, Cloud& cloud) {
    test_output t1("testing custom worldobject");
    t1.set_cout_to_terminal();
    cloud.set_debug(false);
    auto o1 =std::shared_ptr<ScalarResult<>>(new ScalarResult(universe));
    auto o5 =std::shared_ptr<ScalarResult<>>(new ScalarResult(universe));
    *o1=1.2;
    if (universe.rank() == 0) gaxpy(1.0,*o1,2.0,2.8);

    auto adrecords = cloud.store(universe, o1);
    MacroTaskQ::set_pmap(subworld);
    print("world constructible",is_world_constructible<ScalarResult<>>::value);

    cloud.set_force_load_from_cache(false);
    auto o2 = cloud.load<std::shared_ptr<ScalarResult<>>>(subworld, adrecords);
    cloud.set_force_load_from_cache(true);
    auto o3 = cloud.load<std::shared_ptr<ScalarResult<>>>(subworld, adrecords);
    double d1=o1->get_local();
    double d2=o2->get_local();
    double d3=o3->get_local();
    std::cout << "pointer  " << o1->id() << " " << o2->id() << " " << o3->id() <<  " other: " << o5->id() << std::endl;
    std::cout << "numerics (plain)" << d1 << " " << d2 << " " << d3 << std::endl;
    std::cout << "numerics (get)  " << o1->get() << " " << o2->get() << " " << o3->get() << std::endl;
    double error=d1-d2;
    cloud.set_force_load_from_cache(false);
    return t1.end(error < 1.e-10 );
}

int test_custom_serialization(World& universe, Cloud& cloud) {
    test_output t1("testing custom serialization");
    t1.set_cout_to_terminal();
    cloud.set_debug(true);
    custom_serialize_tester cst;
    cst.i=1;
    cst.d=2.0;
    static_assert(Cloud::has_cloud_serialize<custom_serialize_tester>::value,"custom_serialize_tester must have a cloud_serialize method");
    {
        auto records = cloud.store(universe, cst);
        auto cst2=cloud.load<custom_serialize_tester>(universe, records);
        t1.checkpoint(cst==cst2,"custom serialization");
    }

    // test being part of a tuple
    typedef std::tuple<int,double,custom_serialize_tester> tupleT;
    tupleT tuple1=std::make_tuple(1,2.0,cst);
    cloud.clear();
    {
        auto records = cloud.store(universe, tuple1);
        auto tuple2=cloud.load<tupleT>(universe, records);

        t1.checkpoint(tuple1==tuple2,"custom serialization with tuple");
    }

    return t1.end();


}

int main(int argc, char **argv) {

    madness::World &universe = madness::initialize(argc, argv);
    startup(universe, argc, argv);

    chunk_example(universe);
    simple_example(universe);
    int success = 0;
    {
        Cloud cloud(universe);
//        cloud.set_debug(true);

        auto subworld_ptr = MacroTaskQ::create_worlds(universe, universe.size());
        World& subworld = *subworld_ptr;

        // test storing custom WorldObject
        success += test_custom_worldobject(universe, subworld, cloud);
        success += test_custom_serialization(universe, cloud);

        if (universe.rank() == 0) print("entering test_cloud");
        print("my world: universe_rank, subworld_id", universe.rank(), subworld.id());

        auto dotest = [&](auto &arg) {

            test_output test_p("testing cloud/shared_ptr<Function> in world "
                               + std::to_string(subworld.id()) + " " + typeid(std::get<0>(arg)).name());
            MacroTaskQ::set_pmap(subworld);

            typedef std::remove_reference_t<decltype(std::get<0>(arg))> argT;
            auto records = std::get<1>(arg);
            double universe_norm = std::get<2>(arg);

            // the first time we load from the cloud's distributed container
            auto copy_of_arg = cloud.load<argT>(subworld, records);
            double n = norm(copy_of_arg);
            double error = n - universe_norm;
            test_p.logger << "error(container)" << error << std::endl;
            if (error > 1.e-10) success++;

            // the second time we load from the cloud's world-local cache
            cloud.set_force_load_from_cache(true);
            auto cached_copy_of_arg = cloud.load<argT>(subworld, records);
            double n_cached = norm(cached_copy_of_arg);
            double error_cached = n_cached - universe_norm;
            test_p.logger << "error(cache)    " << error_cached << std::endl;
            success += test_p.end(error_cached < 1.e-10 && error < 1.e-10);
            cloud.set_force_load_from_cache(false);
            subworld.gop.fence();
        };
        auto tester = [&](auto &&... args) { (dotest(args), ...); };

        // test some standard objects
        real_function_3d f1 = real_factory_3d(universe).functor(gaussian(1.0));
        real_function_3d f2 = real_factory_3d(universe).functor(gaussian(2.0));
        real_function_3d f3 = real_factory_3d(universe).functor(gaussian(3.0));
        int i = 3;
        long l = 4l;
        Tensor<double> t(3, 3);
        t.fillrandom();
        std::vector<real_function_3d> vf{f2, f3};
        std::vector<double> vd{2.0, 3.0};


        auto ipair = std::make_tuple(i, cloud.store(universe, i), norm(i));
        auto lpair = std::make_tuple(l, cloud.store(universe, l), norm(l));
        auto fpair = std::make_tuple(f1, cloud.store(universe, f1), norm(f1));
        auto vpair = std::make_tuple(vf, cloud.store(universe, vf), norm(vf));
        auto tpair = std::make_tuple(t, cloud.store(universe, t), norm(t));
        auto vdpair = std::make_tuple(vd, cloud.store(universe, vd), norm(vd));

        tester(ipair, lpair, fpair, vpair, tpair, vdpair);
        universe.gop.fence();

        MacroTaskQ::set_pmap(universe);
        universe.gop.fence();
        universe.gop.fence();

        // test pointer to FunctionImpl
        typedef std::shared_ptr<Function<double, 3>::implT> impl_ptrT;
        Function<double, 3> ff = real_factory_3d(universe).functor(gaussian(1.5));
        impl_ptrT p1 = ff.get_impl();
        auto precords = cloud.store(universe, p1);

        {
            test_output test_ptr("testing cloud/shared_ptr<Function> in world " + std::to_string(subworld.id()));
            MacroTaskQ::set_pmap(subworld);

            auto p3 = cloud.load<impl_ptrT>(subworld, precords);
            auto p4 = cloud.load<impl_ptrT>(subworld, precords);
            auto p5 = cloud.load<impl_ptrT>(subworld, precords);
            test_ptr.logger << "p1/p2/p3/p4" << " " << p1.get() << " " << p3.get() << " " << p4.get() << " " << p5.get()
                            << std::endl;
            test_ptr.end(p1 == p3 && p1 == p4 && p1 == p5
                         && p1->get_world().id() == p3->get_world().id()
                         && p1->get_world().id() == p4->get_world().id()
                         && p1->get_world().id() == p5->get_world().id());
            Function<double, 3> fff;
            fff.set_impl(p3);
            Function<double, 3> ffsub = real_factory_3d(subworld).functor(gaussian(1.5));
            fff -= ffsub * (1.0 / universe.size());
            MacroTaskQ::set_pmap(universe);
            cloud.clear_cache(subworld);
        }
        subworld.gop.fence();
        universe.gop.fence();
        test_output test_ptr("testing cloud/shared_ptr<Function> numerics in universe");
        double ffnorm = ff.norm2();
        test_ptr.end((ffnorm < 1.e-10));
        universe.gop.fence();


        // test storing tuple
        test_output test_tuple("testing tuple");
        cloud.set_debug(false);
        typedef std::tuple<double, int, Function<double, 3>, impl_ptrT> tupleT;
        tupleT t1{1.0, 2, f1, f2.get_impl()};
        std::vector<double> norm1{1.0, 2.0, f1.norm2()};
        auto turecords = cloud.store(universe, t1);
        {
            MacroTaskQ::set_pmap(subworld);

            cloud.set_force_load_from_cache(false);
            auto t2 = cloud.load<tupleT>(subworld, turecords);
            cloud.set_force_load_from_cache(true);
            auto t3 = cloud.load<tupleT>(subworld, turecords);
            std::vector<double> norm2{1.0, 2.0, std::get<2>(t2).norm2()};
            test_tuple.logger << "error double, int, Function " << norm1[0] - norm2[0] << "  "
                              << norm1[1] - norm2[1] << " " << norm1[2] - norm2[2];
            std::vector<double> norm3{1.0, 2.0, std::get<2>(t3).norm2()};
            test_tuple.logger << "error double, int, Function " << norm1[0] - norm3[0] << " "
                              << norm1[1] - norm3[1] << " " << norm1[2] - norm3[2];
            double error = std::max({norm1[0] - norm2[0], norm1[1] - norm2[1], norm1[2] - norm2[2], norm1[0] - norm3[0],
                                     norm1[1] - norm3[1], norm1[2] - norm3[2]});
            double error1 = std::min(
                    {norm1[0] - norm2[0], norm1[1] - norm2[1], norm1[2] - norm2[2], norm1[0] - norm3[0],
                     norm1[1] - norm3[1], norm1[2] - norm3[2]});
            success += test_tuple.end(error < 1.e-10 && error1 > -1.e-10);
            cloud.set_force_load_from_cache(false);
        }

        // test storing twice (using cache)
        {
            cloud.clear_timings();
            cloud.store(universe, vd);
            auto recordlist = cloud.store(universe, vd);
            auto vd1 = cloud.load<std::vector<double>>(universe, recordlist);
            vd1 = cloud.load<std::vector<double>>(universe, recordlist);
            cloud.print_timings(universe);
            cloud.clear_cache(subworld);
        }
    }
    universe.gop.fence();
    madness::finalize();

    return success;
}

template<> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template<> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(
        0);
