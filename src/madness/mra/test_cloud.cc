/*
 * test_cloud.cc
 *
 *  Created on: Aug 28, 2020
 *      Author: fbischoff
 */

#include<madness/mra/mra.h>
#include<madness/world/cloud.h>
#include<madness/mra/macrotaskq.h>

using namespace madness;


struct gaussian {
    double a;
    gaussian() : a() {};
    gaussian(double aa) : a(aa) {}
    double operator()(const coord_4d& r) const {
        double x=r[0], y=r[1], z=r[2], aa=r[3];
        return exp(-a*(x*x + y*y + z*z * aa*aa));//*abs(sin(abs(2.0*x))) *cos(y);
    }
    double operator()(const coord_3d& r) const {
        double x=r[0], y=r[1], z=r[2];
        return exp(-a*(x*x + y*y + z*z ));//*abs(sin(abs(2.0*x))) *cos(y);
    }
};

template<typename T>
class variantT {
   std::variant<Tensor<T>> v;
};

template <typename> struct is_vector_of_world_objects: std::false_type {};
template <typename T> struct is_vector_of_world_objects<std::vector<std::is_constructible<T, World&>>> : std::true_type {};

int main(int argc, char** argv) {

    madness::World& universe = madness::initialize(argc, argv);
    startup(universe,argc,argv);
    int success=0;
    {
        Cloud cloud(universe);
        cloud.set_debug(true);
        real_function_3d f = real_factory_3d(universe).functor(gaussian(1.0));
        double n1=f.norm2();
        auto records = cloud.store(universe, f);
        int i=3;
        auto irecords = cloud.store(universe, i);
        long l=4;
        auto lrecords = cloud.store(universe, l);

        std::vector<Function<double,3>> v;
        v.push_back(real_factory_3d(universe).functor(gaussian(2.0)));
        v.push_back(real_factory_3d(universe).functor(gaussian(3.0)));
        auto vrecords = cloud.store(universe, v);

//        print("vector of world objects",is_madness_function_vector<std::vector<Function<double,3>>>::value);


        auto subworld_ptr = MacroTaskQ::create_worlds(universe, universe.size());
        World& subworld=*subworld_ptr;
        print("my world: universe_rank, subworld_id", universe.rank(), subworld.id());

        {
            MacroTaskQ::set_pmap(subworld);
            auto f2=cloud.load1<real_function_3d>(subworld,records);
            double n2=f2.norm2();
            print("sizeof(real_function_3d)",sizeof(real_function_3d));

            auto f3=cloud.load1<real_function_3d>(subworld,records);
            double n3=f3.norm2();
            auto f4=cloud.load1<real_function_3d>(subworld,records);
            double n4=f4.norm2();
            print("n1/n2/n3/n4",n1,n2,n3,n4);


            auto i3=cloud.load1<int>(subworld,irecords);
            auto i4=cloud.load1<int>(subworld,irecords);
            auto i5=cloud.load1<int>(subworld,irecords);
            print("i1/i3/i4/i5",i,i3,i4,i5);

            auto l3=cloud.load1<long>(subworld,lrecords);
            auto l4=cloud.load1<long>(subworld,lrecords);
            auto l5=cloud.load1<long>(subworld,lrecords);
            print("l1/l3/l4/l5",l,l3,l4,l5);

            auto v3=cloud.load1<std::vector<Function<double,3>>>(subworld,vrecords);
            auto v4=cloud.load1<std::vector<Function<double,3>>>(subworld,vrecords);
            auto v5=cloud.load1<std::vector<Function<double,3>>>(subworld,vrecords);

            MacroTaskQ::set_pmap(universe);

            cloud.clear_cache();
        }
        subworld.gop.fence();
    }
    universe.gop.fence();
    madness::finalize();

    return success;
}

template <> volatile std::list<detail::PendingMsg> WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending = std::list<detail::PendingMsg>();
template <> Spinlock WorldObject<WorldContainerImpl<long, std::vector<unsigned char>, madness::Hash<long> > >::pending_mutex(0);
