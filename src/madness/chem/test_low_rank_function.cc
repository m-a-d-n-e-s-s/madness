//
// Created by Florian Bischoff on 4/24/23.
//


#include<madness.h>
#include<madness/chem/lowrankfunction.h>
#include<madness/chem/electronic_correlation_factor.h>


#include<madness/world/test_utilities.h>
#include<madness/world/timing_utilities.h>
#include <random>


using namespace madness;

struct LowRankFunctionParameters : QCCalculationParametersBase {

    LowRankFunctionParameters() : QCCalculationParametersBase() {

        // initialize with: key, value, comment (optional), allowed values (optional)
        initialize<double>("radius",2.0,"the radius");
        initialize<double>("volume_element",0.1,"volume covered by each grid point");
        initialize<long>("rank",500,"the number of grid points in random grids");
        initialize<bool>("stable_power_iteration",true,"use stable power iteration algorithm (orthonormalize)");
        initialize<double>("tol",1.e-12,"rank-reduced choleski tolerance");
        initialize<std::string>("gridtype","random","the grid type",{"random","cartesian"});
        initialize<std::string>("rhsfunctiontype","exponential","the type of function",{"delta","exponential"});
        initialize<int>("optimize",1,"number of optimization iterations");
    }

    void read_and_set_derived_values(World& world, const commandlineparser& parser, std::string tag) {
        read_input_and_commandline_options(world,parser,tag);
    }

    double radius() const {return get<double>("radius");}
    double volume_element() const {return get<double>("volume_element");}
    double tol() const {return get<double>("tol");}
    long rank() const {return get<long>("rank");}
    bool stable_power_iteration() const {return get<bool>("stable_power_iteration");}
    int optimize() const {return get<int>("optimize");}
    std::string gridtype() const {return get<std::string>("gridtype");}
    std::string rhsfunctiontype() const {return get<std::string>("rhsfunctiontype");}
};



int test_lowrank_function(World& world, LowRankFunctionParameters& parameters) {
    test_output t1("CCPairFunction::low rank function");
    t1.set_cout_to_terminal();
    madness::default_random_generator.setstate(int(cpu_time())%4149);
    madness::default_random_generator.setstate(int(cpu_time())%4149);

    constexpr std::size_t LDIM=3;
    constexpr std::size_t NDIM=2*LDIM;
    print("eps, k, NDIM",FunctionDefaults<NDIM>::get_thresh(),FunctionDefaults<NDIM>::get_k(),NDIM);

    parameters.print("grid","end");

    Vector<double,LDIM> offset;
    offset.fill(0.0);
//    Function<double,LDIM> phi1=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r)
//            { return exp(-r.normf());});
    Function<double,LDIM> phi1=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r)
                                                                           { return 1.0;});
    Function<double,LDIM> phi2=FunctionFactory<double,LDIM>(world).functor([&offset](const Vector<double,LDIM>& r)
                                                                           { return exp(-1.0*(r-offset).normf());});
    Function<double,LDIM> one=FunctionFactory<double,LDIM>(world)
            .functor([](const Vector<double,LDIM>& r) { return 1.0;});

    std::shared_ptr<real_convolution_3d> f12(SlaterOperatorPtr(world,1.0,1.e-6,FunctionDefaults<LDIM>::get_thresh()));

//    auto f = [](const coord_6d& r) {
//        coord_3d r1={r[0],r[1],r[2]};
//        coord_3d r2={r[3],r[4],r[5]};
//        return exp(-(r1-r2).normf() -r2.normf());
//    };

    auto compute_result = [&world, &one](const auto& lrf) {
        real_function_3d result=real_factory_3d(world);
        for (int r=0; r<lrf.rank(); ++r) result+=lrf.g[r]*inner(one,lrf.h[r]);
        return result;
    };
    auto compute_relative_error = [&world,&parameters](const auto reference, const auto result, const auto lrf) {
        auto diff=reference-result;
        double refnorm=reference.norm2();
        double resultnorm=result.norm2();
        double error=diff.norm2();
        return error/refnorm;
    };

    auto output=[&parameters] (const double projection_error, const long projection_rank, const double projection_time,
                               const double optimized_error, const long optimized_rank, const double optimized_time) {
        print("error",parameters.radius(),parameters.gridtype(),parameters.rhsfunctiontype(),parameters.volume_element(),
              parameters.tol(),
              projection_error,projection_rank,projection_time,
              optimized_error,optimized_rank,optimized_time);
    };

    // \phi(1) \bar \phi(1) = \intn phi(1) \phi(2) f(1,2) d2
    auto reference = phi1* (*f12)(phi2);
    output(0.0,0.0,0.0,0.0,0.0,0.0);

    LowRank<double,6> lrf(f12,copy(phi1),copy(phi2));
    lrf.stable_power_iteration=parameters.stable_power_iteration();
//    plot_plane<6>(world,lrf.lrfunctor,"plot_f12_r2",PlotParameters(world).set_plane({"x1","x4"}));
//    lrf.project(parameters.rank(),parameters.radius(),parameters.gridtype(),parameters.rhsfunctiontype());
    double cpu0=cpu_time();
    lrf.project(parameters.volume_element(),parameters.radius(),parameters.gridtype(),parameters.rhsfunctiontype(),parameters.tol());
    double cpu1=cpu_time();

    // compare
    // \phi(1) \bar \phi(1) = \int phi(1) \phi(2) f(1,2) d2
    //       = \int \sum_r g_r(1) h_r(2)  d2
    //       = \sum_r g_r(1) <\phi|h_r>
    real_function_3d result=compute_result(lrf);
    double projection_error=compute_relative_error(reference,result,lrf);
    output(projection_error,lrf.rank(),cpu1-cpu0,0.0,0.0,0.0);

    double cpu2=cpu_time();
    lrf.optimize(parameters.optimize());
    double cpu3=cpu_time();
    result=compute_result(lrf);
    double optimization_error=compute_relative_error(reference,result,lrf);
    output(projection_error,lrf.rank(),cpu1-cpu0,optimization_error,lrf.rank(),cpu3-cpu2);
    bool success=(projection_error<5.e-2) and (optimization_error<1.e-2);


    return t1.end(success);

}

int main(int argc, char **argv) {

    madness::World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);
    int k = parser.key_exists("k") ? std::atoi(parser.value("k").c_str()) : 6;
    double thresh  = parser.key_exists("thresh") ? std::stod(parser.value("thresh")) : 1.e-4;
    FunctionDefaults<6>::set_tensor_type(TT_2D);

    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<2>::set_thresh(thresh);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<4>::set_thresh(thresh);
    FunctionDefaults<5>::set_thresh(thresh);
    FunctionDefaults<6>::set_thresh(thresh);

    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<2>::set_k(k);
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<4>::set_k(k);
    FunctionDefaults<5>::set_k(k);
    FunctionDefaults<6>::set_k(k);

    FunctionDefaults<1>::set_cubic_cell(-10.,10.);
    FunctionDefaults<2>::set_cubic_cell(-10.,10.);
    FunctionDefaults<3>::set_cubic_cell(-10.,10.);
    FunctionDefaults<4>::set_cubic_cell(-10.,10.);
    FunctionDefaults<5>::set_cubic_cell(-10.,10.);
    FunctionDefaults<6>::set_cubic_cell(-10.,10.);

    FunctionDefaults<2>::set_tensor_type(TT_FULL);
    print("numerical parameters: k, eps(3D), eps(6D)", FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(),
          FunctionDefaults<6>::get_thresh());
    LowRankFunctionParameters parameters;
    parameters.read_and_set_derived_values(world,parser,"grid");
    int isuccess=0;
#ifdef USE_GENTENSOR

    try {
        parser.set_keyval("geometry", "he");
        parser.print_map();

        isuccess+=test_lowrank_function(world,parameters);
    } catch (std::exception& e) {
        madness::print("an error occured");
        madness::print(e.what());
    }
#else
    print("could not run test_ccpairfunction: U need to compile with ENABLE_GENTENSOR=1");
#endif
    finalize();

    return isuccess;
}