//
// Created by Florian Bischoff on 4/24/23.
//


#include<madness.h>
#include<madness/chem/lowrankfunction.h>


#include<madness/world/test_utilities.h>
#include<madness/world/timing_utilities.h>
#include <random>


using namespace madness;



int test_lowrank_function(World& world, LowRankFunctionParameters& parameters) {
    test_output t1("CCPairFunction::low rank function");
    t1.set_cout_to_terminal();
    madness::default_random_generator.setstate(int(cpu_time())%4149);
    madness::default_random_generator.setstate(int(cpu_time())%4149);
    std::string id=unique_fileid();

    constexpr std::size_t LDIM=3;
    constexpr std::size_t NDIM=2*LDIM;
    print("eps, k, NDIM, id",FunctionDefaults<NDIM>::get_thresh(),FunctionDefaults<NDIM>::get_k(),NDIM,id);

    parameters.print("grid","end");
    std::string f12type=parameters.f12type();
    std::string transpose=parameters.get<std::string>("transpose");

    json j;
    std::string jsonfilename="test_low_rank_function."+id+".json";
    j["radius"]=parameters.radius();
    j["f12type"]=parameters.f12type();
    j["gamma"]=parameters.gamma();
    j["volume_element"]=parameters.volume_element();
    j["tol"]=parameters.tol();
    j["hard_shell"]=parameters.hard_shell();
    j["transpose"]=transpose;
    j["orthomethod"]=parameters.orthomethod();
    j["gridtype"]=parameters.gridtype();
    j["rhsfunctiontype"]=parameters.rhsfunctiontype();
    j["optimize"]=parameters.optimize();
    std::ofstream of(jsonfilename,std::ios::out);
    of<<j;
    of.close();

    Vector<double,LDIM> offset;
    offset.fill(0.0);
    Function<double,LDIM> phi1=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r)
                                                                           { return 1.0;});
    Function<double,LDIM> phi2=FunctionFactory<double,LDIM>(world).functor([&offset](const Vector<double,LDIM>& r)
                                                                           { return exp(-1.0*(r-offset).normf());});
    if (transpose=="slater1") std::swap(phi1,phi2);
    {
        double n1 = phi1.norm2();
        double n2 = phi2.norm2();
        bool first_one = (fabs(phi1({1, 1, 1}) - 1.0) < 1.e-6);
        if (world.rank() == 0) {
            if (first_one) print("1(1) phi(2)");
            else print("phi(1) 1(2)");
            print("norms", n1, n2);
        }
    }

    Function<double,LDIM> one=FunctionFactory<double,LDIM>(world)
            .functor([](const Vector<double,LDIM>& r) { return 1.0;});


    std::shared_ptr<real_convolution_3d> f12;
    if (f12type=="slaterf12") f12.reset(SlaterF12OperatorPtr(world,parameters.gamma(),1.e-6,FunctionDefaults<LDIM>::get_thresh()));
    else if (f12type=="slater") f12.reset(SlaterOperatorPtr(world,parameters.gamma(),1.e-6,FunctionDefaults<LDIM>::get_thresh()));
    else {
        MADNESS_EXCEPTION(std::string("unknown f12type"+f12type).c_str(),1);
    }


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
    plot_plane<3>(world,reference,"reference."+id,PlotParameters(world).set_plane({"x1","x2"}));
    double n2=reference.norm2();
    print("reference.norm2() = int f12 phi2 d2",n2);
    output(0.0,0.0,0.0,0.0,0.0,0.0);

    LowRankFunction<double,6> lrf(f12, copy(phi1), copy(phi2));
    plot_plane<6>(world,lrf.lrfunctor,"plot_original."+id,PlotParameters(world).set_plane({"x1","x4"}));
    double cpu0=cpu_time();
    lrf.project(parameters);
    double cpu1=cpu_time();
    double error1=lrf.l2error();
    print("l2error projection",error1);
//    plot_plane<6>(world,lrf,"plot_lrf_projection."+id,PlotParameters(world).set_plane({"x1","x4"}));

    // compare
    // \phi(1) \bar \phi(1) = \int phi(1) \phi(2) f(1,2) d2
    //       = \int \sum_r g_r(1) h_r(2)  d2
    //       = \sum_r g_r(1) <\phi|h_r>
    real_function_3d result=compute_result(lrf);
    double projection_error=compute_relative_error(reference,result,lrf);
    auto diff=reference-result;
//    plot_plane<3>(world,diff,"plot_diff_int_projection."+id,PlotParameters(world).set_plane({"x1","x2"}));
//    plot_plane<3>(world,result,"plot_lrf_int_projection."+id,PlotParameters(world).set_plane({"x1","x2"}));
    output(projection_error,lrf.rank(),cpu1-cpu0,0.0,0.0,0.0);
    j["projection_error_integrated"]=projection_error;
    j["projection_error_l2"]=error1;
    j["projection_time"]=cpu1-cpu0;
    std::ofstream of1(jsonfilename,std::ios::out);
    of1<<j;
    of1.close();

    double cpu2=cpu_time();
    lrf.optimize(parameters.optimize());
    double error2=lrf.l2error();
    print("l2error optimization",error2);
    double cpu3=cpu_time();
    result=compute_result(lrf);
    diff=reference-result;
    plot_plane<3>(world,diff,"plot_diff_int_optimization."+id,PlotParameters(world).set_plane({"x1","x2"}));
    plot_plane<3>(world,result,"plot_lrf_int_optimization."+id,PlotParameters(world).set_plane({"x1","x2"}));
    double optimization_error=compute_relative_error(reference,result,lrf);
    output(projection_error,lrf.rank(),cpu1-cpu0,optimization_error,lrf.rank(),cpu3-cpu2);
    bool success=(projection_error<5.e-2) and (optimization_error<1.e-2);

    j["optimization_error_integrated"]=optimization_error;
    j["optimization_error_l2"]=error2;
    j["optimization_time"]=cpu3-cpu2;

    std::ofstream of2(jsonfilename,std::ios::out);
    of2<<j;
    of2.close();

    return t1.end(success);

}

/// test the K commutator of the He atom

/// < ij | K f12 | ij >
int test_Kcommutator(World& world, LowRankFunctionParameters& parameters) {
    test_output t1("CCPairFunction::low exchange commutator");
    t1.set_cout_to_terminal();
    madness::default_random_generator.setstate(int(cpu_time())%4149);
    std::string id=unique_fileid();

    constexpr std::size_t LDIM=3;
    constexpr std::size_t NDIM=2*LDIM;
    print("eps, k, NDIM, id",FunctionDefaults<NDIM>::get_thresh(),FunctionDefaults<NDIM>::get_k(),NDIM,id);
    parameters.print("grid");

    real_convolution_3d g12=(CoulombOperator(world,1.e-6,FunctionDefaults<LDIM>::get_thresh()));
    std::shared_ptr<real_convolution_3d> f12ptr;
    std::string f12type=parameters.f12type();
    if (f12type=="slaterf12") f12ptr.reset(SlaterF12OperatorPtr(world,parameters.gamma(),1.e-6,FunctionDefaults<LDIM>::get_thresh()));
    else if (f12type=="slater") f12ptr.reset(SlaterOperatorPtr(world,parameters.gamma(),1.e-6,FunctionDefaults<LDIM>::get_thresh()));
    else {
        MADNESS_EXCEPTION(std::string("unknown f12type"+f12type).c_str(),1);
    }
    real_convolution_3d& f12=*f12ptr;

    real_function_3d phi=real_factory_3d(world).f([](const coord_3d& r){return exp(-r.normf());});
    double n=phi.norm2();
    phi.scale(1/n);
    real_function_3d phi_k=phi; // lookys silly, helps reading.


    // reference term ( < ij | K(1) )  = <Ki(1) j(2) |
    real_function_3d Ki=phi_k*g12(phi*phi_k);
    // < ij | K(1) f12 | ij > = < Ki(1) j(2) | f12 | i(1) j(2) > = <i* Ki(1) f12( j*j(2))
    double reference=inner(Ki*phi,f12(phi*phi));
    print("reference <ij | K(1) f12 | ij>",reference);

    json j;
    std::string jsonfilename="test_kcommuntator."+id+".json";
    j["radius"]=parameters.radius();
    j["f12type"]=parameters.f12type();
    j["gamma"]=parameters.gamma();
    j["thresh"]=FunctionDefaults<3>::get_thresh();
    j["volume_element"]=parameters.volume_element();
    j["tol"]=parameters.tol();
    j["hard_shell"]=parameters.hard_shell();
    j["orthomethod"]=parameters.orthomethod();
    j["gridtype"]=parameters.gridtype();
    j["rhsfunctiontype"]=parameters.rhsfunctiontype();
    j["optimize"]=parameters.optimize();
    j["reference"]=reference;

    auto json2file= [](const json& j, const std::string& jsonfilename) {
        std::ofstream of(jsonfilename, std::ios::out);
        of << j;
        of.close();
    };

    json2file(j,jsonfilename);
    timer t(world);
    auto compute_error = [&](const std::string msg, const LowRankFunction<double,6>& lrf) {
        auto gk = mul(world, phi_k, g12(lrf.g * phi_k)); // function of 1
        auto hj = lrf.h * phi; // function of 2
        Tensor<double> j_hj = inner(world, phi, hj);
        Tensor<double> i_gk = inner(world, phi, gk);
        double result_right = j_hj.trace(i_gk);
        print(msg, result_right);
        j[msg]=result_right-reference;
        j[msg+"_rank"]=lrf.rank();
        j[msg+"_compute_time"]=t.tag(msg+"_compute_time");
        json2file(j,jsonfilename);
    };


    {
        // lowrankfunction left phi: lrf(1',2) = f12(1',2) i(1')
        // K f12 ij = \sum_k k(1) \int g(1,1') f12(1'2) i(1') j(2) k(1') d1'
        //          = \sum_kr k(1) j(2) \int g(1,1') g_r(1') h_r(2) k(1') d1'
        //          = \sum_r j(2) h_r(2) \sum_k k(1) \int g(1,1') g_r(1') k(1') d1'
        real_function_3d one = real_factory_3d(world).f([](const coord_3d& r) { return 1.0; });
        LowRankFunction<double, 6> fi_one(f12ptr, copy(phi), copy(one));
        fi_one.project(parameters);
        double l2error=fi_one.l2error();
        print("left_project_l2error",l2error);

        j["left_project_time"]=t.tag("left_project_time");
        json2file(j,jsonfilename);
        compute_error("left_project",fi_one);

        fi_one.optimize();
        l2error=fi_one.l2error();
        print("left_optimize_l2error",l2error);
        j["left_optimize_time"]=t.tag("left_optimize_time");
        json2file(j,jsonfilename);
        compute_error("left_optimize",fi_one);

        fi_one.reorthonormalize();
        j["left_reorthonormalize"]=t.tag("left_reorthonormalize");
        json2file(j,jsonfilename);
        compute_error("left_reorthonormalize",fi_one);
    }

//    // lowrankfunction right phi: lrf(1',2) = f12(1',2) i(1')
//    {
//        real_function_3d one = real_factory_3d(world).f([](const coord_3d &r) { return 1.0; });
//        LowRankFunction<double, 6> fi_one(f12ptr, copy(one), copy(phi));
//        fi_one.project(parameters);
//        std::swap(fi_one.g,fi_one.h);
//        j["right_project_time"]=t.tag("right_project_time");
//        json2file(j,jsonfilename);
//
//
//        {
//            auto gk = mul(world, phi_k, g12(fi_one.g * phi_k)); // function of 1
//            auto hj = fi_one.h * phi; // function of 2
//            Tensor<double> j_hj = inner(world, phi, hj);
//            Tensor<double> i_gk = inner(world, phi, gk);
//            double result_right = j_hj.trace(i_gk);
//            print("result_right, project only", result_right);
//            j["right_project"]=result_right-reference;
//            j["right_project_rank"]=fi_one.rank();
//            j["left_optimize_compute_time"]=t.tag("left_optimize_compute_time");
//            j["right_project_compute_time"]=t.tag("right_project_compute_time");
//        }
//        json2file(j,jsonfilename);
//        std::swap(fi_one.g,fi_one.h);
//        fi_one.optimize();
//        std::swap(fi_one.g,fi_one.h);
//        {
//            auto gk = mul(world, phi_k, g12(fi_one.g * phi_k)); // function of 1
//            auto hj = fi_one.h * phi; // function of 2
//            Tensor<double> j_hj = inner(world, phi, hj);
//            Tensor<double> i_gk = inner(world, phi, gk);
//            double result_right = j_hj.trace(i_gk);
//            print("result_right, optimize", result_right);
//            j["right_optimize"]=result_right-reference;
//            j["right_optimize_rank"]=fi_one.rank();
//            j["right_optimize_compute_time"]=t.tag("right_optimize_compute_time");
//        }
//        json2file(j,jsonfilename);
//
//    }

    return 0;

}

template<std::size_t LDIM>
int test_grids(World& world, LowRankFunctionParameters& parameters) {
    randomgrid<LDIM> g(parameters.volume_element(),parameters.radius());
    g.get_grid();

    return 0;
}

template<std::size_t LDIM>
int test_construction_optimization(World& world, LowRankFunctionParameters parameters) {
    parameters.set_user_defined_value("volume_element",0.05);
    constexpr std::size_t NDIM=2*LDIM;
    test_output t1("LowRankFunction::construction/optimization in dimension "+std::to_string(NDIM));
    t1.set_cout_to_terminal();
    OperatorInfo info(1.0,1.e-6,FunctionDefaults<LDIM>::get_thresh(),OT_SLATER);
    auto slater=std::shared_ptr<SeparatedConvolution<double,LDIM> >(new SeparatedConvolution<double,LDIM>(world,info));
    Function<double,LDIM> one=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r){return exp(-0.2*inner(r,r));});

    LowRankFunction<double,NDIM> lrf(slater,one,one);
    lrf.project(parameters);
    double error=lrf.l2error();
    t1.checkpoint(error<2.e-2,"l2 error in projection "+std::to_string(error));
    print("l2 error project ",error);
    lrf.optimize();
    error=lrf.l2error();
    print("l2 error optimize",error);
    t1.checkpoint(error<1.e-2,"l2 error in optimization "+std::to_string(error));
    return t1.end();
}

int main(int argc, char **argv) {

    madness::World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);
    int k = parser.key_exists("k") ? std::atoi(parser.value("k").c_str()) : 6;
    double thresh  = parser.key_exists("thresh") ? std::stod(parser.value("thresh")) : 1.e-5;
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
    parameters.print("grid");
    int isuccess=0;
#ifdef USE_GENTENSOR


    try {

        isuccess+=test_grids<1>(world,parameters);
        isuccess+=test_grids<2>(world,parameters);
        isuccess+=test_grids<3>(world,parameters);
        isuccess+=test_construction_optimization<1>(world,parameters);
        isuccess+=test_construction_optimization<2>(world,parameters);
//        isuccess+=test_arithmetic<1>(world,parameters);
//        isuccess+=test_arithmetic<2>(world,parameters);

//        isuccess+=test_lowrank_function(world,parameters);
//       isuccess+=test_Kcommutator(world,parameters);
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