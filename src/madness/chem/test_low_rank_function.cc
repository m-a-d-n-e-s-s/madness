//
// Created by Florian Bischoff on 4/24/23.
//


#include<madness.h>
#include<madness/chem/lowrankfunction.h>


#include<madness/world/test_utilities.h>
#include<madness/world/timing_utilities.h>
#include <random>


using namespace madness;



int test_lowrank_function(World& world, LowRankFunctionParameters parameters) {
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
    auto compute_relative_error = [](const auto reference, const auto result, const auto lrf) {
        auto diff=reference-result;
        double refnorm=reference.norm2();
        // double resultnorm=result.norm2();
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
//    plot_plane<3>(world,reference,"reference."+id,PlotParameters(world).set_plane({"x1","x2"}));
    double n2=reference.norm2();
    print("reference.norm2() = int f12 phi2 d2",n2);
    output(0.0,0.0,0.0,0.0,0.0,0.0);

    LRFunctorF12<double,6> lrfunctor(f12,phi1,phi1);
    double cpu0=cpu_time();
    auto lrf=LowRankFunctionFactory<double,6>(parameters).project(lrfunctor);
    lrf.do_print=true;
//    plot_plane<6>(world,lrfunctor,"plot_original."+id,PlotParameters(world).set_plane({"x1","x4"}));
    double cpu1=cpu_time();
    double error1=lrf.l2error(lrfunctor);
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
    lrf.optimize(lrfunctor,parameters.optimize());
    double error2=lrf.l2error(lrfunctor);
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
    g12.particle()=1;
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
    real_function_3d phi_k=phi; // looks silly, helps reading.


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
        print(msg,"norm ", result_right);
        print(msg,"error", result_right-reference);
        print(msg,"rank ", lrf.rank());
        j[msg]=result_right-reference;
        j[msg+"_rank"]=lrf.rank();
        j[msg+"_compute_time"]=t.tag(msg+"_compute_time");
        json2file(j,jsonfilename);
    };


    if (1) {
        // lowrankfunction left phi: lrf(1',2) = f12(1',2) i(1')
        // K f12 ij = \sum_k k(1) \int g(1,1') f12(1'2) i(1') j(2) k(1') d1'
        //          = \sum_kr k(1) j(2) \int g(1,1') g_r(1') h_r(2) k(1') d1'
        //          = \sum_r j(2) h_r(2) \sum_k k(1) \int g(1,1') g_r(1') k(1') d1'
        real_function_3d one = real_factory_3d(world).f([](const coord_3d& r) { return 1.0; });
        LRFunctorF12<double,6> lrfunctor(f12ptr,phi,one);
//        LowRankFunction<double, 6> fi_one(f12ptr, copy(phi), copy(one));
        auto fi_one=LowRankFunctionFactory<double,6>(parameters).project(lrfunctor);
//        fi_one.project(parameters);
        double l2error=fi_one.l2error(lrfunctor);
        print("left_project_l2error",l2error);

        j["left_project_time"]=t.tag("left_project_time");
        json2file(j,jsonfilename);
        compute_error("left_project",fi_one);

        fi_one.optimize(lrfunctor);
        l2error=fi_one.l2error(lrfunctor);
        print("left_optimize_l2error",l2error);
        j["left_optimize_time"]=t.tag("left_optimize_time");
        json2file(j,jsonfilename);
        compute_error("left_optimize",fi_one);

//        fi_one.reorthonormalize();
//        j["left_reorthonormalize"]=t.tag("left_reorthonormalize");
//        json2file(j,jsonfilename);
//        compute_error("left_reorthonormalize",fi_one);

        LowRankFunction<double,6> phi0(world);
        phi0.g={phi};
        phi0.h={phi};

        timer t2(world);
        // this is f12|ij>
        auto f12ij=copy(fi_one);
        f12ij.h=f12ij.h*phi;
        double result1=inner(phi0,f12ij);
        print("<ij | f12 | ij>",result1);
        t2.tag("multiply 1(1)* phi(1))");

        // this is f12|(ki) j>
        auto f12kij=copy(f12ij);
        f12kij.g=f12kij.g*phi;
        double result3=inner(phi0,f12kij);
        print("<ij | f12 | (ki) j>",result3);
        t2.tag("multiply 1");

        // this is g(f12|(ki) j>);
        auto gf12kij=copy(f12kij);
        gf12kij.g=g12(f12kij.g);
        double result2=inner(phi0,gf12kij);
        print("<ij | g(f12 | (ki) j>)",result2);
        t2.tag("apply g ");

        // this is kg(f12|(ki) j>);
        auto kgf12kij=copy(gf12kij);
        kgf12kij.g=gf12kij.g * phi;
        double result4=inner(phi0,kgf12kij);
        print("<ij | k g(f12 | (ki) j>)",result4);
        t2.tag("apply g ");
    }


    // apply exchange operator in 6d
//    if (f12type=="slaterf12") {
    if (1) {
//        FunctionDefaults<3>::print();
//        FunctionDefaults<6>::print();
        real_function_6d phi0=CompositeFactory<double,6,3>(world).particle1(phi).particle2(phi);

        double thresh=FunctionDefaults<3>::get_thresh();
        double dcut=1.e-6;
        real_function_6d tmp=TwoElectronFactory(world).dcut(dcut).gamma(parameters.gamma()).f12().thresh(thresh);
        real_function_6d f12ij=CompositeFactory<double,6,3>(world).g12(tmp).particle1(copy(phi)).particle2(copy(phi));

        timer t(world);
        f12ij.fill_tree();
        t.tag("exchange: fill_tree");
        f12ij.print_size("f12ij");

        auto result1=madness::inner(phi0,f12ij);
        print("<ij | f12 | ij>", result1);
        double reference1=inner(phi*phi,f12(phi*phi));
        print("reference <ij |f12 | ij>",reference1);


        real_function_6d kf12ij=multiply(f12ij,copy(phi_k),1);
        kf12ij.print_size("kf12ij");
        t.tag("exchange: multiply 1");

        auto result2=madness::inner(phi0,kf12ij);
        print("<ij | f12 | (ki) j>", result2);
        double reference2=inner(phi*phi,f12(phi*phi*phi));
        print("reference <ij |f12 | (ki) j>",reference2);

        kf12ij.change_tree_state(reconstructed);
        real_function_6d gkf12ij=g12(kf12ij).truncate();
        gkf12ij.print_size("gkf12ij");
        t.tag("exchange: apply g");

        auto result3=madness::inner(phi0,gkf12ij);
        print("<ij | g(1'1) f12 | (ki) j>", result3);
        double reference3=inner(phi*phi,f12(phi*phi*g12(copy(phi))));
        print("reference <ij | g(1'1) f12 | (ki) j>",reference3);


        auto exf12ij=multiply(gkf12ij,copy(phi_k),1).truncate();
        exf12ij.print_size("exf12ij");
        t.tag("exchange: multiply 2");

        auto result=madness::inner(phi0,exf12ij);
        print("<ij | K1 f12 | ij>", result);
        print("error",result-reference);


    }

    return t1.end();

}

template<std::size_t LDIM>
int test_full_rank_functor(World& world, LowRankFunctionParameters parameters) {

    test_output t1("test_full_rank_functor");
    t1.set_cout_to_terminal();
    print_header2("entering test_full_rank_functor");
    constexpr int NDIM=2*LDIM;
//    FunctionDefaults<LDIM>::set_thresh(1.e-6);
//    FunctionDefaults<NDIM>::set_thresh(1.e-6);
    double tol=6.e-3;
    double gaussexponent=2.0;

    // const particle<LDIM> p1=particle<LDIM>::particle1();
    // const particle<LDIM> p2=particle<LDIM>::particle2();

    LRFunctorPure<double,NDIM> functorpure;
    Function<double,NDIM> gauss=FunctionFactory<double,NDIM>(world)
            .functor([&gaussexponent](const Vector<double,NDIM>& r){
                Vector<double,LDIM> a,b;
                for (int i=0; i<LDIM; ++i) {
                    a[i]=r[i];
                    b[i]=r[i+LDIM];
                }
                return exp(-gaussexponent*inner(a-b,a-b));
            });
    functorpure.f=gauss;
    t1.checkpoint(true,"prep");

    auto gaussop=std::shared_ptr<SeparatedConvolution<double,LDIM>>(GaussOperatorPtr<LDIM>(world,gaussexponent));
    LRFunctorF12<double,NDIM> functorf12(gaussop,std::vector<Function<double,LDIM>>({}),{});
//    functorf12.f12.reset(GaussOperatorPtr<LDIM>(world,gaussexponent));

    auto builder= LowRankFunctionFactory<double,NDIM>(parameters).set_radius(8)
            .set_volume_element(0.1).set_rank_revealing_tol(1.e-10).set_orthomethod("canonical");

    auto lrfunction1=builder.project(functorf12);
    t1.checkpoint(true,"construction f12 functor");
    auto lrfunction2=builder.project(functorpure);
    t1.checkpoint(true,"construction full rank functor");
    lrfunction1.optimize(functorf12);
    t1.checkpoint(true,"optimization f12 functor");
    lrfunction2.optimize(functorpure);
    t1.checkpoint(true,"optimization full rank functor");


    try {
        double error1=lrfunction1.l2error(functorf12);
        t1.checkpoint(error1,tol,"f12 functor, f12 l2 error: ");
    } catch (...) {
        print("l2 error negative 1");
    }
    try {
        double error2=lrfunction2.l2error(functorpure);
        t1.checkpoint(error2,tol,"full rank functor, full rank l2 error");
    } catch (...) {
        print("l2 error negative 2");
    }
    try {
        double error3=lrfunction2.l2error(functorf12);
        t1.checkpoint(error3,tol,"full rank functor, f12 l2 error:");
    } catch (...) {
        print("l2 error negative 3");
    }
    try {
        double error4=lrfunction1.l2error(functorpure);
        t1.checkpoint(error4,tol,"f12 functor, full rank l2 error: ");
    } catch (...) {
        print("l2 error negative 4");
    }
//    print("errors",error2,error4);

    print_header2("leaving test_full_rank_functor");
    return t1.end();
}


template<std::size_t LDIM>
int test_arithmetic(World& world, LowRankFunctionParameters parameters) {
    constexpr std::size_t NDIM = 2 * LDIM;
    test_output t1("LowRankFunction::arithmetic in dimension " + std::to_string(NDIM));
    t1.set_cout_to_terminal();
    double thresh=FunctionDefaults<LDIM>::get_thresh();
    double thresh_ndim=FunctionDefaults<LDIM>::get_thresh();
    print("thresh ldim/ndim",thresh,thresh_ndim);
    Function<double,LDIM> phi=FunctionFactory<double,LDIM>(world)
            .functor([](const Vector<double,LDIM>& r){return exp(-4.0*inner(r,r));});

    auto gauss1=std::shared_ptr<SeparatedConvolution<double,LDIM>>(GaussOperatorPtr<LDIM>(world,1.0));
    LRFunctorF12<double,NDIM> functor1(gauss1,{phi},{});
    auto gauss2=std::shared_ptr<SeparatedConvolution<double,LDIM>>(GaussOperatorPtr<LDIM>(world,2.0));
    LRFunctorF12<double,NDIM> functor2(gauss2,{phi},{});

    // auto p1=particle<LDIM>::particle1();
    // auto p2=particle<LDIM>::particle2();

    auto builder= LowRankFunctionFactory<double,NDIM>(parameters).set_radius(4)
            .set_volume_element(0.1).set_rank_revealing_tol(1.e-10).set_orthomethod("canonical");
    auto lrf1=builder.project(functor1);
    auto lrf2=builder.project(functor2);

    Vector<double,NDIM> r;
    r.fill(0.2);

    // addition/subtraction
    {
        auto l1=lrf1+lrf1;
        t1.checkpoint(fabs(l1(r)-2.0*lrf1(r))<thresh,"addition - value");
        t1.checkpoint(l1.rank()==lrf1.rank()*2,"addition - rank");
        t1.checkpoint(&l1.get_g().front()!=&lrf1.get_g().front(),"addition - deep copy");

        auto l2=l1-lrf1;
        t1.checkpoint(fabs(l2(r)-lrf1(r))<thresh,"subtraction - value");
        t1.checkpoint(l2.rank()==lrf1.rank()*3,"subtraction - rank");
        t1.checkpoint(&l2.get_g().front()!=&lrf1.get_g().front(),"subtraction - deep copy");

        l2+=lrf1;
        t1.checkpoint(fabs(l2(r)-2.0*lrf1(r))<thresh,"in-place-addition - value");
        t1.checkpoint(l2.rank()==lrf1.rank()*4,"in-place-addition - rank");

        l2-=lrf1;
        t1.checkpoint(fabs(l2(r)-lrf1(r))<thresh,"in-place-subtraction - value");
        t1.checkpoint(l2.rank()==lrf1.rank()*5,"in-place-subtraction - rank");
    }


    // norm
    {
        double n1=lrf1.norm2();
        double refn1=functor1.norm2();
        t1.checkpoint(fabs(n1-refn1)<thresh,"norm2 computation");
        double n2=lrf2.norm2();
        double refn2=functor2.norm2();
        t1.checkpoint(fabs(n2-refn2)<thresh,"norm2 computation");
    }


    // scalar multiplication
    {
        auto l1=2.0*lrf1;
        t1.checkpoint(fabs(l1(r)-2.0*lrf1(r))<thresh,"oop-place multiplication");
        t1.checkpoint(&l1.get_g().front()!=&lrf1.get_g().front(),"subtraction - deep copy");
        l1*=0.5;
        t1.checkpoint(fabs(l1(r)-lrf1(r))<thresh,"in-place multiplication");

    }

    return t1.end();
}

/// test partial inner products

/// with
/// f1(x,y) = exp(-a*x^2) * exp(-(x-y)^2)
/// f2(x,y) = exp(-a*x^2) * exp(-g (x-y)^2)
/// test the inner products
///   a) inner(f1,f2,i,j) for i,j=0,1
///   b) inner(f1,
template<std::size_t LDIM>
int test_inner(World& world, LowRankFunctionParameters parameters) {

    static_assert(LDIM==1 or LDIM==2);
    constexpr std::size_t NDIM = 2 * LDIM;
    test_output t1("LowRankFunction::test_inner in dimension " + std::to_string(NDIM));
//    t1.set_cout_to_terminal();
    double thresh=FunctionDefaults<LDIM>::get_thresh();
    double thresh_ndim=FunctionDefaults<LDIM>::get_thresh();
    print("thresh ldim/ndim",thresh,thresh_ndim);
    Function<double,LDIM> phi=FunctionFactory<double,LDIM>(world)
            .functor([](const Vector<double,LDIM>& r){return exp(-4.0*inner(r,r));});

    auto gauss1=std::shared_ptr<SeparatedConvolution<double,LDIM>>(GaussOperatorPtr<LDIM>(world,1.0));
    LRFunctorF12<double,NDIM> functor1(gauss1,{phi},{});
    auto gauss2=std::shared_ptr<SeparatedConvolution<double,LDIM>>(GaussOperatorPtr<LDIM>(world,2.0));
    LRFunctorF12<double,NDIM> functor2(gauss2,{phi},{});
//   functor2.a={phi};

    auto p1=particle<LDIM>::particle1();
    auto p2=particle<LDIM>::particle2();

    auto builder= LowRankFunctionFactory<double,NDIM>(parameters).set_radius(4)
            .set_volume_element(0.1).set_rank_revealing_tol(1.e-10).set_orthomethod("canonical");
    auto lrf1=builder.project(functor1);
    auto lrf2=builder.project(functor2);

    // reference numbers: (by mathematica)
    // f1(x,y) = exp(-a*x^2) * exp(-(x-y)^2)
    // f2(x,y) = exp(-a*x^2) * exp(-g (x-y)^2)
    // with a=4, g=2
    // int f1(x,y),f2(x,z) dx = inner(f1,f2,0,0) : norm^2 = Pi^2/(2 Sqrt[2] Sqrt[a gamma] Sqrt[1 + 2 a + gamma]) = 0.37197471167788324677
    // int f1(x,y),f2(z,x) dx = inner(f1,f2,0,1) : norm^2 = Pi^2/(2 Sqrt[a (1 + a + gamma) (a + 2 gamma)]) = 0.32972034117743393239
    // int f1(y,x),f2(x,z) dx = inner(f1,f2,1,0) : norm^2 = 0.26921553123369812300
    // int f1(y,x),f2(z,x) dx = inner(f1,f2,1,1) : norm^2 = 0.35613867236025352322

    // inner f(1,2) f(2,3)
    auto fullrank1=lrf1.reconstruct();
    auto fullrank2=lrf2.reconstruct();
    t1.checkpoint(true,"prep inner");
    {
        std::vector<double> reference={ 0.37197471167788324677, 0.32972034117743393239, 0.26921553123369812300, 0.35613867236025352322};
        int counter=0;
        for (auto p11 : {p1,p2}) {
            for (auto p22 : {p1,p2}) {
                double ref=reference[counter];
                if (LDIM==2) ref*=ref;

                // full/full
                auto lhs1=inner(fullrank1,fullrank2,p11.get_tuple(),p22.get_tuple());
                const double l1=lhs1.norm2();
                print("l1",l1,l1*l1,ref);
                t1.checkpoint(l1*l1,ref,thresh,"inner(full,full,"+p11.str()+","+p22.str()+")");
                double asymmetric_ref=inner(fullrank1,lhs1);
                double asymmetric1=inner(fullrank1,lhs1);
                t1.checkpoint(asymmetric1,asymmetric_ref,thresh,"asymmetric(full,full,"+p11.str()+","+p22.str()+")");

                // low rank/full
                auto lhs2=inner(lrf1,fullrank2,p11,p22);
                double l2=lhs2.norm2();
                t1.checkpoint(l2*l2,ref,thresh,"inner(lrf,full,"+p11.str()+","+p22.str()+")");
                double asymmetric2=inner(lrf1,lhs2);
                t1.checkpoint(asymmetric2,asymmetric_ref,thresh,"asymmetric(lrf,full,"+p11.str()+","+p22.str()+")");

                // full/low rank
                auto lhs3=inner(fullrank1,lrf2,p11,p22);
                double l3=lhs3.norm2();
                t1.checkpoint(l3*l3,ref,thresh,"inner(full,lrf,"+p11.str()+","+p22.str()+")");
                double asymmetric3=inner(lrf1,lhs3);
                t1.checkpoint(asymmetric3,asymmetric_ref,thresh,"asymmetric(full,lrf,"+p11.str()+","+p22.str()+")");


                // low rank/low rank
                auto lhs4=inner(lrf1,lrf2,p11,p22);
                double l4=lhs4.norm2();
                t1.checkpoint(l4*l4,ref,thresh,"inner(lrf,lrf,"+p11.str()+","+p22.str()+")");
                double asymmetric4=inner(lrf1,lhs4);
                t1.checkpoint(asymmetric4,asymmetric_ref,thresh,"asymmetric(lrf,lrf,"+p11.str()+","+p22.str()+")");

                print("result norm",p11,p22,"), reference ",l1*l1,l2*l2,l3*l3,l4*l4,ref);
                // l2 norm cannot distinguish between f(1,2) and f(2,1)
                print("result asym",p11,p22,")",asymmetric1,asymmetric2,asymmetric3,asymmetric4);

                counter++;
            }
        }
    }

    // inner f(1,2) g(2)
    // this is surprisingly inaccurate, but algorithm is correct, the error can be systematically decreased
    {
        thresh=FunctionDefaults<LDIM>::get_thresh()*50.0;
        // fresh start
        lrf1=builder.project(functor1);
        fullrank1=FunctionFactory<double,NDIM>(world).functor(functor1);

        std::vector<Function<double,LDIM>> arg(3);
        for (int i=0; i<3; ++i) arg[i]=FunctionFactory<double,LDIM>(world)
                    .functor([&i](const Vector<double,LDIM>& r)
                             {return exp(-(i+1)*r.normf());});

        std::vector<Function<double,LDIM>> lhs_full1, lhs_full2,lhs_func1,lhs_func2;
        for (auto& a : arg) {
            lhs_full1.push_back(inner(fullrank1,a,p1.get_tuple(),p1.get_tuple()));
            lhs_full2.push_back(inner(fullrank1,a,p2.get_tuple(),p1.get_tuple()));

            lhs_func1.push_back(inner(functor1,a,p1,p1));
            lhs_func2.push_back(inner(functor1,a,p2,p1));
        }
        auto lhs_lrf1=inner(lrf1,arg,p1,p1);
        auto lhs_lrf2=inner(lrf1,arg,p2,p1);

        double norm_func1=norm2(world,lhs_func1);
        double norm_func2=norm2(world,lhs_func2);
        double norm_full1=norm2(world,lhs_full1);
        double norm_full2=norm2(world,lhs_full2);
        double norm_lrf1=norm2(world,lhs_lrf1);
        double norm_lrf2=norm2(world,lhs_lrf2);
        print("norms 1",norm_func1,norm_full1,norm_lrf1);
        print("norms 2",norm_func2,norm_full2,norm_lrf2);

        double error1=norm2(world,lhs_full1-lhs_func1);
        double error2=norm2(world,lhs_full2-lhs_func2);
        double error3=norm2(world,lhs_lrf1 -lhs_func1);
        double error4=norm2(world,lhs_lrf2 -lhs_func2);

        print("error1/2",error1,error2,error3,error4);
//        t1.checkpoint(error1<thresh,"inner(full(1,2),g(1)");      // particularly inaccurate, but we're not testing this
//        t1.checkpoint(error2<thresh,"inner(full(1,2),g(2)");      // particularly inaccurate, but we're not testing this
        t1.checkpoint(error3<thresh,"inner(lrf(1,2),g(1)");
        t1.checkpoint(error4<thresh,"inner(lrf(1,2),g(2)");
    }
    return t1.end();
}

template<std::size_t LDIM>
int test_construction_optimization(World& world, LowRankFunctionParameters parameters) {
    constexpr std::size_t NDIM=2*LDIM;
    test_output t1("LowRankFunction::construction/optimization in dimension "+std::to_string(NDIM));
    t1.set_cout_to_terminal();
    OperatorInfo info(1.0,1.e-6,FunctionDefaults<LDIM>::get_thresh(),OT_SLATER);
    auto slater=std::shared_ptr<SeparatedConvolution<double,LDIM> >(new SeparatedConvolution<double,LDIM>(world,info));
    Function<double,LDIM> one=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r){return exp(-0.4*inner(r,r));});
    Function<double,LDIM> half=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r){return sqrt(0.5)*exp(-0.4*inner(r,r));});

    LRFunctorF12<double,NDIM> lrfunctor1(slater,one,one);
    LRFunctorF12<double,NDIM> lrfunctor2(slater,{half,half},{half,half});

    for (auto& lrfunctor : {lrfunctor1,lrfunctor2}) {
        LowRankFunctionFactory<double, NDIM> builder(parameters);
        auto lrf = builder.project(lrfunctor);
        t1.checkpoint(lrf.rank() > 0, "construction");

        // with Slater tol must be relaxed
        double tol = 2.e-2;

        double error = lrf.l2error(lrfunctor);
        double norm = lrf.norm2();
        print("lrf.norm", norm);
        print("l2 error project ", error);
        t1.checkpoint(error, tol, "l2 error in projection");

        auto lrf2(lrf);
        error = lrf2.l2error(lrfunctor);
        print("l2 error copy ctor  ", error);
        MADNESS_CHECK(lrf.rank() == lrf2.rank());
        MADNESS_CHECK(&(lrf.g[0]) != &(lrf2.g[0]));  // deep copy
        t1.checkpoint(error, tol, "l2 error in copy ctor");

        lrf.optimize(lrfunctor);
        error = lrf.l2error(lrfunctor);
        print("l2 error optimize", error);
        t1.checkpoint(error, tol, "l2 error in optimization");

        print_header2(lrf.orthomethod);
        lrf.reorthonormalize();
        error = lrf.l2error(lrfunctor);
        print("l2 error reorthonormalize", error);
        t1.checkpoint(error, tol, "l2 error in reorthonormalization");

        lrf+=lrf;
        lrf*=0.5;
        lrf.reorthonormalize();
        error = lrf.l2error(lrfunctor);
        print("l2 error reorthonormalize with lindep", error);
        t1.checkpoint(error, tol, "l2 error in reorthonormalization with lindep");

    }
    return t1.end();
}

template<std::size_t LDIM>
int test_molecular_grid(World& world, LowRankFunctionParameters parameters) {
    constexpr std::size_t NDIM=2*LDIM;
    test_output t1("LowRankFunction::molecular_grid in dimension "+std::to_string(NDIM));
    t1.set_cout_to_terminal();
    OperatorInfo info(1.0,1.e-6,FunctionDefaults<LDIM>::get_thresh(),OT_SLATER);
    auto slater=std::shared_ptr<SeparatedConvolution<double,LDIM> >(new SeparatedConvolution<double,LDIM>(world,info));

    /*
     * we test the accuracy of the matrix element <dens(1) | f12 | dens(2)>_12
     */

    // a molecule of hydrogen atoms in the xy plane
    std::vector<Vector<double,LDIM>> atomic_sites;
    atomic_sites.push_back(Vector<double,LDIM>( 0.0));
    atomic_sites.push_back(Vector<double,LDIM>( 2.0));
    atomic_sites.push_back(Vector<double,LDIM>(-4.0));
    if (LDIM>1) {
        atomic_sites.back()[0]=0.0;
        atomic_sites.back()[1]=-4.0;
    }

    parameters.print("grid parameters in test_molecular_grid");
    molecular_grid<LDIM> mgrid(atomic_sites,parameters);
    auto grid=mgrid.get_grid();
    mgrid.visualize("grid",grid);

    // something like a density
    auto dens=[&atomic_sites](const Vector<double,LDIM>& r) {
        double result=0.0;
        for (auto& c : atomic_sites) result+=exp(-0.4*inner(r-c,r-c));
        return result;
    };

    Function<double,LDIM> density=FunctionFactory<double,LDIM>(world).functor(dens);
    Function<double,LDIM> one=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r){return exp(-0.4*inner(r,r));});
    Function<double,LDIM> half=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r){return sqrt(0.5)*exp(-0.4*inner(r,r));});

    LRFunctorF12<double,NDIM> lrfunctor(slater,density,density);
    {

//        atomic_sites.erase(atomic_sites.begin()+1, atomic_sites.end());
        LowRankFunctionFactory<double, NDIM> builder(parameters, atomic_sites);
        auto lrf = builder.project(lrfunctor);
        t1.checkpoint(lrf.rank() > 0, "construction");

        // with Slater tol must be relaxed
        double tol = 1.e-2;

        double error = lrf.l2error(lrfunctor);
        double norm = lrf.norm2();
        print("lrf.norm", norm);
        print("l2 error project ", error);
        t1.checkpoint(error, tol, "l2 error in projection");

        auto lrf2(lrf);
        error = lrf2.l2error(lrfunctor);
        print("l2 error copy ctor  ", error);
        MADNESS_CHECK(lrf.rank() == lrf2.rank());
        MADNESS_CHECK(&(lrf.g[0]) != &(lrf2.g[0]));  // deep copy
        t1.checkpoint(error, tol, "l2 error in copy ctor");

        lrf.optimize(lrfunctor,parameters.optimize());
        error = lrf.l2error(lrfunctor);
        print("l2 error optimize", error);
        t1.checkpoint(error, tol, "l2 error in optimization");

        lrf.reorthonormalize();
        error = lrf.l2error(lrfunctor);
        print("l2 error reorthonormalize", error);
        t1.checkpoint(error, tol, "l2 error in reorthonormalization");

        lrf+=lrf;
        lrf*=0.5;
        lrf.reorthonormalize();
        error = lrf.l2error(lrfunctor);
        print("l2 error reorthonormalize with lindep", error);
        t1.checkpoint(error, tol, "l2 error in reorthonormalization with lindep");

    }
    return t1.end();
}



int main(int argc, char **argv) {

    madness::World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);
    bool long_test = parser.key_exists("long_test");
    int k = parser.key_exists("k") ? std::atoi(parser.value("k").c_str()) : 6;
    double thresh  = parser.key_exists("thresh") ? std::stod(parser.value("thresh")) : 3.e-5;
    FunctionDefaults<6>::set_tensor_type(TT_2D);


    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<6>::set_truncate_mode(1);

    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<2>::set_thresh(thresh);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<4>::set_thresh(thresh);
    FunctionDefaults<5>::set_thresh(thresh);
    FunctionDefaults<6>::set_thresh(thresh);
    FunctionDefaults<6>::set_thresh(1.e-3);

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
    parameters.set_derived_value("f12type",std::string("slaterf12"));
    parameters.read_and_set_derived_values(world,parser,"grid");
    parameters.set_user_defined_value("radius",2.5);
    parameters.set_user_defined_value("volume_element",1.e-1);
//    parameters.set_user_defined_value("tol",1.0e-10);
    parameters.print("grid");
    int isuccess=0;
#ifdef USE_GENTENSOR

    // parameters.set_user_defined_value("volume_element",3.e-1);
    // parameters.set_user_defined_value("gridtype",std::string("random"));
    // isuccess+=test_molecular_grid<3>(world,parameters);
    // parameters.set_user_defined_value("gridtype",std::string("dftgrid"));
    // isuccess+=test_molecular_grid<3>(world,parameters);
    // parameters.set_user_defined_value("gridtype",std::string("random"));

    try {

        // isuccess+=test_full_rank_functor<1>(world, parameters);
        // parameters.set_user_defined_value("orthomethod",std::string("canonical"));
        // isuccess+=test_construction_optimization<1>(world,parameters);
        parameters.set_user_defined_value("orthomethod",std::string("cholesky"));
        isuccess+=test_construction_optimization<1>(world,parameters);
        // isuccess+=test_arithmetic<1>(world,parameters);
        // isuccess+=test_inner<1>(world,parameters);

        if (0) {
        // if (long_test) {
            isuccess+=test_construction_optimization<2>(world,parameters);
            isuccess+=test_arithmetic<2>(world,parameters);
            isuccess+=test_inner<2>(world,parameters);
            isuccess+=test_molecular_grid<2>(world,parameters);
        }

//        parameters.set_user_defined_value("volume_element",1.e-1);
//        isuccess+=test_lowrank_function(world,parameters);
//        isuccess+=test_Kcommutator(world,parameters);
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