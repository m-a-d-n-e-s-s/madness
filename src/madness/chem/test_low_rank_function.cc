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
    print("testing the representation of the function:   phi(1) phi(2) f(1,2)");
    print("by computing the norm of the projection : r(1) =  phi(1) int one(1) phi(2) f(1,2) d2");
    madness::default_random_generator.setstate(int(cpu_time())%4149);
    madness::default_random_generator.setstate(int(cpu_time())%4149);
    std::string id=unique_fileid();

    constexpr std::size_t LDIM=2;
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
    j["optimize"]=parameters.optimize();
    std::ofstream of(jsonfilename,std::ios::out);
    of<<j;
    of.close();

    Vector<double,LDIM> offset(0.0);
    Function<double,LDIM> phi1=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r)
                                                                           { return 1.0;});
    Function<double,LDIM> phi2=FunctionFactory<double,LDIM>(world).functor([&offset](const Vector<double,LDIM>& r)
                                                                           { return exp(-1.0*(r-offset).normf());});
    if (transpose=="slater1") std::swap(phi1,phi2);
    {
        double n1 = phi1.norm2();
        double n2 = phi2.norm2();
        Vector<double,LDIM> point(1.0);
        bool first_one = (fabs(phi1(point) - 1.0) < 1.e-6);
        if (world.rank() == 0) {
            if (first_one) print("1(1) phi(2)");
            else print("phi(1) 1(2)");
            print("norms", n1, n2);
        }
    }

    Function<double,LDIM> one=FunctionFactory<double,LDIM>(world)
            .functor([](const Vector<double,LDIM>& r) { return 1.0;});


    std::shared_ptr<SeparatedConvolution<double,LDIM>> f12;
    if (f12type=="slaterf12") f12.reset(SlaterF12OperatorPtr_ND<LDIM>(world,parameters.gamma(),1.e-6,FunctionDefaults<LDIM>::get_thresh()));
    else if (f12type=="slater") f12.reset(SlaterOperatorPtr_ND<LDIM>(world,parameters.gamma(),1.e-6,FunctionDefaults<LDIM>::get_thresh()));
    else {
        MADNESS_EXCEPTION(std::string("unknown f12type"+f12type).c_str(),1);
    }

    // compute relative l2 error between reference and result
    auto compute_relative_error = [](const auto reference, const auto result) {
        auto diff=reference-result;
        double refnorm=reference.norm2();
        // double resultnorm=result.norm2();
        double error=diff.norm2();
        return error/refnorm;
    };

    auto output=[&parameters] (const double projection_error, const Tensor<long> projection_rank, const double projection_time,
                               const double optimized_error, const Tensor<long> optimized_rank, const double optimized_time) {
        print("error",parameters.radius(),parameters.gridtype(),parameters.volume_element(),
              parameters.tol(),
              projection_error,projection_rank,projection_time,
              optimized_error,optimized_rank,optimized_time);
    };

    // \phi(1) \bar \phi(1) = \int phi(1) \phi(2) f(1,2) d2
    auto reference = phi1* (*f12)(phi2);
    plot_plane<LDIM>(world,reference,"reference."+id,PlotParameters(world).set_plane({"x1","x2"}));
    double n2=reference.norm2();
    print("reference.norm2() = int f12 phi2 d2",n2);
    output(0.0,Tensor<long>(2),0.0,0.0,Tensor<long>(2),0.0);

    LRFunctorF12<double,2*LDIM> lrfunctor(f12,phi1,phi2);
    double cpu0=cpu_time();
    auto lrf=LowRankFunctionFactory<double,2*LDIM>(parameters).project(lrfunctor);
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
    // Function<double,LDIM> result=compute_result(lrf);
    Function<double,LDIM> result=inner(lrf,one,particle<LDIM>::particle2(),particle<LDIM>::particle1());
    double projection_error=compute_relative_error(reference,result);
    auto diff=reference-result;
//    plot_plane<3>(world,diff,"plot_diff_int_projection."+id,PlotParameters(world).set_plane({"x1","x2"}));
//    plot_plane<3>(world,result,"plot_lrf_int_projection."+id,PlotParameters(world).set_plane({"x1","x2"}));
    output(projection_error,lrf.rank(),cpu1-cpu0,0.0,Tensor<long>(2),0.0);
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
    // result=compute_result(lrf);
    result=inner(lrf,one,particle<LDIM>::particle2(),particle<LDIM>::particle1());
    diff=reference-result;
    plot_plane<LDIM>(world,diff,"plot_diff_int_optimization."+id,PlotParameters(world).set_plane({"x1","x2"}));
    plot_plane<LDIM>(world,result,"plot_lrf_int_optimization."+id,PlotParameters(world).set_plane({"x1","x2"}));
    double optimization_error=compute_relative_error(reference,result);
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

int test_stuff(World& world, LowRankFunctionParameters parameters) {
    test_output t1("test_stuff");
    t1.set_cout_to_terminal();
    print("this is a placeholder for testing stuff");

    return t1.end();
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
    f12ptr.reset(SlaterF12OperatorPtr(world,parameters.gamma(),1.e-6,FunctionDefaults<LDIM>::get_thresh()));
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
    j["optimize"]=parameters.optimize();
    j["reference"]=reference;

    auto json2file= [](const json& j, const std::string& jsonfilename) {
        std::ofstream of(jsonfilename, std::ios::out);
        of << j;
        of.close();
    };

    json2file(j,jsonfilename);
    timer t(world);
    auto compute_error = [&](const std::string& msg, const LowRankFunction<double,6>& lrf) {
        auto gk = mul(world, phi_k, g12(lrf.g * phi_k)); // function of 1
        auto hj = lrf.h * phi; // function of 2
        Tensor<double> j_hj = inner(world, phi, hj);
        Tensor<double> i_gk = inner(world, phi, gk);

        // Contract through metric for the general representation: i_gk^T * M * j_hj.
        double result_right = lrf.is_canonical()
                ? i_gk.trace(j_hj)
                : i_gk.trace(inner(lrf.metric, j_hj, 1, 0));
        print(msg,"norm ", result_right);
        print(msg,"error", result_right-reference);
        print(msg,"rank ", lrf.rank());
        j[msg]=result_right-reference;
        j[msg+"_rank"]=std::to_string(lrf.rank()(0l))+" "+std::to_string(lrf.rank()(1));
        j[msg+"_error_compute_time"]=t.tag(msg+"_error_compute_time");
        json2file(j,jsonfilename);
    };


    if (true) {
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
        compute_error("exchange error",fi_one);
        t1.checkpoint(l2error,"l2error");


        // compute the exchange term twice, first time with no optimizatio of the low-rank function, second time with
        // optimization. The optimization should reduce the error significantly.
        for (int i=0; i<1; ++i) {
            if (i==1) {
                fi_one.optimize(lrfunctor);
                l2error=fi_one.l2error(lrfunctor);
                print("left_optimize_l2error",l2error);
                j["left_optimize_time"]=t.tag("left_optimize_time");
                json2file(j,jsonfilename);
                compute_error("left_optimize",fi_one);
            }

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
            t2.tag("multiply 2 ");
            std::string msg=(i==0) ? "LRF exchange commutator no optimization" : "LRF exchange commutator optimization";
            t1.checkpoint(result4,reference,1.e-4,msg);
        }
    }


    // apply exchange operator in 6d
//    if (f12type=="slaterf12") {
    if (true) {
//        FunctionDefaults<3>::print();
//        FunctionDefaults<6>::print();
        real_function_6d phi0=CompositeFactory<double,6,3>(world).particle1(phi).particle2(phi);

        double thresh=FunctionDefaults<3>::get_thresh();
        double dcut=1.e-6;
        real_function_6d tmp=TwoElectronFactory(world).dcut(dcut).gamma(parameters.gamma()).f12().thresh(thresh);
        real_function_6d f12ij=CompositeFactory<double,6,3>(world).g12(tmp).particle1(copy(phi)).particle2(copy(phi));

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
        double error=fabs(result-reference);
        print("error 6D <ij | K1 f12 | ij>", error);
        t1.checkpoint(error,1.e-4,"6D-function: exchange commutator");

    }

    return t1.end();

}

template<std::size_t LDIM>
int test_construction(World& world, LowRankFunctionParameters parameters) {

    test_output t1("test_construction");
    t1.set_cout_to_terminal();
    constexpr int NDIM=2*LDIM;
    double gaussexponent=2.0;
    double gauss1=2.0;
    double gauss2=3.0;
    parameters.set_derived_value("lmax",3);
    parameters.set_derived_value("tol",1.e-4);
    parameters.set_derived_value("radius",3.0);

    // make a functor that is not separable, but has a simple low-rank representation
    std::vector<std::shared_ptr<LRFunctorBase<double,NDIM>>> functors;
    Function<double,NDIM> gauss=FunctionFactory<double,NDIM>(world)
            .functor([&gaussexponent,&gauss1,&gauss2](const Vector<double,NDIM>& r){
                Vector<double,LDIM> a,b;
                for (int i=0; i<LDIM; ++i) {
                    a[i]=r[i];
                    b[i]=r[i+LDIM];
                }
                return exp(-gaussexponent*inner(a-b,a-b))
                        *exp(-gauss1*inner(a,a)) * exp(-gauss2*inner(b,b));
            });
    functors.push_back(std::shared_ptr<LRFunctorBase<double,NDIM>>(new LRFunctorPure<double,NDIM>(gauss)));
    t1.checkpoint(true,"prep");

    // same functor, different implementation
    auto gaussop=std::shared_ptr<SeparatedConvolution<double,LDIM>>(GaussOperatorPtr<LDIM>(world,gaussexponent));
    Function<double,LDIM> phi1=FunctionFactory<double,LDIM>(world)
            .functor([&gauss1](const Vector<double,LDIM>& r) {return exp(-gauss1*inner(r,r));});
    Function<double,LDIM> phi2=FunctionFactory<double,LDIM>(world)
            .functor([&gauss2](const Vector<double,LDIM>& r) {return exp(-gauss2*inner(r,r));});
    functors.push_back(std::shared_ptr<LRFunctorBase<double,NDIM>>(
        new LRFunctorF12<double,NDIM>(gaussop,phi1,phi2)));

    plot_plane<NDIM,LRFunctorBase<double,NDIM>>(world,*functors[0],"pure",PlotParameters(world).set_plane({"x1","x2"}));
    plot_plane<NDIM,LRFunctorBase<double,NDIM>>(world,*functors[1],"f12",PlotParameters(world).set_plane({"x1","x2"}));


    for (auto canonicalize : {true,false}) {
        for (auto gridtype : {"random","harmonics","adaptive"}) {
            for (auto ortho : {"cholesky","canonical"}) {
                for (auto& functor : functors) {
                    parameters.set_derived_value("gridtype",std::string(gridtype));
                    parameters.set_derived_value("orthomethod",std::string(ortho));
                    parameters.set_derived_value("canonicalize",canonicalize);
                    parameters.print("grid");

                    print("working with functor,",functor->type());

                    auto builder= LowRankFunctionFactory<double,NDIM>(parameters);
                    auto lrfunction1=builder.project(*functor);

                    // check the accuracy of the lrf projection and the l2error. The lrf projection is notoriously
                    // inaccurate, cf R12 error convergence, so we use a very loose tolerance here.
                    double tol_lrf=2.e-2;
                    double tol_lrf_l2=5.e-2;

                    // turn tol into a string with scientic notation and 2 digits, for better readability of the output
                    char buf[20];
                    snprintf(buf,20,"%.2e",parameters.tol());
                    std::string stol=std::string(buf);

                    std::string description=std::string(parameters.gridtype())+
                                        ", canon="+std::to_string(canonicalize)+
                                        ", "+functor->type()+
                                        ", tol="+stol;
                    t1.checkpoint(true,"lrf projection");
                    { // check function evaluation
                        Vector<double,NDIM> a;
                        a.fill(0.25);
                        double val=lrfunction1(a);
                        double ref1=(*functor)(a);
                        // double ref2=functorpure(a);
                        print("lr function evaluation, val, reference",val,ref1);
                        t1.checkpoint(val,ref1,tol_lrf,"lrf eval, "+description);
                    }
                    { // test l2 error
                        double error1=lrfunction1.l2error(*functor);
                        print("error, tol", error1,tol_lrf);
                        t1.checkpoint(error1,tol_lrf_l2,"lrf l2 err, "+description);
                    }
                }
            }
        }
    }
    return t1.end();
}


template<std::size_t LDIM>
int test_arithmetic(World& world, LowRankFunctionParameters parameters) {
    constexpr std::size_t NDIM = 2 * LDIM;
    test_output t1("LowRankFunction::arithmetic in dimension " + std::to_string(NDIM));
    t1.set_cout_to_terminal();
    double thresh=FunctionDefaults<LDIM>::get_thresh()*10;
    double thresh_ndim=FunctionDefaults<LDIM>::get_thresh();
    print("thresh ldim/ndim",thresh,thresh_ndim);
    Function<double,LDIM> phi=FunctionFactory<double,LDIM>(world)
            .functor([](const Vector<double,LDIM>& r){return exp(-4.0*inner(r,r));});

    auto gauss1=std::shared_ptr<SeparatedConvolution<double,LDIM>>(GaussOperatorPtr<LDIM>(world,1.0));
    LRFunctorF12<double,NDIM> functor1(gauss1,{phi},{});
    auto gauss2=std::shared_ptr<SeparatedConvolution<double,LDIM>>(GaussOperatorPtr<LDIM>(world,2.0));
    LRFunctorF12<double,NDIM> functor2(gauss2,{phi},{});

    for (auto gridtype : {"random","harmonics"}) {
        for (auto canonicalize : {true,false}) {
            parameters.set_derived_value("gridtype",std::string(gridtype));
            parameters.set_derived_value("canonicalize",canonicalize);
            parameters.print("grid");

            auto builder= LowRankFunctionFactory<double,NDIM>(parameters).set_radius(4)
                    .set_volume_element(0.1).set_rank_revealing_tol(1.e-5).set_orthomethod("canonical");
            auto lrf1=builder.project(functor1);
            auto lrf2=builder.project(functor2);

            Vector<double,NDIM> r;
            r.fill(0.2);

            // addition/subtraction
            {
                auto l1=lrf1+lrf1;
                t1.checkpoint(fabs(l1(r)-2.0*lrf1(r)),thresh,"addition - value");
                t1.checkpoint((l1.rank()-lrf1.rank()*2).sumsq()==0,"addition - rank");
                t1.checkpoint(&l1.get_g().front()!=&lrf1.get_g().front(),"addition - deep copy");

                auto l2=l1-lrf1;
                t1.checkpoint(fabs(l2(r)-lrf1(r)),thresh,"subtraction - value");
                t1.checkpoint((l2.rank()-lrf1.rank()*3).sumsq()==0,"subtraction - rank");
                t1.checkpoint(&l2.get_g().front()!=&lrf1.get_g().front(),"subtraction - deep copy");

                l2+=lrf1;
                t1.checkpoint(fabs(l2(r)-2.0*lrf1(r)),thresh,"in-place-addition - value");
                t1.checkpoint((l2.rank()-lrf1.rank()*4).sumsq()==0,"in-place-addition - rank");

                l2-=lrf1;
                t1.checkpoint(fabs(l2(r)-lrf1(r))<thresh,"in-place-subtraction - value");
                t1.checkpoint((l2.rank()-lrf1.rank()*5).sumsq()==0,"in-place-subtraction - rank");
            }


            // norm
            {
                double n1=lrf1.norm2();
                double refn1=functor1.norm2();
                t1.checkpoint(fabs(n1-refn1),thresh,"norm2 computation");
                double n2=lrf2.norm2();
                double refn2=functor2.norm2();
                t1.checkpoint(fabs(n2-refn2),thresh,"norm2 computation");
                auto lrf3=copy(lrf2);
                lrf3.canonicalize();
                double n3=lrf3.norm2();
                t1.checkpoint(fabs(n3-refn2),thresh,"norm2 computation - canon");
            }


            // scalar multiplication
            {
                auto l1=2.0*lrf1;
                t1.checkpoint(fabs(l1(r)-2.0*lrf1(r)),thresh,"oop-place multiplication");
                t1.checkpoint(&l1.get_g().front()!=&lrf1.get_g().front(),"subtraction - deep copy");
                l1*=0.5;
                t1.checkpoint(fabs(l1(r)-lrf1(r)),thresh,"in-place multiplication");

            }
        }
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
//   functor2.a={phi};

    auto p1=particle<LDIM>::particle1();
    auto p2=particle<LDIM>::particle2();

    auto builder= LowRankFunctionFactory<double,NDIM>(parameters).set_radius(4)
            .set_volume_element(0.1).set_rank_revealing_tol(1.e-6).set_orthomethod("canonical");
    auto lrf1=builder.project(functor1);
    auto lrf2=builder.project(functor2);

    // reference numbers: (by mathematica)
    // f1(x,y) = exp(-a*x^2) * exp(-(x-y)^2)
    // f2(x,y) = exp(-a*x^2) * exp(-g (x-y)^2)
    // with a=4, g=2
    // int f1(x,y),f2(z,x) dx = inner(f1,f2,0,0) : norm^2 = Pi^2/(2 Sqrt[2] Sqrt[a gamma] Sqrt[1 + 2 a + gamma]) = 0.37197471167788324677
    // int f1(x,y),f2(z,x) dx = inner(f1,f2,0,1) : norm^2 = Pi^2/(2 Sqrt[a (1 + a + gamma) (a + 2 gamma)]) = 0.32972034117743393239
    // int f1(y,x),f2(z,x) dx = inner(f1,f2,1,0) : norm^2 = 0.26921553123369812300
    // int f1(y,x),f2(z,x) dx = inner(f1,f2,1,1) : norm^2 = 0.35613867236025352322

    // inner f(1,2) f(2,3)
    auto fullrank1=lrf1.reconstruct();
    auto fullrank2=lrf2.reconstruct();
    double fnorm1=fullrank1.norm2();
    double fnorm2=fullrank2.norm2();
    print("full rank norms",fnorm1,fnorm2);
    double lrfnorm1=lrf1.norm2();
    double lrfnorm2=lrf2.norm2();
    print("low rank norms",lrfnorm1,lrfnorm2);
    double in1=inner(lrf1,lrf1);
    double in2=inner(lrf2,lrf2);
    print("inner 1,2",sqrt(in1),sqrt(in2));

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
int test_remove_lindep(World& world, LowRankFunctionParameters parameters) {
    constexpr std::size_t NDIM=2*LDIM;
    test_output t1("LowRankFunction::remove_lindep in dimension "+std::to_string(NDIM));
    t1.set_cout_to_terminal();
    OperatorInfo info(1.0,1.e-6,FunctionDefaults<LDIM>::get_thresh(),OT_SLATER);
    auto slater=std::shared_ptr<SeparatedConvolution<double,LDIM> >(new SeparatedConvolution<double,LDIM>(world,info));
    Function<double,LDIM> one=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r){return exp(-0.4*inner(r,r));});
    Function<double,LDIM> half=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r){return sqrt(0.5)*exp(-0.4*inner(r,r));});

    LRFunctorF12<double,NDIM> lrfunctor1(slater,one,one);
    LRFunctorF12<double,NDIM> lrfunctor2(slater,{half,half},{half,half});

    parameters.set_derived_value("volume_element",3.e-2);
    parameters.set_derived_value("tol",1.e-6);
    parameters.set_derived_value("lmax",3);
    parameters.print("grid");

    // for (auto& lrfunctor : {lrfunctor1,lrfunctor2}) {
    for (auto& lrfunctor : {lrfunctor2}) {
        for (auto& canonicalize : {true,false}) {
            for (auto gridtype : {"random","harmonics"}) {
                parameters.set_derived_value("gridtype",std::string(gridtype));
                print_header2("in functor loop");
                parameters.set_derived_value("canonicalize",canonicalize);
                MADNESS_CHECK_THROW(parameters.canonicalize() == canonicalize,"incorrect setting of canonicalize"); // make sure it isn't overridden
                std::string description=std::string(parameters.gridtype())+", canon="+std::to_string(canonicalize);

                LowRankFunctionFactory<double, NDIM> builder(parameters);
                auto lrf = builder.project(lrfunctor);

                // with Slater tol must be relaxed
                double tol = 2.e-2;
                if (not canonicalize) tol=5.e-2;        // sad..

                double error = lrf.l2error(lrfunctor);
                t1.checkpoint(error, tol, "l2 error in projection "+description);

                auto lrf2(lrf);
                if (canonicalize) MADNESS_CHECK((lrf.rank()-lrf2.rank()).sumsq()==0);
                else {
                    MADNESS_CHECK(lrf.metric.dim(0) == lrf.metric.dim(0));
                    MADNESS_CHECK(lrf.metric.dim(1) == lrf.metric.dim(1));
                }
                MADNESS_CHECK(&(lrf.g[0]) != &(lrf2.g[0]));  // deep copy
                error = lrf2.l2error(lrfunctor);
                t1.checkpoint(error, tol, "l2 error in copy ctor "+description);

                lrf.remove_linear_dependencies();
                error = lrf.l2error(lrfunctor);
                t1.checkpoint(error, tol, "l2 error in remove_lindep "+description);

                lrf+=lrf;
                lrf*=0.5;
                lrf.remove_linear_dependencies();
                error = lrf.l2error(lrfunctor);
                t1.checkpoint(error, tol, "l2 error in remove_lindep with lindep "+description);
            }
        }

    }
    return t1.end();
}

template<std::size_t LDIM>
int test_molecular_grid(World& world, LowRankFunctionParameters parameters) {
    constexpr std::size_t NDIM=2*LDIM;
    test_output t1("LowRankFunction::molecular_grid in dimension "+std::to_string(NDIM));
    t1.set_cout_to_terminal();

    // prepare a small set of atomic centers
    std::vector<Vector<double,LDIM>> atomic_sites;
    atomic_sites.push_back(Vector<double,LDIM>(0.0));
    atomic_sites.push_back(Vector<double,LDIM>(2.0));
    atomic_sites.push_back(Vector<double,LDIM>(-4.0));
    if (LDIM>1) {
        atomic_sites[0][0]=0.0;
        atomic_sites[0][1]=0.0;
        atomic_sites[1][0]=2.0;
        atomic_sites[1][1]=0.0;
        atomic_sites[2][0]=0.0;
        atomic_sites[2][1]=-4.0;
    }

    parameters.print("grid parameters in test_molecular_grid");

    // construct molecular grid
    molecular_grid<LDIM> mgrid(atomic_sites, parameters);
    auto grid = mgrid.get_grid();
    mgrid.visualize("grid", grid);

    // 1) construction: we must have produced some grid points
    t1.checkpoint(grid.size()>0, "construction: grid has points");

    // 2) randomness: if requested, two subsequent grids should not be identical
    if (parameters.gridtype() == "random") {
        auto grid2 = mgrid.get_grid();
        // check that at least one point differs by more than a tiny tolerance
        const double eps = 1e-12;
        bool different = false;
        size_t nmin = std::min(grid.size(), grid2.size());
        for (size_t i=0;i<nmin && !different;++i) {
            double d = (grid[i]-grid2[i]).normf();
            if (d > eps) different = true;
        }
        // if sizes differ that's also evidence of randomness
        if (grid.size() != grid2.size()) different = true;
        t1.checkpoint(different, "randomness: consecutive grids differ for gridtype=random");
    }

    // 3) locality: points should be local to their centers
    double radius = parameters.radius();
    std::vector<size_t> count_per_center(atomic_sites.size(), 0);

    for (const auto &p : grid) {
        // find nearest center
        double mind = std::numeric_limits<double>::infinity();
        size_t best = 0;
        for (size_t ic=0; ic<atomic_sites.size(); ++ic) {
            double d = (p - atomic_sites[ic]).normf();
            if (d < mind) { mind = d; best = ic; }
        }
        // count it
        if (mind <= 1.5*radius) ++count_per_center[best];
    }

    bool any_local = false;
    for (size_t ic=0; ic<count_per_center.size(); ++ic) {
        print("center", ic, "points within 1.5*radius", count_per_center[ic]);
        if (count_per_center[ic] > 0) any_local = true;
    }
    t1.checkpoint(any_local, "locality: at least one center has local points");

    // 4) density: compare number of points per center to expected value based on volume_element
    double expected_per_center = 1.0;
    if constexpr (LDIM==1) expected_per_center = 2.0 * radius / parameters.volume_element();
    if constexpr (LDIM==2) expected_per_center = constants::pi * radius * radius / parameters.volume_element();
    if constexpr (LDIM==3) expected_per_center = 4.0/3.0 * constants::pi * radius * radius * radius / parameters.volume_element();

    // be permissive: accept a broad range due to randomness and boundary effects
    double lower = std::max(1.0, 0.2 * expected_per_center);
    double upper = std::max(1.0, 5.0 * expected_per_center);
    bool density_ok = true;
    for (size_t ic=0; ic<count_per_center.size(); ++ic) {
        double c = double(count_per_center[ic]);
        if (c < lower || c > upper) {
            density_ok = false;
            print("density warning: center",ic,"count",c,"expected~",expected_per_center);
        }
    }
    t1.checkpoint(density_ok, "density: per-center counts within expected (permissive) range");

    return t1.end();
}

/// test norm2 with an asymmetric (non-symmetric) coupling matrix
///
/// Build a LowRankFunction in canonical form, compute its norm as reference,
/// then transform into a general representation with an asymmetric metric
/// and verify that norm2() still returns the correct value.
template<std::size_t LDIM>
int test_norm2_asymmetric_metric(World& world, LowRankFunctionParameters parameters) {
    constexpr std::size_t NDIM = 2 * LDIM;
    test_output t1("LowRankFunction::norm2 with asymmetric metric in dimension " + std::to_string(NDIM));
    t1.set_cout_to_terminal();
    double thresh = FunctionDefaults<LDIM>::get_thresh();

    // build a few basis functions for g and h
    const int nfunc = 3;
    std::vector<Function<double,LDIM>> gfuncs(nfunc), hfuncs(nfunc);
    for (int i = 0; i < nfunc; ++i) {
        double alpha_g = 1.0 + 0.5 * i;
        double alpha_h = 2.0 + 0.3 * i;
        gfuncs[i] = FunctionFactory<double,LDIM>(world)
                .functor([alpha_g](const Vector<double,LDIM>& r) { return exp(-alpha_g * inner(r,r)); });
        hfuncs[i] = FunctionFactory<double,LDIM>(world)
                .functor([alpha_h](const Vector<double,LDIM>& r) { return exp(-alpha_h * inner(r,r)); });
    }

    // 1) canonical LRF: f(1,2) = sum_i g_i(1) h_i(2)
    LowRankFunction<double,NDIM> lrf_canon(gfuncs, hfuncs, 1.e-8, "canonical");
    double norm_canon = lrf_canon.norm2();
    double norm_reconstruct = lrf_canon.reconstruct().norm2();
    print("canonical norm2, reconstruct norm2", norm_canon, norm_reconstruct);
    t1.checkpoint(fabs(norm_canon - norm_reconstruct) < 10.0 * thresh,
                  "canonical norm2 matches reconstruct");

    // 2) general LRF with symmetric metric: M = S (overlap of g)
    //    f(1,2) = sum_{ij} g_i(1) S^{-1}_{ij} S_{jk} h_k(2)  [just identity transform, same function]
    //    Instead, define: new_g = g, new_h = h, metric = I (trivially symmetric) -- boring.
    //    Better: transform g -> g' = g * A, metric = A^{-1}
    //    with A being an arbitrary invertible matrix (asymmetric).
    {
        // build an asymmetric, invertible matrix A
        Tensor<double> A(nfunc, nfunc);
        A(0,0) = 1.0; A(0,1) = 0.3; A(0,2) = 0.0;
        A(1,0) = 0.0; A(1,1) = 1.0; A(1,2) = 0.2;
        A(2,0) = 0.1; A(2,1) = 0.0; A(2,2) = 1.0;
        // A is not symmetric: A(0,1)=0.3 != A(1,0)=0.0

        // compute A^{-1} by solving A * Ainv = I
        Tensor<double> Ainv;
        {
            Tensor<double> I = Tensor<double>(nfunc, nfunc);
            for (int i = 0; i < nfunc; ++i) I(i,i) = 1.0;
            gesv(A, I, Ainv);
        }

        // transform: g' = g * A  (i.e. g'_j = sum_i g_i A_{ij})
        auto gprime = transform(world, gfuncs, A);

        // the function is the same: f = g' A^{-1} h = g * A * A^{-1} * h = g * h
        // so metric = A^{-1}, which is asymmetric
        LowRankFunction<double,NDIM> lrf_asym(gprime, hfuncs, 1.e-8, "canonical", Ainv);
        MADNESS_CHECK(!lrf_asym.is_canonical()); // metric is set

        double norm_asym = lrf_asym.norm2();
        double norm_asym_reconstruct = lrf_asym.reconstruct().norm2();
        print("asymmetric metric norm2, reconstruct norm2", norm_asym, norm_asym_reconstruct);
        print("asymmetric metric norm2 vs canonical norm2", norm_asym, norm_canon);

        // the function hasn't changed, so norm must match the canonical reference
        t1.checkpoint(fabs(norm_asym - norm_canon), 10.0 * thresh,
                      "asymmetric metric: norm2 matches canonical norm2");
        t1.checkpoint(fabs(norm_asym - norm_asym_reconstruct), 10.0 * thresh,
                      "asymmetric metric: norm2 matches reconstruct norm2");

        // 2b) function evaluation with asymmetric metric
        {
            Vector<double,NDIM> r;
            r.fill(0.2);
            double val_canon = lrf_canon(r);
            double val_asym = lrf_asym(r);
            print("eval canonical, asymmetric", val_canon, val_asym);
            t1.checkpoint(fabs(val_canon - val_asym), 10.0 * thresh,
                          "asymmetric metric: operator() matches canonical");
        }

        // 2c) inner product: inner(canonical, general_asymmetric) should equal inner(canonical, canonical)
        {
            double ip_cc = inner(lrf_canon, lrf_canon);
            double ip_ca = inner(lrf_canon, lrf_asym);
            double ip_ac = inner(lrf_asym, lrf_canon);
            double ip_aa = inner(lrf_asym, lrf_asym);
            print("inner cc, ca, ac, aa", ip_cc, ip_ca, ip_ac, ip_aa);
            t1.checkpoint(fabs(ip_cc - ip_ca), 10.0 * thresh,
                          "asymmetric metric: inner(canon,asym) matches inner(canon,canon)");
            t1.checkpoint(fabs(ip_cc - ip_ac), 10.0 * thresh,
                          "asymmetric metric: inner(asym,canon) matches inner(canon,canon)");
            t1.checkpoint(fabs(ip_cc - ip_aa), 10.0 * thresh,
                          "asymmetric metric: inner(asym,asym) matches inner(canon,canon)");
        }

        // 2d) canonicalize from asymmetric metric
        {
            auto lrf_recanonicalized = copy(lrf_asym);
            lrf_recanonicalized.canonicalize();
            t1.checkpoint(lrf_recanonicalized.is_canonical(),
                          "asymmetric metric: canonicalize sets is_canonical()");
            double norm_recanon = lrf_recanonicalized.norm2();
            t1.checkpoint(fabs(norm_recanon - norm_canon), 10.0 * thresh,
                          "asymmetric metric: norm2 after canonicalize matches");
            Vector<double,NDIM> r;
            r.fill(0.2);
            double val_canon = lrf_canon(r);
            double val_recanon = lrf_recanonicalized(r);
            t1.checkpoint(fabs(val_canon - val_recanon), 10.0 * thresh,
                          "asymmetric metric: eval after canonicalize matches");
        }
    }

    // 3) another asymmetric metric obtained via remove_linear_dependencies
    //    after lrf += lrf, the metric becomes block-diagonal, and remove_linear_dependencies
    //    can produce an asymmetric metric
    {
        auto lrf_dup = copy(lrf_canon);
        lrf_dup += lrf_canon;
        lrf_dup *= 0.5;
        // at this point lrf_dup has a block-diagonal metric (from +=) but represents the same function
        double norm_dup = lrf_dup.norm2();
        print("duplicated norm2", norm_dup);
        t1.checkpoint(fabs(norm_dup - norm_canon), 10.0 * thresh,
                      "duplicated: norm2 matches canonical");

        lrf_dup.remove_linear_dependencies();
        // after remove_linear_dependencies the metric is generally asymmetric
        double norm_lindep = lrf_dup.norm2();
        double norm_lindep_reconstruct = lrf_dup.reconstruct().norm2();
        print("after remove_lindep: norm2, reconstruct norm2, is_canonical",
              norm_lindep, norm_lindep_reconstruct, lrf_dup.is_canonical());
        t1.checkpoint(fabs(norm_lindep - norm_canon), 10.0 * thresh,
                      "remove_lindep: norm2 matches canonical");
        t1.checkpoint(fabs(norm_lindep - norm_lindep_reconstruct), 10.0 * thresh,
                      "remove_lindep: norm2 matches reconstruct");
    }

    return t1.end();
}

/// make an RI basis for an atom, namely Gaussian functions
template<std::size_t LDIM>
int make_ri_basis(World& world, LowRankFunctionParameters parameters) {

    constexpr std::size_t NDIM=2*LDIM;

    // set up possible loops
    auto v_rhs = std::vector<std::string>({"exponential"});
    auto v_tol = std::vector<double>({1.e-4,1.e-6,1.e-8});
    auto v_lmax=std::vector<int>({2});
    auto v_canonicalize=std::vector<bool>({true,false});


    for (auto canonicalize : v_canonicalize) {
        for (auto gridtype : v_rhs) {
            for (auto tol : v_tol) {
            // for (auto lmax : v_lmax) {
                 timer t_total(world);
                 print("canonicalize",canonicalize);
                 print("gridtype",gridtype);
                 parameters.set_derived_value("gridtype",std::string(gridtype));
                 parameters.set_derived_value("canonicalize",bool(canonicalize));
                 parameters.set_derived_value("tol",tol);
                 parameters.set_derived_value("radius",1.8);
                 // parameters.set_derived_value("lmax",lmax);
                 print("parameters.lmax",parameters.lmax());

                 parameters.print("in make_ri_basis");
                 auto builder= LowRankFunctionFactory<double,NDIM>(parameters);

                 std::shared_ptr<SeparatedConvolution<double,LDIM>> f12;
                 std::string f12type=parameters.f12type();
                 if (f12type=="slaterf12") {
                     f12.reset(SlaterF12OperatorPtr_ND<LDIM>(world,parameters.gamma(),1.e-6,FunctionDefaults<LDIM>::get_thresh()));
                 } else if (f12type=="slater") {
                     f12.reset(SlaterOperatorPtr_ND<LDIM>(world,parameters.gamma(),1.e-6,FunctionDefaults<LDIM>::get_thresh()));
                 } else {
                     MADNESS_EXCEPTION(std::string("unknown f12type"+f12type).c_str(),1);
                 }
                 // a trial function
                 Function<double,3> phi=FunctionFactory<double,3>(world)
                                 .functor([](const Vector<double,3>& r){return exp(-1.2*inner(r,r));});
                 LRFunctorF12<double,NDIM> functorf12(f12,std::vector<Function<double,LDIM>>({phi}),{phi});

                 auto lrfunction1=builder.project(functorf12);
                 double error=lrfunction1.l2error(functorf12);
                 print("l2error(lrfunction1)", error);


                 world.gop.fence();
                 timer t1(world);
                 auto trial=(*f12)(phi*phi);
                 double trialnorm=inner(phi,trial);
                 print("norm( f12(phi))",trialnorm);

                 auto p1=particle<LDIM>::particle1();
                 auto p2=particle<LDIM>::particle2();
                 auto result1=inner(lrfunction1,phi,p2,p1).trace();
                 t1.tag("inner(f12,phi,1,1)");
                 print("norm(result)",result1);
                 print("relative error",(trialnorm-result1)/trialnorm);

                 world.gop.fence();
                 auto zeta=parameters.tempered();
                 double elapsed = t_total.end("end loop");
                 printf("result: tol, vol, radius, canonicalize lmax, zeta, elapsed, l2error, rel. error: "
                        "%.1e, %.1e, %.1e, %d, %d %.1e, %.1e, %.1e, %e, %e, %e\n",
                        parameters.tol(),parameters.volume_element(),parameters.radius(), bool(canonicalize),
                        parameters.lmax(),zeta[0],zeta[1],zeta[2],elapsed, error,fabs(trialnorm-result1)/trialnorm);
            }
        }
    }
    return 0;
}

template<std::size_t LDIM>
int test_adaptive_grid_projection(World& world, LowRankFunctionParameters parameters) {
    constexpr std::size_t NDIM = 2 * LDIM;
    test_output t1("LowRankFunction::adaptive grid projection in dimension " + std::to_string(NDIM));
    t1.set_cout_to_terminal();

    parameters.set_derived_value("gridtype",std::string("adaptive"));
    parameters.set_derived_value("canonicalize",true);
    parameters.set_derived_value("radius",3.0);
    // parameters.set_derived_value("volume_element",0.08);
    parameters.set_derived_value("tol",1.e-4);
    // parameters.set_derived_value("adaptive_coarse_factor",8.0);
    // parameters.set_derived_value("adaptive_refine_radius",0.5);
    // parameters.set_derived_value("adaptive_significance_ratio",0.2);
    // parameters.set_derived_value("adaptive_max_centers",8);
    // parameters.set_derived_value("adaptive_min_centers",2);
    parameters.print("grid");

    // offset to avoid sampling the function at the origin where it is largest and easiest to approximate
    // fill offset with random numbers in the range [3,3]
    std::vector<Vector<double,LDIM>> offsets(4);
    offsets[0].fill(-3.0);
    offsets[1].fill(2.0);
    offsets[2].fill(6.0);
    offsets[4].fill(-16.0);


    print("set offsets to", offsets);


    auto op = std::shared_ptr<SeparatedConvolution<double,LDIM>>(GaussOperatorPtr<LDIM>(world,2.0));

    Function<double,LDIM> phi1 = FunctionFactory<double,LDIM>(world)
            .functor([&offsets](const Vector<double,LDIM>& r){ return exp(-1.0*inner(r-offsets[0],r-offsets[0])); });
    Function<double,LDIM> phi2 = FunctionFactory<double,LDIM>(world)
            .functor([&offsets](const Vector<double,LDIM>& r){ return 0.2*exp(-2.0*inner(r-offsets[1],r-offsets[1])); });
    Function<double,LDIM> phi3 = FunctionFactory<double,LDIM>(world)
            .functor([&offsets](const Vector<double,LDIM>& r){ return exp(-3.0*inner(r-offsets[2],r-offsets[2])); });
    Function<double,LDIM> one = FunctionFactory<double,LDIM>(world)
            .functor([&offsets](const Vector<double,LDIM>& r){ return 1.0; });
    LRFunctorF12<double,NDIM> functor(op,{phi1,phi2,phi3},{one,one,one});

    for (std::string gridtype : {"adaptive", "random","twostage"}) {
        parameters.set_derived_value("gridtype",std::string(gridtype));
        print_header2("testing gridtype="+parameters.gridtype());
        auto lrf = LowRankFunctionFactory<double,NDIM>(parameters, std::vector<Vector<double,LDIM>>(offsets)).project(functor);
        double error = lrf.l2error(functor);
        print("rank", lrf.rank(), "l2error", error);

        auto rank=lrf.rank();
        t1.checkpoint(rank(0l)>0 and rank(1l)>0, gridtype+" projection has non-empty rank");
        t1.checkpoint(error,1.e-2, gridtype+" projection yields bounded l2error");
    }
    return t1.end();
}

int main(int argc, char **argv) {

    madness::World& world = madness::initialize(argc, argv);
    startup(world, argc, argv);
    commandlineparser parser(argc, argv);
    // bool long_test = parser.key_exists("long_test");
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
    print("cell[1]",FunctionDefaults<1>::get_cell());

    FunctionDefaults<2>::set_tensor_type(TT_FULL);
    print("numerical parameters: k, eps(3D), eps(6D)", FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(),
          FunctionDefaults<6>::get_thresh());
    LowRankFunctionParameters parameters;
    parameters.set_derived_value("f12type",std::string("slater"));
    parameters.read_and_set_derived_values(world,parser,"grid");
    parameters.set_derived_value("radius",2.5);
    parameters.set_derived_value("volume_element",0.08);
    parameters.set_derived_value("tol",1.e-5);
    parameters.set_derived_value("tempered",std::vector<double>({1.e-3,1.e2,4.0}));
    parameters.print("grid");

    bool long_test=false;
    int isuccess=0;
    // isuccess+=test_Kcommutator(world,parameters);
    isuccess+=test_stuff(world,parameters);

    // parameters.set_user_defined_value("volume_element",3.e-1);
    parameters.set_derived_value("gridtype",std::string("random"));
//    isuccess+=test_molecular_grid<1>(world,parameters);
//    isuccess+=test_molecular_grid<2>(world,parameters);
//    isuccess+=test_molecular_grid<3>(world,parameters);

    try {

        // make_ri_basis<3>(world, parameters);
        isuccess+=test_construction<1>(world, parameters);
        isuccess+=test_adaptive_grid_projection<1>(world, parameters);
//        isuccess+=test_norm2_asymmetric_metric<1>(world, parameters);
//        isuccess+=test_remove_lindep<1>(world,parameters);
//        isuccess+=test_arithmetic<1>(world,parameters);
//        isuccess+=test_inner<1>(world,parameters);

        if (long_test) {
            isuccess+=test_remove_lindep<2>(world,parameters);
            isuccess+=test_arithmetic<2>(world,parameters);
            isuccess+=test_inner<2>(world,parameters);
            isuccess+=test_molecular_grid<2>(world,parameters);
        }

    } catch (std::exception& e) {
        madness::print("an error occured");
        madness::print(e.what());
        isuccess+=1;
    }
    finalize();

    return isuccess;
}

















