//
// Created by Florian Bischoff on 6/27/22.
//


#include<madness/mra/mra.h>
#include<madness/chem/ccpairfunction.h>
#include<madness/chem/correlationfactor.h>
#include<madness/chem/electronic_correlation_factor.h>
#include<madness/chem/CCStructures.h>
#include<madness/chem/CCPotentials.h>
#include<madness/chem/projector.h>
#include<madness/chem/lowrankfunction.h>

#include<madness/world/test_utilities.h>

using namespace madness;

struct XParameters : QCCalculationParametersBase {

    XParameters() : QCCalculationParametersBase() {

        // initialize with: key, value, comment (optional), allowed values (optional)
        initialize<double>("radius",5.0,"the radius");
        initialize<int>("rank",500,"the number of grid points in random grids");
        initialize<std::string>("gridtype","random","the grid type",{"random","cartesian"});
    }

    void read_and_set_derived_values(World& world, const commandlineparser& parser, std::string tag) {
        read_input_and_commandline_options(world,parser,tag);
    }

    double radius() const {return get<double>("radius");}
    int rank() const {return get<int>("rank");}
    std::string gridtype() const {return get<std::string>("gridtype");}
};


bool longtest=false;
struct data {
    real_function_3d f1,f2,f3,f4,f5;
    real_function_6d f12,f23;

    std::shared_ptr<CCConvolutionOperator> f12_op;
    data() {}

    data(World& world, CCParameters& parameters) {
        f12_op.reset(new CCConvolutionOperator(world,OT_F12,parameters));
        if (not f1.is_initialized()) initialize(world);
    }

    void initialize(World& world) {
        auto g1 = [](const coord_3d& r) { return exp(-1.0 * inner(r, r)); };
        auto g2 = [](const coord_3d& r) { return exp(-2.0 * inner(r, r)); };
        auto g3 = [](const coord_3d& r) { return exp(-3.0 * inner(r, r)); };
        auto g4 = [](const coord_3d& r) { return exp(-4.0 * inner(r, r)); };
        auto g5 = [](const coord_3d& r) { return exp(-5.0 * inner(r, r)); };
        f1=real_factory_3d(world).f(g1);
        f2=real_factory_3d(world).f(g2);
        f3=real_factory_3d(world).f(g3);
        f4=real_factory_3d(world).f(g4);
        f5=real_factory_3d(world).f(g5);

        auto g = [](const coord_6d& r) {
            double r1=r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
            double r2=r[3]*r[3] + r[4]*r[4] + r[5]*r[5];
            return exp(-1.0*r1 - 2.0*r2);
        };
        auto g23 = [](const coord_6d& r) {
            double r1=r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
            double r2=r[3]*r[3] + r[4]*r[4] + r[5]*r[5];
            return exp(-1.0*r1 - 2.0*r2) + exp(-2.0*r1 - 3.0*r2);
        };
        f12 = real_factory_6d(world).f(g);
        f23 = real_factory_6d(world).f(g23);

    }
    void clear() {
        f12_op.reset();
        f1.clear();
        f2.clear();
        f3.clear();
        f4.clear();
        f5.clear();
        f12.clear();
        f23.clear();
    }

    /// get some standard functions

    /// f1: exp(-1.0 r^2)
    /// f2: exp(-2.0 r^2)
    /// f3: exp(-3.0 r^2)
    /// f4: exp(-4.0 r^2)
    /// f5: exp(-5.0 r^2)
    /// f12: exp(-r_1^2 - 2 r_2^2)
    /// f23: exp(-r_1^2 - 2 r_2^2) + exp(-2 r_1^2 - 3 r_2^2)
    auto get_functions() const {
        return std::make_tuple(f1,f2,f3,f4,f5,f12);
    }

    /// get some standard ccpairfunctions

    /// p1: pure, corresponds to f12
    /// p2: dec, corresponds to f23
    /// p3: op_dec, corresponds to f23
    /// p4: pure, corresponds to f23
    /// p5: op_pure, corresponds to f23
    auto get_ccpairfunctions() const {
        CCPairFunction p1(f12);
        CCPairFunction p2({f1,f2},{f2,f3});
        CCPairFunction p3(f12_op,{f1,f2},{f2,f3});
        CCPairFunction p4(f23); // two-term, corresponds to p2
        CCPairFunction p5(f12_op,f23); // two-term, corresponds to p2
        return std::make_tuple(p1,p2,p3,p4,p5);
    }

};

data data1;

/// Computes the electrostatic potential due to a Gaussian charge distribution

/// stolen from testsuite.cc
class GaussianPotential : public FunctionFunctorInterface<double,3> {
public:
    typedef Vector<double,3> coordT;
    const coordT center;
    const double exponent;
    const double coefficient;

    GaussianPotential(const coordT& center, double expnt, double coefficient)
            : center(center)
            , exponent(sqrt(expnt))
            , coefficient(coefficient*pow(constants::pi/exponent,1.5)*pow(expnt,-0.75)) {}

    double operator()(const coordT& x) const {
        double sum = 00;
        for (int i=0; i<3; ++i) {
            double xx = center[i]-x[i];
            sum += xx*xx;
        };
        double r = sqrt(sum);
        if (r<1.e-4) {	// correct thru order r^3
            const double sqrtpi=sqrt(constants::pi);
            const double a=exponent;
            return coefficient*(2.0*a/sqrtpi - 2.0*a*a*a*r*r/(3.0*sqrtpi));
        } else {
            return coefficient*erf(exponent*r)/r;
        }
    }
};

int test_lowrank_function3(World& world, XParameters& parameters) {
    test_output t1("CCPairFunction::low rank function");
    t1.set_cout_to_terminal();
    madness::default_random_generator.setstate(int(cpu_time())%4149);

    constexpr std::size_t LDIM=3;
    constexpr std::size_t NDIM=2*LDIM;
    print("eps, k, NDIM",FunctionDefaults<NDIM>::get_thresh(),FunctionDefaults<NDIM>::get_k(),NDIM);

    Vector<double,LDIM> offset;
    offset.fill(0.0);
//    Function<double,LDIM> phi1=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r)
//            { return exp(-r.normf());});
    Function<double,LDIM> phi1=FunctionFactory<double,LDIM>(world).functor([](const Vector<double,LDIM>& r)
                                                                           { return 1.0;});
    Function<double,LDIM> phi2=FunctionFactory<double,LDIM>(world).functor([&offset](const Vector<double,LDIM>& r)
                                                                           { return exp(-1.0*(r-offset).normf());});

    std::shared_ptr<real_convolution_3d> f12(SlaterOperatorPtr(world,1.0,1.e-6,FunctionDefaults<LDIM>::get_thresh()));

    auto f = [](const coord_6d& r) {
        coord_3d r1={r[0],r[1],r[2]};
        coord_3d r2={r[3],r[4],r[5]};
        return exp(-(r1-r2).normf() -r2.normf());
    };

    LowRank<double,6> lrf(f12,copy(phi1),copy(phi2));
    lrf.project(parameters.rank(),parameters.radius(),parameters.gridtype());
//    lrf.optimize(1);
    print("lrf.rank()",lrf.rank());

    plot_plane<6>(world,lrf.lrfunctor,"plot_f12_r2");
    plot_plane<6>(world,lrf,"lrf_6d");

    // compare
    // \phi(1) \bar \phi(1) = \int phi(1) \phi(2) f(1,2) d2
    //       = \int \sum_r\phi(1) g_r(1) h_r(2) \phi(2) d2
    //       = \phi(1) \sum_r g_r(1) <\phi|h_r>
    auto reference = phi1* (*f12)(phi2);
    real_function_3d result=real_factory_3d(world);
    for (int r=0; r<lrf.rank(); ++r) result+=lrf.g[r]*lrf.h[r].trace();
    auto diff=reference-result;

    double refnorm=reference.norm2();
    double resultnorm=result.norm2();
    double error=diff.norm2();
    print("refnorm, resultnorm, abs. error, rel. error",refnorm, resultnorm, error, error/refnorm);

    plot<LDIM>({reference, result, diff}, "f_and_approx", std::vector<std::string>({"adsf", "asdf", "diff"}));



    return t1.end();

}

int test_lowrank_function(World& world) {
    test_output t1("CCPairFunction::low rank function");
    t1.set_cout_to_terminal();
    madness::default_random_generator.setstate(int(cpu_time())%4149);

    constexpr std::size_t LDIM=1;
    constexpr std::size_t NDIM=2*LDIM;
    print("eps, k, NDIM",FunctionDefaults<NDIM>::get_thresh(),FunctionDefaults<NDIM>::get_k(),NDIM);

    Function<double,NDIM> f12=FunctionFactory<double,NDIM>(world).functor([&LDIM](const Vector<double,NDIM>& r)
        {
            Vector<double,LDIM> r1,r2;
            for (int i=0; i<LDIM; ++i) {
                r1[i]=r[i];
                r2[i]=r[i+LDIM];
            }
            return exp(-(r1-r2).normf());//* exp(-0.2*inner(r1,r1));
        });
    Function<double,NDIM> phi0=FunctionFactory<double,NDIM>(world).functor([&LDIM](const Vector<double,NDIM>& r)
          {
              Vector<double,LDIM> r1,r2;
              for (int i=0; i<LDIM; ++i) {
                r1[i]=r[i];
                r2[i]=r[i+LDIM];
              }
              return exp(-(r1).normf());//* exp(-0.2*inner(r1,r1));
          });
    Function<double,NDIM> phi1=FunctionFactory<double,NDIM>(world).functor([&LDIM](const Vector<double,NDIM>& r)
          {
            Vector<double,LDIM> r1,r2;
            for (int i=0; i<LDIM; ++i) {
                r1[i]=r[i];
                r2[i]=r[i+LDIM];
            }
            return exp(-(r2).normf());//* exp(-0.2*inner(r1,r1));
        });
    Function<double,NDIM> one=FunctionFactory<double,NDIM>(world)
            .functor([&LDIM](const Vector<double,NDIM>& r){ return 1.0; });

    // f(1,2) = f12 * phi(1) * phi(2)
    Function<double,NDIM> f = phi0*phi1 *f12;
    print("lhs = phi0 * phi1; rhs=phi0*f12");
    auto lhs=phi0*phi1;
    auto rhs=phi0*f12;
    auto rhsnorm=rhs.norm2();


    phi0.reconstruct();
    f12.reconstruct();
    Function<double,NDIM> reference=inner(lhs,rhs,{0},{0});
    LowRank<double,NDIM> lrf(rhs);
    lrf.project(50,3.0,"random");
    lrf.optimize();

    double ferror=lrf.error();
    print("fnorm, error, rel. error",rhsnorm,ferror,ferror/rhsnorm);


    auto result=inner(lhs,lrf,{0},{0});
    double refnorm=reference.norm2();

    auto reconstruct=result.reconstruct();
    auto diff=reference-reconstruct;

    double resultnorm=reconstruct.norm2();
    double diffnorm=diff.norm2();
    print("refnorm, resultnorm ",refnorm, resultnorm);
    print("resultnorm, error, rel. error",resultnorm,diffnorm,diffnorm/resultnorm);
    print("final - one: k,thresh ",FunctionDefaults<2>::get_k(),FunctionDefaults<2>::get_thresh(),
          ferror/rhsnorm, diffnorm/resultnorm);
    print("\n");

    plot<NDIM>({reference, reconstruct, diff}, "f_and_approx", std::vector<std::string>({"adsf", "asdf", "diff"}));

//    return 0;
    if (0) {
        Function<double,NDIM> f;
        double fnorm = f.norm2();
        print("norm(2D-f)", fnorm);
        t1.checkpoint(true, "hi-dim projection");

        long rank = 50;
        double radius = 5.0;

        LowRank<double, NDIM> lr(copy(f));
        lr.project(rank, radius,"random");
        double err = lr.error();
        print("error in f_approx", err);
        lr.orthonormalize(true);
        t1.checkpoint(true, "start stuff");

        err = lr.error();
        print("error in f_approx", err);
        t1.checkpoint(true, "compute error");
        if (LDIM < 3) {
            auto fapprox = lr.reconstruct();
            plot<NDIM>({f, fapprox, f - fapprox}, "f_and_approx", std::vector<std::string>({"adsf", "asdf", "diff"}));
            t1.checkpoint(true, "plotting");
        }

        /*
         * optimize
         */

        lr.optimize(2);
        t1.checkpoint(true, "optimization");


        err = lr.error();
        print("error in f_approx_opt", err);

        if (LDIM < 3) {
            auto fapprox = lr.reconstruct();
            plot<NDIM>({f, fapprox, f - fapprox}, "f_and_approx_opt",
                       std::vector<std::string>({"adsf", "asdf", "diff"}));
            t1.checkpoint(true, "plotting");
        }


        /// int f(1,2) f12(2,3) d2
        f12.print_tree();
        auto reference=inner(f,f12,{1},{0});
        auto hbar=copy(world,lr.h);
        for (auto& hh : hbar) hh=inner(f12,hh,{0},{0});
        auto result=lr;
        result.h=hbar;
        auto reference_approx=result.reconstruct();
        result.lrfunctor.f=reference;
        double error=result.error();
        print("error ",error);
        plot<NDIM>({reference, reference_approx, reference-reference_approx}, "contraction", std::vector<std::string>({"adsf", "asdf", "diff"}));
    }

//    auto result=




    return t1.end();
}

int test_constructor(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {
    test_output t1("CCPairFunction::constructor");

    real_function_6d f=real_factory_6d(world);
    auto [f1,f2,f3,f4,f5,ff]=data1.get_functions();

    vector_real_function_3d  a= zero_functions<double,3>(world,3);
    vector_real_function_3d  b= zero_functions<double,3>(world,3);
    auto f12=CCConvolutionOperatorPtr(world, OT_F12, parameter);
    t1.checkpoint(true,"preparation");

    CCPairFunction p1;
    CCPairFunction p2(f);
    CCPairFunction p3({f1,f2},{f1,f3});
    CCPairFunction p4(f12,{f1,f2},{f2,f3});
    t1.checkpoint(true,"construction");

    {
        MADNESS_CHECK(p2.is_pure());
        MADNESS_CHECK(!p2.is_decomposed());
        MADNESS_CHECK(!p2.is_op_decomposed());
        auto p = p2.pure();
        auto ff = p2.pure().get_function();
        MADNESS_CHECK(!(ff.get_impl()==f.get_impl())); // deep copy of f
    }
    t1.checkpoint(true,"checks on pure");

    {
        MADNESS_CHECK(!p3.is_pure());
        MADNESS_CHECK(p3.is_decomposed());
        MADNESS_CHECK(!p3.is_op_decomposed());
        auto a1=p3.get_a();
        auto b1=p3.get_b();
    }
    t1.checkpoint(true,"checks on decomposed");

    {
        MADNESS_CHECK(!p4.is_pure());
        MADNESS_CHECK(p4.is_decomposed());
        MADNESS_CHECK(p4.is_op_decomposed());
        auto a1=p4.get_a();
        auto b1=p4.get_b();
        auto op=p4.get_operator();
    }
    t1.checkpoint(true,"checks on op_decomposed");

    {
        CCPairFunction tmp;
        tmp=p2;
        MADNESS_CHECK(tmp.is_pure());
        MADNESS_CHECK(tmp.component==p2.component);

        tmp=p3;
        MADNESS_CHECK(!tmp.is_pure());
        MADNESS_CHECK(tmp.component==p3.component);

    }
    t1.checkpoint(true,"checks on assignment");

    return t1.end();
}

int test_operator_apply(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                         const CCParameters& parameter) {
    test_output t1("CCPairFunction::test_operator_apply");
//    t1.set_cout_to_terminal();

    double exponent=1.0; // corresponds to the exponent of data::f1 and data::ff
//    double coefficient=pow(1.0/constants::pi*exponent,0.5*3);
    double coefficient=1.0;
    const coord_3d center={0.0,0.0,-0.0};

    auto Gaussian = [&center, &exponent, &coefficient](const coord_3d& r) {
        return coefficient * exp(-exponent*inner(r-center,r-center));
    };

    // this won't work, the simulation cell is only [-1,1]..
    // Normalized Gaussian exponent a produces potential erf(sqrt(a)*r)/r
//    GaussianPotential gpot(center, exponent, coefficient);
//    real_function_3d op_a=real_factory_3d(world).functor(gpot);

    real_function_3d a=real_factory_3d(world).functor(Gaussian);
    exponent=2.0;
    real_function_3d b=real_factory_3d(world).functor(Gaussian);

    auto [f1,f2,f3,f4,f5,ff]=data1.get_functions();
    CCPairFunction c1(a,b);
    CCPairFunction c2(f1,f2);
//    CCPairFunction ref(op_a,b);
    auto op= CoulombOperator(world,1.e-5,FunctionDefaults<3>::get_thresh());
    op.print_timings=false;
    op.particle()=1;

    auto op_c1=op(c1);
    auto op_c2=op(c2);
//    double a1=inner(ref,ref);
    double norm1=inner(op_c1,op_c1);
    double norm2=inner(op_c2,op_c2);
    print("norm1,norm2",norm1,norm2);
    bool good=fabs(norm1-norm2)<FunctionDefaults<3>::get_thresh();
    t1.checkpoint(good,"op(xx)");
    return t1.end();
}


int test_transformations(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                     const CCParameters& parameter) {
    test_output t1("CCPairFunction::test_transformations");

    real_function_6d f=real_factory_6d(world);
    auto [f1,f2,f3,f4,f5,ff]=data1.get_functions();
    auto f12=CCConvolutionOperatorPtr(world, OT_F12, parameter);
    auto g12=CCConvolutionOperatorPtr(world, OT_G12, parameter);


    CCPairFunction p1(ff);
    t1.checkpoint(p1.is_pure(),"is_pure");
    t1.checkpoint(p1.is_pure_no_op(),"is_pure_no_op");

    CCPairFunction p2(f12,ff);
    t1.checkpoint(p2.is_pure(),"is_pure");
    t1.checkpoint(p2.is_op_pure(),"is_op_pure");
    t1.checkpoint(p2.is_convertible_to_pure_no_op(),"is_convertible_to_pure_no_op");
    CCPairFunction p3=copy(p2);
    p3.convert_to_pure_no_op_inplace();
    t1.checkpoint(p2.is_op_pure(),"is_op_pure");
    t1.checkpoint(p3.is_pure_no_op(),"is_pure_no_op");

    CCPairFunction p4(g12,ff);
    t1.checkpoint(p4.is_pure(),"is_pure");
    t1.checkpoint(p4.is_op_pure(),"is_op_pure");
    t1.checkpoint(not p4.is_convertible_to_pure_no_op(),"not is_convertible_to_pure_no_op");

    CCPairFunction p5(f12,f1,f2);
    t1.checkpoint(not p5.is_pure(),"is_pure");
    t1.checkpoint(p5.is_op_decomposed(),"is_op_decomposed");
    t1.checkpoint(p5.is_convertible_to_pure_no_op(),"is_convertible_to_pure_no_op");
    CCPairFunction p6=copy(p5);
    p6.convert_to_pure_no_op_inplace();
    t1.checkpoint(p6.is_pure_no_op(),"is_pure_no_op");

    return t1.end();
}

int test_multiply_with_f12(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                  const CCParameters& parameters) {
    test_output t1("CCPairFunction::test_multiply_with_f12");

    // p1: pure, corresponds to f12
    // p2: dec, corresponds to f23
    // p3: op_dec, corresponds to f23
    // p4: pure, corresponds to f23
    // p5: op_pure, corresponds to f23
    auto [p1,p2,p3,p4,p5]=data1.get_ccpairfunctions();  // p2-p5 correspond to f23
    auto f12=data1.f12_op;

    double thresh=FunctionDefaults<3>::get_thresh();

    // decomposed
    CCPairFunction tmp1=f12*p2;         // should now be identical to p3
    CCPairFunction tmp2=p2*f12;         // should now be identical to p3
    double ref=inner(p2,p3);

    double r1=inner(p2,tmp1);
    bool good=(fabs(ref-r1)<thresh);
    t1.checkpoint(good,"f(1,2)*"+p2.name());

    double r2=inner(p2,tmp2);
    good=(fabs(ref-r2)<thresh);
    t1.checkpoint(good,p2.name() + "f(1,2)");

    // pure
    tmp1=f12*p4;         // should now be identical to p5
    tmp2=p4*f12;         // should now be identical to p5
    ref=inner(p2,p5);

    r1=inner(p2,tmp1);
    good=(fabs(ref-r1)<thresh);
    t1.checkpoint(good,"f(1,2)*"+p4.name());

    r2=inner(p2,tmp2);
    good=(fabs(ref-r2)<thresh);
    t1.checkpoint(good,p4.name() + "f(1,2)");
    return t1.end();
}

int test_multiply(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameters) {
    test_output t1("CCPairFunction::test_multiply");

    // consistency check, relies on CCPairFunction::inner to work correctly

    double thresh=FunctionDefaults<3>::get_thresh();
    auto [p1,p2,p3,p4,p5]=data1.get_ccpairfunctions();  // p2-p5 correspond to f23
    auto [f1,f2,f3,f4,f5,f]=data1.get_functions();

    // reference value is <bra | f(1/2) p>  = <f1 f | p>
    CCPairFunction bra(f1,f2);
    CCPairFunction bra1(f1*f2,f2);
    CCPairFunction bra2(f1,f2*f2);
    for (auto& p : {p2,p3,p4,p5}) {

        auto tmp1=multiply(p,f2,{0,1,2});
        double ovlp1=inner(bra,tmp1);
        double ref1=p.has_operator() ? inner(bra1,p3) : inner(bra1,p2);

        bool good=(fabs(ovlp1-ref1)<thresh);
        t1.checkpoint(good,"f(1)*"+p.name());

        auto tmp2=multiply(p,f2,{3,4,5});
        double ovlp2=inner(bra,tmp2);
        double ref2=p.has_operator() ? inner(bra2,p3) : inner(bra2,p2);

        good=(fabs(ovlp2-ref2)<thresh);
        t1.checkpoint(good,"f(2)*"+p.name());
    }

    return t1.end();
}
int test_inner(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                 const CCParameters& parameters) {
    test_output t1("CCPairFunction::test_inner");
    t1.set_cout_to_terminal();

    /// f1: exp(-1.0 r^2)
    /// f2: exp(-2.0 r^2)
    /// f3: exp(-3.0 r^2)
    /// f4: exp(-4.0 r^2)
    /// f5: exp(-5.0 r^2)
    /// f: exp(-r_1^2 - 2 r_2^2)
    /// f23: exp(-r_1^2 - 2 r_2^2) + exp(-2 r_1^2 - 3 r_2^2)
    auto [f1,f2,f3,f4,f5,f]=data1.get_functions();
    /// p1: pure, corresponds to f12
    /// p2: dec, corresponds to f23
    /// p3: op_dec, corresponds to f23
    /// p4: pure, corresponds to f23
    /// p5: op_pure, corresponds to f23
    auto [p1,p2,p3,p4,p5]=data1.get_ccpairfunctions();
    auto f12 = *(data1.f12_op);

    /// results
    auto a=std::vector<real_function_3d>({f1,f2});
    auto b=std::vector<real_function_3d>({f2,f3});
    std::vector<real_function_3d> a_ij_functions, b_ij_functions;
    for (int i=0; i<a.size(); ++i) {
        for (int j=0; j<a.size(); ++j) {
            a_ij_functions.push_back(a[i]*a[j]);
        }
    }
    for (int i=0; i<b.size(); ++i) {
        for (int j=0; j<b.size(); ++j) {
            b_ij_functions.push_back(b[i]*b[j]);
        }
    }

    auto aij=matrix_inner(world,a,a);
    auto bij=matrix_inner(world,b,b);

    // <a_ib_i | a_jb_j> = \sum_{ij} <a_i|a_j> <b_i|b_j>
    double ab_ab=aij.trace(bij);

    // <a_ib_i | f | a_jb_j> = \sum_{ij}  < <a_i|f|a_j>_1(2) | b_ib_j(2) >
    double ab_f_ab=dot(world,f12(a_ij_functions),b_ij_functions).trace();

    // <a_ib_i | f2 | a_jb_j> = \sum_{ij}  < <a_i|f^2|a_j>_2 | b_ib_j(2) >
    // f^2 = 1/(4y^2)(1 - 2*f(y) + f2(2y)) , f2(2y) =f2(y)^2
    // operator apply of SlaterF12Operator includes a factor of 1/(2 gamma) and the identity
    // operator apply of SlaterOperator has no further terms
    const double y = parameters.gamma();
    SeparatedConvolution<double, 3> fop= SlaterOperator(world, y, parameters.lo(), parameters.thresh_bsh_3D());
    SeparatedConvolution<double, 3> fsq = SlaterOperator(world, 2.0 * y, parameters.lo(), parameters.thresh_bsh_3D());

    const double prefactor = 1.0 / (4 * y * y);
    const double ab_f2_ab = prefactor*( ab_ab
                            -2.0*dot(world,apply(world,fop,a_ij_functions),b_ij_functions).trace()
                            +dot(world,apply(world,fsq,a_ij_functions),b_ij_functions).trace() );


    for (auto& ket : {p2, p3, p4, p5}) {
        for (auto& bra : {p2, p3, p4, p5}) {
            double ref=0.0;
            if (bra.has_operator() and ket.has_operator()) ref=ab_f2_ab;
            if (bra.has_operator() and (not ket.has_operator())) ref=ab_f_ab;
            if ((not bra.has_operator()) and ket.has_operator()) ref=ab_f_ab;
            if ((not bra.has_operator()) and (not ket.has_operator())) ref=ab_ab;
            double result=inner(bra,ket);

            print(bra.name(true)+ket.name(),"ref, result, diff", ref, result, ref-result);
            double thresh=FunctionDefaults<3>::get_thresh();
            bool good=(fabs(result-ref)<thresh);
            t1.checkpoint(good,bra.name(true)+ket.name());
        }
    }
    return t1.end();
}


int test_partial_inner_6d(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                       const CCParameters& parameter) {

    CCTimer timer(world, "testing");
    test_output t1("CCPairFunction::test_partial_inner_6d");
    t1.set_cout_to_terminal();
    auto [f1,f2,f3,f4,f5,f] = data1.get_functions();

    std::vector<real_function_3d> a = {f1, f2};
    std::vector<real_function_3d> b = {f2, f3};

    auto f12=CCConvolutionOperatorPtr(world, OT_F12, parameter);

    auto [p1,p2,p3,p4,nil]=data1.get_ccpairfunctions();
    CCPairFunction p11({f1},{f1});
    CCPairFunction p12({f1},{f2});

    double g11=inner(f1,f1);
    double g22=inner(f2,f2);
    double g12=inner(f1,f2);
    double g13=inner(f1,f3);
    double g23=inner(f2,f3);
    double g33=inner(f3,f3);
    real_function_3d gf11=(*f12)(f1*f1);
    real_function_3d gf12=(*f12)(f1*f2);
    real_function_3d gf22=(*f12)(f2*f2);
    real_function_3d gf13=(*f12)(f1*f3);
    real_function_3d gf23=(*f12)(f2*f3);

    t1.checkpoint(true,"prep",timer.reset());

    // p1 = p12 = e(-r1 -2r2)
    // p2 = e(-r1) * e(-2r2) + e(-2r1) * e(-3e2)  separated
    // p4 = e(-r1) * e(-2r2) + e(-2r1) * e(-3e2)  6d
    for (auto test_p1 : {p2,p4}) {
        for (auto test_p2 : {p2,p4}) {
            CCPairFunction r1=inner(test_p1,test_p2,{0,1,2},{0,1,2});
            CCPairFunction r2=inner(test_p1,test_p2,{0,1,2},{3,4,5});
            CCPairFunction r3=inner(test_p1,test_p2,{3,4,5},{0,1,2});
            CCPairFunction r4=inner(test_p1,test_p2,{3,4,5},{3,4,5});

            double n1=inner(r1,p11);
            double n2=inner(r2,p11);
            double n3=inner(r3,p11);
            double n4=inner(r4,p11);

            double ref_n1=g11*g12*g12 + g12*g12*g13 + g12*g13*g12 + g22*g13*g13;
            double ref_n2=g12*g12*g11 + g13*g12*g12 + g22*g13*g11 + g23*g13*g12;
            double ref_n3=ref_n2;
            double ref_n4=g22*g11*g11 + g23*g11*g12 + g23*g12*g11 + g33*g12*g12;

            bool good=fabs(n1-ref_n1)<FunctionDefaults<3>::get_thresh();
            t1.checkpoint(good,test_p1.name(true)+test_p2.name()+" -- 1");
            good=fabs(n2-ref_n2)<FunctionDefaults<3>::get_thresh();
            t1.checkpoint(good,test_p1.name(true)+test_p2.name()+" -- 2");
            good=fabs(n3-ref_n3)<FunctionDefaults<3>::get_thresh();
            t1.checkpoint(good,test_p1.name(true)+test_p2.name()+" -- 3");
            good=fabs(n4-ref_n4)<FunctionDefaults<3>::get_thresh();
            t1.checkpoint(good,test_p1.name(true)+test_p2.name()+" -- 4",timer.reset());

        }
    }

    // test < sth | f(1,2) a(1)b(2) >
    // CCPairFunction p2({f1,f2},{f2,f3});
    // CCPairFunction p4(f23); // two-term, corresponds to p2
    CCPairFunction p5(f12,std::vector<real_function_3d>({f1}),std::vector<real_function_3d>({f2}));
    for (auto& test_p1 : {p2, p4}) {
        CCPairFunction r1=inner(test_p1,p5,{0,1,2},{0,1,2});
        CCPairFunction r2=inner(test_p1,p5,{0,1,2},{3,4,5});
        CCPairFunction r3=inner(test_p1,p5,{3,4,5},{0,1,2});
        CCPairFunction r4=inner(test_p1,p5,{3,4,5},{3,4,5});

        double ref_n1=inner(gf11,f2*f1) * g12 + inner(gf12,f1*f2) * g13;
        double ref_n2=inner(gf12,f1*f1) * g12 + inner(gf22,f1*f1) * g13;
        double ref_n3=inner(gf12,f2*f1) * g11 + inner(gf13,f1*f2) * g12;
        double ref_n4=inner(gf22,f1*f1) * g11 + inner(gf23,f1*f1) * g12;

        double n1=inner(r1,p11);
        double n2=inner(r2,p11);
        double n3=inner(r3,p11);
        double n4=inner(r4,p11);
//        print("n1, ref_n1",n1,ref_n1, n1-ref_n1);
//        print("n2, ref_n2",n2,ref_n2, n2-ref_n2);
//        print("n3, ref_n3",n3,ref_n3, n3-ref_n3);
//        print("n4, ref_n4",n4,ref_n4, n4-ref_n4);

        bool good=fabs(n1-ref_n1)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,test_p1.name(true)+p5.name()+" -- 1");
        good=fabs(n2-ref_n2)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,test_p1.name(true)+p5.name()+" -- 2");
        good=fabs(n3-ref_n3)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,test_p1.name(true)+p5.name()+" -- 3");
        good=fabs(n4-ref_n4)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,test_p1.name(true)+p5.name()+" -- 4",timer.reset());

    }

    // test < a(1) b(2) f(1,2) | f(1,3) c(1)d(3) >
    // CCPairFunction p3(f12_op.get(),{f1,f2},{f2,f3});
    // CCPairFunction p5(&f12,{f1},{f2});
    if (longtest) {
        CCPairFunction r1=inner(p3,p5,{0,1,2},{0,1,2});
        print("time after r1 ", timer.reset());
        double n1=inner(r1,p11);
        p3.convert_to_pure_no_op_inplace();
        p5.convert_to_pure_no_op_inplace();
        print("n1",n1);
        print("time after r2a ", timer.reset());
        CCPairFunction r1a=inner(p3,p5,{0,1,2},{0,1,2});
        double n1a=inner(r1a,p11);
        print("n1a",n1a);
        print("diff",n1-n1a);
        CCPairFunction r2=inner(p3,p5,{0,1,2},{3,4,5});
        print("time after r2 ", timer.reset());
        CCPairFunction r3=inner(p3,p5,{3,4,5},{0,1,2});
        print("time after r3 ", timer.reset());
        CCPairFunction r4=inner(p3,p5,{3,4,5},{3,4,5});
        print("time after r4 ", timer.reset());

    }


    return  (t1.get_final_success()) ? 0 : 1;
}
int test_partial_inner_3d(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {

    CCTimer timer(world, "testing");
    test_output t1("CCPairFunction::test_partial_inner");
    auto [f1,f2,f3,f4,f5,f] = data1.get_functions();

    std::vector<real_function_3d> a = {f1, f2};
    std::vector<real_function_3d> b = {f2, f3};

    auto f12=CCConvolutionOperatorPtr(world, OT_F12, parameter);

    CCPairFunction p1(f);   // e(-r1 - 2r2)
    CCPairFunction p2(a,b);
    CCPairFunction p3(f12,a,b);
    CCPairFunction p11({f1},{f1});
    CCPairFunction p12({f1},{f2});

    double g11=inner(f1,f1);
    double g22=inner(f2,f2);
    double g12=inner(f1,f2);
    double g13=inner(f1,f3);
    double g23=inner(f2,f3);
    real_function_3d gf11=(*f12)(f1*f1);
    real_function_3d gf12=(*f12)(f1*f2);
    real_function_3d gf13=(*f12)(f1*f3);
    print("g11, g22",g11,g22);

    print("time in preparation",timer.reset());
    t1.checkpoint(true,"prep");

    // test pure/3d
    {
        real_function_3d r=inner(p1,f1,{0,1,2},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=g11 * g22;
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,"pure -- 1");
    }
    // test pure/3d
    {
        real_function_3d r=inner(p1,f1,{3,4,5},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=g12 * g12;
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,"pure -- 2");
    }
    // test decomposed
    {
        real_function_3d r=inner(p2,f1,{0,1,2},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=g11 * g22 + g12 * g23;
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,"decomposed -- 1");
    }
    // test decomposed
    {
        real_function_3d r=inner(p2,f1,{3,4,5},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=g12 * g12 + g13 * g22;
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,"decomposed -- 2");
    }
    // test op_decomposed
    {
        // < f1 f2 | f | f1 f2 > + < f1 f2 | f | f2 f3>
        real_function_3d r=inner(p3,f1,{0,1,2},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=inner(gf11*f2,f2) + inner(gf12*f2,f3);
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,"op_decomposed -- 1");
    }
    // test op_decomposed
    {
        // < f1 f2 | f | f2 f1 > + < f1 f2 | f | f3 f2>
        real_function_3d r=inner(p3,f1,{3,4,5},{0,1,2});
        double norm=inner(r,f2);
        double ref_norm=inner(gf12*f1,f2) + inner(gf13*f2,f2);
        print("norm, ref_norm",norm,ref_norm);
        bool good=fabs(norm-ref_norm)<FunctionDefaults<3>::get_thresh();
        t1.checkpoint(good,"op_decomposed -- 2");
    }


    return  (t1.get_final_success()) ? 0 : 1;
}

int test_apply(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {
    test_output t1("CCPairFunction::test_apply");

    return  (t1.get_final_success()) ? 0 : 1;
}

int test_scalar_multiplication(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
               const CCParameters& parameter) {
    CCTimer timer(world, "testing");
    test_output t1("CCPairFunction::test_scalar_multiplication");

    auto [f1,f2,f3,f4,f5,f] = data1.get_functions();

    std::vector<real_function_3d> a = {f1, f2};
    std::vector<real_function_3d> b = {f3, f1};

    print("time in preparation",timer.reset());
    t1.checkpoint(true,"prep");

    CCPairFunction p(f);
    CCPairFunction p1(a,b);
    double norm1=inner(p,p1);
    double pnorm=inner(p,p);
    double p1norm=inner(p1,p1);
    print("pnorm,p1norm",pnorm,p1norm);
    double anorm=norm2(world,p1.get_a());
    double bnorm=norm2(world,p1.get_b());
    print("anorm,bornm",anorm,bnorm);

    p*=2.0;
    p1*=2.0;
    anorm=norm2(world,p1.get_a());
    bnorm=norm2(world,p1.get_b());
    print("anorm,bornm",anorm,bnorm);
    pnorm=inner(p,p);
    p1norm=inner(p1,p1);
    print("pnorm,p1norm",pnorm,p1norm);
    double norm2=inner(p,p1);
    print("norm1,norm2",norm1,norm2);

    print("time in scaling and inner",timer.reset());
    bool bsuccess=fabs(4.0*norm1-norm2)<FunctionDefaults<3>::get_thresh();
    t1.checkpoint(bsuccess,"scaling");

    t1.end();
    return  (t1.get_final_success()) ? 0 : 1;
}

int test_swap_particles(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                        const CCParameters& parameters) {
    test_output t1("CCPairFunction::swap_particles");
    CCTimer timer(world, "testing swap_particles");

    // prepare
    auto one = [](const coord_3d& r) { return 1.0; };
    real_function_3d R2 = real_factory_3d(world).f(one);

    auto [f1,f2,f3,f4,f5,f] = data1.get_functions();
    std::vector<real_function_3d> a = {f1, f2};
    std::vector<real_function_3d> b = {f3, f1};

    // test decomposed
    {
        CCPairFunction p1(a, b);
        CCPairFunction p2(b, a);

        double norm1 = inner(p1, p2.swap_particles(), R2);
        double norm1a = inner(p1, p1, R2);
        // <p1 | p2> = \sum_ij <a_i b_i | a_j b_j> = \sum_ij <a_i|a_j> <b_i|b_j>
        double norm2 = matrix_inner(world, a, a).emul(matrix_inner(world, b, b)).sum();
        print("norm1 ", norm1);
        print("norm1a", norm1a);
        print("norm2 ", norm2);
        t1.checkpoint(std::abs(norm1 - norm2) < FunctionDefaults<3>::get_thresh(), "swap_particles a,b");
    }

    // test pure
    {
        CCPairFunction p(f);
        CCPairFunction p_swapped=p.swap_particles();
        double pnorm=p.get_function().norm2();
        double psnorm=p_swapped.get_function().norm2();
        print("p/s norm",pnorm,psnorm);

        CCPairFunction p1({f1}, {f2});
        CCPairFunction p2({f2}, {f1});
        double ref1=inner(f1,f1)*inner(f2,f2);
        double ref2=inner(f1,f2)*inner(f2,f1);
        print("ref1/2",ref1,ref2);
        print("pref1/2",inner(p1,p1),inner(p1,p2));

        double norm12_12=inner(p,p1);
        double norm12_21=inner(p,p1.swap_particles());
        double norm12_12_again=inner(p,p1);
        double norm21_12=inner(p_swapped,p1);
        double norm21_21=inner(p_swapped,p1.swap_particles());

        print("norms in exp(-12 - 12):",norm12_12);
        print("norms in exp(-12 - 21):",norm12_21);
        print("norms in exp(-12 - 12) again:",norm12_12_again);
        print("norms in exp(-21 - 12):",norm21_12);
        print("norms in exp(-21 - 21):",norm21_21);

        double ref_12_12=inner(p1,p1);
        double ref_12_21=inner(p1,p2);

        print("difference norms in exp(-12 - 12):",norm12_12,ref_12_12);
        print("difference norms in exp(-12 - 21):",norm12_21,ref_12_21);
        print("difference norms in exp(-21 - 12):",norm21_12,ref_12_21);
        print("difference norms in exp(-21 - 21):",norm21_21,ref_12_12);

        double total_error= fabs(norm12_12-ref_12_12)+ fabs(norm12_21-ref_12_21)
                + fabs(norm21_12-ref_12_21)+ fabs(norm21_21-ref_12_12);


        t1.checkpoint(total_error < FunctionDefaults<3>::get_thresh(), "swap_particles u");
    };
    t1.end();

    return  (t1.get_final_success()) ? 0 : 1;
}

int test_projector(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                           const CCParameters& parameter) {
    test_output t1("CCPairFunction::test_projector");
    CCTimer timer(world, "testing");

    t1.set_cout_to_logger();
    auto [f1,f2,f3,f4,f5,f] = data1.get_functions();
    double nf1=f1.norm2();
    double nf2=f2.norm2();
    f1.scale(1.0/nf1);
    f2.scale(1.0/nf1);
    std::vector<real_function_3d> a = {f1+f2, f2+f3, f3+f4};
    std::vector<real_function_3d> b = {f3, f1, f2};
    std::vector<real_function_3d> o = orthonormalize_cd<double,3>({f1-f3, f5});  // projects on an orthonormal basis
    a = orthonormalize_cd<double,3>(a);  // projects on an orthonormal basis
    auto f12=CCConvolutionOperatorPtr(world, OT_F12, parameter);

    {
        double n1=inner(f1,o[0]);
        double n2=inner(f2,o[0]);
        auto ovlp=matrix_inner(world,o,o);
        print("<o | o>", ovlp);
    }

    CCPairFunction p1(a,b);
    CCPairFunction p2(f12,a,b);
    CCPairFunction p3(f); // outer (f1,f2)

    std::vector<CCPairFunction> vp1({p1});
    std::vector<CCPairFunction> vp2({p2});
    std::vector<CCPairFunction> vp3({p3});

    Projector<double,3> O(o,o);
    QProjector<double,3> Q(world,o,o);
    StrongOrthogonalityProjector<double,3> Q12(world);
    Q12.set_spaces(o);

    // some hardwire test
    {
        // <ab|k ><k |ab>  = <a|k><k|a><b|b>
        auto ak_tensor=matrix_inner(world,a,o);
        double ak=ak_tensor.trace(ak_tensor);
        auto aa=matrix_inner(world,a,a);
        auto bb=matrix_inner(world,b,b);
        Projector<double,3> O1(a,a);
        O.set_particle(0);
        O1.set_particle(0);
        Projector<double,3> O2(a,a);
        O2.set_particle(0);
        double n1=inner(vp1,O1(vp2));
        double n1a=inner(O1(vp1),vp2);
        print("n1,n1a",n1,n1a);
        double n2=inner(vp1,O2(vp2));
        double n2a=inner(O2(vp1),vp2);
        print("n2,n2a",n1,n1a);

    }
    timer.reset();


    // test {O,Q} with particles 1,2 on all three CCPairFunction types
    {
        int i=0;
        for (auto p : {p1,p2,p3}) {
            std::string s="CCPairFunction p"+std::to_string(++i);
            std::vector<CCPairFunction> vp({p});
            for (int particle : {0,1}) {
                O.set_particle(particle);
                Q.set_particle(particle);
                auto Op = O(vp);
                auto Qp = Q(vp);

                double n1 = inner(Op, vp3) + inner(Qp, vp3);
                double n2 = inner(vp, vp3);
                print("n1 (RI), n2", n1, n2, fabs(n1 - n2));
                t1.checkpoint(fabs(n1 - n2) < FunctionDefaults<3>::get_thresh(), "RI with particle "+std::to_string(particle)+" on "+s ,timer.reset());

                auto Op1 = O(Op);
                double n3=inner(Op, vp3);
                double n4=inner(Op1, vp3);
                print("n3, n4", n3, n4);
                t1.checkpoint(fabs(n3-n4) < FunctionDefaults<3>::get_thresh(), "idempotency with particle "+std::to_string(particle)+" on "+s,timer.reset());

            }
            // testing SO projector
            auto SOp=Q12(vp);
            double n1,n2,n3;
            {
                Q.set_particle(0);
                auto tmp = Q(vp);
                Q.set_particle(1);
                auto tmp1 = Q(tmp);
                n2=inner(vp3,tmp1);
            }
            {
                Q.set_particle(1);
                auto tmp = Q(vp);
                Q.set_particle(0);
                auto tmp1 = Q(tmp);
                n3=inner(vp3,tmp1);
            }
            n1 = inner(SOp, vp3);
            print("SO: n1, n2, n3", n1, n2, n3, fabs(n1 - n2));
            double zero=fabs(n1-n2) + fabs(n1-n3) + fabs(n2-n3);
            t1.checkpoint(zero < FunctionDefaults<3>::get_thresh(), "SO operator on "+s,timer.reset() );
        }
    }
    t1.end();

    return  (t1.get_final_success()) ? 0 : 1;
}

int test_dirac_convolution(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                           const CCParameters& parameter) {
    test_output t1("CCPairFunction::test_dirac_convolution");

    return  (t1.get_final_success()) ? 0 : 1;
}

/// testing <ij | g Q f | ij> = 0.032 mEh
int test_helium(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, const Molecule& molecule,
                   const CCParameters& parameters) {

    test_output t1("CCPairFunction::test_helium");
    CCTimer timer(world, "testing");

    real_function_3d Vnuc = real_factory_3d(world).f([](const coord_3d& r) {return -2.0/(r.normf()+1.e-8);});
    real_function_3d psi  = real_factory_3d(world).f([](const coord_3d& r) {return exp(-r.normf());});

    auto iterate=[&world](const real_function_3d& potential, real_function_3d& psi, double& eps) {
        real_convolution_3d op = BSHOperator3D(world, sqrt(-2*eps), 0.001, 1e-6);
        real_function_3d Vpsi = (potential*psi);
        Vpsi.scale(-2.0).truncate();
        real_function_3d tmp = op(Vpsi).truncate();
        double norm = tmp.norm2();
        real_function_3d r = tmp-psi;
        double rnorm = r.norm2();
        double eps_new = eps - 0.5*inner(Vpsi,r)/(norm*norm);
        if (world.rank() == 0) {
            print("norm=",norm," eps=",eps," err(psi)=",rnorm," err(eps)=",eps_new-eps);
        }
        psi = tmp.scale(1.0/norm);
        eps = eps_new;
    };
    psi.truncate();
    psi.scale(1.0/psi.norm2());
    double eps = -0.6;
    real_convolution_3d op = CoulombOperator(world, 0.001, 1e-6);
    for (int iter=0; iter<10; iter++) {
        real_function_3d rho = square(psi).truncate();
        real_function_3d potential = Vnuc + op(rho).truncate();
        iterate(potential, psi, eps);
    }
    double kinetic_energy = 0.0;
    for (int axis=0; axis<3; axis++) {
        real_derivative_3d D = free_space_derivative<double,3>(world, axis);
        real_function_3d dpsi = D(psi);
        kinetic_energy += inner(dpsi,dpsi);
    }
    real_function_3d rho = square(psi);
    double two_electron_energy = inner(op(rho),rho);
    double nuclear_attraction_energy = 2.0*inner(Vnuc,rho);
    double total_energy = kinetic_energy + two_electron_energy + nuclear_attraction_energy;

    t1.checkpoint(fabs(total_energy+2.8616533)<1.e-4,"helium iterations",timer.reset());
    print("ke, total", kinetic_energy, total_energy);


    CCConvolutionOperator::Parameters param;
    auto f=CCConvolutionOperatorPtr(world,OT_F12,param);
    auto g=CCConvolutionOperatorPtr(world,OT_G12,param);

    CCPairFunction fij(f,psi,psi);
    CCPairFunction gij(g,psi,psi);
    CCPairFunction ij({psi},{psi});
    std::vector<CCPairFunction> vfij={fij};
    std::vector<CCPairFunction> vgij={gij};
    std::vector<CCPairFunction> vij={ij};

    StrongOrthogonalityProjector<double,3> SO(world);
    SO.set_spaces({psi},{psi},{psi},{psi});
    std::vector<CCPairFunction> Qfij=SO(vfij);
    std::vector<CCPairFunction> Qgij=SO(vgij);

    double result1=inner(vgij,Qfij);
    print("<ij | g (Q f | ij>)",result1);
    double result2=inner(vfij,Qgij);
    print("(<ij | g Q) f | ij>",result2);
    bool good=fabs(result1-result2)<1.e-5;
    good=good and fabs(result1+3.2624783e-02)<1.e-5;
    t1.checkpoint(good,"V matrix element",timer.reset());

    return  (t1.get_final_success()) ? 0 : 1;

}

/** functionality
 *
 *  - ctor                      OK
 *  - assignment                OK
 *  - add
 *  - scalar multiplication     OK
 *  - inner                     OK
 *  - inner_partial             OK
 *  - swap_particles            OK
 *  - apply
 *  - apply_partial (i.e. exchange)
 *  - serialize
 *  - collapse_to_pure (excl g!)
 *  - mul_partial
 */

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
    FunctionDefaults<6>::set_thresh(1.e-4);

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
    FunctionDefaults<6>::set_cubic_cell(-10.,10.);

    FunctionDefaults<2>::set_tensor_type(TT_FULL);
    print("numerical parameters: k, eps(3D), eps(6D)", FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(),
          FunctionDefaults<6>::get_thresh());
    XParameters parameters;
    parameters.read_and_set_derived_values(world,parser,"grid");
    parameters.print("grid parameters");
    int isuccess=0;
#ifdef USE_GENTENSOR

    try {
        parser.set_keyval("geometry", "he");
        parser.print_map();
        Molecule mol(world, parser);
        mol.print();
        CCParameters ccparam(world, parser);

//        data1=data(world,ccparam);

        std::shared_ptr<NuclearCorrelationFactor> ncf = create_nuclear_correlation_factor(world,
                         mol, nullptr, std::make_pair("slater", 2.0));

        isuccess+=test_lowrank_function3(world,parameters);
//        isuccess+=test_constructor(world, ncf, mol, ccparam);
//        isuccess+=test_operator_apply(world, ncf, mol, ccparam);
//        isuccess+=test_transformations(world, ncf, mol, ccparam);
//        isuccess+=test_inner(world, ncf, mol, ccparam);
//        isuccess+=test_multiply(world, ncf, mol, ccparam);
//        isuccess+=test_multiply_with_f12(world, ncf, mol, ccparam);
//        isuccess+=test_swap_particles(world, ncf, mol, ccparam);
//        isuccess+=test_scalar_multiplication(world, ncf, mol, ccparam);
//        isuccess+=test_partial_inner_3d(world, ncf, mol, ccparam);
//        isuccess+=test_partial_inner_6d(world, ncf, mol, ccparam);
//        isuccess+=test_projector(world, ncf, mol, ccparam);
//        FunctionDefaults<3>::set_cubic_cell(-10,10);
//        isuccess+=test_helium(world,ncf,mol,ccparam);
        data1.clear();
    } catch (std::exception& e) {
        madness::print("an error occured");
        madness::print(e.what());
        data1.clear();
    }
#else
    print("could not run test_ccpairfunction: U need to compile with ENABLE_GENTENSOR=1");
#endif
    finalize();

    return isuccess;
}