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

#include<madness/world/test_utilities.h>
#include<random>

using namespace madness;


bool longtest=false;
template<typename T, std::size_t NDIM>
struct data {
    static constexpr std::size_t LDIM=NDIM/2;
    Function<T,LDIM> f1,f2,f3,f4,f5;
    Function<T,NDIM> f12,f23;
    World& world;
    CCParameters parameters;

    std::shared_ptr<CCConvolutionOperator<T,LDIM>> f12_op;

    data(World& world, const CCParameters& parameters) : world(world), parameters(parameters){}

    bool is_initialized() const {
        return f1.is_initialized();
    }

    void initialize() {
        print("initializing data for NDIM=",NDIM);
        f12_op.reset(new CCConvolutionOperator<T,LDIM>(world,OT_F12,parameters));
        auto g1 = [](const Vector<double,LDIM>& r) { return exp(-1.0 * inner(r, r)); };
        auto g2 = [](const Vector<double,LDIM>& r) { return exp(-2.0 * inner(r, r)); };
        auto g3 = [](const Vector<double,LDIM>& r) { return exp(-3.0 * inner(r, r)); };
        auto g4 = [](const Vector<double,LDIM>& r) { return exp(-4.0 * inner(r, r)); };
        auto g5 = [](const Vector<double,LDIM>& r) { return exp(-5.0 * inner(r, r)); };
        f1=FunctionFactory<T,LDIM>(world).f(g1);
        f2=FunctionFactory<T,LDIM>(world).f(g2);
        f3=FunctionFactory<T,LDIM>(world).f(g3);
        f4=FunctionFactory<T,LDIM>(world).f(g4);
        f5=FunctionFactory<T,LDIM>(world).f(g5);

        auto g = [](const Vector<double,NDIM>& r) {
            double r1=0.0, r2=0.0;
            for (int i=0; i<LDIM; ++i) {
                r1+=r[i]*r[i];
                r2+=r[i+LDIM]*r[i+LDIM];
            }
            return exp(-1.0*r1 - 2.0*r2);
        };
        auto g23 = [](const Vector<double,NDIM>& r) {
            double r1=0.0, r2=0.0;
            for (int i=0; i<LDIM; ++i) {
                r1+=r[i]*r[i];
                r2+=r[i+LDIM]*r[i+LDIM];
            }
            return exp(-1.0*r1 - 2.0*r2) + exp(-2.0*r1 - 3.0*r2);
        };

        f12=FunctionFactory<T,NDIM>(world);
        f23=FunctionFactory<T,NDIM>(world);
        std::string name_f12="test_ccpairfunction_f12_ndim_"+std::to_string(NDIM);
        std::string name_f23="test_ccpairfunction_f23_ndim_"+std::to_string(NDIM);
//        try {
//            load(f12,name_f12);
//        } catch (...) {
            f12 = FunctionFactory<T,NDIM>(world).f(g);
//            save(f12,name_f12);
//        }
//        try {
//            load(f23,name_f23);
//        } catch (...) {
            f23 = FunctionFactory<T,NDIM>(world).f(g23);
//            save(f23,name_f23);
//        }

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
    auto get_functions() {
        if (not is_initialized()) initialize();
        return std::make_tuple(f1,f2,f3,f4,f5,f12);
    }

    /// get some standard ccpairfunctions

    /// p1: pure, corresponds to f12
    /// p2: dec, corresponds to f23
    /// p3: op_dec, corresponds to f23
    /// p4: pure, corresponds to f23
    /// p5: op_pure, corresponds to f23
    auto get_ccpairfunctions() {
        if (not is_initialized()) initialize();
        CCPairFunction<T,NDIM> p1(copy(f12));
        CCPairFunction<T,NDIM> p2({f1,f2},{f2,f3});
        CCPairFunction<T,NDIM> p3(f12_op,{f1,f2},{f2,f3});
        CCPairFunction<T,NDIM> p4(copy(f23)); // two-term, corresponds to p2
        CCPairFunction<T,NDIM> p5(f12_op,copy(f23)); // two-term, corresponds to p2
        return std::make_tuple(p1,p2,p3,p4,p5);
    }

};


template<typename T, std::size_t NDIM>
int test_constructor(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                     const CCParameters& parameter) {
    test_output t1("constructor<T,"+std::to_string(NDIM)+">");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;

    Function<T,NDIM> f=FunctionFactory<T,NDIM>(world);
    auto [f1,f2,f3,f4,f5,ff]=data.get_functions();

    std::vector<Function<T,LDIM>>  a= zero_functions<T,LDIM>(world,LDIM);
    std::vector<Function<T,LDIM>>  b= zero_functions<T,LDIM>(world,LDIM);
    auto f12=CCConvolutionOperatorPtr<double,LDIM>(world, OT_F12, parameter);
    t1.checkpoint(true,"preparation");

    auto f_copy=copy(f);
    CCPairFunction<T,NDIM> p1;
    CCPairFunction<T,NDIM> p2(f_copy);
    CCPairFunction<T,NDIM> p3({f1,f2},{f1,f3});
    CCPairFunction<T,NDIM> p4(f12,{f1,f2},{f2,f3});
    t1.checkpoint(true,"construction");

    {
        MADNESS_CHECK(p2.is_pure());
        MADNESS_CHECK(!p2.is_decomposed());
        MADNESS_CHECK(!p2.is_op_decomposed());
        auto p = p2.pure();
        auto ff = p2.pure().get_function();
        MADNESS_CHECK((ff.get_impl()==f_copy.get_impl())); // shallow copy of f
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
        CCPairFunction<T,NDIM> tmp;
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

template<typename T, std::size_t NDIM>
int test_norm(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                     const CCParameters& parameter) {
    test_output t1("norm of <T,"+std::to_string(NDIM)+">");

    auto [p1,p2,p3,p4,p5]=data.get_ccpairfunctions();  // p2-p5 correspond to f230
    for (const CCPairFunction<T,NDIM>& p : {p2,p3,p4,p5}) {
        double n=p.norm2();
        print("norm of ",p.name(),n);
        double n1=sqrt(inner(p,p));
        print("inner",n1);
        t1.checkpoint(n,n1,FunctionDefaults<NDIM>::get_thresh(),"norm of p");
    }

    double n2=p2.norm2();
    double n3=p3.norm2();
    double n4=p4.norm2();
    double n5=p5.norm2();
    t1.checkpoint(n2,n4,FunctionDefaults<NDIM>::get_thresh(),"norm of p2/4");
    t1.checkpoint(n3,n5,FunctionDefaults<NDIM>::get_thresh(),"norm of p3/5");

    return t1.end();

}

template<typename T, std::size_t NDIM>
int test_load_store(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                     const CCParameters& parameter) {
    test_output t1("load/store of <T,"+std::to_string(NDIM)+">");

//    t1.set_cout_to_terminal();
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;

    auto [p1,p2,p3,p4,p5]=data.get_ccpairfunctions();  // p2-p5 correspond to f23
    auto [f1,f2,f3,f4,f5,ff]=data.get_functions();

    auto compute_diff_norm = [](const CCPairFunction<T,NDIM>& f1, const CCPairFunction<T,NDIM> f2) {
        std::vector<CCPairFunction<T,NDIM>> diff;
        diff+={f1};
        diff-={f2};
        double in = inner(diff,diff);
        if (in<0) return -sqrt(-in);
        return sqrt(in);
    };

    std::string fname="ccpairfunction_test";
    {
        archive::ParallelOutputArchive<archive::BinaryFstreamOutputArchive> ar(world, fname, 1);
        ar & f1;
        ar & p2;
        ar & p3;
        ar & p4;
        ar & p5;
    }
    {
        archive::ParallelInputArchive<archive::BinaryFstreamInputArchive> ar(world, fname, 1);
        CCPairFunction<T,NDIM> g2, g3, g4, g5;
        ar & f1;
        ar & g2;
        ar & g3;
        ar & g4;
        ar & g5;

        double n2=compute_diff_norm(g2,p2);
        double n3=compute_diff_norm(g3,p3);
        double n4=compute_diff_norm(g4,p4);
        double n5=compute_diff_norm(g5,p5);
        t1.checkpoint(n2,FunctionDefaults<LDIM>::get_thresh(),"store/load "+p2.name());
        t1.checkpoint(n3,FunctionDefaults<LDIM>::get_thresh(),"store/load "+p3.name());
        t1.checkpoint(n4,FunctionDefaults<LDIM>::get_thresh(),"store/load "+p4.name());
        t1.checkpoint(n5,FunctionDefaults<LDIM>::get_thresh(),"store/load "+p5.name());
    }
    t1.checkpoint(true,"checks on load/store");
    return t1.end();
}


template<typename T, std::size_t NDIM>
int test_operator_apply(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                        const CCParameters& parameter) {
    test_output t1("test_operator_apply<double,"+std::to_string(NDIM)+">");
//    t1.set_cout_to_terminal();

    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;

    double exponent=1.0; // corresponds to the exponent of data::f1 and data::ff
//    double coefficient=pow(1.0/constants::pi*exponent,0.5*3);
    double coefficient=1.0;
    const Vector<double,LDIM> center(0.0);

    auto Gaussian = [&center, &exponent, &coefficient](const Vector<double,LDIM>& r) {
        return coefficient * exp(-exponent*inner(r-center,r-center));
    };

    exponent=1.5;
    Function<T,LDIM> i=FunctionFactory<T,LDIM>(world).functor(Gaussian);
    exponent=2.5;
    Function<T,LDIM> j=FunctionFactory<T,LDIM>(world).functor(Gaussian);

    auto gop= BSHOperator<LDIM>(world,0.5,1.e-5,std::min(1.e-4,FunctionDefaults<LDIM>::get_thresh()));

    auto [f1,f2,f3,f4,f5,ff]=data.get_functions();
    auto [p1,p2,p3,p4,p5]=data.get_ccpairfunctions();  // p2-p5 correspond to f23
    // p2: dec, corresponds to f23
    // p3: op_dec, corresponds to f23
    // p4: pure, corresponds to f23
    // p5: op_pure, corresponds to f23

    // integrate over particle 2: | u(1,3)> = g(1,2) uu(2,3), project u(1,3) onto bra <i(1) j(3) |
    // reference for p2, p4:
    // <i(1) j(3) | g(1,2) f(2,3) | k(2) l(3)> = <g_i(2) | k(2) f_lj(2)>
    // reference for p1, p3:
    // <i(1) j(3) | g(1,2) | k(2) l(3)> = <g_i(2) *j(3) | k(2) l(3)> = <g_i(2) | k(2)> <j(3)|l(3)>
    //
    // note that the ccpairfunctions always have 2 terms
    // CCPairFunction<T,NDIM> p2({f1,f2},{f2,f3});
    // with l=f1, k=f2; and l=f2, k=f3

    auto& fop=*(data.f12_op->get_op());
    Function<T,LDIM> gi=gop(i);

    Function<T,LDIM> k_1=f1;
    Function<T,LDIM> l_1=f2;
    Function<T,LDIM> f_lj_1=fop(l_1*j);
    Function<T,LDIM> f_kj_1=fop(k_1*j);

    Function<T,LDIM> k_2=f2;
    Function<T,LDIM> l_2=f3;
    Function<T,LDIM> f_lj_2=fop(l_2*j);
    Function<T,LDIM> f_kj_2=fop(k_2*j);

    CCPairFunction<T,NDIM> bra(i,j);

    double ref_p2p4=(gi*k_1*f_lj_1).trace() + (gi*k_2*f_lj_2).trace();
    double ref_p1p3=(gi*k_1).trace() * (j*l_1).trace() + (gi*k_2).trace() * (j*l_2).trace();
    std::vector<double> reference({ref_p1p3,ref_p2p4,ref_p1p3,ref_p2p4});

    int counter=0;

    for (auto& p : {p2,p3,p4,p5}) {
        gop.particle() = 1;
        auto tmp1=gop(p);
        double result=inner(bra,tmp1);
        double ref=reference[counter++];
        print("p=",p.name(),"result=",result,"reference=",ref);
        t1.checkpoint(result,ref,FunctionDefaults<LDIM>::get_thresh(),"op(1) to "+p.name());
    }

    // integrate over particle 3: | u(1,2)> = g(1,3) uu(2,3), project u(1,2) onto bra <i(1) j(2) |
    // reference for p2, p4:
    // <i(1) j(2) | g(1,3) f(2,3) | k(2) l(3)> = <g_i(3) | l(3) f_kj(3)>
    // reference for p1, p3:
    // <i(1) j(2) | g(1,3) | k(2) l(3)> = <g_i(3) j(2) | k(2) l(3)> = <j(2) | k(2)> <g_i(3)|l(3)>
    // note the ordering of the indices in the ket! the above line is equivalent to
    // <i(1) j(2) | g(1,3) | k(2) l(3)> = <i(1) j(2) | g_l(1) k(2)>
    ref_p2p4=(gi*l_1*f_kj_1).trace() + (gi*l_2*f_kj_2).trace();
    ref_p1p3=(gi*l_1).trace() * (j*k_1).trace() + (gi*l_2).trace() * (j*k_2).trace();
    reference=std::vector<double>({ref_p1p3,ref_p2p4,ref_p1p3,ref_p2p4});

    counter=0;
    CCPairFunction<T,NDIM> bra_ji(j,i);
    for (auto& p : {p2,p3,p4,p5}) {
        gop.particle() = 2;
        auto tmp1=gop(p);
        double result=inner(bra,tmp1);
        double ref=reference[counter++];
        print("p=",p.name(),"result=",result,"reference=",ref);
        t1.checkpoint(result,ref,FunctionDefaults<LDIM>::get_thresh(),"op(2) to "+p.name());
    }



    return t1.end();
}


template<typename T, std::size_t NDIM>
int test_transformations(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                         const CCParameters& parameter) {
    test_output t1("test_transformations<double,"+std::to_string(NDIM)+">");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;

//    auto data=get_data<NDIM>(world,parameter);
    auto [f1,f2,f3,f4,f5,ff]=data.get_functions();
    auto f12=CCConvolutionOperatorPtr<double,LDIM>(world, OT_F12, parameter);
    auto g12=CCConvolutionOperatorPtr<double,LDIM>(world, OT_G12, parameter);

    auto compute_diff_norm = [](const CCPairFunction<T,NDIM>& f1, const CCPairFunction<T,NDIM> f2) {
        std::vector<CCPairFunction<T,NDIM>> diff;
        diff+={f1};
        diff-={f2};
        return sqrt(inner(diff,diff));
    };

    CCPairFunction<T,NDIM> p1(copy(ff));
    t1.checkpoint(p1.is_pure(),"is_pure");
    t1.checkpoint(p1.is_pure_no_op(),"is_pure_no_op");

    CCPairFunction<T,NDIM> p2(f12,ff);
    t1.checkpoint(p2.is_pure(),"is_pure");
    t1.checkpoint(p2.is_op_pure(),"is_op_pure");
    t1.checkpoint(p2.is_convertible_to_pure_no_op(),"is_convertible_to_pure_no_op");
    CCPairFunction<T,NDIM> p3=copy(p2);
    p3.convert_to_pure_no_op_inplace();
    t1.checkpoint(p2.is_op_pure(),"is_op_pure");
    t1.checkpoint(p3.is_pure_no_op(),"is_pure_no_op");

    CCPairFunction<T,NDIM> p4(g12,copy(ff));
    t1.checkpoint(p4.is_pure(),"is_pure");
    t1.checkpoint(p4.is_op_pure(),"is_op_pure");
    t1.checkpoint(not p4.is_convertible_to_pure_no_op(),"not is_convertible_to_pure_no_op");

    // convert f12 f1 f2 to pure_op_op
    CCPairFunction<T,NDIM> p5(f12,f1,f2);
    t1.checkpoint(not p5.is_pure(),"is_pure");
    t1.checkpoint(p5.is_op_decomposed(),"is_op_decomposed");
    t1.checkpoint(p5.is_convertible_to_pure_no_op(),"is_convertible_to_pure_no_op");
    CCPairFunction<T,NDIM> p6=copy(p5);
    p6.convert_to_pure_no_op_inplace();
    t1.checkpoint(p6.is_pure_no_op(),"is_pure_no_op");
    double d6=compute_diff_norm(p5,p6);
    t1.checkpoint(d6,FunctionDefaults<NDIM>::get_thresh()*50,"numerics");

    // convert \sum_i f12 f1_i f2_i to pure_op_op
    CCPairFunction<T,NDIM> p7(f12,{f1,f2,f3},{f1,f2,f3});
    t1.checkpoint(not p7.is_pure(),"is_pure");
    t1.checkpoint(p7.is_op_decomposed(),"is_op_decomposed");
    t1.checkpoint(p7.is_convertible_to_pure_no_op(),"is_convertible_to_pure_no_op");
    CCPairFunction<T,NDIM> p8=copy(p7);
    p8.convert_to_pure_no_op_inplace();
    t1.checkpoint(p8.is_pure_no_op(),"is_pure_no_op");
    double d8=compute_diff_norm(p7,p8);
    t1.checkpoint(d8,FunctionDefaults<NDIM>::get_thresh()*50,"numerics");



    return t1.end();
}

template<typename T, std::size_t NDIM>
int test_multiply_with_f12(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                           const CCParameters& parameters) {
    test_output t1("test_multiply_with_f12<double,"+std::to_string(NDIM)+">");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;

    // p1: pure, corresponds to f12
    // p2: dec, corresponds to f23
    // p3: op_dec, corresponds to f23
    // p4: pure, corresponds to f23
    // p5: op_pure, corresponds to f23
//    auto data=get_data<NDIM>(world,parameters);
    auto [p1,p2,p3,p4,p5]=data.get_ccpairfunctions();  // p2-p5 correspond to f23
    auto f12=data.f12_op;

    double thresh=FunctionDefaults<LDIM>::get_thresh();

    // decomposed
    CCPairFunction<double,NDIM> tmp1=f12*p2;         // should now be identical to p3
    CCPairFunction<double,NDIM> tmp2=p2*f12;         // should now be identical to p3
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

template<typename T, std::size_t NDIM>
int test_multiply(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                  const CCParameters& parameters) {
    test_output t1("test_multiply<"+std::to_string(NDIM)+">");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;

    // consistency check, relies on CCPairFunction<double,6>::inner to work correctly
    double thresh=FunctionDefaults<LDIM>::get_thresh();
//    auto data=get_data<NDIM>(world,parameters);
    auto [p1,p2,p3,p4,p5]=data.get_ccpairfunctions();  // p2-p5 correspond to f23
    auto [f1,f2,f3,f4,f5,f]=data.get_functions();

    auto particle1=std::array<int,LDIM>();
    auto particle2=std::array<int,LDIM>();
    for (int i=0; i<LDIM; ++i) {
        particle1[i]=i;
        particle2[i]=i+LDIM;
    }

    // reference value is <bra | f(1/2) p>  = <f1 f | p>
    CCPairFunction<T,NDIM> bra(f1,f2);
    CCPairFunction<T,NDIM> bra1(f1*f2,f2);
    CCPairFunction<T,NDIM> bra2(f1,f2*f2);
    for (auto& p : {p2,p3,p4,p5}) {

        auto tmp1=multiply(p,f2,particle1);
        double ovlp1=inner(bra,tmp1);
        double ref1=p.has_operator() ? inner(bra1,p3) : inner(bra1,p2);

        bool good=(fabs(ovlp1-ref1)<thresh);
        t1.checkpoint(good,"f(1)*"+p.name());

        auto tmp2=multiply(p,f2,particle2);
        double ovlp2=inner(bra,tmp2);
        double ref2=p.has_operator() ? inner(bra2,p3) : inner(bra2,p2);

        good=(fabs(ovlp2-ref2)<thresh);
        t1.checkpoint(good,"f(2)*"+p.name());
    }

    return t1.end();
}

template<typename T, std::size_t NDIM>
int test_inner(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
               const CCParameters& parameters) {
    test_output t1("test_inner<"+std::to_string(NDIM)+">");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;
//    t1.set_cout_to_terminal();

    /// f1: exp(-1.0 r^2)
    /// f2: exp(-2.0 r^2)
    /// f3: exp(-3.0 r^2)
    /// f4: exp(-4.0 r^2)
    /// f5: exp(-5.0 r^2)
    /// f: exp(-r_1^2 - 2 r_2^2)
    /// f23: exp(-r_1^2 - 2 r_2^2) + exp(-2 r_1^2 - 3 r_2^2)
//    auto data=get_data<NDIM>(world,parameters);
    auto [f1,f2,f3,f4,f5,f]=data.get_functions();
    /// p1: pure, corresponds to f12
    /// p2: dec, corresponds to f23
    /// p3: op_dec, corresponds to f23
    /// p4: pure, corresponds to f23
    /// p5: op_pure, corresponds to f23
    auto [p1,p2,p3,p4,p5]=data.get_ccpairfunctions();
    auto f12 = *(data.f12_op);

    /// results
    auto a=std::vector<Function<T,LDIM>>({f1,f2});
    auto b=std::vector<Function<T,LDIM>>({f2,f3});
    std::vector<Function<T,LDIM>> a_ij_functions, b_ij_functions;
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
    SeparatedConvolution<double, LDIM> fop= SlaterOperator<LDIM>(world, y, parameters.lo(), parameters.thresh_bsh_3D());
    SeparatedConvolution<double, LDIM> fsq = SlaterOperator<LDIM>(world, 2.0 * y, parameters.lo(), parameters.thresh_bsh_3D());

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
            double thresh=FunctionDefaults<LDIM>::get_thresh();
            t1.checkpoint(result,ref,thresh,bra.name(true)+ket.name());
        }
    }
    return t1.end();
}


template<typename T, std::size_t NDIM>
int test_partial_inner_6d(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                          const CCParameters& parameter) {

    test_output t1("test_partial_inner6d<"+std::to_string(NDIM)+">");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;
//    t1.set_cout_to_terminal();
//    auto data=get_data<6>(world,parameter);
    auto [f1,f2,f3,f4,f5,f] = data.get_functions();

    std::vector<Function<T,LDIM>> a = {f1, f2};
    std::vector<Function<T,LDIM>> b = {f2, f3};

    auto f12=CCConvolutionOperatorPtr<double,LDIM>(world, OT_F12, parameter);

    auto [p1,p2,p3,p4,nil]=data.get_ccpairfunctions();
    CCPairFunction<T,NDIM> p11({f1},{f1});
    CCPairFunction<T,NDIM> p12({f1},{f2});

    double g11=inner(f1,f1);
    double g22=inner(f2,f2);
    double g12=inner(f1,f2);
    double g13=inner(f1,f3);
    double g23=inner(f2,f3);
    double g33=inner(f3,f3);
    Function<T,LDIM> gf11=(*f12)(f1*f1);
    Function<T,LDIM> gf12=(*f12)(f1*f2);
    Function<T,LDIM> gf22=(*f12)(f2*f2);
    Function<T,LDIM> gf13=(*f12)(f1*f3);
    Function<T,LDIM> gf23=(*f12)(f2*f3);

    t1.checkpoint(true,"prep");

    auto particle1=std::array<int,LDIM>();
    auto particle2=std::array<int,LDIM>();
    for (int i=0; i<LDIM; ++i) {
        particle1[i]=i;
        particle2[i]=i+LDIM;
    }

    // p1 = p12 = e(-r1 -2r2)
    // p2 = e(-r1) * e(-2r2) + e(-2r1) * e(-3e2)  separated
    // p4 = e(-r1) * e(-2r2) + e(-2r1) * e(-3e2)  6d
    for (auto test_p1 : {p2,p4}) {
        for (auto test_p2 : {p2,p4}) {
            CCPairFunction<T,NDIM> r1=inner(test_p1,test_p2,particle1,particle1);
            CCPairFunction<T,NDIM> r2=inner(test_p1,test_p2,particle1,particle2);
            CCPairFunction<T,NDIM> r3=inner(test_p1,test_p2,particle2,particle1);
            CCPairFunction<T,NDIM> r4=inner(test_p1,test_p2,particle2,particle2);

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
            t1.checkpoint(good,test_p1.name(true)+test_p2.name()+" -- 4");

        }
    }

    // test < sth | f(1,2) a(1)b(2) >
    // CCPairFunction<double,6> p2({f1,f2},{f2,f3});
    // CCPairFunction<double,6> p4(f23); // two-term, corresponds to p2
    CCPairFunction<T,NDIM> p5(f12,std::vector<Function<T,LDIM>>({f1}),std::vector<Function<T,LDIM>>({f2}));
    for (auto& test_p1 : {p2, p4}) {
        CCPairFunction<T,NDIM> r1=inner(test_p1,p5,particle1,particle1);
        CCPairFunction<T,NDIM> r2=inner(test_p1,p5,particle1,particle2);
        CCPairFunction<T,NDIM> r3=inner(test_p1,p5,particle2,particle1);
        CCPairFunction<T,NDIM> r4=inner(test_p1,p5,particle2,particle2);

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
        t1.checkpoint(good,test_p1.name(true)+p5.name()+" -- 4");

    }

    // test < a(1) b(2) f(1,2) | f(1,3) c(1)d(3) >
    // CCPairFunction<double,6> p3(f12_op.get(),{f1,f2},{f2,f3});
    // CCPairFunction<double,6> p5(&f12,{f1},{f2});
    if (longtest) {
        CCPairFunction<T,NDIM> r1=inner(p3,p5,particle1,particle1);
        double n1=inner(r1,p11);
        p3.convert_to_pure_no_op_inplace();
        p5.convert_to_pure_no_op_inplace();
        print("n1",n1);
        CCPairFunction<T,NDIM> r1a=inner(p3,p5,particle1,particle1);
        double n1a=inner(r1a,p11);
        print("n1a",n1a);
        print("diff",n1-n1a);
        CCPairFunction<T,NDIM> r2=inner(p3,p5,particle1,particle2);
        CCPairFunction<T,NDIM> r3=inner(p3,p5,particle2,particle1);
        CCPairFunction<T,NDIM> r4=inner(p3,p5,particle2,particle2);

    }
    return t1.end();
}


template<typename T, std::size_t NDIM>
int test_partial_inner_3d(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                          const CCParameters& parameter) {

    test_output t1("test_partial_inner<"+std::to_string(NDIM)+">");
//    t1.set_cout_to_terminal();
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;
    auto [f1,f2,f3,f4,f5,f] = data.get_functions();

    auto particle1=std::array<int,LDIM>();
    auto particle2=std::array<int,LDIM>();
    for (int i=0; i<LDIM; ++i) {
        particle1[i]=i;
        particle2[i]=i+LDIM;
    }

    std::vector<Function<T,LDIM>> a = {f1, f2};
    std::vector<Function<T,LDIM>> b = {f2, f3};

    auto f12=CCConvolutionOperatorPtr<double,LDIM>(world, OT_F12, parameter);

    CCPairFunction<T,NDIM> p1(copy(f));   // e(-r1 - 2r2)
    CCPairFunction<T,NDIM> p2(a,b);
    CCPairFunction<T,NDIM> p3(f12,a,b);
    CCPairFunction<T,NDIM> p11({f1},{f1});
    CCPairFunction<T,NDIM> p12({f1},{f2});

    double g11=inner(f1,f1);
    double g22=inner(f2,f2);
    double g12=inner(f1,f2);
    double g13=inner(f1,f3);
    double g23=inner(f2,f3);
    Function<T,LDIM> gf11=(*f12)(f1*f1);
    Function<T,LDIM> gf12=(*f12)(f1*f2);
    Function<T,LDIM> gf13=(*f12)(f1*f3);
    print("g11, g22",g11,g22);

    double thresh=FunctionDefaults<LDIM>::get_thresh();

    t1.checkpoint(true,"prep");

    // test pure/3d
    {
        double aa=inner(p1,p1);
        print("aa",aa);
        // < e1 e2 | e1 >_1 = | e2 > ||e1||
        Function<T,LDIM> r=inner(p1,f1,particle1,particle1);
        double norm=inner(r,f2);
        double ref_norm=g11 * g22;
        print("norm    ",norm);
        print("ref_norm", ref_norm);
        Function<T,LDIM> r_swap=inner(p1,f1,particle2,particle1);
        double norm_swap=inner(r_swap,f2);
        print("norm1 swap",norm_swap);
        t1.checkpoint(norm,ref_norm,thresh,"pure -- 1");
    }
    // test pure/3d
    {
        Function<T,LDIM> r=inner(p1,f1,particle2,particle1);
        double norm=inner(r,f2);
        double ref_norm=g12 * g12;
        t1.checkpoint(norm, ref_norm, thresh,"pure -- 2");
    }
    // test decomposed
    {
        Function<T,LDIM> r=inner(p2,f1,particle1,particle1);
        double norm=inner(r,f2);
        double ref_norm=g11 * g22 + g12 * g23;
        t1.checkpoint(norm, ref_norm, thresh,"decomposed -- 1");
    }
    // test decomposed
    {
        Function<T,LDIM> r=inner(p2,f1,particle2,particle1);
        double norm=inner(r,f2);
        double ref_norm=g12 * g12 + g13 * g22;
        t1.checkpoint(norm, ref_norm, thresh,"decomposed -- 2");
    }
    // test op_decomposed
    {
        // < f1 f2 | f | f1 f2 > + < f1 f2 | f | f2 f3>
        Function<T,LDIM> r=inner(p3,f1,particle1,particle1);
        double norm=inner(r,f2);
        double ref_norm=inner(gf11*f2,f2) + inner(gf12*f2,f3);
        t1.checkpoint(norm, ref_norm, thresh,"op_decomposed -- 1");
    }
    // test op_decomposed
    {
        // < f1 f2 | f | f2 f1 > + < f1 f2 | f | f3 f2>
        Function<T,LDIM> r=inner(p3,f1,particle2,particle1);
        double norm=inner(r,f2);
        double ref_norm=inner(gf12*f1,f2) + inner(gf13*f2,f2);
        t1.checkpoint(norm, ref_norm, thresh,"op_decomposed -- 2");
    }

    return t1.end();
}


template<typename T, std::size_t NDIM>
int test_consolidate(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
               const CCParameters& parameter) {
    test_output t1("CCPairFunction::test_consolidate");
    static_assert(NDIM % 2 == 0, "NDIM must be even");
    constexpr std::size_t LDIM = NDIM / 2;
//    t1.set_cout_to_terminal();

    /// f12: exp(-r_1^2 - 2 r_2^2)
    /// f23: exp(-r_1^2 - 2 r_2^2) + exp(-2 r_1^2 - 3 r_2^2)
    /// p1: pure, corresponds to f12
    /// p2: dec, corresponds to f23
    /// p3: op_dec, corresponds to f23
    /// p4: pure, corresponds to f23
    /// p5: op_pure, corresponds to f23
    auto [p1,p2,p3,p4,p5]=data.get_ccpairfunctions();

    // collect all terms of similar type, no conversions
    for (const auto& p : {p1,p4,p5}) {
        auto tmp=std::vector<CCPairFunction<T,NDIM>>({p,p});
        double r0=inner(p,{p1});
        double r1=inner(tmp,{p1});
        auto tmp1=consolidate(tmp,{});
        double r2=inner(tmp1,{p1});
        t1.checkpoint(tmp1.size()==1 && tmp.size()==2,"vector size");
        t1.checkpoint(2.0*r0,r1,FunctionDefaults<LDIM>::get_thresh(),"consolidate");
        t1.checkpoint(2.0*r0,r2,FunctionDefaults<LDIM>::get_thresh(),"consolidate");
    }

    // convert op_pure to pure
    for (const auto& p : {p5}) {
        auto tmp=std::vector<CCPairFunction<T,NDIM>>({p});
        t1.checkpoint(tmp.front().is_op_pure(),"correct initial type: op_pure");
        auto tmp1=consolidate(tmp,{"op_pure_to_pure"});
        t1.checkpoint(tmp1.front().is_pure_no_op(),"correct final type: pure");
        t1.checkpoint(is_collected(tmp1),"is_collected");

        double r0=inner(p,{p1});
        double r1=inner(tmp,{p1});
        t1.checkpoint(r0,r1,FunctionDefaults<LDIM>::get_thresh(),"correct numbers");
    }

    // convert op_decomposed to decomposed
    for (const auto& p : {p3}) {
        auto tmp=std::vector<CCPairFunction<T,NDIM>>({p});
        t1.checkpoint(tmp.front().is_op_decomposed(),"correct initial type: op_decomposed");
        auto tmp1=consolidate(tmp,{"op_dec_to_dec"});
        t1.checkpoint(tmp1.front().is_decomposed_no_op(),"correct final type: decomposed");

        double r0=inner(p,{p1});
        double r1=inner(tmp,{p1});
        t1.checkpoint(r0,r1,FunctionDefaults<LDIM>::get_thresh(),"correct numbers");
    }

    // remove linear dependencies
    for (const auto& p : {p3}) {
        auto tmp=std::vector<CCPairFunction<T,NDIM>>({p});
        tmp+=tmp;
        t1.checkpoint(tmp.size()==2,"correct number of terms");
        t1.checkpoint(tmp.front().get_a().size()==2,"correct number of vectors in a");
        t1.checkpoint(tmp.front().is_op_decomposed(),"correct initial type: op_decomposed");
        auto tmp1=consolidate(tmp,{"remove_lindep"});
        t1.checkpoint(tmp1.front().is_decomposed(),"correct final type: decomposed");
        t1.checkpoint(tmp1.size()==1,"correct number of terms");
        t1.checkpoint(tmp1.front().get_a().size()==2,"correct number of vectors in a");

        double r0=2*inner(p,{p1});
        double r1=inner(tmp,{p1});
        t1.checkpoint(r0,r1,FunctionDefaults<LDIM>::get_thresh(),"correct numbers");
    }



    // some random composition of the above
    // a vector of numerically identical terms is created, then a random permutation is applied, and the result is checked
    for (int i=0; i<5; ++i) {
        std::vector<CCPairFunction<T,NDIM>> pvec={p2,p3,p4,p5};     // numerically all the same

        std::random_device rd;  // a seed source for the random number engine
        std::mt19937 gen(rd()); // mersenne_twister_engine seeded with rd()
        std::uniform_int_distribution<> distrib(0, 3);
        std::vector<int> mapping(5);
        for (int i=0; i<5; ++i) mapping[i]=distrib(gen);
        print("mapping",mapping);

        std::vector<CCPairFunction<T,NDIM>> tmp;
        for (auto m: mapping) tmp.push_back(pvec[m]);

        double n0=inner({p1},tmp);

        auto tmp1=consolidate(tmp,{}); // collect only, no conversions
        print("tmp");
        for (auto& c : tmp) c.print_size();
        print("tmp1");
        for (auto& c : tmp1) c.print_size();

        t1.checkpoint(is_collected(tmp1),"random is_collected");
        double n1=inner({p1},tmp1);
        t1.checkpoint(n0,n1,FunctionDefaults<LDIM>::get_thresh(),"random collect numerics");

        tmp1=consolidate(tmp1,{"op_pure_to_pure"}); //
        print("tmp1 after op_pure_to_pure");
        for (auto& c : tmp1) c.print_size();
        t1.checkpoint(is_collected(tmp1),"random is_collected");
        double n2=inner({p1},tmp1);
        t1.checkpoint(n0,n2,FunctionDefaults<LDIM>::get_thresh(),"random op_pure_to_pure numerics");

        tmp1=consolidate(tmp1,{"op_dec_to_dec"}); //
        print("tmp1 after op_dec_to_dec");
        for (auto& c : tmp1) c.print_size();
        t1.checkpoint(is_collected(tmp1),"random is_collected");
        double n3=inner({p1},tmp1);
        t1.checkpoint(n0,n3,FunctionDefaults<LDIM>::get_thresh(),"random op_dec_to_dec numerics");
        t1.checkpoint(tmp1.size()<=2,"only max. two types of terms after consolidate");
    }

    return t1.end();
}


template<typename T, std::size_t NDIM>
int test_apply(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
               const CCParameters& parameter) {
    test_output t1("CCPairFunction::test_apply");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;
//    t1.set_cout_to_terminal();

    /// f12: exp(-r_1^2 - 2 r_2^2)
    /// f23: exp(-r_1^2 - 2 r_2^2) + exp(-2 r_1^2 - 3 r_2^2)
    /// p1: pure, corresponds to f12
    /// p2: dec, corresponds to f23
    /// p3: op_dec, corresponds to f23
    /// p4: pure, corresponds to f23
    /// p5: op_pure, corresponds to f23
    auto [p1,p2,p3,p4,p5]=data.get_ccpairfunctions();
    auto [f1,f2,f3,f4,f5,f]=data.get_functions();

    auto f12=CCConvolutionOperatorPtr<T,LDIM>(world, OT_F12, parameter);
    auto& op=*(f12->get_op());
    std::vector<CCPairFunction<T,NDIM>> vp2({p2});
    std::vector<CCPairFunction<T,NDIM>> vp2ex({p2.swap_particles()});

    // tmp(2) = \int a(1)b(2') f(1,2) d1
    // result=inner(tmp,f1);
    for (auto& p : {p1,p2,p3,p4,p5}) {
        print("working on ",p.name());
        std::vector<CCPairFunction<T,NDIM>> vp({p});
        op.set_particle(1);
        auto op1_p=op(vp);
        double r1=inner(vp2,op1_p);
        double r1ex=inner(vp2ex,op1_p);
        printf("r1   %12.8f\n",r1);
        printf("r1ex %12.8f\n",r1ex);
        op.set_particle(2);
        auto op2_p=op(vp);
        double r2=inner(vp2,op2_p);
        double r2ex=inner(vp2ex,op2_p);
        printf("r2   %12.8f\n",r2);
        printf("r2ex %12.8f\n",r2ex);
    }

    return t1.end();
}

template<typename T, std::size_t NDIM>
int test_scalar_multiplication(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                               const CCParameters& parameter) {
    CCTimer timer(world, "testing");
    test_output t1("CCPairFunction<double,6>::test_scalar_multiplication");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;

//    auto data=get_data<NDIM>(world,parameter);
    auto [f1,f2,f3,f4,f5,f] = data.get_functions();

    std::vector<Function<T,LDIM>> a = {f1, f2};
    std::vector<Function<T,LDIM>> b = {f3, f1};

    print("time in preparation",timer.reset());
    t1.checkpoint(true,"prep");

    CCPairFunction<T,NDIM> p(copy(f));
    CCPairFunction<T,NDIM> p1(a,b);
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

    bool bsuccess=fabs(4.0*norm1-norm2)<FunctionDefaults<LDIM>::get_thresh();
    t1.checkpoint(bsuccess,"scaling");

    return t1.end();
}

template<typename T, std::size_t NDIM>
int test_swap_particles(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                        const CCParameters& parameters) {
    test_output t1("swap_particles<"+std::to_string(NDIM)+">");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;

    // prepare
    auto one = [](const Vector<double,LDIM>& r) { return 1.0; };
    Function<T,LDIM> R2 = FunctionFactory<T,LDIM>(world).f(one);

//    auto data=get_data<NDIM>(world,parameters);
    auto [f1,f2,f3,f4,f5,f] = data.get_functions();
    std::vector<Function<T,LDIM>> a = {f1, f2};
    std::vector<Function<T,LDIM>> b = {f3, f1};

    // test decomposed
    {
        CCPairFunction<T,NDIM> p1(a, b);
        CCPairFunction<T,NDIM> p2(b, a);

        double norm1 = inner(p1, p2.swap_particles(), R2);
        double norm1a = inner(p1, p1, R2);
        // <p1 | p2> = \sum_ij <a_i b_i | a_j b_j> = \sum_ij <a_i|a_j> <b_i|b_j>
        double norm2 = matrix_inner(world, a, a).emul(matrix_inner(world, b, b)).sum();
        print("norm1 ", norm1);
        print("norm1a", norm1a);
        print("norm2 ", norm2);
        t1.checkpoint(std::abs(norm1 - norm2) < FunctionDefaults<LDIM>::get_thresh(), "swap_particles a,b");
    }

    // test pure
    {
        CCPairFunction<T,NDIM> p(copy(f));
        CCPairFunction<T,NDIM> p_swapped=p.swap_particles();
        double pnorm=p.get_function().norm2();
        double psnorm=p_swapped.get_function().norm2();
        print("p/s norm",pnorm,psnorm);

        CCPairFunction<T,NDIM> p1({f1}, {f2});
        CCPairFunction<T,NDIM> p2({f2}, {f1});
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

    return t1.end();
}


template<typename T, std::size_t NDIM>
int test_projector(World& world, std::shared_ptr<NuclearCorrelationFactor> ncf, data<T,NDIM>& data,
                   const CCParameters& parameter) {
    test_output t1("test_projector<"+std::to_string(NDIM)+">");
    static_assert(NDIM%2==0, "NDIM must be even");
    constexpr std::size_t LDIM=NDIM/2;

//    t1.set_cout_to_logger();
//    t1.set_cout_to_terminal();
    const auto [f11,f22,f3,f4,f5,f] = data.get_functions();
    auto f1=copy(f11);  // keep f1, f2 constant for use in other tests!
    auto f2=copy(f22);
    double nf1=f1.norm2();
    f1.scale(1.0/nf1);
    f2.scale(1.0/nf1);
    std::vector<Function<T,LDIM>> a = {f1+f2, f2+f3, f3+f4};
    std::vector<Function<T,LDIM>> b = {f3, f1, f2};
    std::vector<Function<T,LDIM>> o = orthonormalize_cd<T,LDIM>({f1-f3, f5});  // projects on an orthonormal basis
    a = orthonormalize_cd<T,LDIM>(a);  // projects on an orthonormal basis
    auto f12=CCConvolutionOperatorPtr<double,LDIM>(world, OT_F12, parameter);

    {
        auto ovlp=matrix_inner(world,o,o);
        print("<o | o>", ovlp);
    }

    CCPairFunction<T,NDIM> p1(a,b);
    CCPairFunction<T,NDIM> p2(f12,a,b);
    CCPairFunction<T,NDIM> p3(copy(f)); // outer (f1,f2)

    std::vector<CCPairFunction<T,NDIM>> vp1({p1});
    std::vector<CCPairFunction<T,NDIM>> vp2({p2});
    std::vector<CCPairFunction<T,NDIM>> vp3({p3});

    Projector<T,LDIM> O(o,o);
    QProjector<T,LDIM> Q(o,o);
    StrongOrthogonalityProjector<T,LDIM> Q12(world);
    Q12.set_spaces(o);

    double thresh=FunctionDefaults<LDIM>::get_thresh();

    // compute reference values as: (<f1 f2 | projector) | px >
    // compute result values as:    <f1 f2 | (projector | px >)
    Function<T,LDIM> of1=O(f1);
    Function<T,LDIM> of2=O(f2);
    Function<T,LDIM> qf1=Q(f1);
    Function<T,LDIM> qf2=Q(f2);
    std::vector<std::vector<CCPairFunction<T,NDIM>>> vp({vp1,vp2,vp3});
    for (int i=0; i<3; ++i) {
        // O1
        O.set_particle(0);
        {
            double ref=inner({CCPairFunction<T,NDIM>({of1},{f2})},vp[i]);
            auto tmp=O(vp[i]);
            t1.checkpoint(tmp.size()==1,"vector size correct");
            t1.checkpoint(tmp[0].is_decomposed(),"O(argument) is decomposed");
            double result=inner({CCPairFunction<T,NDIM>({f1},{f2})},tmp);
            t1.checkpoint(result,ref,thresh,"O1 p"+std::to_string(i));
        }

        // O2
        O.set_particle(1);
        {
            double ref=inner({CCPairFunction<T,NDIM>({f1},{of2})},vp[i]);
            auto tmp=O(vp[i]);
            t1.checkpoint(tmp.size()==1,"vector size correct");
            t1.checkpoint(tmp[0].is_decomposed(),"O(argument) is decomposed");
            double result=inner({CCPairFunction<T,NDIM>({f1},{f2})},tmp);
            t1.checkpoint(result,ref,thresh,"O2 p"+std::to_string(i));
        }
        // Q1
        Q.set_particle(0);
        {
            double ref=inner({CCPairFunction<T,NDIM>({qf1},{f2})},vp[i]);
            double result=inner({CCPairFunction<T,NDIM>({f1},{f2})},Q(vp[i]));
            t1.checkpoint(result,ref,thresh,"Q1 p"+std::to_string(i));
        }
        // Q2
        Q.set_particle(1);
        {
            double ref=inner({CCPairFunction<T,NDIM>({f1},{qf2})},vp[i]);
            double result=inner({CCPairFunction<T,NDIM>({f1},{f2})},Q(vp[i]));
            t1.checkpoint(result,ref,thresh,"Q2 p"+std::to_string(i));
        }
    }

    // some hardwire test
    {
        // <ab|k ><k |ab>  = <a|k><k|a><b|b>
        auto ak_tensor=matrix_inner(world,a,o);
        auto aa=matrix_inner(world,a,a);
        auto bb=matrix_inner(world,b,b);
        Projector<double,LDIM> O1(a,a);
        O.set_particle(0);
        O1.set_particle(0);
        Projector<double,LDIM> O2(a,a);
        O2.set_particle(1);
        double n1=inner(vp1,O1(vp2));
        double n1a=inner(O1(vp1),vp2);
        t1.checkpoint(fabs(n1-n1a)<thresh,"O1");
        print("n1,n1a",n1,n1a);
        double n2=inner(vp1,O2(vp2));
        double n2a=inner(O2(vp1),vp2);
        t1.checkpoint(fabs(n2-n2a)<thresh,"O2");
        print("n2,n2a",n2,n2a);

    }


    // test {O,Q} with particles 1,2 on all three CCPairFunction<double,6> types
    {
        int i=0;
        for (auto p : {p1,p2,p3}) {
            std::string s="CCPairFunction<double,"+std::to_string(NDIM)+"> p"+std::to_string(++i);
            std::vector<CCPairFunction<T,NDIM>> vp({p});
            for (int particle : {0,1}) {
                O.set_particle(particle);
                Q.set_particle(particle);
                auto Op = O(vp);
                auto Qp = Q(vp);

                double n1 = inner(Op, vp3) + inner(Qp, vp3);
                double n2 = inner(vp, vp3);
                print("n1 (RI), n2", n1, n2, fabs(n1 - n2));
                t1.checkpoint(fabs(n1 - n2) < FunctionDefaults<LDIM>::get_thresh(), "RI with particle "+std::to_string(particle)+" on "+s);

                auto Op1 = O(Op);
                double n3=inner(Op, vp3);
                double n4=inner(Op1, vp3);
                print("n3, n4", n3, n4);
                t1.checkpoint(fabs(n3-n4) < FunctionDefaults<LDIM>::get_thresh(), "idempotency with particle "+std::to_string(particle)+" on "+s);

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
            t1.checkpoint(zero, FunctionDefaults<3>::get_thresh(), "SO operator on "+s);
        }
    }
    return t1.end();
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
    int k=5;
    double thresh=1.e-4;
    if (parser.key_exists("k")) k=std::stoi(parser.value("k"));
    if (parser.key_exists("thresh")) thresh=std::stod(parser.value("thresh"));
    FunctionDefaults<1>::set_k(k);
    FunctionDefaults<2>::set_k(k);
    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<4>::set_k(k);
    FunctionDefaults<5>::set_k(k);
    FunctionDefaults<6>::set_k(k);
    FunctionDefaults<1>::set_thresh(thresh);
    FunctionDefaults<2>::set_thresh(thresh);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<4>::set_thresh(thresh);
    FunctionDefaults<5>::set_thresh(thresh);
    FunctionDefaults<6>::set_thresh(thresh);
    FunctionDefaults<6>::set_tensor_type(TT_2D);
    FunctionDefaults<1>::set_cubic_cell(-10.,10.);
    FunctionDefaults<2>::set_cubic_cell(-10.,10.);
    FunctionDefaults<3>::set_cubic_cell(-10.,10.);
    FunctionDefaults<4>::set_cubic_cell(-10.,10.);
    FunctionDefaults<5>::set_cubic_cell(-10.,10.);
    FunctionDefaults<6>::set_cubic_cell(-10.,10.);
    print("numerical parameters: k, eps(3D), eps(6D)", FunctionDefaults<3>::get_k(), FunctionDefaults<3>::get_thresh(),
          FunctionDefaults<6>::get_thresh());
    int isuccess=0;
#ifdef USE_GENTENSOR

    CCParameters ccparam(world, parser);
    try {
        parser.set_keyval("geometry", "he");
        parser.print_map();
        Molecule mol(world, parser);
        mol.print();
        std::shared_ptr<NuclearCorrelationFactor> ncf = create_nuclear_correlation_factor(world,
                                             mol, nullptr, std::make_pair("slater", 2.0));

        auto data2=data<double,2>(world,ccparam);
        auto data4=data<double,4>(world,ccparam);
        auto data6=data<double,6>(world,ccparam);

        isuccess+=test_constructor<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_load_store<double,2>(world,ncf,data2,ccparam);
        isuccess+=test_operator_apply<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_norm<double,2>(world,ncf,data2,ccparam);
        isuccess+=test_transformations<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_multiply_with_f12<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_inner<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_multiply<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_swap_particles<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_scalar_multiplication<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_projector<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_partial_inner_3d<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_partial_inner_6d<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_apply<double,2>(world, ncf, data2, ccparam);
        isuccess+=test_consolidate<double,2>(world, ncf, data2, ccparam);


//        isuccess+=test_constructor<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_load_store<double,4>(world,ncf,data2,ccparam);
//        isuccess+=test_operator_apply<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_transformations<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_multiply_with_f12<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_inner<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_multiply<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_swap_particles<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_scalar_multiplication<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_projector<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_partial_inner_3d<double,4>(world, ncf, data4, ccparam);
//        isuccess+=test_partial_inner_6d(world, ncf, data4, ccparam);
//        isuccess+=test_apply(world, ncf, data, ccparam);


//        isuccess+=test_projector<double,6>(world, ncf, data6, ccparam);

        data2.clear();
        data4.clear();
        data6.clear();
        world.gop.fence();
    } catch (std::exception& e) {
        madness::print("an error occurred");
        madness::print(e.what());
        world.gop.fence();
    }
    world.gop.fence();
#else
    print("could not run test_ccpairfunction: U need to compile with ENABLE_GENTENSOR=1");
#endif
    finalize();

    return isuccess;
}