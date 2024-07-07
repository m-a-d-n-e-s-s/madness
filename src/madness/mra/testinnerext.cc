

#include <madness/mra/mra.h>
#include <array>
using namespace madness;

bool smalltest = true;

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %-20.20s %8.2fs %8.2fs\n", msg, sss, ttt)

typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

static const double R = 1.4;    // bond length
static const double L = 5.0*R; // box size
static const long k = 6;        // wavelet order
static const double thresh = 1e-6; // precision

static double alpha_func(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return ((x*x + y*y + z*z) * sin(x*x + y*y + z*z));
};

static double beta_func(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (exp(- x*x - y*y - z*z));
};

class alpha_functor : public FunctionFunctorInterface<double,3> {
private:
    double coeff;
public:
    alpha_functor(double coeff=1.0) : coeff(coeff) {}

    virtual double operator()(const coord_3d& r) const {
        const double x=r[0], y=r[1], z=r[2];
        return (coeff * (x*x + y*y + z*z) * sin(x*x + y*y + z*z));
    }
};

class beta_functor : public FunctionFunctorInterface<double,3> {
private:
    double coeff;
public:
    beta_functor(double coeff=1.0) : coeff(coeff) {}

    virtual double operator()(const coord_3d& r) const {
        return (coeff * exp(- r[0]*r[0] - r[1]*r[1] - r[2]*r[2]));
    }
};

struct gauss_1d {
    gauss_1d(double s) : s(s) {}
    double s=1.0;
    double operator()(const double x) const {return 1.0/(s * std::pow(constants::pi,0.25))* exp(-0.5*x*x/(s*s));}
};
template<typename T, std::size_t NDIM>
struct gauss {
    std::array<double,NDIM> e;
    gauss(const std::array<double,NDIM> exponent) : e(exponent) {};

    T operator()(const Vector<double,NDIM>& coord) const {
        double result=1.0;
        for (std::size_t i=0; i<NDIM; ++i) result*=gauss_1d(e[i])(coord[i]);
        return result;
    }
};

bool is_like(double a, double b, double tol) {
    return (std::abs((a - b)/a) <= tol);
}
bool test(std::string msg, double a, double b, double tol=thresh) {
    int len=msg.size();
    std::string padding(40-len,' ');
    bool ok=(std::abs((a - b)/a) <= tol);
    print(msg+padding,a,b,std::abs((a-b)/a),ok);
    return ok;
}
bool test_loose1(std::string msg, double a, double b, double tol=thresh) {
    return test(msg,a,b,tol*10.0);
}


int test_tight_diffuse(World& world) {
    FunctionDefaults<4>::set_thresh(1.e-5);
    double a=1.e2;
    double b=1.e-2;
    for (int i=0; i<4; ++i) {
        a=std::pow(10.0,double(i));
        b=std::pow(0.1,double(i));
        print("a,b",a,b);

        real_function_2d aa=real_factory_2d(world).functor([&](const coord_2d& r){return exp(-a*inner(r,r));});
        real_function_2d bb=real_factory_2d(world).functor([&](const coord_2d& r){return exp(-b*inner(r,r));});
        real_function_4d cc=real_factory_4d(world).functor([&](const coord_4d& r){return exp(-b*inner(r,r));});
        real_function_4d dd=CompositeFactory<double,4,2>(world).particle1(aa).particle2(copy(aa));
        aa.print_size("exp(-1000 r^2");
        bb.print_size("exp(-0.001 r^2");
        double result=inner(cc,dd);
        double refresult=std::pow(constants::pi/(a+b),2.0);
        print("result,refresult,error",result,refresult,result-refresult);
        MADNESS_CHECK(test(" inner(exp(-a r^2 , exp(-b r^2))  ", result,refresult));
    }




    return 0;


}

int test_partial_inner(World& world) {
    bool do_low_rank=false;
#if HAVE_GENTENSOR
    do_low_rank=true;
#endif
    print("\ntesting partial inner; low rank: ",do_low_rank,"\n");

    real_function_1d one_1d=real_factory_1d(world).functor([](const coord_1d& r){return 1.0;});
    real_function_2d one_2d=real_factory_2d(world).functor([](const coord_2d& r){return 1.0;});
    real_function_1d g1=real_factory_1d(world).functor(gauss<double,1>({1.0}));
    real_function_1d g2=real_factory_1d(world).functor(gauss<double,1>({2.0}));
    real_function_1d g3=real_factory_1d(world).functor(gauss<double,1>({3.0}));
//    real_function_1d g3=real_factory_1d(world).functor([](const Vector<double,1>& r) {return exp(r[0])*cos(r[0]);} );
    real_function_1d g4=real_factory_1d(world).functor(gauss<double,1>({4.0}));

    real_function_2d f2=real_factory_2d(world).functor(gauss<double,2>({1.0,2.0}));
    if (do_low_rank) FunctionDefaults<2>::set_tensor_type(TT_2D);
    real_function_2d f2_svd=real_factory_2d(world).functor(gauss<double,2>({1.0,2.0}));
    FunctionDefaults<2>::set_tensor_type(TT_FULL);
    real_function_2d f2_swap=real_factory_2d(world).functor(gauss<double,2>({2.0,1.0}));
    real_function_2d f2_tight=real_factory_2d(world).functor(gauss<double,2>({3.0,4.0}));
    real_function_3d f3=real_factory_3d(world).functor(gauss<double,3>({1.0,2.0,3.0}));
//    real_function_3d f3=real_factory_3d(world).functor([](const Vector<double,3>& r) {return gauss<double,2>({1.0,2.0})({r[0],r[1]}) * exp(r[2])*cos(r[2]);});

    double g11=inner(g1,g1);
    double g12=inner(g1,g2);
    double g13=inner(g1,g3);
    double g14=inner(g1,g4);
    double g22=inner(g2,g2);
    double g23=inner(g2,g3);
    double g24=inner(g2,g4);
    double g33=inner(g3,g3);
    // double g34=inner(g3,g4);
    // double g44=inner(g4,g4);

    {   // test unevenly refined functions
        real_function_2d f12=real_factory_2d(world)
                .functor([](const coord_2d& r) {
                    return exp(-abs(r[0]-r[1]));
                });

        print("done with projection");
        real_function_2d r = inner(f2, f12, {0}, {0});
    }
    {
        real_function_2d r = inner(f2, f2, {0}, {1});
        double n=inner(f2,r);
//        MADNESS_CHECK(test(" int f2(1,2)*f2(2,1) d1 (full)", n,g12*g12*g12));
        test(" int f2(1,2)*f2(2,1) d1 (full)", n,g12*g12*g12);


        FunctionDefaults<2>::set_tensor_type(TT_2D);
        real_function_2d r_svd = inner(f2_svd, f2_svd, {0}, {1});
        FunctionDefaults<2>::set_tensor_type(TT_FULL);
        // double n_svd=inner(f2_svd,r_svd);
        MADNESS_CHECK(test(" int f2(1,2)*f2(2,1) d1 (svd)", n,g12*g12*g12));
    }
    {
        real_function_1d r=inner(f2_svd,g1,{0},{0});
        double n=inner(g1,r);
        MADNESS_CHECK(test(" int f2(1,2)*g1(1) d1  ", n,g11*g12));
    }
    if (do_low_rank) {
        FunctionDefaults<6>::set_thresh(1.e-4);
        FunctionDefaults<6>::set_tensor_type(TT_2D);
//        real_function_6d f6=real_factory_6d(world).functor(gauss<double,6>({1.0,2.0,3.0,1.0,2.0,3.0}));
        real_function_6d f6=hartree_product(f3,f3);
        double cpu0=cpu_time();
        double p1 = inner(f6.project_out(f3, 1), f3);
        double cpu1=cpu_time();
        double p2 = inner(inner(f6,f3, {0,1,2},{0,1,2}), f3);
        double cpu2=cpu_time();
        MADNESS_CHECK(test_loose1("project_out 1 ", p1, g11*g22*g33*g11*g22*g33));
        MADNESS_CHECK(test_loose1("project_out 2 ", p2, g11*g22*g33*g11*g22*g33));
        print("timings project_out, partial_inner",cpu1-cpu0, cpu2-cpu1);
    }
    {
        real_function_1d r=inner(f2_svd,g1,{1},{0});
        double n=inner(g3,r);
        MADNESS_CHECK(test(" int f2(1,2)*g1(2) d2  ", n,g12*g13));
    }
    {
        real_function_1d r=inner(g1,f2_svd,{0},{1});
        double n=inner(g3,r);
        MADNESS_CHECK(test(" int g1(2)*f2(1,2) d2  ", n,g12*g13));
    }
    {
        real_function_1d r=inner(g1,f2_svd,{0},{0});
        double n=inner(g3,r);
        MADNESS_CHECK(test(" int g1(1)*f2(1,2) d1 (svd) ", n,g11*g23));
    }
    {
        real_function_1d r=inner(g1,f2,{0},{0});
        double n=inner(g3,r);
        MADNESS_CHECK(test(" int g1(1)*f2(1,2) d1 (full) ", n,g11*g23));
    }

    {

        double c1=inner(inner(f2,g2,{0},{0}),g3);
        MADNESS_CHECK(test(" int f(1,2)*g2(1) d1 ", c1, g12 * g23));

        double c2=inner(inner(f2,g2,{1},{0}),g3);
        MADNESS_CHECK(test("result 1 - 2", c2, g22 * g13));

        double c3=inner(inner(g1,f2,{0},{0}),g3);
        MADNESS_CHECK(test("result 1 - 3", c3, g11 * g23));

        double c4=inner(inner(g1,f2,{0},{1}),g3);
        MADNESS_CHECK(test("result 1 - 4", c4, g12 * g13));
    }
    {
        double c1 = inner(inner(f2, f2, {0}, {0}), f2_tight);
        MADNESS_CHECK(test(" f2 f2 0 0", c1, g11 * g23 * g24));

        double c2 = inner(inner(f2, f2, {0}, {1}), f2_tight);
        MADNESS_CHECK(test(" f2 f2 0 1", c2, g12 * g23 * g14));

        double c3 = inner(inner(f2, f2, {1}, {0}), f2_tight);
        MADNESS_CHECK(test(" f2 f2 1 0 ", c3, g12 * g13 * g24));

        double c4 = inner(inner(f2, f2, {1}, {1}), f2_tight);
        MADNESS_CHECK(test(" f2 f2 1 1", c4, g22 * g13 * g14));
    }

    {
        double c2 = inner(inner(f3, f2, {0}, {0}), f3);
        MADNESS_CHECK(test(" f3 f2 0 0", c2, g11 * g12 * g23 * g23));

        double c3 = inner(inner(f3, f2, {0}, {1}), f3);
        MADNESS_CHECK(test(" f3 f2 0 1", c3, g12 * g12 * g13 * g23));

        double c4 = inner(inner(f3, f2, {2}, {0}), f3);
        MADNESS_CHECK(test(" f3 f2 2 0", c4, g11 * g22 * g13 * g23));

        double c5 = inner(inner(f3, f2, {2}, {1}), f3);
        MADNESS_CHECK(test(" f3 f2 2 1", c5, g11 * g22 * g13 * g23));

        double c6 = inner(inner(f2, f3, {1}, {2}), f3);
        MADNESS_CHECK(test(" f2 f3 1 2", c6, g11 * g12 * g23 * g23));

    }
    {
        double c1=inner(inner(f2,f3,{0,1},{0,1}),g3);
        MADNESS_CHECK(test(" f2 f3 {01} {01}", c1, g11 * g22 * g33));

        double c2=inner(inner(f3,f2,{0,1},{0,1}),g3);
        MADNESS_CHECK(test(" f3 f2 {01} {01}", c2, g11 * g22 * g33));

        double c3=inner(inner(f3,f2,{1,2},{0,1}),g3);
        MADNESS_CHECK(test(" f3 f1 {12} {01}", c3, g13 * g12 * g23));

    }

    return 0;
}

template<std::size_t NDIM>
void initialize(World& world) {
    FunctionDefaults<NDIM>::set_defaults(world);
    FunctionDefaults<NDIM>::set_k(k);
    FunctionDefaults<NDIM>::set_thresh(thresh);
    FunctionDefaults<NDIM>::set_refine(true);
    FunctionDefaults<NDIM>::set_initial_level(5);
    FunctionDefaults<NDIM>::set_truncate_mode(1);
    FunctionDefaults<NDIM>::set_cubic_cell(-L/2, L/2);
}
int main(int argc, char** argv) {
    World& world=initialize(argc, argv,false);
//    World world(SafeMPI::COMM_WORLD);

    int success = 0;

    startup(world,argc,argv,true);
    std::cout.precision(6);

    if (getenv("MAD_SMALL_TESTS")) smalltest=true;
    for (int iarg=1; iarg<argc; iarg++) if (strcmp(argv[iarg],"--small")==0) smalltest=true;
    std::cout << "small test : " << smalltest << std::endl;

    initialize<1>(world);
    initialize<2>(world);
    initialize<3>(world);
    initialize<4>(world);
    initialize<5>(world);
    initialize<6>(world);

    test_partial_inner(world);
    test_tight_diffuse(world);

    if (!smalltest) {
        real_function_3d alpha1 = real_factory_3d(world).f(alpha_func);
        real_functor_3d alpha1_ffi = real_functor_3d(new alpha_functor());
        
        if (world.rank() == 0) {
            print("***************************************************************************");
            print("alpha is a highly oscillatory function (x^2 sin(x^2)).");
            print("beta is a pretty boring Gaussian.");
            print("For the first two timers, the cost of computing alpha in");
            print("the numerical basis is not included. We see that since beta");
            print("is a simple Gaussian, it is cheaper to use the inner() method.\n");
        }
        
        START_TIMER;
        real_function_3d beta = real_factory_3d(world).f(beta_func);
        double ab = alpha1.inner(beta);
        END_TIMER("1. < a | b >)");
        
        START_TIMER;
        real_functor_3d beta_ffi = real_functor_3d(new beta_functor());
        double ab_ffi = alpha1.inner_ext(beta_ffi);
        END_TIMER("3. < a | b_ffi >");
        
        if (world.rank() == 0) {
            print("\n***************************************************************************");
            print("For the next two timers, the cost of computing beta in");
            print("the numerical basis is not included. We see that since alpha");
            print("is complicated, it is cheaper to use the inner_ext() method.\n");
        }
        
        START_TIMER;
        real_function_3d alpha = real_factory_3d(world).f(alpha_func);
        double ba = beta.inner(alpha);
        END_TIMER("4. < b | a >");
        
        START_TIMER;
        real_functor_3d alpha_ffi = real_functor_3d(new alpha_functor());
        double ba_ffi = beta.inner_ext(alpha_ffi);
        END_TIMER("6. < b | a_ffi >");
        
        double aa = alpha.inner(alpha);
        double bb = beta.inner(beta);
        
        double aa_ffi = alpha.inner_ext(alpha_ffi);
        double bb_ffi = beta.inner_ext(beta_ffi);
        
        if (world.rank() == 0) {
            print("\nTest for correctness");
            print("***************************************************************************");
            printf("<a|a> (using inner() with Function) =                      %7.10f\n", aa);
            printf("<a|a> (using inner_ext() with FunctionFunctor Interface) = %7.10f\n", aa_ffi);
            print("***************************************************************************");
            printf("<b|b> (using inner() with Function) =                      %7.10f\n", bb);
            printf("<b|b> (using inner_ext() with FunctionFunctor Interface) = %7.10f\n", bb_ffi);
            print("***************************************************************************");
            printf("<a|b> (using inner() with Function) =                      %7.10f\n", ab);
            printf("<a|b> (using inner_ext() with FunctionFunctor Interface) = %7.10f\n", ab_ffi);
            printf("<b|a> (using inner() with Function) =                      %7.10f\n", ba);
            printf("<b|a> (using inner_ext() with FunctionFunctor Interface) = %7.10f\n", ba_ffi);
            print("***************************************************************************");
        }
        
        real_function_3d alphabeta = real_factory_3d(world);
        alphabeta = alpha + beta;
        double aba = alphabeta.inner(alpha);
        double aba_ffi = alphabeta.inner_ext(alpha_ffi);
        
        if (world.rank() == 0) {
            print("\nCheck that inner_ext works for Function that lacks a functor");
            print("***************************************************************************");
            printf("<a+b|a> (using inner() with Function) =                      %7.10f\n", aba);
            printf("<a+b|a> (using inner_ext() with FunctionFunctor Interface) = %7.10f\n", aba_ffi);
            print("***************************************************************************");
        }
        
        if (not is_like(aa, aa_ffi, thresh)) ++success;
        if (not is_like(bb, bb_ffi, thresh)) ++success;
        if (not is_like(ab, ab_ffi, thresh)) ++success;
        if (not is_like(ab, ba_ffi, thresh)) ++success;
    }

    world.gop.fence();

    finalize();
    return success;
}
