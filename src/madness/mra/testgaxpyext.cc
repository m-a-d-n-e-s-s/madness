
#include <madness/mra/mra.h>

using namespace madness;

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %-20.20s %8.2fs %8.2fs\n", msg, sss, ttt)

typedef std::shared_ptr< FunctionFunctorInterface<double,3> > functorT;

static const double R = 1.4;    // bond length
static const double L = 5.0*R;  // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision

static double alpha_func(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return ((x*x + y*y + z*z) * sin(x*x + y*y + z*z));
};

static double beta_func(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double x0=-1.0, y0=-1.0, z0=-1.0;
    const double a=1.0;
    return (a * exp(- (x-x0)*(x-x0) - (y-y0)*(y-y0) - (z-z0)*(z-z0)));
};

/* static double gamma_func(const coord_3d& r) { */
/*    const double x=r[0], y=r[1], z=r[2]; */
/*    const double x0=2.0, y0=2.0, z0=2.0; */
/*    const double a=0.5, sig2=0.0625; */
/*    return (a * exp(- sig2*(x-x0)*(x-x0) - sig2*(y-y0)*(y-y0) - sig2*(z-z0)*(z-z0))); */
/* }; */

/* static double alpha_p_beta_func(const coord_3d& r) { */
/*     const double x=r[0], y=r[1], z=r[2]; */
/*     const double x0=-1.0, y0=-1.0, z0=-1.0; */
/*     const double a=1.0; */
/*     return ((x*x + y*y + z*z) * sin(x*x + y*y + z*z) + a * exp(- (x-x0)*(x-x0) - (y-y0)*(y-y0) - (z-z0)*(z-z0))); */
/* } */

/* static double beta_p_gamma_func(const coord_3d& r) { */
/*     const double x=r[0], y=r[1], z=r[2]; */
/*     const double x0=-1.0, y0=-1.0, z0=-1.0; */
/*     const double x1=2.0, y1=2.0, z1=2.0; */
/*     const double a=1.0, b=0.5, sig2=0.0625; */
/*     return (a * exp(- (x-x0)*(x-x0) - (y-y0)*(y-y0) - (z-z0)*(z-z0)) */
/*           + b * exp(- sig2*(x-x1)*(x-x1) - sig2*(y-y1)*(y-y1) - sig2*(z-z1)*(z-z1))); */
/* } */

bool is_like(double a, double b, double tol) {
    return (std::abs((a - b)/a) <= tol);
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    int success = 0;

    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_defaults(world);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_truncate_on_project(true);
    FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);

    real_function_3d alpha1 = real_factory_3d(world).f(alpha_func);
    /* real_function_3d bpc = real_factory_3d(world).f(beta_p_gamma_func); */

    if (world.rank() == 0) {
        print("***************************************************************************");
        print("alpha is a highly oscillatory function (x^2 sin(x^2)).");
        print("beta is a pretty boring Gaussian.");
        print("For the first two timers, the cost of computing alpha in");
        print("the numerical basis is not included.\n");
    }

    START_TIMER;
    real_function_3d beta = real_factory_3d(world).f(beta_func);
    real_function_3d apb1 = real_factory_3d(world);
    apb1 = alpha1 + beta;
    END_TIMER("1. alpha + beta");

    START_TIMER;
    real_function_3d apb_ext1 = real_factory_3d(world);
    apb_ext1.gaxpy_ext<double>(alpha1, beta_func, 1.0, 1.0, 1e-06, true);
    END_TIMER("2. alpha + beta_func");

    if (world.rank() == 0) {
        print("***************************************************************************");
        print("alpha is a highly oscillatory function (x^2 sin(x^2)).");
        print("beta is a pretty boring Gaussian.");
        print("For the next two timers, the cost of computing beta in");
        print("the numerical basis is not included.\n");
    }
    
    START_TIMER;
    real_function_3d alpha2 = real_factory_3d(world).f(alpha_func);
    real_function_3d apb2 = real_factory_3d(world);
    apb2 = alpha2 + beta;
    END_TIMER("3. beta + alpha");

    START_TIMER;
    real_function_3d apb_ext2 = real_factory_3d(world);
    apb_ext2.gaxpy_ext<double>(beta, alpha_func, 1.0, 1.0, 1e-06, true);
    END_TIMER("4. beta + alpha_func");

    double norm_apb1 = apb1.norm2();
    double norm_apb2 = apb2.norm2();
    double norm_apbext1 = apb_ext1.norm2();
    double norm_apbext2 = apb_ext2.norm2();

    if (world.rank() == 0) {
        print("***************************************************************************");
        print("Test for correctness by comparing the 2-norms.");
        printf("(1) norm of apb =     %20.16f\n", norm_apb1);
        printf("(2) norm of apb_ext = %20.16f\n", norm_apbext1);
        printf("(3) norm of apb =     %20.16f\n", norm_apb2);
        printf("(4) norm of apb_ext = %20.16f\n", norm_apbext2);
    }

    /* START_TIMER; */
    /* real_function_3d gamma = real_factory_3d(world).f(gamma_func); */
    /* real_function_3d bpc = real_factory_3d(world); */
    /* bpc = alpha1 + beta; */
    /* END_TIMER("1. bpc = beta + gamma"); */

    /* START_TIMER; */
    /* real_function_3d apb_ext1 = real_factory_3d(world); */
    /* apb_ext1.gaxpy_ext<double>(alpha, beta_func, 1.0, 1.0, 1e-06, true); */
    /* END_TIMER("2. apb_ext = alpha + beta_func"); */

    /* START_TIMER; */
    /* real_function_3d apb_ext2 = real_factory_3d(world); */
    /* apb_ext2.gaxpy_ext<double>(beta, alpha_func, 1.0, 1.0, 1e-06, true); */
    /* END_TIMER("3. apb_ext = alpha_func + beta"); */

    if (not is_like(norm_apb1, norm_apbext1, thresh)) ++success;
    if (not is_like(norm_apb2, norm_apbext2, thresh)) ++success;

    world.gop.fence();

    finalize();
    return success;
}
