

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>

using namespace madness;

double ttt, sss;
#define START_TIMER world.gop.fence(); ttt=wall_time(); sss=cpu_time()
#define END_TIMER(msg) ttt=wall_time()-ttt; sss=cpu_time()-sss; if (world.rank()==0) printf("timer: %-20.20s %8.2fs %8.2fs\n", msg, sss, ttt)

static const double R = 1.4;    // bond length
static const double L = 5.0*R;  // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision

/* static double alpha_func(const coord_3d& r) { */
/*     const double x=r[0], y=r[1], z=r[2]; */
/*     return ((x*x + y*y + z*z) * sin(x*x + y*y + z*z)); */
/* }; */

static double beta_func(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double x0=-1.0, y0=-1.0, z0=-1.0;
    const double a=1.0;
    return (a * exp(- (x-x0)*(x-x0) - (y-y0)*(y-y0) - (z-z0)*(z-z0)));
};

static double gamma_func(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double x0=2.0, y0=2.0, z0=2.0;
    const double a=0.5, sig2=0.0625;
    return (a * exp(- sig2*(x-x0)*(x-x0) - sig2*(y-y0)*(y-y0) - sig2*(z-z0)*(z-z0)));
};

/* static double alpha_p_beta_func(const coord_3d& r) { */
/*     const double x=r[0], y=r[1], z=r[2]; */
/*     const double x0=-1.0, y0=-1.0, z0=-1.0; */
/*     const double a=1.0; */
/*     return ((x*x + y*y + z*z) * sin(x*x + y*y + z*z) + a * exp(- (x-x0)*(x-x0) - (y-y0)*(y-y0) - (z-z0)*(z-z0))); */
/* } */

static double beta_p_gamma_func(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    const double x0=-1.0, y0=-1.0, z0=-1.0;
    const double x1=2.0, y1=2.0, z1=2.0;
    const double a=1.0, b=0.5, sig2=0.0625;
    return (a * exp(- (x-x0)*(x-x0) - (y-y0)*(y-y0) - (z-z0)*(z-z0))
          + b * exp(- sig2*(x-x1)*(x-x1) - sig2*(y-y1)*(y-y1) - sig2*(z-z1)*(z-z1)));
}

static double onefn(const coord_3d& r) {
    return 1.0;
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    startup(world,argc,argv);
    std::cout.precision(6);

    FunctionDefaults<3>::set_defaults(world);

    FunctionDefaults<3>::set_k(k);
    FunctionDefaults<3>::set_thresh(thresh);
    FunctionDefaults<3>::set_refine(true);
    FunctionDefaults<3>::set_initial_level(5);
    FunctionDefaults<3>::set_truncate_mode(1);
    FunctionDefaults<3>::set_cubic_cell(-L/2, L/2);

    // Create numerical functions for alpha, beta, gamma
    // and also for alpha + beta and beta + gamma
    print("***************************************************");
    print("Creating madness functions");
    /* real_function_3d alpha = real_factory_3d(world).f(alpha_func); */
    real_function_3d beta = real_factory_3d(world).f(beta_func);
    real_function_3d gamma = real_factory_3d(world).f(gamma_func);
    /* real_function_3d apb = real_factory_3d(world).f(alpha_p_beta_func); */
    real_function_3d bpc = real_factory_3d(world).f(beta_p_gamma_func);

    // Create empty numerical functions to store
    // the results of gaxpy_ext()
    print("Creating empty functions to store results.");
    /* real_function_3d apb_ext = real_factory_3d(world); */
    real_function_3d bpc_ext = real_factory_3d(world);

    // Create empty numerical functions to store
    // the results of gaxpy computed with the 
    // overloaded operators + and *
    /* real_function_3d apb_gax = real_factory_3d(world); */
    real_function_3d bpc_gax = real_factory_3d(world);

    print("Adding madness functions using overloaded + and * operators");
    bpc_gax = beta + gamma;

    /* bpc_gax.compress(); */
    /* bpc_gax.reconstruct(); */
    /* Translation l1,l2,l3; */
    /* l1 = 3; */
    /* l2 = 3; */
    /* l3 = 2; */
    /* const Vector<Translation,3> ll=vec(l1,l2,l3); */
    /* Key<3> test_key(4, ll); */ 
    /* print(test_key); */
    /* print(bpc.get_impl().get()->project(test_key)); */

    print("Adding using gaxpy_ext_local");
    bpc_ext.gaxpy_ext_local<double>(beta, gamma_func, 1.0, 1.0, 1e-02);

    bpc_gax.print_tree();
    bpc_ext.print_tree();

    /* if (world.rank() == 0) {} */

    world.gop.fence();

    finalize();
    return 0;
}
