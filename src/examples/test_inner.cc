

#define WORLD_INSTANTIATE_STATIC_TEMPLATES
#include <madness/mra/mra.h>

using namespace madness;

static const double R = 1.4;    // bond length
static const double L = 5.0*R; // box size
static const long k = 8;        // wavelet order
static const double thresh = 1e-6; // precision

static double easyfunc1(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (sin(x*x + y*y + z*z));
}

static double easyfunc2(const coord_3d& r) {
    const double x=r[0], y=r[1], z=r[2];
    return (exp(- x*x - y*y - z*z));
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

    real_function_3d beta = real_factory_3d(world).f(easyfunc1);
    real_function_3d zeta = real_factory_3d(world).f(easyfunc2);

    double norm = zeta.norm2();
    zeta.compress();
    double normc = zeta.norm2();
    zeta.reconstruct();
    double normr = zeta.norm2();

    printf("%7.10f\n", zeta.inner_ext(easyfunc1));
    printf("%7.10f\n", zeta.inner_ext(easyfunc2));
    printf("%7.10f\n", beta.inner_ext(easyfunc1));
    printf("%7.10f\n", beta.inner_ext(easyfunc2));

    if (world.rank() == 0) {
    print("                    Norm is ", norm, normc, normr);
    }

    world.gop.fence();

    finalize();
    return 0;
}
