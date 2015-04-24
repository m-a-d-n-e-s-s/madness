#include <madness/mra/mra.h>
#include <madness/mra/operator.h>
using namespace madness;

double gaussian(const coord_3d& r) {
    double x=r[0], y=r[1], z=r[2];
    return exp(-(x*x + y*y + z*z));
}

int main(int argc, char** argv) {
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);
    startup(world,argc,argv);

    FunctionDefaults<3>::set_cubic_cell(-6.0,6.0);
    real_convolution_3d op = CoulombOperator(world, 1e-4, 1e-6);
    real_function_3d g = real_factory_3d(world).f(gaussian);

    print(g.trace()); // exact trace is Pi^(3/2)=5.56832799683175
    print(g.norm2()); // exact norm2 is (Pi/2)^3/4=1.40310414553422
    print(g.inner(op(g))); // exact answer is Pi^(5/2)*sqrt(2) = 24.7394294511936

    finalize();
    return 0;
}
