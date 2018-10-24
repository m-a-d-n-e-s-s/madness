#include <madness/mra/mra.h>
//For basic computing with MADNESS functions, the only header file i need
// to include is <mra/mra.h>

using namespace madness;

double myf(const coor_1_d &r)
{
    return std::sin(r[0]);
} //translate the numerical representation from math

int main(int argc, char **argv)
{
    initialize(argc, argv);           //initialize the madness runtime in parallel
    World world(SafeMPI::COMM_WORLD); //initialzing the MADNESS numerical enivronment

    startup(world, argc, argv);

    FunctionDefaults<1>::set_cubic_cell(0, 10);

    real_function_1d f = real_factory_1d(world).f(myf);

    double integral = f.trace();
    if (world.rank() == 0)
        print("The result is ", integral);

    finalize();
    return 0;
}
