#include <madness/mra/mra.h>
#include <iostream>
#include <cmath>

using namespace madness;

// Define a test mathematical function: f(x) = sin(x)
double my_func(const coord_1d& r) {
    return std::sin(r[0]);
}

int main(int argc, char** argv) {
    // 1. Initialize the MADNESS parallel runtime. initialize() constructs the
    //    default World on MPI_COMM_WORLD and returns it -- do NOT construct
    //    another World on the same communicator, and do not let any World or
    //    Function outlive finalize() (see the "Runtime" invariants in CLAUDE.md).
    World& world = initialize(argc, argv);

    // 2. Initialize the numerical environment (MRA basis, FunctionDefaults, etc.)
    startup(world, argc, argv);

    if (world.rank() == 0) {
        print("=================================================");
        print("   Hello World -- MADNESS MRA Application        ");
        print("=================================================");
        print("World rank:", world.rank(), "of size:", world.size());
    }

    // 3. Define the 1D computation box: interval [0, pi]. FunctionDefaults must
    //    be set before any Function of that dimension is constructed.
    FunctionDefaults<1>::set_cubic_cell(0.0, M_PI);

    // 4. Scope every Function so it is destroyed before finalize() tears the
    //    runtime down.
    {
        // Project the analytic function into an adaptive multiresolution
        // wavelet representation
        real_function_1d f = real_factory_1d(world).f(my_func);

        // Compute the definite integral over the cell:
        // integral_0^pi sin(x) dx = 2.0. trace() is an observing operation and
        // therefore fences, so no explicit world.gop.fence() is needed here.
        const double integral = f.trace();

        if (world.rank() == 0) {
            print("Computed integral of sin(x) on [0, pi]:", integral);
            print("Expected exact mathematical value     :", 2.0);
            print("Difference from exact                 :", std::abs(integral - 2.0));
            print("Calculation completed successfully.");
        }
    }

    // 5. Clean up and finalize the runtime
    finalize();
    return 0;
}
