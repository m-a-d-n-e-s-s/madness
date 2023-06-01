# Numerical library

MADNESS provides a high-level environment for the solution of integral and differential equations 
in many dimensions using adaptive, fast methods with guaranteed precision based on multi-resolution 
analysis and novel separated representations. 

Useful examples can be found in the `src/examples` directory, where the use of the numerical library can be 
practiced by example.
 * `sininteg.cc`: Create a function and integrate
 * `hatom_energy.cc`: Compute the energy of the hydrogen atom 
 * `heat.cc`: apply the Greens function to the heat equation
 * `3dharmonic.cc`: solve the 3D harmonic oscillator
 * `h2.cc`: solve the Hartree-Fock equations for the H2 molecule
 * `nonlinschro.cc`: use the KAIN solver for accelerating the solution of a system of non-linear equations
 * `hedft.cc`: use density functional theory 
 
Some (dated) documentation for the examples and the API is [here](https://m-a-d-n-e-s-s.github.io/madness/api-doc/modules.html).


## Example for MADNESS as an external library
To use MADNESS as an external library in a code we recommend cmake. Build and install MADNESS to 
an install directory. Set

`export MADNESS_DIR=/path/to/madness/install/directory/`

In file `CMakeLists.txt`:

````
cmake_minimum_required(VERSION 3.22)
project(yourbinary)
set(CMAKE_CXX_STANDARD 17)
find_package(MADNESS CONFIG REQUIRED)
add_executable(yourbinary main.cc)
target_link_libraries(yourbinary madness)
````

In file `main.cc`:
````c
#include <madness.h>
using namespace madness;
int main(int argc, char* argv[]) {
    World& world=initialize(argc,argv);
    startup(world,argc,argv,true);
    FunctionDefaults<1>::set_cubic_cell(-10,10);
    FunctionDefaults<1>::set_k(8);
    FunctionDefaults<1>::set_thresh(1.e-6);
    try {
        auto gaussian = [](const Vector<double,1>& r){return std::exp(-r[0]*r[0]);};
        Function<double,1> f = FunctionFactory<double,1>(world).f(gaussian);
        double I = f.trace();
        double value = f(0.7);
        if (world.rank() == 0) {
          print("trace(exp(-r^2)",I,"error",I-std::sqrt(M_PI));
          print("f(0.7)",value,"error",value-exp(-0.7*0.7));
        }
    } catch (...) {
         std::cout << "caught an error " << std::endl;
    } 
    finalize();
    return 0;
}
}
````

Going through the code line by line:

````c
#include <madness.h>
````

Includes the MADNESS library.

````c
using namespace madness;
int main(int argc, char* argv[]) {
    World& world=initialize(argc,argv)
````
 
World contains the MPI communicator and must always be set, even if the code will run only 
on one node.

````c
    startup(world,argc,argv,true)
````

Reads all necessary numerical data (e.g. the wavelet twoscale coefficients) and prints out
compiler flags.

````c
    FunctionDefaults<1>::set_cubic_cell(-10,10)
````

Defines the computation cell/intervall, here in 1 dimension. MADNESS is templated with respect
to the dimensions.

````
    FunctionDefaults<1>::set_k(8);
````

Sets the wavelet order to 8. Anything between 2 and 30 will work, the optimal choice depends 
on the specific problem, 8 is usually a good guess.

````
    FunctionDefaults<1>::set_thresh(1.e-7);
````

Sets the precision threshold to $\epsilon=10^{-6}$ in the $L^2$ norm.
Unless the function is singular or has other pathological features the threshold will be met. Up to the limits of numerical precision,
it is possible to tighten the threshold to secure more digits.  There are also different truncation modes.  The default is `mode=0` which is designed to yield accurate function values and norms, whereas `mode=1` aims to also yield accurate derivatives.

````
    try {
        auto functor=[](const Vector<double,1>& r){return std::exp(-r[0]*r[0]);}
````

Defines the function to be represented in MRA. `Vector` is a MADNESS class defining coordinates.

````c
        Function<double,1> f=FunctionFactory<double,1>(world).f(functor);
````

Projects the function $e^{-r^2}$ onto the MRA representation. A `FunctionFactory` is used to 
define various properties of the `Function`, it requires `world` as input, a function (or functor)
to compute the value.  Additional optional arguments can be provided using the [named-parameter idiom](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__getting__started.html).

````c
        double I=f.trace();
````
Integrates the function $\int_{-10}^{-10} e^{-r^2}\mathrm dx$.  It is a collective operation.
````c
        double value = f(0.7);
````
evalues the numerical function at point `x=0.7.  This is executed by a single process but might involve remote communication.  If you wish to evaluate at many points (e.g., a line or a cube, there are more efficient interfaces).

````c
        if (world.rank() == 0) {
          print("trace(exp(-r^2)",I,"error",I-std::sqrt(M_PI));
          print("f(0.7)",value,"error",value-exp(-0.7*0.7));
        }
````
Prints out the values using the MADNESS Python-like `print` function.   Be sure to add the `if` block to avoid verbose output if 
run on many MPI ranks. Also make sure that the `trace` or any other collective operation is *not* called inside this
`if` block.  A collective operation called by only one rank will cause the program to hang.  This is a common error.

````c
    } catch (...) {
         std::cout << "caught an error " << std::endl;
    } 
    finalize();
````
 
Finalizes the communicator.
It is important that all MRA objects (e.g. `Function<double,1>`) are destructed before
`finalize()` is called, otherwise segmentation faults might occur since the destructor for these objects will erroneously be called at the very end of the program *after* the runtime has been dismantled.
Thus, for simple programs enclose all MRA code after `initialize` in a sub-scope (e.g., using braces or inside a `try/catch` block), or after obtaining `World` pass it (by reference) into another procedure
that contains your MRA code.

````c
    return 0;
}
````

