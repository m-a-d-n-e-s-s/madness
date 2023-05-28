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

In file `main.cpp`:
````c
#include <madness.h>
using namespace madness;
int main(int argc, char* argv[]) {
    World& world=initialize(argc,argv);
    startup(world,argc,argv,true);
    FunctionDefaults<1>::set_cubic_cell(-10,10);
    FunctionDefaults<1>::set_k(8);
    FunctionDefaults<1>::set_thresh(1.e-5);
    try {
        auto gaussian=[](const Vector<double,1>& r){return exp(-r[0]*r[0]);};
        Function<double,1> f=FunctionFactory<double,1>(world).f(gaussian);
        double I=f.trace();
        std::cout << "trace(exp(-r^2) " << I << std::endl;
    } catch (...) {
         std::cout << "caught an error " << std::endl;
    } 
    finalize();
    return 0;
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
    FunctionDefaults<1>::set_thresh(1.e-5);
````

Sets the precision threshold to $\epsilon=10^{-5}$ in the $L^2$ norm.
Unless the function is singular or has other pathological features the threshold will be met. It is always 
possible to tighten the threshold to secure more digits.

````
    try {
        auto functor=[](const Vector<double,1>& r){return exp(-r[0]*r[0]);}
````

Defines the function to be represented in MRA. `Vector` is a MADNESS class defining coordinates.

````c
        Function<double,1> f=FunctionFactory<double,1>(world).f(functor);
````

Projects the function $e^{-r^2}$ onto the MRA representation. A `FunctionFactory` is used to 
define various properties of the `Function`, it requires `world` as input, optionally a function/lambda
defining the mathematical function. 

````c
        double I=f.trace();
````

Integrates the function $\int_{-10}^{-10} e^{-r^2}\mathrm dx$.

````c
        if (world.rank()==0) std::cout << "trace(exp(-r^2) " << I << std::endl;
````
 
Prints out the value of the integral. Be sure to add the `if` block to avoid verbose output if 
run on many MPI ranks. Also make sure that the `trace` operation is not called inside this
`if` block, because it is a collective operation of all MPI ranks and having it executed only
by one rank will cause the program to hang. This is a common error.

````c
    } catch (...) {
         std::cout << "caught an error " << std::endl;
    } 
    finalize();
````
 
Finalizes the communicator.
It is important that all MRA objects (e.g. `Function<double,1>`) are destructed before
`finalize()` is called, otherwise segmentation faults will occur,
so best enclose all MRA code after `startup` inside a `try/catch` block.

````c
    return 0;
}
````
