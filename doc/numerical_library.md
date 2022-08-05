# MADNESS numerical library

MADNESS provides a high-level environment for the solution of integral and differential equations 
in many dimensions using adaptive, fast methods with guaranteed precision based on multi-resolution 
analysis and novel separated representations. 

## Example for MADNESS as an external library
To use MADNESS as an external library in a code we recommend cmake. Build and install MADNESS to 
an install directory. Set

`export MADNESS_DIR=/path/to/madness/install/directory/`


In file `CMakeLists.txt`:
> `cmake_minimum_required(VERSION 3.22)`\
> `project(yourbinary)`\
> `set(CMAKE_CXX_STANDARD 17)`\
> `find_package(MADNESS CONFIG REQUIRED)`\
> `add_executable(yourbinary main.cpp)`\
> `target_link_libraries(yourbinary madness)`

In file `main.cpp`:
>`#define USE_GENTENSOR 1` // only needed if madness was configured with `-D ENABLE_GENTENSOR=1`\
>`#include <madness.h>`\
>`using namespace madness;`\
>`int main(int argc, char* argv[]) {`\
>`  World& world=initialize(argc,argv);`\
>`  startup(world,argc,argv,true);`\
>`  FunctionDefaults<1>::set_cubic_cell(-10,10);`\
>`  FunctionDefaults<1>::set_k(8);`\
>`  try {`\
>`    auto functor=[](const Vector<double,1>& r){return exp(-r[0]*r[0]);};`\
>`    Function<double,1> f=FunctionFactory<double,1>(world).f(functor);`\
>`    double I=f.trace();`\
>`    std::cout << "trace(exp(-r^2) " << I << std::endl;`\
>`  } catch (...) {`\
>`     std::cout << "caught an error " << std::endl;`\
>`  } `\
>`  finalize();`\
>`  return 0;`\
>`}`

Going through the code line by line:
>`#define USE_GENTENSOR 1` // only needed if madness was configured with `-D ENABLE_GENTENSOR=1`\
 
If MADNESS is to be used for high-dimensional problems (d>3) you will probably use low-rank tensor 
approximations, which must be set at MADNESS configure time with the `-D ENABLE_GENTENSOR=1` flag, 
and then again in your code. GenTensor and complex Functions are currently exlusive.

>`#include <madness.h>`

Includes the MADNESS library.

>`using namespace madness;`\
>`int main(int argc, char* argv[]) {`\
>`  World& world=initialize(argc,argv);`
 
World contains the MPI communicator and must always be set, even if the code will run only 
on one node.
 
>`  startup(world,argc,argv,true);`

Reads all necessary numerical data (e.g. the wavelet twoscale coefficients) and prints out
compiler flags.

>`  FunctionDefaults<1>::set_cubic_cell(-10,10);`

Defines the computation cell/intervall, here in 1 dimension. MADNESS is templated with respect
to the dimensions.
 
>`  FunctionDefaults<1>::set_k(8);`

Sets the wavelet order to 8. Anything between 2 and 30 will work, the optimal choice depends 
on the specific problem, 8 is usually a good guess.


>`  try {`\
>`    auto functor=[](const Vector<double,1>& r){return exp(-r[0]*r[0]);};`

Defines the function to be represented in MRA. `Vector` is a MADNESS class defining coordinates.

>`    Function<double,1> f=FunctionFactory<double,1>(world).f(functor);`

Projects the function $e^{-r^2}$ onto the MRA representation. A `FunctionFactory` is used to 
define various properties of the `Function`, it requires `world` as input, optionally a function/lambda
defining the mathematical function. 

>`    double I=f.trace();`

Integrate the function $\int_{-10}^{-10} e^{-r^2}\mathrm dx$..

>`    if (world.rank()==0) std::cout << "trace(exp(-r^2) " << I << std::endl;`
 
Prints out the value of the integral. Be sure to add the `if` block to avoid verbose output if 
run on many MPI ranks. Also make sure that the `trace` operation is not called inside this
`if` block, because it is a collective operation of all MPI ranks and having it executed only
by one rank will cause the program to hang. This is a common error.

>`  } catch (...) {`\
>`     std::cout << "caught an error " << std::endl;`\
>`  } `\
>`  finalize();`\
 
Finalizes the communicator.
It is important that all MRA objects (e.g. Function<double,1>) are destructed before
`finalize()` is called, otherwise segmentation faults will occur,
so best enclose all MRA code after `startup` inside a `try/catch` block.

>`  return 0;`\
>`}`
