# MADNESS tutorial

## Table of contents

1. Building MADNESS
   1. [Building from source](#building-from-source)
       1. [Downloading](#downloading)
       1. [Building and installing](#building-and-installing)
   1. [Use from a Docker container](#use-from-a-docker-container)
1. [Chemistry applications](#chemistry)
1. [Parallel runtime --- basic concepts and architecture](#Parallel-runtime)
1. [Numerical API and example application](#Numerical-API)
1. [Chemical API and example Hartree-Fock program](#Chemical-API)
1. [Exercises](#Exercises)
  
## Building MADNESS

There are two ways to obtain/build MADNESS:
- download the MADNESS source code from the Github repo and build, or
- use a Docker image containing pre-built MADNESS.

### Building from source

#### Downloading

From the command line, clone the [MADNESS GitHub repository](https://github.com/m-a-d-n-e-s-s/madness) using one of the below
* HTTPS
```
    git clone https://github.com/m-a-d-n-e-s-s/madness.git
```
* SSH (for which you will need to have a GitHub account with SSH key)
```
    git clone git@github.com:m-a-d-n-e-s-s/madness.git
```
This will make a new directory `madness` --- if you wish use a different name, append it to the command.

The head of the master branch usually builds and runs correctly due to the continuous integration.  If you are looking for the version tested for the Tromso June 2023 workshop, after cloning 
```
cd madness
checkout XXXXXXXXXXXXXXXXXXX
```

#### Building and installing

Please refer to the [installation instructions](https://madness.readthedocs.io/en/latest/INSTALL.html).

Minimal recipe without MPI and assuming Intel MKL is installed in a standard location
```
mkdir build # CANNOT be in the madness source tree
cd build
cmake -DENABLE_MPI=OFF -DCMAKE_INSTALL_PREFIX=/home/me/madinstall /path/to/madness/source
make applications
make install  # optional
```

#### Modules to load to build on FRAM for the Trømso tutorial

```
    module purge
    module load Emacs/27.2-GCCcore-11.2.0
    module load OpenMPI/4.1.4-GCC-11.3.0
    module load imkl/2022.2.1
    module load tbb/2021.5.0-GCCcore-11.3.0
    module load CMake/3.23.1-GCCcore-11.3.0
    module load git/2.36.0-GCCcore-11.3.0-nodocs
    module load Python/3.10.4-GCCcore-11.3.0
```
Ignore any warnings or informational messages.

### Use from a Docker container

It is also possible to use a Docker image with pre-built MADNESS. Although this method is not recommended
for optimal use of MADNESS on distributed-memory machines, it is sufficient to run MADNESS applications on
a single multicore machine.

To get started, make sure your machine has [a Docker engine](https://docs.docker.com/engine/install/). Once the engine is installed (and, on MacOS, running)
you can start a container containing the pre-built MADNESS:
```shell
docker run -it --rm rjharrison/ubuntu:22.04 bash
```
This will start a shell within the container and put you in the directory containing the MADNESS source (`source`),
the build directory with MADNESS libraries and applications  already configured and built (`build`), and
the directory with MADNESS libraries and applications installed (`install`).
E.g., to run a `moldft` calculation use the installed binary:
```shell
MAD_NUM_THREADS=2 ./install/bin/moldft --geometry=he
```

You can also run other commands directly, e.g. to run tests execute:
```shell
docker run -it --rm rjharrison/ubuntu:22.04 cmake --build build --target check-madness
```
It is convenient to execute MADNESS applications this way using input files located on the host file system;
to run `moldft` with input file `my.input` located in the current working directory
do
```shell
docker run -it --rm -v `pwd`:/pwd rjharrison/ubuntu:22.04 install/bin/moldft /pwd/my.input
```
P.S. Note that unfortunately it is not possible to stop this process using Ctrl+C, see some workarounds [here](https://forums.docker.com/t/docker-run-cannot-be-killed-with-ctrl-c/13108/10).

## Chemistry

* [Please look here](chemistry.md)

## Parallel runtime

Please refer to slides 166-189 in [MADNESSeverything4.pdf](https://github.com/m-a-d-n-e-s-s/madness/blob/tutorial/doc/MADNESSeverything4.pdf).
* Since the file is large, it is probably easier to look at your local version using your system PDF viewer

Hello world in MADNESS (runtime only --- numerical library is discussed below)
```
    #include <madness/world/MADworld.h>
    using namespace madness;

    int main(int argc, char** argv) {
        initialize(argc,argv);
        World world(SafeMPI::COMM_WORLD);

        print("Hello from processor",world.rank());

        finalize();
        return 0;
    }
```

More detailed (but dated) documentation
* https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__parallel__runtime.html


## Numerical API

Please refer to the [Numerical Library User Documentation](https://madness.readthedocs.io/en/latest/numerical_library.html) for a worked example.

* Setup up your project files 
  * Make a new directory (`test`) with two new sub-directories `src` and `build`
    * `mkdir test`
    * `cd test`
    * `mkdir src build`
  * Copy and paste the example cmake code into `src/CMakeLists.txt`
  * Copy and paste the example C++ into `src/main.cc`
* Setup your environment so cmake can find MADNES
  * `export MADNESS_DIR=/path/to/madness/install/directory/`
* Build the code
  * `cd build`
  * `cmake ../src`
  * `make`
* Execute
  * `MAD_NUM_THREADS=2 ./yourbinary`
  
[Developer documentation](https://m-a-d-n-e-s-s.github.io/madness/api-doc/) generated by doxygen
* Needs extensive updating and will (eventually) be mostly superseded by readthedocs

If you are going to program more extensively in MADNESS, then also worth looking at 
* [Getting started](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__getting__started.html)
* [Tensor library](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__tensor.html)
* [Linear algebra](https://m-a-d-n-e-s-s.github.io/madness/api-doc/tensor__lapack_8h.html)
* [MRA functions](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__gstart__functions.html) and [here](https://m-a-d-n-e-s-s.github.io/madness/api-doc/mra_8h.html)
  * [Function factory](https://m-a-d-n-e-s-s.github.io/madness/api-doc/classmadness_1_1FunctionFactory.html#details)
  * [Function defaults](https://m-a-d-n-e-s-s.github.io/madness/api-doc/classmadness_1_1FunctionDefaults.html)
  * [Operations on vectors of functions](https://m-a-d-n-e-s-s.github.io/madness/api-doc/vmra_8h.html)
* [Serialization](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__serialization.html) and [here]()

## Chemical API

* [Please look here](API.md)

## Exercises

### Chemistry examples
* Run some calculations on small molecules of either your choosing or from the [structure library](https://github.com/m-a-d-n-e-s-s/madness/blob/master/src/madness/chem/structure_library).  E.g.,
```shell
    MAD_NUM_THREADS=10 moldft --geometry="water" --dft="xc lda; maxiter 5"
```
  * Without linking with LIBXC, with `moldft` you can run either LDA or HF calculations.
  * Play with some of the `moldft` DFT parameters
    * E.g., the sequence of tolerances used for DFT solution `protocol [1e-4,1e-6,1e-8]` and the error in the energy.
    * E.g., the interaction between convergence requested for the density (e.g., `dconv 1e-6`), the error in the energy, and the required tolerances (you'll find tighter convergence will need tighter tolerance to get convergence).
    * E.g., the smoothing of the nuclear potential (e.g., `eprec 1e-4` in the `geometry` block) and its impact on the cost of the calculation, the error in the energy, or the error in an optimized geometry.
  * Compare results with your favorite other code using a different basis (e.g., Gaussian functions).

### Explore the example codes
* The source code is [here](https://github.com/m-a-d-n-e-s-s/madness/tree/master/src/examples)
* The doxygen generated documentation
  * [Solves the 3D harmonic oscillator](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__example3dharm.html)
  * [Illustrates general composition of two functions](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__examplebinop.html)
  * [Data and load balancing](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__loadbaleg.html)
  * [Poisson's equation in a dielectric medium](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__exampledielectric.html)
  * [Laplace's equations for dielectric sphere in an external field](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__exampledielectricfield.html)
  * [Example of function I/O from getting started guide](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__functionioeg.html)
  * [Compute the dielectric cavity and the electrostatic potential of hydrogen atom in water](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__examplegygi.html)
  * [Hartree-Fock equations for the hydrogen molecule](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__examplesh2hf.html)
  * [Energy of the hydrogen atom ground state](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__hatom__energy.html)
  * [Solves heat equation using the Green's function](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__exampleheat.html)
  * [Evolve in time 3D heat equation with a linear term](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__heatex2.html)
  * [Hartree-Fock equations for the helium atom](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__examplehehf.html)
  * [Solves the two-particle system exactly](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__helium__exact.html)
  * [Hello world MADNESS style](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__hellowworldmad.html)
  * [Solves a Navier-Stokes equation](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__examplense.html)
  * [Solves a 1D nonlinear Schrödinger equation](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__examplenonlinsc.html)
  * [Demonstrates/tests use of 3D shape functions](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__shape__tester.html)
  * [First example from getting started guide](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__sininteg.html)
  * [Spectral propagator in time using semigroup approach](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__spectralprop.html)
  * [Solves a 1D time-dependent Schrödinger equation using splitting and semi-group approaches with the free-particle propagator](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__exampletdse1d.html)
 
### Write your own simple test
 
Starting from the simple `cmake` and C++ file from the [Numerical API](#Numerical-API) section, make a code to solve a simple problem you are interested in
* E.g., project and do arithmetic on a function of your selection.
* E.g., plot a function (look at some of the examples and also the [intrductory documentation](https://m-a-d-n-e-s-s.github.io/madness/api-doc/group__gstart__io.html)

### MADNESS detailed presentation

Have a skim through [MADNESSeverything4.pdf](https://github.com/m-a-d-n-e-s-s/madness/blob/tutorial/doc/MADNESSeverything4.pdf).
* Since the file is large, it is probably easier to look at your local version using your system PDF viewer.
* It goes through all the key concepts and is hopefully useful for self study.




 
 
