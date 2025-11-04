# Agent Instructions for MADNESS


This document provides instructions for agents working with the MADNESS codebase from the github repo `git@github.com:m-a-d-n-e-s-s/madness.git`

## Dependencies

The following packages are required to build and test the project. They should be installed using `apt-get`:
- `build-essential`
- `gcc`
- `gfortran`
- `libopenblas-serial-dev`
- `cmake`
- `git`
- `ninja-build`
- `openmpi-bin`
- `openmpi-common`
- `libopenmpi-dev`

Example installation command:
```
sudo apt-get update && sudo apt-get install -y build-essential gcc gfortran libopenblas-serial-dev cmake git ninja-build openmpi-bin openmpi-common libopenmpi-dev
```

## Configuration

To configure the project for fast builds and debugging, create a `build` directory in the root of the repository, change into it, and then run `cmake` with the following flags:
```
mkdir build
cd build
cmake .. -G ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG="-O0 -g -Wall" -DENABLE_NEVER_SPIN=ON -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=ON -DLAPACK_LIBRARIES="-L/usr/lib/x86_64-linux-gnu/openblas-serial -lopenblas -llapack"
```

## Building

After configuring the project, you can build the targets using `ninja`. For example, to build the `testsuite` target:
```
ninja testsuite
```

## Testing

To build and run the short test suite, from the build directory execute
```
ninja check-short-madness
```

To build the detailed tests of the numerical library, from the build directory execute
```
ninja testsute
```
After a successful build, run these tests with
```
./src/madness/src/testsuite
```
Success will print the message "testsuite passed:  true" at the end and will have exit code zero.

To run the same tests in parallel using MPI and two processes with two threads each, use the following command
```
MAD_NUM_THREADS=2 mpiexc -np 2 ./src/madness/src/testsuite
```
Success will print the message "testsuite passed:  true" at the end and will have exit code zero.

To build the molecular density functional (DFT) code, from the build directory execute
```
ninja moldft
```
After a successful build, from the build directory test the DFT code with the command
```
./src/apps/moldft/moldft --geometry=water --dft="xc=lda"
```
Success will have exit code zero.

