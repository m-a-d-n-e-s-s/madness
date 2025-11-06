# Getting Started with MADNESS

This document provides instructions for configuring, building, and testing the MADNESS project.

## Introduction

MADNESS (Multiresolution Adaptive Numerical Environment for Scientific Simulation) provides a high-level environment for the solution of integral and differential equations in many dimensions using adaptive, fast methods with guaranteed precision based on multi-resolution analysis and novel separated representations.

## Quick Start (Ubuntu)

These instructions are for a fast setup on an Ubuntu system.

### 1. Install Dependencies

Install the required packages using `apt-get`:

```bash
sudo apt-get update && sudo apt-get install -y build-essential gcc gfortran libopenblas-serial-dev cmake git ninja-build openmpi-bin openmpi-common libopenmpi-dev
```

### 2. Configure the Build

Create a build directory and run CMake with the recommended flags for a debug build:

```bash
mkdir build
cd build
cmake .. -G ninja -DCMAKE_BUILD_TYPE=Debug -DCMAKE_CXX_FLAGS_DEBUG="-O0 -g -Wall" -DENABLE_NEVER_SPIN=ON -DBUILD_SHARED_LIBS=OFF -DBUILD_TESTING=ON -DLAPACK_LIBRARIES="-L/usr/lib/x86_64-linux-gnu/openblas-serial -lopenblas -llapack"
```

### 3. Build the Code

Build the desired targets using `ninja`. To build all targets:

```bash
ninja
```

To build a specific target, like the test suite:

```bash
ninja testsuite
```

### 4. Run Tests

To run the short test suite:

```bash
ninja check-short-madness
```

To run the more detailed numerical library tests:

```bash
./src/madness/src/testsuite
```

To run the tests in parallel with MPI:

```bash
MAD_NUM_THREADS=2 mpiexec -np 2 ./src/madness/src/testsuite
```

## Detailed Configuration

For more advanced configuration, you can use various CMake variables to control the build.

### General Build Process

The general process is to create a build directory and run `cmake` to configure the project, followed by `make` or `ninja` to build it.

```bash
mkdir build
cd build
cmake /path/to/madness/source
make applications
```

### Important CMake Options

*   `CMAKE_BUILD_TYPE`: Set the build type (e.g., `Debug`, `Release`).
*   `ENABLE_MPI`: Enable or disable MPI (`ON`/`OFF`).
*   `CMAKE_INSTALL_PREFIX`: Set the installation directory.
*   `BUILD_TESTING`: Enable or disable the build of tests (`ON`/`OFF`).
*   `BUILD_SHARED_LIBS`: Build shared libraries instead of static (`ON`/`OFF`).

For a full list of options, see the `INSTALL.md` file.

### External Libraries

MADNESS can use several external libraries for enhanced performance. You can specify their locations using CMake variables.

*   **Intel MKL:** Set the `MKLROOT` environment variable or the `MKL_ROOT_DIR` CMake variable.
*   **Intel TBB:** Set the `TBBROOT` environment variable or the `TBB_ROOT_DIR` CMake variable.
*   **Other libraries:** See `INSTALL.md` for details on linking libraries like `LibXC`, `Gperftools`, `PAPI`, and others.

## Troubleshooting

### Address Space Layout Randomization (ASLR)

If you experience segmentation faults when running with multiple MPI processes, it might be due to ASLR. On platforms with ASLR, it is recommended to build MADNESS with static libraries. By default, MADNESS assumes ASLR is enabled and builds static libraries.

If you need to build shared libraries on a platform with ASLR, you may need to adjust the `MADNESS_ASSUMES_ASLR_DISABLED` and `CMAKE_POSITION_INDEPENDENT_CODE` CMake variables. See the `INSTALL.md` file for more details.
