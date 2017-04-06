#! /bin/sh

# Exit on error
set -ev

# Configure MADNESS
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE=RelWithDebInfo \
    -D ENABLE_UNITTESTS=ON \
    -D ENABLE_NEVER_SPIN=ON \
    ..

if [ "$RUN_TEST" = "buildonly" ]; then
    # Build all libraries, examples, and applications
    make -j2 all
else
    # Run unit tests
    export MAD_NUM_THREADS=2
    export CTEST_OUTPUT_ON_FAILURE=1
    make -C src/madness/$RUN_TEST -j2 test
fi
