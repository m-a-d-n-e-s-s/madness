#! /bin/sh

# OSX tests presently not supported since Travis OSX environment is so problematic
exit 1

# Exit on error
set -ev

CC=gcc
CXX=clang

$CC --version
$CXX --version

# Configure MADNESS
mkdir build
cd build
export LDFLAGS=`cmake --find-package -DNAME=clang -DCOMPILER_ID=GNU -DLANGUAGE=CXX -DMODE=LINK`
echo $LDFLAGS
cmake \
    -D CMAKE_BUILD_TYPE=RelWithDebInfo \
    -D ENABLE_UNITTESTS=ON \
    -D ENABLE_NEVER_SPIN=ON \
    -D CMAKE_C_COMPILER=$CC \
    -D CMAKE_CXX_COMPILER=$CXX \
    $CMAKE_EXTRA_OPTIONS \
    ..

if [ "$RUN_TEST" = "deponly" ]; then
    echo "Build dependencies only --- no source compiled or tests being run"
elif [ "$RUN_TEST" = "buildonly" ]; then
    # Build all libraries, examples, and applications
    make -j2 all
else
    # Run unit tests
    export MAD_NUM_THREADS=2
    export CTEST_OUTPUT_ON_FAILURE=1
    # Verbose output to reassure Travis that stuff is happening
    make -C src/madness/$RUN_TEST -j2 test ARGS="-V"
fi
