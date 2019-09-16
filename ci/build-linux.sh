#! /bin/sh

# Exit on error 
set -ev

echo " IN BUILD SCRIPT "

pwd
ls -ltr
echo ~
ls -ltr ~

# Set up paths to stuff we installed
export MPIDIR=$HOME/mpich
export LIBXCDIR=$HOME/libxc
export PATH=${HOME}/ccache/bin:$PATH

export PATH=$MPIDIR/bin:$PATH
export CC=mpicc
export CXX=mpicxx
export FC=gfortran-8

echo $CC
which $CC
$CC --version
echo $CXX
which $CXX
$CXX --version
echo $FC
which $FC
$FC --version
echo $LIBXCDIR
ls $LIBXCDIR
ls $LIBXCDIR/lib
ls -l $LIBXCDIR/lib/libxc.a
ls -l /opt/intel
which ccache
ccache --version

# Configure MADNESS 
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE=MinSizeRel \
    -D ENABLE_UNITTESTS=ON \
    -D ENABLE_NEVER_SPIN=ON \
    -D BUILD_SHARED_LIBS=OFF \
    -D ENABLE_GPERFTOOLS=OFF \
    -D ENABLE_MKL=ON \
    -D CMAKE_C_COMPILER=$CC \
    -D CMAKE_CXX_COMPILER=$CXX \
    -DLIBXC_LIBRARIES=$LIBXCDIR/lib/libxc.a \
    -DLIBXC_INCLUDE_DIRS=$LIBXCDIR/include \
    $CMAKE_EXTRA_OPTIONS \
    ..

export MAD_SMALL_TESTS=1
export MAD_NUM_THREADS=2
export CTEST_OUTPUT_ON_FAILURE=1

if [ "$RUN_TEST" = "buildonly" ]; then
    echo "Build all libraries, examples, and applications ($CMAKE_EXTRA_OPTIONS)"
    make -j2 all VERBOSE=1
elif [ "$RUN_TEST" = "all" ]; then
    echo "Build all libraries, examples, and applications; then run all test programs ($CMAKE_EXTRA_OPTIONS)"
    make -j2 all  VERBOSE=1
    make -j2 test ARGS="-V"
else
    echo "Running test --- $RUN_TEST ($CMAKE_EXTRA_OPTIONS)"
    make -C src/madness/$RUN_TEST -j2 test ARGS="-V"
fi
ccache -s
