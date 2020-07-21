#! /bin/sh

# Exit on error 
set -ev

echo " IN BUILD SCRIPT "

pwd
ls -ltr
echo ~
ls -ltr ~

# enable MKL direct call for clang only
if [ "$CXX" = "clang++" ]; then
    export EXTRACXXFLAGS="-DMKL_DIRECT_CALL"
fi

# add cereal header files to include dirs
if [ -f "${HOME}/cereal/include/cereal/cereal.hpp" ]; then
    export EXTRACXXFLAGS="${EXTRACXXFLAGS} -I${HOME}/cereal/include"
fi

# Set up paths to stuff we installed
export MPIDIR=$HOME/mpich
export LIBXCDIR=$HOME/libxc
export PATH=${HOME}/ccache/bin:$PATH

export PATH=$MPIDIR/bin:$PATH
export CC=mpicc
export CXX=mpicxx
export FC=gfortran-8

if [ "X${BUILD_SHARED}" = "X1" ]; then
  LIBEXT="so"
else
  LIBEXT="a"
fi

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
ls -l $LIBXCDIR/lib/libxc.${LIBEXT}
ls -l /opt/intel
which ccache
ccache --version
ls -l /opt/intel/tbb/bin/*
#source /opt/intel/tbb/bin/tbbvars.sh intel64
echo ${TBBROOT}

# Configure MADNESS
mkdir build
cd build
cmake \
    -D CMAKE_BUILD_TYPE=MinSizeRel \
    -D ENABLE_UNITTESTS=ON \
    -D ENABLE_NEVER_SPIN=ON \
    -D ENABLE_GPERFTOOLS=OFF \
    -D ENABLE_MKL=ON \
    -D ENABLE_TBB=ON \
    -D CMAKE_C_COMPILER=$CC \
    -D CMAKE_CXX_COMPILER=$CXX \
    -D CMAKE_CXX_FLAGS="${EXTRACXXFLAGS}" \
    -D LIBXC_LIBRARIES=$LIBXCDIR/lib/libxc.${LIBEXT} \
    -D LIBXC_INCLUDE_DIRS=$LIBXCDIR/include \
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
