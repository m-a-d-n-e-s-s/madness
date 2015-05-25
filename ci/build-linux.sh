#! /bin/sh

# Exit on error
set -ev

# Environment variables
export CXXFLAGS="-std=c++11 -mno-avx"
export CPPFLAGS=-DDISABLE_SSE3
if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
fi
export F77=/usr/bin/gfortran-$GCC_VERSION
export MPICH_CC=$CC
export MPICH_CXX=$CXX
export MPICC=/usr/bin/mpicc.mpich2
export MPICXX=/usr/bin/mpicxx.mpich2
export LD_LIBRARY_PATH=/usr/lib/lapack:/usr/lib/openblas-base:$LD_LIBRARY_PATH

# Configure and build MADNESS
./autogen.sh 
./configure \
    --enable-debugging --disable-optimization --enable-warning --disable-optimal \
    --with-google-test \
    --enable-never-spin \
    LIBS="-L/usr/lib/lapack -L/usr/lib/libblas -llapack -lblas -lpthread"
    
make -j2 libraries

# Run unit tests
if [ "$RUN_TEST" != "buildonly" ]; then
    export MAD_NUM_THREADS=4
    make -C src/madness/$RUN_TEST -j2 check # run all tests, even if some fail
fi
