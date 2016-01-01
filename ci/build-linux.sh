#! /bin/sh

# Exit on error
set -ev

# Environment variables
export CXXFLAGS="-std=c++11 -mno-avx"
export CPPFLAGS=-DDISABLE_SSE3
if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
else
    # Assume CXX = clang
    export CC=/usr/bin/clang-3.6
    export CXX=/usr/bin/clang-3.6
    export CXXFLAGS="-stdlib=libc++"
    #export LDFLAGS="-fdefine-sized-deallocation"
    # clang-3.6 gives this error:
    #       clang: error: unknown argument: '-fdefine-sized-deallocation'
fi
export F77=/usr/bin/gfortran-$GCC_VERSION
# Jeff: these are unnecessary because MPICH has been compiler
#       from source with these compilers.
#export MPICH_CC=$CC
#export MPICH_CXX=$CXX
export MPICC=$HOME/mpich/bin/mpicc
export MPICXX=$HOME/mpich/bin/mpicxx
export LD_LIBRARY_PATH=/usr/lib/lapack:/usr/lib/openblas-base:$LD_LIBRARY_PATH

# Configure and build MADNESS
./autogen.sh 
mkdir build
cd build
../configure \
    --enable-debugging --enable-optimization --enable-warning --disable-optimal --disable-static \
    --with-google-test \
    --enable-never-spin \
    --with-libxc=${HOME}/libxc \
    LIBS="-L/usr/lib/lapack -L/usr/lib/libblas -llapack -lblas -lpthread"

if [ "$RUN_TEST" = "buildonly" ]; then
    # Build all libraries, examples, and applications
    make -j2 all
else
    # Run unit tests
    export MAD_NUM_THREADS=2
    make -C src/madness/$RUN_TEST -j2 check
fi
