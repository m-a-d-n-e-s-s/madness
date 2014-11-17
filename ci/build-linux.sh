#! /bin/sh

# Exit on error
set -ev

# Environment variables
if [ "$CXX" = "g++" ]; then
    export CC=gcc-4.7
    export CXX=g++-4.7
    export F77=gfortran-4.7
fi
export MPICH_CC=$CC
export MPICH_CXX=$CXX
export MPICC=/usr/bin/mpicc.mpich2
export MPICXX=/usr/bin/mpicxx.mpich2
export LD_LIBRARY_PATH=/usr/lib/lapack:/usr/lib/openblas-base:$LD_LIBRARY_PATH

# Configure MADNESS
./autogen.sh 
./configure \
    --enable-debugging --disable-optimization --enable-warning --disable-optimal \
    --enable-never-spin \
    CXXFLAGS="-std=c++11" \
    LIBS="-L/usr/lib/lapack -L/usr/lib/openblas-base -llapack -lopenblas -lpthread"
    
make

# Run unit tests
make check