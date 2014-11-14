#! /bin/sh

# Exit on error
set -e

# Environment variables
export MPICH_CC=$CC
export MPICH_CXX=$CXX
export MPICC=/usr/bin/mpicc.mpich2
export MPICXX=/usr/bin/mpicxx.mpich2

# Configure MADNESS
./autogen.sh 
./configure \
    --enable-debugging --disable-optimization --enable-warning --disable-optimal \
    --enable-never-spin \
    LIBS="-L/usr/lib -L/usr/lib/lapack -L/usr/lib/openblas-base -llapack -lopenblas -ltcmalloc_minimal -ltbb -lpthread"
    
make

# Run unit tests
make check