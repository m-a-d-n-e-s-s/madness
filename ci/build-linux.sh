#! /bin/sh

# Exit on error
set -e

# Configure MADNESS
./autogen.sh 
./configure \
    --enable-debugging --disable-optimization --enable-warning --disable-optimal \
    --enable-never-spin \
    LIBS="-llapack -lopenblas" \
    MPICC=/usr/bin/mpicc.mpich2 MPICXX=/usr/bin/mpicxx.mpich2
make

# Run unit tests
make check