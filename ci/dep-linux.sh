#! /bin/sh

# Exit on error
set -ev

# Install packages

if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
else
    # Assume CXX = clang
    export CC=/usr/bin/clang-3.6
    export CXX=/usr/bin/clang-3.6
    export CXXFLAGS="-std=c++11"
    export LDFLAGS="-std=c++11"
fi
export FC=/usr/bin/gfortran-$GCC_VERSION

# Print compiler information
$CC --version
$CXX --version
$FC --version

# Install libxc
wget -O libxc-2.2.1.tar.gz "http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.2.1.tar.gz"
tar -xzf libxc-2.2.1.tar.gz
cd libxc-2.2.1
autoreconf -i
./configure --prefix=${HOME}/libxc --enable-shared --disable-static CFLAGS="-mno-avx" FCFLAGS="-mno-avx"
make -j2
make install

# Install CMake 3
cmake --version

# Install MPICH
if [ ! -d "${HOME}/mpich" ]; then
    wget --no-check-certificate -q http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
    tar -xzf mpich-3.2.tar.gz
    cd mpich-3.2
    mkdir build && cd build
    ../configure CC=$CC CXX=$CXX --disable-fortran --disable-romio --prefix=${HOME}/mpich
    make -j2
    make install
    ${HOME}/mpich/bin/mpichversion
    ${HOME}/mpich/bin/mpicc -show
    ${HOME}/mpich/bin/mpicxx -show
else
    echo "MPICH installed..."
    find ${HOME}/mpich -name mpiexec
    find ${HOME}/mpich -name mpicc
    find ${HOME}/mpich -name mpicxx
fi
