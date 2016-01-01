#! /bin/sh

# Exit on error
set -ev

# Add repository for libxc
#sudo add-apt-repository ppa:hogliux/misstep -y

# Add PPA for a newer version GCC
#sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y
# Add PPA for newer cmake (3.2.3)
#sudo add-apt-repository ppa:george-edison55/precise-backports -y

# Update package list
#sudo apt-get update -qq

# Install packages

#sudo apt-get install -qq -y gcc-$GCC_VERSION g++-$GCC_VERSION gfortran-$GCC_VERSION
if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
else
    # Assume CXX = clang
    export CC=/usr/bin/clang-3.6
    export CXX=/usr/bin/clang-3.6
    export LDFLAGS="-fdefine-sized-deallocation"
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
#sudo apt-get -y -qq --no-install-suggests --no-install-recommends --force-yes install cmake cmake-data
cmake --version

# Install the rest
#sudo apt-get install -qq -y libblas-dev liblapack-dev libgoogle-perftools-dev mpich2 libtbb-dev

if [ ! -d "${HOME}/mpich" ]; then
    wget --no-check-certificate -q http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
    tar -xzf mpich-3.2.tar.gz
    cd mpich-3.2
    mkdir build && cd build
    ../configure CC=$CC CXX=$CXX --disable-fortran --disable-romio --prefix=${HOME}/mpich
    make -j2
    make install
else
    echo "MPICH installed..."
    find ${HOME}/mpich -name mpiexec
    find ${HOME}/mpich -name mpicc
fi
