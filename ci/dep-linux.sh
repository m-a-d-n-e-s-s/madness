#! /bin/sh

# Exit on error
set -ev

# Add repository for libxc
#sudo add-apt-repository ppa:hogliux/misstep -y

# Add repository for a newer version GCC
sudo add-apt-repository ppa:ubuntu-toolchain-r/test -y

# Update package list
sudo apt-get update -qq

# Install packages

sudo apt-get install -qq -y gcc-$GCC_VERSION g++-$GCC_VERSION gfortran-$GCC_VERSION
if [ "$CXX" = "g++" ]; then
    export CC=/usr/bin/gcc-$GCC_VERSION
    export CXX=/usr/bin/g++-$GCC_VERSION
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
./configure --prefix=/usr/local CFLAGS="-mno-avx" FCFLAGS="-mno-avx"
make -j2
sudo make install

sudo apt-get install -qq -y cmake libopenblas-dev liblapack-dev libgoogle-perftools-dev mpich2 libtbb-dev
