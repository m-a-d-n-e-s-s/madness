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

if [ "$CXX" = "g++" ]; then
    sudo apt-get install -qq -y gcc-4.7 g++-4.7 gfortran-4.7
    export CC=gcc-4.7
    export CXX=g++-4.7
    export FC=gfortran-4.7
else
    sudo apt-get install -qq -y gfortran
    export FC=gfortran
fi

wget -O libxc-2.2.1.tar.gz "http://www.tddft.org/programs/octopus/down.php?file=libxc/libxc-2.2.1.tar.gz"
tar -xzf libxc-2.2.1.tar.gz
cd libxc-2.2.1
autoreconf -i
./configure --prefix=/usr/local
make
sudo make install

sudo apt-get install -qq -y cmake libopenblas-dev liblapack-dev libgoogle-perftools-dev mpich2 libtbb-dev
