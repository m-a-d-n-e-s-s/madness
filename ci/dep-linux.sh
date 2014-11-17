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
#sudo apt-get install -qq -y libxc-devel

if [ "$CXX" = "g++" ]; then
   sudo apt-get install -qq gcc-4.7 g++-4.7
fi

sudo apt-get install -qq -y cmake libopenblas-dev liblapack-dev  libgoogle-perftools0 mpich2 libtbb-dev
