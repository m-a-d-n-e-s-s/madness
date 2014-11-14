#! /bin/sh

# Exit on error
set -e

sudo add-apt-repository ppa:hogliux/misstep -y
sudo apt-get update -qq
sudo apt-get install -qq -y cmake libopenblas-dev liblapack-dev  libgoogle-perftools0 mpich2 libxc-devel libtbb-dev