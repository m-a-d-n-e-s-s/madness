#! /bin/sh

# Exit on error
set -e

brew update
brew install autoconf automake libtool cmake libxc mpich2 tbb gperftools
