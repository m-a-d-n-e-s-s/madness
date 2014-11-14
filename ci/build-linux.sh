#! /bin/sh

# Exit on error
set -e

# Configure MADNESS
./autogen.sh 
./configure --enable-never-spin
make

# Run unit tests
make check