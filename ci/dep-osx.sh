#! /bin/sh

# Exit on error
set -e

brew update

#
# Workaround Homebrew issues in Travis CI...
#
#==> Installing gcc
#==> Downloading https://homebrew.bintray.com/bottles/gcc-7.2.0.sierra.bottle.tar.gz
#==> Pouring gcc-7.2.0.sierra.bottle.tar.gz
#Error: The `brew link` step did not complete successfully
#The formula built, but is not symlinked into /usr/local
#Could not symlink include/c++
#Target /usr/local/include/c++
#already exists. You may want to remove it:
#  rm '/usr/local/include/c++'
brew upgrade gcc || brew install gcc || true
brew link --overwrite gcc || true

# Proceed assuming GCC is installed properly.
brew install libxc mpich tbb
