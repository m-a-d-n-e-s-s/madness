#! /bin/sh

# OSX tests presently not supported since Travis OSX environment is so problematic
exit 1

# Exit on error
set -e

# Configure ccache
if [ ! -f "${HOME}/.ccache/ccache.conf" ]; then
    mkdir ${HOME}/.ccache
    cat <<EOF > ${HOME}/.ccache/ccache.conf
cache_dir_levels = 8
compiler_check = %compiler% --version
compression = true
run_second_cpp = true
EOF
fi

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

brew install --HEAD ccache || true
ccache --version || true

# TBB requires python@2 but it's plagued by same issues as gcc
# use similare workaround
brew upgrade python@2 || brew install python@2 || true
brew link --overwrite python@2 || true

# Proceed assuming GCC is installed properly.
brew install libxc mpich tbb cereal
