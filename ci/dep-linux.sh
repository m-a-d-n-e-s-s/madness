#! /bin/sh

# Exit on error
set -ev

# Install packages
echo " IN DEP SCRIPT "
env
pwd
ls -ltr
echo ~
ls -ltr ~

# Verify CMake 3
cmake --version

# Install dependencies (libxc, mpich)

export FC=/usr/bin/gfortran-8

# Set up CC and CXX commands
export COMPILER=$CXX
case "$COMPILER" in
    g++)
        echo "setting up for GCC version 8"
        export CC=/usr/bin/gcc-8
        export CXX=/usr/bin/g++-8
        ;;
    clang++)
        echo "setting up for CLANG version 7"

        which clang || true
        which clang++ || true
        which clang++-7.0 || true
        export CC=clang
        export CXX=clang++
        ;;
    *)
        echo "Unknown C++ compiler:"
        echo "$CXX"
        exit 1
        ;;
esac

# Confirm compiler information
$CC --version
$CXX --version
$FC --version


# Install ccache
if [ ! -f "${HOME}/ccache/bin/ccache" ]; then
    wget https://www.samba.org/ftp/ccache/ccache-3.5.tar.bz2
    tar xf ccache-3.5.tar.bz2 
    cd ccache-3.5/
    ./configure --prefix=${HOME}/ccache
    make -j 2
    make install
    cd ..
else
    echo "ccache already installed"
    ls -l ${HOME}/ccache/bin/
fi
export PATH=${HOME}/ccache/bin:$PATH
ccache --version || true

# Configure ccache
if [ ! -f "${HOME}/.ccache/ccache.conf" ]; then
    mkdir ${HOME}/.ccache
    cat <<EOF > ${HOME}/.ccache/ccache.conf
hash_dir = false
cache_dir_levels = 8
compiler_check = %compiler% --version
compression = true
run_second_cpp = true
EOF
fi

# Install libxc
if [ ! -f "${HOME}/libxc/lib/libxc.a" ]; then
    export LIBXC_VERSION=4.3.4
    wget -O libxc-${LIBXC_VERSION}.tar.gz "https://gitlab.com/libxc/libxc/-/archive/${LIBXC_VERSION}/libxc-${LIBXC_VERSION}.tar.gz"
    tar -xzf libxc-${LIBXC_VERSION}.tar.gz
    ls -l
    cd libxc-${LIBXC_VERSION}
    autoreconf -i
    ./configure --prefix=${HOME}/libxc --enable-static --disable-fortran CFLAGS="-mno-avx -O1" CXXFLAGS="-mno-avx -O1" FCFLAGS="-mno-avx -O1"
    make -j2
    make install
    cd ..
    rm -rf libxc-${LIBXC_VERSION}
else
    echo "LIBXC installed..."
    ls -l ${HOME}/libxc
fi

# Install MPICH
if [ ! -f "${HOME}/mpich/bin/mpicc" ]; then
    wget --no-check-certificate -q http://www.mpich.org/static/downloads/3.2/mpich-3.2.tar.gz
    tar -xzf mpich-3.2.tar.gz
    cd mpich-3.2
    ./configure CC=$CC CXX=$CXX --disable-fortran --disable-romio --with-pm=gforker --prefix=${HOME}/mpich CFLAGS="-mno-avx -O0"
    make -j2
    make install
    ${HOME}/mpich/bin/mpichversion
    ${HOME}/mpich/bin/mpicc -show
    ${HOME}/mpich/bin/mpicxx -show
    cd ..
    rm -rf mpich-3.2
else
    echo "MPICH installed..."
    find ${HOME}/mpich -name mpiexec
    find ${HOME}/mpich -name mpicc
    find ${HOME}/mpich -name mpicxx
fi

# Install Cereal (headers only)
if [ ! -f "${HOME}/cereal/include/cereal/cereal.hpp" ]; then
    cd
    git clone https://github.com/USCiLab/cereal cereal_repo
    cd cereal_repo
    git checkout a5a30953125e70b115a2
    mkdir build
    cd build
    cmake -D JUST_INSTALL_CEREAL=ON -D CMAKE_INSTALL_PREFIX=${HOME}/cereal ..
    make install
    cd 
    rm -rf cereal_repo
else
    echo "Cereal installed..."
    find ${HOME}/cereal/include/cereal -name "cereal.hpp"
fi

# Do not exit on error because MKL is optional
set +e

# Install MKL
wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
sudo apt-get update
sudo apt-get install intel-mkl-2019.4-070

