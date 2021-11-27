# check if conda is there
if ! command -v conda &> /dev/null
then
    echo "conda is not installed"
    exit 1
fi

set -e

# install compilers (linux only: others need to install manually)
echo building on $(uname -a)
os=$(uname -a | cut -d ' ' -f 1)
if [ "$os" = "Linux" ]; then
    echo "This is a Linux system, we can install compilers with conda"
    conda install gxx_linux-64 -y
else
    echo "please check and install compilers g++ >8 and mpich manually"
fi


# define your target directories
# to override: just call the script like e.g. "CC=whatever bash build.sh"
CC=${CC-mpicc}
MPI_CC=${MPI_CC-mpicc}
CXX=${CXX-mpicxx}
MPI_CXX=${MPI_CXX-mpicxx}
MAD_SRC_DIR=${MAD_SRC_DIR-madness_src}
MAD_ROOT_DIR=${MAD_ROOT_DIR-madness}
NUMCPP_SRC_DIR=${NUMCPP_SRC_DIR-numcpp}

mkdir ${MAD_ROOT_DIR}

export MAD_SRC_DIR=${MAD_SRC_DIR} # add path to directory where you want the source code
export MAD_ROOT_DIR=${MAD_ROOT_DIR} # add path to directory where you want the compiled code
export NUMCPP_SRC_DIR=${NUMCPP_SRC_DIR} # add path to directory where you want the numcpp dependency

# make sure that tequila will find it later
echo "added by madness install script" >> ~/.bashrc
echo "export MAD_ROOT_DIR=${MAD_ROOT_DIR}" >> ~/.bashrc

# get the sources
git clone https://github.com/kottmanj/madness.git $MAD_SRC_DIR
git clone https://github.com/dpilger26/numcpp $NUMCPP_SRC_DIR

# install dependencies
conda install cmake mkl mpich boost -y

# export paths to dependencies
export CPLUS_INCLUDE_PATH=$(realpath $NUMCPP_SRC_DIR/include):$CPLUS_INCLUDE_PATH
# make sure that the right boost and mkl are found and used
export CPLUS_INCLUDE_PATH=$(realpath $CONDA_PREFIX/include):$CPLUS_INCLUDE_PATH

cmake -D ENABLE_MKL=ON -D CMAKE_CXX_FLAGS='-O3 -DNDEBUG -march=native' -S $MAD_SRC_DIR -B $MAD_ROOT_DIR

# compile
make -j -C $MAD_ROOT_DIR

# add MAD_ROOT_DIR to conda env
conda env config vars set MAD_ROOT_DIR=$(realpath ${MAD_ROOT_DIR})
~                                                        
