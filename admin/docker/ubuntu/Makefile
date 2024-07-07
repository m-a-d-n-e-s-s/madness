ubuntu_codename ?= $(shell lsb_release -sc)
clang_version ?= 13

export DEBIAN_FRONTEND=noninteractive

ubuntu_release = $(shell lsb_release -sr)

intel_mkl_version ?= 2020.4-304
intel_mkl := intel-mkl-\*-${intel_mkl_version}

install/main:
	apt-get update
	apt-get -yq install \
          lsb-release coreutils sudo bash-completion \
          apt-transport-https software-properties-common ca-certificates gnupg \
          linux-tools-common time pciutils \
          build-essential wget curl \
          git make ninja-build \
          gcc g++ gfortran gdb valgrind \
          libeigen3-dev \
          libblas-dev liblapack-dev liblapacke-dev \
          libunwind-dev libtbb-dev libomp-dev \
          libopenmpi-dev openmpi-bin libscalapack-openmpi-dev \
          libxc-dev \
          python3 python3-pip python3-numpy python3-dev python3-pytest \
          vim emacs-nox ccache

install/docker:
	sudo apt-get install docker.io

install/g++-%:
	apt-get -yq install gcc-$* g++-$* gfortran-$*

install/add-apt-repository:
	apt-get update
	apt-get -yq install software-properties-common

# CMake
install/cmake:
	wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc | apt-key add -
	apt-add-repository "deb https://apt.kitware.com/ubuntu/ ${ubuntu_codename} main"
	apt-get -yq update
	apt-get -yq install cmake cmake-curses-gui


# LLVM
install/clang-%:
	wget https://apt.llvm.org/llvm-snapshot.gpg.key -O - | apt-key add -
	add-apt-repository "deb http://apt.llvm.org/${ubuntu_codename}/ llvm-toolchain-${ubuntu_codename}-$* main"
	apt-get -yq update
	apt-get -yq install clang-$* libomp-$*-dev

install/clang: install/clang-${clang_version}


# Intel
apt-add-repository/intel-mkl:
	wget -O - https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB | apt-key add -
	apt-add-repository "deb https://apt.repos.intel.com/mkl all main"
	apt-get -yq update

download/intel-mkl:
	apt-get install -yq --download-only ${intel_mkl}

install/intel-mkl:
	apt-get -yq install ${intel_mkl}
	rm -rf /opt/intel/mkl
	ln -s /opt/intel/compilers_and_libraries_$(subst -,.,${intel_mkl_version})/linux/mkl /opt/intel/mkl
	test -x /opt/intel/mkl/bin/mklvars.sh

# MADNESS
build/madness:
	git clone https://github.com/m-a-d-n-e-s-s/madness /home/m-a-d-n-e-s-s/source
	MKLROOT=/opt/intel/mkl cmake -S /home/m-a-d-n-e-s-s/source -B /home/m-a-d-n-e-s-s/build -GNinja -DCMAKE_INSTALL_PREFIX=/home/m-a-d-n-e-s-s/install
	cmake --build /home/m-a-d-n-e-s-s/build --target applications-madness
	cmake --build /home/m-a-d-n-e-s-s/build --target install

