ubuntu_codename ?= $(shell lsb_release -sc 2>/dev/null || echo jammy)
clang_version ?= 13

export DEBIAN_FRONTEND=noninteractive

ubuntu_release ?= $(shell lsb_release -sr 2>/dev/null || echo 22.04)

HOST_ARCH := $(shell uname -m)
ifeq ($(HOST_ARCH),aarch64)
  DEFAULT_ARCH := arm64
else ifeq ($(HOST_ARCH),arm64)
  DEFAULT_ARCH := arm64
else ifeq ($(HOST_ARCH),x86_64)
  DEFAULT_ARCH := x86_64
else
  DEFAULT_ARCH := $(HOST_ARCH)
endif

ARCH ?= $(DEFAULT_ARCH)

BUILD_TYPE ?= Release
BUILD_PARALLEL_JOBS ?= 2
CMAKE_EXTRA_ARGS ?=
MADNESS_REPO ?= https://github.com/m-a-d-n-e-s-s/madness
ARM64_FLAGS ?= -O3 -march=armv8-a -mtune=generic

intel_mkl_version ?= 2020.4-304
intel_mkl := intel-mkl-\*-${intel_mkl_version}

ifneq ($(BUILD_PARALLEL_JOBS),)
  PARALLEL_OPT := -j$(BUILD_PARALLEL_JOBS)
else
  PARALLEL_OPT :=
endif

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
          liblapack-dev liblapacke-dev \
          libunwind-dev libtbb-dev libomp-dev \
          libopenmpi-dev openmpi-bin libscalapack-openmpi-dev \
          libxc-dev \
          python3 python3-pip python3-numpy python3-dev python3-pytest \
          vim emacs-nox ccache

install/openblas:
	apt-get update
	apt-get -yq install libopenblas-serial-dev liblapack-dev liblapacke-dev

install/docker:
	sudo apt-get install -y docker.io

install/g++-%:
	apt-get -yq install gcc-$* g++-$* gfortran-$*

install/add-apt-repository:
	apt-get update
	apt-get -yq install software-properties-common

# CMake
install/cmake:
	if [ "$$(uname -m)" = "x86_64" ]; then \
		(wget -O - https://apt.kitware.com/keys/kitware-archive-latest.asc 2>/dev/null | apt-key add - 2>/dev/null && \
		 apt-add-repository "deb https://apt.kitware.com/ubuntu/ ${ubuntu_codename} main" 2>/dev/null && \
		 apt-get -yq update) || true; \
	fi
	apt-get -yq install cmake cmake-curses-gui

# LLVM / Clang
install/clang-%:
	(wget https://apt.llvm.org/llvm-snapshot.gpg.key -O - 2>/dev/null | apt-key add - 2>/dev/null && \
	 add-apt-repository "deb http://apt.llvm.org/${ubuntu_codename}/ llvm-toolchain-${ubuntu_codename}-$* main" 2>/dev/null && \
	 apt-get -yq update) || true
	apt-get -yq install clang-$* libomp-$*-dev || apt-get -yq install clang libomp-dev

install/clang: install/clang-${clang_version}

# Intel MKL (x86_64)
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

# MADNESS build targets
build/madness: build/madness/$(ARCH)

strip/madness:
	find /home/m-a-d-n-e-s-s/install/bin -type f -exec strip --strip-all {} + 2>/dev/null || true
	find /home/m-a-d-n-e-s-s/install/lib -type f -name "*.so*" -exec strip --strip-unneeded {} + 2>/dev/null || true
	find /home/m-a-d-n-e-s-s/install/lib -type f -name "*.a" -exec strip --strip-debug {} + 2>/dev/null || true

build/madness/x86_64:
	if [ ! -d /home/m-a-d-n-e-s-s/source ]; then \
		git clone --depth 1 $(MADNESS_REPO) /home/m-a-d-n-e-s-s/source; \
	fi
	MKLROOT=/opt/intel/mkl cmake -S /home/m-a-d-n-e-s-s/source -B /home/m-a-d-n-e-s-s/build \
		-GNinja \
		-DCMAKE_INSTALL_PREFIX=/home/m-a-d-n-e-s-s/install \
		-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
		$(CMAKE_EXTRA_ARGS)
	cmake --build /home/m-a-d-n-e-s-s/build --target applications-madness $(PARALLEL_OPT)
	cmake --build /home/m-a-d-n-e-s-s/build --target install $(PARALLEL_OPT)
	$(MAKE) strip/madness
	MAD_SHARE_DATA=$$(find /home/m-a-d-n-e-s-s/install/share/madness -type d -name "data" | head -n 1); \
	if [ -n "$$MAD_SHARE_DATA" ] && [ "$$MAD_SHARE_DATA" != "/home/m-a-d-n-e-s-s/install/share/madness/data" ]; then \
		ln -sfn "$$MAD_SHARE_DATA" /home/m-a-d-n-e-s-s/install/share/madness/data; \
	fi
	rm -rf /home/m-a-d-n-e-s-s/source /home/m-a-d-n-e-s-s/build

build/madness/amd64: build/madness/x86_64

build/madness/arm64:
	if [ ! -d /home/m-a-d-n-e-s-s/source ]; then \
		git clone --depth 1 $(MADNESS_REPO) /home/m-a-d-n-e-s-s/source; \
	fi
	cmake -S /home/m-a-d-n-e-s-s/source -B /home/m-a-d-n-e-s-s/build \
		-GNinja \
		-DCMAKE_INSTALL_PREFIX=/home/m-a-d-n-e-s-s/install \
		-DCMAKE_BUILD_TYPE=$(BUILD_TYPE) \
		-DCMAKE_C_FLAGS="$(ARM64_FLAGS)" \
		-DCMAKE_CXX_FLAGS="$(ARM64_FLAGS)" \
		-DCMAKE_EXE_LINKER_FLAGS="$(ARM64_FLAGS)" \
		-DCMAKE_SHARED_LINKER_FLAGS="$(ARM64_FLAGS)" \
		-DCMAKE_MODULE_LINKER_FLAGS="$(ARM64_FLAGS)" \
		-DLAPACK_LIBRARIES="openblas" \
		-DENABLE_MPI=OFF \
		-DENABLE_MKL=OFF \
		-DENABLE_ACML=OFF \
		$(CMAKE_EXTRA_ARGS)
	cmake --build /home/m-a-d-n-e-s-s/build --target applications-madness $(PARALLEL_OPT)
	cmake --build /home/m-a-d-n-e-s-s/build --target install $(PARALLEL_OPT)
	$(MAKE) strip/madness
	MAD_SHARE_DATA=$$(find /home/m-a-d-n-e-s-s/install/share/madness -type d -name "data" | head -n 1); \
	if [ -n "$$MAD_SHARE_DATA" ] && [ "$$MAD_SHARE_DATA" != "/home/m-a-d-n-e-s-s/install/share/madness/data" ]; then \
		ln -sfn "$$MAD_SHARE_DATA" /home/m-a-d-n-e-s-s/install/share/madness/data; \
	fi
	rm -rf /home/m-a-d-n-e-s-s/source /home/m-a-d-n-e-s-s/build

build/madness/aarch64: build/madness/arm64
