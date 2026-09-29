# syntax=docker/dockerfile:1
# ------------------------------------------------------------------------------
# MADNESS Dockerfile (MPI with TBB task backend, Intel MKL / ARMPL, TCMalloc-minimal)
# Provides 'madqc' and MADNESS development libraries/headers.
#
# BLAS/LAPACK is chosen from the target architecture, because MADNESS owns
# parallelism through its own task pool and therefore requires a *sequential*
# BLAS (see CLAUDE.md):
#   amd64 -> Intel MKL (cmake/modules/FindMKL.cmake requires mkl_sequential)
#   arm64 -> Arm Performance Libraries (ARMPL, sequential)
#   other -> OpenBLAS-serial + LAPACKE (MKL is x86-64 only)
# ------------------------------------------------------------------------------

ARG UBUNTU_VERSION=24.04

# ==============================================================================
# Stage 1: Build MADNESS and applications
# ==============================================================================
FROM ubuntu:${UBUNTU_VERSION} AS builder

ARG DEBIAN_FRONTEND=noninteractive
# TARGETARCH is populated automatically by BuildKit; if empty (e.g. classic builder),
# we fall back to dpkg --print-architecture so both amd64 and arm64 are detected.
ARG TARGETARCH

# Install build-time dependencies
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
      build-essential \
      cmake \
      ninja-build \
      python3 \
      libgoogle-perftools-dev \
      libxc-dev \
      libfftw3-dev \
      libopenmpi-dev \
      openmpi-bin \
      libtbb-dev \
      libeigen3-dev \
      nlohmann-json3-dev \
      ca-certificates \
      curl \
      git; \
    arch="${TARGETARCH:-$(dpkg --print-architecture)}"; \
    if [ "$arch" = "amd64" ]; then \
      apt-get install -y --no-install-recommends libmkl-dev; \
      find /usr/lib -name "libmkl_*.a" -delete; \
      find /usr/lib -name "*mkl_*ilp64*.so*" -delete; \
      find /usr/lib -name "*mkl_*thread*.so*" -delete; \
      find /usr/lib -name "*mkl_blacs*.so*" -delete; \
      find /usr/lib -name "*mkl_scalapack*.so*" -delete; \
      find /usr/lib -name "*mkl_cdft*.so*" -delete; \
    elif [ "$arch" = "arm64" ]; then \
      . /etc/os-release; \
      curl -sL "https://developer.arm.com/packages/arm-toolchains/ubuntu/pool/arm-toolchains-repository_2-2~${VERSION_CODENAME}_all.deb" -o /tmp/arm-repo.deb; \
      dpkg -i /tmp/arm-repo.deb; \
      rm -f /tmp/arm-repo.deb; \
      apt-get update; \
      apt-get install -y --no-install-recommends arm-performance-libraries; \
      rm -rf /opt/arm/arm-performance-libraries/lib/*ilp64* \
             /opt/arm/arm-performance-libraries/lib/*_mp* \
             /opt/arm/arm-performance-libraries/lib/*int64* \
             /opt/arm/arm-performance-libraries/examples_* \
             /opt/arm/arm-performance-libraries/lib/*.a; \
      echo "/opt/arm/arm-performance-libraries/lib" > /etc/ld.so.conf.d/armpl.conf; \
      ldconfig; \
    else \
      apt-get install -y --no-install-recommends libopenblas-serial-dev liblapacke-dev; \
    fi; \
    rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/madness

# Copy source tree
COPY . /usr/src/madness

# Configure, build and install MADNESS with MPI, TBB, tcmalloc-minimal and a sequential BLAS.
#
# `set -eux` plus `;` separators are load-bearing. With `&&` chaining, the
# trailing `|| true` on the strip steps would also swallow a cmake or ninja
# failure -- && and || are equal-precedence and left-associative in sh -- so a
# broken build would produce an image that looks fine but has no madqc in it.
#
# The checks in this step, in order:
#   * libxc is found in CONFIG mode only (external/libxc.cmake), so a missing
#     LibxcConfig.cmake would silently yield a DFT-less quantum-chemistry
#     image. Assert on the generated config.h instead.
#   * assert tcmalloc-minimal and TBB backend are active in config.h.
#   * likewise assert vendor-tuned sequential BLAS was actually found: MKL on amd64,
#     ARMPL on arm64. If detection fails, external/lapack.cmake quietly falls
#     back to a generic (possibly threaded) BLAS, which this image must not ship.
#   * install-madness-libraries covers every library component: headers, the
#     aggregate libmadness.a, and the mra/chem data files.
#   * The data directory is versioned (share/madness/<version>/data), so it is
#     symlinked to a stable path and the runtime stage does not have to hardcode
#     MADNESS_VERSION. `coeffs` comes from the mra component and `sto-3g` from
#     chem, so testing for both also proves both components installed.
ARG CMAKE_BUILD_TYPE=Release
RUN set -eux; \
    arch="${TARGETARCH:-$(dpkg --print-architecture)}"; \
    if [ "$arch" = "amd64" ]; then \
      blas_args="-DENABLE_MKL=ON"; \
    elif [ "$arch" = "arm64" ]; then \
      blas_args="-DENABLE_MKL=OFF -DENABLE_ARMPL=ON"; \
    else \
      blas_args="-DENABLE_MKL=OFF"; \
    fi; \
    cmake -S /usr/src/madness -B /tmp/build -GNinja \
      -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}" \
      -DMADNESS_TARGET_ARCH=performance \
      -DMADNESS_RELAXED_MATH=fast \
      -DENABLE_MPI=ON \
      -DMPI_CXX_SKIP_MPICXX=ON \
      -DMADNESS_TASK_BACKEND=TBB \
      -DENABLE_TCMALLOC_MINIMAL=ON \
      -DENABLE_LIBXC=ON \
      ${blas_args} \
      -DBUILD_TESTING=OFF \
      -DCMAKE_INSTALL_PREFIX=/usr/local; \
    grep -q '^#define MADNESS_HAS_LIBXC 1' /tmp/build/src/madness/config.h; \
    grep -q '^#define MADNESS_HAS_GOOGLE_PERF_TCMALLOC_MINIMAL 1' /tmp/build/src/madness/config.h; \
    grep -q '^#define HAVE_INTEL_TBB 1' /tmp/build/src/madness/config.h; \
    if [ "$arch" = "amd64" ]; then \
      grep -q '^#define HAVE_INTEL_MKL 1' /tmp/build/src/madness/config.h; \
    elif [ "$arch" = "arm64" ]; then \
      grep -q '^#define HAVE_ARMPL 1' /tmp/build/src/madness/config.h; \
    fi; \
    ninja -C /tmp/build madness madqc libxsmm_gemm_generator libxsmm_binaryexport_generator; \
    ninja -C /tmp/build install-madness-libraries install-madness-config; \
    if [ -d /tmp/build/_deps/libxsmm-build ]; then \
      cmake --install /tmp/build/_deps/libxsmm-build; \
    fi; \
    cmake --install /tmp/build/src/apps/madqc; \
    test -x /usr/local/bin/madqc; \
    datadir="$(dirname "$(find /usr/local/share/madness -name coeffs -print -quit)")"; \
    test -f "${datadir}/coeffs"; \
    test -f "${datadir}/sto-3g"; \
    ln -sfn "${datadir}" /usr/local/share/madness/data; \
    if [ "${CMAKE_BUILD_TYPE}" != "Debug" ]; then \
      strip /usr/local/bin/* || true; \
      strip -g /usr/local/lib/*.a || true; \
    fi; \
    for comp in chem linalg misc mra muparser tensor tinyxml world; do \
      ln -sfn libmadness.a /usr/local/lib/libMAD${comp}.a; \
    done; \
    rm -rf /tmp/build

# ==============================================================================
# Stage 2: Final image (Runtime & Development Environment)
# ==============================================================================
FROM ubuntu:${UBUNTU_VERSION}

ARG DEBIAN_FRONTEND=noninteractive
ARG TARGETARCH

# Install runtime & development dependencies (compiler, cmake, sequential BLAS,
# LibXC, FFTW, OpenMPI, TBB, TCMalloc-minimal, Eigen, JSON).
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
      build-essential \
      cmake \
      ninja-build \
      pkg-config \
      libgoogle-perftools-dev \
      libxc-dev \
      libfftw3-dev \
      libopenmpi-dev \
      openmpi-bin \
      libtbb-dev \
      libeigen3-dev \
      nlohmann-json3-dev \
      ca-certificates \
      curl; \
    arch="${TARGETARCH:-$(dpkg --print-architecture)}"; \
    if [ "$arch" = "amd64" ]; then \
      apt-get install -y --no-install-recommends libmkl-dev; \
      find /usr/lib -name "libmkl_*.a" -delete; \
      find /usr/lib -name "*mkl_*ilp64*.so*" -delete; \
      find /usr/lib -name "*mkl_*thread*.so*" -delete; \
      find /usr/lib -name "*mkl_blacs*.so*" -delete; \
      find /usr/lib -name "*mkl_scalapack*.so*" -delete; \
      find /usr/lib -name "*mkl_cdft*.so*" -delete; \
    elif [ "$arch" = "arm64" ]; then \
      . /etc/os-release; \
      curl -sL "https://developer.arm.com/packages/arm-toolchains/ubuntu/pool/arm-toolchains-repository_2-2~${VERSION_CODENAME}_all.deb" -o /tmp/arm-repo.deb; \
      dpkg -i /tmp/arm-repo.deb; \
      rm -f /tmp/arm-repo.deb; \
      apt-get update; \
      apt-get install -y --no-install-recommends arm-performance-libraries; \
      rm -rf /opt/arm/arm-performance-libraries/lib/*ilp64* \
             /opt/arm/arm-performance-libraries/lib/*_mp* \
             /opt/arm/arm-performance-libraries/lib/*int64* \
             /opt/arm/arm-performance-libraries/examples_* \
             /opt/arm/arm-performance-libraries/lib/*.a; \
      echo "/opt/arm/arm-performance-libraries/lib" > /etc/ld.so.conf.d/armpl.conf; \
      ldconfig; \
    else \
      apt-get install -y --no-install-recommends libopenblas-serial-dev liblapacke-dev; \
    fi; \
    rm -rf /var/lib/apt/lists/*

# Copy installed MADNESS binaries, libraries, headers, CMake configs and data files
COPY --from=builder /usr/local /usr/local

# Copy developer example project
COPY admin/docker/examples/hello_world /usr/local/share/madness/examples/hello_world

# Configure environment variables. N.B. no PKG_CONFIG_PATH: MADNESS has no
# working pkg-config module (config/MADNESS.pc.in is an unmaintained autotools
# leftover whose substitutions come out empty), so find_package(madness CONFIG)
# is the only supported way to consume the installation.
ENV PATH="/usr/local/bin:/opt/arm/arm-performance-libraries/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib:/opt/arm/arm-performance-libraries/lib"
ENV CMAKE_PREFIX_PATH="/usr/local:/opt/arm/arm-performance-libraries"
ENV MRA_DATA_DIR="/usr/local/share/madness/data"
ENV MRA_CHEMDATA_DIR="/usr/local/share/madness/data"
ENV OMPI_ALLOW_RUN_AS_ROOT=1
ENV OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1
ENV OMPI_MCA_btl_vader_single_copy_mechanism=none
ENV OMPI_MCA_hwloc_base_binding_policy=none

# Default working directory for volume mounts
WORKDIR /work

CMD ["/bin/bash"]
