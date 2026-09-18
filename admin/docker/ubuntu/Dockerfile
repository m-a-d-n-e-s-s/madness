# syntax=docker/dockerfile:1
# ------------------------------------------------------------------------------
# MADNESS Dockerfile (Single-Node / Non-MPI, sequential BLAS)
# Provides 'madqc' and MADNESS development libraries/headers.
#
# BLAS/LAPACK is chosen from the target architecture, because MADNESS owns
# parallelism through its own task pool and therefore requires a *sequential*
# BLAS (see CLAUDE.md):
#   amd64 -> Intel MKL (cmake/modules/FindMKL.cmake requires mkl_sequential)
#   other -> OpenBLAS-serial + LAPACKE (MKL is x86-64 only)
# ------------------------------------------------------------------------------

ARG UBUNTU_VERSION=24.04

# ==============================================================================
# Stage 1: Build MADNESS and applications
# ==============================================================================
FROM ubuntu:${UBUNTU_VERSION} AS builder

ARG DEBIAN_FRONTEND=noninteractive
# TARGETARCH is populated automatically by BuildKit; it is empty with the
# classic builder, in which case we assume the x86-64/MKL configuration.
ARG TARGETARCH

# Install build-time dependencies
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
      build-essential \
      cmake \
      ninja-build \
      python3 \
      libxc-dev \
      libeigen3-dev \
      nlohmann-json3-dev \
      ca-certificates \
      git; \
    if [ "${TARGETARCH:-amd64}" = "amd64" ]; then \
      apt-get install -y --no-install-recommends libmkl-dev; \
    else \
      apt-get install -y --no-install-recommends libopenblas-serial-dev liblapacke-dev; \
    fi; \
    rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/madness

# Copy source tree
COPY . /usr/src/madness

# Configure, build and install MADNESS without MPI and with a sequential BLAS.
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
#   * likewise assert MKL was actually found on amd64: if FindMKL fails,
#     external/lapack.cmake quietly falls back to a generic (possibly threaded)
#     BLAS, which is exactly what this image must not ship.
#   * install-madness-libraries covers every library component: headers, the
#     aggregate libmadness.a, and the mra/chem data files.
#   * The data directory is versioned (share/madness/<version>/data), so it is
#     symlinked to a stable path and the runtime stage does not have to hardcode
#     MADNESS_VERSION. `coeffs` comes from the mra component and `sto-3g` from
#     chem, so testing for both also proves both components installed.
ARG CMAKE_BUILD_TYPE=Release
RUN set -eux; \
    if [ "${TARGETARCH:-amd64}" = "amd64" ]; then blas_args="-DENABLE_MKL=ON"; \
    else blas_args="-DENABLE_MKL=OFF"; fi; \
    cmake -S /usr/src/madness -B /tmp/build -GNinja \
      -DCMAKE_BUILD_TYPE="${CMAKE_BUILD_TYPE}" \
      -DENABLE_MPI=OFF \
      -DENABLE_LIBXC=ON \
      ${blas_args} \
      -DBUILD_TESTING=OFF \
      -DCMAKE_INSTALL_PREFIX=/usr/local; \
    grep -q '^#define MADNESS_HAS_LIBXC 1' /tmp/build/src/madness/config.h; \
    if [ "${TARGETARCH:-amd64}" = "amd64" ]; then \
      grep -q '^#define HAVE_INTEL_MKL 1' /tmp/build/src/madness/config.h; \
    fi; \
    ninja -C /tmp/build madness madqc; \
    ninja -C /tmp/build install-madness-libraries install-madness-config; \
    cmake --install /tmp/build/src/apps/madqc; \
    test -x /usr/local/bin/madqc; \
    datadir="$(dirname "$(find /usr/local/share/madness -name coeffs -print -quit)")"; \
    test -f "${datadir}/coeffs"; \
    test -f "${datadir}/sto-3g"; \
    ln -sfn "${datadir}" /usr/local/share/madness/data; \
    strip /usr/local/bin/* || true; \
    strip -g /usr/local/lib/*.a || true; \
    rm -rf /tmp/build

# ==============================================================================
# Stage 2: Final image (Runtime & Development Environment)
# ==============================================================================
FROM ubuntu:${UBUNTU_VERSION}

ARG DEBIAN_FRONTEND=noninteractive
ARG TARGETARCH

# Install runtime & development dependencies (compiler, cmake, sequential BLAS,
# LibXC, Eigen, JSON). The BLAS, libxc and nlohmann_json packages must match the
# builder's: the installed madness-config.cmake points consumers back at the
# very same system packages.
RUN set -eux; \
    apt-get update; \
    apt-get install -y --no-install-recommends \
      build-essential \
      cmake \
      ninja-build \
      pkg-config \
      libxc-dev \
      libeigen3-dev \
      nlohmann-json3-dev \
      ca-certificates; \
    if [ "${TARGETARCH:-amd64}" = "amd64" ]; then \
      apt-get install -y --no-install-recommends libmkl-dev; \
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
ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib"
ENV CMAKE_PREFIX_PATH="/usr/local"
ENV MRA_DATA_DIR="/usr/local/share/madness/data"
ENV MRA_CHEMDATA_DIR="/usr/local/share/madness/data"

# Default working directory for volume mounts
WORKDIR /work

CMD ["/bin/bash"]
