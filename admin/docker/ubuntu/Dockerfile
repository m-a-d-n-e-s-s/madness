# syntax=docker/dockerfile:1
# ------------------------------------------------------------------------------
# MADNESS Dockerfile (Single-Node / Non-MPI with Sequential Intel MKL)
# Provides 'madqc' and MADNESS development libraries/headers.
# ------------------------------------------------------------------------------

ARG UBUNTU_VERSION=24.04

# ==============================================================================
# Stage 1: Build MADNESS and applications
# ==============================================================================
FROM ubuntu:${UBUNTU_VERSION} AS builder

ENV DEBIAN_FRONTEND=noninteractive

# Install build-time dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    python3 \
    libmkl-dev \
    libxc-dev \
    libeigen3-dev \
    nlohmann-json3-dev \
    ca-certificates \
    git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /usr/src/madness

# Copy source tree
COPY . /usr/src/madness

# Configure, build and install MADNESS without MPI and with sequential MKL
ARG CMAKE_BUILD_TYPE=Release
RUN cmake -S /usr/src/madness -B /tmp/build -GNinja \
    -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE} \
    -DENABLE_MPI=OFF \
    -DENABLE_MKL=ON \
    -DBUILD_TESTING=OFF \
    -DCMAKE_INSTALL_PREFIX=/usr/local && \
    ninja -C /tmp/build madness madqc && \
    ninja -C /tmp/build install-madness-libraries install-madness-config && \
    cmake --install /tmp/build/src/apps/madqc && \
    cmake --install /tmp/build --component chem && \
    strip /usr/local/bin/* || true && \
    strip -g /usr/local/lib/*.a || true && \
    rm -rf /tmp/build

# ==============================================================================
# Stage 2: Final image (Runtime & Development Environment)
# ==============================================================================
FROM ubuntu:${UBUNTU_VERSION}

ENV DEBIAN_FRONTEND=noninteractive

# Install runtime & development dependencies (compiler, cmake, MKL sequential, LibXC, Eigen, JSON)
RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    cmake \
    ninja-build \
    pkg-config \
    libmkl-dev \
    libxc-dev \
    libeigen3-dev \
    nlohmann-json3-dev \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy installed MADNESS binaries, libraries, headers, CMake configs and data files
COPY --from=builder /usr/local /usr/local

# Copy developer example project
COPY admin/docker/examples/hello_world /usr/local/share/madness/examples/hello_world

# Configure environment variables
ENV PATH="/usr/local/bin:${PATH}"
ENV LD_LIBRARY_PATH="/usr/local/lib"
ENV PKG_CONFIG_PATH="/usr/local/lib/pkgconfig"
ENV CMAKE_PREFIX_PATH="/usr/local"
ENV MRA_DATA_DIR="/usr/local/share/madness/0.10.1/data"
ENV MRA_CHEMDATA_DIR="/usr/local/share/madness/0.10.1/data"

# Default working directory for volume mounts
WORKDIR /work

CMD ["/bin/bash"]
