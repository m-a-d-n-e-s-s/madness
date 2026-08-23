# MADNESS Docker Containers

This directory contains Docker configurations and build scripts for building MADNESS container images supporting application users (e.g., running `madqc`) and developers linking against the MADNESS numerical libraries.

## Supported Architectures and Linear Algebra Backends

- **ARM64 (`linux/arm64` / `aarch64`)**: Uses single-threaded OpenBLAS (installed via `libopenblas-serial-dev` and linked against `-lopenblas`; compiled and linked with `-O3 -march=armv8-a -mtune=generic` for generic ARMv8-A portability, `-DENABLE_MPI=OFF`, and symbol stripping to minimize container image size).
- **x86_64 (`linux/amd64`)**: Uses the Intel Math Kernel Library (MKL).

## Directory Structure

```
madness/admin/docker/
├── README.md               # Documentation and usage instructions
├── images/
│   └── Makefile            # Orchestrates building, tagging, and pushing Docker images
└── ubuntu/
    ├── Dockerfile          # Multi-architecture Dockerfile (auto-detects ARM64 / x86_64)
    ├── Dockerfile.arm64    # Dedicated ARM64 Dockerfile using OpenBLAS
    ├── Dockerfile.x86_64   # Dedicated x86_64 Dockerfile using Intel MKL
    ├── Makefile            # Build targets for packages, linear algebra, and MADNESS compilation
    └── sudoers             # Sudoers configuration for container users
```

## How to Build

### Building with Make (from `images/` directory)

Build the ARM64 image (defaults to 2 parallel jobs for Raspberry Pi 5 / memory-constrained hosts):
```bash
cd madness/admin/docker/images
make arm64
```

Build the x86_64 image:
```bash
cd madness/admin/docker/images
make x86_64
```

Build with custom compilation options (e.g., Debug / `-O0` optimization and concurrency limit):
```bash
cd madness/admin/docker/images
make arm64 BUILD_TYPE=Debug CMAKE_EXTRA_ARGS="-DCMAKE_CXX_FLAGS=-O0" BUILD_PARALLEL_JOBS=2
```

### Building Directly with Docker

For ARM64:
```bash
docker build \
    --build-arg BUILD_TYPE=Release \
    --build-arg BUILD_PARALLEL_JOBS=2 \
    -f madness/admin/docker/ubuntu/Dockerfile.arm64 \
    -t rjharrison/ubuntu-arm64:22.04 \
    madness/admin/docker/ubuntu/
```

For x86_64:
```bash
docker build \
    --build-arg BUILD_TYPE=Release \
    -f madness/admin/docker/ubuntu/Dockerfile.x86_64 \
    -t rjharrison/ubuntu-x86_64:22.04 \
    madness/admin/docker/ubuntu/
```

### Pushing to Docker Hub

```bash
docker login -u <docker.com username>
cd madness/admin/docker/images
make push/arm64/22.04
make push/x86_64/22.04
# or push all:
make push/all
```

## User Guide: Running the `madqc` Application

The container includes `madqc` and MADNESS application executables pre-installed under `/home/m-a-d-n-e-s-s/install/bin` which is automatically included in `$PATH`.

### Running `madqc` on an Input File

Mount your current working directory containing the input file (e.g. `scf_he_hf.in`):
```bash
docker run --rm -it \
    -v $(pwd):/work \
    -w /work \
    rjharrison/ubuntu-arm64:22.04 \
    madqc scf_he_hf.in
```

### Interactive Shell

```bash
docker run --rm -it \
    -v $(pwd):/work \
    -w /work \
    rjharrison/ubuntu-arm64:22.04 \
    bash
```

## Developer Guide: Linking Against MADNESS Numerical Libraries

The container provides all development headers, static/shared libraries, and CMake configuration packages:
- **Headers**: `/home/m-a-d-n-e-s-s/install/include`
- **Libraries**: `/home/m-a-d-n-e-s-s/install/lib`
- **CMake Config**: `/home/m-a-d-n-e-s-s/install/lib/cmake/madness`
- **Preconfigured Environment**: `CMAKE_PREFIX_PATH`, `PATH`, and `LD_LIBRARY_PATH` are pre-set.

### Example: Building a C++ Program with CMake

Create a `main.cpp`:
```cpp
#include <madness/world/world.h>
#include <iostream>

int main(int argc, char** argv) {
    madness::initialize(argc, argv);
    madness::World& world = madness::World::find_instance();
    std::cout << "Hello from MADNESS World! Rank: " << world.rank()
              << " Size: " << world.size() << std::endl;
    madness::finalize();
    return 0;
}
```

Create a `CMakeLists.txt`:
```cmake
cmake_minimum_required(VERSION 3.12)
project(myapp CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(madness REQUIRED)

add_executable(myapp main.cpp)
target_link_libraries(myapp PRIVATE madness)
```

Compile inside the container:
```bash
docker run --rm -it \
    -v $(pwd):/project \
    -w /project \
    rjharrison/ubuntu-arm64:22.04 \
    bash -c "cmake -B build -S . && cmake --build build && ./build/myapp"
```
