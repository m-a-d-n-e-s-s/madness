# MADNESS Container Image Guide

This directory provides the Docker container configuration for **MADNESS** (Multiresolution Adaptive Numerical Environment for Scientific Simulation) and the **`madqc`** quantum chemistry application.

The container is configured for single-node / non-MPI execution with stubbed-out MPI and uses **Intel MKL** for single-threaded (sequential) BLAS and LAPACK. It provides both:
1. Ready-to-use binaries (including **`madqc`** and chemistry datasets).
2. A complete C++20 development environment (headers, static libraries, CMake configuration files, and `pkg-config`) for building and running new MADNESS-based applications.

---

## 1. Building the Docker Image

From the root of the MADNESS source tree (or from `madrel/madness`):

```bash
docker build -t madness:latest -f admin/docker/ubuntu/Dockerfile .
```

### Build Options
- **`UBUNTU_VERSION`**: Base Ubuntu release (default: `24.04`).
- **`CMAKE_BUILD_TYPE`**: CMake build configuration (default: `Release`, or `Debug` for fast non-optimized test builds).

Example building with custom arguments:
```bash
docker build --build-arg CMAKE_BUILD_TYPE=Release -t madness:latest -f admin/docker/ubuntu/Dockerfile .
```

### Multi-Architecture & ARM64 Notes
The Dockerfile is structured to support multi-platform builds (e.g. `linux/amd64` and `linux/arm64`). To build for multiple architectures with `docker buildx`:
```bash
docker buildx build --platform linux/amd64 -t madness:amd64 -f admin/docker/ubuntu/Dockerfile . --load
```

---

## 2. Instructions for `madqc` Users

When running quantum chemistry calculations with `madqc`, you want all generated output files (`.out` files, checkpoint JSON files, moldft/response archives, cube files, etc.) stored directly in a directory on the host filesystem.

### Basic Invocations

Mount your current host directory containing your input file into the container's `/work` directory:

```bash
docker run --rm \
  -v "$(pwd)":/work \
  -w /work \
  --user "$(id -u):$(id -g)" \
  madness:latest \
  madqc input.in
```

### Explanation of Flags:
- `-v "$(pwd)":/work`: Mounts the current host working directory into the container at `/work`.
- `-w /work`: Sets the container's working directory to `/work` so relative paths and output files resolve to the mounted directory.
- `--user "$(id -u):$(id -g)"`: Runs the container process with your host user and group ID so generated output files are owned by your user account rather than `root`.
- `--rm`: Automatically removes the container instance after execution completes.

### Controlling Threading
MADNESS uses its internal task scheduler with a thread pool. By default, it detects and utilizes available hardware threads. You can restrict CPU cores using Docker's `--cpus` flag or set the `MAD_NUM_THREADS` environment variable:

```bash
# Limit to 4 CPU threads
docker run --rm --cpus=4 -v "$(pwd)":/work -w /work --user "$(id -u):$(id -g)" madness:latest madqc input.in
```

### Recommended Shell Wrapper / Alias

To use `madqc` seamlessly as if it were installed natively on your host system, add the following shell function to your `~/.bashrc` or `~/.zshrc`:

```bash
madqc() {
    docker run --rm -i \
      -v "$(pwd)":/work \
      -w /work \
      --user "$(id -u):$(id -g)" \
      madness:latest \
      madqc "$@"
}
```

Then you can simply run:
```bash
madqc my_calculation.in
```

---

## 3. Instructions for Application Developers

The container provides the full MADNESS header tree in `/usr/local/include`, static libraries in `/usr/local/lib`, and CMake config files in `/usr/local/lib/cmake/madness`.

An example "Hello World" application is included in the container at `/usr/local/share/madness/examples/hello_world` and in the source repository at `admin/docker/examples/hello_world`.

### Option A: Interactive Development Shell

Launch an interactive shell inside the container with your current project directory mounted:

```bash
docker run -it --rm \
  -v "$(pwd)":/work \
  -w /work \
  --user "$(id -u):$(id -g)" \
  madness:latest /bin/bash
```

Inside the container, you can use `g++`, `make`, and `cmake` directly.

---

### Option B: Building via Make

Create a `Makefile` for your application:

```makefile
CXX ?= g++
CXXFLAGS ?= -std=c++20 -O2
CPPFLAGS ?= -I/usr/local/include
LDFLAGS ?= -L/usr/local/lib
LIBS ?= -lmadness -lxc -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl

TARGET = my_app
SRCS = main.cc

all: $(TARGET)

$(TARGET): $(SRCS)
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $< $(LDFLAGS) $(LIBS) -o $@

clean:
	rm -f $(TARGET) *.o
```

Build and run from host:
```bash
docker run --rm -v "$(pwd)":/work -w /work --user "$(id -u):$(id -g)" madness:latest make
docker run --rm -v "$(pwd)":/work -w /work --user "$(id -u):$(id -g)" madness:latest ./my_app
```

---

### Option C: Building via CMake

Create a `CMakeLists.txt` for your application:

```cmake
cmake_minimum_required(VERSION 3.20)
project(MyApp CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find the installed MADNESS package
find_package(madness REQUIRED CONFIG)

add_executable(my_app main.cc)
target_link_libraries(my_app PRIVATE madness)
```

Configure, build, and run from host:
```bash
docker run --rm -v "$(pwd)":/work -w /work --user "$(id -u):$(id -g)" madness:latest \
  bash -c "cmake -B build -S . && cmake --build build && ./build/my_app"
```

---

## 4. Hello World Developer Example Walkthrough

A complete example source file (`hello.cc`):

```cpp
#include <madness/mra/mra.h>
#include <iostream>
#include <cmath>

using namespace madness;

double my_func(const coord_1d& r) {
    return std::sin(r[0]);
}

int main(int argc, char** argv) {
    // 1. Initialize MADNESS runtime & World communicator
    initialize(argc, argv);
    World world(SafeMPI::COMM_WORLD);

    // 2. Initialize numerical environment
    startup(world, argc, argv);

    if (world.rank() == 0) {
        print("Hello from MADNESS on rank", world.rank(), "of size", world.size());
    }

    // 3. Set 1D computation domain: [0, pi]
    FunctionDefaults<1>::set_cubic_cell(0.0, M_PI);

    // 4. Project f(x) = sin(x) into adaptive wavelet basis
    real_function_1d f = real_factory_1d(world).f(my_func);

    // 5. Compute integral: int_0^pi sin(x) dx = 2.0
    double integral = f.trace();

    if (world.rank() == 0) {
        print("Computed integral of sin(x) on [0, pi]:", integral);
        print("Expected exact value                  :", 2.0);
    }

    // 6. Finalize runtime
    finalize();
    return 0;
}
```

To run this example using the bundled files:
```bash
cd admin/docker/examples/hello_world
docker run --rm -v "$(pwd)":/work -w /work madness:latest make run
```
