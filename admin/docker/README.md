# MADNESS Container Image Guide

This directory provides the Docker container configuration for **MADNESS** (Multiresolution Adaptive Numerical Environment for Scientific Simulation) and the **`madqc`** quantum chemistry application.

The container is configured for single-node / non-MPI execution with stubbed-out MPI and a **sequential** (single-threaded) BLAS/LAPACK — MADNESS owns parallelism through its own task pool, so a threaded BLAS would oversubscribe the cores. It provides both:
1. Ready-to-use binaries (including **`madqc`** and chemistry datasets).
2. A complete C++20 development environment (headers, static libraries and CMake configuration files) for building and running new MADNESS-based applications.

Consume the installation with `find_package(madness CONFIG)`. MADNESS has no
working `pkg-config` module — `config/MADNESS.pc.in` is an unmaintained
autotools leftover whose substitutions come out empty — so the image
deliberately does not ship or advertise one.

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

### Architecture & the BLAS choice
Intel MKL is x86-64 only, so the BLAS/LAPACK dependency is selected from
BuildKit's `TARGETARCH`:

| Target | Packages | CMake |
| --- | --- | --- |
| `linux/amd64` | `libmkl-dev` | `-DENABLE_MKL=ON` (`FindMKL.cmake` requires `mkl_sequential`) |
| anything else | `libopenblas-serial-dev`, `liblapacke-dev` | `-DENABLE_MKL=OFF` |

`linux/amd64` is the configuration this image is built and used in. The
non-x86-64 path is provided so that `docker buildx build --platform linux/arm64`
has a sequential BLAS to fall back on, but it is **not** covered by CI — treat it
as best-effort and expect to adjust package names per Ubuntu release.

```bash
docker buildx build --platform linux/amd64 -t madness:amd64 -f admin/docker/ubuntu/Dockerfile . --load
```

With the classic (non-BuildKit) builder `TARGETARCH` is empty and the x86-64/MKL
configuration is assumed.

### Publishing images
`admin/docker/images/Makefile` builds and pushes the tagged images. It invokes
the Dockerfile with the repository root as the build context, which is what the
`COPY . /usr/src/madness` step requires:

```bash
docker login -u <docker.com username>
cd admin/docker/images
make all         # build an image per $(versions), tag the newest as :latest
make push/all    # build, tag and push everything
```

Override `repo`, `versions` or `latest` on the command line to build elsewhere,
e.g. `make repo=myorg versions="24.04 22.04" all`.

### CI coverage
`.github/workflows/docker.yml` builds the image and exercises it: it reproduces
the `src/examples/qc/scf_he_hf` reference *through the container* (by pointing
the qctest harness's `MADQC` at a `docker run` wrapper), runs an LDA deck to
prove libxc is live at runtime, and builds the developer example against the
install tree with CMake. It is path-filtered to `admin/docker/**` plus a weekly
schedule, since a full in-container MADNESS build is too expensive for every
push. The header of that file records what it deliberately does not cover.

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

Option C (CMake) is the supported route — it picks up MADNESS's public compile
definitions and BLAS include paths automatically. The hand-written form below
hardcodes the amd64/MKL link line and must be adjusted for a non-x86-64 image
(`-lopenblas -llapacke` in place of the `mkl_*` libraries).

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

The example lives at `admin/docker/examples/hello_world/hello.cc` in the source
tree and at `/usr/local/share/madness/examples/hello_world/hello.cc` in the
image. Read it there rather than from a copy in this file — a duplicated listing
drifts, and this one previously drifted into misusing the runtime.

It projects `f(x) = sin(x)` onto an adaptive wavelet basis on `[0, pi]` and
computes `int_0^pi sin(x) dx = 2.0`. Two runtime invariants it exists to
demonstrate (see the "Runtime" section of `CLAUDE.md`):

- `World& world = initialize(argc, argv);` — `initialize()` already constructs
  the default `World` on `MPI_COMM_WORLD` and returns a reference to it.
  Constructing a second `World` on the same communicator makes `finalize()`
  report `MADNESS runtime finalized but 1 world still exists`, and leaves that
  `World`'s destructor running after MPI has been torn down.
- Every `Function` is scoped so that it is destroyed *before* `finalize()`.

To build and run it with the bundled files:

```bash
cd admin/docker/examples/hello_world
docker run --rm -v "$(pwd)":/work -w /work --user "$(id -u):$(id -g)" \
  madness:latest bash -c "cmake -B build -S . && cmake --build build && ./build/hello_madness"
```

Expected output ends with:

```
Computed integral of sin(x) on [0, pi]: 2.000000e+00
Difference from exact                 : 2.220446e-16
Calculation completed successfully.
```
