# MADNESS Developer Example: Hello World

This example demonstrates how to build and execute a C++20 numerical application based on the MADNESS library.

The application initializes the parallel runtime, sets up the 1D Multiresolution Analysis (MRA) numerical environment, projects an analytic function $f(x) = \sin(x)$ into an adaptive wavelet basis, and computes its definite integral $\int_0^\pi \sin(x)\,dx = 2.0$.

---

## Building with Make

To compile and link using the provided `Makefile`:

```bash
make
./hello_madness
```

To clean build artifacts:
```bash
make clean
```

---

## Building with CMake

To configure and compile using CMake:

```bash
mkdir -p build && cd build
cmake ..
cmake --build .
./hello_madness
```

---

## Running via Docker from the Host

You can compile and run this example directly using the MADNESS container without entering an interactive shell:

### Using Make:
```bash
docker run --rm -v "$(pwd)":/work -w /work madness:latest bash -c "make && ./hello_madness"
```

### Using CMake:
```bash
docker run --rm -v "$(pwd)":/work -w /work madness:latest bash -c "cmake -B build -S . && cmake --build build && ./build/hello_madness"
```
