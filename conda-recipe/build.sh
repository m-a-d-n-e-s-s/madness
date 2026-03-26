#!/bin/bash
set -ex

mkdir -p build
cd build

cmake ${CMAKE_ARGS} \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DBUILD_SHARED_LIBS=OFF \
  -DBUILD_TESTING=OFF \
  -DMADNESS_TASK_BACKEND=Pthreads \
  -DMPI_CXX_SKIP_MPICXX=ON \
  -DMPIEXEC_PREFLAGS='--allow-run-as-root' \
  ..

ninja install
