#!/bin/bash
set -ex

MPI_FLAG="-DENABLE_MPI=OFF"
if [ "${MPI_VARIANT}" = "openmpi" ]; then
  MPI_FLAG="-DENABLE_MPI=ON"
  # Point OpenMPI wrappers at the host-env libraries/headers
  export OPAL_PREFIX="${PREFIX}"
  # Use the build-env wrappers (they have the compiler toolchain)
  MPI_FLAG="${MPI_FLAG} -DMPI_C_COMPILER=${BUILD_PREFIX}/bin/mpicc"
  MPI_FLAG="${MPI_FLAG} -DMPI_CXX_COMPILER=${BUILD_PREFIX}/bin/mpicxx"
fi

mkdir -p build
cd build

# The static MADNESS libraries get linked into the _pymadness shared
# extension module, so everything must be compiled position-independent.
cmake ${CMAKE_ARGS} \
  -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${PREFIX}" \
  -DCMAKE_PREFIX_PATH="${PREFIX}" \
  -DBUILD_SHARED_LIBS=OFF \
  -DCMAKE_POSITION_INDEPENDENT_CODE=ON \
  -DBUILD_TESTING=OFF \
  ${MPI_FLAG} \
  -DENABLE_LIBXC=OFF \
  -DMADNESS_TASK_BACKEND=Pthreads \
  -DMPI_CXX_SKIP_MPICXX=ON \
  -DMPIEXEC_PREFLAGS='--allow-run-as-root' \
  -DENABLE_PYTHON=ON \
  -DPython3_EXECUTABLE="${PYTHON}" \
  -DPython3_ROOT_DIR="${PREFIX}" \
  ..

ninja install

# --- Relocate pymadness into site-packages -------------------------------
# CMake installs the Python package and the _pymadness extension into
# ${PREFIX}/lib/pymadness, which is not on Python's import path.  The
# package imports the extension as a *top-level* module (`import _pymadness`),
# so the layout must be:
#   ${SP_DIR}/pymadness/{__init__,plotting}.py   <- the package
#   ${SP_DIR}/_pymadness*.so                      <- sibling extension module
PYM_SRC="${PREFIX}/lib/pymadness"
if [ -d "${PYM_SRC}" ]; then
  : "${SP_DIR:=$("${PYTHON}" -c 'import site; print(site.getsitepackages()[0])')}"
  mkdir -p "${SP_DIR}/pymadness"
  # pure-Python package files
  find "${PYM_SRC}" -maxdepth 1 -name '*.py' -exec mv {} "${SP_DIR}/pymadness/" \;
  # compiled extension module as a sibling of the package
  mv "${PYM_SRC}"/_pymadness*.so "${SP_DIR}/"
  rm -rf "${PYM_SRC}"
fi

# Set up environment activation for MADNESS data directories
mkdir -p "${PREFIX}/etc/conda/activate.d"
mkdir -p "${PREFIX}/etc/conda/deactivate.d"

cat > "${PREFIX}/etc/conda/activate.d/activate_madness.sh" <<'EOF'
# MADNESS installs data into a versioned subdirectory; find the first one.
# Only one version of MADNESS is expected per conda environment.
MADNESS_DATA_BASE="${CONDA_PREFIX}/share/madness"
if [ -d "${MADNESS_DATA_BASE}" ]; then
  for d in "${MADNESS_DATA_BASE}"/*; do
    if [ -d "$d/data" ]; then
      export MRA_DATA_DIR="$d/data"
      export MRA_CHEMDATA_DIR="$d/data"
      break
    fi
  done
fi
EOF

cat > "${PREFIX}/etc/conda/deactivate.d/deactivate_madness.sh" <<'EOF'
unset MRA_DATA_DIR
unset MRA_CHEMDATA_DIR
EOF
