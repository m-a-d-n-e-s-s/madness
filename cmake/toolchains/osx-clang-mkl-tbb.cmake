#
# Generic Toolchain for OS X + MPI + MKL + TBB
#
# REQUIREMENTS:
# - in PATH:
#   mpicc and mpicxx, both using clang as the backend
#
# OPTIONS:
# - environment variables:
#   * MKLROOT: the MKL root directory; if not set, will use /opt/intel/mkl
#   * TBBROOT: the TBB root directory; if not set, will use /opt/intel/tbb
#

# this is key to avoid RPATH problems
set(CMAKE_SYSTEM_NAME Darwin)

# Compilers
set(CMAKE_C_COMPILER clang)
set(CMAKE_CXX_COMPILER clang++)
set(MPI_C_COMPILER mpicc)
set(MPI_CXX_COMPILER mpicxx)

# Compile flags
set(CMAKE_C_FLAGS_INIT             "-std=c99" CACHE STRING "Inital C compile flags")
set(CMAKE_C_FLAGS_DEBUG            "-g -Wall" CACHE STRING "Inital C debug compile flags")
set(CMAKE_C_FLAGS_MINSIZEREL       "-Os -march=native -DNDEBUG" CACHE STRING "Inital C minimum size release compile flags")
set(CMAKE_C_FLAGS_RELEASE          "-O3 -march=native -DNDEBUG" CACHE STRING "Inital C release compile flags")
set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g -Wall" CACHE STRING "Inital C release with debug info compile flags")
set(CMAKE_CXX_FLAGS_INIT           "-I/usr/local/include -stdlib=libc++" CACHE STRING "Inital C++ compile flags")
set(CMAKE_CXX_FLAGS_DEBUG          "-g -Wall" CACHE STRING "Inital C++ debug compile flags")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -march=native -DNDEBUG" CACHE STRING "Inital C++ minimum size release compile flags")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG" CACHE STRING "Inital C++ release compile flags")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -Wall" CACHE STRING "Inital C++ release with debug info compile flags")

# Libraries
if(EXISTS $ENV{MKLROOT})
  set(MKL_ROOT_DIR "$ENV{MKLROOT}" CACHE PATH "MKL root directory")
else()
  set(MKL_ROOT_DIR "/opt/intel/mkl" CACHE PATH "MKL root directory")
endif()
if(EXISTS $ENV{TBBROOT})
  set(TBB_ROOT_DIR "$ENV{TBBROOT}" CACHE PATH "TBB root directory")
else()
  set(TBB_ROOT_DIR "/opt/intel/tbb" CACHE PATH "TBB root directory")
endif()

# Flags
set(LAPACK_LIBRARIES "-L${MKL_ROOT_DIR}/lib" "-lmkl_intel_lp64" "-lmkl_core"
    "-lmkl_sequential" "-lm" CACHE STRING "LAPACK linker flags")
set(LAPACK_INCLUDE_DIRS ${MKL_ROOT_DIR}/include CACHE STRING "LAPACK include directories")
set(LAPACK_COMPILE_DEFINITIONS MKL_DIRECT_CALL CACHE STRING "LAPACK compile definitions")
set(INTEGER4 TRUE CACHE BOOL "Set Fortran integer size to 4 bytes")