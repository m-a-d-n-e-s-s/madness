#
# Generic Toolchain for MPI + MKL + TBB
#
# REQUIREMENTS:
# - in PATH:
#   mpicc and mpicxx
#
# OPTIONS:
# - environment variables:
#   * MKLROOT: the MKL root directory; if not set, will use /opt/intel/mkl
#   * TBBROOT: the TBB root directory; if not set, will use /opt/intel/tbb
#

# Compilers

# These seem to break findmpi module
#if(NOT CMAKE_C_COMPILER)
#  set(CMAKE_C_COMPILER mpicc)
#endif()
#if(NOT CMAKE_CXX_COMPILER)
#  set(CMAKE_CXX_COMPILER mpicxx)
#endif()
#if (NOT MPI_C_COMPILER)
#  set(MPI_C_COMPILER mpicc)
#endif()
#if (NOT MPI_CXX_COMPILER)
#  set(MPI_CXX_COMPILER mpicxx)
#endif()

# Specify the GNU compilers to use with Intel C/C++
set(GCC_ROOT_DIR "$ENV{GCC_DIR}")
set(GCC_SUFFIX "$ENV{GCC_SUFFIX}")
if (CMAKE_CXX_COMPILER_ID EQUAL "Intel" AND GCC_ROOT_DIR)
  set(CMAKE_C_FLAGS_INIT             "-std=c99 -gcc-name=${GCC_ROOT_DIR}/bin/gcc${GCC_SUFFIX}" CACHE STRING "Inital C compile flags")
  set(CMAKE_CXX_FLAGS_INIT           "-gxx-name=${GCC_ROOT_DIR}/bin/g++${GCC_SUFFIX}" CACHE STRING "Initial C++ compile flags")
else()
  set(CMAKE_C_FLAGS_INIT             "-std=c99" CACHE STRING "Inital C compile flags")
  set(CMAKE_CXX_FLAGS_INIT           "" CACHE STRING "Initial C++ compile flags")
endif()

set(CMAKE_C_FLAGS_DEBUG            "-g -Wall -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C debug compile flags")
set(CMAKE_C_FLAGS_MINSIZEREL       "-O1 -DNDEBUG -march=native" CACHE STRING "Inital C minimum size release compile flags")
set(CMAKE_C_FLAGS_RELEASE          "-O3 -DNDEBUG -march=native" CACHE STRING "Inital C release compile flags")
set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -DNDEBUG -g -Wall -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C release with debug info compile flags")
set(CMAKE_CXX_FLAGS_DEBUG          "-g -Wall -Wno-class-memaccess -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C++ debug compile flags")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-O1 -DNDEBUG -march=native -Wall -Wno-class-memaccess -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C++ minimum size release compile flags")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG -march=native -Wall -Wno-class-memaccess -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C++ release compile flags")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -DNDEBUG -g -Wall -Wno-class-memaccess -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C++ release with debug info compile flags")

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
set(LAPACK_LIBRARIES "-L${MKL_ROOT_DIR}/lib" "-lmkl_intel_ilp64" "-lmkl_core"
        "-lmkl_sequential" "-lm" CACHE STRING "LAPACK linker flags")
set(LAPACK_INCLUDE_DIRS ${MKL_ROOT_DIR}/include CACHE STRING "LAPACK include directories")
set(LAPACK_COMPILE_DEFINITIONS MKL_DIRECT_CALL;MKL_INT=long CACHE STRING "LAPACK compile definitions")
set(FORTRAN_INTEGER_SIZE "8" CACHE STRING "Set Fortran integer size to 8 bytes")
