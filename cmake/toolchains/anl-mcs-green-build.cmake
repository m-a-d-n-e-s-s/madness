#
# Toolchain for Argonne MCS Green Build Linux: Ubuntu + MPI + MKL
#
# REQUIREMENTS:
# - softenv keys in ~/.soft
#     +mkl
#     +gcc-4.7.3
# - Will need to build CMake 3.5.1 and MPICH on your own
# - Note that MKL root directory is stored in $MKL_HOME 
#
# OPTIONS:
# - environment variables:
#   * MKLROOT: the MKL root directory; if not set, will use /opt/intel/mkl
#   * TBBROOT: the TBB root directory; if not set, will use /opt/intel/tbb
#

# Compilers
set(CMAKE_C_COMPILER gcc)
set(CMAKE_CXX_COMPILER g++)
set(MPI_C_COMPILER mpicc)
set(MPI_CXX_COMPILER mpicxx)

# Compile flags
set(CMAKE_C_FLAGS_INIT             "-std=c99" CACHE STRING "Inital C compile flags")
set(CMAKE_C_FLAGS_DEBUG            "-g -Wall" CACHE STRING "Inital C debug compile flags")
set(CMAKE_C_FLAGS_MINSIZEREL       "-Os -march=native -DNDEBUG" CACHE STRING "Inital C minimum size release compile flags")
set(CMAKE_C_FLAGS_RELEASE          "-O3 -march=native -DNDEBUG" CACHE STRING "Inital C release compile flags")
set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2 -g -Wall" CACHE STRING "Inital C release with debug info compile flags")
set(CMAKE_CXX_FLAGS_INIT           "-std=c++11 " CACHE STRING "Inital C++ compile flags")
set(CMAKE_CXX_FLAGS_DEBUG          "-g -Wall" CACHE STRING "Inital C++ debug compile flags")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -march=native -DNDEBUG" CACHE STRING "Inital C++ minimum size release compile flags")
set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -march=native -DNDEBUG" CACHE STRING "Inital C++ release compile flags")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g -Wall" CACHE STRING "Inital C++ release with debug info compile flags")

# Libraries
if(EXISTS $ENV{MKL_HOME})
  set(MKL_ROOT_DIR "$ENV{MKL_HOME}" CACHE PATH "MKL root directory")
else()
  set(MKL_ROOT_DIR "/opt/intel/mkl" CACHE PATH "MKL root directory")
endif()
if(EXISTS $ENV{TBBROOT})
  set(TBB_ROOT_DIR "$ENV{TBBROOT}" CACHE PATH "TBB root directory")
else()
  set(TBB_ROOT_DIR "/opt/intel/tbb" CACHE PATH "TBB root directory")
endif()

# Flags
set(LAPACK_LIBRARIES "-L${MKL_ROOT_DIR}/lib/intel64 -lmkl_intel_lp64 -lmkl_core -lmkl_sequential" "-lm" CACHE STRING "LAPACK linker flags")
set(INTEGER4 TRUE CACHE BOOL "Set Fortran integer size to 4 bytes")
