#
# ARM Toolchain for MPI + TBB
#
# REQUIREMENTS:
# - in PATH:
#   mpicc and mpicxx
#
# OPTIONS:
# - environment variables:
#

# Compilers

if(NOT CMAKE_C_COMPILER)
  set(CMAKE_C_COMPILER mpicc)
endif()
if(NOT CMAKE_CXX_COMPILER)
  set(CMAKE_CXX_COMPILER mpicxx)
endif()
if (NOT MPI_C_COMPILER)
  set(MPI_C_COMPILER mpicc)
endif()
if (NOT MPI_CXX_COMPILER)
  set(MPI_CXX_COMPILER mpicxx)
endif()

set(CMAKE_C_FLAGS_DEBUG            "-g -Wall -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C debug compile flags")
set(CMAKE_C_FLAGS_MINSIZEREL       "-O1 -DNDEBUG -march=native" CACHE STRING "Inital C minimum size release compile flags")
set(CMAKE_C_FLAGS_RELEASE          "-Ofast  -floop-optimize -falign-loops -falign-labels -falign-functions -falign-jumps -fomit-frame-pointer -DNDEBUG -march=native" CACHE STRING "Inital C release compile flags")
set(CMAKE_C_FLAGS_RELWITHDEBINFO   "-O2  -floop-optimize -falign-loops -falign-labels -falign-functions -falign-jumps -DNDEBUG -g -Wall -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C release with debug info compile flags")
set(CMAKE_CXX_FLAGS_DEBUG          "-g -Wall -Wno-class-memaccess -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C++ debug compile flags")
set(CMAKE_CXX_FLAGS_MINSIZEREL     "-O1 -DNDEBUG -march=native -Wall -Wno-class-memaccess -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C++ minimum size release compile flags")
set(CMAKE_CXX_FLAGS_RELEASE        "-Ofast  -floop-optimize -falign-loops -falign-labels -falign-functions -falign-jumps -fomit-frame-pointer -DNDEBUG -march=native -Wall -Wno-class-memaccess -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C++ release compile flags")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2  -floop-optimize -falign-loops -falign-labels -falign-functions -falign-jumps -DNDEBUG -g -Wall -Wno-class-memaccess -Wno-deprecated-declarations -Wno-comment" CACHE STRING "Inital C++ release with debug info compile flags")

# Libraries
if(EXISTS $ENV{TBBROOT})
  set(TBB_ROOT_DIR "$ENV{TBBROOT}" CACHE PATH "TBB root directory")
else()
  set(TBB_ROOT_DIR "/opt/intel/tbb" CACHE PATH "TBB root directory")
endif()

# Flags
set(LAPACK_LIBRARIES "-L/usr/lib/aarch64-linux-gnu/openblas" "-llapacke" "-llapack" "-lopenblas")
#set(LAPACK_LIBRARIES "-L/opt/arm/armpl-19.2.0_Cortex-A72_Ubuntu-16.04_arm-hpc-compiler_19.2_aarch64-linux/lib" "-larmpl" "-lamath")

#set(LAPACK_COMPILE_DEFINITIONS MADNESS_LINALG_USE_LAPACKE CACHE STRING "LAPACK preprocessor definitions")
#set(LAPACK_INCLUDE_DIRS ${MKL_ROOT_DIR}/include CACHE STRING "LAPACK include directories")

set(INTEGER4 TRUE CACHE BOOL "Set Fortran integer size to 4 bytes")
set(FORTRAN_INTEGER_SIZE "4" CACHE STRING "Set Fortran integer size to 4 bytes")
