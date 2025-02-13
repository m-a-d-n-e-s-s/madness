cmake_minimum_required (VERSION 3.12.0)  # for FindPython

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

## I know this is obnoxious but I don't know any other way to keep
## things clean.  Make the build verbose so that people see the warning
## messages and fix the issues --- some of them are real
#set(CMAKE_VERBOSE_MAKEFILE ON CACHE BOOL "ON")

# safety net for dev workflow: accidental install will not affect FindOrFetch*
if (NOT DEFINED CACHE{CMAKE_FIND_NO_INSTALL_PREFIX})
  set(CMAKE_FIND_NO_INSTALL_PREFIX ON CACHE BOOL "Whether find_* commands will search CMAKE_INSTALL_PREFIX and CMAKE_STAGING_PREFIX; see https://cmake.org/cmake/help/latest/variable/CMAKE_FIND_NO_INSTALL_PREFIX.html#variable:CMAKE_FIND_NO_INSTALL_PREFIX")
endif()

# Set MADNESS version
set(MADNESS_MAJOR_VERSION 0)
set(MADNESS_MINOR_VERSION 10)
set(MADNESS_MICRO_VERSION 1)
set(MADNESS_VERSION "${MADNESS_MAJOR_VERSION}.${MADNESS_MINOR_VERSION}.${MADNESS_MICRO_VERSION}")

# Declare ourselves ============================================================

project(MADNESS
        VERSION ${MADNESS_VERSION}
        DESCRIPTION "MADNESS: fast, adaptive multiresolution integrodifferential calculus in 1..6 dimensions"
        LANGUAGES CXX C ASM
        HOMEPAGE_URL "https://madness.readthedocs.io/")

# Add source directory =========================================================

add_definitions(-DMAD_ROOT_DIR="${PROJECT_SOURCE_DIR}")

# Add module directory and modules =============================================

list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake/modules/")
include(AddCustomTargetSubproject)
include(AddOptions)
include(AppendFlags)
include(CheckIncludeFile)
include(CheckTypeSize)
include(CheckCXXSourceCompiles)
include(CheckFunctionExists)
include(CMakeDependentOption)
include(AddMADLibrary)
include(AddMADExecutable)
include(AddUnittests)
include(AddScriptedTests)
include(CMakePackageConfigHelpers)
include(CopyTargetProperties)
include(FeatureSummary)
include(RedefaultableOption)
include(FetchContent)

# Check compiler feature support ===============================================

set(CMAKE_CXX_STANDARD 17 CACHE STRING  "C++ ISO Standard version")
if (NOT(CMAKE_CXX_STANDARD EQUAL 17 OR CMAKE_CXX_STANDARD EQUAL 20))
  message(FATAL_ERROR "C++ 2017 ISO Standard or higher is required to compile MADNESS")
endif()
# C++20 is only configurable via compile features with cmake 3.12 and older
if (CMAKE_CXX_STANDARD EQUAL 20 AND CMAKE_VERSION VERSION_LESS 3.12.0)
  message(FATAL_ERROR "C++ 2020 ISO Standard requires CMake version 3.12 or higher")
endif()
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF CACHE BOOL "Whether to use extensions of C++ ISO Standard version")

# workaround for cmake bug: CheckCXXSourceCompiles does not respect CXX_STANDARD
# see https://github.com/OPM/opm-common/blob/master/cmake/Modules/FindCXX11Features.cmake
add_options(CXX ALL_BUILDS ${CMAKE_CXX${CMAKE_CXX_STANDARD}_STANDARD_COMPILE_OPTION})

# Set install paths ============================================================
include(GNUInstallDirs)
set(MADNESS_INSTALL_BINDIR "${CMAKE_INSTALL_BINDIR}"
    CACHE PATH "MADNESS BIN install directory")
set(MADNESS_INSTALL_INCLUDEDIR "${CMAKE_INSTALL_INCLUDEDIR}"
    CACHE PATH "MADNESS INCLUDE install directory")
set(MADNESS_INSTALL_LIBDIR "${CMAKE_INSTALL_LIBDIR}"
    CACHE PATH "MADNESS LIB install directory")
set(MADNESS_INSTALL_DATADIR "${CMAKE_INSTALL_DATAROOTDIR}/madness/${MADNESS_VERSION}/data"
    CACHE PATH "MADNESS DATA install directory")
set(MADNESS_INSTALL_DOCDIR "${CMAKE_INSTALL_DATAROOTDIR}/madness/${MADNESS_VERSION}/doc"
    CACHE PATH "MADNESS DOC install directory")
set(MADNESS_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake/madness"
    CACHE PATH "MADNESS CMAKE install directory")

# Build options  ==============================================================

redefaultable_option(MADNESS_BUILD_LIBRARIES_ONLY "Build only MADNESS libraries" OFF)
if (NOT MADNESS_BUILD_LIBRARIES_ONLY)
  redefaultable_option(MADNESS_BUILD_MADWORLD_ONLY "Build only MADWorld runtime (excluding BLAS/LAPACK interfaces)" OFF)
endif()
redefaultable_option(MADNESS_ENABLE_CEREAL "Support use of Cereal archives as backends for MADNESS serialization" OFF)

# Enable optional libraries ====================================================

# Preload versions/tags of all dependencies ====================================
include(external/versions.cmake)

# Always search for libraries that are 'ON' by default. If a library is not
# found, it is disabled without error. If the library is 'OFF' by default,
# failure to find the library is an error. 

option(ENABLE_MPI "Enable Message Passing Interface (MPI) Library" ON)
#redefaultable_option(ENABLE_MKL "Search for Intel Math Kernel Library (MKL) for BLAS and LAPACK support" ON)
option(ENABLE_MKL "Search for Intel Math Kernel Library (MKL) for BLAS and LAPACK support" ON)
option(ENABLE_ACML "Search for AMD Core Math Library (ACML) for BLAS and LAPACK support" ON)
option(ENABLE_TCMALLOC_MINIMAL "Enable use of tcmalloc_minimal library from Google Performance Tools (gperftools)" OFF)
option(ENABLE_GPERFTOOLS "Enable use of full Google Performance Tools (gperftools)" OFF)
option(ENABLE_LIBUNWIND
       "Forces detection (and possible use) of Libunwind (normally Libunwind will be sought automatically if needed)"
        OFF)
option(ENABLE_PAPI "Enables use of PAPI" OFF)
option(ENABLE_LIBXC "Enables use of the libxc library of density functionals" ON)
option(ENABLE_PCM "Enables use of the polarizable continuum model library" ON)
option(ENABLE_INTEGRATORXX "Enables use of the IntegratorXX for building DFT quadrature grids" ON)

set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES})
check_cxx_source_compiles(
    "
    namespace madness {
      template <typename T>
      void mTxmq(long dimi, long dimj, long dimk,
		 T* __restrict__ c, const T* a, const T* b, long ldb=-1);

      template <>
      void mTxmq(long dimi, long dimj, long dimk, double* __restrict__ c, const  double* a, const double* b, long ldb);
    }
    int main() { double a[1], b[1], c[1]; madness::mTxmq(1,1,1,c,a,b,1); return 0; }
    " HAVE_MTXMQ)

  ####### preferred task backend
set(MADNESS_TASK_VALID_BACKENDS Pthreads TBB PaRSEC)
set(MADNESS_TASK_BACKEND Pthreads CACHE STRING "The backend to use for dispatching tasks")
set_property(CACHE MADNESS_TASK_BACKEND PROPERTY STRINGS ${MADNESS_TASK_VALID_BACKENDS})

if (DEFINED MADNESS_TASK_BACKEND)
  if(NOT MADNESS_TASK_BACKEND IN_LIST MADNESS_TASK_VALID_BACKENDS)
    message(FATAL_ERROR "MADNESS_TASK_BACKEND must be one of ${MADNESS_TASK_VALID_BACKENDS}")
  endif()
  if (MADNESS_TASK_BACKEND STREQUAL TBB)
    set(ENABLE_TBB ON)
  else()
    set(ENABLE_TBB OFF)
  endif()
  if (MADNESS_TASK_BACKEND STREQUAL PaRSEC)
    set(ENABLE_PARSEC ON)
  else()
    set(ENABLE_PARSEC OFF)
  endif()
else (DEFINED MADNESS_TASK_BACKEND)
  redefaultable_option(ENABLE_TBB "Enables use of Intel Thread Building Blocks (TBB) as the task scheduler" OFF)
  option(ENABLE_PARSEC "Enables use of PaRSEC as the task scheduler" OFF)
  if(ENABLE_PARSEC AND ENABLE_TBB)
    message(FATAL_ERROR "TBB and PaRSEC cannot be both enabled.")
  endif()
endif(DEFINED MADNESS_TASK_BACKEND)
if (ENABLE_TBB)
  redefaultable_option(MADNESS_EXPLOIT_TBB_PRIORITY "Enables use of Intel TBB task priorities" OFF)
endif(ENABLE_TBB)

option(ENABLE_ELEMENTAL "Enable Elemental library for distributed-memory linear algebra" OFF)
if (ENABLE_ELEMENTAL)
  if (DEFINED ELEMENTAL_TAG)
    if("${ELEMENTAL_TAG}" STREQUAL "")
      message(FATAL_ERROR "Invalid value given for ELEMENTAL_TAG; specify a valid hash or tag.")
    endif()
    message(STATUS "Will compile Elemental source (tag=${ELEMENTAL_TAG}) in the build directory")
    set(ENABLE_ELEMENTAL_EMBEDDED OFF)
    set(MADNESS_HAS_ELEMENTAL_EMBEDDED 0)
  else (DEFINED ELEMENTAL_TAG)
    message(STATUS "Will use the embedded Elemental source(tag=v0.84)")
    set(ENABLE_ELEMENTAL_EMBEDDED ON)
    set(MADNESS_HAS_ELEMENTAL_EMBEDDED 1)
  endif (DEFINED ELEMENTAL_TAG)
  set(MADNESS_HAS_ELEMENTAL 1)
else (ENABLE_ELEMENTAL)
  set(MADNESS_HAS_ELEMENTAL 0)
endif (ENABLE_ELEMENTAL)
add_feature_info(Elemental ENABLE_ELEMENTAL "compiles parallel linear-algebra library Elemental as part of MADNESS")

option(ENABLE_BOOST "Enable C++ Boost Libraries" OFF)

# Configure options ============================================================

option(ENABLE_GENTENSOR "Enable generic tensor support" OFF)
add_feature_info(GENTENSOR ENABLE_GENTENSOR "enables use of MRA Function compression that allows computing in 6 dimensions")
set(MADNESS_USE_GENTENSOR ${ENABLE_GENTENSOR} CACHE BOOL
    "enables use of MRA Function compression that allows computing in 6 dimensions")

if(ENABLE_GENTENSOR)
  add_definitions(-DUSE_GENTENSOR)
endif()
set(TENSOR_USE_GENTENSOR ${ENABLE_GENTENSOR} CACHE BOOL "enable generic tensors")

option(ENABLE_TASK_PROFILER
    "Enable task profiler that collects per-task start and stop times." OFF)
add_feature_info(TASK_PROFILER ENABLE_TASK_PROFILER "supports task-level tracing of program execution")
set(MADNESS_TASK_PROFILING ${ENABLE_TASK_PROFILER} CACHE BOOL
    "Enable task profiler that collects per-task start and stop times.")

option(ENABLE_WORLD_PROFILE "Enables profiling" OFF)
add_feature_info(WORLD_PROFILE ENABLE_WORLD_PROFILE "supports simple profiling of MADworld runtime")
set(WORLD_PROFILE_ENABLE ${ENABLE_WORLD_PROFILE} CACHE BOOL 
    "Enables world profiling")

option(ENABLE_MEM_STATS "Gather fine-grained memory statistics (expensive)" OFF)
add_feature_info(MEM_STATS ENABLE_MEM_STATS "gather fine-grained memory statistics (expensive)")
set(WORLD_GATHER_MEM_STATS ${ENABLE_MEM_STATS} CACHE BOOL "Gather fine-grained memory statistics (expensive)")

option(ENABLE_MEM_PROFILE "Turn on instrumented aggregate memory profiling (print_meminfo)" OFF)
add_feature_info(MEM_PROFILE ENABLE_MEM_PROFILE "instrumented aggregate memory profiling")
set(WORLD_MEM_PROFILE_ENABLE ${ENABLE_MEM_PROFILE} CACHE BOOL "Turn on instrumented aggregate memory profiling (print_meminfo)")

option(ENABLE_TENSOR_BOUNDS_CHECKING
    "Enable checking of bounds in tensors ... slow but useful for debugging" OFF)
add_feature_info(TENSOR_BOUNDS_CHECKING ENABLE_TENSOR_BOUNDS_CHECKING
    "Enable checking of bounds in tensors ... slow but useful for debugging")
set(TENSOR_BOUNDS_CHECKING ${ENABLE_TENSOR_BOUNDS_CHECKING} CACHE BOOL
    "Enable checking of bounds in tensors ... slow but useful for debugging")

option(ENABLE_TENSOR_INSTANCE_COUNT
    "Enable counting of allocated tensors for memory leak detection" OFF)
add_feature_info(TENSOR_INSTANCE_COUNT ENABLE_TENSOR_INSTANCE_COUNT
    "Enable counting of allocated tensors for memory leak detection")
set(TENSOR_INSTANCE_COUNT CACHE BOOL
    "Enable counting of allocated tensors for memory leak detection")

option(ENABLE_SPINLOCKS
    "Enables use of spinlocks instead of mutexes (faster unless over subscribing processors)" ON)
add_feature_info(SPINLOCKS ENABLE_SPINLOCKS
    "Enables use of spinlocks instead of mutexes (faster unless over subscribing processors)")
set(USE_SPINLOCKS ${ENABLE_SPINLOCKS} CACHE BOOL
    "Enables use of spinlocks instead of mutexes (faster unless over subscribing processors)")

option(ENABLE_NEVER_SPIN
    "Disables use of spinlocks (notably for use inside virtual machines)" OFF)
add_feature_info(NEVER_SPIN ENABLE_NEVER_SPIN
    "Disables use of spinlocks (notably for use inside virtual machines)")
set(NEVER_SPIN ${ENABLE_NEVER_SPIN} CACHE BOOL
    "Disables use of spinlocks (notably for use inside virtual machines)")

option(ENABLE_DQ_PREBUF
    "Enables thread-local buffer for task aggregation to reduce lock contention" ON)
add_feature_info(DQ_PREBUF ENABLE_DQ_PREBUF
    "Enables thread-local buffer for task aggregation to reduce lock contention")
set(MADNESS_DQ_USE_PREBUF ${ENABLE_DQ_PREBUF} CACHE BOOL
    "Enables thread-local buffer for task aggregation to reduce lock contention")

set(MADNESS_DQ_PREBUF_SIZE 20 CACHE STRING "Number of entries in the thread-pool prebuffer for task aggregation to reduce lock contention")

option(ENABLE_BSEND_ACKS 
    "Use MPI Send instead of MPI Bsend for huge message acknowledgements" ON)
add_feature_info(BSEND_ACKS ENABLE_BSEND_ACKS
    "Use MPI Send instead of MPI Bsend for huge message acknowledgements")
set(MADNESS_USE_BSEND_ACKS ${ENABLE_BSEND_ACKS} CACHE BOOL
    "Use MPI Send instead of MPI Bsend for huge message acknowledgements")

option(DISABLE_WORLD_GET_DEFAULT "Disables World::get_default()" OFF)
add_feature_info(WORLD_GET_DEFAULT_DISABLE DISABLE_WORLD_GET_DEFAULT "Disables World::get_default()")
set(WORLD_GET_DEFAULT_DISABLED ${DISABLE_WORLD_GET_DEFAULT} CACHE BOOL 
    "Disables World::get_default()")

option(ENABLE_TASK_DEBUG_TRACE
        "Enable task debug tracing." OFF)
add_feature_info(TASK_DEBUG_TRACE ENABLE_TASK_DEBUG_TRACE "supports debug trace of task engine")
set(MADNESS_TASK_DEBUG_TRACE ${ENABLE_TASK_DEBUG_TRACE} CACHE BOOL
        "Enable task debug tracing.")

option(ENABLE_WORLDOBJECT_FUTURE_TRACE
        "Enable tracing of futures assicuated with WorldObjects." OFF)
add_feature_info(WORLDOBJECT_FUTURE_TRACE ENABLE_WORLDOBJECT_FUTURE_TRACE "supports tracing of futures associated with WorldObjects")
set(MADNESS_WORLDOBJECT_FUTURE_TRACE ${ENABLE_WORLDOBJECT_FUTURE_TRACE} CACHE BOOL
        "Enable tracing of futures assicuated with WorldObjects.")

set(FORTRAN_INTEGER_SIZE 4 CACHE STRING "The fortran integer size (4 or 8 bytes) used for BLAS and LAPACK function calls")
if(NOT (FORTRAN_INTEGER_SIZE EQUAL 4 OR FORTRAN_INTEGER_SIZE EQUAL 8))
  message(FATAL_ERROR "Incorrect fortran integer size '${FORTRAN_INTEGER_SIZE}'\n"
                       "FORTRAN_INTEGER_SIZE must be equal to 4 or 8")
endif()
add_feature_info("FORTRAN_INTEGER_SIZE=${FORTRAN_INTEGER_SIZE}" TRUE "assumes Fortran integers to be ${FORTRAN_INTEGER_SIZE} bytes long")

if (CMAKE_BUILD_TYPE STREQUAL Release OR
    CMAKE_BUILD_TYPE STREQUAL MinSizeRel)
  set (DEFAULT_ASSERTION_TYPE "disable")
else ()
  set (DEFAULT_ASSERTION_TYPE "assert")
endif()
set(ASSERTION_TYPE ${DEFAULT_ASSERTION_TYPE} CACHE STRING "Define MADNESS assertion behavior (abort|assert|disable|throw)")
if(ASSERTION_TYPE STREQUAL "abort")
  set(MADNESS_ASSERTIONS_ABORT 1)
elseif(ASSERTION_TYPE STREQUAL "assert")
  set(MADNESS_ASSERTIONS_ASSERT 1)
elseif(ASSERTION_TYPE STREQUAL "disable")
  set(MADNESS_ASSERTIONS_DISABLE 1)
elseif(ASSERTION_TYPE STREQUAL "throw")
  set(MADNESS_ASSERTIONS_THROW 1)
else()
  message(WARNING "Unsupported ASSERTION_TYPE '${ASSERTION_TYPE}'")
  set(ASSERTION_TYPE "throw")
  set(MADNESS_ASSERTIONS_THROW 1)
endif()
message(STATUS "Assertion type: ${ASSERTION_TYPE}")
add_feature_info("ASSERTION_TYPE=${ASSERTION_TYPE}" TRUE "controls how MADNESS assertions (MADNESS_ASSERT) are handled")

set(MPI_THREAD "multiple" CACHE STRING "Thread level for MPI (multiple|serialized)")
if(MPI_THREAD STREQUAL "multiple")
  set(MADNESS_MPI_THREAD_LEVEL "MPI_THREAD_MULTIPLE")
elseif(MPI_THREAD STREQUAL "serialized")
  set(MADNESS_MPI_THREAD_LEVEL "MPI_THREAD_SERIALIZED")
  if (ENABLE_ELEMENTAL)
    message (FATAL_ERROR "Use of Elemental precludes MPI_THREAD=serialized, must use MPI_THREAD=multiple")
  endif (ENABLE_ELEMENTAL)
else()
  message(FATAL_ERROR "Invalid value for MPI_THREAD '${MPI_THREAD}'; valid values are 'multiple' or 'serialized'")
endif()
add_feature_info("MPI_THREAD=${MPI_THREAD}" TRUE "controls the level of thread-safety support in MPI")

# Enable support for shared libraries ==========================================

redefaultable_option(MADNESS_ASSUMES_ASLR_DISABLED "MADNESS runtime assumes the Address Space Layout Randomization (ASLR) to be disabled" OFF)
add_feature_info(ASSUMES_ASLR_DISABLED MADNESS_ASSUMES_ASLR_DISABLED
    "MADNESS runtime assumes the Address Space Layout Randomization (ASLR) to be disabled")
get_property(SUPPORTS_SHARED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
if (SUPPORTS_SHARED AND MADNESS_ASSUMES_ASLR_DISABLED)
  set(default_BUILD_SHARED_LIBS ON)
else ()
  set(default_BUILD_SHARED_LIBS OFF)
endif()
redefaultable_option(BUILD_SHARED_LIBS "Enable shared libraries" default_BUILD_SHARED_LIBS)
if (BUILD_SHARED_LIBS OR NOT MADNESS_ASSUMES_ASLR_DISABLED)
  set(default_CMAKE_POSITION_INDEPENDENT_CODE ON)
else ()
  set(default_CMAKE_POSITION_INDEPENDENT_CODE OFF)
endif()
redefaultable_option(CMAKE_POSITION_INDEPENDENT_CODE "Default value for POSITION_INDEPENDENT_CODE of targets" default_CMAKE_POSITION_INDEPENDENT_CODE)

if(BUILD_SHARED_LIBS)
  set(BLA_STATIC_DEFAULT FALSE)
  set(CMAKE_MACOSX_RPATH TRUE)
else(BUILD_SHARED_LIBS)
  set(BLA_STATIC_DEFAULT TRUE)
  set(CMAKE_MACOSX_RPATH FALSE)
endif(BUILD_SHARED_LIBS)
if (NOT DEFINED BLA_STATIC)
  set(BLA_STATIC ${BLA_STATIC_DEFAULT})
endif()

set(CMAKE_SKIP_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH  FALSE)
set(CMAKE_SKIP_INSTALL_RPATH FALSE)

# Get host and platform information ============================================

if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  set(ON_A_MAC 1)
  # look for frameworks and appbundles last
  set(CMAKE_FIND_FRAMEWORK LAST)
  set(CMAKE_FIND_APPBUNDLE LAST)
endif()

if (NOT DEFINED MADNESS_CONFIGURATION_USER)
  if($ENV{USERNAME})
    set(MADNESS_CONFIGURATION_USER "$ENV{USERNAME}" CACHE INTERNAL "")
  elseif($ENV{USER})
    set(MADNESS_CONFIGURATION_USER "$ENV{USER}" CACHE INTERNAL "")
  else()
    set(MADNESS_CONFIGURATION_USER "$ENV{USER}" CACHE INTERNAL "")
  endif()
endif (NOT DEFINED MADNESS_CONFIGURATION_USER)
if (NOT DEFINED MADNESS_CONFIGURATION_HOST)
  cmake_host_system_information(RESULT madness_configuration_host QUERY HOSTNAME)
  set(MADNESS_CONFIGURATION_HOST "${madness_configuration_host}" CACHE INTERNAL "")
endif (NOT DEFINED MADNESS_CONFIGURATION_HOST)
if (NOT DEFINED MADNESS_CONFIGURATION_DATE)
  string(TIMESTAMP madness_configuration_date)
  set(MADNESS_CONFIGURATION_DATE "${madness_configuration_date}" CACHE INTERNAL "")
endif (NOT DEFINED MADNESS_CONFIGURATION_DATE)

set(MAD_BIND_DEFAULT "OFF" CACHE STRING "The default binding for threads")

# Check if the target platform is CRAY XE
check_cxx_source_compiles(
    "
    #ifndef __CRAYXE
    #error choke me
    #endif
    int main() { return 0; }
    " HAVE_CRAYXE)
# Check if the target platform is CRAY XE
check_cxx_source_compiles(
    "
    #ifndef __CRAYXT
    #error choke me
    #endif
    int main() { return 0; }
    " HAVE_CRAYXT)
if(HAVE_CRAYXE)
  set(AMD_QUADCORE_TUNE ON CACHE BOOL "Target for tuning mtxmq kernels")
  set(USE_SPINLOCKS ON CACHE BOOL 
      "Enables use of spinlocks instead of mutexs (faster unless over subscribing processors)" FORCE)
  set(MAD_BIND_DEFAULT "1 0 2" CACHE STRING "The default binding for threads" FORCE)
  set(MPI_C_COMPILER cc CACHE STRING "CRAY MPI C compiler")
  set(MPI_CXX_COMPILER CC CACHE STRING "CRAY MPI C++ compiler")
endif()

# Check if the target platform is BG/P
check_cxx_source_compiles(
    "
    #ifndef __bgp__
    #error choke me
    #endif
    int main() { return 0; }
    " HAVE_IBMBGP)

# Check if the target platform is BG/Q
check_cxx_source_compiles(
    "
    #ifndef __bgq__
    #error choke me
    #endif
    int main() { return 0; }
    " HAVE_IBMBGQ)

if(HAVE_IBMBGQ OR HAVE_IBMBGP)
  set(USE_SPINLOCKS ON CACHE BOOL 
      "Enables use of spinlocks instead of mutexs (faster unless over subscribing processors)" FORCE)
endif()


# Check if the target is x86_64
check_cxx_source_compiles(
    "
    #if !(defined(__x86_64__) || defined(_M_X64))
    #error Not x86_64
    #endif
    int main() { return 0; }
    " USE_X86_64_ASM)

if(NOT USE_X86_64_ASM)
  # Check if the target is x86
  check_cxx_source_compiles(
      "
      #if !(defined(__i386) || defined(_M_IX86))
      #error Not x86
      #endif
      int main() { return 0; }
      " USE_X86_32_ASM)
endif()

# (try to) determine C++ ABI
# ABI kinds are named as in https://clang.llvm.org/doxygen/classclang_1_1TargetCXXABI.html
# we only need ABI for serializing member pointers, hence all ARM-based ABIs are represented by same kind
set(MADNESS_CXX_ABI_GenericItanium 1)
set(MADNESS_CXX_ABI_GenericARM 2)
set(MADNESS_CXX_ABI_Microsoft 3)
if (NOT CMAKE_SYSTEM_PROCESSOR)
  set(processor_name ${CMAKE_HOST_SYSTEM_PROCESSOR})
else(NOT CMAKE_SYSTEM_PROCESSOR)
  set(processor_name ${CMAKE_SYSTEM_PROCESSOR})
endif(NOT CMAKE_SYSTEM_PROCESSOR)
if(processor_name MATCHES "^(aarch64.*|AARCH64.*|arm64.*|ARM64.*)" OR processor_name MATCHES "^(arm.*|ARM.*)" OR processor_name MATCHES "^(mips.*|MIPS.*)")
  # on ARM everything, including Windows, targets ARM C++ ABI
  set(MADNESS_CXX_ABI ${MADNESS_CXX_ABI_GenericARM})
elseif(processor_name MATCHES "amd64.*|x86_64.*|AMD64.*" OR processor_name MATCHES "i686.*|i386.*|x86.*" OR processor_name MATCHES "^(powerpc|ppc)64le" OR processor_name MATCHES "^(powerpc|ppc)64")
  if (CMAKE_SYSTEM_NAME STREQUAL Windows) # on x86 Windows all compilers (clang, intel) target Windows C++ ABI
    set(MADNESS_CXX_ABI ${MADNESS_CXX_ABI_Microsoft})
  else (CMAKE_SYSTEM_NAME STREQUAL Windows)  # everyone else uses Itanium C++ ABI
    set(MADNESS_CXX_ABI ${MADNESS_CXX_ABI_GenericItanium})
  endif (CMAKE_SYSTEM_NAME STREQUAL Windows)
else()
  message(FATAL_ERROR "Cannot determine C++ ABI for CMAKE_SYSTEM_PROCESSOR=${CMAKE_SYSTEM_PROCESSOR}, please file an issue")
endif()

# Check compiler feature support ===============================================

# Check for system include files
check_include_file(sys/stat.h HAVE_SYS_STAT_H)
check_include_file(sys/types.h HAVE_SYS_TYPES_H)
check_include_file(unistd.h HAVE_UNISTD_H)
check_include_file(sys/resource.h HAVE_RESOURCE_H)
if(MADNESS_TASK_PROFILING)
  check_include_file(execinfo.h HAVE_EXECINFO_H)
  check_include_file(cxxabi.h HAVE_CXXABI_H)
  if(NOT (HAVE_EXECINFO_H AND HAVE_CXXABI_H))
    message(FATAL_ERROR "Unable to find required header files execinfo.h and/or cxxabi.h")
  endif()
endif()

# Check type support
check_type_size("int64_t" HAVE_INT64_T)
check_type_size("long double" HAVE_LONG_DOUBLE)
check_type_size("long long" HAVE_LONG_LONG)
check_cxx_source_compiles(
    "
    #include <sys/types.h>
    int main() { typedef pid_t test_t; return 0; }
    " SYS_TYPES_H_HAS_PID_T)

# Check function support
check_function_exists(fork HAVE_FORK)
check_function_exists(memset HAVE_MEMSET)
check_function_exists(posix_memalign HAVE_POSIX_MEMALIGN)
if(HAVE_POSIX_MEMALIGN)
  # look for both version of posix_memalign, with and without throw()
  check_cxx_source_compiles(
      "
      #include <stddef.h>
      #include <stdlib.h>
      extern \"C\"  int posix_memalign(void **memptr, size_t alignment, size_t size) throw();
      int main() { void *m; posix_memalign(&m, 16, 1024); }
      " CHECK_STDLIB_H_HAS_POSIX_MEMALIGN_THROW)
  if(NOT CHECK_STDLIB_H_HAS_POSIX_MEMALIGN_THROW)
    check_cxx_source_compiles(
        "
        #include <stddef.h>
        #include <stdlib.h>
        extern \"C\"  int posix_memalign(void **memptr, size_t alignment, size_t size);
        int main() { void *m; posix_memalign(&m, 16, 1024); }
        " CHECK_STDLIB_H_HAS_POSIX_MEMALIGN)
  endif()
  if(NOT CHECK_STDLIB_H_HAS_POSIX_MEMALIGN_THROW AND NOT CHECK_STDLIB_H_HAS_POSIX_MEMALIGN)
    set(MISSING_POSIX_MEMALIGN_PROTO 1)
  endif()

else()
  message(WARNING "posix_memalign NOT FOUND ... enabling override of new/delete ... THIS WILL BE SLOW")
endif()
check_function_exists(pow HAVE_POW)
check_function_exists(random HAVE_RANDOM)
check_function_exists(sleep HAVE_SLEEP)
check_function_exists(strchr HAVE_STRCHR)
# look for both version of posix_memalign, with and without throw()
check_cxx_source_compiles(
    "
    #include <cmath>
    #include <cstdlib>
    long (*absptr)(long) = &std::abs; 
    long a = -1;  
    long b = std::abs(a);
    int main() { return 0; }
    " HAVE_STD_ABS_LONG)
if(NOT HAVE_STD_ABS_LONG)
  check_cxx_source_compiles(
      "
      #include <cmath>
      #include <cstdlib>
      long (*labsptr)(long) = &std::labs; 
      long a = -1l;  
      long b = labs(a);
      int main() { return 0; }
      " HAVE_LABS)
endif()
if(NOT CHECK_STDLIB_H_HAS_POSIX_MEMALIGN_THROW AND NOT CHECK_STDLIB_H_HAS_POSIX_MEMALIGN)
  set(MISSING_POSIX_MEMALIGN_PROTO 1)
endif()

# Check for thread local storage keyword support.
# thread_local, __thread , __thread_local, or __declspec(thread)
if(NOT DEFINED THREAD_LOCAL_KEYWORD)
  foreach(_thread_local_keyword thread_local __thread __thread_local)
    check_cxx_source_compiles(
        "
        ${_thread_local_keyword} int i = 0;
        int main() { i = 1; return 0; }
        " THREAD_LOCAL_SUPPORT)
    if(THREAD_LOCAL_SUPPORT AND _thread_local_keyword STREQUAL "thread_local")
      message(STATUS "Thread local keyword: thread_local")
      unset(THREAD_LOCAL_KEYWORD CACHE)
      break()
    elseif(THREAD_LOCAL_SUPPORT)
      message(STATUS "Thread local keyword: ${_thread_local_keyword}")
      set(THREAD_LOCAL_KEYWORD "${_thread_local_keyword}"
          CACHE STRING "thread local storage keyword, 'thread_local' in C++>11")
      break()
    else()
      unset(THREAD_LOCAL_SUPPORT CACHE)
    endif()
  endforeach()
  
  if(NOT DEFINED THREAD_LOCAL_SUPPORT)
    message(FATAL_ERROR "Unable to detect mandatory support for thread-local storage")
  endif()
endif()


# Check for restrict keyword support
# restrict, __restrict, __restrict__, or _Restrict
if(NOT DEFINED RESTRICT_KEYWORD)
  foreach(_restrict_keyword restrict __restrict __restrict__ _Restrict)
    check_cxx_source_compiles(
        "
        int*  ${_restrict_keyword} i = nullptr;
        int main() { *i = 1; return 0; }
        " RESTRICT_SUPPORT)
    if(RESTRICT_SUPPORT AND _restrict_keyword STREQUAL "restrict")
      unset(RESTRICT_KEYWORD CACHE)
      break()
    elseif(RESTRICT_SUPPORT)
      set(RESTRICT_KEYWORD "${_restrict_keyword}"
          CACHE STRING "C++ equivialent of the C 'restrict' keyword")
      break()
    else()
      unset(RESTRICT_SUPPORT CACHE)
    endif()
  endforeach()

  if(NOT DEFINED RESTRICT_SUPPORT)
    # Set the restrict keyword to nothing so that it is not used
    set(RESTRICT_KEYWORD ""
        CACHE STRING "C++ equivialent of the C 'restrict' keyword")
  endif()
  
  # Print restrict keyword search results
  message(STATUS "Restrict keyword: ${RESTRICT_KEYWORD}")
endif()

check_cxx_source_compiles(
    "
    template <typename T>
    static inline void f(T* a) {};
    template <typename T> void g(T* a) { f(a); }
    template void g(int* a);
    int main() { return 0; }
    " HAVE_UNQUALIFIED_STATIC_DECL)

# Check linker feature support ===============================================

include(CheckDisablePIESupport)
check_disablepie_support(LINKER_HAS_DISABLEPIE_SUPPORT DISABLEPIE_LINKER_FLAG)

# Check for applications =======================================================

find_package(Doxygen)
find_package(LATEX)
find_program(XTERM_EXECUTABLE xterm)
if(XTERM_EXECUTABLE)
  set(HAVE_XTERM 1)
  message(STATUS "Found xterm: ${XTERM_EXECUTABLE}")
endif()
find_program(GDB_EXECUTABLE gdb)
if(GDB_EXECUTABLE)
  set(HAVE_GDB 1)
  message(STATUS "Found gdb: ${GDB_EXECUTABLE}")
endif()
find_program(LLDB_EXECUTABLE lldb)
if(LLDB_EXECUTABLE)
  set(HAVE_LLDB 1)
  message(STATUS "Found lldb: ${LLDB_EXECUTABLE}")
endif()

# ccache for caching builds
find_program(CCACHE ccache)
if(CCACHE)
  mark_as_advanced(CCACHE)
  message (STATUS "Found ccache: ${CCACHE}")
  set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE}" CACHE STRING "Compiler launcher to use for compiling C++")
  set(CMAKE_C_COMPILER_LAUNCHER "${CCACHE}" CACHE STRING "Compiler launcher to use for compiling C")
endif(CCACHE)

# Check for external dependencies ==============================================
include(external/pthread.cmake)
include(external/mpi.cmake)
include(external/papi.cmake)
include(external/libunwind.cmake)
include(external/gperftools.cmake)
include(external/tbb.cmake)
include(external/parsec.cmake)
include(external/boost.cmake)
if (NOT MADNESS_BUILD_MADWORLD_ONLY OR MADNESS_BUILD_LIBRARIES_ONLY)
  include(external/lapack.cmake)
  include(external/libxc.cmake)
  include(external/pcm.cmake)
  include(external/integratorxx.cmake)
endif ()

if (DEFINED ELEMENTAL_TAG)
  include(external/elemental.cmake)
endif (DEFINED ELEMENTAL_TAG)

# Add project subdirectories ===================================================

# predefine targets
add_custom_target(everything)
add_dependencies(everything madness)

if(NOT (MADNESS_BUILD_MADWORLD_ONLY OR MADNESS_BUILD_LIBRARIES_ONLY) AND LAPACK_FOUND)
	add_custom_target_subproject(madness numerical-examples)
	add_custom_target_subproject(madness applications)
endif()

# Create build and install libraries libraries
add_custom_target(madness-libraries)
add_custom_target(install-madness-libraries)
add_dependencies(everything madness-libraries)

include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_BINARY_DIR}/src)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE TRUE)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)
include(CTest)
if(BUILD_TESTING)
  set(MADNESS_HAS_GOOGLE_TEST 1)
  add_custom_target_subproject(madness unittests)
  add_custom_target_subproject(madness scripted_tests)
  add_custom_target_subproject(madness check COMMAND ${CMAKE_CTEST_COMMAND} -V -R "madness/test/" USES_TERMINAL) # to be GNU compatible
  add_custom_target_subproject(madness check-short COMMAND ${CMAKE_CTEST_COMMAND} -L \"short|medium\" -V -R "madness/test/" USES_TERMINAL) # to be GNU compatible
  add_dependencies(everything unittests)
endif()
add_subdirectory(src)
add_subdirectory(doc)

# Generate configure files =====================================================

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/config.h.in
  ${PROJECT_BINARY_DIR}/src/madness/config.h
)

configure_file(
  ${PROJECT_SOURCE_DIR}/cmake/doxygen.cfg.in
  ${PROJECT_BINARY_DIR}/doc/doxygen.cfg
  @ONLY
)

CONFIGURE_FILE(
  ${PROJECT_SOURCE_DIR}/config/MADNESS.pc.in
  ${PROJECT_BINARY_DIR}/config/MADNESS.pc
)

CONFIGURE_FILE( 
  ${PROJECT_SOURCE_DIR}/admin/conda/Makefile.in
  ${PROJECT_BINARY_DIR}/admin/conda/Makefile
  @ONLY
)

CONFIGURE_FILE( 
  ${PROJECT_SOURCE_DIR}/admin/conda/recipe/build.sh.in
  ${PROJECT_BINARY_DIR}/admin/conda/recipe/build.sh
  @ONLY
)

CONFIGURE_FILE( 
  ${PROJECT_SOURCE_DIR}/admin/conda/recipe/meta.yaml.in
  ${PROJECT_BINARY_DIR}/admin/conda/recipe/meta.yaml
  @ONLY
)


# install config files
install(FILES ${PROJECT_BINARY_DIR}/config/MADNESS.pc
    DESTINATION lib/pkgconfig)

# Install some CMake modules
install(FILES
        "${PROJECT_SOURCE_DIR}/cmake/modules/FindTBB.cmake"
        DESTINATION "${MADNESS_INSTALL_CMAKEDIR}/modules"
        COMPONENT madness-config)

# Create the version file
write_basic_package_version_file(madness-config-version.cmake
  VERSION ${MADNESS_VERSION} COMPATIBILITY AnyNewerVersion)

# Create the targets file for madness export set
# N.B. ALL targets end up in madness export set, even if building only madworld
export(EXPORT madness FILE "${PROJECT_BINARY_DIR}/madness-targets.cmake")

# Create the configure file
configure_package_config_file(cmake/madness-config.cmake.in
    "${MADNESS_BINARY_DIR}/madness-config.cmake"
  INSTALL_DESTINATION "${MADNESS_INSTALL_CMAKEDIR}"
  PATH_VARS CMAKE_INSTALL_PREFIX MADNESS_INSTALL_BINDIR 
            MADNESS_INSTALL_INCLUDEDIR MADNESS_INSTALL_LIBDIR
            MADNESS_INSTALL_DATADIR MADNESS_INSTALL_DOCDIR 
            MADNESS_INSTALL_CMAKEDIR)

# Install config, version, and target files
install(EXPORT madness
   	FILE "madness-targets.cmake"
   	DESTINATION "${MADNESS_INSTALL_CMAKEDIR}" 
   	COMPONENT madness-config)
install(FILES
    "${MADNESS_BINARY_DIR}/madness-config.cmake"
    "${MADNESS_BINARY_DIR}/madness-config-version.cmake"
    DESTINATION "${MADNESS_INSTALL_CMAKEDIR}" 
    COMPONENT madness-config)
add_custom_target(install-madness-config
     COMMAND ${CMAKE_COMMAND} -DCOMPONENT=madness-config -P ${PROJECT_BINARY_DIR}/cmake_install.cmake
     COMMENT "Installing MADNESS config components"
     USES_TERMINAL)

feature_summary(WHAT ALL
                DESCRIPTION "=== MADNESS Package/Feature Info ===")

