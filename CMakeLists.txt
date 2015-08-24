cmake_minimum_required (VERSION 2.8.12)
project(MADNESS CXX C ASM)

# Set MADNESS version
set(MADNESS_MAJOR_VERSION 0)
set(MADNESS_MINOR_VERSION 10)
set(MADNESS_MICRO_VERSION 0)
set(MADNESS_VERSION "${MADNESS_MAJOR_VERSION}.${MADNESS_MINOR_VERSION}.${MADNESS_MICRO_VERSION}")

# Add module directory and modules =============================================

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules/)
include(AppendFlags)
include(CheckCXX11Support)
include(CheckIncludeFile)
include(CheckTypeSize)
include(CheckCXXSourceCompiles)
include(CheckFunctionExists)
include(CMakeDependentOption)
include(AddUnittests)
include(AddMADLibrary)
include(CMakePackageConfigHelpers)
include(CopyTargetProperties)

# Set install paths ============================================================

set(MADNESS_INSTALL_BINDIR "bin" 
    CACHE PATH "MADNESS binary install directory")
set(MADNESS_INSTALL_INCLUDEDIR "include" 
    CACHE PATH "MADNESS INCLUDE install directory")
set(MADNESS_INSTALL_LIBDIR "lib"
    CACHE PATH "MADNESS LIB install directory")
set(MADNESS_INSTALL_DATADIR "share"
    CACHE PATH "MADNESS DATA install directory")
set(MADNESS_INSTALL_DOCDIR "${MADNESS_INSTALL_DATADIR}/doc"
    CACHE PATH "MADNESS DOC install directory")
set(MADNESS_INSTALL_CMAKEDIR "lib/cmake/Madness"
    CACHE PATH "MADNESS CMAKE install directory")

# Enable optional libraries ====================================================

# Always search for libraries that are 'ON' by default. If a library is not
# found, it is disable without error. If the library is 'OFF' by default,
# failure to find the library is an error. 

option(WITH_MPI "Enable MPI" ON)
option(WITH_TBB "Enables use of Intel TBB as the task scheduler" ON)
option(WITH_MKL "Search for Intel MKL for BLAS and LAPACK support" ON)
option(WITH_ACML "Search for AMD math library for BLAS and LAPACK support" ON)
option(ENABLE_ELEMENTAL "Enable Elemental" OFF)
option(WITH_GPERFTOOLS "Enable use of gperftools" ON)
option(WITH_LIBXC "Enables use of the libxc library of density functionals" ON)
option(WITH_PAPI "Enables use of PAPI" OFF)

# Configure options ============================================================

option(ENABLE_GENTENSOR "Enable gentensor C and C++ compilers." OFF)
if(ENABLE_GENTENSOR)
  add_definitions(-DUSE_GENTENSOR)
endif()

option(ENABLE_MEM_STATS "Enable generic tensors" OFF)
set(WORLD_GATHER_MEM_STATS ${ENABLE_MEM_STATS} CACHE BOOL "Enable generic tensors")

option(ENABLE_TASK_PROFILER
    "Enable task profiler that collects per-task start and stop times." OFF)
set(MADNESS_TASK_PROFILING ${ENABLE_TASK_PROFILER} CACHE BOOL
    "Enable task profiler that collects per-task start and stop times.")

option(ENABLE_WORLD_PROFILE "Enables profiling" OFF)
set(WORLD_PROFILE_ENABLE ${ENABLE_WORLD_PROFILE} CACHE BOOL 
    "Enables world profiling")

option(ENABLE_MEMSTATS "Gather memory statistics (expensive)" OFF)
set(WORLD_GATHER_MEM_STATS ${WORLD_GATHER_MEM_STATS})

option(ENABLE_TENSOR_BOUNDS_CHECKING
    "Enable checking of bounds in tensors ... slow but useful for debugging" OFF)
set(TENSOR_BOUNDS_CHECKING ${ENABLE_TENSOR_BOUNDS_CHECKING} CACHE BOOL
    "Enable checking of bounds in tensors ... slow but useful for debugging")

option(ENABLE_TENSOR_INSTANCE_COUNT
    "Enable counting of allocated tensors for memory leak detection" OFF)
set(TENSOR_INSTANCE_COUNT CACHE BOOL
    "Enable counting of allocated tensors for memory leak detection")

option(ENABLE_SPINLOCKS
    "Enables use of spinlocks instead of mutexes (faster unless over subscribing processors)" ON)
set(USE_SPINLOCKS ${ENABLE_SPINLOCKS} CACHE BOOL
    "Enables use of spinlocks instead of mutexs (faster unless over subscribing processors)")

option(ENABLE_NEVER_SPIN
    "Disables use of spinlocks (notably for use inside virtual machines" OFF)
set(NEVER_SPIN ${ENABLE_NEVER_SPIN} CACHE BOOL
    "Disables use of spinlocks (notably for use inside virtual machines")

option(ENABLE_BSEND_ACKS 
    "Use MPI Send instead of MPI Bsend for huge message acknowledgements" ON)
set(MADNESS_USE_BSEND_ACKS ${ENABLE_BSEND_ACKS} CACHE BOOL
    "Use MPI Send instead of MPI Bsend for huge message acknowledgements")

option(ENABLE_UNITTESTS "Enables unit tests targets" ON)

set(FORTRAN_INTEGER_SIZE 4 CACHE STRING "The fortran integer size (4 or 8 bytes) used for BLAS and LAPACK function calls")
if(NOT (FORTRAN_INTEGER_SIZE EQUAL 4 OR FORTRAN_INTEGER_SIZE EQUAL 8))
  message(FATAL_ERROR "Incorrect fortran integer size '${FORTRAN_INTEGER_SIZE}'\n"
                       "FORTRAN_INTEGER_SIZE must be equal to 4 or 8")
endif()

set(ASSERTION_TYPE "throw" CACHE STRING "Define MADNESS assertion behavior (abort|assert|disable|throw)")
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

set(MPI_THREAD "multiple" CACHE STRING "Thread level for MPI (multiple|serialized)")
if(MPI_THREAD STREQUAL "multiple")
  set(MADNESS_MPI_THREAD_LEVEL "MPI_THREAD_MULTIPLE")
elseif(MPI_THREAD STREQUAL "serialized")
  set(MADNESS_MPI_THREAD_LEVEL "MPI_THREAD_SERIALIZED")
else()
  message(FATAL_ERROR "Invalid value for MPI_THREAD '${MPI_THREAD}'; valid valuse are 'multiple' or 'serialized'")
endif()

# Enable support for shared libraries ==========================================

get_property(SUPPORTS_SHARED GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS)
cmake_dependent_option(BUILD_SHARED_LIBS "Enable shared libraries" ON
      "SUPPORTS_SHARED" OFF)
if(BUILD_SHARED_LIBS)
  set(BLA_STATIC FALSE)
  set(CMAKE_MACOSX_RPATH TRUE)
else()
  set(BLA_STATIC TURE)
  set(CMAKE_MACOSX_RPATH FALSE)
endif()

set(CMAKE_SKIP_RPATH FALSE)
set(CMAKE_SKIP_BUILD_RPATH  FALSE)

# Get host and platform information ============================================

if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
  set(ON_A_MAC 1)
endif()

if($ENV{USERNAME})
set(MADNESS_CONFIGURATION_USER "$ENV{USERNAME}")
elseif($ENV{USER})
set(MADNESS_CONFIGURATION_USER "$ENV{USER}")
else()
set(MADNESS_CONFIGURATION_USER "$ENV{USER}")
endif()
cmake_host_system_information(RESULT MADNESS_CONFIGURATION_HOST QUERY HOSTNAME)
string(TIMESTAMP MADNESS_CONFIGURATION_DATE)

set(MAD_BIND_DEFAULT "-1 -1 -1" CACHE STRING "The default binding for threads")

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

# Check compiler feature support ===============================================

check_cxx11_support(COMPILER_HAS_CXX11_SUPPORT CXX11_COMPILE_FLAG)
if(COMPILER_HAS_CXX11_SUPPORT)
  if(CXX11_COMPILE_FLAG)
    append_flags(CMAKE_CXX_FLAGS "${CXX11_COMPILE_FLAG}")
  endif()
else()
  message(FATAL_ERROR "${CMAKE_CXX_COMPILER} does not support C++11.")
endif()

# Check for system include files
check_include_file(sys/stat.h HAVE_SYS_STAT_H)
check_include_file(sys/types.h HAVE_SYS_TYPES_H)
check_include_file(unistd.h HAVE_UNISTD_H)
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
    " CHECK_SYS_TYPES_H_HAS_PID_T)
if(NOT CHECK_SYS_TYPES_H_HAS_PID_T)
  set(pid_t 1)
endif()

# Check function support
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
      unset(THREAD_LOCAL_KEYWORD CACHE)
      break()
    elseif(THREAD_LOCAL_SUPPORT)
      set(THREAD_LOCAL_KEYWORD "${_thread_local_keyword}"
          CACHE STRING "thread local storage keyword, 'thread_local' in C++11")
      break()
    else()
      unset(THREAD_LOCAL_SUPPORT CACHE)
    endif()
  endforeach()
  
  if(NOT DEFINED THREAD_LOCAL_SUPPORT)
      set(THREAD_LOCAL_KEYWORD ""
          CACHE STRING "thread local storage keyword, 'thread_local' in C++11")
  endif()
  
  # Print thread_local keyword search results
  message(STATUS "Thread local keyword: ${THREAD_LOCAL_KEYWORD}")
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
find_package(PythonInterp)

# Check for external dependencies ==============================================

include(external/pthread.cmake)
include(external/mpi.cmake)
include(external/papi.cmake)
include(external/lapack.cmake)
include(external/libxc.cmake)
include(external/libunwind.cmake)
include(external/gperftools.cmake)
include(external/tbb.cmake)
#include(external/elemental.cmake)

# Add project subdirectories ===================================================

# Create build and install libraries libraries
add_custom_target(libraries)
add_custom_target(install-libraries)

include_directories(${CMAKE_SOURCE_DIR}/src ${CMAKE_BINARY_DIR}/src)
set(CMAKE_INCLUDE_CURRENT_DIR_IN_INTERFACE TRUE)
set(CMAKE_INCLUDE_CURRENT_DIR TRUE)
if(ENABLE_UNITTESTS)
  include(CTest)
  enable_testing()
  set(MADNESS_HAS_GOOGLE_TEST 1)
  add_custom_target(unittests)
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

# Create the version file
write_basic_package_version_file(madness-config-version.cmake
  VERSION ${MADNESS_VERSION} COMPATIBILITY AnyNewerVersion)

# Create the targets file
export(EXPORT madness
  FILE "${PROJECT_BINARY_DIR}/madness-targets.cmake")

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
    COMPONENT madness)
install(FILES
    "${MADNESS_BINARY_DIR}/madness-config.cmake"
    "${MADNESS_BINARY_DIR}/madness-config-version.cmake"
    DESTINATION "${MADNESS_INSTALL_CMAKEDIR}" 
    COMPONENT madness)
