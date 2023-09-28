# Find BLAS and LAPACK.

set(lapack_is_optional 0)
set(missing_lapack_message_level "FATAL_ERROR")

# if building just the runtime, included this by mistake, warn and make optional
if (MADNESS_BUILD_MADWORLD_ONLY)
  message(WARNING "MADNESS_BUILD_MADWORLD_ONLY=ON, but included external/lapack.cmake; must be error in CMakeLists.txt")
  set(lapack_is_optional 1)
  set(missing_lapack_message_level "STATUS")
endif (MADNESS_BUILD_MADWORLD_ONLY)

include(CheckCFortranFunctionExists)
include(CMakePushCheckState)
include(CheckCXXSourceCompiles)

if(NOT LAPACK_LIBRARIES)
  set(USER_LAPACK_LIBRARIES FALSE)

  if(ENABLE_MKL)
    find_package(MKL)
    
    if(MKL_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
      set(HAVE_INTEL_MKL 1)
      set(LAPACK_COMPILE_DEFINITIONS MADNESS_LINALG_USE_LAPACKE)
      set(LAPACK_INCLUDE_DIRS ${MKL_INCLUDE_DIRS})
    endif()
  endif()
  
  if(ENABLE_ACML AND NOT LAPACK_FOUND)
    find_package(ACML)
    
    if(ACML_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${ACML_LIBRARIES})
      set(HAVE_ACML 1)
    endif()
  endif()
  
  # Search for system specific BLAS/LAPACK checks
  if(NOT LAPACK_FOUND AND CMAKE_SYSTEM_NAME MATCHES "Darwin")
    # Accelerate is always present, so no need to search
    set(LAPACK_LIBRARIES "-framework Accelerate")
    set(LAPACK_FOUND TRUE)
  endif()
  
  # Search for netlib lapack and blas libraries
  if(NOT LAPACK_FOUND)
    find_library(LAPACK_lapack_LIBRARY lapack)
    find_library(LAPACK_blas_LIBRARY blas)
    
    if(LAPACK_lapack_LIBRARY AND LAPACK_blas_LIBRARY)
      set(LAPACK_LIBRARIES ${LAPACK_lapack_LIBRARY} ${LAPACK_blas_LIBRARY})
      set(LAPACK_FOUND TRUE)
    endif()
  endif()  

else()
  set(USER_LAPACK_LIBRARIES TRUE)
endif()

cmake_push_check_state()

# process LAPACK_LIBRARIES for CMAKE_REQUIRED_LIBRARIES (this is likely only to work with Makefile generator):
# 1. get rid of the surrounding quotes
string(REGEX REPLACE "\"" "" PROCESSED_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
# 2. convert a space-separated string of libs into a list
string(REGEX REPLACE " " ";" PROCESSED_LAPACK_LIBRARIES "${PROCESSED_LAPACK_LIBRARIES}")
# 3. restore (and protect!) the space in "-framework X"
string(REGEX REPLACE "-framework;(.*)" "-framework\\\\ \\1" PROCESSED_LAPACK_LIBRARIES "${PROCESSED_LAPACK_LIBRARIES}")
#message(STATUS "PROCESSED_LAPACK_LIBRARIES=${PROCESSED_LAPACK_LIBRARIES}")
set(CMAKE_REQUIRED_LIBRARIES ${PROCESSED_LAPACK_LIBRARIES} ${CMAKE_REQUIRED_LIBRARIES}
        Threads::Threads)

# Verify that we can link against BLAS
check_c_fortran_function_exists(sgemm BLAS_WORKS)

if(BLAS_WORKS)
  message(STATUS "A library with BLAS API found.")
else()
  message("${missing_lapack_message_level}" "Unable to link against BLAS function. Specify the BLAS library in LAPACK_LIBRARIES.")
endif()

# Verify that we can link against LAPACK
check_c_fortran_function_exists(cheev LAPACK_WORKS)

if(LAPACK_WORKS)
  message(STATUS "A library with LAPACK API found.")
else()
  message("${missing_lapack_message_level}" "Unable to link against LAPACK function. Specify the LAPACK library in LAPACK_LIBRARIES.")
endif()

set(LAPACK_FOUND TRUE)
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

# introspect LAPACK_LIBRARIES given by the user
if (USER_LAPACK_LIBRARIES)

  # check for MKL
  check_function_exists(mkl_get_version USER_LAPACK_LIBRARIES_IS_MKL)
  if(USER_LAPACK_LIBRARIES_IS_MKL)
    message(STATUS "User-defined LAPACK_LIBRARIES provides an MKL library")
    set(HAVE_INTEL_MKL 1)
    # ensure that MADNESS_LINALG_USE_LAPACKE is defined
    list(APPEND LAPACK_COMPILE_DEFINITIONS MADNESS_LINALG_USE_LAPACKE)
    list(REMOVE_DUPLICATES LAPACK_COMPILE_DEFINITIONS)

  else(USER_LAPACK_LIBRARIES_IS_MKL)
    # check for ACML
    check_function_exists(acmlversion USER_LAPACK_LIBRARIES_IS_ACML)

    if(USER_LAPACK_LIBRARIES_IS_ACML)
      message(STATUS "User-defined LAPACK_LIBRARIES provides an ACML library")
      set(HAVE_ACML 1)
    endif(USER_LAPACK_LIBRARIES_IS_ACML)
  endif(USER_LAPACK_LIBRARIES_IS_MKL)

  # LAPACK_LIBRARIES might include IMPORTED targets that are not globally available
  if (LAPACK_LIBRARIES MATCHES OpenMP::OpenMP_C AND NOT TARGET OpenMP::OpenMP_C)
    find_package(OpenMP REQUIRED COMPONENTS C)
  endif()
  if (LAPACK_LIBRARIES MATCHES Threads::Threads AND NOT TARGET Threads::Threads)
    find_package(Threads REQUIRED)
  endif()

endif(USER_LAPACK_LIBRARIES)

cmake_pop_check_state()

# Set the fortran mangling scheme.
if(LAPACK_WORKS STREQUAL "cheev_")
  set(FORTRAN_LINKAGE_LCU 1)
elseif(LAPACK_WORKS STREQUAL "cheev")
  set(FORTRAN_LINKAGE_LC 1)
elseif(LAPACK_WORKS STREQUAL "cheev__")
  set(FORTRAN_LINKAGE_LCUU 1)
elseif(LAPACK_WORKS STREQUAL "CHEEV")
  set(FORTRAN_LINKAGE_UC 1)
elseif(LAPACK_WORKS STREQUAL "CHEEV_")
  set(FORTRAN_LINKAGE_UCU 1)
endif()

# unquote LAPACK_COMPILE_OPTIONS, LAPACK_INCLUDE_DIRS, and LAPACK_COMPILE_DEFINITIONS also
string(REGEX REPLACE "\"" "" LAPACK_COMPILE_OPTIONS "${LAPACK_COMPILE_OPTIONS}")
string(REGEX REPLACE "\"" "" LAPACK_INCLUDE_DIRS "${LAPACK_INCLUDE_DIRS}")
string(REGEX REPLACE "\"" "" LAPACK_COMPILE_DEFINITIONS "${LAPACK_COMPILE_DEFINITIONS}")

# epilogue: final sanity checks
if(USER_LAPACK_LIBRARIES_IS_MKL)
  cmake_push_check_state()
  # ensure that can include mkl.h
  set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIRS})
  foreach(def ${LAPACK_COMPILE_DEFINITIONS})
    set(CMAKE_REQUIRED_DEFINITIONS "${CMAKE_REQUIRED_DEFINITIONS};-D${def}")
  endforeach()
  set(CMAKE_REQUIRED_FLAGS ${LAPACK_COMPILE_OPTIONS})
  check_cxx_source_compiles(
          "
  #include <mkl.h>
  int main(int argc, char** argv) {
    return 0;
  }
  "  MADNESS_CAN_INCLUDE_MKL_H)
  if (NOT MADNESS_CAN_INCLUDE_MKL_H)
    message("${missing_lapack_message_level}" "User-provided LAPACK_LIBRARIES provides MKL but cannot include its headers; ensure that corresponding LAPACK_INCLUDE_DIRS, LAPACK_COMPILE_DEFINITIONS, or LAPACK_COMPILE_OPTIONS were provided")
  endif()

  cmake_pop_check_state()
endif(USER_LAPACK_LIBRARIES_IS_MKL)

