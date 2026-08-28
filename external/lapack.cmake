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
include(CheckFunctionExists)

# Helper function: Check if a library or list of libraries is OpenBLAS
function(madness_check_is_openblas _lib_list _out_var)
  set(${_out_var} FALSE PARENT_SCOPE)
  # 1. Path check (including resolved symlinks)
  foreach(_lib ${_lib_list})
    if(EXISTS "${_lib}")
      get_filename_component(_real_path "${_lib}" REALPATH)
      if(_real_path MATCHES "openblas" OR _lib MATCHES "openblas")
        set(${_out_var} TRUE PARENT_SCOPE)
        return()
      endif()
    elseif(_lib MATCHES "openblas")
      set(${_out_var} TRUE PARENT_SCOPE)
      return()
    endif()
  endforeach()
  # 2. Symbol check by test compilation/linking
  cmake_push_check_state()
  set(CMAKE_REQUIRED_LIBRARIES ${_lib_list})
  check_function_exists(openblas_get_config _has_openblas_symbol)
  cmake_pop_check_state()
  if(_has_openblas_symbol)
    set(${_out_var} TRUE PARENT_SCOPE)
  endif()
endfunction()

if(NOT LAPACK_LIBRARIES)
  set(USER_LAPACK_LIBRARIES FALSE)

  # 1. Search for Intel MKL
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

  # 2. Search for Arm Performance Libraries (ARMPL)
  if(ENABLE_ARMPL AND NOT LAPACK_FOUND)
    find_package(ARMPL)
    if(ARMPL_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${ARMPL_LIBRARIES})
      set(LAPACK_INCLUDE_DIRS ${ARMPL_INCLUDE_DIRS})
      set(HAVE_ARMPL 1)
    endif()
  endif()

  # 3. Search for BLIS (serial) + LAPACK (Netlib)
  if(ENABLE_BLIS AND NOT LAPACK_FOUND)
    find_package(BLIS)
    if(BLIS_FOUND)
      # BLIS provides BLAS; find a compatible LAPACK library (e.g. Netlib LAPACK)
      # Search specifically in Netlib directories first to avoid picking up OpenBLAS symlinks
      find_library(BLIS_LAPACK_LIBRARY NAMES lapack
        HINTS
          /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}/lapack
          /usr/lib/aarch64-linux-gnu/lapack
          /usr/lib/x86_64-linux-gnu/lapack
          /usr/lib/lapack
          /usr/lib64/lapack
          /usr/local/lib/lapack
          /usr/local/lib
          /usr/lib
        PATH_SUFFIXES lapack
      )
      if(BLIS_LAPACK_LIBRARY AND NOT ENABLE_OPENBLAS)
        madness_check_is_openblas("${BLIS_LAPACK_LIBRARY}" _lapack_is_openblas)
        if(_lapack_is_openblas)
          message(STATUS "Rejecting OpenBLAS library found for LAPACK: ${BLIS_LAPACK_LIBRARY}")
          unset(BLIS_LAPACK_LIBRARY CACHE)
          unset(BLIS_LAPACK_LIBRARY)
        endif()
      endif()

      if(BLIS_LAPACK_LIBRARY)
        set(LAPACK_FOUND TRUE)
        set(LAPACK_LIBRARIES ${BLIS_LIBRARIES} ${BLIS_LAPACK_LIBRARY})
        set(LAPACK_INCLUDE_DIRS ${BLIS_INCLUDE_DIRS})
        set(HAVE_BLIS 1)
      endif()
    endif()
  endif()

  # 4. Search for ACML
  if(ENABLE_ACML AND NOT LAPACK_FOUND)
    find_package(ACML)
    if(ACML_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${ACML_LIBRARIES})
      set(HAVE_ACML 1)
    endif()
  endif()

  # 5. Search for system specific BLAS/LAPACK checks (Darwin / Accelerate)
  if(NOT LAPACK_FOUND AND CMAKE_SYSTEM_NAME MATCHES "Darwin")
    # Accelerate is always present, so no need to search
    set(LAPACK_LIBRARIES "-framework Accelerate")
    set(LAPACK_FOUND TRUE)
  endif()

  # 6. Search for Netlib lapack and blas libraries
  if(NOT LAPACK_FOUND)
    find_library(LAPACK_lapack_LIBRARY NAMES lapack
      HINTS
        /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}/lapack
        /usr/lib/aarch64-linux-gnu/lapack
        /usr/lib/x86_64-linux-gnu/lapack
        /usr/lib/lapack
        /usr/lib64/lapack
        /usr/local/lib/lapack
        /usr/local/lib
        /usr/lib
      PATH_SUFFIXES lapack
    )
    find_library(LAPACK_blas_LIBRARY NAMES blas
      HINTS
        /usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}/blas
        /usr/lib/aarch64-linux-gnu/blas
        /usr/lib/x86_64-linux-gnu/blas
        /usr/lib/blas
        /usr/lib64/blas
        /usr/local/lib/blas
        /usr/local/lib
        /usr/lib
      PATH_SUFFIXES blas
    )

    if(LAPACK_lapack_LIBRARY AND NOT ENABLE_OPENBLAS)
      madness_check_is_openblas("${LAPACK_lapack_LIBRARY}" _lapack_is_ob)
      if(_lapack_is_ob)
        message(STATUS "Rejecting OpenBLAS library found for LAPACK: ${LAPACK_lapack_LIBRARY}")
        unset(LAPACK_lapack_LIBRARY CACHE)
        unset(LAPACK_lapack_LIBRARY)
      endif()
    endif()

    if(LAPACK_blas_LIBRARY AND NOT ENABLE_OPENBLAS)
      madness_check_is_openblas("${LAPACK_blas_LIBRARY}" _blas_is_ob)
      if(_blas_is_ob)
        message(STATUS "Rejecting OpenBLAS library found for BLAS: ${LAPACK_blas_LIBRARY}")
        unset(LAPACK_blas_LIBRARY CACHE)
        unset(LAPACK_blas_LIBRARY)
      endif()
    endif()

    if(LAPACK_lapack_LIBRARY AND LAPACK_blas_LIBRARY)
      set(LAPACK_LIBRARIES ${LAPACK_lapack_LIBRARY} ${LAPACK_blas_LIBRARY})
      set(LAPACK_FOUND TRUE)
    endif()
  endif()

  # 7. Fallback to OpenBLAS only if explicitly requested via ENABLE_OPENBLAS=ON
  if(ENABLE_OPENBLAS AND NOT LAPACK_FOUND)
    find_library(LAPACK_openblas_LIBRARY NAMES openblas)
    if(LAPACK_openblas_LIBRARY)
      set(LAPACK_LIBRARIES ${LAPACK_openblas_LIBRARY})
      set(LAPACK_FOUND TRUE)
    endif()
  endif()

else()
  set(USER_LAPACK_LIBRARIES TRUE)
endif()

cmake_push_check_state()

# process LAPACK_LIBRARIES for CMAKE_REQUIRED_LIBRARIES
string(REGEX REPLACE "\"" "" PROCESSED_LAPACK_LIBRARIES "${LAPACK_LIBRARIES}")
string(REGEX REPLACE " " ";" PROCESSED_LAPACK_LIBRARIES "${PROCESSED_LAPACK_LIBRARIES}")
string(REGEX REPLACE "-framework;(.*)" "-framework\\\\ \\1" PROCESSED_LAPACK_LIBRARIES "${PROCESSED_LAPACK_LIBRARIES}")

set(CMAKE_REQUIRED_LIBRARIES ${PROCESSED_LAPACK_LIBRARIES} ${CMAKE_REQUIRED_LIBRARIES}
        Threads::Threads)
if(LAPACK_INCLUDE_DIRS)
  set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIRS} ${CMAKE_REQUIRED_INCLUDES})
endif()

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

# Check for OpenBLAS rejection/warning
check_function_exists(openblas_get_config LAPACK_IS_OPENBLAS)
if(LAPACK_IS_OPENBLAS)
  if(USER_LAPACK_LIBRARIES OR ENABLE_OPENBLAS)
    message(WARNING "OpenBLAS detected in LAPACK_LIBRARIES. Note: OpenBLAS is not safe for multithreaded concurrent calls on ARM.")
  else()
    message("${missing_lapack_message_level}" "Auto-detected LAPACK/BLAS library is OpenBLAS, which is not thread-safe on ARM. Please install ARMPL, BLIS serial, or Netlib LAPACK, or configure with -DENABLE_OPENBLAS=ON to override.")
  endif()
endif()

set(LAPACK_FOUND TRUE)
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

# Introspect LAPACK_LIBRARIES (both user-specified and auto-detected)
check_function_exists(mkl_get_version USER_LAPACK_LIBRARIES_IS_MKL)
if(USER_LAPACK_LIBRARIES_IS_MKL)
  message(STATUS "LAPACK provides an MKL library")
  set(HAVE_INTEL_MKL 1)
  list(APPEND LAPACK_COMPILE_DEFINITIONS MADNESS_LINALG_USE_LAPACKE)
  list(REMOVE_DUPLICATES LAPACK_COMPILE_DEFINITIONS)
endif()

check_function_exists(armplversion USER_LAPACK_LIBRARIES_IS_ARMPL)
if(USER_LAPACK_LIBRARIES_IS_ARMPL)
  message(STATUS "LAPACK provides an ARMPL library")
  set(HAVE_ARMPL 1)
endif()

check_function_exists(bli_info_get_version_str USER_LAPACK_LIBRARIES_IS_BLIS)
if(USER_LAPACK_LIBRARIES_IS_BLIS)
  message(STATUS "LAPACK provides a BLIS library")
  set(HAVE_BLIS 1)
endif()

check_function_exists(acmlversion USER_LAPACK_LIBRARIES_IS_ACML)
if(USER_LAPACK_LIBRARIES_IS_ACML)
  message(STATUS "LAPACK provides an ACML library")
  set(HAVE_ACML 1)
endif()

# LAPACK_LIBRARIES might include IMPORTED targets that are not globally available
if(LAPACK_LIBRARIES MATCHES OpenMP::OpenMP_C AND NOT TARGET OpenMP::OpenMP_C)
  find_package(OpenMP REQUIRED COMPONENTS C)
endif()
if(LAPACK_LIBRARIES MATCHES Threads::Threads AND NOT TARGET Threads::Threads)
  find_package(Threads REQUIRED)
endif()

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

# unquote LAPACK_COMPILE_OPTIONS, LAPACK_INCLUDE_DIRS, and LAPACK_COMPILE_DEFINITIONS
string(REGEX REPLACE "\"" "" LAPACK_COMPILE_OPTIONS "${LAPACK_COMPILE_OPTIONS}")
string(REGEX REPLACE "\"" "" LAPACK_INCLUDE_DIRS "${LAPACK_INCLUDE_DIRS}")
string(REGEX REPLACE "\"" "" LAPACK_COMPILE_DEFINITIONS "${LAPACK_COMPILE_DEFINITIONS}")

# epilogue: final sanity checks
if(USER_LAPACK_LIBRARIES_IS_MKL)
  cmake_push_check_state()
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
    " MADNESS_CAN_INCLUDE_MKL_H)
  if(NOT MADNESS_CAN_INCLUDE_MKL_H)
    message("${missing_lapack_message_level}" "LAPACK provides MKL but cannot include its headers; ensure that corresponding LAPACK_INCLUDE_DIRS, LAPACK_COMPILE_DEFINITIONS, or LAPACK_COMPILE_OPTIONS were provided")
  endif()
  cmake_pop_check_state()
endif(USER_LAPACK_LIBRARIES_IS_MKL)

if(USER_LAPACK_LIBRARIES_IS_ARMPL AND LAPACK_INCLUDE_DIRS)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_FLAGS ${LAPACK_COMPILE_OPTIONS})
  check_cxx_source_compiles(
    "
    #include <armpl.h>
    int main(int argc, char** argv) {
      return 0;
    }
    " MADNESS_CAN_INCLUDE_ARMPL_H)
  cmake_pop_check_state()
endif()

if(USER_LAPACK_LIBRARIES_IS_BLIS AND LAPACK_INCLUDE_DIRS)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_INCLUDES ${LAPACK_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_FLAGS ${LAPACK_COMPILE_OPTIONS})
  check_cxx_source_compiles(
    "
    #include <blis.h>
    int main(int argc, char** argv) {
      return 0;
    }
    " MADNESS_CAN_INCLUDE_BLIS_H)
  cmake_pop_check_state()
endif()
