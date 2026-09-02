# - Try to find NVPL (NVIDIA Performance Libraries) BLAS & LAPACK for ARM
#
# This module defines:
#  NVPL_FOUND - System has NVPL
#  NVPL_INCLUDE_DIRS - The NVPL include directories
#  NVPL_LIBRARIES - The libraries needed to use NVPL
#  NVPL_VERSION - The version string for NVPL
#  NVPL::NVPL - Imported target for NVPL
#
# Input variables:
#  NVPL_ROOT_DIR, NVPL_DIR, NVPL_ROOT - The NVPL install root directory
#  NVPL_INCLUDE_DIR - The NVPL include directory
#  NVPL_LIBRARY_DIR - The NVPL library directory
#  FORTRAN_INTEGER_SIZE - Integer size (4 for lp64, 8 for ilp64, default 4)
#  INTEGER4, INTEGER8 - Boolean indicators for integer size
#  NVPL_THREADING - "seq" (default) or "omp"/"mp"
#  BLA_STATIC, NVPL_STATIC - Prefer static libraries

include(FindPackageHandleStandardArgs)

if(NOT NVPL_FOUND)

  # Collect candidate root directories
  set(_nvpl_roots)
  if(NVPL_ROOT_DIR)
    list(APPEND _nvpl_roots "${NVPL_ROOT_DIR}")
  endif()
  if(NVPL_DIR)
    list(APPEND _nvpl_roots "${NVPL_DIR}")
  endif()
  if(NVPL_ROOT)
    list(APPEND _nvpl_roots "${NVPL_ROOT}")
  endif()
  if(DEFINED ENV{NVPL_DIR})
    list(APPEND _nvpl_roots "$ENV{NVPL_DIR}")
  endif()
  if(DEFINED ENV{NVPL_ROOT})
    list(APPEND _nvpl_roots "$ENV{NVPL_ROOT}")
  endif()
  if(DEFINED ENV{NVPL_ROOT_DIR})
    list(APPEND _nvpl_roots "$ENV{NVPL_ROOT_DIR}")
  endif()
  if(DEFINED ENV{NVPL_HOME})
    list(APPEND _nvpl_roots "$ENV{NVPL_HOME}")
  endif()
  if(DEFINED ENV{NVHPC_ROOT})
    list(APPEND _nvpl_roots "$ENV{NVHPC_ROOT}")
  endif()
  if(DEFINED ENV{HPCX_DIR})
    list(APPEND _nvpl_roots "$ENV{HPCX_DIR}")
  endif()

  # Standard installation paths on Linux / AArch64
  list(APPEND _nvpl_roots
    "/opt/nvidia/nvpl"
    "/opt/nvidia/hpc_sdk"
    "/opt/nvidia"
    "/usr/local/cuda"
    "/usr/local"
    "/usr/lib/aarch64-linux-gnu"
    "/usr/lib64"
    "/usr/lib"
    "/usr"
  )
  file(GLOB _opt_nvpl_dirs "/opt/nvidia/nvpl*")
  file(GLOB _opt_nvhpc_dirs "/opt/nvidia/hpc_sdk/Linux_aarch64/*/math_libs")
  list(APPEND _nvpl_roots ${_opt_nvpl_dirs} ${_opt_nvhpc_dirs})
  if(_nvpl_roots)
    list(REMOVE_DUPLICATES _nvpl_roots)
  endif()

  # Determine integer size
  if((DEFINED FORTRAN_INTEGER_SIZE AND FORTRAN_INTEGER_SIZE EQUAL 8) OR INTEGER8)
    set(_nvpl_int "ilp64")
  else()
    set(_nvpl_int "lp64")
  endif()

  # Determine threading model (default sequential for thread-safe concurrent tasks)
  if(NVPL_THREADING MATCHES "omp|mp|openmp")
    set(_nvpl_thread "omp")
  else()
    set(_nvpl_thread "seq")
  endif()

  # Candidate names for NVPL BLAS & LAPACK libraries
  set(_nvpl_blas_names
    "nvpl_blas_${_nvpl_int}_${_nvpl_thread}"
    "nvpl_blas_core_${_nvpl_thread}"
    "nvpl_blas_${_nvpl_thread}"
    "nvpl_blas_${_nvpl_int}"
    "nvpl_blas"
  )
  set(_nvpl_lapack_names
    "nvpl_lapack_${_nvpl_int}_${_nvpl_thread}"
    "nvpl_lapack_core_${_nvpl_thread}"
    "nvpl_lapack_${_nvpl_thread}"
    "nvpl_lapack_${_nvpl_int}"
    "nvpl_lapack"
  )

  # Search for NVPL include directory
  find_path(NVPL_INCLUDE_DIR
    NAMES nvpl_blas.h nvpl_lapack.h
    HINTS ${NVPL_INCLUDE_DIR} ${_nvpl_roots}
    PATH_SUFFIXES
      include
      nvpl_include
      include/nvpl
      math_libs/include
  )

  # Search for NVPL library directory and libraries
  set(_nvpl_lib_hints)
  if(NVPL_LIBRARY_DIR)
    list(APPEND _nvpl_lib_hints "${NVPL_LIBRARY_DIR}")
  endif()
  list(APPEND _nvpl_lib_hints ${_nvpl_roots})

  set(_nvpl_suffixes
    lib
    lib64
    lib/aarch64-linux-gnu
    math_libs/lib64
    math_libs/lib
  )

  find_library(NVPL_BLAS_LIBRARY
    NAMES ${_nvpl_blas_names}
    HINTS ${_nvpl_lib_hints}
    PATH_SUFFIXES ${_nvpl_suffixes}
  )

  find_library(NVPL_LAPACK_LIBRARY
    NAMES ${_nvpl_lapack_names}
    HINTS ${_nvpl_lib_hints}
    PATH_SUFFIXES ${_nvpl_suffixes}
  )

  if(NVPL_BLAS_LIBRARY AND NVPL_LAPACK_LIBRARY)
    set(NVPL_LIBRARIES ${NVPL_LAPACK_LIBRARY} ${NVPL_BLAS_LIBRARY})
    set(NVPL_INCLUDE_DIRS ${NVPL_INCLUDE_DIR})

    find_package_handle_standard_args(NVPL
      REQUIRED_VARS NVPL_LIBRARIES
    )

    if(NVPL_FOUND AND NOT TARGET NVPL::NVPL)
      add_library(NVPL::NVPL INTERFACE IMPORTED)
      set_target_properties(NVPL::NVPL PROPERTIES
        INTERFACE_LINK_LIBRARIES "${NVPL_LIBRARIES}"
      )
      if(NVPL_INCLUDE_DIRS)
        set_target_properties(NVPL::NVPL PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${NVPL_INCLUDE_DIRS}"
        )
      endif()
    endif()
  endif()

endif()
