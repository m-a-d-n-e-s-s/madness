# - Find or fetch the LIBXSMM library
# This module defines the following targets:
#   libxsmm::libxsmm  - The LIBXSMM target (or alias to xsmm)
#   xsmm              - The LIBXSMM library target
#
# Variables defined:
#   MADNESS_HAS_LIBXSMM - TRUE if LIBXSMM was found or fetched

if (NOT TARGET libxsmm::libxsmm AND NOT TARGET xsmm)
  find_package(libxsmm QUIET CONFIG)
  if (TARGET libxsmm::libxsmm)
    message(STATUS "Found LIBXSMM CONFIG at ${libxsmm_DIR}")
  endif ()
endif ()

if (NOT TARGET libxsmm::libxsmm AND NOT TARGET xsmm)
  cmake_minimum_required(VERSION 3.14.0)
  include(FetchContent)

  find_package(Python3 COMPONENTS Interpreter QUIET)
  if (NOT Python3_EXECUTABLE AND TARGET Python::Interpreter)
    get_target_property(Python3_EXECUTABLE Python::Interpreter LOCATION)
  endif ()
  if (NOT Python3_EXECUTABLE AND Python_EXECUTABLE)
    set(Python3_EXECUTABLE "${Python_EXECUTABLE}")
  endif ()
  if (NOT Python3_EXECUTABLE)
    set(Python3_EXECUTABLE "/usr/bin/python3")
  endif ()

  set(CMAKE_Fortran_COMPILER NOTFOUND CACHE FILEPATH "Fortran compiler" FORCE)
  set(LIBXSMM_FORTRAN OFF CACHE BOOL "Disable Fortran support in LIBXSMM" FORCE)

  set(MADNESS_TRACKED_LIBXSMM_TAG "main" CACHE STRING "The tag/branch of LIBXSMM repository to track")

  FetchContent_Declare(
    libxsmm
    GIT_REPOSITORY https://github.com/libxsmm/libxsmm.git
    GIT_TAG        ${MADNESS_TRACKED_LIBXSMM_TAG}
    GIT_SHALLOW    TRUE
  )
  FetchContent_MakeAvailable(libxsmm)

  # Export libxsmm targets for build tree so MADNESS can export targets that depend on it
  export(EXPORT libxsmm-targets FILE "${PROJECT_BINARY_DIR}/libxsmm-targets.cmake")
endif ()

if (TARGET xsmm AND NOT TARGET libxsmm::libxsmm)
  add_library(libxsmm::libxsmm ALIAS xsmm)
endif ()

if (TARGET libxsmm::libxsmm OR TARGET xsmm)
  set(MADNESS_HAS_LIBXSMM ON CACHE BOOL "MADNESS has access to LIBXSMM")
  set(HAVE_MTXMQ ON CACHE BOOL "MADNESS has mTxmq enabled" FORCE)
else ()
  message(FATAL_ERROR "FindOrFetchLibxsmm could not make libxsmm::libxsmm target available")
endif ()
