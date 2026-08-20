# - Try to find simple-dftd3 (s-dftd3), Grimme's D3 dispersion correction
#
# Only needed for installations that do not ship a CMake package config; the
# conda-forge binary is built with meson and provides pkg-config only.  A
# source build of simple-dftd3 exports `s-dftd3::s-dftd3` and is picked up by
# find_package(s-dftd3 CONFIG) in external/dftd3.cmake before this module runs.
#
# Search order for hints: DFTD3_ROOT_DIR, then pkg-config, then $CONDA_PREFIX.
# The conda fallback exists because `conda install -c conda-forge simple-dftd3`
# is the documented way to get the library and leaves nothing on the default
# CMake or pkg-config search paths -- without it, configuring inside an
# activated environment reports the package as missing for no visible reason.
#
# Input variables:
#  DFTD3_ROOT_DIR     - The s-dftd3 install directory (e.g. $CONDA_PREFIX)
#  DFTD3_INCLUDE_DIR  - The s-dftd3 include directory | optional, else determined from DFTD3_ROOT_DIR
#  DFTD3_LIBRARY      - The s-dftd3 library directory | optional, else determined from DFTD3_ROOT_DIR
# Output variables:
#  DFTD3_FOUND        - System has s-dftd3
#  DFTD3_INCLUDE_DIRS - The s-dftd3 include directories
#  DFTD3_LIBRARIES    - The libraries needed to use s-dftd3
#  DFTD3_VERSION      - The version string for s-dftd3, if pkg-config supplied one

include(FindPackageHandleStandardArgs)

if(NOT DFTD3_FOUND)

  # define include and library directories based on root directory
  if(DFTD3_ROOT_DIR)
    set(DFTD3_INCLUDE_DIR ${DFTD3_ROOT_DIR}/include CACHE PATH "The include directory for s-dftd3")
    if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set(DFTD3_LIBRARY ${DFTD3_ROOT_DIR}/lib;${DFTD3_ROOT_DIR}/lib64 CACHE PATH "Linker flags for the s-dftd3 library")
    else()
      set(DFTD3_LIBRARY ${DFTD3_ROOT_DIR}/lib CACHE PATH "Linker flags for the s-dftd3 library")
    endif()
  endif()

  # pkg-config only supplies hints here -- an explicit DFTD3_ROOT_DIR wins, and a
  # missing/broken .pc file (s-dftd3.pc declares `Requires.private: mctc-lib`, so
  # it fails outright when mctc-lib.pc is absent) must not make the find fail.
  if(NOT DFTD3_ROOT_DIR)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      # An activated conda environment does not put its .pc files on
      # PKG_CONFIG_PATH, so add them for the duration of this query.
      set(_dftd3_saved_pkg_config_path "$ENV{PKG_CONFIG_PATH}")
      if(DEFINED ENV{CONDA_PREFIX})
        set(ENV{PKG_CONFIG_PATH} "$ENV{CONDA_PREFIX}/lib/pkgconfig:$ENV{PKG_CONFIG_PATH}")
      endif()
      pkg_check_modules(PC_SDFTD3 QUIET s-dftd3)
      set(ENV{PKG_CONFIG_PATH} "${_dftd3_saved_pkg_config_path}")
      unset(_dftd3_saved_pkg_config_path)
      if(PC_SDFTD3_VERSION)
        set(DFTD3_VERSION ${PC_SDFTD3_VERSION})
      endif()
    endif()
  endif()

  # last-resort hint: the active conda environment
  if(NOT DFTD3_ROOT_DIR AND DEFINED ENV{CONDA_PREFIX})
    set(_dftd3_conda_include "$ENV{CONDA_PREFIX}/include")
    set(_dftd3_conda_lib "$ENV{CONDA_PREFIX}/lib")
  endif()

  # s-dftd3.h, not dftd3.h: the latter is a one-line wrapper around the former
  # and shares its name with Grimme's older standalone dftd3 program.
  find_path(DFTD3_INCLUDE_DIRS NAMES s-dftd3.h
      HINTS ${DFTD3_INCLUDE_DIR} ${PC_SDFTD3_INCLUDE_DIRS} ${_dftd3_conda_include})

  find_library(DFTD3_LIBRARIES NAMES s-dftd3
      HINTS ${DFTD3_LIBRARY} ${PC_SDFTD3_LIBRARY_DIRS} ${_dftd3_conda_lib})

  unset(_dftd3_conda_include)
  unset(_dftd3_conda_lib)

  # handle the QUIETLY and REQUIRED arguments and set DFTD3_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(DFTD3
      FOUND_VAR DFTD3_FOUND
      VERSION_VAR DFTD3_VERSION
      REQUIRED_VARS DFTD3_LIBRARIES DFTD3_INCLUDE_DIRS)

  if(NOT DFTD3_FOUND)
    message(STATUS "simple-dftd3 not found; the `dispersion` keyword will be unavailable. "
                   "Install it (conda install -c conda-forge simple-dftd3) and point "
                   "-DDFTD3_ROOT_DIR at the prefix, or configure -DENABLE_DFTD3=OFF to stop looking.")
  endif()

  mark_as_advanced(DFTD3_INCLUDE_DIR DFTD3_LIBRARY
      DFTD3_INCLUDE_DIRS DFTD3_LIBRARIES)

endif()
