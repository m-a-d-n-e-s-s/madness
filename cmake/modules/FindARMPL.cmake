# - Try to find ARMPL (Arm Performance Libraries)
#
# This module defines:
#  ARMPL_FOUND - System has ARMPL
#  ARMPL_INCLUDE_DIRS - The ARMPL include directories
#  ARMPL_LIBRARIES - The libraries needed to use ARMPL
#  ARMPL_VERSION - The version string for ARMPL
#  ARMPL::ARMPL - Imported target for ARMPL
#
# Input variables:
#  ARMPL_ROOT_DIR, ARMPL_DIR, ARMPL_ROOT - The ARMPL install root directory
#  ARMPL_INCLUDE_DIR - The ARMPL include directory
#  ARMPL_LIBRARY_DIR - The ARMPL library directory
#  FORTRAN_INTEGER_SIZE - Integer size (4 for lp64, 8 for ilp64, default 4)
#  INTEGER4, INTEGER8 - Boolean indicators for integer size
#  ARMPL_THREADING - "seq" (default) or "omp"/"mp"
#  BLA_STATIC, ARMPL_STATIC - Prefer static libraries

include(FindPackageHandleStandardArgs)

if(NOT ARMPL_FOUND)

  # Collect candidate root directories
  set(_armpl_roots)
  if(ARMPL_ROOT_DIR)
    list(APPEND _armpl_roots "${ARMPL_ROOT_DIR}")
  endif()
  if(ARMPL_DIR)
    list(APPEND _armpl_roots "${ARMPL_DIR}")
  endif()
  if(ARMPL_ROOT)
    list(APPEND _armpl_roots "${ARMPL_ROOT}")
  endif()
  if(DEFINED ENV{ARMPL_DIR})
    list(APPEND _armpl_roots "$ENV{ARMPL_DIR}")
  endif()
  if(DEFINED ENV{ARMPL_ROOT})
    list(APPEND _armpl_roots "$ENV{ARMPL_ROOT}")
  endif()
  if(DEFINED ENV{ARMPL_ROOT_DIR})
    list(APPEND _armpl_roots "$ENV{ARMPL_ROOT_DIR}")
  endif()
  if(DEFINED ENV{ARMPL_HOME})
    list(APPEND _armpl_roots "$ENV{ARMPL_HOME}")
  endif()
  if(DEFINED ENV{ARM_DIR})
    list(APPEND _armpl_roots "$ENV{ARM_DIR}")
  endif()
  if(DEFINED ENV{ARM_ROOT})
    list(APPEND _armpl_roots "$ENV{ARM_ROOT}")
  endif()

  # Standard installation paths on Linux / AArch64
  list(APPEND _armpl_roots
    "/opt/arm/arm-performance-libraries"
    "/opt/arm"
    "/usr/local"
    "/usr"
  )
  file(GLOB _opt_arm_pl_dirs "/opt/arm/arm-performance-libraries*")
  file(GLOB _opt_armpl_dirs "/opt/arm/armpl*")
  list(APPEND _armpl_roots ${_opt_arm_pl_dirs} ${_opt_armpl_dirs})
  if(_armpl_roots)
    list(REMOVE_DUPLICATES _armpl_roots)
  endif()

  # Determine integer size
  if((DEFINED FORTRAN_INTEGER_SIZE AND FORTRAN_INTEGER_SIZE EQUAL 8) OR INTEGER8)
    set(_armpl_int "ilp64")
  else()
    set(_armpl_int "lp64")
  endif()

  # Determine threading model (default sequential for thread-safe concurrent tasks)
  if(ARMPL_THREADING MATCHES "omp|mp|openmp")
    set(_armpl_thread "_mp")
  else()
    set(_armpl_thread "")
  endif()

  # Candidate names for the main ARMPL library
  set(_armpl_lib_names
    "armpl_${_armpl_int}${_armpl_thread}"
    "armpl_${_armpl_int}"
    "armpl"
  )

  # Search for ARMPL include directory
  find_path(ARMPL_INCLUDE_DIR
    NAMES armpl.h
    HINTS ${ARMPL_INCLUDE_DIR} ${_armpl_roots}
    PATH_SUFFIXES
      include
      armpl_include
      include/armpl
  )

  # Search for ARMPL library directory and libraries
  set(_armpl_lib_hints)
  if(ARMPL_LIBRARY_DIR)
    list(APPEND _armpl_lib_hints "${ARMPL_LIBRARY_DIR}")
  endif()
  foreach(_root ${_armpl_roots})
    list(APPEND _armpl_lib_hints "${_root}/lib" "${_root}/lib64" "${_root}")
  endforeach()

  if(BLA_STATIC OR ARMPL_STATIC)
    set(_static_names)
    foreach(_name ${_armpl_lib_names})
      list(APPEND _static_names "${CMAKE_STATIC_LIBRARY_PREFIX}${_name}${CMAKE_STATIC_LIBRARY_SUFFIX}" "${_name}")
    endforeach()
    find_library(ARMPL_armpl_LIBRARY NAMES ${_static_names} HINTS ${_armpl_lib_hints})
    find_library(ARMPL_amath_LIBRARY NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}amath${CMAKE_STATIC_LIBRARY_SUFFIX}" amath amath_repro HINTS ${_armpl_lib_hints})
    find_library(ARMPL_astring_LIBRARY NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}astring${CMAKE_STATIC_LIBRARY_SUFFIX}" astring HINTS ${_armpl_lib_hints})
  else()
    find_library(ARMPL_armpl_LIBRARY NAMES ${_armpl_lib_names} HINTS ${_armpl_lib_hints})
    find_library(ARMPL_amath_LIBRARY NAMES amath amath_repro HINTS ${_armpl_lib_hints})
    find_library(ARMPL_astring_LIBRARY NAMES astring HINTS ${_armpl_lib_hints})
  endif()

  find_library(ARMPL_m_LIBRARY NAMES m)

  # Extract version if headers found
  if(ARMPL_INCLUDE_DIR AND EXISTS "${ARMPL_INCLUDE_DIR}/armpl.h")
    file(READ "${ARMPL_INCLUDE_DIR}/armpl.h" _armpl_h_content)
    if(_armpl_h_content MATCHES "Arm Performance Libraries version ([0-9]+(\\.[0-9]+)+)")
      set(ARMPL_VERSION "${CMAKE_MATCH_1}")
    endif()
  endif()

  # If version not found in header, check pkgconfig if available
  if(NOT ARMPL_VERSION)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      foreach(_root ${_armpl_roots})
        if(EXISTS "${_root}/lib/pkgconfig/armpl-${_armpl_int}-seq.pc")
          file(READ "${_root}/lib/pkgconfig/armpl-${_armpl_int}-seq.pc" _pc_content)
          if(_pc_content MATCHES "Version:[ \t]*([0-9]+(\\.[0-9]+)+)")
            set(ARMPL_VERSION "${CMAKE_MATCH_1}")
            break()
          endif()
        endif()
      endforeach()
    endif()
  endif()

  # Build list of libraries
  if(ARMPL_armpl_LIBRARY)
    set(ARMPL_LIBRARIES "${ARMPL_armpl_LIBRARY}")
    if(ARMPL_amath_LIBRARY)
      list(APPEND ARMPL_LIBRARIES "${ARMPL_amath_LIBRARY}")
    endif()
    if(ARMPL_astring_LIBRARY)
      list(APPEND ARMPL_LIBRARIES "${ARMPL_astring_LIBRARY}")
    endif()
    if(ARMPL_m_LIBRARY)
      list(APPEND ARMPL_LIBRARIES "${ARMPL_m_LIBRARY}")
    else()
      list(APPEND ARMPL_LIBRARIES "-lm")
    endif()
  endif()

  if(ARMPL_INCLUDE_DIR)
    set(ARMPL_INCLUDE_DIRS "${ARMPL_INCLUDE_DIR}")
  endif()

  find_package_handle_standard_args(ARMPL
    FOUND_VAR ARMPL_FOUND
    REQUIRED_VARS ARMPL_LIBRARIES ARMPL_INCLUDE_DIRS
    VERSION_VAR ARMPL_VERSION
  )

  if(ARMPL_FOUND)
    set(HAVE_ARMPL 1)
    if(NOT TARGET ARMPL::ARMPL)
      add_library(ARMPL::ARMPL UNKNOWN IMPORTED)
      set_target_properties(ARMPL::ARMPL PROPERTIES
        IMPORTED_LOCATION "${ARMPL_armpl_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${ARMPL_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${ARMPL_LIBRARIES}"
      )
    endif()
  endif()

  mark_as_advanced(ARMPL_INCLUDE_DIR ARMPL_armpl_LIBRARY ARMPL_amath_LIBRARY ARMPL_astring_LIBRARY)

endif()
