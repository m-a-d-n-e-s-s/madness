# - Try to find Libxc
# Input variables:
#  LIBXC_ROOT_DIR     - The libxc install directory
#  LIBXC_INCLUDE_DIR  - The libxc include directory
#  LIBXC_LIBRARY      - The libxc library directory
# Output variables:
#  LIBXC_FOUND        - System has libxc
#  LIBXC_INCLUDE_DIRS - The libxc include directories
#  LIBXC_LIBRARIES    - The libraries needed to use libxc
#  LIBXC_VERSION      - The version string for libxc

include(FindPackageHandleStandardArgs)

if(NOT LIBXC_FOUND)

  # Set default search paths for libxc
  if(LIBXC_ROOT_DIR)
    set(LIBXC_INCLUDE_DIR ${LIBXC_ROOT_DIR}/include CACHE PATH "The include directory for libxc")
    if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set(LIBXC_LIBRARY ${LIBXC_ROOT_DIR}/lib64;${LIBXC_ROOT_DIR}/lib CACHE PATH "The library directory for libxc")
    else()
      set(LIBXC_LIBRARY ${LIBXC_ROOT_DIR}/lib CACHE PATH "The library directory for libxc")
    endif()
  endif()
  
  find_path(LIBXC_INCLUDE_DIRS NAMES xc.h xc_funcs.h
      HINTS ${LIBXC_INCLUDE_DIR})
  
  find_library(LIBXC_LIBRARIES xc 
      HINTS ${LIBXC_LIBRARY})
  
  # Get libxc version
  if(LIBXC_INCLUDE_DIRS)
    file(READ "${LIBXC_INCLUDE_DIRS}/xc_version.h" _libxc_version_header)
    string(REGEX MATCH "define[ \t]+XC_VERSION[ \t]+\\\"([0-9\\.]+)\\\"" 
        LIBXC_VERSION "${_libxc_version_header}")
    string(REGEX MATCH "([0-9\\.]+)" LIBXC_VERSION "${LIBXC_VERSION}")
    unset(_libxc_version_header)
  endif()

  # handle the QUIETLY and REQUIRED arguments and set LIBXC_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(Libxc
      FOUND_VAR LIBXC_FOUND
      VERSION_VAR LIBXC_VERSION 
      REQUIRED_VARS LIBXC_LIBRARIES LIBXC_INCLUDE_DIRS)

  mark_as_advanced(LIBXC_INCLUDE_DIR LIBXC_LIBRARY 
      LIBXC_INCLUDE_DIRS LIBXC_LIBRARIES)

endif()
