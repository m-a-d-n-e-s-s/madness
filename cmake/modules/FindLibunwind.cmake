# - Try to find Libunwind
# Input variables:
#  LIBUNWIND_DIR     - The libunwind install directory;
#                      if not set the LIBUNWIND_DIR environment variable will be used
# Output variables:
#  LIBUNWIND_FOUND        - System has libunwind
#  LIBUNWIND_INCLUDE_DIR  - The libunwind include directories
#  LIBUNWIND_LIBRARIES    - The libraries needed to use libunwind
#  LIBUNWIND_VERSION      - The version string for libunwind

include(FindPackageHandleStandardArgs)
  
if(NOT DEFINED LIBUNWIND_FOUND)

  # if not set already, set LIBUNWIND_DIR from environment
  if (DEFINED ENV{LIBUNWIND_DIR} AND NOT DEFINED LIBUNWIND_DIR)
    set(LIBUNWIND_DIR $ENV{LIBUNWIND_DIR})
  endif()

  # Set default search paths for libunwind
  if(LIBUNWIND_DIR)
    set(LIBUNWIND_INCLUDE_DIR_HINT ${LIBUNWIND_DIR}/include CACHE PATH "The include directory for libunwind")
    if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set(LIBUNWIND_LIBRARY_DIR_HINT ${LIBUNWIND_DIR}/lib64;${LIBUNWIND_DIR}/lib CACHE PATH "The library directory for libunwind")
    else()
      set(LIBUNWIND_LIBRARY_DIR_HINT ${LIBUNWIND_DIR}/lib CACHE PATH "The library directory for libunwind")
    endif()
  endif()

  find_path(LIBUNWIND_INCLUDE_DIR NAMES libunwind.h
      HINTS ${LIBUNWIND_INCLUDE_DIR_HINT})
  
  find_library(LIBUNWIND_LIBRARIES unwind 
      HINTS ${LIBUNWIND_LIBRARY_DIR_HINT})
  
  # Get libunwind version
  if(EXISTS "${LIBUNWIND_INCLUDE_DIR}/libunwind-common.h")
    file(READ "${LIBUNWIND_INCLUDE_DIR}/libunwind-common.h" _libunwind_version_header)
    string(REGEX REPLACE ".*define[ \t]+UNW_VERSION_MAJOR[ \t]+([0-9]+).*" "\\1" 
        LIBUNWIND_MAJOR_VERSION "${_libunwind_version_header}")
    string(REGEX REPLACE ".*define[ \t]+UNW_VERSION_MINOR[ \t]+([0-9]+).*" "\\1"
        LIBUNWIND_MINOR_VERSION "${_libunwind_version_header}")
    string(REGEX REPLACE ".*define[ \t]+UNW_VERSION_EXTRA[ \t]+([0-9]+).*" "\\1"
        LIBUNWIND_MICRO_VERSION "${_libunwind_version_header}")
    set(LIBUNWIND_VERSION "${LIBUNWIND_MAJOR_VERSION}.${LIBUNWIND_MINOR_VERSION}.${LIBUNWIND_MICRO_VERSION}")
    unset(_libunwind_version_header)
  endif()

  # handle the QUIETLY and REQUIRED arguments and set LIBUNWIND_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(Libunwind
      FOUND_VAR LIBUNWIND_FOUND
      VERSION_VAR LIBUNWIND_VERSION 
      REQUIRED_VARS LIBUNWIND_LIBRARIES LIBUNWIND_INCLUDE_DIR)

  mark_as_advanced(LIBUNWIND_INCLUDE_DIR_HINT LIBUNWIND_LIBRARY_DIR_HINT
      LIBUNWIND_INCLUDE_DIR LIBUNWIND_LIBRARIES)

endif()