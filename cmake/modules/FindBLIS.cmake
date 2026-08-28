# - Try to find BLIS (BLAS-like Library Instantiation Software)
#
# This module defines:
#  BLIS_FOUND - System has BLIS
#  BLIS_INCLUDE_DIRS - The BLIS include directories
#  BLIS_LIBRARIES - The libraries needed to use BLIS
#  BLIS_VERSION - The version string for BLIS
#  BLIS::BLIS - Imported target for BLIS
#
# Input variables:
#  BLIS_ROOT_DIR, BLIS_DIR, BLIS_ROOT - The BLIS install root directory
#  BLIS_INCLUDE_DIR - The BLIS include directory
#  BLIS_LIBRARY_DIR - The BLIS library directory
#  BLIS_SERIAL_ONLY - Require sequential/serial BLIS (default TRUE)
#  BLA_STATIC, BLIS_STATIC - Prefer static libraries

include(FindPackageHandleStandardArgs)

if(NOT DEFINED BLIS_SERIAL_ONLY)
  set(BLIS_SERIAL_ONLY TRUE)
endif()

if(NOT BLIS_FOUND)

  # Collect candidate root directories
  set(_blis_roots)
  if(BLIS_ROOT_DIR)
    list(APPEND _blis_roots "${BLIS_ROOT_DIR}")
  endif()
  if(BLIS_DIR)
    list(APPEND _blis_roots "${BLIS_DIR}")
  endif()
  if(BLIS_ROOT)
    list(APPEND _blis_roots "${BLIS_ROOT}")
  endif()
  if(DEFINED ENV{BLIS_DIR})
    list(APPEND _blis_roots "$ENV{BLIS_DIR}")
  endif()
  if(DEFINED ENV{BLIS_ROOT})
    list(APPEND _blis_roots "$ENV{BLIS_ROOT}")
  endif()
  if(DEFINED ENV{BLIS_ROOT_DIR})
    list(APPEND _blis_roots "$ENV{BLIS_ROOT_DIR}")
  endif()
  if(DEFINED ENV{BLIS_HOME})
    list(APPEND _blis_roots "$ENV{BLIS_HOME}")
  endif()
  if(DEFINED ENV{AOCL_ROOT})
    list(APPEND _blis_roots "$ENV{AOCL_ROOT}")
  endif()
  if(DEFINED ENV{AOCL_DIR})
    list(APPEND _blis_roots "$ENV{AOCL_DIR}")
  endif()
  if(DEFINED ENV{AMD_ROOT})
    list(APPEND _blis_roots "$ENV{AMD_ROOT}")
  endif()

  # System paths including Debian/Ubuntu multiarch directories, Fedora/RHEL lib64, and AOCL paths on x86
  list(APPEND _blis_roots
    ${CMAKE_PREFIX_PATH}
    "/usr/include/${CMAKE_LIBRARY_ARCHITECTURE}"
    "/usr/lib/${CMAKE_LIBRARY_ARCHITECTURE}"
    "/usr/include/aarch64-linux-gnu"
    "/usr/lib/aarch64-linux-gnu"
    "/usr/include/x86_64-linux-gnu"
    "/usr/lib/x86_64-linux-gnu"
    "/usr/local/lib64"
    "/usr/local/lib"
    "/usr/local"
    "/usr/lib64"
    "/usr/lib"
    "/usr"
    "/opt/blis"
    "/opt/local"
    "/opt/AMD/aocl"
    "/opt/aocl"
  )
  file(GLOB _aocl_dirs "/opt/AMD/aocl*" "/opt/aocl*")
  list(APPEND _blis_roots ${_aocl_dirs})
  if(_blis_roots)
    list(REMOVE_DUPLICATES _blis_roots)
  endif()

  # Search for BLIS include directory
  # If BLIS_SERIAL_ONLY is set, search serial-specific subdirectories first to avoid
  # picking up symlinks in /usr/include/aarch64-linux-gnu/ or /usr/include/ that point to blis-openmp
  if(BLIS_SERIAL_ONLY)
    find_path(BLIS_INCLUDE_DIR
      NAMES blis.h
      HINTS ${BLIS_INCLUDE_DIR} ${_blis_roots}
      PATH_SUFFIXES
        blis-serial
        include/blis-serial
        blis-serial/include
        blis
        include/blis
        include
      NO_DEFAULT_PATH
    )
    if(NOT BLIS_INCLUDE_DIR)
      find_path(BLIS_INCLUDE_DIR
        NAMES blis.h
        PATH_SUFFIXES
          blis-serial
          include/blis-serial
          blis
          include/blis
          include
      )
    endif()
  else()
    find_path(BLIS_INCLUDE_DIR
      NAMES blis.h
      HINTS ${BLIS_INCLUDE_DIR} ${_blis_roots}
      PATH_SUFFIXES
        blis
        blis-serial
        blis-openmp
        bliso
        blisp
        include/blis
        include
    )
  endif()

  # Check that include path is actually serial if BLIS_SERIAL_ONLY is ON
  if(BLIS_SERIAL_ONLY AND BLIS_INCLUDE_DIR)
    get_filename_component(_blis_inc_real "${BLIS_INCLUDE_DIR}/blis.h" REALPATH)
    if(_blis_inc_real MATCHES "openmp|pthread|bliso|blisp" AND NOT _blis_inc_real MATCHES "serial")
      message(STATUS "Rejecting non-serial BLIS include directory: ${BLIS_INCLUDE_DIR} (resolves to ${_blis_inc_real})")
      unset(BLIS_INCLUDE_DIR CACHE)
      unset(BLIS_INCLUDE_DIR)
    endif()
  endif()

  # Search for BLIS library
  set(_blis_lib_hints)
  if(BLIS_LIBRARY_DIR)
    list(APPEND _blis_lib_hints "${BLIS_LIBRARY_DIR}")
  endif()
  foreach(_root ${_blis_roots})
    list(APPEND _blis_lib_hints
      "${_root}/blis-serial" "${_root}/lib/blis-serial" "${_root}/lib64/blis-serial"
      "${_root}/lib64" "${_root}/lib" "${_root}"
    )
  endforeach()

  if(BLIS_SERIAL_ONLY)
    if(BLA_STATIC OR BLIS_STATIC)
      find_library(BLIS_LIBRARY
        NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}blis${CMAKE_STATIC_LIBRARY_SUFFIX}"
              "${CMAKE_STATIC_LIBRARY_PREFIX}blis-serial${CMAKE_STATIC_LIBRARY_SUFFIX}"
              blis blis-serial
        HINTS ${_blis_lib_hints}
        PATH_SUFFIXES blis-serial lib/blis-serial lib64/blis-serial
        NO_DEFAULT_PATH
      )
      if(NOT BLIS_LIBRARY)
        find_library(BLIS_LIBRARY
          NAMES "${CMAKE_STATIC_LIBRARY_PREFIX}blis${CMAKE_STATIC_LIBRARY_SUFFIX}"
                "${CMAKE_STATIC_LIBRARY_PREFIX}blis-serial${CMAKE_STATIC_LIBRARY_SUFFIX}"
                blis blis-serial
          HINTS ${_blis_lib_hints}
          PATH_SUFFIXES blis-serial lib/blis-serial lib64/blis-serial lib64 lib
        )
      endif()
    else()
      find_library(BLIS_LIBRARY
        NAMES blis blis-serial
        HINTS ${_blis_lib_hints}
        PATH_SUFFIXES blis-serial lib/blis-serial lib64/blis-serial
        NO_DEFAULT_PATH
      )
      if(NOT BLIS_LIBRARY)
        find_library(BLIS_LIBRARY
          NAMES blis blis-serial
          HINTS ${_blis_lib_hints}
          PATH_SUFFIXES blis-serial lib/blis-serial lib64/blis-serial lib64 lib
        )
      endif()
    endif()
  else()
    find_library(BLIS_LIBRARY
      NAMES blis blis-serial blis-openmp blis-threads bliso blisp blis-mt
      HINTS ${_blis_lib_hints}
      PATH_SUFFIXES lib64 lib blis-serial blis-openmp
    )
  endif()

  # Check that library is actually serial if BLIS_SERIAL_ONLY is ON
  if(BLIS_SERIAL_ONLY AND BLIS_LIBRARY)
    get_filename_component(_blis_lib_real "${BLIS_LIBRARY}" REALPATH)
    if(_blis_lib_real MATCHES "(openmp|pthread|bliso|blisp|-mt|_mt)" AND NOT _blis_lib_real MATCHES "serial")
      message(STATUS "Rejecting non-serial BLIS library: ${BLIS_LIBRARY} (resolves to ${_blis_lib_real})")
      unset(BLIS_LIBRARY CACHE)
      unset(BLIS_LIBRARY)
    endif()
  endif()

  # Runtime/compile check to ensure BLIS is single-threaded when BLIS_SERIAL_ONLY is ON
  if(BLIS_SERIAL_ONLY AND BLIS_LIBRARY AND BLIS_INCLUDE_DIR)
    include(CheckCSourceRuns)
    include(CMakePushCheckState)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_LIBRARIES "${BLIS_LIBRARY}")
    find_library(_blis_m_check NAMES m)
    if(_blis_m_check)
      list(APPEND CMAKE_REQUIRED_LIBRARIES "${_blis_m_check}")
    endif()
    find_package(Threads QUIET)
    if(TARGET Threads::Threads)
      list(APPEND CMAKE_REQUIRED_LIBRARIES Threads::Threads)
    endif()
    set(CMAKE_REQUIRED_INCLUDES "${BLIS_INCLUDE_DIR}")
    if(NOT CMAKE_CROSSCOMPILING)
      check_c_source_runs(
        "
        #include <blis.h>
        int main(void) {
          return (bli_info_get_enable_threading() == 0) ? 0 : 1;
        }
        " _BLIS_RUN_CHECK_SERIAL)
      if(NOT _BLIS_RUN_CHECK_SERIAL)
        message(STATUS "Rejecting multithreaded BLIS library: ${BLIS_LIBRARY} (bli_info_get_enable_threading() != 0)")
        unset(BLIS_LIBRARY CACHE)
        unset(BLIS_LIBRARY)
        unset(_BLIS_RUN_CHECK_SERIAL CACHE)
      endif()
    endif()
    cmake_pop_check_state()
  endif()

  # Extract version if header found
  if(BLIS_INCLUDE_DIR AND EXISTS "${BLIS_INCLUDE_DIR}/blis.h")
    file(READ "${BLIS_INCLUDE_DIR}/blis.h" _blis_h_content)
    if(_blis_h_content MATCHES "BLIS_VERSION_STRING[ \t]+\"([^\"]+)\"")
      set(BLIS_VERSION "${CMAKE_MATCH_1}")
    endif()
  endif()

  if(NOT BLIS_VERSION)
    find_package(PkgConfig QUIET)
    if(PKG_CONFIG_FOUND)
      foreach(_root ${_blis_roots})
        if(EXISTS "${_root}/blis-serial/pkgconfig/blis.pc" OR EXISTS "${_root}/lib/pkgconfig/blis.pc")
          file(GLOB _blis_pc_files "${_root}/blis-serial/pkgconfig/blis.pc" "${_root}/lib/pkgconfig/blis.pc")
          foreach(_pc ${_blis_pc_files})
            file(READ "${_pc}" _pc_content)
            if(_pc_content MATCHES "Version:[ \t]*([0-9]+(\\.[0-9]+)+)")
              set(BLIS_VERSION "${CMAKE_MATCH_1}")
              break()
            endif()
          endforeach()
        endif()
      endforeach()
    endif()
  endif()

  if(BLIS_LIBRARY)
    set(BLIS_LIBRARIES "${BLIS_LIBRARY}")
    find_library(BLIS_m_LIBRARY NAMES m)
    if(BLIS_m_LIBRARY)
      list(APPEND BLIS_LIBRARIES "${BLIS_m_LIBRARY}")
    else()
      list(APPEND BLIS_LIBRARIES "-lm")
    endif()
    # Also add pthread if needed for thread safety
    find_package(Threads QUIET)
    if(TARGET Threads::Threads)
      list(APPEND BLIS_LIBRARIES Threads::Threads)
    endif()
  endif()

  if(BLIS_INCLUDE_DIR)
    set(BLIS_INCLUDE_DIRS "${BLIS_INCLUDE_DIR}")
  endif()

  find_package_handle_standard_args(BLIS
    FOUND_VAR BLIS_FOUND
    REQUIRED_VARS BLIS_LIBRARIES BLIS_INCLUDE_DIRS
    VERSION_VAR BLIS_VERSION
  )

  if(BLIS_FOUND)
    set(HAVE_BLIS 1)
    if(NOT TARGET BLIS::BLIS)
      add_library(BLIS::BLIS UNKNOWN IMPORTED)
      set_target_properties(BLIS::BLIS PROPERTIES
        IMPORTED_LOCATION "${BLIS_LIBRARY}"
        INTERFACE_INCLUDE_DIRECTORIES "${BLIS_INCLUDE_DIRS}"
        INTERFACE_LINK_LIBRARIES "${BLIS_LIBRARIES}"
      )
    endif()
  endif()

  mark_as_advanced(BLIS_INCLUDE_DIR BLIS_LIBRARY)

endif()
