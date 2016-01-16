# - Try to find Intel TBB
# Input variables:
#  TBB_ROOT_DIR     - The TBB install directory
#  TBB_INCLUDE_DIR  - The TBB include directory
#  TBB_LIBRARY      - The TBB library directory
# Output variables:
#  TBB_FOUND        - System has TBB
#  TBB_INCLUDE_DIRS - The tbb include directories
#  TBB_LIBRARIES    - The libraries needed to use TBB
#  TBB_VERSION      - The version string for TBB

include(FindPackageHandleStandardArgs)

if(NOT TBB_FOUND)

  # Set default sarch paths for TBB
  if(NOT TBB_ROOT_DIR AND NOT DEFINED TBB_ROOT_DIR)
    if(EXISTS $ENV{TBBROOT})
      set(TBB_ROOT_DIR "$ENV{TBBROOT}")
    elseif(EXISTS /opt/intel/tbb)
      set(TBB_ROOT_DIR /opt/intel/tbb)
    endif()
  endif()
  if(TBB_ROOT_DIR)
    # NOTE: Will not overwrite user defined include and library directory variables
    set(TBB_INCLUDE_DIR ${TBB_ROOT_DIR}/include 
        CACHE PATH "The include directory for TBB")
    if(CMAKE_SYSTEM_NAME MATCHES "Darwin")
      set(TBB_LIBRARY ${TBB_ROOT_DIR}/lib/libc++;${TBB_ROOT_DIR}/lib 
          CACHE PATH "The library directory for TBB")
    elseif(CMAKE_SYSTEM_NAME MATCHES "Linux")
      set(TBB_LIBRARY ${TBB_ROOT_DIR}/lib/intel64/gcc4.4;${TBB_ROOT_DIR}/lib 
          CACHE PATH "The library directory for TBB")
    else()
      set(TBB_LIBRARY ${TBB_ROOT_DIR}/lib
          CACHE PATH "The library directory for TBB")
    endif()
  endif()
  
  if(CMAKE_BUILD_TYPE STREQUAL "Debug" OR CMAKE_BUILD_TYPE STREQUAL "RelWithDebInfo")
    set(TBB_USE_DEBUG TRUE)
  else()
    set(TBB_USE_DEBUG FALSE)
  endif()
  
  # Search for TBB include directory
  find_path(TBB_INCLUDE_DIRS NAMES tbb/tbb.h
      HINTS ${TBB_INCLUDE_DIR})
  
  # Search for TBB libraries
  find_library(TBB_tbb_LIBRARY tbb 
      HINTS ${TBB_LIBRARY})
  if(TBB_tbb_LIBRARY)
    get_filename_component(TBB_tbb_LIBRARY_DIR "${TBB_tbb_LIBRARY}" PATH)
    find_library(TBB_tbb_debug_LIBRARY tbb_debug 
        HINTS ${TBB_tbb_LIBRARY_DIR}
        NO_DEFAULT_PATH)
      
    foreach(_comp tbb_preview tbbmalloc tbbmalloc_proxy) 
      find_library(TBB_${_comp}_LIBRARY ${_comp} 
          HINTS ${TBB_tbb_LIBRARY_DIR}
          NO_DEFAULT_PATH)
      find_library(TBB_${_comp}_debug_LIBRARY ${_comp}_debug 
          HINTS ${TBB_tbb_LIBRARY_DIR}
          NO_DEFAULT_PATH)
    endforeach()
  endif()
  
  # Process TBB libaraies
  foreach(_lib tbb tbb_preview tbbmalloc tbbmalloc_proxy)
    # Set library found variables
    if(TBB_${_lib}_LIBRARY)
      set(TBB_${_lib}_FOUND TRUE)
    else()
      set(TBB_${_lib}_FOUND FALSE)
    endif()
    if(TBB_${_lib}_debug_LIBRARY)
      set(TBB_${_lib}_debug_FOUND TRUE)
    else()
      set(TBB_${_lib}_debug_FOUND FALSE)
    endif()

    # Set the build type TBB library variables
    if(_lib STREQUAL "tbb" OR ";${TBB_FIND_COMPONENTS};" MATCHES ";${_lib};")
      if(TBB_${_lib}_FOUND)
        set(TBB_LIBRARIES_RELEASE ${TBB_${_lib}_LIBRARY} ${TBB_LIBRARIES_RELEASE})
      endif()
      if(TBB_${_lib}_debug_FOUND)
        set(TBB_LIBRARIES_DEBUG ${TBB_${_lib}_debug_LIBRARY} ${TBB_LIBRARIES_DEBUG})
      endif()
    endif()
  endforeach()
  
  # Set the TBB_LIBRARIES variable
  if(TBB_USE_DEBUG AND TBB_LIBRARIES_DEBUG)
    set(TBB_LIBRARIES ${TBB_LIBRARIES_DEBUG})
  else()
    set(TBB_LIBRARIES ${TBB_LIBRARIES_RELEASE})
  endif()
  
  # Get TBB version
  if(TBB_INCLUDE_DIRS)
    file(READ "${TBB_INCLUDE_DIRS}/tbb/tbb_stddef.h" _tbb_version_file)
    string(REGEX REPLACE ".*#define TBB_VERSION_MAJOR ([0-9]+).*" "\\1"
            TBB_VERSION_MAJOR "${_tbb_version_file}")
    string(REGEX REPLACE ".*#define TBB_VERSION_MINOR ([0-9]+).*" "\\1"
            TBB_VERSION_MINOR "${_tbb_version_file}")
    string(REGEX REPLACE ".*#define TBB_INTERFACE_VERSION ([0-9]+).*" "\\1"
            TBB_INTERFACE_VERSION "${_tbb_version_file}")
    set(TBB_VERSION "${TBB_VERSION_MAJOR}.${TBB_VERSION_MINOR}")
    unset(_tbb_version_header)
  endif()

  # handle the QUIETLY and REQUIRED arguments and set TBB_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(TBB
      FOUND_VAR TBB_FOUND
      VERSION_VAR TBB_VERSION 
      REQUIRED_VARS TBB_LIBRARIES TBB_INCLUDE_DIRS
      HANDLE_COMPONENTS)

  if(TBB_LIBRARIES_DEBUG)
    set(TBB_COMPILE_FLAGS_DEBUG "-DTBB_USE_DEBUG=1")
  endif()

  mark_as_advanced(TBB_INCLUDE_DIR TBB_LIBRARY TBB_INCLUDE_DIRS TBB_LIBRARIES)

endif()