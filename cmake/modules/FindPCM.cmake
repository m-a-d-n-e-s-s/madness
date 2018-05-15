# - Try to find PCM
# Input variables:
#  PCM_ROOT_DIR     - The pcm  install directory
#  PCM_INCLUDE_DIR  - The pcm  include directory | optional, else determined from PCM_ROOT_DIR, dont set PCM_ROOT_DIR if you want to use this
#  PCM_LIBRARY      - The pcm  library directory | optional, else determined from PCM_ROOT_DIR, ...
# Output variables:
#  PCM_FOUND        - System has pcm
#  PCM_INCLUDE_DIRS - The pcm include directories
#  PCM_LIBRARIES    - The libraries needed to use pcm
#  PCM_VERSION      - The version string for pcm | currently unused

include(FindPackageHandleStandardArgs)

if(NOT PCM_FOUND)

  # define include and library directories based on root directory
  if(PCM_ROOT_DIR) 
    set(PCM_INCLUDE_DIR ${PCM_ROOT_DIR}/include/ CACHE PATH "The include directory for PCM")
    if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set(PCM_LIBRARY ${PCM_ROOT_DIR}/lib;${PCM_ROOT_DIR}/lib64 CACHE PATH "Linker Flags for PCM Library")
    else()
      set(PCM_LIBRARY ${PCM_ROOT_DIR}/lib CACHE PATH "Linker Flags for PCM Library")
    endif()
  endif()
  
  find_path(PCM_INCLUDE_DIRS NAMES PCMSolver/pcmsolver.h 
      HINTS ${PCM_INCLUDE_DIR})
  
  find_library(PCM_LIBRARIES NAMES pcm 
      HINTS ${PCM_LIBRARY})
  
  # Get PCM version
#  if(PCM_INCLUDE_DIRS)
#    file(READ "${PCM_INCLUDE_DIRS}/GitInfo.hpp" _PCM_version_header)
#    string(REGEX MATCH "define[ \t]+PCM_VERSION[ \t]+\\\"([0-9\\.]+)\\\"" 
#        PCM_VERSION "${_PCM_version_header}")
#    string(REGEX MATCH "([0-9\\.]+)" PCM_VERSION "${PCM_VERSION}")
#    unset(_PCM_version_header)
#  endif()


  # handle the QUIETLY and REQUIRED arguments and set PCM_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(PCM
      FOUND_VAR PCM_FOUND
      VERSION_VAR PCM_VERSION 
      REQUIRED_VARS PCM_LIBRARIES PCM_INCLUDE_DIRS)

  mark_as_advanced(PCM_INCLUDE_DIR PCM_LIBRARY 
      PCM_INCLUDE_DIRS PCM_LIBRARIES)

endif()
