# - Try to find Libunwind
# Input variables:
#  ELEMENTAL_ROOT_DIR     - The Elemental install directory
#  ELEMENTAL_INCLUDE_DIR  - The Elemental include directory
#  ELEMENTAL_LIBRARY      - The Elemental library directory
# Output variables:
#  ELEMENTAL_FOUND        - System has Elemental
#  ELEMENTAL_INCLUDE_DIRS - The Elemental include directories
#  ELEMENTAL_LIBRARIES    - The libraries needed to use Elemental
#  ELEMENTAL_VERSION      - The version string for Elemental

include(FindPackageHandleStandardArgs)

if(NOT ELEMENTAL_FOUND)

  # Set default sarch paths for Elemental
  if(ELEMENTAL_ROOT_DIR)
    set(ELEMENTAL_INCLUDE_DIR ${ELEMENTAL_ROOT_DIR}/include CACHE PATH "The include directory for Elemental")
    set(ELEMENTAL_LIBRARY ${ELEMENTAL_ROOT_DIR}/lib CACHE PATH "The library directory for Elemental")
  endif()
  
  find_path(ELEMENTAL_INCLUDE_DIRS NAMES El.hpp elemental.hpp
      HINTS ${ELEMENTAL_INCLUDE_DIR})
  
  find_library(ELEMENTAL_LIBRARIES El elemental 
      HINTS ${ELEMENTAL_LIBRARY})
  if(ELEMENTAL_LIBRARIES)
    get_filename_component(ELEMENTAL_LIB_SEARCH_DIR "${ELEMENTAL_LIBRARIES}" PATH)
    foreach(_comp pmrrr lapack-addons kiss_fft metis)
      find_library(ELEMENTAL_${_comp}_LIBRARY ${_comp} 
          HINTS ${ELEMENTAL_LIB_SEARCH_DIR}
          NO_DEFAULT_PATH)
          
      if(ELEMENTAL_${_comp}_LIBRARY)
        set(Elemental_${_comp}_FOUND TRUE)
        set(ELEMENTAL_LIBRARIES ${ELEMENTAL_${_comp}_LIBRARY} ${ELEMENTAL_LIBRARIES})
      else()
        set(Elemental_${_comp}_FOUND FALSE)
      endif()
    endforeach()
  endif()
  
  # Get Elemental version
  if(ELEMENTAL_INCLUDE_DIRS)
    if(EXISTS ${ELEMENTAL_INCLUDE_DIRS}/El)
      file(READ "${ELEMENTAL_INCLUDE_DIRS}/El/config.h" _elemental_version_file)
    else()
      file(READ "${ELEMENTAL_INCLUDE_DIRS}/elemental/config.h" _elemental_version_file)
    endif()
    string(REGEX REPLACE ".*#define[ \t]+(EL|Elemental)_VERSION_MAJOR \\\"([0-9]+)\\\".*" "\\2"
            ELEMENTAL_VERSION_MAJOR "${_elemental_version_file}")
    string(REGEX REPLACE ".*#define[ \t]+(EL|Elemental)_VERSION_MINOR \\\"([0-9]+)(\\-dev|)\\\".*" "\\2"
            ELEMENTAL_VERSION_MINOR "${_elemental_version_file}")
    set(ELEMENTAL_VERSION "${ELEMENTAL_VERSION_MAJOR}.${ELEMENTAL_VERSION_MINOR}")
    unset(_elemental_version_file)
  endif()

  # handle the QUIETLY and REQUIRED arguments and set ELEMENTAL_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(Elemental
      FOUND_VAR ELEMENTAL_FOUND
      VERSION_VAR ELEMENTAL_VERSION
      HANDLE_COMPONENTS
      REQUIRED_VARS ELEMENTAL_LIBRARIES ELEMENTAL_INCLUDE_DIRS)

  mark_as_advanced(ELEMENTAL_INCLUDE_DIR ELEMENTAL_LIBRARY 
      ELEMENTAL_INCLUDE_DIRS ELEMENTAL_LIBRARIES)

endif()