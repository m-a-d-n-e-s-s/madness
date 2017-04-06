# - Try to find Papi
# Input variables:
#  PAPI_ROOT_DIR     - The papi install directory
#  PAPI_INCLUDE_DIR  - The papi include directory
#  PAPI_LIBRARY      - The papi library directory
# Output variables:
#  PAPI_FOUND        - System has papi
#  PAPI_INCLUDE_DIRS - The papi include directories
#  PAPI_LIBRARIES    - The libraries needed to use papi

include(FindPackageHandleStandardArgs)
  
if(NOT PAPI_FOUND)

  # Set default sarch paths for papi
  if(PAPI_ROOT_DIR)
    set(PAPI_INCLUDE_DIR ${PAPI_ROOT_DIR}/include CACHE PATH "The include directory for papi")
    if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set(PAPI_LIBRARY ${PAPI_ROOT_DIR}/lib64;${PAPI_ROOT_DIR}/lib CACHE PATH "The library directory for papi")
    else()
      set(PAPI_LIBRARY ${PAPI_ROOT_DIR}/lib CACHE PATH "The library directory for papi")
    endif()
  endif()
  
  find_path(PAPI_INCLUDE_DIRS NAMES papi.h
      HINTS ${PAPI_INCLUDE_DIR})

  find_library(PAPI_papi_LIBRARIES papi 
      HINTS ${PAPI_LIBRARY})

  if(PAPI_INCLUDE_DIRS AND PAPI_papi_LIBRARY)
    add_library(papi UNKNOWN IMPORTED)
    set_target_properties(papi PROPERTIES
        IMPORTED_LOCATION "${PAPI_papi_LIBRARIES}"
        INTERFACE_INCLUDE_DIRECTORIES "${PAPI_INCLUDE_DIRS}")
  endif()
  
  set(PAPI_LIBRARIES papi)
  
  # handle the QUIETLY and REQUIRED arguments and set PAPI_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(Papi
      FOUND_VAR PAPI_FOUND
      REQUIRED_VARS PAPI_LIBRARIES PAPI_INCLUDE_DIRS)

  mark_as_advanced(PAPI_INCLUDE_DIR PAPI_LIBRARY 
      PAPI_INCLUDE_DIRS PAPI_LIBRARIES)

endif()