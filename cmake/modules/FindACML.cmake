# - Try to find acml
# Input variables:
#  ACML_ROOT_DIR     - The acml install directory
#  ACML_LIBRARY      - The acml library directory
# Output variables:
#  ACML_FOUND        - System has acml
#  ACML_LIBRARIES    - The libraries needed to use acml

include(FindPackageHandleStandardArgs)

if(NOT ACML_FOUND)

  # Set default sarch paths for acml
  if(ACML_ROOT_DIR)
    if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set(ACML_LIBRARY ${ACML_ROOT_DIR}/lib64;${ACML_ROOT_DIR}/lib CACHE PATH "The library directory for acml")
    else()
      set(ACML_LIBRARY ${ACML_ROOT_DIR}/lib CACHE PATH "The library directory for acml")
    endif()
  endif()
  
  find_library(ACML_LIBRARIES acml 
      HINTS ${ACML_LIBRARY})

  # handle the QUIETLY and REQUIRED arguments and set ACML_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(ACML
      FOUND_VAR ACML_FOUND
      REQUIRED_VARS ACML_LIBRARIES)

  mark_as_advanced(ACML_INCLUDE_DIR ACML_LIBRARY 
      ACML_INCLUDE_DIRS ACML_LIBRARIES)

endif()