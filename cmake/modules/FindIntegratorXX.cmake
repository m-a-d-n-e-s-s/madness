# - Try to find INTEGRATORXX
# Input variables:
#  INTEGRATORXX_ROOT_DIR     - The INTEGRATORXX  install directory
#  INTEGRATORXX_INCLUDE_DIR  - The INTEGRATORXX  include directory | optional, else determined from INTEGRATORXX_ROOT_DIR, dont set INTEGRATORXX_ROOT_DIR if you want to use this
#  INTEGRATORXX_LIBRARY      - The INTEGRATORXX  library directory | optional, else determined from INTEGRATORXX_ROOT_DIR, ...
# Output variables:
#  INTEGRATORXX_FOUND        - System has INTEGRATORXX
#  INTEGRATORXX_INCLUDE_DIRS - The INTEGRATORXX include directories
#  INTEGRATORXX_LIBRARIES    - The libraries needed to use INTEGRATORXX
#  INTEGRATORXX_VERSION      - The version string for INTEGRATORXX | currently unused

include(FindPackageHandleStandardArgs)

if(NOT INTEGRATORXX_FOUND)

  # define include and library directories based on root directory
  if(INTEGRATORXX_ROOT_DIR)
    set(INTEGRATORXX_INCLUDE_DIR ${INTEGRATORXX_ROOT_DIR}/include/ CACHE PATH "The include directory for INTEGRATORXX")
    if(CMAKE_SIZEOF_VOID_P EQUAL 8 AND CMAKE_SYSTEM_NAME STREQUAL "Linux")
      set(INTEGRATORXX_LIBRARY ${INTEGRATORXX_ROOT_DIR}/lib;${INTEGRATORXX_ROOT_DIR}/lib64 CACHE PATH "Linker Flags for IntegratorXX Library")
    else()
      set(INTEGRATORXX_LIBRARY ${INTEGRATORXX_ROOT_DIR}/lib CACHE PATH "Linker Flags for IntegratorXX Library")
    endif()
  endif()
  
  find_path(INTEGRATORXX_INCLUDE_DIRS NAMES integratorxx/generators/impl/impl.hpp
      HINTS ${INTEGRATORXX_INCLUDE_DIR})
  

  # handle the QUIETLY and REQUIRED arguments and set INTEGRATORXX_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(IntegratorXX
      FOUND_VAR INTEGRATORXX_FOUND
      VERSION_VAR INTEGRATORXX_VERSION
      REQUIRED_VARS INTEGRATORXX_INCLUDE_DIRS)

  mark_as_advanced(INTEGRATORXX_INCLUDE_DIR INTEGRATORXX_INCLUDE_DIRS)

endif()
