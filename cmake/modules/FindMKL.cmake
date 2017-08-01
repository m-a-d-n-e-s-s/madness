# - Try to find MKL
# Input variables:
#  MKL_ROOT_DIR         - The MKL install directory
#  MKL_INCLUDE_DIR      - The MKL include directory
#  MKL_LIBRARY          - The MKL library directory
#  FORTRAN_INTEGER_SIZE - Set the integer size (4 or 8 bytes)
#  BLA_STATIC           - Prefer to link against static lbraries
# Output variables:
#  MKL_FOUND            - System has MKL
#  MKL_INCLUDE_DIRS     - The tbb include directories
#  MKL_LIBRARIES        - The libraries needed to use MKL
#  MKL_VERSION          - The version string for MKL

if(NOT MKL_FOUND)

  # If the user did not specify an MKL root directory, check for the MKLROOT
  # environment variable or the presence of /opt/intel/mkl.
  if(NOT MKL_ROOT_DIR OR NOT DEFINED MKL_ROOT_DIR)
    if(EXISTS $ENV{MKLROOT})
      set(MKL_ROOT_DIR "$ENV{MKLROOT}")
    elseif(EXISTS /opt/intel/mkl)
      set(MKL_ROOT_DIR /opt/intel/mkl)
    endif()
  endif()

  if(MKL_ROOT_DIR)
    set(MKL_INCLUDE_DIR ${MKL_ROOT_DIR}/include 
        CACHE PATH "The include directory for MKL")
    # Set the MKL library directory if not specified by the user.
    if(CMAKE_SYSTEM_NAME MATCHES "Linux")
      set(MKL_LIBRARY ${MKL_ROOT_DIR}/lib/intel64 
          CACHE PATH "The library directory for MKL")
    else()
      set(MKL_LIBRARY ${MKL_ROOT_DIR}/lib
          CACHE PATH "The library directory for MKL")
    endif()
  endif()
  
  if(FORTRAN_INTEGER_SIZE EQUAL 4)
    set(MKL_INT_TYPE "lp64")
  elseif(FORTRAN_INTEGER_SIZE EQUAL 8)
    set(MKL_INT_TYPE "ilp64")
  else()
    set(MKL_INT_TYPE "lp64")
  endif()
  
  # There are no user specified components, but we use the component handling
  # mechanism to make sure we find all the required libraries.
  set(MKL_FIND_COMPONENTS mkl_intel_${MKL_INT_TYPE} mkl_core mkl_sequential)
  set(MKL_FIND_REQUIRED_mkl_intel_${MKL_INT_TYPE} TRUE)
  set(MKL_FIND_REQUIRED_mkl_core TRUE)
  set(MKL_FIND_REQUIRED_mkl_sequential TRUE)
  
  # Search for MKL header files
  find_path(MKL_INCLUDE_DIRS mkl.h
      HINTS ${MKL_INCLUDE_DIR})
      
  # Get MKL version
  if(MKL_INCLUDE_DIRS)
    file(READ "${MKL_INCLUDE_DIRS}/mkl_version.h" _mkl_version_file)
    string(REGEX REPLACE ".*#define __INTEL_MKL__ ([0-9]+).*" "\\1"
            MKL_VERSION_MAJOR "${_mkl_version_file}")
    string(REGEX REPLACE ".*#define __INTEL_MKL_MINOR__ ([0-9]+).*" "\\1"
            MKL_VERSION_MINOR "${_mkl_version_file}")
    string(REGEX REPLACE ".*#define __INTEL_MKL_UPDATE__ ([0-9]+).*" "\\1"
            MKL_VERSION_UPDATE "${_mkl_version_file}")
    set(MKL_VERSION "${MKL_VERSION_MAJOR}.${MKL_VERSION_MINOR}.${MKL_VERSION_UPDATE}")
    unset(_mkl_version_header)
  endif()
  
  # Search for MKL libraries
  foreach(_lib ${MKL_FIND_COMPONENTS})
    if(BLA_STATIC)
      find_library(MKL_${_lib}_LIBRARY ${CMAKE_STATIC_LIBRARY_PREFIX}${_lib}${CMAKE_STATIC_LIBRARY_SUFFIX} ${_lib}
          HINTS ${MKL_LIBRARY})
    else()
      find_library(MKL_${_lib}_LIBRARY ${_lib}
          HINTS ${MKL_LIBRARY})
    endif()
    if(MKL_${_lib}_LIBRARY)
      set(MKL_${_lib}_FOUND TRUE)
      list(APPEND MKL_LIBRARIES ${MKL_${_lib}_LIBRARY})
    else()
      set(MKL_${_lib}_FOUND FALSE)
    endif()
  endforeach()
  
  # Set LAPACK_LIBRARIES variable if MKL was found
  if(MKL_mkl_core_FOUND)
    set(MKL_FOUND TRUE)
    if(UNIX AND BLA_STATIC)
      set(MKL_LIBRARIES -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group -lm -ldl
          CACHE STRING "The Intel MKL libraries")
    else()
      set(MKL_LIBRARIES ${MKL_LIBRARIES} -lm -ldl
          CACHE STRING "The Intel MKL libraries")
    endif()
  endif()
  
  # handle the QUIETLY and REQUIRED arguments and set MKL_FOUND to TRUE
  # if all listed variables are TRUE
  find_package_handle_standard_args(MKL
      FOUND_VAR MKL_FOUND
      VERSION_VAR MKL_VERSION 
      REQUIRED_VARS MKL_LIBRARIES MKL_INCLUDE_DIRS
      HANDLE_COMPONENTS)

endif()