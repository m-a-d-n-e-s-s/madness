# - CMAKE Project file for the MADNESS package
# The following variables are defined:
#  MADNESS_INCLUDE_DIRS         - The MADNESS include directory
#  MADNESS_LIBRARIES            - The MADNESS libraries and their dependencies
#  MADNESS_<COMPONENT>_LIBRARY  - The MADNESS COMPONENT library
#  MADNESS_COMPILE_FLAGS        - Compile flags required to build with MADNESS
#  MADNESS_LINKER_FLAGS         - Linker flags required to link with MADNESS
#  MADNESS_VERSION              - MADNESS version number
#  MADNESS_F77_INTEGER_SIZE     - The default F77 integer size used for BLAS calls 

# Compute paths
set(MADNESS_SOURCE_DIR /home/kimjins/madness/madness/build/..)
set(MADNESS_BUILD_DIR /home/kimjins/madness/madness/build)
set(MADNESS_INCLUDE_DIRS ${MADNESS_SOURCE_DIR}/src ${MADNESS_BUILD_DIR}/src)

# Set package version
set(MADNESS_VERSION "0.10.1")

# Set compile and link flags, and remove optimization and debug flags
string(REGEX REPLACE "-(O[0-9s]|g[0-9]?)([ ]+|$)" "" MADNESS_COMPILE_FLAGS " -DMPICH_SKIP_MPICXX=1 -DOMPI_SKIP_MPICXX=1  -std=c++11 -O3 -Wall -Wno-strict-aliasing -Wno-deprecated -Wno-unused-local-typedefs  -ffast-math -march=native")
string(STRIP "${MADNESS_COMPILE_FLAGS}" MADNESS_COMPILE_FLAGS)
string(REGEX REPLACE "-(O[0-9s]|g[0-9]?)([ ]+|$)" "" MADNESS_LINKER_FLAGS " -std=c++11 -O3 -Wall -Wno-strict-aliasing -Wno-deprecated -Wno-unused-local-typedefs  -ffast-math -march=native ")
string(STRIP "${MADNESS_LINKER_FLAGS}" MADNESS_LINKER_FLAGS)

set(MADNESS_BUILD_SHARED yes)
set(MADNESS_BUILD_STATIC yes)

# Set variables for shared library paths
if(MADNESS_BUILD_SHARED)
  set(MADNESS_MADworld_SHARED_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/world/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}MADworld${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(MADNESS_MADmisc_SHARED_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/misc/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}MADmisc${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(MADNESS_MADtensor_SHARED_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/tensor/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}MADtensor${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(MADNESS_MADlinalg_SHARED_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/tensor/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}MADlinalg${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(MADNESS_MADmra_SHARED_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/mra/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}MADmra${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(MADNESS_MADmuparser_SHARED_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/external/muparser/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}MADmuparser${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(MADNESS_MADtinyxml_SHARED_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/external/tinyxml/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}MADtinyxml${CMAKE_SHARED_LIBRARY_SUFFIX})
  set(MADNESS_MADchem_SHARED_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/apps/chem/.libs/${CMAKE_SHARED_LIBRARY_PREFIX}MADchem${CMAKE_SHARED_LIBRARY_SUFFIX})
endif()

# Set variables for static library paths
if(MADNESS_BUILD_STATIC)
  set(MADNESS_MADworld_STATIC_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/world/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}MADworld${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(MADNESS_MADmisc_STATIC_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/misc/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}MADmisc${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(MADNESS_MADtensor_STATIC_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/tensor/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}MADtensor${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(MADNESS_MADlinalg_STATIC_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/tensor/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}MADlinalg${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(MADNESS_MADmra_STATIC_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/mra/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}MADmra${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(MADNESS_MADmuparser_STATIC_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/external/muparser/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}MADmuparser${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(MADNESS_MADtinyxml_STATIC_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/madness/external/tinyxml/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}MADtinyxml${CMAKE_STATIC_LIBRARY_SUFFIX})
  set(MADNESS_MADchem_STATIC_LIBRARY 
      ${MADNESS_BUILD_DIR}/src/apps/chem/.libs/${CMAKE_STATIC_LIBRARY_PREFIX}MADchem${CMAKE_STATIC_LIBRARY_SUFFIX})
endif()

# Set default libraries
foreach(_lib MADchem MADmra MADtinyxml MADmuparser MADlinalg MADtensor MADmisc MADworld)
  if(MADNESS_BUILD_SHARED)
    set(MADNESS_${_lib}_LIBRARY ${MADNESS_${_lib}_SHARED_LIBRARY})
  else()
    set(MADNESS_${_lib}_LIBRARY ${MADNESS_${_lib}_STATIC_LIBRARY})
  endif()
  
  list(APPEND MADNESS_LIBRARIES ${MADNESS_${_lib}_LIBRARY})
endforeach()

list(APPEND MADNESS_LIBRARIES " -L/opt/intel/mkl/lib/intel64 -Wl,--start-group -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -Wl,--end-group -lpthread -lm -ldl  ")

# Set Fortran 77 integer size used by MADNESS
set(MADNESS_F77_INTEGER_SIZE 4)

unset(MAD_LIBRARY_PREFIX)
unset(MAD_LIBRARY_SUFFIX)
