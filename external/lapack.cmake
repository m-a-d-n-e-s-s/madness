# Find BLAS and LAPACK.
include(CheckCFortranFunctionExists)
include(CMakePushCheckState)

if(NOT LAPACK_LIBRARIES)

  if(ENABLE_MKL)
    find_package(MKL)
    
    if(MKL_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${MKL_LIBRARIES})
      set(HAVE_INTEL_MKL 1)
    endif()
  endif()
  
  if(ENABLE_ACML AND NOT LAPACK_FOUND)
    find_package(ACML)
    
    if(ACML_FOUND)
      set(LAPACK_FOUND TRUE)
      set(LAPACK_LIBRARIES ${ACML_LIBRARIES})
      set(HAVE_ACML 1)
    endif()
  endif()
  
  # Search for system specific BLAS/LAPACK checks
  if(NOT LAPACK_FOUND AND CMAKE_SYSTEM_NAME MATCHES "Darwin")
    # Accelerate is always present, so no need to search
    set(LAPACK_LIBRARIES "-framework Accelerate")
    set(LAPACK_FOUND TRUE)
  endif()
  
  # Search for netlib lapack and blas libraries
  if(NOT LAPACK_FOUND)
    find_library(LAPACK_lapack_LIBRARY lapack)
    find_library(LAPACK_blas_LIBRARY blas)
    
    if(LAPACK_lapack_LIBRARY AND LAPACK_blas_LIBRARY)
      set(LAPACK_LIBRARIES ${LAPACK_lapack_LIBRARY} ${LAPACK_blas_LIBRARY})
      set(LAPACK_FOUND TRUE)
    endif()
  endif()  

endif()

cmake_push_check_state()

set(CMAKE_REQUIRED_LIBRARIES ${LAPACK_LIBRARIES} ${CMAKE_REQUIRED_LIBRARIES} 
    ${CMAKE_THREAD_LIBS_INIT})

# Verify that we can link against BLAS
check_c_fortran_function_exists(sgemm BLAS_WORKS)

if(BLAS_WORKS)
  message(STATUS "A library with BLAS API found.")
else()
  message(FATAL_ERROR "Uable to link against BLAS function. Specify the BLAS library in LAPACK_LIBRARIES.")
endif()

# Verify that we can link against LAPACK
check_c_fortran_function_exists(cheev LAPACK_WORKS)

if(LAPACK_WORKS)
  message(STATUS "A library with LAPACK API found.")
else()
  message(FATAL_ERROR "Uable to link against LAPACK function. Specify the LAPACK library in LAPACK_LIBRARIES.")
endif()

set(LAPACK_FOUND TRUE)
message(STATUS "Found LAPACK: ${LAPACK_LIBRARIES}")

cmake_pop_check_state()

# Set the fortran mangling scheme.
if(LAPACK_WORKS STREQUAL "cheev_")
  set(FORTRAN_LINKAGE_LCU 1)
elseif(LAPACK_WORKS STREQUAL "cheev")
  set(FORTRAN_LINKAGE_LC 1)
elseif(LAPACK_WORKS STREQUAL "cheev__")
  set(FORTRAN_LINKAGE_LCUU 1)
elseif(LAPACK_WORKS STREQUAL "CHEEV")
  set(FORTRAN_LINKAGE_UC 1)
elseif(LAPACK_WORKS STREQUAL "CHEEV_")
  set(FORTRAN_LINKAGE_UCU 1)
endif()
