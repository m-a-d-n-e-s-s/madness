# src/madness/tensor

# Add targets for generated header files
find_package(Python COMPONENTS Interpreter)
if(TARGET Python::Interpreter)
  add_custom_command(
      OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/tensoriter_spec.h
             ${CMAKE_CURRENT_SOURCE_DIR}/tensor_spec.h
      COMMAND ${Python_EXECUTABLE} tempspec.py
      MAIN_DEPENDENCY tempspec.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating tensoriter_spec.h and tensor_spec.h")
elseif(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tensoriter_spec.h.h
    OR NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tensor_spec.h)
  message(FATAL_ERROR "Python executable is not available to generate required source files")
endif()

# Source lists for MADtensor
set(MADTENSOR_HEADERS 
    aligned.h mxm.h tensorexcept.h tensoriter_spec.h type_data.h basetensor.h
    tensor.h tensor_macros.h vector_factory.h slice.h tensoriter.h
    tensor_spec.h vmath.h systolic.h gentensor.h srconf.h distributed_matrix.h
    tensortrain.h SVDTensor.h tensor_json.hpp)
set(MADTENSOR_SOURCES tensor.cc tensoriter.cc basetensor.cc vmath.cc)

# logically these headers should be part of their own library (MADclapack)
# however CMake right now does not support a mechanism to properly handle header-only libs.
# so will keep this a part of MADlinalg, add an install rule for these header only
# and will remember that these files depend on LAPACK
#
# Source lists for MADclapack
set(MADCLAPACK_HEADERS cblas.h cblas_types.h
        clapack.h clapack_fortran.h
        lapacke_types.h linalg_wrappers.h ) # this part of MADlinalg is purely independent of MADtensor

# Source lists for MADlinalg
set(MADLINALG_HEADERS ${MADCLAPACK_HEADERS} tensor_lapack.h solvers.h elem.h jacobi.h)
set(MADLINALG_SOURCES lapack.cc solvers.cc elem.cc SVDTensor.cc RandomizedMatrixDecomposition.cc linalg_wrappers.cc )

set(MADGENTENSOR_HEADERS gentensor.h srconf.h distributed_matrix.h tensortrain.h SVDTensor.h
        lowranktensor.h RandomizedMatrixDecomposition.h)
set(MADGENTENSOR_SOURCES SVDTensor.cc RandomizedMatrixDecomposition.cc)

# Create libraries MADtensor and MADlinalg
add_mad_library(tensor MADTENSOR_SOURCES MADTENSOR_HEADERS "misc" "madness/tensor")
#add_mad_library(clapack MADLINALG_HEADERS "" "madness/tensor")
add_mad_library(linalg MADLINALG_SOURCES MADLINALG_HEADERS "tensor" "madness/tensor")

if(ENABLE_GENTENSOR)
  add_mad_library(gentensor MADGENTENSOR_SOURCES MADGENTENSOR_HEADERS "linalg" "madness/tensor")
endif()

  set(l_target MADlinalg)
  set(t_target MADtensor)

  target_link_libraries(${l_target} PUBLIC ${LAPACK_LIBRARIES}) # logically it's MADclapack that depends on LAPACK
  target_link_libraries(${t_target} PUBLIC ${LAPACK_LIBRARIES}) # some of the tests depend on MKL/Blas
  if (LAPACK_INCLUDE_DIRS)
    target_include_directories(${l_target} PUBLIC ${LAPACK_INCLUDE_DIRS})
    target_include_directories(${t_target} PUBLIC ${LAPACK_INCLUDE_DIRS})
  endif (LAPACK_INCLUDE_DIRS)
  if (LAPACK_COMPILE_OPTIONS)
    target_compile_options(${l_target} PUBLIC ${LAPACK_COMPILE_OPTIONS})
    target_compile_options(${t_target} PUBLIC ${LAPACK_COMPILE_OPTIONS})
  endif (LAPACK_COMPILE_OPTIONS)
  if (LAPACK_COMPILE_DEFINITIONS)
    target_compile_definitions(${l_target} PUBLIC ${LAPACK_COMPILE_DEFINITIONS})
    target_compile_definitions(${t_target} PUBLIC ${LAPACK_COMPILE_DEFINITIONS})
  endif (LAPACK_COMPILE_DEFINITIONS)

  if (HAVE_IBMBGQ)
    message(FATAL_ERROR "IBM BG/Q support ended")
  elseif (HAVE_IBMBGP)
    message(FATAL_ERROR "IBM BG/P support ended")
  endif ()

  if (ELEMENTAL_FOUND AND MADNESS_HAS_ELEMENTAL_EMBEDDED)
    target_link_libraries(${l_target} PUBLIC ${ELEMENTAL_PACKAGE_NAME})
  endif ()

# Add unit tests
if(BUILD_TESTING)
  
  # The list of unit test source files
  set(TENSOR_TEST_SOURCES test_tensor.cc oldtest.cc test_mtxmq.cc
      jimkernel.cc test_distributed_matrix.cc test_Zmtxmq.cc test_systolic.cc)
  set(LINALG_TEST_SOURCES test_linalg.cc test_solvers.cc testseprep.cc test_jacobi.cc)

  if(ENABLE_GENTENSOR)
    list(APPEND LINALG_TEST_SOURCES test_gentensor.cc test_lowranktensor.cc)
  endif()

  add_unittests(tensor "${TENSOR_TEST_SOURCES}" "MADtensor;MADgtest" "unittests;short")
  add_unittests(linalg "${LINALG_TEST_SOURCES}" "MADlinalg;MADgtest" "unittests;short")
  
endif()
