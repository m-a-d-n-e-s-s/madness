# src/madness/tensor

# Add targets for generated header files
if(PYTHONINTERP_FOUND)
  add_custom_command(
      OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/tensoriter_spec.h
             ${CMAKE_CURRENT_SOURCE_DIR}/tensor_spec.h
      COMMAND ${PYTHON_EXECUTABLE} tempspec.py
      MAIN_DEPENDENCY tempspec.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating tensoriter_spec.h and tensor_spec.h")
elseif(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tensoriter_spec.h.h
    OR NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tensor_spec.h)
  message(FATAL_ERROR "PYTHON executable is not available to generate required source files")
endif()

# Source lists for MADtensor
set(MADTENSOR_HEADERS 
    aligned.h mxm.h tensorexcept.h tensoriter_spec.h type_data.h basetensor.h
    tensor.h tensor_macros.h vector_factory.h slice.h tensoriter.h
    tensor_spec.h vmath.h systolic.h gentensor.h srconf.h distributed_matrix.h
    tensortrain.h)
set(MADTENSOR_SOURCES tensor.cc tensoriter.cc basetensor.cc vmath.cc)

# logically these headers should be part of their own library (MADclapack)
# however CMake right now does not support a mechanism to properly handle header-only libs.
# so will keep this a part of MADlinalg, add an install rule for these header only
# and will remember that these files depend on LAPACK
#
# Source lists for MADclapack
set(MADCLAPACK_HEADERS cblas.h clapack.h)  # this part of MADlinalg is purely independent of MADtensor

# Source lists for MADlinalg
set(MADLINALG_HEADERS ${MADCLAPACK_HEADERS} tensor_lapack.h solvers.h elem.h)
set(MADLINALG_SOURCES lapack.cc solvers.cc elem.cc)

# Create libraries MADtensor and MADlinalg
add_mad_library(tensor MADTENSOR_SOURCES MADTENSOR_HEADERS "misc" "madness/tensor")
#add_mad_library(clapack MADLINALG_HEADERS "" "madness/tensor")
add_mad_library(linalg MADLINALG_SOURCES MADLINALG_HEADERS "tensor" "madness/tensor")
target_link_libraries(MADlinalg PUBLIC ${LAPACK_LIBRARIES}) # logically it's MADclapack that depends on LAPACK

if(HAVE_IBMBGQ)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/mtxmq_bg/libMADMTXM.a
      COMMAND TARGET=BGQ ${CMAKE_MAKE_PROGRAM} 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mtxmq_bg
      COMMENT "Building libMADMTXM.a")
  add_custom_target (libmadmtxm DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mtxmq_bg/libMADMTXM.a)
  add_dependencies(MADtensor libmadmtxm)
  target_link_libraries(MADtensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mtxmq_bg/libMADMTXM.a)
elseif(HAVE_IBMBGP)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/mtxmq_bg/libMADMTXM.a
      COMMAND TARGET=BGP ${CMAKE_MAKE_PROGRAM} 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/mtxmq_bg
      COMMENT "Building libMADMTXM.a")
  add_custom_target (libmadmtxm DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/mtxmq_bg/libMADMTXM.a)
  add_dependencies(MADtensor libmadmtxm)
  target_link_libraries(MADtensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/mtxmq_bg/libMADMTXM.a)
endif()

if(ELEMENTAL_FOUND AND MADNESS_HAS_ELEMENTAL_EMBEDDED)
  target_link_libraries(MADlinalg PUBLIC ${ELEMENTAL_PACKAGE_NAME})
endif()

# add custom target for "installing" MADclapack headers
install(FILES ${MADCLAPACK_HEADERS}
        DESTINATION "${MADNESS_INSTALL_INCLUDEDIR}/madness/tensor"
        COMPONENT clapack)
add_custom_target(install-clapack
    COMMAND ${CMAKE_COMMAND} -DCOMPONENT=clapack -P ${CMAKE_BINARY_DIR}/cmake_install.cmake
    COMMENT "Installing clapack headers")

# Add unit tests
if(ENABLE_UNITTESTS)
  
  # The list of unit test source files
  set(TENSOR_TEST_SOURCES test_tensor.cc oldtest.cc test_mtxmq.cc
      jimkernel.cc test_distributed_matrix.cc test_Zmtxmq.cc test_systolic.cc)
  if(ENABLE_GENTENSOR)
    list(APPEND TENSOR_TEST_SOURCES test_gentensor.cc)
  endif()
  set(LINALG_TEST_SOURCES test_linalg.cc test_solvers.cc testseprep.cc)

  add_unittests(tensor TENSOR_TEST_SOURCES "MADtensor;MADgtest")
  add_unittests(linalg LINALG_TEST_SOURCES "MADlinalg;MADgtest")
  
endif()
