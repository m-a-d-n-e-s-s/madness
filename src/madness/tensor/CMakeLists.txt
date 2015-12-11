# src/madness/tensor

# Add targets for generated header files
if(PYTHONINTERP_FOUND)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/mtxm_gen.h
      COMMAND ${PYTHON_EXECUTABLE} genmtxm.py > mtxm_gen.h
      MAIN_DEPENDENCY genmtxm.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating mtxm_gen.h")
  add_custom_command(
      OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/tensoriter_spec.h
             ${CMAKE_CURRENT_SOURCE_DIR}/tensor_spec.h
      COMMAND ${PYTHON_EXECUTABLE} tempspec.py
      MAIN_DEPENDENCY tempspec.py
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
      COMMENT "Generating tensoriter_spec.h and tensor_spec.h")
elseif(NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/mtxm_gen.h
    OR NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tensoriter_spec.h.h
    OR NOT EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/tensor_spec.h)
  message(FATAL_ERROR "PYTHON executable is not available to generate required source files")
endif()

# Source lists for MADtensor
set(MADTENSOR_HEADERS 
    aligned.h mxm.h tensorexcept.h tensoriter_spec.h type_data.h basetensor.h
    tensor.h tensor_macros.h vector_factory.h mtxmq.h slice.h tensoriter.h
    tensor_spec.h vmath.h systolic.h gentensor.h srconf.h distributed_matrix.h
    tensortrain.h)
set(MADTENSOR_SOURCES tensor.cc tensoriter.cc basetensor.cc mtxmq.cc vmath.cc)
if(USE_X86_64_ASM OR USE_X86_32_ASM)
  list(APPEND MADTENSOR_SOURCES mtxmq_asm.S)
endif()
  
# Source lists for MADlinalg
set(MADLINALG_HEADERS cblas.h tensor_lapack.h clapack.h  lapack_functions.h solvers.h)
set(MADLINALG_SOURCES lapack.cc solvers.cc)
if(ELEMENTAL_FOUND)
  list(APPEND MADLINALG_HEADERS elem.h)
  list(APPEND MADLINALG_SOURCE  elem.cc)
endif()

# Create MADtensor and MADlinalg libraries
add_mad_library(tensor MADTENSOR_SOURCES MADTENSOR_HEADERS "misc" "madness/tensor")
add_mad_library(linalg MADLINALG_SOURCES MADLINALG_HEADERS "tensor" "madness/tensor")
target_link_libraries(MADlinalg PUBLIC ${LAPACK_LIBRARIES})
if(USE_X86_64_ASM)
  target_compile_definitions(MADtensor-obj PRIVATE -DX86_64)
elseif(USE_X86_32_ASM)
  target_compile_definitions(MADtensor-obj PRIVATE -DX86_32)
endif()
if(HAVE_IBMBGQ)
  add_custom_target(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/new_mtxmq/bests/libMADMTXM.a
      COMMAND TARGET=BGQ ${CMAKE_MAKE_PROGRAM} 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/new_mtxmq/bests
      COMMENT "Building libMADMTXM.a")
  target_link_libraries(MADtensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/new_mtxmq/bests/libMADMTXM.a)
elseif(HAVE_IBMBGP)
  add_custom_target(OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/new_mtxmq/bests/libMADMTXM.a
      COMMAND TARGET=BGP ${CMAKE_MAKE_PROGRAM} 
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/new_mtxmq/bests
      COMMENT "Building libMADMTXM.a")
  target_link_libraries(MADtensor PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/new_mtxmq/bests/libMADMTXM.a)
endif()

# Add unit tests
if(ENABLE_UNITTESTS)
  
  # The list of unit test source files
  set(TENSOR_TEST_SOURCES test_tensor.cc oldtest.cc test_scott.cc test_mtxmq.cc
      jimkernel.cc test_distributed_matrix.cc test_Zmtxmq.cc test_systolic.cc)
  if(ENABLE_GENTENSOR)
    list(APPEND TENSOR_TEST_SOURCES test_gentensor.cc)
  endif()
  set(LINALG_TEST_SOURCES test_linalg.cc test_solvers.cc testseprep.cc)

  add_unittests(tensor TENSOR_TEST_SOURCES "MADtensor;MADgtest")
  add_unittests(linalg LINALG_TEST_SOURCES "MADlinalg;MADgtest")
  
endif()