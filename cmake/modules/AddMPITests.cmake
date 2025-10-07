# AddMPITests.cmake
# 
# This module provides a macro for adding MPI-based tests with different
# process counts. All MPI options are controlled via runtime environment variable.
#
# add_mpi_tests(_component _test_name _nprocs _libs _labels)
#   - Adds MPI tests for a single test executable
#   - Tests can run on single or multiple nodes depending on MPI configuration
#
# Environment Variables (runtime - can be set when running ctest):
#   MADNESS_MPI_NODE_OPTIONS - Specify MPI execution options
#                              Example: MADNESS_MPI_NODE_OPTIONS="--bind-to none" ctest -R mpi
#                              Example: MADNESS_MPI_NODE_OPTIONS="--hostfile /path/to/hostfile --map-by node" ctest -R mpi
#                              If not set, mpiexec will use its default behavior

# Add MPI tests for a single test
# Usage: add_mpi_tests(component test_name "2;4;8" "libs" "labels")
macro(add_mpi_tests _component _test_name _nprocs _libs _labels)
  
  # Track that MPI tests have been added
  set(MADNESS_HAS_MPI_TESTS TRUE CACHE INTERNAL "MPI tests have been configured")
  
  if(NOT ENABLE_MPI OR NOT MPIEXEC_EXECUTABLE)
    message(STATUS "MPI not enabled or MPIEXEC not found, skipping MPI tests for ${_test_name}")
    # Create placeholder tests that will skip at runtime with explanation
    foreach(NPROC ${_nprocs})
      set(_mpi_test_name "${_test_name}_mpi${NPROC}")
      add_test(NAME madness/test/${_component}/${_mpi_test_name}/run
               COMMAND ${CMAKE_COMMAND} -E echo "SKIPPED: MPI not enabled or mpiexec not found for ${_test_name}_mpi${NPROC}")
      set_tests_properties(madness/test/${_component}/${_mpi_test_name}/run
                           PROPERTIES SKIP_RETURN_CODE 0
                           LABELS "${_labels};mpi")
    endforeach()
    return()
  endif()
  
  # Ensure the test executable exists or will be created
  if(NOT TARGET ${_test_name})
    message(WARNING "Test target ${_test_name} does not exist. Make sure it is created before calling add_mpi_tests.")
  endif()
  
  foreach(NPROC ${_nprocs})
    # Create test name
    set(_mpi_test_name "${_test_name}_mpi${NPROC}")
    
    # Create a CMake wrapper script that will check for MADNESS_MPI_NODE_OPTIONS at runtime
    set(_wrapper_script_template "${CMAKE_CURRENT_BINARY_DIR}/run_${_mpi_test_name}.cmake.in")
    set(_wrapper_script "${CMAKE_CURRENT_BINARY_DIR}/run_${_mpi_test_name}.cmake")
    
    # Convert flags to space-separated strings
    string(REPLACE ";" " " MPIEXEC_PREFLAGS_STR "${MPIEXEC_PREFLAGS}")
    string(REPLACE ";" " " MPIEXEC_POSTFLAGS_STR "${MPIEXEC_POSTFLAGS}")
    
    # Write CMake script template with placeholder for target file
    file(WRITE ${_wrapper_script_template} "# Auto-generated wrapper script for MPI test\n")
    file(APPEND ${_wrapper_script_template} "# Check if MADNESS_MPI_NODE_OPTIONS is set in environment\n")
    file(APPEND ${_wrapper_script_template} "if(DEFINED ENV{MADNESS_MPI_NODE_OPTIONS})\n")
    file(APPEND ${_wrapper_script_template} "  set(MPI_OPTIONS \"\$ENV{MADNESS_MPI_NODE_OPTIONS}\")\n")
    file(APPEND ${_wrapper_script_template} "  # Convert MPI_OPTIONS to list\n")
    file(APPEND ${_wrapper_script_template} "  separate_arguments(MPI_OPTIONS_LIST UNIX_COMMAND \"\${MPI_OPTIONS}\")\n")
    file(APPEND ${_wrapper_script_template} "else()\n")
    file(APPEND ${_wrapper_script_template} "  set(MPI_OPTIONS_LIST \"\")\n")
    file(APPEND ${_wrapper_script_template} "endif()\n")
    file(APPEND ${_wrapper_script_template} "# Execute MPI command\n")
    file(APPEND ${_wrapper_script_template} "message(STATUS \"Running: ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} \${MPI_OPTIONS_LIST} ${MPIEXEC_PREFLAGS_STR} \\\"\$<TARGET_FILE:${_test_name}>\\\" ${MPIEXEC_POSTFLAGS_STR}\")\n")
    file(APPEND ${_wrapper_script_template} "execute_process(\n")
    file(APPEND ${_wrapper_script_template} "  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} \${MPI_OPTIONS_LIST} ${MPIEXEC_PREFLAGS_STR} \"$<TARGET_FILE:${_test_name}>\" ${MPIEXEC_POSTFLAGS_STR}\n")
    file(APPEND ${_wrapper_script_template} "  RESULT_VARIABLE _result\n")
    file(APPEND ${_wrapper_script_template} ")\n")
    file(APPEND ${_wrapper_script_template} "if(NOT _result EQUAL 0)\n")
    file(APPEND ${_wrapper_script_template} "  message(FATAL_ERROR \"Test failed with exit code \${_result}\")\n")
    file(APPEND ${_wrapper_script_template} "endif()\n")

    # Use file(GENERATE) to resolve generator expressions at build time
    file(GENERATE OUTPUT ${_wrapper_script} INPUT ${_wrapper_script_template})
    
    # Add the MPI test using the wrapper script
    add_test(NAME madness/test/${_component}/${_mpi_test_name}/run
             COMMAND ${CMAKE_COMMAND} -P ${_wrapper_script})
    
    # Set test properties
    set_tests_properties(madness/test/${_component}/${_mpi_test_name}/run
                         PROPERTIES DEPENDS madness/test/${_component}/build 
                         LABELS "${_labels};mpi")
    
    # Add dependency to component unittests
    if(TARGET ${_component}_unittests-madness)
      add_dependencies(${_component}_unittests-madness ${_test_name})
    endif()
  endforeach()
  
endmacro()
