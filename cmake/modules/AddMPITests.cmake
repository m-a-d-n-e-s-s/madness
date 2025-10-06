# AddMPITests.cmake
# 
# This module provides macros for adding MPI-based tests with different
# process counts and node configurations.
#
# add_mpi_tests(_component _test_name _nprocs _libs _labels)
#   - Adds multi-rank MPI tests for a single test executable
#   - Tests run on a single node with multiple processes (oversubscribed if needed)
#
# add_multinode_tests(_component _test_name _configs _libs _labels)
#   - Adds multi-node MPI tests for a single test executable
#   - Tests run across multiple nodes using a hostfile
#   - Only created if MADNESS_HOSTFILE or HOSTFILE environment variable is set
#
# Environment Variables (configure-time):
#   MADNESS_HOSTFILE or HOSTFILE - Path to MPI hostfile for multi-node tests
#
# Environment Variables (runtime - can be set when running ctest):
#   MADNESS_MPI_NODE_OPTIONS - Override default node mapping options for multi-node tests
#                              Example: MADNESS_MPI_NODE_OPTIONS="--bind-to none" ctest -R multinode
#                              Default: "--map-by node -N {procs_per_node}"

# Add multi-rank MPI tests for a single test
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
    
    # Add the MPI test
    add_test(NAME madness/test/${_component}/${_mpi_test_name}/run
             COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} 
                     ${MPIEXEC_PREFLAGS} --oversubscribe 
                     $<TARGET_FILE:${_test_name}> ${MPIEXEC_POSTFLAGS})
    
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

# Add multi-node MPI tests for a single test
# Usage: add_multinode_tests(component test_name "config1;config2" "libs" "labels")
# Each config should be in format "name:nprocs:procs_per_node"
# Example: "2nodes_4procs:4:2" means 4 total procs, 2 per node
macro(add_multinode_tests _component _test_name _configs _libs _labels)
  
  # Track that multi-node tests have been added
  set(MADNESS_HAS_MULTINODE_TESTS TRUE CACHE INTERNAL "Multi-node tests have been configured")
  
  if(NOT ENABLE_MPI OR NOT MPIEXEC_EXECUTABLE)
    message(STATUS "MPI not enabled or MPIEXEC not found, skipping multi-node tests for ${_test_name}")
    # Create placeholder tests that will skip at runtime with explanation
    foreach(CONFIG ${_configs})
      string(REPLACE ":" ";" CONFIG_LIST ${CONFIG})
      list(LENGTH CONFIG_LIST CONFIG_LENGTH)
      if(CONFIG_LENGTH EQUAL 3)
        list(GET CONFIG_LIST 0 CONFIG_NAME)
        set(_multinode_test_name "${_test_name}_${CONFIG_NAME}")
        add_test(NAME madness/test/${_component}/${_multinode_test_name}/run
                 COMMAND ${CMAKE_COMMAND} -E echo "SKIPPED: MPI not enabled or mpiexec not found for ${_test_name}_${CONFIG_NAME}")
        set_tests_properties(madness/test/${_component}/${_multinode_test_name}/run
                             PROPERTIES SKIP_RETURN_CODE 0
                             LABELS "${_labels};mpi;multinode")
      endif()
    endforeach()
    return()
  endif()
  
  # Ensure the test executable exists or will be created
  if(NOT TARGET ${_test_name})
    message(WARNING "Test target ${_test_name} does not exist. Make sure it is created before calling add_multinode_tests.")
  endif()
  
  # Check if a hostfile is specified via environment variable
  if(DEFINED ENV{MADNESS_HOSTFILE} AND EXISTS "$ENV{MADNESS_HOSTFILE}")
    set(HOSTFILE_OPTION "--hostfile" "$ENV{MADNESS_HOSTFILE}")
  elseif(DEFINED ENV{HOSTFILE} AND EXISTS "$ENV{HOSTFILE}")
    set(HOSTFILE_OPTION "--hostfile" "$ENV{HOSTFILE}")
  else()
    set(HOSTFILE_OPTION "")
  endif()
  
  # Only create tests if hostfile is available
  if(NOT HOSTFILE_OPTION)
    message(STATUS "No hostfile found (MADNESS_HOSTFILE or HOSTFILE), skipping multi-node tests for ${_test_name}")
    # Create placeholder tests that will skip at runtime with explanation
    foreach(CONFIG ${_configs})
      string(REPLACE ":" ";" CONFIG_LIST ${CONFIG})
      list(LENGTH CONFIG_LIST CONFIG_LENGTH)
      if(CONFIG_LENGTH EQUAL 3)
        list(GET CONFIG_LIST 0 CONFIG_NAME)
        set(_multinode_test_name "${_test_name}_${CONFIG_NAME}")
        add_test(NAME madness/test/${_component}/${_multinode_test_name}/run
                 COMMAND ${CMAKE_COMMAND} -E echo "SKIPPED: No hostfile specified (set MADNESS_HOSTFILE or HOSTFILE) for ${_test_name}_${CONFIG_NAME}")
        set_tests_properties(madness/test/${_component}/${_multinode_test_name}/run
                             PROPERTIES SKIP_RETURN_CODE 0
                             LABELS "${_labels};mpi;multinode")
      endif()
    endforeach()
    return()
  endif()
  
  foreach(CONFIG ${_configs})
    # Parse config string: "name:nprocs:procs_per_node"
    string(REPLACE ":" ";" CONFIG_LIST ${CONFIG})
    list(LENGTH CONFIG_LIST CONFIG_LENGTH)
    
    if(CONFIG_LENGTH EQUAL 3)
      list(GET CONFIG_LIST 0 CONFIG_NAME)
      list(GET CONFIG_LIST 1 NPROC)
      list(GET CONFIG_LIST 2 PROCS_PER_NODE)
      
      # Create test name
      set(_multinode_test_name "${_test_name}_${CONFIG_NAME}")
      
      # Create a CMake wrapper script that will check for MADNESS_MPI_NODE_OPTIONS at runtime
      set(_wrapper_script_template "${CMAKE_CURRENT_BINARY_DIR}/run_${_multinode_test_name}.cmake.in")
      set(_wrapper_script "${CMAKE_CURRENT_BINARY_DIR}/run_${_multinode_test_name}.cmake")
      
      # Convert HOSTFILE_OPTION list to space-separated strings for command
      string(REPLACE ";" " " HOSTFILE_OPTION_STR "${HOSTFILE_OPTION}")
      string(REPLACE ";" " " MPIEXEC_PREFLAGS_STR "${MPIEXEC_PREFLAGS}")
      string(REPLACE ";" " " MPIEXEC_POSTFLAGS_STR "${MPIEXEC_POSTFLAGS}")
      
      # Write CMake script template with placeholder for target file
      file(WRITE ${_wrapper_script_template} "# Auto-generated wrapper script for MPI test\n")
      file(APPEND ${_wrapper_script_template} "# Check if MADNESS_MPI_NODE_OPTIONS is set in environment\n")
      file(APPEND ${_wrapper_script_template} "if(DEFINED ENV{MADNESS_MPI_NODE_OPTIONS})\n")
      file(APPEND ${_wrapper_script_template} "  set(NODE_OPTIONS \"\$ENV{MADNESS_MPI_NODE_OPTIONS}\")\n")
      file(APPEND ${_wrapper_script_template} "else()\n")
      file(APPEND ${_wrapper_script_template} "  set(NODE_OPTIONS \"--map-by node -N ${PROCS_PER_NODE}\")\n")
      file(APPEND ${_wrapper_script_template} "endif()\n")
      file(APPEND ${_wrapper_script_template} "# Convert NODE_OPTIONS to list\n")
      file(APPEND ${_wrapper_script_template} "separate_arguments(NODE_OPTIONS_LIST UNIX_COMMAND \"\${NODE_OPTIONS}\")\n")
      file(APPEND ${_wrapper_script_template} "# Execute MPI command\n")
      file(APPEND ${_wrapper_script_template} "execute_process(\n")
      file(APPEND ${_wrapper_script_template} "  COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${HOSTFILE_OPTION_STR} \${NODE_OPTIONS_LIST} ${MPIEXEC_PREFLAGS_STR} \"$<TARGET_FILE:${_test_name}>\" ${MPIEXEC_POSTFLAGS_STR}\n")
      file(APPEND ${_wrapper_script_template} "  RESULT_VARIABLE _result\n")
      file(APPEND ${_wrapper_script_template} ")\n")
      file(APPEND ${_wrapper_script_template} "if(NOT _result EQUAL 0)\n")
      file(APPEND ${_wrapper_script_template} "  message(FATAL_ERROR \"Test failed with exit code \${_result}\")\n")
      file(APPEND ${_wrapper_script_template} "endif()\n")
      
      # Use file(GENERATE) to resolve generator expressions at build time
      file(GENERATE OUTPUT ${_wrapper_script} INPUT ${_wrapper_script_template})
      
      # Add the multi-node test using the wrapper script
      add_test(NAME madness/test/${_component}/${_multinode_test_name}/run
               COMMAND ${CMAKE_COMMAND} -P ${_wrapper_script})
      
      # Set test properties
      set_tests_properties(madness/test/${_component}/${_multinode_test_name}/run
                           PROPERTIES DEPENDS madness/test/${_component}/build 
                           LABELS "${_labels};mpi;multinode")
      
      # Add dependency to component unittests
      if(TARGET ${_component}_unittests-madness)
        add_dependencies(${_component}_unittests-madness ${_test_name})
      endif()
    else()
      message(WARNING "Invalid multi-node config format '${CONFIG}'. Expected 'name:nprocs:procs_per_node'")
    endif()
  endforeach()
  
endmacro()
