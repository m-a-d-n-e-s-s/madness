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
# Environment Variables:
#   MADNESS_MPI_NODE_OPTIONS - Override default node mapping options for multi-node tests
#                              Example: "--map-by node -N 2" or "--bind-to none"
#   MADNESS_HOSTFILE or HOSTFILE - Path to MPI hostfile for multi-node tests

# Add multi-rank MPI tests for a single test
# Usage: add_mpi_tests(component test_name "2;4;8" "libs" "labels")
macro(add_mpi_tests _component _test_name _nprocs _libs _labels)
  
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
  
  # Check for custom node options from environment variable
  if(DEFINED ENV{MADNESS_MPI_NODE_OPTIONS})
    # Parse the environment variable into a list
    string(REPLACE " " ";" NODE_OPTIONS_LIST "$ENV{MADNESS_MPI_NODE_OPTIONS}")
    message(STATUS "Using custom MPI node options from MADNESS_MPI_NODE_OPTIONS: $ENV{MADNESS_MPI_NODE_OPTIONS}")
  else()
    # Use default node options (empty list - will be set per config)
    set(NODE_OPTIONS_LIST "")
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
      
      # Set node mapping options - use custom if provided, otherwise use defaults
      if(NODE_OPTIONS_LIST)
        set(NODE_OPTIONS ${NODE_OPTIONS_LIST})
      else()
        # Default OpenMPI node mapping options
        set(NODE_OPTIONS "--map-by" "node" "-N" "${PROCS_PER_NODE}")
      endif()
      
      # Add the multi-node test
      add_test(NAME madness/test/${_component}/${_multinode_test_name}/run
               COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} 
                       ${HOSTFILE_OPTION} ${NODE_OPTIONS} ${MPIEXEC_PREFLAGS}
                       $<TARGET_FILE:${_test_name}> ${MPIEXEC_POSTFLAGS})
      
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
