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
#   MADNESS_MPI_NODE_OPTIONS - Specify MPI execution options; replaces the default
#                              options entirely (see below).
#                              Example: MADNESS_MPI_NODE_OPTIONS="--bind-to none" ctest -R mpi
#                              Example: MADNESS_MPI_NODE_OPTIONS="--hostfile /path/to/hostfile --map-by node" ctest -R mpi
#   MAD_CHECK_BINDING        - Defaulted to OFF for these tests, see below.  Set it
#                              explicitly to keep MADNESS' start-up binding check enabled.
#   MAD_NUM_THREADS          - Defaulted to a per-rank budget for these tests, see below.
#                              Set it explicitly to pick the pool size yourself.
#
# Default launch options
#   mpiexec pins each rank to a single core when a machine has more cores than ranks,
#   so the 2-rank unit tests end up running all their MADNESS threads on one core.
#   That is not just slow, it is crippling: on a 96-core node test_vectormacrotask
#   takes 206 s pinned versus 14 s unpinned, and test_eval does not finish at all.
#   These tests are about correctness, not placement, so unless the user overrides
#   MADNESS_MPI_NODE_OPTIONS we launch them unbound.  The flag spelling is
#   vendor-specific, hence the version-string sniffing; unknown vendors get no flag.
#
#   Unbinding alone still leaves the ranks fighting each other, because each one
#   sizes its pool independently: MADNESS defaults MAD_NUM_THREADS to
#   (logical cores - 1), so an N-rank test asks for N*(cores-1) threads on one
#   host.  On a 12-core arm64 Mac the 2-rank tests ran ~2.5x oversubscribed, and
#   it cost more than an order of magnitude in wall time: test_eval_mpi2 301 s
#   and test_halo_mpi2 15.6 s at the old pool size, versus 11.7 s and 2.2 s once
#   the pool was budgeted.  So the wrapper also defaults MAD_NUM_THREADS to
#   (cores / nprocs - 2), the two spare cores per rank being the MADNESS
#   communication thread and the MPI implementation's own progress thread.
#
#   That budget is only applied to the default, single-host launch.  The wrapper
#   runs on the host ctest runs on, so the core count it probes is that host's --
#   fine when the ranks stay there, wrong the moment they do not.  Setting
#   MADNESS_MPI_NODE_OPTIONS is how this module spells "I am choosing the layout"
#   (a hostfile, --map-by node, a batch allocation), so when it is set the budget
#   is skipped and MADNESS' own per-rank default applies instead.  That default is
#   (cores - 1) measured on the host the rank actually lands on, which for one
#   rank per node is already right; a global number derived here would not be.
#   Set MAD_NUM_THREADS explicitly to budget such a run.
#
#   The matching CTest PROCESSORS property tells 'ctest -j' what one run of the
#   test actually costs, so it does not schedule several of them side by side and
#   put the oversubscription straight back.  CTest fixes that reservation at
#   configure time, so it describes the default budget and nothing else: override
#   MAD_NUM_THREADS at run time and the reservation no longer matches what the
#   test asks for -- larger overrides under-reserve, smaller ones over-reserve.
#   Run overridden tests without 'ctest -j', or reserve by hand.
#
#   MAD_CHECK_BINDING is defaulted to OFF on top of all that.  The budget above is
#   what its aggregate-demand test asks for, so it would now mostly pass -- but
#   cmake_host_system_information reports the host's cores, not the cpuset a
#   container or batch allocation has confined the job to.  Where those differ the
#   budget overshoots, and with the check on that turns a slow test into a hard
#   abort.  'MAD_CHECK_BINDING=ON ctest -L mpi' still exercises it.

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

  # Host size, probed once per configure.  Used only for the PROCESSORS property,
  # which CTest reads at configure time; the wrapper re-probes at run time for the
  # thread budget, where the machine running the test is the one that matters.
  if(NOT DEFINED MADNESS_TEST_HOST_CORES)
    cmake_host_system_information(RESULT MADNESS_TEST_HOST_CORES
                                  QUERY NUMBER_OF_LOGICAL_CORES)
  endif()

  # Default launch options: run unbound, so the ranks' threads get the whole node.
  # Only Open MPI and MPICH-family launchers are spelled out; anything else (Intel
  # MPI, Cray, vendor wrappers) keeps mpiexec's own defaults rather than risking an
  # unrecognized flag, and can be tuned through MADNESS_MPI_NODE_OPTIONS.
  # Probed once per configure, not once per test.
  if(NOT DEFINED MADNESS_MPI_DEFAULT_OPTIONS)
    # Ask the launcher itself rather than trusting MPI_<lang>_LIBRARY_VERSION_STRING,
    # which FindMPI leaves empty in some configurations.  Open MPI's mpiexec reports
    # "mpiexec (OpenRTE) 4.1.8" or "(Open MPI) 5.x"; MPICH/Hydra reports "HYDRA".
    execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version
                    OUTPUT_VARIABLE _mad_mpiexec_version
                    ERROR_VARIABLE _mad_mpiexec_version
                    OUTPUT_STRIP_TRAILING_WHITESPACE
                    ERROR_STRIP_TRAILING_WHITESPACE)
    set(_mad_mpi_version "${_mad_mpiexec_version} ${MPI_C_LIBRARY_VERSION_STRING} ${MPI_CXX_LIBRARY_VERSION_STRING}")
    string(REGEX REPLACE "[\r\n]+" " " _mad_mpi_version "${_mad_mpi_version}")
    if(_mad_mpi_version MATCHES "Intel")
      # Intel MPI is Hydra-based but does not take -bind-to; pinning there is driven
      # by I_MPI_PIN / I_MPI_PIN_DOMAIN, which we leave to the user.
      set(MADNESS_MPI_DEFAULT_OPTIONS "")
      message(STATUS "add_mpi_tests: Intel MPI detected; launching MPI tests without "
                     "binding flags -- export I_MPI_PIN=0 if they run pinned to one core")
    elseif(_mad_mpi_version MATCHES "Open MPI|OpenRTE|Open RTE")
      set(MADNESS_MPI_DEFAULT_OPTIONS "--bind-to none")
    elseif(_mad_mpi_version MATCHES "MPICH|HYDRA|Hydra")
      set(MADNESS_MPI_DEFAULT_OPTIONS "-bind-to none")
    else()
      set(MADNESS_MPI_DEFAULT_OPTIONS "")
      message(STATUS "add_mpi_tests: unrecognized MPI (\"${_mad_mpi_version}\"), "
                     "launching MPI tests without binding flags; "
                     "set MADNESS_MPI_NODE_OPTIONS if the tests run pinned to one core")
    endif()
    message(STATUS "add_mpi_tests: default mpiexec options: \"${MADNESS_MPI_DEFAULT_OPTIONS}\"")
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
    file(APPEND ${_wrapper_script_template} "  # Default: launch unbound, otherwise every rank's threads share a single core\n")
    file(APPEND ${_wrapper_script_template} "  separate_arguments(MPI_OPTIONS_LIST UNIX_COMMAND \"${MADNESS_MPI_DEFAULT_OPTIONS}\")\n")
    file(APPEND ${_wrapper_script_template} "endif()\n")
    file(APPEND ${_wrapper_script_template} "# Disable the start-up CPU-binding check (see thread_binding.h).  On a machine\n")
    file(APPEND ${_wrapper_script_template} "# with more cores than ranks mpiexec pins each rank to a single core by default,\n")
    file(APPEND ${_wrapper_script_template} "# which the check rejects -- correctly for production, but these unit tests are\n")
    file(APPEND ${_wrapper_script_template} "# about correctness, not throughput.  An explicit MAD_CHECK_BINDING in the\n")
    file(APPEND ${_wrapper_script_template} "# environment wins, so 'MAD_CHECK_BINDING=ON ctest -L mpi' still exercises it.\n")
    file(APPEND ${_wrapper_script_template} "if(NOT DEFINED ENV{MAD_CHECK_BINDING})\n")
    file(APPEND ${_wrapper_script_template} "  set(ENV{MAD_CHECK_BINDING} \"OFF\")\n")
    file(APPEND ${_wrapper_script_template} "endif()\n")
    file(APPEND ${_wrapper_script_template} "# Budget the per-rank thread pool.  Left alone every rank sizes its own pool at\n")
    file(APPEND ${_wrapper_script_template} "# (cores - 1), and the ranks then time-share the host several times over.  Leave\n")
    file(APPEND ${_wrapper_script_template} "# two cores per rank for the MADNESS communication thread and the MPI progress\n")
    file(APPEND ${_wrapper_script_template} "# thread.  An explicit MAD_NUM_THREADS in the environment wins.\n")
    file(APPEND ${_wrapper_script_template} "# Only for the default layout: this script runs on the host ctest runs on, so\n")
    file(APPEND ${_wrapper_script_template} "# the core count below is that host's.  Once MADNESS_MPI_NODE_OPTIONS puts the\n")
    file(APPEND ${_wrapper_script_template} "# ranks somewhere else -- a hostfile, --map-by node, a batch allocation -- a\n")
    file(APPEND ${_wrapper_script_template} "# number derived here describes the wrong machine.  MADNESS' own default,\n")
    file(APPEND ${_wrapper_script_template} "# (cores - 1) measured wherever the rank lands, is the better answer there.\n")
    file(APPEND ${_wrapper_script_template} "if(NOT DEFINED ENV{MAD_NUM_THREADS} AND NOT DEFINED ENV{MADNESS_MPI_NODE_OPTIONS})\n")
    file(APPEND ${_wrapper_script_template} "  cmake_host_system_information(RESULT _ncores QUERY NUMBER_OF_LOGICAL_CORES)\n")
    file(APPEND ${_wrapper_script_template} "  math(EXPR _nthreads \"\${_ncores} / ${NPROC} - 2\")\n")
    file(APPEND ${_wrapper_script_template} "  if(_nthreads LESS 1)\n")
    file(APPEND ${_wrapper_script_template} "    set(_nthreads 1)\n")
    file(APPEND ${_wrapper_script_template} "  endif()\n")
    file(APPEND ${_wrapper_script_template} "  set(ENV{MAD_NUM_THREADS} \"\${_nthreads}\")\n")
    file(APPEND ${_wrapper_script_template} "  message(STATUS \"MAD_NUM_THREADS=\${_nthreads} (${NPROC} rank(s) on \${_ncores} logical cores)\")\n")
    file(APPEND ${_wrapper_script_template} "elseif(NOT DEFINED ENV{MAD_NUM_THREADS})\n")
    file(APPEND ${_wrapper_script_template} "  message(STATUS \"MADNESS_MPI_NODE_OPTIONS is set: leaving the thread count to MADNESS' \"\n")
    file(APPEND ${_wrapper_script_template} "                 \"per-host default.  Set MAD_NUM_THREADS to budget this run.\")\n")
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
    
    # What one run costs the host, so 'ctest -j' does not start several of these
    # side by side and undo the budget.  Mirrors the wrapper's arithmetic, with the
    # two spare cores per rank added back in.  CTest fixes this at configure time,
    # so it describes the default budget only: a run-time MAD_NUM_THREADS override
    # or a MADNESS_MPI_NODE_OPTIONS layout will not match it (see the header).
    math(EXPR _mad_rank_threads "${MADNESS_TEST_HOST_CORES} / ${NPROC} - 2")
    if(_mad_rank_threads LESS 1)
      set(_mad_rank_threads 1)
    endif()
    math(EXPR _mad_test_procs "${NPROC} * (${_mad_rank_threads} + 2)")

    # Set test properties
    set_tests_properties(madness/test/${_component}/${_mpi_test_name}/run
                         PROPERTIES DEPENDS madness/test/${_component}/build 
                         LABELS "${_labels};mpi"
                         PROCESSORS ${_mad_test_procs})
    
    # Add dependency to component unittests
    if(TARGET ${_component}_unittests-madness)
      add_dependencies(${_component}_unittests-madness ${_test_name})
    endif()
  endforeach()
  
endmacro()
