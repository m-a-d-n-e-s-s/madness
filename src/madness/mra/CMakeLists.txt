# src/madness/mra

add_definitions(-DMRA_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Set the MRA sources and header files
set(MADMRA_HEADERS
    adquad.h  funcimpl.h  indexit.h  legendre.h  operator.h  vmra.h
    funcdefaults.h  key.h  mra.h  power.h  qmprop.h  twoscale.h lbdeux.h
    mraimpl.h  funcplot.h  function_common_data.h function_factory.h
    function_interface.h gfit.h convolution1d.h simplecache.h derivative.h
    displacements.h functypedefs.h sdf_shape_3D.h sdf_domainmask.h vmra1.h
    leafop.h nonlinsol.h macrotaskq.h macrotaskpartitioner.h QCCalculationParametersBase.h
    commandlineparser.h operatorinfo.h bc.h kernelrange.h mw.h memory_measurement.h)
set(MADMRA_SOURCES
    mra1.cc mra2.cc mra3.cc mra4.cc mra5.cc mra6.cc startup.cc legendre.cc 
    twoscale.cc qmprop.cc QCCalculationParametersBase.cc)

# Create the MADmra library
add_mad_library(mra MADMRA_SOURCES MADMRA_HEADERS "linalg;tinyxml;muparser" "madness/mra")
install(FILES autocorr coeffs gaussleg ble-first.txt ble-second.txt b-spline-deriv1.txt b-spline-deriv2.txt b-spline-deriv3.txt
    DESTINATION "${MADNESS_INSTALL_DATADIR}"
    COMPONENT mra)

# Create executables
if (NOT MADNESS_BUILD_LIBRARIES_ONLY)
  add_mad_executable(mraplot "mraplot.cc" "MADmra") # installation fails with gnu-8 and gnu-9
  install(TARGETS mraplot DESTINATION "${MADNESS_INSTALL_BINDIR}")
endif()

# Add unit tests    
if(BUILD_TESTING)

  # The list of unit test source files
  add_library(libtest_sepop EXCLUDE_FROM_ALL test_sepop.cc)
  target_link_libraries(libtest_sepop MADmra)
  
  set(MRA_TEST_SOURCES testbsh.cc testproj.cc 
      testpdiff.cc testdiff1Db.cc testgconv.cc testopdir.cc testinnerext.cc 
      testgaxpyext.cc testvmra.cc, test_vectormacrotask.cc test_cloud.cc test_tree_state.cc testsolver.cc
      test_macrotaskpartitioner.cc test_QCCalculationParametersBase.cc test_memory_measurement.cc)
  add_unittests(mra "${MRA_TEST_SOURCES}" "MADmra;MADgtest" "unittests;short")
  set(MRA_SEPOP_TEST_SOURCES testsuite.cc
      testper.cc)
  add_unittests(mra_sepop "${MRA_SEPOP_TEST_SOURCES}" "libtest_sepop;MADgtest" "unittests;short")
  
  # Add MPI-specific tests for test_memory_measurement
  if(ENABLE_MPI AND MPIEXEC_EXECUTABLE)
    # Add MPI tests for different process counts and configurations
    set(MPI_PROCESS_COUNTS 2 4 8)
    foreach(NPROC ${MPI_PROCESS_COUNTS})
      add_test(NAME madness/test/mra/test_memory_measurement_mpi${NPROC}/run
               COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${NPROC} ${MPIEXEC_PREFLAGS} --oversubscribe 
                       $<TARGET_FILE:test_memory_measurement> ${MPIEXEC_POSTFLAGS})
      set_tests_properties(madness/test/mra/test_memory_measurement_mpi${NPROC}/run
                           PROPERTIES DEPENDS madness/test/mra/build LABELS "unittests;mpi;short")
      
      # Add dependency to build step
      add_dependencies(mra_unittests-madness test_memory_measurement)
    endforeach()
  endif()

  # Test executables that are not run with unit tests ... consider these executables (unlike unit tests)
  if (NOT MADNESS_BUILD_LIBRARIES_ONLY)
    set(MRA_OTHER_TESTS testperiodic testbc testqm test6
        testdiff1D testdiff2D testdiff3D)
  
    foreach(_test ${MRA_OTHER_TESTS})
      add_mad_executable(${_test} "${_test}.cc" "MADmra")
    endforeach()
  endif()
  
endif()
