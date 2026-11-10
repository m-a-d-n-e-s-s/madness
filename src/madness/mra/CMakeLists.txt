# src/madness/mra

add_definitions(-DMRA_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Set the MRA sources and header files
set(MADMRA_HEADERS
    adquad.h  funcimpl.h  indexit.h  legendre.h  operator.h  vmra.h
    funcdefaults.h  key.h  mra.h  power.h  qmprop.h  twoscale.h lbdeux.h
    mraimpl.h  funcplot.h  function_common_data.h function_factory.h
    function_interface.h gfit.h convolution1d.h simplecache.h derivative.h
    displacements.h functypedefs.h sdf_shape_3D.h sdf_domainmask.h
    leafop.h nonlinsol.h macrotaskq.h macrotaskpartitioner.h QCCalculationParametersBase.h
    commandlineparser.h operatorinfo.h bc.h kernelrange.h mw.h memory_measurement.h composite_product_op.h
    interpolation_1d.h phandler.h)
# mra1.cc also defines dimension-independent statics, so it is always compiled;
# mra<N>.cc (N>1) holds nothing but the explicit instantiations for dimension N.
set(MADMRA_SOURCES
    mra1.cc startup.cc legendre.cc
    twoscale.cc qmprop.cc QCCalculationParametersBase.cc)
foreach(_dim IN LISTS MADNESS_MRA_DIMENSIONS)
  if (NOT _dim EQUAL 1)
    list(APPEND MADMRA_SOURCES mra${_dim}.cc)
  endif()
endforeach()

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

  # MRA dimensions each test needs; see madness_mra_filter_targets() (unlisted tests need none)
  set(mra_dims_testbsh 3)
  set(mra_dims_testproj 1 2 3)
  set(mra_dims_testpdiff 3)
  set(mra_dims_testdiff1Db 1)
  set(mra_dims_testgconv 1 2)
  set(mra_dims_testopdir 3)
  set(mra_dims_testinnerext 1 2 3 4 5 6)
  set(mra_dims_testgaxpyext 3)
  set(mra_dims_testvmra 1 2 3)
  set(mra_dims_test_vectormacrotask 2 3)
  set(mra_dims_test_cloud 3)
  set(mra_dims_test_tree_state 2)
  set(mra_dims_testsolver 1 2)
  set(mra_dims_test_macrotaskpartitioner 1)
  set(mra_dims_test_memory_measurement 1 2 3 4)
  set(mra_dims_test_eval 1 2 3 6)
  set(mra_dims_test_mul_sparse 3)
  set(mra_dims_test_halo 3)
  set(mra_dims_test_displacements 3)
  set(mra_dims_test_truncate_tol 3)
  set(mra_dims_test_periodic_truncation 3)
  set(mra_dims_test_images_operator 3)
  set(mra_dims_test_anisotropic_reach 3)
  set(mra_dims_testsuite 1 2 3 4)
  set(mra_dims_testper 1)
  set(mra_dims_testperiodic 1 3)
  set(mra_dims_testbc 3)
  set(mra_dims_testqm 1)
  set(mra_dims_test6 1 2 3 4 5 6)
  set(mra_dims_testdiff1D 1)
  set(mra_dims_testdiff2D 2)
  set(mra_dims_testdiff3D 3)
  set(mra_dims_test_composite_product 1 2 3 4 5 6)

  # The list of unit test source files
  add_library(libtest_sepop EXCLUDE_FROM_ALL test_sepop.cc)
  target_link_libraries(libtest_sepop MADmra)

  set(MRA_TEST_SOURCES testbsh.cc testproj.cc
      testpdiff.cc testdiff1Db.cc testgconv.cc testopdir.cc testinnerext.cc
      testgaxpyext.cc testvmra.cc test_vectormacrotask.cc test_cloud.cc test_tree_state.cc testsolver.cc
      test_macrotaskpartitioner.cc test_QCCalculationParametersBase.cc test_memory_measurement.cc
      test_keybox.cc test_eval.cc test_mul_sparse.cc test_halo.cc interp3.cc test_displacements.cc
      test_truncate_tol.cc test_periodic_truncation.cc test_images_operator.cc
      test_anisotropic_reach.cc)
  madness_mra_filter_targets(MRA_TEST_SOURCES)
  add_unittests(mra "${MRA_TEST_SOURCES}" "MADmra;MADgtest" "unittests;short")
  set(MRA_SEPOP_TEST_SOURCES testsuite.cc
      testper.cc)
  madness_mra_filter_targets(MRA_SEPOP_TEST_SOURCES)
  if (MRA_SEPOP_TEST_SOURCES)
    add_unittests(mra_sepop "${MRA_SEPOP_TEST_SOURCES}" "libtest_sepop;MADgtest" "unittests;short")
  endif()

  # Add MPI tests for test_memory_measurement
  # MPI options can be set at runtime via MADNESS_MPI_NODE_OPTIONS environment variable
  # Example: MADNESS_MPI_NODE_OPTIONS="--hostfile /path/to/hostfile --map-by node" ctest -R mpi
  # use with more processes as:
  # add_mpi_tests(mra test_memory_measurement "2;4;8" "MADmra;MADgtest" "unittests;short")
  # (test_mul_sparse: the tree-state and norm checkpoints use gop.max/gop.sum and branch on a
  #  replicated state, so only a multi-rank run exercises them;
  #  test_halo: the halo only carries data between ranks, so the single-rank run cannot exercise it)
  set(MRA_MPI_TESTS test_memory_measurement test_cloud test_vectormacrotask test_eval
      test_mul_sparse test_halo)
  madness_mra_filter_targets(MRA_MPI_TESTS)
  foreach(_test IN LISTS MRA_MPI_TESTS)
    add_mpi_tests(mra ${_test} "2" "MADmra;MADgtest" "unittests;short")
  endforeach()

  # Test executables that are not run with unit tests ... consider these executables (unlike unit tests)
  if (NOT MADNESS_BUILD_LIBRARIES_ONLY)
    set(MRA_OTHER_TESTS testperiodic testbc testqm test6
        testdiff1D testdiff2D testdiff3D test_composite_product)
    madness_mra_filter_targets(MRA_OTHER_TESTS)

    foreach(_test ${MRA_OTHER_TESTS})
      add_mad_executable(${_test} "${_test}.cc" "MADmra")
    endforeach()
  endif()

endif()
