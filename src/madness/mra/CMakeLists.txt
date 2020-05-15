# src/madness/mra

add_definitions(-DMRA_DATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Set the MRA sources and header files
set(MADMRA_HEADERS
    adquad.h  funcimpl.h  indexit.h  legendre.h  operator.h  vmra.h
    funcdefaults.h  key.h  mra.h  power.h  qmprop.h  twoscale.h lbdeux.h
    mraimpl.h  funcplot.h  function_common_data.h function_factory.h
    function_interface.h gfit.h convolution1d.h simplecache.h derivative.h
    displacements.h functypedefs.h sdf_shape_3D.h sdf_domainmask.h vmra1.h
    leafop.h nonlinsol.h)
set(MADMRA_SOURCES
    mra1.cc mra2.cc mra3.cc mra4.cc mra5.cc mra6.cc startup.cc legendre.cc 
    twoscale.cc qmprop.cc)

# Create the MADmra library
add_mad_library(mra MADMRA_SOURCES MADMRA_HEADERS "linalg;tinyxml;muparser" "madness/mra")

# Create executables
add_mad_executable(mraplot "mraplot.cc" "MADmra")

# Install the MADmra library
install(TARGETS mraplot DESTINATION "${MADNESS_INSTALL_BINDIR}")
install(FILES autocorr coeffs gaussleg
    DESTINATION "${MADNESS_INSTALL_DATADIR}"
    COMPONENT mra)

# Add unit tests    
if(ENABLE_UNITTESTS)

  # The list of unit test source files
  add_library(libtest_sepop EXCLUDE_FROM_ALL test_sepop.cc)
  target_link_libraries(libtest_sepop MADmra)
  
  set(MRA_TEST_SOURCES testbsh.cc testproj.cc 
      testpdiff.cc testdiff1Db.cc testgconv.cc testopdir.cc testinnerext.cc 
      testgaxpyext.cc testvmra.cc)
  add_unittests(mra "${MRA_TEST_SOURCES}" "MADmra;MADgtest")
  set(MRA_SEPOP_TEST_SOURCES testsuite.cc
      testper.cc)
  add_unittests(mra_sepop "${MRA_SEPOP_TEST_SOURCES}" "libtest_sepop;MADgtest")
  
  # Test executables that are not run with unit tests
  set(MRA_OTHER_TESTS testperiodic testbc testqm test6
      testdiff1D testdiff2D testdiff3D)
  
  foreach(_test ${MRA_OTHER_TESTS})  
    add_mad_executable(${_test} "${_test}.cc" "MADmra")
  endforeach()
  
endif()