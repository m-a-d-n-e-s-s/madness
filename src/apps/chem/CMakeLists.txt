# src/apps/chem

add_definitions(-DMRA_CHEMDATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Set the CHEM sources and header files
set(MADCHEM_HEADERS
    correlationfactor.h molecule.h molecularbasis.h
    corepotential.h atomutil.h SCF.h xcfunctional.h
    mp2.h nemo.h potentialmanager.h gth_pseudopotential.h
    molecular_optimizer.h projector.h
    SCFOperators.h CCStructures.h CalculationParameters.h
    electronic_correlation_factor.h vibanal.h molopt.h TDHF.h CC2.h CCPotentials.h
    pcm.h SCFProtocol.h AC.h pointgroupoperator.h pointgroupsymmetry.h MolecularOrbitals.h 
	znemo.h zcis.h  diamagneticpotentialfactor.h oep.h
    PNO.h PNOF12Potentials.h PNOGuessFunctions.h PNOParameters.h PNOStructures.h PNOTensors.h)
set(MADCHEM_SOURCES
    correlationfactor.cc molecule.cc molecularbasis.cc vibanal.cc
    corepotential.cc atomutil.cc lda.cc
    distpm.cc SCF.cc gth_pseudopotential.cc nemo.cc mp2.cc pcm.cc
    SCFOperators.cc TDHF.cc GuessFactory.cc CCStructures.cc CC2.cc CCPotentials.cc 
    AC.cc pointgroupsymmetry.cc 
    znemo.cc zcis.cc diamagneticpotentialfactor.cc oep.cc
    PNO.cpp PNOF12Potentials.cpp PNOGuessFunctions.cpp PNOParameters.cpp PNOStructures.cpp
    QCCalculationParametersBase.cc MolecularOrbitals.cc)
if(LIBXC_FOUND)
  list(APPEND MADCHEM_SOURCES xcfunctional_libxc.cc)
else()
  list(APPEND MADCHEM_SOURCES xcfunctional_ldaonly.cc)
endif()

# Create the MADchem library
add_mad_library(chem MADCHEM_SOURCES MADCHEM_HEADERS "mra" "chem/")

  set(targetname MADchem)

  # Add dependencies for MADchem
  target_include_directories(${targetname} PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/apps/>)
  if (LIBXC_FOUND)
    target_include_directories(${targetname} PUBLIC ${LIBXC_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${LIBXC_LIBRARIES})
  endif ()

  if (PCM_FOUND)
    target_include_directories(${targetname} PUBLIC ${PCM_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${PCM_LIBRARIES})
  endif ()

  if (BOOST_FOUND)
    target_include_directories(${targetname} PUBLIC ${BOOST_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${BOOST_LIBRARIES})
  endif ()

# Add executables
add_mad_executable(plotxc "plotxc.cc;xcfunctional.h" "MADchem")

# Install files
install(FILES sto-3g sto-6g 6-31g coredata/mcp coredata/mcp2 coredata/mcp_guess
    DESTINATION "${MADNESS_INSTALL_DATADIR}"
    COMPONENT chem)

# Add unit tests    
if(ENABLE_UNITTESTS)
  
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  # The list of unit test source files
  set(CHEM_TEST_SOURCES test_SCFOperators.cc test_dft.cc test_pointgroupsymmetry.cc test_masks_and_boxes.cc
  test_QCCalculationParametersBase.cc test_qc.cc test_MolecularOrbitals.cc)
  
  add_unittests(chem "${CHEM_TEST_SOURCES}" "MADchem;MADgtest;MADgit")
  
  # Create other test executables not included in the unit tests
  set(CHEM_OTHER_TESTS testxc)
  foreach(_test ${CHEM_OTHER_TESTS})  
    add_mad_executable(${_test} "${_test}.cc" "MADchem;MADgit")
  endforeach()
  
endif()
