# src/apps/chem

add_definitions(-DMRA_CHEMDATA_DIR=\"${CMAKE_CURRENT_SOURCE_DIR}\")

# Set the CHEM sources and header files
set(MADCHEM_HEADERS
    AC.h
    atomutil.h
    basis.h
    BSHApply.h
    CalculationParameters.h
    CC2.h
    mp3.h
    ccpairfunction.h
    CCPotentials.h
    CCStructures.h
    corepotential.h
    correlationfactor.h
    diamagneticpotentialfactor.h
    distpm.h
    electronic_correlation_factor.h
    ESInterface.h
    exchangeoperator.h
    gaussian.h
    gth_pseudopotential.h
    GuessFactory.h
    IntegratorXX.h
    localizer.h
    lowrankfunction.h
    masks_and_boxes.h
    molecularbasis.h
    molecular_functors.h
    molecular_optimizer.h
    MolecularOrbitals.h
    molecule.h
    molopt.h
    mp2.h
    nemo.h
    NWChem.h
    oep.h
    pcm.h
    PNOF12Potentials.h
    PNOGuessFunctions.h
    PNO.h
    PNOParameters.h
    PNOStructures.h
    PNOTensors.h
    pointgroupoperator.h
    pointgroupsymmetry.h
    polynomial.h
    potentialmanager.h
    projector.h
    QCPropertyInterface.h
    SAP.h
    SCF.h
    SCFOperators.h
    SCFProtocol.h
    TDHF.h
    vibanal.h
    write_test_input.h
    xcfunctional.h
    zcis.h
    znemo.h
)
set(MADCHEM_SOURCES
    AC.cc
    atomutil.cc
    CC2.cc
    mp3.cc
    ccpairfunction.cc
    CCPotentials.cc
    CCStructures.cc
    corepotential.cc
    correlationfactor.cc
    diamagneticpotentialfactor.cc
    distpm.cc
    exchangeoperator.cc
    gaussian.cc
    gth_pseudopotential.cc
    GuessFactory.cc
    lda.cc
    localizer.cc
    molecularbasis.cc
    MolecularOrbitals.cc
    molecule.cc
    mp2.cc
    nemo.cc
    NWChem.cc
    oep.cc
    pcm.cc
    pointgroupsymmetry.cc
    potentialmanager.cc
    polynomial.cc
    SAP.cc
    SCF.cc
    SCFOperators.cc
    TDHF.cc
    vibanal.cc
    zcis.cc
    znemo.cc
    PNO.cpp  PNOF12Potentials.cpp  PNOGuessFunctions.cpp  PNOParameters.cpp  PNOStructures.cpp
)
if(TARGET Libxc::xc)
  list(APPEND MADCHEM_SOURCES xcfunctional_libxc.cc)
else()
  list(APPEND MADCHEM_SOURCES xcfunctional_ldaonly.cc)
endif()

# Create the MADchem library
add_mad_library(chem MADCHEM_SOURCES MADCHEM_HEADERS "mra" "madness/chem/")

  # N.B. MADchem too flaky for unity builds
  set_target_properties(MADchem-obj PROPERTIES UNITY_BUILD FALSE)

  set(targetname MADchem)

  # Add dependencies for MADchem
  target_include_directories(${targetname} PUBLIC
      $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/src/apps/>)
  if (TARGET Libxc::xc)
    target_link_libraries(${targetname} PUBLIC Libxc::xc)
  endif ()

  if (PCM_FOUND)
    target_include_directories(${targetname} PUBLIC ${PCM_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${PCM_LIBRARIES})
  endif ()

  if (INTEGRATORXX_FOUND)
    target_include_directories(${targetname} PUBLIC ${INTEGRATORXX_INCLUDE_DIRS})
    target_link_libraries(${targetname} PUBLIC ${INTEGRATORXX_LIBRARIES})
  endif ()

  if (BOOST_FOUND)
    target_link_libraries(${targetname} PUBLIC Boost::headers)
  endif ()

# Add executables
if (NOT MADNESS_BUILD_LIBRARIES_ONLY)
  add_mad_executable(plotxc "plotxc.cc;xcfunctional.h" "MADchem")
endif()

# Install files
install(FILES sto-3g sto-6g 6-31g coredata/mcp coredata/mcp2 coredata/mcp_guess structure_library
    DESTINATION "${MADNESS_INSTALL_DATADIR}"
    COMPONENT chem)

# Add unit tests    
if(BUILD_TESTING)
  
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  # The list of unit test source files
  set(CHEM_TEST_SOURCES_SHORT test_pointgroupsymmetry.cc test_masks_and_boxes.cc
  test_qc.cc test_MolecularOrbitals.cc test_BSHApply.cc test_projector.cc)
  set(CHEM_TEST_SOURCES_LONG test_localizer.cc test_ccpairfunction.cc test_low_rank_function.cc)
  if (TARGET Libxc::xc)
    list(APPEND CHEM_TEST_SOURCES_SHORT test_dft.cc )
    list(APPEND CHEM_TEST_SOURCES_LONG test_SCFOperators.cc)
  endif(TARGET Libxc::xc)
  if (INTEGRATORXX_FOUND)
    list(APPEND CHEM_TEST_SOURCES_SHORT test_IntegratorXX.cc )
  endif(INTEGRATORXX_FOUND)

  add_unittests(chem "${CHEM_TEST_SOURCES_SHORT}" "MADchem;MADgtest" "unittests;short")
  add_unittests(chem "${CHEM_TEST_SOURCES_LONG}" "MADchem;MADgtest" "unittests;long")

  # Create other test executables not included in the unit tests ... consider these executables (unlike unit tests)
  if (NOT MADNESS_BUILD_LIBRARIES_ONLY)
    set(CHEM_OTHER_TESTS testxc benchmark_exchange_operator)
    foreach(_test ${CHEM_OTHER_TESTS})
      add_mad_executable(${_test} "${_test}.cc" "MADchem")
    endforeach()
  endif()

endif()
