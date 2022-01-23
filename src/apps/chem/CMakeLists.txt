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
    CCPotentials.h
    CCStructures.h
    commandlineparser.h
    corepotential.h
    correlationfactor.h
    diamagneticpotentialfactor.h
    electronic_correlation_factor.h
    ESInterface.h
    exchangeoperator.h
    gaussian.h
    gth_pseudopotential.h
    GuessFactory.h
    localizer.h
    masks_and_boxes.h
    molecularbasis.h
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
    QCCalculationParametersBase.h
    SCF.h
    SCFOperators.h
    SCFProtocol.h
    TDHF.h
    test_correlation_factor.h
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
    polynomial.cc
    QCCalculationParametersBase.cc
    SCF.cc
    SCFOperators.cc
    TDHF.cc
    vibanal.cc
    zcis.cc
    znemo.cc
    PNO.cpp  PNOF12Potentials.cpp  PNOGuessFunctions.cpp  PNOParameters.cpp  PNOStructures.cpp
)
# these source files require specific header order inclusion ... just skip them from unity build for now
set_source_files_properties(PNO.cpp  PNOF12Potentials.cpp  PNOGuessFunctions.cpp  PNOParameters.cpp  PNOStructures.cpp
        PROPERTIES SKIP_UNITY_BUILD_INCLUSION TRUE)
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
if(BUILD_TESTING)
  
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  # The list of unit test source files
  set(CHEM_TEST_SOURCES test_pointgroupsymmetry.cc test_masks_and_boxes.cc test_localizer.cc 
  test_QCCalculationParametersBase.cc test_qc.cc test_MolecularOrbitals.cc test_BSHApply.cc)
  if (LIBXC_FOUND)
    list(APPEND CHEM_TEST_SOURCES test_dft.cc test_SCFOperators.cc)
  endif(LIBXC_FOUND)
  
  add_unittests(chem "${CHEM_TEST_SOURCES}" "MADchem;MADgtest;MADgit")
  
  # Create other test executables not included in the unit tests
  set(CHEM_OTHER_TESTS testxc benchmark_exchange_operator)
  foreach(_test ${CHEM_OTHER_TESTS})  
    add_mad_executable(${_test} "${_test}.cc" "MADchem;MADgit")
  endforeach()
  
endif()
