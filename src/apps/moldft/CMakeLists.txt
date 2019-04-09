# src/apps/moldft

add_executable(moldft moldft.cc)
target_link_libraries(moldft MADchem)

add_executable(mcpfit mcpfit.cc)
target_link_libraries(mcpfit MADchem)

#add_executable(preal preal.cc)
#target_link_libraries(preal MADchem)

add_executable(testpg testpg.cc)
target_link_libraries(testpg MADchem)

add_executable(testperiodic_moldft testperiodic.cc)
target_link_libraries(testperiodic_moldft MADchem)

add_executable(DFdriver DFdriver.cc DF.cc fcwf.cc gaussian.cc NWChem.cc polynomial.cc)
target_link_libraries(DFdriver MADchem MADmra)

add_executable(DFdriver2 DFdriver.cc DF.cc fcwf.cc gaussian.cc NWChem.cc polynomial.cc)
target_link_libraries(DFdriver2 MADchem MADmra)

add_executable(rk rk-3.cc)
target_link_libraries(rk MADchem MADmra -lgmp -lmpfr)

add_executable(dk DKhydrogen3.cc)
target_link_libraries(dk MADchem MADmra -lgmp -lmpfr)

install(TARGETS moldft mcpfit DESTINATION "${MADNESS_INSTALL_BINDIR}")
