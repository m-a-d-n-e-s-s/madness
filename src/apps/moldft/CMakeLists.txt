# src/apps/moldft

add_mad_executable(moldft moldft.cc MADchem)

add_executable(mcpfit mcpfit.cc)
target_link_libraries(mcpfit MADchem ${DISABLEPIE_LINKER_FLAG})

#add_executable(preal preal.cc)
#target_link_libraries(preal MADchem ${DISABLEPIE_LINKER_FLAG})

add_executable(testpg testpg.cc)
target_link_libraries(testpg MADchem ${DISABLEPIE_LINKER_FLAG})

add_executable(testperiodic_moldft testperiodic.cc)
target_link_libraries(testperiodic_moldft MADchem ${DISABLEPIE_LINKER_FLAG})

install(TARGETS moldft mcpfit DESTINATION "${MADNESS_INSTALL_BINDIR}")