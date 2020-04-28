# src/apps/moldft

add_mad_executable(moldft moldft.cc MADchem)
<<<<<<< Updated upstream
=======
target_link_libraries(moldft MADchem ${DISABLEPIE_LINKER_FLAG})
>>>>>>> Stashed changes

add_mad_executable(mcpfit mcpfit.cc MADchem)

#add_mad_executable(preal preal.cc MADchem)

add_mad_executable(testpg testpg.cc MADchem)

add_mad_executable(testperiodic_moldft testperiodic.cc MADChem)

install(TARGETS moldft mcpfit DESTINATION "${MADNESS_INSTALL_BINDIR}")
