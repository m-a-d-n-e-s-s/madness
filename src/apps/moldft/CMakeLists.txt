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

install(TARGETS moldft mcpfit DESTINATION "${MADNESS_INSTALL_BINDIR}")