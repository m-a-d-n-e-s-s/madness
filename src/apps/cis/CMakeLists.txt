# src/apps/moldft

add_executable(cis cis.cpp)
target_link_libraries(cis MADchem ${DISABLEPIE_LINKER_FLAG})

install(TARGETS cis DESTINATION "${MADNESS_INSTALL_BINDIR}")