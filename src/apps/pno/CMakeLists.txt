# src/apps/moldft

add_executable(pno pno.cpp)
target_link_libraries(pno MADchem ${DISABLEPIE_LINKER_FLAG})
add_executable(make_canonical_integrals make_canonical_integrals.cpp)
target_link_libraries(make_canonical_integrals MADchem ${DISABLEPIE_LINKER_FLAG})

install(TARGETS pno DESTINATION "${MADNESS_INSTALL_BINDIR}")
