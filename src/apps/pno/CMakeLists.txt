# src/apps/moldft

add_executable(pno pno.cpp)
target_link_libraries(pno MADchem ${DISABLEPIE_LINKER_FLAG})
add_executable(pno_integrals pno_integrals.cpp)
target_link_libraries(pno_integrals MADchem ${DISABLEPIE_LINKER_FLAG})

install(TARGETS pno DESTINATION "${MADNESS_INSTALL_BINDIR}")
install(TARGETS pno_integrals DESTINATION "${MADNESS_INSTALL_BINDIR}")
