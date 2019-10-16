# src/apps/moldft

add_executable(pno pno.cpp)
target_link_libraries(pno MADchem ${DISABLEPIE_LINKER_FLAG})

install(TARGETS pno DESTINATION "${MADNESS_INSTALL_BINDIR}")