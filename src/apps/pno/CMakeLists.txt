# src/apps/moldft

add_executable(pno pno.cpp)
target_link_libraries(pno-mp2 MADchem ${DISABLEPIE_LINKER_FLAG})

install(TARGETS pno DESTINATION "${MADNESS_INSTALL_BINDIR}")