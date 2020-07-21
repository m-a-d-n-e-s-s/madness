# src/apps/moldft

add_mad_executable(pno pno.cpp MADchem)

install(TARGETS pno DESTINATION "${MADNESS_INSTALL_BINDIR}")
