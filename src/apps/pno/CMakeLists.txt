# src/apps/pno

add_mad_executable(pno pno.cpp MADchem)
add_dependencies(applications-madness pno)

install(TARGETS pno DESTINATION "${MADNESS_INSTALL_BINDIR}")
