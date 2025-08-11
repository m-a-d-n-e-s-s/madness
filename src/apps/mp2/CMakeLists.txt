# src/apps/mp2

add_mad_executable(mp2 mp2.cc "MADchem")
add_dependencies(applications-madness mp2)

install(TARGETS mp2 DESTINATION "${MADNESS_INSTALL_BINDIR}")
