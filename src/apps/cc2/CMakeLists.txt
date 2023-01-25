# src/apps/moldft

add_mad_executable(cc2 cc2.cc "MADchem")
add_dependencies(applications-madness cc2)


install(TARGETS cc2 DESTINATION "${MADNESS_INSTALL_BINDIR}")
