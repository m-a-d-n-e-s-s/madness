# src/apps/moldft

add_mad_executable(nemo nemo.cc "MADchem")
add_dependencies(applications-madness nemo)

install(TARGETS nemo DESTINATION "${MADNESS_INSTALL_BINDIR}")
