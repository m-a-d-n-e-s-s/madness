# src/apps/moldft

add_mad_executable(oep oep.cc "MADchem;MADgit")
add_dependencies(applications-madness oep)

install(TARGETS oep DESTINATION "${MADNESS_INSTALL_BINDIR}")
