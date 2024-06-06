# src/apps/moldft

add_mad_executable(madqc madqc.cc "MADchem")
add_dependencies(applications-madness madqc)

install(TARGETS madqc DESTINATION "${MADNESS_INSTALL_BINDIR}")
