# src/apps/moldft

add_mad_executable(madqc macqc.cc "MADchem")
add_dependencies(applications-madness madcq)

install(TARGETS madqc DESTINATION "${MADNESS_INSTALL_BINDIR}")
