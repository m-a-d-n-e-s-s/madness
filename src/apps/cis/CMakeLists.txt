# src/apps/moldft

add_mad_executable(cis cis.cpp "MADchem;MADgit")
add_dependencies(applications-madness cis)


install(TARGETS cis DESTINATION "${MADNESS_INSTALL_BINDIR}")
