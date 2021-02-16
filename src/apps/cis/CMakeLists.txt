# src/apps/moldft

add_mad_executable(cis cis.cpp "MADchem;MADgit")

install(TARGETS cis DESTINATION "${MADNESS_INSTALL_BINDIR}")
