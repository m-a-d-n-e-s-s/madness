# src/apps/moldft

add_mad_executable(znemo znemo.cc "MADchem;MADgit")
add_dependencies(applications-madness znemo)


install(TARGETS znemo DESTINATION "${MADNESS_INSTALL_BINDIR}")
