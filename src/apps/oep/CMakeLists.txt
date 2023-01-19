# src/apps/moldft

add_mad_executable(oep oep.cc "MADchem")
add_dependencies(applications-madness oep)
add_scripted_tests(test_energy.py oep "medium;application")

install(TARGETS oep DESTINATION "${MADNESS_INSTALL_BINDIR}")
