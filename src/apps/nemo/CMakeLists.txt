# src/apps/moldft

add_mad_executable(nemo nemo.cc "MADchem")
add_dependencies(applications-madness nemo)

add_scripted_tests(test_energy.py nemo "short;applications")
add_scripted_tests(test_localization.py nemo "long;applications")

install(TARGETS nemo DESTINATION "${MADNESS_INSTALL_BINDIR}")
