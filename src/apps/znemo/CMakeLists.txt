# src/apps/moldft

add_mad_executable(znemo znemo.cc "MADchem")
add_dependencies(applications-madness znemo)
add_scripted_tests(test_energy.py znemo "medium;application")


install(TARGETS znemo DESTINATION "${MADNESS_INSTALL_BINDIR}")
