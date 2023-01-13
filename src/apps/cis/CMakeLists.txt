# src/apps/moldft

add_mad_executable(cis cis.cpp "MADchem")
add_dependencies(applications-madness cis)

add_scripted_tests(test_callable.py cis "short;application")
add_scripted_tests(test_energy_he.py cis "medium;application")
add_scripted_tests(test_symmetry_h2o.py cis "long;application")

install(TARGETS cis DESTINATION "${MADNESS_INSTALL_BINDIR}")
