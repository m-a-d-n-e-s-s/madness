# src/apps/cis

add_mad_executable(cis cis.cpp "MADchem")
add_dependencies(applications-madness cis)

if(BUILD_TESTING)
  add_scripted_tests(test_callable.py cis "short;applications")
  add_scripted_tests(test_energy_he.py cis "medium;applications")
  add_scripted_tests(test_symmetry_h2o.py cis "long;applications")
endif()

install(TARGETS cis DESTINATION "${MADNESS_INSTALL_BINDIR}")
