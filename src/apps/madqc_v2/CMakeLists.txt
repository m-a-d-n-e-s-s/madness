# src/apps/moldft

add_mad_executable(madqc2 madqc2.cpp "MADchem;MADresponse2")
add_dependencies(applications-madness madqc2)

# if(BUILD_TESTING)
# 	add_scripted_tests(test_callable.py cis "short;applications")
# 	add_scripted_tests(test_energy_he.py cis "medium;applications")
# 	add_scripted_tests(test_symmetry_h2o.py cis "long;applications")
# endif()

install(TARGETS madqc2 DESTINATION "${MADNESS_INSTALL_BINDIR}")
