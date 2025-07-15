# src/apps/moldft

add_mad_executable(madqc madqc.cpp "MADchem;MADresponse2")
add_dependencies(applications-madness madqc)


# integration tests
#add_scripted_tests(test_nemo_energy.py madqc "short;applications")
add_scripted_tests(test_nemo_localization.py madqc "short;applications")
add_scripted_tests(test_moldft_energy.py madqc "short;applications")
add_scripted_tests(test_moldft_energy_mpi_parallel.py madqc "short;applications")
add_scripted_tests(test_cc2_callable.py madqc "short;applications")
add_scripted_tests(test_oep_energy.py madqc "short;applications")
#add_scripted_tests(test_cis_energy_he.py madqc "medium;applications")
#add_scripted_tests(test_cis_symmetry_h2o.py madqc "long;applications")
#add_scripted_tests(test_cis_callable.py madqc "short;applications")
# add_scripted_tests(test_response.py madqc "short;applications")
# add_scripted_tests(test_znemo_energy.py madqc "short;applications")


install(TARGETS madqc DESTINATION "${MADNESS_INSTALL_BINDIR}")
