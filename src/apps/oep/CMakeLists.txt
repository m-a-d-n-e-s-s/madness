# src/apps/moldft

add_mad_executable(oep oep.cc "MADchem")
add_dependencies(applications-madness oep)

#if(BUILD_TESTING)
#  add_scripted_tests(../madqc_v2/test_energy.py oep "medium;applications")
#endif()

install(TARGETS oep DESTINATION "${MADNESS_INSTALL_BINDIR}")
