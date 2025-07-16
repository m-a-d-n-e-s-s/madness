# src/apps/cc2

add_mad_executable(cc2 cc2.cc "MADchem")
add_dependencies(applications-madness cc2)

#if (ENABLE_GENTENSOR AND BUILD_TESTING)
#  add_scripted_tests(../madqc_v2/test_callable.py cc2 "short;applications") # runs he with one frozen orbital..
#endif()

install(TARGETS cc2 DESTINATION "${MADNESS_INSTALL_BINDIR}")
