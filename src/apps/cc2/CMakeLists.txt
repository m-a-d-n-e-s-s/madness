# src/apps/moldft

add_mad_executable(cc2 cc2.cc "MADchem")
add_dependencies(applications-madness cc2)

if (ENABLE_GENTENSOR)
  add_scripted_tests(test_callable.py cc2 "short;application") # runs he with one frozen orbital..
endif()

install(TARGETS cc2 DESTINATION "${MADNESS_INSTALL_BINDIR}")
