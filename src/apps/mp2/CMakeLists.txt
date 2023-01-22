# src/apps/moldft

add_mad_executable(mp2 mp2.cc "MADchem")
add_dependencies(applications-madness mp2)

if (ENABLE_GENTENSOR)
   add_scripted_tests(test_callable.py mp2 "short;applications")
endif()


install(TARGETS mp2 DESTINATION "${MADNESS_INSTALL_BINDIR}")
