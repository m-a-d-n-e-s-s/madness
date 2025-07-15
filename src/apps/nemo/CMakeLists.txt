# src/apps/moldft

add_mad_executable(nemo nemo.cc "MADchem")
add_dependencies(applications-madness nemo)

if(BUILD_TESTING)
  add_scripted_tests(test_localization.py nemo "long;applications")
endif()

install(TARGETS nemo DESTINATION "${MADNESS_INSTALL_BINDIR}")
