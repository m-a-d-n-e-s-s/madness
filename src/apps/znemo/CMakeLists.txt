# src/apps/moldft

add_mad_executable(znemo znemo.cc "MADchem")
add_dependencies(applications-madness znemo)

if(BUILD_TESTING)
  if (NOT ENABLE_GENTENSOR)
    add_scripted_tests(test_energy.py znemo "medium;applications")
  endif()
endif()


install(TARGETS znemo DESTINATION "${MADNESS_INSTALL_BINDIR}")
