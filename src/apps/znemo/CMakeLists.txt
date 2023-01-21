# src/apps/moldft

add_mad_executable(znemo znemo.cc "MADchem")
add_dependencies(applications-madness znemo)

if (!ENABLE_GENTENSOR)
    add_scripted_tests(test_energy.py znemo "medium;application")
endif()


install(TARGETS znemo DESTINATION "${MADNESS_INSTALL_BINDIR}")
