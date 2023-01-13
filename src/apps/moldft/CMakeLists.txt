# src/apps/moldft

add_mad_executable(moldft moldft.cc MADchem)
add_dependencies(applications-madness moldft)

add_mad_executable(mcpfit mcpfit.cc MADchem)

#add_mad_executable(preal preal.cc MADchem)

add_mad_executable(testpg testpg.cc MADchem)

add_mad_executable(testperiodicdft testperiodicdft.cc MADchem)

add_scripted_tests(test_energy.py moldft "short;application")


# removed mcpfit since it is likely only of historic interest
install(TARGETS moldft DESTINATION "${MADNESS_INSTALL_BINDIR}")
