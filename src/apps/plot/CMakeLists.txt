# src/apps/moldft

add_mad_executable(plot2cube plot2cube.cpp MADchem)
add_dependencies(applications-madness plot2cube)
add_mad_executable(plot2plane plot2plane.cpp MADchem)
add_dependencies(applications-madness plot2plane)

install(TARGETS plot2cube DESTINATION "${MADNESS_INSTALL_BINDIR}")
install(TARGETS plot2plane DESTINATION "${MADNESS_INSTALL_BINDIR}")
