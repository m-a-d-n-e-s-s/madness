# src/apps/moldft

add_mad_executable(plot2cube plot2cube.cpp MADchem)
add_dependencies(applications-madness plot2cube)
install(TARGETS plot2cube DESTINATION "${MADNESS_INSTALL_BINDIR}")

# plot2plane also plots 6D functions (see MADNESS_MRA_DIMENSIONS)
if (6 IN_LIST MADNESS_MRA_DIMENSIONS)
  add_mad_executable(plot2plane plot2plane.cpp MADchem)
  add_dependencies(applications-madness plot2plane)
  install(TARGETS plot2plane DESTINATION "${MADNESS_INSTALL_BINDIR}")
endif()
