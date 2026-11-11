# src/apps/tdse

add_mad_executable(tdse tdse.cc MADmra)

if (4 IN_LIST MADNESS_MRA_DIMENSIONS)
  add_mad_executable(tdse4 tdse4.cc MADmra)
endif()

#install(TARGETS tdse tdse4 DESTINATION "${MADNESS_INSTALL_BINDIR}")
