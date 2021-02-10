# src/apps/tdse

add_mad_executable(tdse tdse.cc MADmra)

add_mad_executable(tdse4 tdse4.cc MADmra)

install(TARGETS tdse tdse4 DESTINATION "${MADNESS_INSTALL_BINDIR}")
