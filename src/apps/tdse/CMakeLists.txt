# src/apps/tdse

add_executable(tdse tdse.cc)
target_link_libraries(tdse MADmra)

add_executable(tdse4 tdse4.cc)
target_link_libraries(tdse4 MADmra)

install(TARGETS tdse tdse4 DESTINATION "${MADNESS_INSTALL_BINDIR}")