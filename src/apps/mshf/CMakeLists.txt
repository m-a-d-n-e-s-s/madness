# src/apps/mshf

add_executable(mshf mshf.cc densities.cc loadbalance.cc states.cc potential.cc energies.cc output.cc orthonorm.cc iterate.cc ground.cc)
target_link_libraries(mshf MADmra)

install(TARGETS mshf
    DESTINATION "${MADNESS_INSTALL_BINDIR}")
