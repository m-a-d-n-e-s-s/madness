# src/apps/interior_bc

add_mad_executable(embedded_dirichlet embedded_dirichlet.cc embedded_dirichlet MADmra)

#install(TARGETS embedded_dirichlet
#    DESTINATION "${MADNESS_INSTALL_BINDIR}")
