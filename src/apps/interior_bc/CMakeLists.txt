# src/apps/interior_bc

add_executable(embedded_dirichlet embedded_dirichlet.cc)
target_link_libraries(embedded_dirichlet MADmra)

#install(TARGETS embedded_dirichlet
#    DESTINATION "${MADNESS_INSTALL_BINDIR}")