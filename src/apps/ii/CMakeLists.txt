# src/apps/ii

add_mad_executable(testsystolic testsystolic.cc MADmra)

add_mad_executable(eigen_solver eigen_solver.cc MADmra)

add_mad_executable(localize_boys localize_boys.cc MADmra)

add_mad_executable(tensor_training tensor_training.cc MADmra)

#install(TARGETS testsystolic eigen_solver localize_boys tensor_training
#    DESTINATION "${MADNESS_INSTALL_BINDIR}")
