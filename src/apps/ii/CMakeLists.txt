# src/apps/ii

add_executable(testsystolic testsystolic.cc)
target_link_libraries(testsystolic MADmra)

add_executable(eigen_solver eigen_solver.cc)
target_link_libraries(eigen_solver MADmra)

add_executable(localize_boys localize_boys.cc)
target_link_libraries(localize_boys MADmra)

add_executable(tensor_training tensor_training.cc)
target_link_libraries(tensor_training MADmra)

#install(TARGETS testsystolic eigen_solver localize_boys tensor_training
#    DESTINATION "${MADNESS_INSTALL_BINDIR}")