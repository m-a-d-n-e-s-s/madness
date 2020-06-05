# src/apps/moldft

add_executable(tequila_test tequila_test.cpp)
target_link_libraries(tequila_test MADchem ${DISABLEPIE_LINKER_FLAG})

find_library(CNPY_LOCATION cnpy HINTS /home/jsk/install/cnpy/lib/)

add_executable(tequila tequila.cpp)
target_link_libraries(tequila MADchem ${DISABLEPIE_LINKER_FLAG})
target_link_libraries(tequila ${CNPY_LOCATION})

add_executable(two_electron_example two_electron_example.cpp)
target_link_libraries(two_electron_example MADchem ${DISABLEPIE_LINKER_FLAG})
target_link_libraries(two_electron_example ${CNPY_LOCATION})

add_executable(compute_two_electron_density_matrix compute_two_electron_density_matrix.cpp)
target_link_libraries(compute_two_electron_density_matrix MADchem ${DISABLEPIE_LINKER_FLAG})
target_link_libraries(compute_two_electron_density_matrix ${CNPY_LOCATION})
