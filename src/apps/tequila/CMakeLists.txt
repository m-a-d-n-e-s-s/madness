# src/apps/moldft

add_executable(tequila_test tequila_test.cpp)
target_link_libraries(tequila_test MADchem ${DISABLEPIE_LINKER_FLAG})

find_library(CNPY_LOCATION cnpy HINTS /home/jsk/install/cnpy/lib/)

add_executable(tequila tequila.cpp)
target_link_libraries(tequila MADchem ${DISABLEPIE_LINKER_FLAG})
target_link_libraries(tequila ${CNPY_LOCATION})

