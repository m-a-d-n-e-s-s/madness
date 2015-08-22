# src/apps/hf

set(HARTREEFOCK_HEADERS 
    dft.h eigsolver.h poperator.h util.h lda.h electronicstructureparams.h)
set(HARTREEFOCK_SOURCES eigsolver.cc dft.cc)

# Create HartreeFock library
add_library(HartreeFock 
    ${HARTREEFOCK_SOURCES} ${HARTREEFOCK_HEADERS})
target_link_libraries(HartreeFock PUBLIC MADchem)

# Create executables

add_executable(test_hf test_hf.cc)
target_link_libraries(test_hf HartreeFock)

add_executable(test_he test_he.cc)
target_link_libraries(test_he HartreeFock)

add_executable(test_be test_be.cc)
target_link_libraries(test_he HartreeFock)

add_executable(test_hydro test_hydro.cc)
target_link_libraries(test_hydro HartreeFock)

add_executable(test_lattice test_lattice.cc)
target_link_libraries(test_lattice HartreeFock)

add_executable(test_coulomb test_coulomb.cc)
target_link_libraries(test_coulomb HartreeFock)

add_executable(test_comm test_comm.cc)
target_link_libraries(test_comm HartreeFock)

add_executable(testconv testconv.cc)
target_link_libraries(testconv HartreeFock)

add_executable(test_xc test_xc.cc)
target_link_libraries(test_xc HartreeFock)

add_executable(esolver solver_driver.cc mentity.cc)
target_link_libraries(esolver HartreeFock)
