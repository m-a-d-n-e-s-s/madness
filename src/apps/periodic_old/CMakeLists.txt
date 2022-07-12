# src/apps/hf

set(HARTREEFOCK_HEADERS 
    dft.h eigsolver.h poperator.h util.h lda.h electronicstructureparams.h)
set(HARTREEFOCK_SOURCES eigsolver.cc dft.cc)

# Create HartreeFock library
add_library(HartreeFock 
    ${HARTREEFOCK_SOURCES} ${HARTREEFOCK_HEADERS})
target_link_libraries(HartreeFock PUBLIC MADchem)

# Create executables

add_mad_executable(test_hf test_hf.cc HartreeFock)

add_mad_executable(test_he test_he.cc HartreeFock)

add_mad_executable(test_be test_be.cc HartreeFock)

add_mad_executable(test_hydro test_hydro.cc HartreeFock)

add_mad_executable(test_lattice test_lattice.cc HartreeFock)

add_mad_executable(test_coulomb test_coulomb.cc HartreeFock)

add_mad_executable(test_comm test_comm.cc HartreeFock)

add_mad_executable(testconv testconv.cc HartreeFock)

add_mad_executable(test_xc test_xc.cc HartreeFock)

add_mad_executable(esolver solver_driver.cc mentity.cc HartreeFock)
