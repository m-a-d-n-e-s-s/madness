# src/apps/moldft

add_executable(DFdriver DFdriver.cc DF.cc fcwf.cc ../chem/gaussian.cc ../chem/NWChem.cc ../chem/polynomial.cc)
target_link_libraries(DFdriver MADchem MADmra)

add_executable(rk rk.cc)
target_link_libraries(rk MADchem MADmra -lgmp -lmpfr)

