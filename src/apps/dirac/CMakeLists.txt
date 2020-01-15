# src/apps/moldft

add_executable(DFdriverJ DFdriver.cc DF.cc fcwf.cc gaussian.cc NWChem.cc polynomial.cc)
target_link_libraries(DFdriverJ MADchem MADmra)

add_executable(rkJ rk.cc)
target_link_libraries(rkJ MADchem MADmra -lgmp -lmpfr)

