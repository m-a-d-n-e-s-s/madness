# src/apps/moldft

add_executable(driver5 driver5.cc TDHF2.cc polynomial.cc gaussian.cc NWChem.cc)
target_link_libraries(driver5 MADchem)

